//! Join allele-level outlier calls to locus-level evidence and rank per sample.
//!
//! Neither input can answer the question on its own. The outlier scan knows which alleles are
//! extreme but nothing about which loci matter; the locus table knows which loci are unstable
//! but nothing about any individual. Joining them is what turns a flat list of thousands of
//! outliers into a short, ordered list of candidates for each genome.
//!
//! It reads only the small derived tables, never the combined call file, so re-ranking at a
//! different cutoff is cheap.

use crate::locus_stats::{LocusKey, LocusStatsRow};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

/// Configuration for the prioritisation step
#[derive(Debug, Clone)]
pub struct PrioritizeConfig {
    pub calls: PathBuf,
    pub index: PathBuf,
    pub annotation: Option<PathBuf>,
    pub output: Option<PathBuf>,
    /// Keep loci in the top N percent by variability index
    pub top_index_pct: f64,
    /// Keep at most this many candidates per sample; zero keeps all
    pub top_n: usize,
    /// Drop calls whose allele exceeds this multiple of the locus reference length.
    /// Zero disables the filter.
    pub max_fold: f64,
    /// Require the locus to have at least this many called alleles in the reference
    pub min_locus_alleles: usize,
    /// Flag samples producing more than this multiple of the cohort median call count
    pub sample_outlier_factor: f64,
    pub funnel: Option<PathBuf>,
}

/// One outlier call read back from the scan
#[derive(Clone, Debug)]
struct Call {
    key: LocusKey,
    chrom: String,
    motif: String,
    sample: String,
    event: String,
    allele: String,
    length: f64,
    alleles_passing: usize,
    score: f64,
    significance: f64,
    evidence: String,
}

/// A call that survived filtering, carrying the locus evidence it was joined to
struct Candidate {
    call: Call,
    /// Variability index of the locus, NaN for a dropout retained without a score
    index: f64,
    /// Called alleles behind the locus statistics
    locus_alleles: usize,
}

/// Counts at each stage of filtering, reported so a user can see whether their thresholds do
/// anything. This is the single most informative diagnostic the pipeline produces.
#[derive(Default, Debug, Clone)]
pub struct Funnel {
    pub calls_in: usize,
    pub dropouts_in: usize,
    pub unmapped: usize,
    pub failed_locus_alleles: usize,
    pub failed_index: usize,
    pub failed_fold: usize,
    pub kept: usize,
    pub kept_after_top_n: usize,
    pub samples: usize,
}

impl Funnel {
    fn write<W: Write>(&self, w: &mut W, per_sample_median: f64, allele_sum: usize) {
        writeln!(w, "stage\tcalls\tper_sample_median").unwrap();
        writeln!(w, "outlier_calls_in\t{}\t", self.calls_in).unwrap();
        writeln!(w, "dropout_events_in\t{}\t", self.dropouts_in).unwrap();
        writeln!(w, "unmapped_to_index\t{}\t", self.unmapped).unwrap();
        writeln!(w, "removed_locus_call_rate\t{}\t", self.failed_locus_alleles).unwrap();
        writeln!(w, "removed_implausible_length\t{}\t", self.failed_fold).unwrap();
        writeln!(w, "removed_by_index_filter\t{}\t", self.failed_index).unwrap();
        writeln!(w, "kept\t{}\t{:.0}", self.kept, per_sample_median).unwrap();
        writeln!(w, "final_after_top_n\t{}\t", self.kept_after_top_n).unwrap();
        writeln!(
            w,
            "# per-allele sum (the published figures count alleles, not samples): {}",
            allele_sum
        )
        .unwrap();
        writeln!(w, "# samples: {}", self.samples).unwrap();
    }
}

/// Read the scan output
fn read_calls(path: &Path) -> (Vec<Call>, usize) {
    let reader = crate::utils::reader(&path.to_string_lossy());
    let mut calls = Vec::new();
    let mut dropouts = 0usize;
    let mut columns: HashMap<String, usize> = HashMap::new();

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("ERROR: Failed reading {}: {}", path.display(), e);
                std::process::exit(1);
            }
        };
        if line.starts_with('#') {
            continue;
        }
        if columns.is_empty() {
            for (i, name) in line.split('\t').enumerate() {
                columns.insert(name.to_string(), i);
            }
            for required in [
                "chromosome",
                "begin",
                "end",
                "sample",
                "event",
                "length",
                "score",
            ]
            .iter()
            {
                if !columns.contains_key(*required) {
                    eprintln!(
                        "ERROR: {} is missing the '{}' column; it does not look like \
                         inquiSTR outlier output.",
                        path.display(),
                        required
                    );
                    std::process::exit(1);
                }
            }
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let get = |name: &str| -> &str { columns.get(name).and_then(|i| f.get(*i)).unwrap_or(&"") };
        let event = get("event").to_string();
        if event == "dropout" {
            dropouts += 1;
        }
        let chrom = get("chromosome").to_string();
        let (begin, end) = match (get("begin").parse::<u32>(), get("end").parse::<u32>()) {
            (Ok(b), Ok(e)) => (b, e),
            _ => continue,
        };
        calls.push(Call {
            key: LocusKey::new(&chrom, begin, end),
            chrom,
            motif: get("motif").to_string(),
            sample: get("sample").to_string(),
            event,
            allele: get("allele").to_string(),
            length: get("length").parse().unwrap_or(f64::NAN),
            alleles_passing: get("alleles_passing").parse().unwrap_or(1),
            score: get("score").parse().unwrap_or(f64::NAN),
            significance: get("significance").parse().unwrap_or(f64::NAN),
            evidence: get("evidence").to_string(),
        });
    }
    (calls, dropouts)
}

/// Read an annotation mapping produced by `annotate`, keyed by locus
fn read_annotation(path: &Path) -> HashMap<LocusKey, (f64, String)> {
    let (header, reader) = crate::locus_stats::open(path);
    let mut map = HashMap::new();
    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };
        if let Some(row) = crate::locus_stats::parse_row(&header, &line)
            && let Some(idx) = row.index
        {
            let column = header
                .index_column
                .clone()
                .unwrap_or_else(|| "external_index".to_string());
            map.insert(row.key.clone(), (idx, column));
        }
    }
    eprintln!("Loaded {} annotated loci from {}", map.len(), path.display());
    map
}

/// Rank candidates within each sample and write the shortlist
pub fn prioritize(config: PrioritizeConfig) {
    let (calls, dropouts_in) = read_calls(&config.calls);
    eprintln!("Read {} rows from {}", calls.len(), config.calls.display());
    if calls.is_empty() {
        eprintln!("No calls to prioritise.");
        return;
    }

    let wanted: HashSet<LocusKey> = calls.iter().map(|c| c.key.clone()).collect();
    eprintln!("Looking up {} distinct loci in {}", wanted.len(), config.index.display());
    let (stats_header, stats) = crate::locus_stats::load_subset(&config.index, &wanted);
    eprintln!("Matched {} of {} loci", stats.len(), wanted.len());
    if stats.is_empty() {
        eprintln!(
            "ERROR: no locus in the call file was found in the index. Check that both were \
             produced from the same catalog."
        );
        std::process::exit(1);
    }
    let match_rate = stats.len() as f64 / wanted.len() as f64;
    if match_rate < 0.9 {
        eprintln!(
            "Warning: only {:.1}% of loci matched the index. Rankings below cover the matched \
             loci only.",
            100.0 * match_rate
        );
    }
    if stats_header.index_column.is_none() {
        eprintln!(
            "ERROR: {} carries no variability index column, so loci cannot be ranked by it. \
             Use a table produced by `inquiSTR varindex`.",
            config.index.display()
        );
        std::process::exit(1);
    }

    let annotation = config.annotation.as_deref().map(read_annotation);

    // The index cutoff is a percentile of the index itself, which is already a rank in [0, 1],
    // so the threshold is a direct comparison rather than a re-ranking.
    let index_cutoff = 1.0 - config.top_index_pct / 100.0;
    let mut funnel =
        Funnel { calls_in: calls.len() - dropouts_in, dropouts_in, ..Funnel::default() };

    let mut kept: Vec<Candidate> = Vec::new();
    for call in calls.into_iter() {
        if call.event == "dropout" {
            // Dropout events carry no length to rank, but are retained so that a failed
            // genotype at a high-index locus is visible rather than silently absent.
            match stats.get(&call.key) {
                Some(row) if row.index.is_some_and(|i| i >= index_cutoff) => {
                    kept.push(Candidate { call, index: f64::NAN, locus_alleles: row.n });
                }
                _ => {}
            }
            continue;
        }
        let Some(row): Option<&LocusStatsRow> = stats.get(&call.key) else {
            funnel.unmapped += 1;
            continue;
        };
        if row.n < config.min_locus_alleles {
            funnel.failed_locus_alleles += 1;
            continue;
        }
        if config.max_fold > 0.0 {
            let reflen = (call.key.end.saturating_sub(call.key.begin)).max(1) as f64;
            if call.length / reflen > config.max_fold {
                funnel.failed_fold += 1;
                continue;
            }
        }
        match row.index {
            Some(index) if index >= index_cutoff => {
                kept.push(Candidate { call, index, locus_alleles: row.n });
            }
            _ => funnel.failed_index += 1,
        }
    }
    funnel.kept = kept.len();

    // Rank within each sample. The score is in repeat units and comparable across loci, which
    // is what makes ordering a whole genome meaningful.
    let mut by_sample: HashMap<String, Vec<Candidate>> = HashMap::new();
    for item in kept {
        by_sample
            .entry(item.call.sample.clone())
            .or_default()
            .push(item);
    }
    funnel.samples = by_sample.len();

    let mut counts: Vec<usize> = by_sample.values().map(|v| v.len()).collect();
    counts.sort_unstable();
    let median = if counts.is_empty() {
        0.0
    } else if counts.len().is_multiple_of(2) {
        (counts[counts.len() / 2 - 1] + counts[counts.len() / 2]) as f64 / 2.0
    } else {
        counts[counts.len() / 2] as f64
    };

    let mut out: Box<dyn Write> = match &config.output {
        Some(p) => {
            let f = std::fs::File::create(p).unwrap_or_else(|e| {
                eprintln!("ERROR: Failed to create {}: {}", p.display(), e);
                std::process::exit(1);
            });
            if p.extension().is_some_and(|e| e == "gz") {
                Box::new(std::io::BufWriter::new(flate2::write::GzEncoder::new(
                    f,
                    flate2::Compression::default(),
                )))
            } else {
                Box::new(std::io::BufWriter::new(f))
            }
        }
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };

    writeln!(out, "# file_type=prioritized_calls").unwrap();
    writeln!(out, "# version={}", crate::VERSION).unwrap();
    writeln!(out, "# command=prioritize").unwrap();
    writeln!(out, "# calls={}", config.calls.display()).unwrap();
    writeln!(out, "# index={}", config.index.display()).unwrap();
    writeln!(
        out,
        "# index_column={}",
        stats_header
            .index_column
            .clone()
            .unwrap_or_else(|| "unknown".to_string())
    )
    .unwrap();
    writeln!(out, "# top_index_pct={}", config.top_index_pct).unwrap();
    writeln!(out, "# index_cutoff={}", index_cutoff).unwrap();
    writeln!(out, "# max_fold={}", config.max_fold).unwrap();
    writeln!(out, "# samples={}", funnel.samples).unwrap();
    writeln!(
        out,
        "sample\trank\tchromosome\tbegin\tend\tmotif\tevent\tallele\tlength\tscore\t\
         significance\tlocus_index\tlocus_alleles\texternal_index\tsample_flag\tevidence"
    )
    .unwrap();

    let mut sample_names: Vec<&String> = by_sample.keys().collect();
    sample_names.sort();
    let mut final_rows = 0usize;

    for name in sample_names {
        let items = &by_sample[name];
        let flagged = config.sample_outlier_factor > 0.0
            && median > 0.0
            && items.len() as f64 > config.sample_outlier_factor * median;
        let mut sorted: Vec<&Candidate> = items.iter().collect();
        // Dropouts carry no score. NaN must be ordered explicitly rather than folded into
        // `Ordering::Equal`: doing that breaks transitivity, and the sort then silently
        // scrambles every rank rather than merely misplacing the scoreless rows.
        sorted.sort_by(|a, b| {
            let (x, y) = (a.call.score, b.call.score);
            match (x.is_nan(), y.is_nan()) {
                (true, true) => a.call.key.cmp(&b.call.key),
                (true, false) => std::cmp::Ordering::Greater,
                (false, true) => std::cmp::Ordering::Less,
                (false, false) => y
                    .partial_cmp(&x)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.call.key.cmp(&b.call.key)),
            }
        });
        let limit = if config.top_n == 0 {
            sorted.len()
        } else {
            config.top_n.min(sorted.len())
        };
        for (rank, candidate) in sorted.iter().take(limit).enumerate() {
            let call = &candidate.call;
            let index = candidate.index;
            let n_alleles = candidate.locus_alleles;
            let external = annotation
                .as_ref()
                .and_then(|a| a.get(&call.key))
                .map(|(v, _)| format!("{:.6}", v))
                .unwrap_or_else(|| ".".to_string());
            let index_str = if index.is_nan() {
                ".".to_string()
            } else {
                format!("{:.6}", index)
            };
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                call.sample,
                rank + 1,
                call.chrom,
                call.key.begin,
                call.key.end,
                call.motif,
                call.event,
                call.allele,
                if call.length.is_nan() {
                    "NaN".into()
                } else {
                    format!("{}", call.length)
                },
                if call.score.is_nan() {
                    "NaN".into()
                } else {
                    format!("{:.6}", call.score)
                },
                if call.significance.is_nan() {
                    "NaN".into()
                } else {
                    format!("{:.6e}", call.significance)
                },
                index_str,
                n_alleles,
                external,
                if flagged { "excess_calls" } else { "." },
                call.evidence
            )
            .unwrap();
            final_rows += 1;
        }
    }
    out.flush().unwrap();
    funnel.kept_after_top_n = final_rows;

    // The scan reports a sample once per locus but records how many of its alleles passed,
    // so the published per-allele counting can be reproduced rather than approximated.
    // Dropout rows carry no passing allele, and `read_calls` already defaults the column to 1
    // for files written before it existed, so clamping here would only inflate the count.
    let allele_sum: usize = by_sample
        .values()
        .flatten()
        .map(|c| c.call.alleles_passing)
        .sum();
    if let Some(path) = &config.funnel {
        let mut f = std::fs::File::create(path).unwrap_or_else(|e| {
            eprintln!("ERROR: Failed to create {}: {}", path.display(), e);
            std::process::exit(1);
        });
        funnel.write(&mut f, median, allele_sum);
        eprintln!("Wrote funnel report to {}", path.display());
    }
    let mut stderr = std::io::stderr();
    writeln!(stderr, "\nFiltering funnel:").unwrap();
    funnel.write(&mut stderr, median, allele_sum);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    const CALLS: &str = "# file_type=outlier_calls\n\
        chromosome\tbegin\tend\tmotif\tsample\tevent\tallele\tlength\tscore\tsignificance\tmethod\tevidence\n\
        chr1\t100\t130\tCAG\tsA\toutlier\tH1\t300\t90.0\t1.0e-3\tpercentile\tref_n=10\n\
        chr1\t200\t230\tCAG\tsA\toutlier\tH2\t150\t40.0\t1.0e-3\tpercentile\tref_n=10\n\
        chr1\t300\t330\tCAG\tsB\toutlier\tH1\t333\t101.0\t1.0e-3\tpercentile\tref_n=10\n";

    // locus 100 is highly variable, 200 is not, 300 is highly variable
    const INDEX: &str = "# file_type=locus_stats\n# samples=5\n# index_column=lvi_length_internal\n\
        chromosome\tbegin\tend\tmotif\tN\tstd\tlvi_length_internal\n\
        chr1\t100\t130\tCAG\t10\t5.0\t0.99\n\
        chr1\t200\t230\tCAG\t10\t1.0\t0.10\n\
        chr1\t300\t330\tCAG\t10\t9.0\t0.999\n";

    fn run(cfg: PrioritizeConfig) -> String {
        let out = cfg.output.clone().unwrap();
        prioritize(cfg);
        std::fs::read_to_string(out).unwrap()
    }

    fn base(calls: &Path, index: &Path, out: &Path) -> PrioritizeConfig {
        PrioritizeConfig {
            calls: calls.to_path_buf(),
            index: index.to_path_buf(),
            annotation: None,
            output: Some(out.to_path_buf()),
            top_index_pct: 5.0,
            top_n: 0,
            max_fold: 0.0,
            min_locus_alleles: 0,
            sample_outlier_factor: 0.0,
            funnel: None,
        }
    }

    #[test]
    fn test_index_filter_removes_uninteresting_loci() {
        let c = tmp(CALLS);
        let i = tmp(INDEX);
        let o = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        let out = run(base(c.path(), i.path(), o.path()));
        // The low-index locus at 200 must not survive a top-5% filter
        assert!(out.contains("\t100\t130\t"), "high-index locus should be kept");
        assert!(!out.contains("\t200\t230\t"), "low-index locus should be filtered out");
        assert!(out.contains("\t300\t330\t"));
    }

    #[test]
    fn test_ranking_is_per_sample_and_by_score() {
        let calls = "# file_type=outlier_calls\n\
            chromosome\tbegin\tend\tmotif\tsample\tevent\tallele\tlength\tscore\tsignificance\tmethod\tevidence\n\
            chr1\t100\t130\tCAG\tsA\toutlier\tH1\t300\t10.0\t1.0e-3\tpercentile\t.\n\
            chr1\t300\t330\tCAG\tsA\toutlier\tH1\t333\t99.0\t1.0e-3\tpercentile\t.\n";
        let c = tmp(calls);
        let i = tmp(INDEX);
        let o = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        let out = run(base(c.path(), i.path(), o.path()));
        let data: Vec<&str> = out.lines().filter(|l| l.starts_with("sA\t")).collect();
        assert_eq!(data.len(), 2);
        assert!(data[0].contains("\t1\t"), "rank 1 first");
        assert!(data[0].contains("\t300\t330\t"), "the higher score must rank first");
    }

    #[test]
    fn test_implausible_lengths_are_removed_and_counted() {
        // 3000 bp on a 30 bp locus is a hundredfold expansion
        let calls = "# file_type=outlier_calls\n\
            chromosome\tbegin\tend\tmotif\tsample\tevent\tallele\tlength\tscore\tsignificance\tmethod\tevidence\n\
            chr1\t100\t130\tCAG\tsA\toutlier\tH1\t3000\t990.0\t1.0e-3\tpercentile\t.\n";
        let c = tmp(calls);
        let i = tmp(INDEX);
        let o = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        let mut cfg = base(c.path(), i.path(), o.path());
        cfg.max_fold = 50.0;
        let out = run(cfg);
        assert!(!out.contains("\tsA\t") && !out.lines().any(|l| l.starts_with("sA\t")));
    }

    #[test]
    fn test_dropout_at_a_high_index_locus_survives() {
        let calls = "# file_type=outlier_calls\n\
            chromosome\tbegin\tend\tmotif\tsample\tevent\tallele\tlength\tscore\tsignificance\tmethod\tevidence\n\
            chr1\t100\t130\tCAG\tsA\tdropout\tone\tNaN\tNaN\tNaN\tpercentile\t.\n\
            chr1\t200\t230\tCAG\tsB\tdropout\tone\tNaN\tNaN\tNaN\tpercentile\t.\n";
        let c = tmp(calls);
        let i = tmp(INDEX);
        let o = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        let out = run(base(c.path(), i.path(), o.path()));
        assert!(
            out.lines()
                .any(|l| l.starts_with("sA\t") && l.contains("dropout")),
            "a failed genotype at a high-index locus must remain visible"
        );
        assert!(
            !out.lines().any(|l| l.starts_with("sB\t")),
            "a dropout at an uninteresting locus is not a candidate"
        );
    }

    #[test]
    fn test_unmatched_loci_are_counted_not_silently_dropped() {
        let calls = "# file_type=outlier_calls\n\
            chromosome\tbegin\tend\tmotif\tsample\tevent\tallele\tlength\tscore\tsignificance\tmethod\tevidence\n\
            chr9\t999\t1099\tCAG\tsA\toutlier\tH1\t300\t90.0\t1.0e-3\tpercentile\t.\n\
            chr1\t100\t130\tCAG\tsA\toutlier\tH1\t300\t90.0\t1.0e-3\tpercentile\t.\n";
        let c = tmp(calls);
        let i = tmp(INDEX);
        let o = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        let out = run(base(c.path(), i.path(), o.path()));
        assert!(!out.contains("chr9"));
        assert!(out.contains("\t100\t130\t"));
    }

    #[test]
    fn test_top_n_truncates_per_sample() {
        let c = tmp(CALLS);
        let i = tmp(INDEX);
        let o = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        let mut cfg = base(c.path(), i.path(), o.path());
        cfg.top_index_pct = 100.0;
        cfg.top_n = 1;
        let out = run(cfg);
        assert_eq!(out.lines().filter(|l| l.starts_with("sA\t")).count(), 1);
    }

    #[test]
    fn test_chromosome_naming_differences_still_join() {
        // Calls use chr-prefixed names, the index does not
        let index = "# file_type=locus_stats\n# index_column=lvi_length_internal\n\
            TRID\tN\tstd\tlvi_length_internal\n\
            1-100-130-CAG\t10\t5.0\t0.99\n";
        let c = tmp(CALLS);
        let i = tmp(index);
        let o = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        let out = run(base(c.path(), i.path(), o.path()));
        assert!(out.contains("\t100\t130\t"), "chr1 must join to 1");
    }
}
