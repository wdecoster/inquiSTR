//! Genome-wide outlier scan over a combined call file.
//!
//! Tests every allele of every sample against a reference distribution at its locus and emits
//! a scored row for each allele extreme enough to pass a deliberately permissive gate. It
//! applies no locus-level knowledge: ranking, index filtering and quality control belong to
//! the prioritisation step, so that re-ranking never requires rescanning the combined file.
//!
//! The default reference is the cohort itself with the tested sample removed. That exclusion
//! is what makes the scan able to see a carrier at all: an individual holding the longest
//! allele at a locus would otherwise be compared against their own allele.

use crate::outlier_methods::{
    Evidence, LocusReference, Method, OutlierResult, dbscan_result, score_exceedance, score_moment,
};
use clap::ValueEnum;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

/// Which allele of a sample to test
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum Zygosity {
    /// Test both alleles independently and report the more extreme, naming which it was
    Both,
    /// Test only the longer allele, appropriate for dominant models
    Max,
    /// Test only the shorter allele. Biallelic disorders such as RFC1 and FXN are missed
    /// entirely by `max`, because there the shorter allele being long is the signal.
    Min,
}

/// Configuration for the outlier scan
#[derive(Debug, Clone)]
pub struct ScanConfig {
    pub combined: PathBuf,
    pub output: Option<PathBuf>,
    pub method: Method,
    pub threads: usize,
    pub zygosity: Zygosity,
    /// Gate for the percentile method: report an allele when at most this many reference
    /// alleles are at least as long. Zero means "longer than everything in the reference".
    pub max_exceedance: usize,
    /// Gate for the moment methods, in standard deviations
    pub zscore: f64,
    /// Minimum score, in repeat units, for a result to be reported
    pub min_units: f64,
    /// Skip loci where no allele reaches this length change relative to the reference
    pub minsize: u32,
    /// Restrict reporting to these samples
    pub samples: Option<Vec<String>>,
    /// Groups of related samples; all members are excluded from each other's reference
    pub related: Option<PathBuf>,
    /// Emit a row for uncalled genotypes as their own event type
    pub report_dropout: bool,
    /// Skip loci called in fewer than this fraction of alleles
    pub min_call_rate: f64,
    /// Restrict the scan to loci overlapping this BED. Required in practice for
    /// `report_dropout`, which would otherwise emit a row per missing genotype genome-wide.
    pub regions: Option<PathBuf>,
}

/// Intervals to restrict the scan to, indexed per chromosome for bounded lookup
struct Regions {
    by_chrom: HashMap<String, (Vec<(u32, u32)>, u32)>,
}

impl Regions {
    fn read(path: &Path) -> Self {
        let reader = crate::utils::reader(&path.to_string_lossy());
        let mut by_chrom: HashMap<String, (Vec<(u32, u32)>, u32)> = HashMap::new();
        for line in reader.lines().map_while(Result::ok) {
            if line.starts_with('#') || line.starts_with("track") || line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 3 {
                continue;
            }
            if let (Ok(b), Ok(e)) = (f[1].parse::<u32>(), f[2].parse::<u32>()) {
                let slot = by_chrom
                    .entry(crate::locus_stats::normalise_chrom(f[0]))
                    .or_insert_with(|| (Vec::new(), 0));
                slot.0.push((b, e));
                slot.1 = slot.1.max(e.saturating_sub(b));
            }
        }
        let mut n = 0;
        for (v, _) in by_chrom.values_mut() {
            v.sort_unstable();
            n += v.len();
        }
        eprintln!("Restricting the scan to {} regions from {}", n, path.display());
        Regions { by_chrom }
    }

    fn overlaps(&self, chrom: &str, begin: u32, end: u32) -> bool {
        let Some((v, max_len)) = self
            .by_chrom
            .get(&crate::locus_stats::normalise_chrom(chrom))
        else {
            return false;
        };
        let from = v.partition_point(|&(b, _)| b < begin.saturating_sub(*max_len));
        v[from..]
            .iter()
            .take_while(|&&(b, _)| b < end)
            .any(|&(b, e)| e > begin && b < end)
    }
}

/// One reported observation
struct Row {
    chrom: String,
    begin: u32,
    end: u32,
    motif: String,
    sample: String,
    event: &'static str,
    allele: &'static str,
    length: f64,
    /// How many of the sample's alleles passed the gate at this locus. A sample is reported
    /// once per locus, so this is what makes the published per-allele counts reconstructable.
    alleles_passing: usize,
    result: Option<OutlierResult>,
}

/// Parse related-sample groups. Each line lists the members of one group, separated by tabs
/// or commas. Returns a map from sample name to group id.
fn read_related(path: &Path) -> HashMap<String, usize> {
    let mut map = HashMap::new();
    let reader = crate::utils::reader(&path.to_string_lossy());
    for (i, line) in reader.lines().enumerate() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("ERROR: Failed reading {}: {}", path.display(), e);
                std::process::exit(1);
            }
        };
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        for name in line.split(['\t', ',']) {
            let name = name.trim();
            if !name.is_empty() {
                map.insert(name.to_string(), i);
            }
        }
    }
    eprintln!("Loaded {} samples in related groups from {}", map.len(), path.display());
    map
}

/// Strip the _H1/_H2 suffix from a column name
fn base_name(column: &str) -> &str {
    column
        .strip_suffix("_H1")
        .or_else(|| column.strip_suffix("_H2"))
        .unwrap_or(column)
}

/// Genotypes of one sample at one locus, as absolute lengths
#[derive(Clone, Copy)]
struct SampleAlleles {
    h1: Option<f64>,
    h2: Option<f64>,
}

impl SampleAlleles {
    fn present(&self) -> Vec<f64> {
        self.h1.into_iter().chain(self.h2).collect()
    }
    fn missing(&self) -> bool {
        self.h1.is_none() || self.h2.is_none()
    }
}

/// Everything the scan needs about one locus
struct Locus {
    chrom: String,
    begin: u32,
    end: u32,
    motif: String,
    motif_len: f64,
    per_sample: Vec<SampleAlleles>,
    sorted: Vec<f64>,
    sum: f64,
    sumsq: f64,
    median: f64,
    mad: f64,
    called: usize,
    total: usize,
}

/// Parse one combined-file line into a locus, converting genotypes to absolute lengths.
/// `NaN` is missing and is never coerced to zero.
fn parse_locus(line: &str, n_samples: usize) -> Option<Locus> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < crate::filetype::STR_FIXED_COLUMNS + 2 {
        return None;
    }
    let begin: u32 = f[1].parse().ok()?;
    let end: u32 = f[2].parse().ok()?;
    let motif = crate::varindex::parse_motif(f[3])
        .unwrap_or(".")
        .to_string();
    let motif_len = if motif == "." {
        1.0
    } else {
        motif.len() as f64
    };
    let reflen = end.saturating_sub(begin) as f64;

    let parse = |s: &str| -> Option<f64> {
        if s.is_empty() || s.eq_ignore_ascii_case("nan") || s == "." {
            return None;
        }
        let v: f64 = s.parse().ok()?;
        if v.is_finite() {
            Some((reflen + v).max(0.0))
        } else {
            None
        }
    };

    let mut per_sample = Vec::with_capacity(n_samples);
    let mut sorted = Vec::with_capacity(n_samples * 2);
    let cols = &f[crate::filetype::STR_FIXED_COLUMNS..];
    for pair in cols.chunks(2) {
        let h1 = parse(pair[0]);
        let h2 = pair.get(1).and_then(|s| parse(s));
        if let Some(v) = h1 {
            sorted.push(v);
        }
        if let Some(v) = h2 {
            sorted.push(v);
        }
        per_sample.push(SampleAlleles { h1, h2 });
    }

    let total = per_sample.len() * 2;
    let called = sorted.len();
    sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let sum: f64 = sorted.iter().sum();
    let sumsq: f64 = sorted.iter().map(|v| v * v).sum();

    let median = if sorted.is_empty() {
        f64::NAN
    } else if sorted.len().is_multiple_of(2) {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };
    let mad = if sorted.is_empty() {
        f64::NAN
    } else {
        let mut d: Vec<f64> = sorted.iter().map(|v| (v - median).abs()).collect();
        d.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        if d.len().is_multiple_of(2) {
            (d[d.len() / 2 - 1] + d[d.len() / 2]) / 2.0
        } else {
            d[d.len() / 2]
        }
    };

    Some(Locus {
        chrom: f[0].to_string(),
        begin,
        end,
        motif,
        motif_len,
        per_sample,
        sorted,
        sum,
        sumsq,
        median,
        mad,
        called,
        total,
    })
}

/// Nearest value in `sorted` that is not the query itself, used to give DBSCAN noise points a
/// magnitude rather than only a label.
fn nearest_other(sorted: &[f64], query: f64, eps: f64) -> f64 {
    let mut best = f64::NAN;
    let mut best_d = f64::INFINITY;
    for &v in sorted {
        let d = (v - query).abs();
        if d > 0.0 && d < best_d && d > eps {
            best_d = d;
            best = v;
        }
    }
    best
}

/// Mode of a sorted slice, used by the DBSCAN radius heuristic
fn mode_of(sorted: &[f64]) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let (mut best, mut best_c) = (sorted[0], 1usize);
    let (mut cur, mut cur_c) = (sorted[0], 1usize);
    for &v in &sorted[1..] {
        if v == cur {
            cur_c += 1;
        } else {
            if cur_c > best_c {
                best_c = cur_c;
                best = cur;
            }
            cur = v;
            cur_c = 1;
        }
    }
    if cur_c > best_c {
        best = cur;
    }
    best
}

/// Collect one allele per sample according to the zygosity mode, ascending.
/// Used to build a reference drawn from the same allele slot that is being tested.
fn collect_by_zygosity(per_sample: &[SampleAlleles], zygosity: Zygosity) -> Vec<f64> {
    let mut v: Vec<f64> = per_sample
        .iter()
        .filter_map(|a| match (a.h1, a.h2) {
            (Some(x), Some(y)) => Some(if zygosity == Zygosity::Min {
                x.min(y)
            } else {
                x.max(y)
            }),
            (Some(x), None) | (None, Some(x)) => Some(x),
            (None, None) => None,
        })
        .collect();
    v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v
}

/// Median and median absolute deviation of the reference, with `excluded` removed.
///
/// Robust statistics resist a single extreme value, but a large related group is not a single
/// value, so the exclusion list is applied here as it is for every other method.
fn robust_center_scale(reference: &[f64], excluded: &[f64]) -> (f64, f64) {
    let mut kept: Vec<f64> = Vec::with_capacity(reference.len());
    let mut drop: Vec<f64> = excluded.to_vec();
    for &v in reference {
        if let Some(pos) = drop.iter().position(|&d| d == v) {
            drop.swap_remove(pos);
            continue;
        }
        kept.push(v);
    }
    if kept.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let median = |v: &mut Vec<f64>| -> f64 {
        v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        if v.len().is_multiple_of(2) {
            (v[v.len() / 2 - 1] + v[v.len() / 2]) / 2.0
        } else {
            v[v.len() / 2]
        }
    };
    let med = median(&mut kept);
    let mut devs: Vec<f64> = kept.iter().map(|v| (v - med).abs()).collect();
    (med, median(&mut devs))
}

/// Score one locus and return the rows that pass the gate
fn scan_locus(
    locus: &Locus,
    config: &ScanConfig,
    sample_names: &[String],
    keep: &Option<Vec<bool>>,
    groups: &Option<HashMap<String, usize>>,
) -> Vec<Row> {
    let mut rows = Vec::new();
    if locus.called == 0 {
        return rows;
    }
    if (locus.called as f64 / locus.total as f64) < config.min_call_rate {
        return rows;
    }

    // The reference must be drawn from the same allele slot being tested, or the comparison is
    // not like with like. Under `min`, a biallelic carrier's shorter allele would otherwise be
    // measured against every sample's *longer* allele and could never stand out — which is
    // exactly the recessive signal the mode exists to find.
    let pooled: Vec<f64>;
    let reference_values: &[f64] = match config.zygosity {
        Zygosity::Both => &locus.sorted,
        Zygosity::Max | Zygosity::Min => {
            pooled = collect_by_zygosity(&locus.per_sample, config.zygosity);
            &pooled
        }
    };
    let sum: f64 = match config.zygosity {
        Zygosity::Both => locus.sum,
        _ => reference_values.iter().sum(),
    };
    let sumsq: f64 = match config.zygosity {
        Zygosity::Both => locus.sumsq,
        _ => reference_values.iter().map(|v| v * v).sum(),
    };
    if reference_values.is_empty() {
        return rows;
    }
    let reference =
        LocusReference { sorted: reference_values, sum, sumsq, motif_len: locus.motif_len };

    // DBSCAN classifies all points jointly, so its clustering is computed once per locus.
    let (eps, dbscan_noise) = if config.method == Method::Dbscan {
        let eps = (2.0 * mode_of(&locus.sorted)).max(10.0);
        let min_points = (locus.sorted.len() as f64).log2().max(2.0) as usize;
        let points: Vec<Vec<f64>> = locus.sorted.iter().map(|v| vec![*v]).collect();
        let model = dbscan::Model::new(eps, min_points);
        let classes = model.run(&points);
        let noise: Vec<f64> = classes
            .iter()
            .zip(locus.sorted.iter())
            .filter(|(c, _)| matches!(c, dbscan::Classification::Noise))
            .map(|(_, v)| *v)
            .collect();
        (eps, Some(noise))
    } else {
        (0.0, None)
    };
    let min_points = (locus.sorted.len() as f64).log2().max(2.0) as usize;

    let p_gate = crate::outlier_methods::normal_sf_public(config.zscore);

    for (idx, alleles) in locus.per_sample.iter().enumerate() {
        // A data line with more genotype columns than the header has is malformed; skip the
        // surplus rather than indexing past the sample list.
        let Some(name) = sample_names.get(idx) else {
            break;
        };
        if let Some(keep) = keep
            && !keep[idx]
        {
            continue;
        }

        if config.report_dropout && alleles.missing() {
            rows.push(Row {
                chrom: locus.chrom.clone(),
                begin: locus.begin,
                end: locus.end,
                motif: locus.motif.clone(),
                sample: name.clone(),
                event: "dropout",
                allele: if alleles.h1.is_none() && alleles.h2.is_none() {
                    "both"
                } else {
                    "one"
                },
                length: f64::NAN,
                alleles_passing: 0,
                result: None,
            });
        }

        // Alleles removed from the reference: the sample's own contribution to whichever
        // slot the reference was built from, plus that of any relatives.
        let own_contribution = |a: &SampleAlleles| -> Vec<f64> {
            match config.zygosity {
                Zygosity::Both => a.present(),
                Zygosity::Max | Zygosity::Min => match (a.h1, a.h2) {
                    (Some(x), Some(y)) => vec![if config.zygosity == Zygosity::Min {
                        x.min(y)
                    } else {
                        x.max(y)
                    }],
                    (Some(x), None) | (None, Some(x)) => vec![x],
                    (None, None) => vec![],
                },
            }
        };
        let mut excluded = own_contribution(alleles);
        if let Some(groups) = groups
            && let Some(g) = groups.get(base_name(name))
        {
            for (j, other) in locus.per_sample.iter().enumerate() {
                // Same guard as the outer loop: a ragged line can hold more genotype pairs
                // than the header names samples.
                let Some(other_name) = sample_names.get(j) else {
                    break;
                };
                if j != idx && groups.get(base_name(other_name)) == Some(g) {
                    // Relatives contribute to the reference through the same allele slot as
                    // everyone else, so they must be removed through that slot too. Removing
                    // both their alleles from a one-allele-per-sample reference deletes
                    // unrelated samples' values by coincidental equality.
                    excluded.extend(own_contribution(other));
                }
            }
        }

        let candidates: Vec<(&'static str, f64)> = match config.zygosity {
            Zygosity::Both => {
                let mut v = Vec::new();
                if let Some(x) = alleles.h1 {
                    v.push(("H1", x));
                }
                if let Some(x) = alleles.h2 {
                    v.push(("H2", x));
                }
                v
            }
            Zygosity::Max => match (alleles.h1, alleles.h2) {
                (Some(a), Some(b)) => {
                    vec![if a >= b { ("H1", a) } else { ("H2", b) }]
                }
                (Some(a), None) => vec![("H1", a)],
                (None, Some(b)) => vec![("H2", b)],
                _ => vec![],
            },
            Zygosity::Min => match (alleles.h1, alleles.h2) {
                (Some(a), Some(b)) => {
                    vec![if a <= b { ("H1", a) } else { ("H2", b) }]
                }
                (Some(a), None) => vec![("H1", a)],
                (None, Some(b)) => vec![("H2", b)],
                _ => vec![],
            },
        };

        // Both alleles are tested, but a sample is reported once per locus: the more extreme
        // allele stands for the sample and is named in the output.
        let mut best: Option<(&'static str, f64, OutlierResult)> = None;
        let mut passing = 0usize;
        for (label, value) in candidates {
            let result = match config.method {
                Method::Percentile => score_exceedance(&reference, value, &excluded),
                Method::Zscore => {
                    score_moment(&reference, value, &excluded, false, locus.median, locus.mad)
                }
                Method::Robustz => {
                    // Centre and scale must come from the same allele slot the reference was
                    // built from, and with relatives removed. Using the pooled both-allele
                    // median under `--zygosity min` would measure a recessive signal against
                    // the wrong null, and would ignore `--related` entirely.
                    let (med, mad) = robust_center_scale(reference_values, &excluded);
                    score_moment(&reference, value, &excluded, true, med, mad)
                }
                Method::Dbscan => {
                    let noise = dbscan_noise.as_ref().is_some_and(|n| n.contains(&value));
                    let nearest = nearest_other(&locus.sorted, value, eps);
                    dbscan_result(&reference, value, eps, min_points, noise, nearest)
                }
            };

            let passes = match config.method {
                Method::Percentile => match &result.evidence {
                    Evidence::Exceedance { at_least, .. } => *at_least <= config.max_exceedance,
                    _ => false,
                },
                Method::Zscore | Method::Robustz => result.significance <= p_gate,
                Method::Dbscan => matches!(result.evidence, Evidence::Cluster { noise: true, .. }),
            } && (config.min_units <= 0.0 || result.score > config.min_units)
                // Nothing to compare against: every allele at this locus belongs to the
                // sample under test, or to a relative of it.
                && result.evidence.ref_n() > 0;

            if !passes {
                continue;
            }
            passing += 1;
            let better = match &best {
                None => true,
                Some((_, _, prev)) => result.score > prev.score,
            };
            if better {
                best = Some((label, value, result));
            }
        }

        if let Some((label, value, result)) = best {
            rows.push(Row {
                chrom: locus.chrom.clone(),
                begin: locus.begin,
                end: locus.end,
                motif: locus.motif.clone(),
                sample: base_name(name).to_string(),
                event: "outlier",
                allele: label,
                length: value,
                alleles_passing: passing,
                result: Some(result),
            });
        }
    }
    rows
}

/// Whether any allele at this locus reaches the minimum expansion size
fn passes_minsize(line: &str, minsize: u32) -> bool {
    if minsize == 0 {
        return true;
    }
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() <= crate::filetype::STR_FIXED_COLUMNS {
        return false;
    }
    f.iter().skip(crate::filetype::STR_FIXED_COLUMNS).any(|s| {
        s.parse::<f64>()
            .map(|v| v >= minsize as f64)
            .unwrap_or(false)
    })
}

/// Run the outlier scan
pub fn scan(config: ScanConfig) {
    if config.threads > 0
        && rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build_global()
            .is_err()
    {
        eprintln!("Warning: thread pool already configured, using existing settings");
    }

    let path = config.combined.display().to_string();
    let file = crate::utils::reader(&path);
    let mut lines = file.lines();
    let header = crate::utils::skip_metadata_lines(&mut lines, &path);
    let columns: Vec<&str> = header.split('\t').collect();
    if columns.len() <= crate::filetype::STR_FIXED_COLUMNS {
        eprintln!("ERROR: {} does not look like a combined call file", path);
        std::process::exit(1);
    }
    // One name per sample, taken from the H1 column of each pair
    let sample_names: Vec<String> = columns[crate::filetype::STR_FIXED_COLUMNS..]
        .chunks(2)
        .map(|c| base_name(c[0]).to_string())
        .collect();
    let n_samples = sample_names.len();
    eprintln!("Scanning {} samples with method '{}'", n_samples, config.method.as_str());

    let keep: Option<Vec<bool>> = config.samples.as_ref().map(|wanted| {
        let set: std::collections::HashSet<&str> = wanted.iter().map(|s| s.as_str()).collect();
        let flags: Vec<bool> = sample_names
            .iter()
            .map(|n| set.contains(n.as_str()))
            .collect();
        let found = flags.iter().filter(|f| **f).count();
        if found != wanted.len() {
            eprintln!(
                "Warning: {} of {} requested samples were not found in the file",
                wanted.len() - found,
                wanted.len()
            );
        }
        flags
    });
    let groups = config.related.as_deref().map(read_related);
    let regions = config.regions.as_deref().map(Regions::read);
    if config.report_dropout && regions.is_none() {
        eprintln!(
            "Warning: --report-dropout without --regions emits one row per missing genotype \
             across the whole catalog. Pass --regions to restrict it to loci of interest."
        );
    }

    let mut out: Box<dyn Write> = match &config.output {
        Some(p) => {
            let f = std::fs::File::create(p).unwrap_or_else(|e| {
                eprintln!("ERROR: Failed to create {}: {}", p.display(), e);
                std::process::exit(1);
            });
            if p.extension().is_some_and(|e| e == "gz") {
                Box::new(std::io::BufWriter::with_capacity(
                    1 << 20,
                    flate2::write::GzEncoder::new(f, flate2::Compression::default()),
                ))
            } else {
                Box::new(std::io::BufWriter::with_capacity(1 << 20, f))
            }
        }
        None => Box::new(std::io::BufWriter::with_capacity(1 << 20, std::io::stdout())),
    };

    writeln!(out, "# file_type=outlier_calls").unwrap();
    writeln!(out, "# version={}", crate::VERSION).unwrap();
    writeln!(out, "# command=outlier").unwrap();
    writeln!(out, "# source={}", path).unwrap();
    writeln!(out, "# method={}", config.method.as_str()).unwrap();
    writeln!(out, "# reference=internal_leave_one_out").unwrap();
    writeln!(out, "# leave_one_out={}", config.method.uses_leave_one_out()).unwrap();
    writeln!(out, "# zygosity={}", format!("{:?}", config.zygosity).to_lowercase()).unwrap();
    writeln!(out, "# length_axis=absolute_bp").unwrap();
    writeln!(out, "# score_units=repeat_units").unwrap();
    writeln!(out, "# max_exceedance={}", config.max_exceedance).unwrap();
    writeln!(out, "# min_units={}", config.min_units).unwrap();
    writeln!(out, "# samples={}", n_samples).unwrap();
    writeln!(out, "# report_dropout={}", config.report_dropout).unwrap();
    writeln!(
        out,
        "chromosome\tbegin\tend\tmotif\tsample\tevent\tallele\tlength\talleles_passing\tscore\tsignificance\tmethod\tevidence"
    )
    .unwrap();

    const CHUNK: usize = 2048;
    let mut buffer: Vec<String> = Vec::with_capacity(CHUNK);
    let mut loci = 0usize;
    let mut reported = 0usize;
    let mut emitted = 0usize;
    let method_name = config.method.as_str();

    let flush = |buffer: &mut Vec<String>, out: &mut Box<dyn Write>, emitted: &mut usize| {
        let rows: Vec<Vec<Row>> = buffer
            .par_iter()
            .filter(|l| match &regions {
                Some(r) => {
                    let mut f = l.split('\t');
                    match (f.next(), f.next(), f.next()) {
                        (Some(c), Some(b), Some(e)) => match (b.parse(), e.parse()) {
                            (Ok(b), Ok(e)) => r.overlaps(c, b, e),
                            _ => false,
                        },
                        _ => false,
                    }
                }
                None => true,
            })
            .filter(|l| passes_minsize(l, config.minsize))
            .filter_map(|l| parse_locus(l, n_samples))
            .map(|locus| scan_locus(&locus, &config, &sample_names, &keep, &groups))
            .collect();
        for group in rows {
            for r in group {
                let (score, sig, ev) = match &r.result {
                    Some(res) => (
                        format!("{:.6}", res.score),
                        format!("{:.6e}", res.significance),
                        res.evidence.render(),
                    ),
                    None => ("NaN".to_string(), "NaN".to_string(), String::from(".")),
                };
                let length = if r.length.is_nan() {
                    "NaN".to_string()
                } else {
                    format!("{}", r.length)
                };
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    r.chrom,
                    r.begin,
                    r.end,
                    r.motif,
                    r.sample,
                    r.event,
                    r.allele,
                    length,
                    r.alleles_passing,
                    score,
                    sig,
                    method_name,
                    ev
                )
                .unwrap();
                *emitted += 1;
            }
        }
        buffer.clear();
    };

    for line in lines {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("ERROR: Failed reading {}: {}", path, e);
                std::process::exit(1);
            }
        };
        buffer.push(line);
        loci += 1;
        if buffer.len() == CHUNK {
            flush(&mut buffer, &mut out, &mut emitted);
            if loci / 1_000_000 != reported {
                reported = loci / 1_000_000;
                eprintln!("  {} loci scanned, {} calls emitted...", loci, emitted);
            }
        }
    }
    if !buffer.is_empty() {
        flush(&mut buffer, &mut out, &mut emitted);
    }
    out.flush().unwrap();
    eprintln!("Scanned {} loci, emitted {} calls", loci, emitted);
    if n_samples > 0 {
        eprintln!("Mean of {} calls per sample", emitted / n_samples.max(1));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg(method: Method) -> ScanConfig {
        ScanConfig {
            combined: PathBuf::from("x"),
            output: None,
            method,
            threads: 1,
            zygosity: Zygosity::Both,
            max_exceedance: 0,
            zscore: 3.0,
            min_units: 0.0,
            minsize: 0,
            samples: None,
            related: None,
            report_dropout: false,
            min_call_rate: 0.0,
            regions: None,
        }
    }

    fn names(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("s{}", i)).collect()
    }

    #[test]
    fn test_parse_locus_converts_and_keeps_missing() {
        // reference length 30; deltas +9, -3, NaN, 0
        let line = "chr1\t100\t130\tCAG\t9\t-3\tNaN\t0";
        let l = parse_locus(line, 2).unwrap();
        assert_eq!(l.motif_len, 3.0);
        assert_eq!(l.called, 3);
        assert_eq!(l.total, 4);
        assert_eq!(l.sorted, vec![27.0, 30.0, 39.0]);
        assert_eq!(l.per_sample[1].h1, None, "NaN must stay missing");
    }

    #[test]
    fn test_unique_maximum_is_reported_once_per_sample() {
        // Four samples; the last carries a large expansion on one allele.
        let line = "chr1\t0\t30\tCAG\t0\t0\t0\t3\t0\t0\t60\t0";
        let l = parse_locus(line, 4).unwrap();
        let rows = scan_locus(&l, &cfg(Method::Percentile), &names(4), &None, &None);
        assert_eq!(rows.len(), 1, "only the unique maximum should pass the default gate");
        assert_eq!(rows[0].sample, "s3");
        assert_eq!(rows[0].allele, "H1");
        assert_eq!(rows[0].length, 90.0);
    }

    #[test]
    fn test_monomorphic_locus_emits_nothing() {
        // Every sample identical: nobody is longer than the reference maximum.
        let line = "chr1\t0\t30\tCAG\t0\t0\t0\t0\t0\t0";
        let l = parse_locus(line, 3).unwrap();
        let rows = scan_locus(&l, &cfg(Method::Percentile), &names(3), &None, &None);
        assert!(rows.is_empty(), "a tie at the top is not an outlier");
    }

    #[test]
    fn test_min_zygosity_finds_biallelic_signal_that_max_misses() {
        // s2 carries one very long allele; s3 carries two moderately long ones. Under `max`
        // the single long allele dominates and s3 is invisible; under `min` the biallelic
        // carrier is the one that stands out. Each mode compares against a reference built
        // from the same allele slot, which is what makes both views possible.
        let line = "chr1\t0\t30\tCAG\t0\t0\t0\t0\t0\t120\t60\t60";
        let l = parse_locus(line, 4).unwrap();
        let nm = names(4);

        let mut c = cfg(Method::Percentile);
        c.zygosity = Zygosity::Max;
        let by_max = scan_locus(&l, &c, &nm, &None, &None);
        assert_eq!(by_max.len(), 1);
        assert_eq!(by_max[0].sample, "s2", "max finds the single very long allele");

        c.zygosity = Zygosity::Min;
        let by_min = scan_locus(&l, &c, &nm, &None, &None);
        assert_eq!(by_min.len(), 1);
        assert_eq!(by_min[0].sample, "s3", "min finds the biallelic carrier that max missed");
    }

    #[test]
    fn test_related_samples_are_excluded_from_each_others_reference() {
        // Two siblings share an expansion. Tested against the whole cohort each masks the
        // other; with the relatedness group applied, both surface.
        let line = "chr1\t0\t30\tCAG\t0\t0\t0\t0\t60\t0\t60\t0";
        let l = parse_locus(line, 4).unwrap();
        let nm = names(4);

        let plain = scan_locus(&l, &cfg(Method::Percentile), &nm, &None, &None);
        assert!(plain.is_empty(), "the sibling pair mutually suppresses without exclusion");

        let mut groups = HashMap::new();
        groups.insert("s2".to_string(), 0usize);
        groups.insert("s3".to_string(), 0usize);
        let with = scan_locus(&l, &cfg(Method::Percentile), &nm, &None, &Some(groups));
        assert_eq!(with.len(), 2, "excluding relatives must reveal both carriers");
    }

    #[test]
    fn test_dropout_is_reported_as_its_own_event() {
        let line = "chr1\t0\t30\tCAG\t0\t0\tNaN\tNaN\t60\t0";
        let l = parse_locus(line, 3).unwrap();
        let mut c = cfg(Method::Percentile);
        c.report_dropout = true;
        let rows = scan_locus(&l, &c, &names(3), &None, &None);
        let dropouts: Vec<_> = rows.iter().filter(|r| r.event == "dropout").collect();
        assert_eq!(dropouts.len(), 1);
        assert_eq!(dropouts[0].sample, "s1");
        assert!(dropouts[0].length.is_nan(), "a dropout has no length, and is not zero");
    }

    #[test]
    fn test_call_rate_filter_skips_poorly_called_loci() {
        // Four samples, two of them entirely uncalled: a 50% call rate.
        let line = "chr1\t0\t30\tCAG\tNaN\tNaN\tNaN\tNaN\t0\t0\t60\t0";
        let l = parse_locus(line, 4).unwrap();
        let mut c = cfg(Method::Percentile);
        c.min_call_rate = 0.6;
        assert!(scan_locus(&l, &c, &names(4), &None, &None).is_empty());
        c.min_call_rate = 0.4;
        assert_eq!(scan_locus(&l, &c, &names(4), &None, &None).len(), 1);
    }

    #[test]
    fn test_locus_with_no_reference_left_yields_nothing() {
        // Only one sample is called, so removing it leaves no reference to compare against.
        // A call here would be an artefact of comparing a sample with itself.
        let line = "chr1\t0\t30\tCAG\tNaN\tNaN\tNaN\tNaN\t60\t0";
        let l = parse_locus(line, 3).unwrap();
        let mut c = cfg(Method::Percentile);
        c.min_call_rate = 0.0;
        assert!(scan_locus(&l, &c, &names(3), &None, &None).is_empty());
    }

    #[test]
    fn test_regions_restrict_the_scan() {
        use std::io::Write;
        let mut bed = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        writeln!(bed, "chr1\t50\t80").unwrap();
        bed.flush().unwrap();
        let r = Regions::read(bed.path());
        assert!(r.overlaps("chr1", 60, 70), "contained interval overlaps");
        assert!(r.overlaps("1", 70, 90), "chr prefix must not matter, partial overlap counts");
        assert!(!r.overlaps("chr1", 80, 90), "half-open: touching at the end is not overlap");
        assert!(!r.overlaps("chr1", 10, 50), "touching at the start is not overlap");
        assert!(!r.overlaps("chr2", 60, 70), "different chromosome");
    }

    #[test]
    fn test_min_units_gate_is_in_repeat_units() {
        // 6 bp longer than the reference maximum is 2 units on a trinucleotide
        let line = "chr1\t0\t30\tCAG\t0\t0\t0\t0\t6\t0";
        let l = parse_locus(line, 3).unwrap();
        let mut c = cfg(Method::Percentile);
        c.min_units = 1.5;
        assert_eq!(scan_locus(&l, &c, &names(3), &None, &None).len(), 1);
        c.min_units = 2.5;
        assert!(scan_locus(&l, &c, &names(3), &None, &None).is_empty());
    }
}
