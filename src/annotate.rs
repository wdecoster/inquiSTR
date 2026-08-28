//! Join a published population table to a local catalog.
//!
//! This is deliberately a separate, inspectable step rather than a lookup hidden inside the
//! outlier or prioritisation commands: the mapping it produces is auditable, and the match rate
//! is reported before any ranking depends on it.
//!
//! The published index is derived from the **longest pure segment**, a quantity inquiSTR does
//! not measure. An imported value is therefore named for its own axis and never substituted for
//! a locally computed length-based index. Where both exist, both are reported.

use crate::locus_stats::{LocusKey, parse_row};
use clap::ValueEnum;
use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

/// Which quantity the imported table describes
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum Axis {
    /// Longest pure segment length: the published PLV/PLVI axis
    Lps,
    /// Total allele length: the same axis inquiSTR measures
    Length,
}

impl Axis {
    /// Column name for the index derived from this axis. The name states both what was
    /// measured and that it came from elsewhere, so it can never be read as a local value.
    fn column(&self) -> &'static str {
        match self {
            Axis::Lps => "plvi_lps_external",
            Axis::Length => "lvi_length_external",
        }
    }
}

#[derive(Debug, Clone)]
pub struct AnnotateConfig {
    /// Published locus statistics table
    pub table: PathBuf,
    /// BED catalog the local calls were made against
    pub catalog: PathBuf,
    pub output: PathBuf,
    pub axis: Option<Axis>,
    /// Genome build of the local catalog
    pub build: String,
    /// Minimum reciprocal overlap for a non-exact match
    pub min_overlap: f64,
    /// Refuse if fewer than this fraction of catalog loci map
    pub min_match_rate: f64,
    /// Minimum loci per motif-length bin when ranking the imported values
    pub min_bin_size: usize,
}

/// How confidently a catalog locus was matched
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Quality {
    Exact,
    HighConfidence,
    Ambiguous,
    Unmapped,
}

impl Quality {
    fn as_str(&self) -> &'static str {
        match self {
            Quality::Exact => "exact",
            Quality::HighConfidence => "high_confidence",
            Quality::Ambiguous => "ambiguous",
            Quality::Unmapped => "unmapped",
        }
    }
}

/// Reverse complement of a motif
fn revcomp(motif: &str) -> String {
    motif
        .chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// Canonical form of a motif: the lexicographically smallest rotation of the motif or its
/// reverse complement. Two motifs describing the same repeat share a canonical form, which is
/// what makes `CAG`, `AGC` and `CTG` compatible.
pub fn canonical_motif(motif: &str) -> String {
    let m = motif.to_ascii_uppercase();
    // Rotation slices by byte offset, so anything that is not a plain ACGT motif has to be
    // rejected here rather than panicking on a character boundary. Imported tables are
    // untrusted input: the motif is whatever followed the third `-` in an identifier.
    if m.is_empty() || m == "." || !m.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
        return String::from(".");
    }
    let mut best: Option<String> = None;
    for candidate in [m.clone(), revcomp(&m)] {
        let n = candidate.len();
        for i in 0..n {
            let rotated: String = candidate[i..]
                .chars()
                .chain(candidate[..i].chars())
                .collect();
            if best.as_ref().is_none_or(|b| rotated < *b) {
                best = Some(rotated);
            }
        }
    }
    best.unwrap_or_else(|| String::from("."))
}

/// Whether two motifs describe the same repeat
fn motifs_compatible(a: &str, b: &str) -> bool {
    if a == "." || b == "." {
        return true;
    }
    a.len() == b.len() && canonical_motif(a) == canonical_motif(b)
}

/// One catalog interval
struct CatalogLocus {
    key: LocusKey,
    motif: String,
}

/// Read a BED catalog, taking column 4 as the motif where it is one
fn read_catalog(path: &Path) -> Vec<CatalogLocus> {
    let reader = crate::utils::reader(&path.to_string_lossy());
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("ERROR: Failed reading {}: {}", path.display(), e);
                std::process::exit(1);
            }
        };
        if line.starts_with('#') || line.starts_with("track") || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let (begin, end) = match (f[1].parse::<u32>(), f[2].parse::<u32>()) {
            (Ok(b), Ok(e)) => (b, e),
            _ => continue,
        };
        let motif = f
            .get(3)
            .and_then(|s| crate::varindex::parse_motif(s))
            .unwrap_or(".")
            .to_string();
        out.push(CatalogLocus { key: LocusKey::new(f[0], begin, end), motif });
    }
    out
}

/// Reciprocal overlap of two intervals
fn reciprocal_overlap(a: &LocusKey, b: &LocusKey) -> f64 {
    let lo = a.begin.max(b.begin);
    let hi = a.end.min(b.end);
    if hi <= lo {
        return 0.0;
    }
    let ov = (hi - lo) as f64;
    let la = (a.end.saturating_sub(a.begin)).max(1) as f64;
    let lb = (b.end.saturating_sub(b.begin)).max(1) as f64;
    (ov / la).min(ov / lb)
}

/// Rank imported values within motif-length bins, reproducing the published index definition.
///
/// The published tables ship the standard deviation, not the index; the index is that value
/// ranked within groups of comparable motif length, so it has to be recomputed here.
fn rank_within_bins(
    entries: &[(LocusKey, String, f64)],
    min_bin_size: usize,
) -> HashMap<LocusKey, f64> {
    let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
    for (_, motif, std) in entries {
        if motif != "." && std.is_finite() {
            *counts.entry(motif.len()).or_insert(0) += 1;
        }
    }
    let bins = crate::varindex::bin_motif_lengths(&counts, min_bin_size);
    let bin_of = |len: usize| bins.iter().position(|&(lo, hi)| len >= lo && len <= hi);

    let mut per_bin: Vec<Vec<f64>> = vec![Vec::new(); bins.len()];
    for (_, motif, std) in entries {
        if motif != "."
            && std.is_finite()
            && let Some(b) = bin_of(motif.len())
        {
            per_bin[b].push(*std);
        }
    }
    for v in per_bin.iter_mut() {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    }

    let mut out = HashMap::with_capacity(entries.len());
    for (key, motif, std) in entries {
        if motif == "." || !std.is_finite() {
            continue;
        }
        if let Some(b) = bin_of(motif.len()) {
            let sorted = &per_bin[b];
            if sorted.is_empty() {
                continue;
            }
            let less = sorted.partition_point(|x| *x < *std);
            out.insert(key.clone(), less as f64 / sorted.len() as f64);
        }
    }
    eprintln!(
        "Ranked imported values within {} motif-length bins (floor {} loci)",
        bins.len(),
        min_bin_size
    );
    out
}

/// Build the mapping from a local catalog to a published table
pub fn annotate(config: AnnotateConfig) {
    if !config.build.eq_ignore_ascii_case("GRCh38") {
        eprintln!(
            "ERROR: the published population tables are GRCh38 / TR-Explorer v1.0.1. \
             Refusing to join them to a '{}' catalog: coordinates would not correspond.",
            config.build
        );
        std::process::exit(1);
    }

    // Read the published table, deciding which axis it describes.
    let (header, reader) = crate::locus_stats::open(&config.table);
    let mut entries: Vec<(LocusKey, String, f64)> = Vec::new();
    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("ERROR: Failed reading {}: {}", config.table.display(), e);
                std::process::exit(1);
            }
        };
        if let Some(row) = parse_row(&header, &line) {
            entries.push((row.key, row.motif, row.std));
        }
    }
    if entries.is_empty() {
        eprintln!("ERROR: {} contained no usable rows", config.table.display());
        std::process::exit(1);
    }
    eprintln!("Read {} loci from {}", entries.len(), config.table.display());

    let axis = match config.axis {
        Some(a) => a,
        None => {
            eprintln!(
                "ERROR: cannot tell whether {} describes longest-pure-segment or total allele \
                 length. Pass --axis lps or --axis length. These are different measurements and \
                 must not be conflated.",
                config.table.display()
            );
            std::process::exit(1);
        }
    };
    eprintln!("Treating imported values as the {:?} axis -> column {}", axis, axis.column());

    let index = rank_within_bins(&entries, config.min_bin_size);

    // Index the published loci by chromosome for interval matching
    // Per chromosome: loci sorted by start, plus the longest interval on it. The longest
    // interval bounds how far back a candidate overlap can begin, which turns the search for
    // overlapping loci into a bounded window instead of a scan from the start of the
    // chromosome. Without that bound this is quadratic and does not finish on a real catalog.
    let mut by_chrom: HashMap<String, (Vec<usize>, u32)> = HashMap::new();
    for (i, (key, _, _)) in entries.iter().enumerate() {
        let slot = by_chrom
            .entry(key.chrom.clone())
            .or_insert_with(|| (Vec::new(), 0));
        slot.0.push(i);
        slot.1 = slot.1.max(key.end.saturating_sub(key.begin));
    }
    for (v, _) in by_chrom.values_mut() {
        v.sort_by_key(|i| entries[*i].0.begin);
    }
    let exact: HashMap<&LocusKey, usize> = entries
        .iter()
        .enumerate()
        .map(|(i, (k, _, _))| (k, i))
        .collect();

    let catalog = read_catalog(&config.catalog);
    eprintln!("Read {} loci from {}", catalog.len(), config.catalog.display());
    if catalog.is_empty() {
        eprintln!("ERROR: catalog {} is empty", config.catalog.display());
        std::process::exit(1);
    }

    let out_file = std::fs::File::create(&config.output).unwrap_or_else(|e| {
        eprintln!("ERROR: Failed to create {}: {}", config.output.display(), e);
        std::process::exit(1);
    });
    let sink: Box<dyn Write> = if config.output.extension().is_some_and(|e| e == "gz") {
        Box::new(flate2::write::GzEncoder::new(out_file, flate2::Compression::default()))
    } else {
        Box::new(out_file)
    };
    let mut out = std::io::BufWriter::with_capacity(1 << 20, sink);

    let mut counts: HashMap<&'static str, usize> = HashMap::new();
    let mut rows: Vec<String> = Vec::with_capacity(catalog.len());

    for locus in &catalog {
        let mut quality = Quality::Unmapped;
        let mut value: Option<f64> = None;
        let mut candidates = 0usize;
        let mut matched_std = f64::NAN;

        if let Some(&i) = exact.get(&locus.key)
            && motifs_compatible(&locus.motif, &entries[i].1)
        {
            quality = Quality::Exact;
            candidates = 1;
            value = index.get(&entries[i].0).copied();
            matched_std = entries[i].2;
        }

        if quality == Quality::Unmapped
            && let Some((idxs, max_len)) = by_chrom.get(&locus.key.chrom)
        {
            // Nothing starting before this point can still reach the query interval
            let window_start = locus.key.begin.saturating_sub(*max_len);
            let first = idxs.partition_point(|&i| entries[i].0.begin < window_start);
            let mut hits: Vec<usize> = Vec::new();
            for &i in &idxs[first..] {
                let other = &entries[i].0;
                if other.begin >= locus.key.end {
                    break;
                }
                if reciprocal_overlap(&locus.key, other) >= config.min_overlap
                    && motifs_compatible(&locus.motif, &entries[i].1)
                {
                    hits.push(i);
                }
            }
            candidates = hits.len();
            match hits.len() {
                0 => {}
                1 => {
                    quality = Quality::HighConfidence;
                    value = index.get(&entries[hits[0]].0).copied();
                    matched_std = entries[hits[0]].2;
                }
                _ => {
                    // Several published loci overlap this one - the compound-locus case.
                    // Resolving it by taking the largest value would systematically inflate
                    // exactly the messiest loci, so no value is assigned at all.
                    quality = Quality::Ambiguous;
                }
            }
        }

        *counts.entry(quality.as_str()).or_insert(0) += 1;
        rows.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            locus.key.chrom,
            locus.key.begin,
            locus.key.end,
            locus.motif,
            quality.as_str(),
            candidates,
            if matched_std.is_finite() {
                format!("{}", matched_std)
            } else {
                ".".to_string()
            },
            match value {
                Some(v) => format!("{:.6}", v),
                None => ".".to_string(),
            }
        ));
    }

    let mapped = counts.get("exact").copied().unwrap_or(0)
        + counts.get("high_confidence").copied().unwrap_or(0);
    let match_rate = mapped as f64 / catalog.len() as f64;

    eprintln!("\nJoin quality:");
    for label in ["exact", "high_confidence", "ambiguous", "unmapped"] {
        let n = counts.get(label).copied().unwrap_or(0);
        eprintln!("  {:<16} {:>9} ({:5.2}%)", label, n, 100.0 * n as f64 / catalog.len() as f64);
    }
    eprintln!("  match rate: {:.2}%", 100.0 * match_rate);

    if match_rate < config.min_match_rate {
        eprintln!(
            "\nERROR: only {:.2}% of catalog loci mapped, below the required {:.2}%. \
             This usually means the catalog is a different version or a different genome build \
             from the published table. Refusing to write a mapping that would mislead.",
            100.0 * match_rate,
            100.0 * config.min_match_rate
        );
        std::process::exit(1);
    }

    writeln!(out, "# file_type=locus_stats").unwrap();
    writeln!(out, "# version={}", crate::VERSION).unwrap();
    writeln!(out, "# command=annotate").unwrap();
    writeln!(out, "# table={}", config.table.display()).unwrap();
    writeln!(out, "# catalog={}", config.catalog.display()).unwrap();
    writeln!(out, "# build={}", config.build).unwrap();
    writeln!(out, "# axis={:?}", axis).unwrap();
    writeln!(out, "# index_column={}", axis.column()).unwrap();
    writeln!(out, "# match_rate={:.6}", match_rate).unwrap();
    writeln!(out, "# min_overlap={}", config.min_overlap).unwrap();
    writeln!(
        out,
        "chromosome\tbegin\tend\tmotif\tjoin_quality\tcandidates\tstd\t{}",
        axis.column()
    )
    .unwrap();
    for row in rows {
        writeln!(out, "{}", row).unwrap();
    }
    out.flush().unwrap();
    eprintln!("Wrote {}", config.output.display());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonical_motif_groups_equivalent_repeats() {
        // Rotations and reverse complements of the same repeat share a canonical form
        assert_eq!(canonical_motif("CAG"), canonical_motif("AGC"));
        assert_eq!(canonical_motif("CAG"), canonical_motif("GCA"));
        assert_eq!(canonical_motif("CAG"), canonical_motif("CTG"), "revcomp of CAG is CTG");
        assert_ne!(canonical_motif("CAG"), canonical_motif("CGG"));
    }

    #[test]
    fn test_non_acgt_motifs_are_rejected_not_sliced() {
        // A motif from an imported table is untrusted; multi-byte characters must not reach
        // the rotation loop, which slices by byte offset.
        assert_eq!(canonical_motif("Aé"), ".");
        assert_eq!(canonical_motif("ID=locus1"), ".");
        assert_eq!(canonical_motif(""), ".");
        // An unresolvable motif must not join to anything. An absent motif (".") is a
        // different case: it carries no claim and so cannot contradict one.
        assert!(!motifs_compatible("Aé", "CAG"), "garbage must not silently match a real motif");
        assert!(motifs_compatible(".", "CAG"), "an absent motif still joins");
    }

    #[test]
    fn test_motif_compatibility_requires_same_length() {
        assert!(motifs_compatible("CAG", "CTG"));
        assert!(!motifs_compatible("CAG", "CAGCAG"), "a doubled motif is not the same locus");
        // An absent motif cannot contradict anything
        assert!(motifs_compatible(".", "CAG"));
    }

    #[test]
    fn test_reciprocal_overlap() {
        let a = LocusKey::new("chr1", 100, 200);
        assert_eq!(reciprocal_overlap(&a, &LocusKey::new("chr1", 100, 200)), 1.0);
        assert_eq!(reciprocal_overlap(&a, &LocusKey::new("chr1", 300, 400)), 0.0);
        // Half of the larger interval, all of the smaller: reciprocal overlap is the smaller
        let partial = reciprocal_overlap(&a, &LocusKey::new("chr1", 150, 250));
        assert!((partial - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_rank_within_bins_reproduces_a_percentile_rank() {
        let entries: Vec<(LocusKey, String, f64)> = (0..100)
            .map(|i| (LocusKey::new("chr1", i * 10, i * 10 + 9), "CAG".to_string(), i as f64))
            .collect();
        let ranked = rank_within_bins(&entries, 10);
        assert_eq!(ranked.len(), 100);
        // The smallest value sits at rank 0, the largest just below 1
        assert_eq!(ranked[&LocusKey::new("chr1", 0, 9)], 0.0);
        assert!((ranked[&LocusKey::new("chr1", 990, 999)] - 0.99).abs() < 1e-9);
    }

    #[test]
    fn test_axis_columns_are_distinct_and_named_for_their_source() {
        assert_eq!(Axis::Lps.column(), "plvi_lps_external");
        assert_eq!(Axis::Length.column(), "lvi_length_external");
        assert_ne!(Axis::Lps.column(), Axis::Length.column());
    }
}
