//! Per-locus summary statistics and the length-variability index.
//!
//! `varindex` collapses a combined call file into one row per locus: the allele-length
//! distribution of the cohort, plus a variability index obtained by ranking loci against
//! others with a comparable motif length.
//!
//! Allele lengths are reported as **absolute lengths in bp**, converted from inquiSTR's
//! reference-relative genotypes as `(end - begin) + value`. This matches the layout of the
//! published population tables so that a locally computed table and an imported one can be
//! read by the same code.
//!
//! See `docs/PLVI_PORT.md` for the design rationale and `docs/PRIORITIZATION.md` for usage.

use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

/// Percentile grid, matching the published locus tables.
/// Stored as (percentile, column label suffix).
pub const PERCENTILES: [(f64, &str); 24] = [
    (0.0, "0th"),
    (1.0, "1st"),
    (5.0, "5th"),
    (10.0, "10th"),
    (15.0, "15th"),
    (20.0, "20th"),
    (25.0, "25th"),
    (30.0, "30th"),
    (35.0, "35th"),
    (40.0, "40th"),
    (45.0, "45th"),
    (50.0, "50th"),
    (55.0, "55th"),
    (60.0, "60th"),
    (65.0, "65th"),
    (70.0, "70th"),
    (75.0, "75th"),
    (80.0, "80th"),
    (85.0, "85th"),
    (90.0, "90th"),
    (95.0, "95th"),
    (99.0, "99th"),
    (99.9, "99.9th"),
    (100.0, "100th"),
];

/// Configuration for the varindex command
#[derive(Debug, Clone)]
pub struct VarIndexConfig {
    pub combined: PathBuf,
    pub output: PathBuf,
    pub threads: usize,
    /// Minimum number of loci a motif-length bin must contain before it can be ranked
    pub min_bin_size: usize,
    /// Minimum number of called alleles a locus needs before it is eligible for ranking
    pub min_alleles: usize,
    pub tmpdir: Option<PathBuf>,
}

/// Per-locus allele-length summary
#[derive(Debug, Clone)]
pub struct LocusStats {
    pub chrom: String,
    pub begin: u32,
    pub end: u32,
    pub motif: String,
    /// Number of called (non-NaN) alleles
    pub n: usize,
    pub percentiles: [f64; PERCENTILES.len()],
    pub mean: f64,
    pub std: f64,
    pub mad: f64,
    pub mode: f64,
    pub unique: usize,
    /// Number of alleles whose absolute length was negative and clamped to zero
    pub clamped: usize,
}

/// Nearest-rank percentile on an ascending-sorted slice.
///
/// Deliberately not interpolated: allele lengths are heavily tied integers, and an
/// interpolated quantile can land on a value no allele actually has. Every value this
/// returns is an observed length.
fn nearest_rank(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    if p <= 0.0 {
        return sorted[0];
    }
    let n = sorted.len() as f64;
    // ceil(p/100 * n) as a 1-based rank, clamped into range
    let rank = (p / 100.0 * n).ceil().max(1.0) as usize;
    sorted[rank.min(sorted.len()) - 1]
}

/// Median of an ascending-sorted slice
fn median_sorted(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n.is_multiple_of(2) {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

/// Most frequent value in an ascending-sorted slice; ties resolve to the smallest value
fn mode_sorted(sorted: &[f64]) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let (mut best_val, mut best_count) = (sorted[0], 1usize);
    let (mut cur_val, mut cur_count) = (sorted[0], 1usize);
    for &v in &sorted[1..] {
        if v == cur_val {
            cur_count += 1;
        } else {
            if cur_count > best_count {
                best_count = cur_count;
                best_val = cur_val;
            }
            cur_val = v;
            cur_count = 1;
        }
    }
    if cur_count > best_count {
        best_val = cur_val;
    }
    best_val
}

/// Count distinct values in an ascending-sorted slice
fn count_unique(sorted: &[f64]) -> usize {
    if sorted.is_empty() {
        return 0;
    }
    let mut n = 1;
    for w in sorted.windows(2) {
        if w[0] != w[1] {
            n += 1;
        }
    }
    n
}

/// Extract a usable motif from the `info` column.
///
/// Returns `None` when the field carries something other than a motif, which is the case for
/// catalogs that use column 4 for an identifier. Comma-separated alternatives (as used by the
/// STRchive catalog) resolve to the first entry; only motif *length* is used downstream, and
/// listed alternatives are near-always the same length.
pub fn parse_motif(info: &str) -> Option<&str> {
    let candidate = info.split(',').next()?.trim();
    if candidate.is_empty() || candidate == "." {
        return None;
    }
    if candidate
        .bytes()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
    {
        Some(candidate)
    } else {
        None
    }
}

/// Compute the summary for one locus from its absolute allele lengths.
/// `values` is sorted in place.
pub fn summarise(
    chrom: &str,
    begin: u32,
    end: u32,
    motif: &str,
    values: &mut [f64],
    clamped: usize,
) -> LocusStats {
    values.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = values.len();

    let mut percentiles = [f64::NAN; PERCENTILES.len()];
    let (mut mean, mut std, mut mad, mut mode) = (f64::NAN, f64::NAN, f64::NAN, f64::NAN);

    if n > 0 {
        for (i, (p, _)) in PERCENTILES.iter().enumerate() {
            percentiles[i] = nearest_rank(values, *p);
        }
        let sum: f64 = values.iter().sum();
        mean = sum / n as f64;
        // Population standard deviation, computed from the centred sum of squares so that
        // heavily tied integer data does not lose precision the way sum-of-squares does.
        let var = values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n as f64;
        std = var.max(0.0).sqrt();

        let med = median_sorted(values);
        let mut devs: Vec<f64> = values.iter().map(|v| (v - med).abs()).collect();
        devs.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        mad = median_sorted(&devs);

        mode = mode_sorted(values);
    }

    LocusStats {
        chrom: chrom.to_string(),
        begin,
        end,
        motif: motif.to_string(),
        n,
        percentiles,
        mean,
        std,
        mad,
        mode,
        unique: count_unique(values),
        clamped,
    }
}

/// Parse one combined-file data line into a locus summary.
///
/// Genotypes are converted from reference-relative to absolute length. `NaN` and empty fields
/// are treated as missing and excluded, never coerced to zero. Absolute lengths below zero
/// (a deletion longer than the reference locus) are clamped to zero and counted.
fn summarise_line(line: &str) -> Option<LocusStats> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() <= crate::filetype::STR_FIXED_COLUMNS {
        return None;
    }
    let chrom = fields[0];
    let begin: u32 = fields[1].parse().ok()?;
    let end: u32 = fields[2].parse().ok()?;
    let motif = parse_motif(fields[3]).unwrap_or(".");
    let reflen = end.saturating_sub(begin) as f64;

    let mut values = Vec::with_capacity(fields.len() - crate::filetype::STR_FIXED_COLUMNS);
    let mut clamped = 0usize;
    for field in fields.iter().skip(crate::filetype::STR_FIXED_COLUMNS) {
        if field.is_empty() || field.eq_ignore_ascii_case("nan") || *field == "." {
            continue;
        }
        match field.parse::<f64>() {
            Ok(v) if v.is_finite() => {
                let abs = reflen + v;
                if abs < 0.0 {
                    clamped += 1;
                    values.push(0.0);
                } else {
                    values.push(abs);
                }
            }
            _ => continue,
        }
    }

    Some(summarise(chrom, begin, end, motif, &mut values, clamped))
}

/// Group adjacent motif lengths into bins holding at least `floor` loci each.
///
/// Grouping by exact motif length is a proxy for comparing like with like, and it breaks down
/// in the tail: most catalogs hold a few loci at each of many long motif lengths, and a locus
/// alone in its group ranks first of one. Merging adjacent lengths keeps those loci rankable
/// instead of discarding them. Returns inclusive `(low, high)` motif-length ranges.
pub fn bin_motif_lengths(counts: &BTreeMap<usize, usize>, floor: usize) -> Vec<(usize, usize)> {
    let mut bins: Vec<(usize, usize)> = Vec::new();
    let (mut low, mut acc) = (None::<usize>, 0usize);

    for (&len, &count) in counts {
        if low.is_none() {
            low = Some(len);
        }
        acc += count;
        if acc >= floor {
            bins.push((low.unwrap(), len));
            low = None;
            acc = 0;
        }
    }
    // A trailing run that never reached the floor folds into the previous bin, so that no
    // locus is left in a bin too small to rank.
    if let Some(start) = low {
        let last = *counts.keys().next_back().unwrap_or(&start);
        match bins.last_mut() {
            Some(prev) => prev.1 = last,
            None => bins.push((start, last)),
        }
    }
    bins
}

/// Map a motif length to the index of its bin
fn bin_of(bins: &[(usize, usize)], motif_len: usize) -> Option<usize> {
    bins.iter()
        .position(|&(lo, hi)| motif_len >= lo && motif_len <= hi)
}

/// Write the column header for a locus_stats file
fn write_header<W: Write>(w: &mut W) -> std::io::Result<()> {
    write!(w, "chromosome\tbegin\tend\tmotif\tN")?;
    for (_, label) in PERCENTILES.iter() {
        write!(w, "\t{}Percentile", label)?;
    }
    writeln!(w, "\tmean\tstd\tmad\tmode\tuniqueLengths\tmotif_bin\tlvi_length_internal")
}

/// Write the statistics portion of a row (everything except bin and index)
fn write_stats_row<W: Write>(w: &mut W, s: &LocusStats) -> std::io::Result<()> {
    write!(w, "{}\t{}\t{}\t{}\t{}", s.chrom, s.begin, s.end, s.motif, s.n)?;
    for v in s.percentiles.iter() {
        write!(w, "\t{}", fmt(*v))?;
    }
    write!(
        w,
        "\t{}\t{}\t{}\t{}\t{}",
        fmt(s.mean),
        fmt(s.std),
        fmt(s.mad),
        fmt(s.mode),
        s.unique
    )
}

/// Format a float compactly, rendering NaN as the string used throughout inquiSTR output
fn fmt(v: f64) -> String {
    if v.is_nan() {
        "NaN".to_string()
    } else {
        format!("{}", v)
    }
}

/// Compute per-locus statistics and the length-variability index for a combined file.
pub fn varindex(config: VarIndexConfig) {
    if config.threads > 0
        && rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build_global()
            .is_err()
    {
        eprintln!("Warning: thread pool already configured, using existing settings");
    }

    let combined_path = config.combined.display().to_string();
    match crate::filetype::read_file_type_metadata(&config.combined) {
        Some(crate::filetype::FileType::CombinedCall) => {}
        Some(other) => {
            eprintln!(
                "ERROR: {} is a {:?} file; varindex requires a combined call file.",
                combined_path, other
            );
            std::process::exit(1);
        }
        None => eprintln!(
            "Warning: {} has no file_type metadata; assuming it is a combined call file.",
            combined_path
        ),
    }

    let file = crate::utils::reader(&combined_path);
    let mut lines = file.lines();
    let header = crate::utils::skip_metadata_lines(&mut lines, &combined_path);
    let n_samples = header
        .split('\t')
        .count()
        .saturating_sub(crate::filetype::STR_FIXED_COLUMNS)
        / 2;
    eprintln!("Reading {} samples from {}", n_samples, combined_path);

    // Pass 1: statistics per locus, streamed to a temporary file. Only the motif length and
    // the standard deviation are retained in memory, since ranking needs nothing else.
    let tmpdir = config.tmpdir.clone().unwrap_or_else(|| {
        config
            .output
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."))
    });
    let mut tmp = match tempfile::NamedTempFile::new_in(&tmpdir) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("ERROR: Failed to create temporary file in {}: {}", tmpdir.display(), e);
            std::process::exit(1);
        }
    };

    let mut motif_len_of: Vec<u32> = Vec::new();
    let mut std_of: Vec<f32> = Vec::new();
    let mut eligible: Vec<bool> = Vec::new();
    let mut total_clamped = 0usize;
    let mut no_motif = 0usize;
    let mut processed = 0usize;
    let mut reported = 0usize;

    const CHUNK: usize = 4096;
    let mut buffer: Vec<String> = Vec::with_capacity(CHUNK);

    let flush = |buffer: &mut Vec<String>,
                 tmp: &mut tempfile::NamedTempFile,
                 motif_len_of: &mut Vec<u32>,
                 std_of: &mut Vec<f32>,
                 eligible: &mut Vec<bool>,
                 total_clamped: &mut usize,
                 no_motif: &mut usize| {
        let stats: Vec<LocusStats> = buffer
            .par_iter()
            .filter_map(|l| summarise_line(l))
            .collect();
        for s in &stats {
            write_stats_row(tmp, s).expect("Failed writing temporary statistics");
            writeln!(tmp).expect("Failed writing temporary statistics");
            *total_clamped += s.clamped;
            let mlen = if s.motif == "." {
                *no_motif += 1;
                0
            } else {
                s.motif.len()
            };
            motif_len_of.push(mlen as u32);
            std_of.push(s.std as f32);
            eligible.push(mlen > 0 && s.n >= config.min_alleles && s.std.is_finite());
        }
        buffer.clear();
    };

    for line in lines {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("ERROR: Failed reading {}: {}", combined_path, e);
                std::process::exit(1);
            }
        };
        buffer.push(line);
        if buffer.len() == CHUNK {
            flush(
                &mut buffer,
                &mut tmp,
                &mut motif_len_of,
                &mut std_of,
                &mut eligible,
                &mut total_clamped,
                &mut no_motif,
            );
            processed += CHUNK;
            if processed / 1_000_000 != reported {
                reported = processed / 1_000_000;
                eprintln!("  {} loci summarised...", processed);
            }
        }
    }
    if !buffer.is_empty() {
        processed += buffer.len();
        flush(
            &mut buffer,
            &mut tmp,
            &mut motif_len_of,
            &mut std_of,
            &mut eligible,
            &mut total_clamped,
            &mut no_motif,
        );
    }
    tmp.flush().expect("Failed flushing temporary statistics");
    eprintln!("Summarised {} loci", processed);
    if no_motif > 0 {
        eprintln!(
            "Warning: {} loci have no usable motif in the info column and receive no index",
            no_motif
        );
    }
    if total_clamped > 0 {
        eprintln!(
            "Note: {} allele calls implied a negative length and were clamped to 0",
            total_clamped
        );
    }

    // Bin motif lengths, then rank loci by standard deviation within each bin.
    let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
    for (i, &mlen) in motif_len_of.iter().enumerate() {
        if eligible[i] {
            *counts.entry(mlen as usize).or_insert(0) += 1;
        }
    }
    let bins = bin_motif_lengths(&counts, config.min_bin_size);
    eprintln!(
        "Ranking within {} motif-length bins (floor {} loci): {}",
        bins.len(),
        config.min_bin_size,
        bins.iter()
            .map(|&(lo, hi)| if lo == hi {
                format!("{}", lo)
            } else {
                format!("{}-{}", lo, hi)
            })
            .collect::<Vec<_>>()
            .join(" ")
    );

    // Sorted standard deviations per bin, so each locus can be placed by how many loci in its
    // bin are strictly less variable.
    let mut per_bin: Vec<Vec<f32>> = vec![Vec::new(); bins.len()];
    let bin_index: Vec<Option<usize>> = motif_len_of
        .iter()
        .enumerate()
        .map(|(i, &mlen)| {
            if eligible[i] {
                bin_of(&bins, mlen as usize)
            } else {
                None
            }
        })
        .collect();
    for (i, b) in bin_index.iter().enumerate() {
        if let Some(b) = b {
            per_bin[*b].push(std_of[i]);
        }
    }
    per_bin.par_iter_mut().for_each(|v| {
        v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    });

    // Pass 2: re-read the temporary file and append the bin label and index.
    let out_file = match std::fs::File::create(&config.output) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("ERROR: Failed to create {}: {}", config.output.display(), e);
            std::process::exit(1);
        }
    };
    // A genome-wide table is large enough that compressing it is usually the right default;
    // inquiSTR's reader decompresses transparently, so a `.gz` output stays readable everywhere.
    let sink: Box<dyn Write> = if config.output.extension().is_some_and(|e| e == "gz") {
        Box::new(flate2::write::GzEncoder::new(out_file, flate2::Compression::default()))
    } else {
        Box::new(out_file)
    };
    let mut out = std::io::BufWriter::with_capacity(1 << 20, sink);

    writeln!(out, "# file_type=locus_stats").unwrap();
    writeln!(out, "# version={}", crate::VERSION).unwrap();
    writeln!(out, "# command=varindex").unwrap();
    writeln!(out, "# source={}", combined_path).unwrap();
    writeln!(out, "# samples={}", n_samples).unwrap();
    writeln!(out, "# loci={}", processed).unwrap();
    writeln!(out, "# length_axis=absolute_bp").unwrap();
    writeln!(out, "# percentile_method=nearest_rank").unwrap();
    writeln!(out, "# index_column=lvi_length_internal").unwrap();
    writeln!(out, "# min_alleles={}", config.min_alleles).unwrap();
    writeln!(out, "# min_bin_size={}", config.min_bin_size).unwrap();
    writeln!(
        out,
        "# motif_bins={}",
        bins.iter()
            .map(|&(lo, hi)| if lo == hi {
                format!("{}", lo)
            } else {
                format!("{}-{}", lo, hi)
            })
            .collect::<Vec<_>>()
            .join(",")
    )
    .unwrap();
    write_header(&mut out).unwrap();

    let tmp_reader = std::io::BufReader::with_capacity(1 << 20, tmp.reopen().unwrap());
    for (i, row) in tmp_reader.lines().enumerate() {
        let row = row.expect("Failed reading temporary statistics");
        match bin_index.get(i).copied().flatten() {
            Some(b) => {
                let (lo, hi) = bins[b];
                let label = if lo == hi {
                    format!("{}", lo)
                } else {
                    format!("{}-{}", lo, hi)
                };
                let sorted = &per_bin[b];
                let value = std_of[i];
                let less = sorted.partition_point(|&x| x < value);
                let index = less as f64 / sorted.len() as f64;
                writeln!(out, "{}\t{}\t{}", row, label, index).unwrap();
            }
            None => writeln!(out, "{}\t.\tNaN", row).unwrap(),
        }
    }
    out.flush().unwrap();
    eprintln!("Wrote {}", config.output.display());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nearest_rank_is_an_observed_value() {
        let v = vec![1.0, 1.0, 2.0, 5.0, 5.0, 5.0, 9.0];
        assert_eq!(nearest_rank(&v, 0.0), 1.0);
        assert_eq!(nearest_rank(&v, 100.0), 9.0);
        assert_eq!(nearest_rank(&v, 50.0), 5.0);
        // Every returned value must exist in the data - no interpolation
        for (p, _) in PERCENTILES.iter() {
            assert!(v.contains(&nearest_rank(&v, *p)), "percentile {} interpolated", p);
        }
    }

    #[test]
    fn test_nearest_rank_empty() {
        assert!(nearest_rank(&[], 50.0).is_nan());
    }

    #[test]
    fn test_mode_and_unique() {
        let v = vec![1.0, 2.0, 2.0, 3.0, 3.0, 3.0];
        assert_eq!(mode_sorted(&v), 3.0);
        assert_eq!(count_unique(&v), 3);
        // Ties resolve to the smallest value
        let tie = vec![1.0, 1.0, 4.0, 4.0];
        assert_eq!(mode_sorted(&tie), 1.0);
    }

    #[test]
    fn test_median_and_mad() {
        let mut v = vec![10.0, 10.0, 12.0, 14.0, 100.0];
        let s = summarise("chr1", 0, 10, "AC", &mut v, 0);
        assert_eq!(s.n, 5);
        assert_eq!(s.percentiles[11], 12.0); // 50th percentile
        assert_eq!(s.mad, 2.0); // deviations 2,2,0,2,88 -> median 2
    }

    #[test]
    fn test_parse_motif() {
        assert_eq!(parse_motif("CAG"), Some("CAG"));
        assert_eq!(parse_motif("GAAAT,AAAAT"), Some("GAAAT"));
        assert_eq!(parse_motif("."), None);
        assert_eq!(parse_motif("1"), None);
        assert_eq!(parse_motif(""), None);
        assert_eq!(parse_motif("ID=locus1"), None);
    }

    #[test]
    fn test_absolute_length_conversion() {
        // reference length 100, genotypes +10 and -5 -> absolute 110 and 95
        let line = "chr1\t1000\t1100\tCAG\t10\t-5\tNaN";
        let s = summarise_line(line).unwrap();
        assert_eq!(s.n, 2, "NaN must be excluded, not coerced to zero");
        assert_eq!(s.percentiles[0], 95.0);
        assert_eq!(s.percentiles[23], 110.0);
    }

    #[test]
    fn test_negative_absolute_length_is_clamped_and_counted() {
        // reference length 10, a deletion of 50 bp cannot leave a negative repeat
        let line = "chr1\t1000\t1010\tCAG\t-50\t0";
        let s = summarise_line(line).unwrap();
        assert_eq!(s.clamped, 1);
        assert_eq!(s.percentiles[0], 0.0);
    }

    #[test]
    fn test_bin_motif_lengths_merges_sparse_tail() {
        let counts: BTreeMap<usize, usize> = [
            (1, 5000),
            (2, 5000),
            (3, 5000),
            (7, 10),
            (9, 10),
            (12, 5),
            (40, 2),
        ]
        .into_iter()
        .collect();
        let bins = bin_motif_lengths(&counts, 1000);
        // Common lengths keep their own bin; the sparse tail merges into the last one
        assert_eq!(bins[0], (1, 1));
        assert_eq!(bins[1], (2, 2));
        assert_eq!(bins.last().unwrap().1, 40);
        // Every motif length present must land in exactly one bin
        for len in counts.keys() {
            assert!(bin_of(&bins, *len).is_some(), "motif length {} unbinned", len);
        }
    }

    #[test]
    fn test_bin_never_leaves_an_unrankable_remainder() {
        let counts: BTreeMap<usize, usize> = [(1, 100), (2, 3)].into_iter().collect();
        let bins = bin_motif_lengths(&counts, 1000);
        // Nothing reaches the floor, so everything must still be covered by one bin
        assert_eq!(bins, vec![(1, 2)]);
    }
}
