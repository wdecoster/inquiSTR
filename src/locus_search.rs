//! # Locus Search Utilities
//!
//! Common functionality for searching and extracting data from specific genomic loci
//! in inquiSTR combined files. Used by query, plot, and histogram subcommands.

use std::io::BufRead;
use std::path::PathBuf;

/// Overlap detection strategy for locus matching
#[derive(Debug, Clone, Copy)]
pub enum OverlapStrategy {
    /// True genomic overlap: max(start1, start2) < min(end1, end2)
    Overlap,
    /// Target region must be fully contained within query region
    Containment,
}

/// Result of finding a locus in the combined file
#[derive(Debug, Clone)]
pub struct LocusMatch {
    /// Chromosome name
    pub chromosome: String,
    /// Start position
    pub start: u32,
    /// End position  
    pub end: u32,
    /// All sample values (as parsed f64)
    pub values: Vec<f64>,
    /// The raw line content
    pub raw_line: String,
}

/// Configuration for locus search
pub struct LocusSearchConfig {
    /// Combined file path
    pub combined_file: PathBuf,
    /// Target region string (chr:start-end)
    pub target_region: String,
    /// How to detect overlaps
    pub overlap_strategy: OverlapStrategy,
}

/// Scan a combined inquiSTR file for loci matching the target region.
///
/// When `stop_at_first` is true the scan returns as soon as the first match is
/// found (cheap lookups); otherwise it returns every matching locus.
fn scan_loci(config: LocusSearchConfig, stop_at_first: bool) -> Vec<LocusMatch> {
    if !config.combined_file.exists() {
        eprintln!("ERROR: Combined file does not exist: {}", config.combined_file.display());
        std::process::exit(1);
    }

    let (target_chrom, target_start, target_end) =
        crate::utils::process_region(config.target_region).unwrap();

    let file = crate::utils::reader(&config.combined_file.to_string_lossy());
    let mut lines = file.lines();
    let reg_chrom = format!("{target_chrom}\t");

    // Read header / first data line to determine the expected column count.
    // Skip metadata lines (those starting with '#') if present.
    let first_line =
        crate::utils::skip_metadata_lines(&mut lines, &config.combined_file.to_string_lossy());
    let expected_columns = first_line.split('\t').count();
    if expected_columns < 4 {
        eprintln!(
            "ERROR: Invalid file format. Expected at least 4 columns, got {}.",
            expected_columns
        );
        eprintln!("First line: '{}'", first_line);
        std::process::exit(1);
    }

    let mut results = Vec::new();

    // Treat the first data line and the remaining lines uniformly.
    let all_lines = std::iter::once(first_line).chain(lines.map(|l| l.unwrap()));
    for (idx, line) in all_lines.enumerate() {
        // Skip lines that don't start with our target chromosome.
        if !line.starts_with(&reg_chrom) {
            continue;
        }

        let splitline: Vec<&str> = line.split('\t').collect();
        if splitline.len() != expected_columns {
            eprintln!("ERROR: Malformed line {} in combined file.", idx + 1);
            eprintln!("Expected {} columns, got {}", expected_columns, splitline.len());
            eprintln!("Line content: '{}'", line);
            std::process::exit(1);
        }

        let begin: u32 = splitline[1].parse().expect("Failed parsing interval start");
        let end: u32 = splitline[2].parse().expect("Failed parsing interval end");

        let matches = match config.overlap_strategy {
            // True genomic overlap.
            OverlapStrategy::Overlap => {
                std::cmp::max(target_start, begin) < std::cmp::min(target_end, end)
            }
            // Target region must fully contain the found region.
            OverlapStrategy::Containment => target_start <= begin && end <= target_end,
        };

        if matches {
            let values: Vec<f64> = splitline
                .iter()
                .skip(crate::filetype::STR_FIXED_COLUMNS) // Skip chr, start, end, info
                .map(|number| number.parse::<f64>().expect("Failed parsing lengths"))
                .collect();

            results.push(LocusMatch {
                chromosome: target_chrom.clone(),
                start: begin,
                end,
                values,
                raw_line: line,
            });

            if stop_at_first {
                break;
            }
        }
    }

    results
}

/// Search for a specific locus in a combined inquiSTR file.
/// Returns the first matching locus, or None if not found.
pub fn find_locus(config: LocusSearchConfig) -> Option<LocusMatch> {
    scan_loci(config, true).into_iter().next()
}

/// Return every locus in the combined file matching the target region.
/// Unlike [`find_locus`], this does not stop at the first match, so callers
/// can detect when a region spans more than one locus.
pub fn find_overlapping_loci(config: LocusSearchConfig) -> Vec<LocusMatch> {
    scan_loci(config, false)
}

/// Resolve the single locus overlapping the requested region for an
/// interactive (single-region) command such as `plot` or `histogram`.
///
/// Exits with a clear message if no locus overlaps the region, and prints a
/// warning (then uses the first locus) if the region spans more than one.
pub fn select_overlapping_locus(config: LocusSearchConfig) -> LocusMatch {
    let region = config.target_region.clone();
    let combined_file = config.combined_file.clone();
    let mut matches = find_overlapping_loci(config);

    match matches.len() {
        0 => {
            eprintln!(
                "ERROR: No locus overlapping region '{}' was found in {}.\n\
                 Check the coordinates (expected chrom:begin-end) and that a locus in this \
                 interval is present in the combined file.",
                region,
                combined_file.display()
            );
            std::process::exit(1);
        }
        1 => matches.swap_remove(0),
        n => {
            let coords: Vec<String> = matches
                .iter()
                .map(|m| format!("{}:{}-{}", m.chromosome, m.start, m.end))
                .collect();
            eprintln!(
                "WARNING: region '{}' overlaps {} loci ({}). Using the first; narrow the \
                 region to select a specific locus.",
                region,
                n,
                coords.join(", ")
            );
            matches.swap_remove(0)
        }
    }
}

/// Search for multiple loci in a combined inquiSTR file
/// Supports both single region strings and files with multiple regions
pub fn find_multiple_loci(
    combined_file: PathBuf,
    region_input: String,
    overlap_strategy: OverlapStrategy,
) -> Vec<LocusMatch> {
    let mut results = Vec::new();

    // Determine if region_input is a file or a region string
    let target_regions = if std::path::Path::new(&region_input).exists() {
        // It's a file, read regions from it
        let mut regions = Vec::new();
        for line in crate::utils::reader(&region_input).lines() {
            regions.push(line.unwrap());
        }
        regions
    } else {
        // It's a single region string
        vec![region_input]
    };

    // Search for each target region
    for target_region in target_regions {
        let config = LocusSearchConfig {
            combined_file: combined_file.clone(),
            target_region,
            overlap_strategy,
        };

        if let Some(locus_match) = find_locus(config) {
            results.push(locus_match);
        }
    }

    results
}

/// Utility to get sample names from header line
pub fn extract_sample_names(header_line: &str) -> Vec<String> {
    header_line
        .split('\t')
        .skip(crate::filetype::STR_FIXED_COLUMNS) // Skip chromosome, begin, end, info
        .map(|s| s.to_string())
        .collect()
}

/// Utility to get clean sample names (removing _H1/_H2 suffixes)  
pub fn extract_clean_sample_names(header_line: &str) -> Vec<String> {
    header_line
        .split('\t')
        .skip(crate::filetype::STR_FIXED_COLUMNS) // Skip chromosome, begin, end, info
        .map(|s| s.replace("_H1", "").replace("_H2", ""))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Write a combined-file body to a temporary `.tsv` and return the handle
    /// (kept alive by the caller so the file isn't deleted while in use).
    fn temp_combined(body: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        f.write_all(body.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    fn config(
        path: &std::path::Path,
        region: &str,
        strategy: OverlapStrategy,
    ) -> LocusSearchConfig {
        LocusSearchConfig {
            combined_file: path.to_path_buf(),
            target_region: region.to_string(),
            overlap_strategy: strategy,
        }
    }

    #[test]
    fn test_find_overlapping_loci_single_multiple_none() {
        let body = "chromosome\tbegin\tend\tinfo\ts1_H1\ts1_H2\n\
                    chr1\t100\t200\t.\t10\t12\n\
                    chr1\t400\t500\t.\t20\t22\n\
                    chr1\t900\t1000\t.\t30\t32\n";
        let f = temp_combined(body);
        let path = f.path();

        // Region overlapping two loci.
        let multi = find_overlapping_loci(config(path, "chr1:150-450", OverlapStrategy::Overlap));
        assert_eq!(multi.len(), 2);
        assert_eq!((multi[0].start, multi[0].end), (100, 200));
        assert_eq!((multi[1].start, multi[1].end), (400, 500));

        // find_locus returns only the first overlapping locus.
        let first = find_locus(config(path, "chr1:150-450", OverlapStrategy::Overlap)).unwrap();
        assert_eq!((first.start, first.end), (100, 200));
        assert_eq!(first.values, vec![10.0, 12.0]);

        // Region overlapping exactly one locus.
        let one = find_overlapping_loci(config(path, "chr1:120-130", OverlapStrategy::Overlap));
        assert_eq!(one.len(), 1);

        // Region overlapping no locus.
        let none = find_overlapping_loci(config(path, "chr1:600-700", OverlapStrategy::Overlap));
        assert!(none.is_empty());

        // Different chromosome.
        let other = find_overlapping_loci(config(path, "chr2:100-200", OverlapStrategy::Overlap));
        assert!(other.is_empty());
    }

    #[test]
    fn test_containment_vs_overlap_strategy() {
        let body = "chromosome\tbegin\tend\tinfo\ts1\nchr1\t100\t200\t.\t10\n";
        let f = temp_combined(body);
        let path = f.path();

        // Partial overlap: Containment misses, Overlap hits.
        assert!(find_locus(config(path, "chr1:150-250", OverlapStrategy::Containment)).is_none());
        assert!(find_locus(config(path, "chr1:150-250", OverlapStrategy::Overlap)).is_some());

        // Region fully containing the locus: both strategies hit.
        assert!(find_locus(config(path, "chr1:50-250", OverlapStrategy::Containment)).is_some());
        assert!(find_locus(config(path, "chr1:50-250", OverlapStrategy::Overlap)).is_some());
    }

    #[test]
    fn test_find_locus_overlap_strategy() {
        // For now, just test that the enum variants work
        let _overlap = OverlapStrategy::Overlap;
        let _containment = OverlapStrategy::Containment;
    }

    #[test]
    fn test_extract_sample_names() {
        let header = "chromosome\tbegin\tend\tinfo\tsample1_H1\tsample1_H2\tsample2_H1\tsample2_H2";
        let samples = extract_sample_names(header);
        assert_eq!(samples, vec!["sample1_H1", "sample1_H2", "sample2_H1", "sample2_H2"]);
    }

    #[test]
    fn test_extract_clean_sample_names() {
        let header = "chromosome\tbegin\tend\tinfo\tsample1_H1\tsample1_H2\tsample2_H1\tsample2_H2";
        let samples = extract_clean_sample_names(header);
        assert_eq!(samples, vec!["sample1", "sample1", "sample2", "sample2"]);
    }
}
