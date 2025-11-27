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

/// Search for a specific locus in a combined inquiSTR file
/// Returns the first matching locus, or None if not found
pub fn find_locus(config: LocusSearchConfig) -> Option<LocusMatch> {
    if !config.combined_file.exists() {
        eprintln!("ERROR: Combined file does not exist: {}", config.combined_file.display());
        std::process::exit(1);
    }

    let (target_chrom, target_start, target_end) =
        crate::utils::process_region(config.target_region).unwrap();

    let file = crate::utils::reader(&config.combined_file.to_string_lossy());
    let mut lines = file.lines();
    let reg_chrom = format!("{target_chrom}\t");

    // Read header to determine expected number of columns
    // Skip metadata lines if present and get actual header/first data line
    let first_line = crate::utils::skip_metadata_lines(&mut lines);
    let first_cols = first_line.split('\t').count();
    if first_cols < 4 {
        eprintln!("ERROR: Invalid file format. Expected at least 4 columns, got {}.", first_cols);
        eprintln!("First line: '{}'", first_line);
        std::process::exit(1);
    }

    // If this line matches our target chromosome, process it
    if first_line.starts_with(&reg_chrom) {
        let splitline: Vec<&str> = first_line.split('\t').collect();
        let begin: u32 = splitline[1].parse().expect("Failed parsing interval start");
        let end: u32 = splitline[2].parse().expect("Failed parsing interval end");

        // Check overlap based on strategy
        let matches = match config.overlap_strategy {
            OverlapStrategy::Overlap => {
                std::cmp::max(target_start, begin) < std::cmp::min(target_end, end)
            }
            OverlapStrategy::Containment => target_start <= begin && end <= target_end,
        };

        if matches {
            let values: Vec<f64> = splitline
                .iter()
                .skip(4) // Skip chr, start, end, info
                .map(|number| number.parse::<f64>().expect("Failed parsing lengths"))
                .collect();

            return Some(LocusMatch {
                chromosome: target_chrom.clone(),
                start: begin,
                end,
                values,
                raw_line: first_line,
            });
        }
    }

    let expected_columns = first_cols;

    for (line_num, line_result) in lines.enumerate() {
        let line = line_result.unwrap();

        // Skip lines that don't start with our target chromosome
        if !line.starts_with(&reg_chrom) {
            continue;
        }

        let splitline: Vec<&str> = line.split('\t').collect();
        if splitline.len() != expected_columns {
            eprintln!("ERROR: Malformed line {} in combined file.", line_num + 2);
            eprintln!("Expected {} columns, got {}", expected_columns, splitline.len());
            eprintln!("Line content: '{}'", line);
            std::process::exit(1);
        }

        let begin: u32 = splitline[1].parse().expect("Failed parsing interval start");
        let end: u32 = splitline[2].parse().expect("Failed parsing interval end");

        // Check overlap based on strategy
        let matches = match config.overlap_strategy {
            OverlapStrategy::Overlap => {
                // True genomic overlap
                std::cmp::max(target_start, begin) < std::cmp::min(target_end, end)
            }
            OverlapStrategy::Containment => {
                // Target region must contain the found region
                target_start <= begin && end <= target_end
            }
        };

        if matches {
            let values: Vec<f64> = splitline
                .iter()
                .skip(4) // Skip chr, start, end, info
                .map(|number| number.parse::<f64>().expect("Failed parsing lengths"))
                .collect();

            return Some(LocusMatch {
                chromosome: target_chrom.clone(),
                start: begin,
                end,
                values,
                raw_line: line,
            });
        }
    }

    None
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
        .skip(3)
        .map(|s| s.to_string())
        .collect()
}

/// Utility to get clean sample names (removing _H1/_H2 suffixes)  
pub fn extract_clean_sample_names(header_line: &str) -> Vec<String> {
    header_line
        .split('\t')
        .skip(3)
        .map(|s| s.replace("_H1", "").replace("_H2", ""))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_find_locus_overlap_strategy() {
        // This test would require actual test data
        // For now, just test that the enum variants work
        let _overlap = OverlapStrategy::Overlap;
        let _containment = OverlapStrategy::Containment;
    }

    #[test]
    fn test_extract_sample_names() {
        let header = "chromosome\tbegin\tend\tsample1_H1\tsample1_H2\tsample2_H1\tsample2_H2";
        let samples = extract_sample_names(header);
        assert_eq!(samples, vec!["sample1_H1", "sample1_H2", "sample2_H1", "sample2_H2"]);
    }

    #[test]
    fn test_extract_clean_sample_names() {
        let header = "chromosome\tbegin\tend\tsample1_H1\tsample1_H2\tsample2_H1\tsample2_H2";
        let samples = extract_clean_sample_names(header);
        assert_eq!(samples, vec!["sample1", "sample1", "sample2", "sample2"]);
    }
}
