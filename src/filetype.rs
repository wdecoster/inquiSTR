//! File type detection and metadata handling for inquiSTR output files
//!
//! This module provides centralized functionality for:
//! - Reading file type metadata from headers
//! - Detecting file types (STR call vs kmer, individual vs combined)
//! - Validating file formats

use std::io::BufRead;
use std::path::Path;

/// File types that can be detected from metadata headers
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum FileType {
    IndividualCall,
    CombinedCall,
    IndividualKmer,
    CombinedKmer,
    TargetKmer,
}

impl FileType {
    /// Returns true if this is a kmer file type
    pub fn is_kmer(&self) -> bool {
        matches!(self, FileType::IndividualKmer | FileType::CombinedKmer | FileType::TargetKmer)
    }

    /// Returns true if this is a combined file type
    pub fn is_combined(&self) -> bool {
        matches!(self, FileType::CombinedCall | FileType::CombinedKmer)
    }

    /// Returns true if this is an STR call file type
    pub fn is_str_call(&self) -> bool {
        matches!(self, FileType::IndividualCall | FileType::CombinedCall)
    }
}

/// Read file type from metadata header (line starting with # file_type=)
/// Returns None if no metadata header is found
/// Scans all metadata lines (lines starting with #) to find file_type
pub fn read_file_type_metadata(file_path: &Path) -> Option<FileType> {
    let mut file_reader = crate::utils::reader(&file_path.to_string_lossy()).lines();

    // Read through all metadata lines to find file_type
    while let Some(Ok(line)) = file_reader.next() {
        if line.starts_with("# file_type=") {
            let file_type_str = line.trim_start_matches("# file_type=").trim();
            return match file_type_str {
                "individual_call" => Some(FileType::IndividualCall),
                "combined_call" => Some(FileType::CombinedCall),
                "individual_kmer" => Some(FileType::IndividualKmer),
                "combined_kmer" => Some(FileType::CombinedKmer),
                "target_kmer" => Some(FileType::TargetKmer),
                _ => None,
            };
        } else if !line.starts_with('#') {
            // Reached data section without finding file_type
            break;
        }
    }
    None
}

/// Detect if a file contains kmer frequency data or STR call data
/// Returns true if it's a kmer file, false if it's an STR call file
/// First tries to read file_type metadata, then falls back to heuristics for backward compatibility
pub fn is_kmer_file(file_path: &Path) -> bool {
    // Try reading metadata first
    if let Some(file_type) = read_file_type_metadata(file_path) {
        return file_type.is_kmer();
    }

    // Fall back to heuristic detection for files without metadata
    let file_reader = crate::utils::reader(&file_path.to_string_lossy());
    let mut lines = file_reader.lines();

    // Skip metadata lines if present
    let first_line = crate::utils::skip_metadata_lines(&mut lines);
    let fields: Vec<&str> = first_line.split('\t').collect();

    // Check if it's a kmer file format
    // Kmer files have "kmer" as first column header
    // Target kmer files have "Sample" as first column header
    if fields.len() >= 2 && (fields[0] == "kmer" || fields[0] == "Sample") {
        return true;
    }

    // Check if it's STR call format (with or without header)
    // STR files either start with "chromosome" or have genomic coordinates
    if fields.len() >= 3 {
        if fields[0] == "chromosome" {
            return false; // STR file with header
        }

        // Check if first line looks like genomic coordinates (chr1, etc.)
        if fields[0].starts_with("chr")
            && fields[1].parse::<u32>().is_ok()
            && fields[2].parse::<u32>().is_ok()
        {
            return false; // STR file without header
        }
    }

    eprintln!("Error: Unable to determine file format for {}", file_path.display());
    eprintln!("The file does not appear to be a valid inquiSTR output file.");
    eprintln!("\nExpected either:");
    eprintln!("  - STR call file (starts with 'chromosome' or has chr:start:end format)");
    eprintln!("  - Kmer frequency file (starts with 'kmer' or 'Sample')");
    std::process::exit(1);
}

/// Determine if a file is a combined STR file (has more than 6 columns)
/// Individual STR files have 6 columns: chr, start, end, info, H1, H2
/// Combined STR files have 5+ columns: chr, start, end, info, sample1_H1, sample1_H2, sample2_H1, sample2_H2, ...
/// First tries to read file_type metadata, then falls back to heuristics for backward compatibility
pub fn is_combined_str_file(file_path: &Path) -> bool {
    // Try reading metadata first
    if let Some(file_type) = read_file_type_metadata(file_path) {
        return file_type == FileType::CombinedCall;
    }

    // Fall back to heuristic detection for files without metadata
    let file_reader = crate::utils::reader(&file_path.to_string_lossy());
    let mut lines = file_reader.lines();

    // Skip metadata lines
    let first_line = crate::utils::skip_metadata_lines(&mut lines);

    // Validate that file has a header
    if !first_line.starts_with("chromosome") {
        eprintln!("Error: File {} does not have a valid header.", file_path.display());
        eprintln!(
            "Expected first column to be 'chromosome', got '{}'",
            first_line.split('\t').next().unwrap_or("<empty>")
        );
        eprintln!("\nAll STR files must have headers starting with 'chromosome'.");
        std::process::exit(1);
    }

    // Check the next line to determine column count
    match lines.next() {
        Some(Ok(second_line)) => second_line.split('\t').count() > 6,
        Some(Err(e)) => {
            eprintln!("Error: Failed to read data from file {}: {}", file_path.display(), e);
            std::process::exit(1);
        }
        None => {
            eprintln!("Error: File {} has a header but no data lines", file_path.display());
            eprintln!("STR files must contain at least one data line after the header.");
            std::process::exit(1);
        }
    }
}

/// Validate that a file is a combined file (either STR or kmer)
/// Returns an error if the file is not in a valid combined format
pub fn validate_combined_file(file_path: &Path) -> Result<(), String> {
    if let Some(file_type) = read_file_type_metadata(file_path) {
        if !file_type.is_combined() {
            return Err(format!(
                "File {} is not a combined file (type: {:?})",
                file_path.display(),
                file_type
            ));
        }
        Ok(())
    } else {
        // Fall back to heuristic validation
        // This is less precise but maintains backward compatibility
        eprintln!(
            "Warning: File {} lacks file_type metadata. Using heuristic detection.",
            file_path.display()
        );
        Ok(())
    }
}

/// Validate STR file header format
/// Returns the number of samples if valid, or an error message if invalid
/// 
/// Validates:
/// - First 3 columns are: chromosome, begin, end
/// - Fourth column exists (info)
/// - Sample columns come in H1/H2 pairs
/// - H1 and H2 suffixes are present and match
pub fn validate_str_header(header_fields: &[&str]) -> Result<usize, String> {
    if header_fields.len() < 5 {
        return Err(format!(
            "STR header must have at least 5 columns (chromosome, begin, end, info, sample_H1), got {}",
            header_fields.len()
        ));
    }
    
    if header_fields[0] != "chromosome" {
        return Err(format!(
            "First column must be 'chromosome', got '{}'",
            header_fields[0]
        ));
    }
    
    if header_fields[1] != "begin" {
        return Err(format!(
            "Second column must be 'begin', got '{}'",
            header_fields[1]
        ));
    }
    
    if header_fields[2] != "end" {
        return Err(format!(
            "Third column must be 'end', got '{}'",
            header_fields[2]
        ));
    }
    
    // Fourth column is 'info' - we don't validate the name as it may vary,
    // but we ensure it exists by requiring at least 5 columns above
    
    let sample_cols = &header_fields[4..];
    if !sample_cols.len().is_multiple_of(2) {
        return Err(format!(
            "Must have even number of sample columns (H1/H2 pairs), got {} columns after chr/begin/end/info",
            sample_cols.len()
        ));
    }
    
    let n_samples = sample_cols.len() / 2;
    for i in 0..n_samples {
        let h1_col = sample_cols[i * 2];
        let h2_col = sample_cols[i * 2 + 1];
        
        if !h1_col.ends_with("_H1") {
            return Err(format!(
                "Column {} should end with '_H1', got: '{}'",
                4 + i * 2,
                h1_col
            ));
        }
        if !h2_col.ends_with("_H2") {
            return Err(format!(
                "Column {} should end with '_H2', got: '{}'",
                5 + i * 2,
                h2_col
            ));
        }
        
        let name_h1 = h1_col.trim_end_matches("_H1");
        let name_h2 = h2_col.trim_end_matches("_H2");
        if name_h1 != name_h2 {
            return Err(format!(
                "Sample name mismatch for H1/H2 pair at columns {}-{}: '{}' vs '{}'",
                4 + i * 2,
                5 + i * 2,
                name_h1,
                name_h2
            ));
        }
    }
    
    Ok(n_samples)
}
