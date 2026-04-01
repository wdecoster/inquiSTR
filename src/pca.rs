//! # Principal Component Analysis (PCA) for STR Genotypes and Kmer Frequencies
//!
//! This module implements Principal Component Analysis for Short Tandem Repeat (STR) genotype data
//! and kmer frequency data from inquiSTR combined files. It provides dimensionality reduction and
//! visualization to identify population structure and relationships between samples.
//!
//! ## Features
//!
//! - **Automated data parsing** from inquiSTR combined files (both STR and kmer formats)
//! - **Simplified PCA implementation** using variance-based feature selection
//! - **Interactive HTML plots** using Plotly for visualization
//! - **Support for missing data** (NaN values handled gracefully)
//!
//! ## Usage
//!
//! ```bash
//! # Generate PCA plot from combined STR data
//! inquiSTR pca combined_strs.tsv --output str_pca.html
//!
//! # Generate PCA plot from combined kmer frequency data
//! inquiSTR pca inquiSTR_unmapped.tsv --output kmer_pca.html
//! ```
//!
//! ## Implementation
//!
//! This is a simplified PCA implementation that uses the top variance features as principal
//! components. For production use with large datasets, consider using more sophisticated
//! eigenvalue decomposition methods.

use kuva::prelude::*;
use ndarray::{Array1, Array2, Axis};
use rayon::prelude::*;
use std::io::BufRead;
use std::path::PathBuf;

/// Parse a numeric field, treating NaN and empty as 0.0.
/// Exits with an error for genuinely non-numeric values.
fn parse_numeric_or_nan(field: &str, col: usize) -> f64 {
    if field.eq_ignore_ascii_case("nan") || field.is_empty() {
        0.0
    } else {
        field.parse().unwrap_or_else(|_| {
            eprintln!("ERROR: Non-numeric value '{}' in column {} of input file", field, col + 1);
            std::process::exit(1);
        })
    }
}

/// Methods for aggregating H1 and H2 allele lengths for PCA analysis
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleAggregation {
    /// Use maximum allele length (default - emphasizes longer alleles)
    Max,
    /// Use minimum allele length (emphasizes shorter alleles)
    Min,
    /// Use sum of allele lengths (emphasizes total repeat content)
    Sum,
}

impl std::fmt::Display for AlleleAggregation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlleleAggregation::Max => write!(f, "max"),
            AlleleAggregation::Min => write!(f, "min"),
            AlleleAggregation::Sum => write!(f, "sum"),
        }
    }
}

impl std::str::FromStr for AlleleAggregation {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "max" => Ok(AlleleAggregation::Max),
            "min" => Ok(AlleleAggregation::Min),
            "sum" => Ok(AlleleAggregation::Sum),
            _ => Err(format!("Invalid aggregation method: '{}'. Valid options: max, min, sum", s)),
        }
    }
}

impl AlleleAggregation {
    /// Apply the aggregation method to H1 and H2 allele lengths
    pub fn aggregate(&self, h1_val: f64, h2_val: f64) -> f64 {
        match self {
            AlleleAggregation::Max => h1_val.max(h2_val),
            AlleleAggregation::Min => h1_val.min(h2_val),
            AlleleAggregation::Sum => h1_val + h2_val,
        }
    }
}

/// Parse a combined kmer file for PCA analysis
fn parse_combined_kmer_file_with_selection(
    combined: &std::path::Path,
    max_features: Option<usize>,
) -> (Array2<f64>, Vec<String>) {
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();

    // Read header line - skip metadata if present
    let header_line = crate::utils::skip_metadata_lines(&mut lines);
    let header_fields: Vec<&str> = header_line.trim().split('\t').collect();

    // Validate kmer file header format
    if header_fields.len() < 2 || header_fields[0] != "kmer" {
        eprintln!("ERROR: Invalid combined kmer file header.");
        eprintln!("Expected format: kmer\tsample1\tsample2\t...");
        std::process::exit(1);
    }

    // Extract sample names (all columns after "kmer")
    let sample_names: Vec<String> = header_fields[1..].iter().map(|s| s.to_string()).collect();
    let num_samples = sample_names.len();

    if num_samples < 2 {
        eprintln!("ERROR: Need at least 2 samples for PCA, found {}", num_samples);
        std::process::exit(1);
    }

    println!("Detected {} samples in kmer file", num_samples);

    // For kmer files, we may have many features (kmers)
    // Use feature selection if requested or if we detect many features
    if let Some(max_features) = max_features {
        return parse_kmer_with_feature_selection(
            combined,
            num_samples,
            sample_names,
            max_features,
        );
    }

    // Single-pass approach for smaller datasets
    let mut data_rows = Vec::new();

    for (line_num, line_result) in lines.enumerate() {
        let line = line_result
            .map_err(|e| format!("Error reading line {}: {}", line_num + 2, e))
            .expect("IO error reading file");

        let fields: Vec<&str> = line.trim().split('\t').collect();

        let expected_cols = 1 + num_samples;
        if fields.len() != expected_cols {
            eprintln!(
                "ERROR: Malformed kmer line {} (expected {} columns, got {})",
                line_num + 2,
                expected_cols,
                fields.len()
            );
            std::process::exit(1);
        }

        // Parse kmer frequencies for this kmer across all samples
        let mut row_data = Vec::with_capacity(num_samples);
        for sample_idx in 0..num_samples {
            let freq_idx = 1 + sample_idx;
            let freq = parse_numeric_or_nan(fields[freq_idx], freq_idx);
            row_data.push(freq);
        }
        data_rows.push(row_data);
    }

    if data_rows.is_empty() {
        eprintln!("ERROR: No data lines found after header in kmer file.");
        std::process::exit(1);
    }

    let num_kmers = data_rows.len();
    let mut data_matrix = Array2::<f64>::zeros((num_samples, num_kmers));

    // Fill matrix from collected row data
    for (kmer_idx, row_data) in data_rows.iter().enumerate() {
        for (sample_idx, &value) in row_data.iter().enumerate() {
            data_matrix[[sample_idx, kmer_idx]] = value;
        }
    }

    println!("Loaded {} kmers for {} samples", num_kmers, num_samples);

    (data_matrix, sample_names)
}

/// Parse kmer file with feature selection for large datasets
fn parse_kmer_with_feature_selection(
    combined: &std::path::Path,
    num_samples: usize,
    sample_names: Vec<String>,
    max_features: usize,
) -> (Array2<f64>, Vec<String>) {
    println!("Using memory-efficient two-pass parsing for large kmer dataset...");

    // PASS 1: Calculate variance scores for each kmer
    println!("Pass 1: Analyzing kmers to find most informative features...");
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();
    let _header = crate::utils::skip_metadata_lines(&mut lines);

    let mut feature_scores: Vec<(usize, f64)> = Vec::new();

    for (kmer_idx, line_result) in lines.enumerate() {
        let line = line_result.unwrap_or_else(|e| {
            eprintln!("ERROR: Failed to read line {} in kmer file: {}", kmer_idx + 2, e);
            std::process::exit(1);
        });

        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 1 + num_samples {
            eprintln!(
                "ERROR: Malformed kmer line {} (expected {} columns, got {})",
                kmer_idx + 2,
                1 + num_samples,
                fields.len()
            );
            std::process::exit(1);
        }

        // Calculate variance for this kmer
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        let mut count = 0;

        for sample_idx in 0..num_samples {
            let freq = parse_numeric_or_nan(fields[1 + sample_idx], 1 + sample_idx);
            sum += freq;
            sum_sq += freq * freq;
            count += 1;
        }

        if count > 1 {
            let variance = (sum_sq - sum * sum / count as f64) / (count as f64 - 1.0);
            feature_scores.push((kmer_idx, variance));
        }
    }

    // Sort by variance and select top features
    feature_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let selected_indices: Vec<usize> = feature_scores
        .into_iter()
        .take(max_features)
        .map(|(idx, _)| idx)
        .collect();

    println!("Pass 2: Loading {} selected kmers into memory...", selected_indices.len());

    // PASS 2: Load only selected features
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();
    let _header = crate::utils::skip_metadata_lines(&mut lines);

    let mut data_matrix = Array2::<f64>::zeros((num_samples, selected_indices.len()));
    let mut selected_idx_map: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();

    for (new_idx, &original_idx) in selected_indices.iter().enumerate() {
        selected_idx_map.insert(original_idx, new_idx);
    }

    for (kmer_idx, line_result) in lines.enumerate() {
        if let Some(&new_idx) = selected_idx_map.get(&kmer_idx) {
            let line = line_result.unwrap_or_else(|e| {
                eprintln!("ERROR: Failed to read line {} in kmer file: {}", kmer_idx + 2, e);
                std::process::exit(1);
            });

            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() != 1 + num_samples {
                eprintln!(
                    "ERROR: Malformed kmer line {} (expected {} columns, got {})",
                    kmer_idx + 2,
                    1 + num_samples,
                    fields.len()
                );
                std::process::exit(1);
            }

            for sample_idx in 0..num_samples {
                let freq = parse_numeric_or_nan(fields[1 + sample_idx], 1 + sample_idx);
                data_matrix[[sample_idx, new_idx]] = freq;
            }
        }
    }

    println!("Loaded {} kmers for {} samples", selected_indices.len(), num_samples);

    (data_matrix, sample_names)
}

/// Parse a combined inquiSTR file with optional feature pre-selection for memory efficiency
/// This is a dispatcher that routes to the appropriate parser based on file type
fn parse_combined_file_with_selection(
    combined: &std::path::Path,
    max_features: Option<usize>,
    aggregation: AlleleAggregation,
) -> (Array2<f64>, Vec<String>) {
    // Detect file type and route to appropriate parser
    let file_type = crate::filetype::read_file_type_metadata(combined);

    match file_type {
        Some(crate::filetype::FileType::CombinedKmer) => {
            parse_combined_kmer_file_with_selection(combined, max_features)
        }
        Some(crate::filetype::FileType::CombinedCall) => {
            parse_combined_str_file_with_selection(combined, max_features, aggregation)
        }
        _ => {
            // Try to auto-detect from header
            let file = crate::utils::reader(&combined.to_string_lossy());
            let mut lines = file.lines();
            let header_line = crate::utils::skip_metadata_lines(&mut lines);
            let header_fields: Vec<&str> = header_line.trim().split('\t').collect();

            // Check if this looks like a kmer file (first column is "kmer")
            if header_fields.len() >= 2 && header_fields[0] == "kmer" {
                parse_combined_kmer_file_with_selection(combined, max_features)
            // Check if this looks like an STR file (first three columns are chromosome, begin, end)
            } else if header_fields.len() >= 3
                && header_fields[0] == "chromosome"
                && header_fields[1] == "begin"
                && header_fields[2] == "end"
            {
                parse_combined_str_file_with_selection(combined, max_features, aggregation)
            } else {
                // Unable to determine file type
                eprintln!("Error: Unable to determine file type from header.");
                eprintln!("Expected either:");
                eprintln!(
                    "  - STR file: chromosome\\tbegin\\tend\\tinfo\\tsample1_H1\\tsample1_H2..."
                );
                eprintln!("  - Kmer file: kmer\\tsample1\\tsample2...");
                eprintln!("\nGot header: {}", header_line);
                std::process::exit(1);
            }
        }
    }
}

/// Parse a combined STR file for PCA analysis
fn parse_combined_str_file_with_selection(
    combined: &std::path::Path,
    max_features: Option<usize>,
    aggregation: AlleleAggregation,
) -> (Array2<f64>, Vec<String>) {
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();

    // Read header line - this should contain sample names
    // Skip metadata lines if present
    let header_line = crate::utils::skip_metadata_lines(&mut lines);

    let header_fields: Vec<&str> = header_line.trim().split('\t').collect();

    // Use centralized validation function
    let num_samples = match crate::filetype::validate_str_header(&header_fields) {
        Ok(n) => n,
        Err(e) => {
            eprintln!("Error: Invalid STR file header format.");
            eprintln!("{}", e);
            eprintln!(
                "\nExpected: chromosome\\tbegin\\tend\\tinfo\\tsample1_H1\\tsample1_H2\\t..."
            );
            eprintln!("Got: {}", header_line);
            eprintln!("\nThis file does not appear to be a valid combined STR file.");
            eprintln!("If this is a kmer file, it should have been auto-detected.");
            std::process::exit(1);
        }
    };

    // Extract sample names from validated header
    let sample_names: Vec<String> = (0..num_samples)
        .map(|i| {
            let h1_col = header_fields[crate::filetype::STR_FIXED_COLUMNS + i * 2];
            h1_col.trim_end_matches("_H1").to_string()
        })
        .collect();

    // Two-pass approach for large datasets with feature selection
    if let Some(max_features) = max_features {
        return parse_with_feature_selection(
            combined,
            num_samples,
            sample_names,
            max_features,
            aggregation,
        );
    }

    // Single-pass approach for smaller datasets (no feature selection)
    let mut data_rows = Vec::new();

    for (line_num, line_result) in lines.enumerate() {
        let line = line_result
            .map_err(|e| format!("Error reading line {}: {}", line_num + 2, e))
            .expect("IO error reading file");

        let fields: Vec<&str> = line.trim().split('\t').collect();

        let expected_cols = crate::filetype::STR_FIXED_COLUMNS + num_samples * 2;
        if fields.len() != expected_cols {
            eprintln!(
                "ERROR: Malformed line {} (expected {} columns, got {})",
                line_num + 2,
                expected_cols,
                fields.len(),
            );
            std::process::exit(1);
        }

        // Skip region name creation for memory efficiency

        // Parse STR lengths for this region across all samples
        let mut row_data = Vec::with_capacity(num_samples);
        for sample_idx in 0..num_samples {
            let h1_idx = crate::filetype::STR_FIXED_COLUMNS + sample_idx * 2;
            let h2_idx = crate::filetype::STR_FIXED_COLUMNS + 1 + sample_idx * 2;

            let h1_val = parse_numeric_or_nan(fields[h1_idx], h1_idx);
            let h2_val = parse_numeric_or_nan(fields[h2_idx], h2_idx);

            // Apply selected aggregation method for H1/H2 allele lengths
            let aggregated_value = aggregation.aggregate(h1_val, h2_val);
            row_data.push(aggregated_value);
        }
        data_rows.push(row_data);
    }

    if data_rows.is_empty() {
        eprintln!("ERROR: No data lines found after header.");
        eprintln!("The input file appears to have a header but no data.");
        std::process::exit(1);
    }

    let num_regions = data_rows.len();
    // NOTE: Using f64 for precision in PCA calculations. For 1000 samples × 5000 features,
    // switching to f32 would save ~20MB but could impact eigenvalue calculation accuracy
    let mut data_matrix = Array2::<f64>::zeros((num_samples, num_regions));

    // Fill matrix from collected row data
    for (region_idx, row_data) in data_rows.iter().enumerate() {
        for (sample_idx, &value) in row_data.iter().enumerate() {
            data_matrix[[sample_idx, region_idx]] = value;
        }
    }

    (data_matrix, sample_names)
}

/// Parse combined file with intelligent feature pre-selection (two-pass approach)
/// First pass: calculate feature scores, Second pass: load only selected features
fn parse_with_feature_selection(
    combined: &std::path::Path,
    num_samples: usize,
    sample_names: Vec<String>,
    max_features: usize,
    aggregation: AlleleAggregation,
) -> (Array2<f64>, Vec<String>) {
    println!("Using memory-efficient two-pass parsing for large dataset...");

    // PASS 1: Calculate feature scores without storing all data
    println!("Pass 1: Analyzing millions of loci to find most informative features...");
    let selected_indices =
        first_pass_feature_selection(combined, num_samples, max_features, aggregation);

    println!(
        "Pass 2: Loading only {} selected features into memory...",
        selected_indices.len()
    );
    // PASS 2: Load only selected features
    let data_matrix =
        second_pass_data_loading(combined, num_samples, &selected_indices, aggregation);

    (data_matrix, sample_names)
}

/// First pass: scan file to calculate feature scores and select best features
fn first_pass_feature_selection(
    combined: &std::path::Path,
    num_samples: usize,
    max_features: usize,
    aggregation: AlleleAggregation,
) -> Vec<usize> {
    let estimated_lines = estimate_feature_count(combined);

    // For large datasets, use parallel processing
    if estimated_lines > 10_000 {
        return parallel_feature_selection(combined, num_samples, max_features, aggregation);
    }

    // Fall back to sequential for smaller datasets
    sequential_feature_selection(combined, num_samples, max_features, aggregation)
}

/// Sequential feature selection for smaller datasets
fn sequential_feature_selection(
    combined: &std::path::Path,
    num_samples: usize,
    max_features: usize,
    aggregation: AlleleAggregation,
) -> Vec<usize> {
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();

    // Skip metadata lines and header
    let _header = crate::utils::skip_metadata_lines(&mut lines);

    let mut feature_scores: Vec<(usize, f64)> = Vec::new();

    // Add progress tracking for large datasets
    let progress_bar = if estimate_feature_count(combined) > 50_000 {
        let pb = indicatif::ProgressBar::new_spinner();
        pb.set_style(indicatif::ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {pos:>7} loci analyzed, {per_sec} loci/sec, {msg}")
            .unwrap());
        pb.set_message("Scoring features for informativeness...");
        Some(pb)
    } else {
        None
    };

    for (region_idx, line_result) in lines.enumerate() {
        let line = line_result.expect("IO error reading file");
        let fields: Vec<&str> = line.trim().split('\t').collect();

        // Quick validation
        let expected_cols = crate::filetype::STR_FIXED_COLUMNS + num_samples * 2;
        if fields.len() != expected_cols {
            eprintln!(
                "ERROR: Malformed data at line {} (expected {} columns, got {})",
                region_idx + 2,
                expected_cols,
                fields.len()
            );
            eprintln!("The combined file may be corrupted or was not generated correctly.");
            std::process::exit(1);
        }

        // Streaming statistics calculation - single pass, no Vec allocation
        let mut count = 0;
        let mut non_zero_count = 0;
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        let mut min_val = f64::INFINITY;
        let mut max_val = f64::NEG_INFINITY;

        // PHASE 1: Quick pre-filter - check data quality before parsing
        // Count missing values with minimal string operations
        let mut missing_samples = 0;
        let max_allowed_missing = num_samples * 20 / 100; // Allow up to 20% missing

        for sample_idx in 0..num_samples {
            let h1_idx = crate::filetype::STR_FIXED_COLUMNS + sample_idx * 2;
            let h2_idx = crate::filetype::STR_FIXED_COLUMNS + 1 + sample_idx * 2;

            // Quick missing data check (both alleles missing)
            let h1_missing =
                fields[h1_idx].eq_ignore_ascii_case("nan") || fields[h1_idx].is_empty();
            let h2_missing =
                fields[h2_idx].eq_ignore_ascii_case("nan") || fields[h2_idx].is_empty();

            if h1_missing && h2_missing {
                missing_samples += 1;
                // Early termination: if too much missing data, skip entirely
                if missing_samples > max_allowed_missing {
                    break;
                }
            }
        }

        // Skip this locus entirely if too much missing data
        if missing_samples > max_allowed_missing {
            continue;
        }

        // PHASE 2: Full parsing and statistics (only for promising loci)
        for sample_idx in 0..num_samples {
            let h1_idx = crate::filetype::STR_FIXED_COLUMNS + sample_idx * 2;
            let h2_idx = crate::filetype::STR_FIXED_COLUMNS + 1 + sample_idx * 2;

            // Optimized parsing with early NaN check
            let h1_val = parse_numeric_or_nan(fields[h1_idx], h1_idx);
            let h2_val = parse_numeric_or_nan(fields[h2_idx], h2_idx);

            let aggregated_value = aggregation.aggregate(h1_val, h2_val);

            // Update streaming statistics in single pass
            count += 1;
            sum += aggregated_value;
            sum_sq += aggregated_value * aggregated_value;

            if aggregated_value > 0.0 {
                non_zero_count += 1;
            }

            if aggregated_value < min_val {
                min_val = aggregated_value;
            }
            if aggregated_value > max_val {
                max_val = aggregated_value;
            }
        }

        // Skip features with too much missing data (early termination caught some)
        if count == 0 || non_zero_count < num_samples / 2 {
            continue;
        }

        // Calculate statistics from streaming accumulators
        let mean = sum / count as f64;
        let variance = if count > 1 {
            (sum_sq - sum * sum / count as f64) / (count as f64 - 1.0)
        } else {
            0.0
        };

        let cv = if mean > 0.0 {
            variance.sqrt() / mean
        } else {
            0.0
        };
        let range = max_val - min_val;

        // Combined score
        let coverage_score = non_zero_count as f64 / num_samples as f64;
        let combined_score = variance * coverage_score + cv * 0.5 + range * 0.1;

        feature_scores.push((region_idx, combined_score));

        // Update progress bar
        if let Some(ref pb) = progress_bar {
            pb.inc(1);
            if region_idx % 10_000 == 0 {
                let msg = format!("Scored {} loci so far...", region_idx + 1);
                pb.set_message(msg);
            }
        }
    }

    // Finish progress bar
    if let Some(pb) = progress_bar {
        let final_msg = format!("Completed analysis of {} loci", feature_scores.len());
        pb.finish_with_message(final_msg);
    }

    // Sort and select top features
    feature_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let selected_indices: Vec<usize> = feature_scores
        .into_iter()
        .take(max_features)
        .map(|(idx, _)| idx)
        .collect();

    selected_indices
}

/// Parallel feature selection for large datasets using chunked streaming
fn parallel_feature_selection(
    combined: &std::path::Path,
    num_samples: usize,
    max_features: usize,
    aggregation: AlleleAggregation,
) -> Vec<usize> {
    let num_cores = rayon::current_num_threads();
    println!(
        "Using streaming parallel processing for feature selection ({} CPU cores)...",
        num_cores
    );

    // Read file in chunks to avoid loading everything into memory
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();

    // Skip metadata lines and header
    let _header = crate::utils::skip_metadata_lines(&mut lines);

    let mut all_scores = Vec::new();
    let chunk_size = 10_000; // Process 10k lines at a time
    let mut chunk_buffer = Vec::with_capacity(chunk_size);
    let mut processed_count = 0;

    // Process file in chunks
    for (idx, line_result) in lines.enumerate() {
        let line = line_result.unwrap_or_else(|e| {
            eprintln!("ERROR: Failed to read line {}: {}", idx + 2, e);
            std::process::exit(1);
        });
        chunk_buffer.push((idx, line));

        // When chunk is full, process it in parallel
        if chunk_buffer.len() == chunk_size {
            let chunk_scores: Vec<(usize, f64)> = chunk_buffer
                .par_iter()
                .filter_map(|(region_idx, line)| {
                    calculate_feature_score(line, num_samples, *region_idx, aggregation)
                })
                .collect();

            all_scores.extend(chunk_scores);
            processed_count += chunk_buffer.len();

            if processed_count % 50_000 == 0 {
                println!("Processed {} loci so far...", processed_count);
            }

            chunk_buffer.clear();
        }
    }

    // Process remaining lines in the last chunk
    if !chunk_buffer.is_empty() {
        let chunk_scores: Vec<(usize, f64)> = chunk_buffer
            .par_iter()
            .filter_map(|(region_idx, line)| {
                calculate_feature_score(line, num_samples, *region_idx, aggregation)
            })
            .collect();

        all_scores.extend(chunk_scores);

        println!("Completed chunked processing. Total features scored: {}", all_scores.len());
    }

    println!(
        "Scored {} features using chunked parallel processing, selecting top {}",
        all_scores.len(),
        max_features
    );

    // Step 3: Sort and select top features
    all_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let selected_indices: Vec<usize> = all_scores
        .into_iter()
        .take(max_features)
        .map(|(idx, _)| idx)
        .collect();

    selected_indices
}

/// Calculate feature score for a single line (extracted for parallel processing)
fn calculate_feature_score(
    line: &str,
    num_samples: usize,
    region_idx: usize,
    aggregation: AlleleAggregation,
) -> Option<(usize, f64)> {
    let fields: Vec<&str> = line.trim().split('\t').collect();

    // Quick validation
    let expected_cols = crate::filetype::STR_FIXED_COLUMNS + num_samples * 2;
    if fields.len() != expected_cols {
        eprintln!(
            "ERROR: Malformed data at line {} (expected {} columns, got {})",
            region_idx + 2,
            expected_cols,
            fields.len()
        );
        eprintln!("The combined file may be corrupted or was not generated correctly.");
        std::process::exit(1);
    }

    // Streaming statistics calculation - single pass, no Vec allocation
    let mut count = 0;
    let mut non_zero_count = 0;
    let mut sum = 0.0;
    let mut sum_sq = 0.0;
    let mut min_val = f64::INFINITY;
    let mut max_val = f64::NEG_INFINITY;

    // PHASE 1: Quick pre-filter - check data quality before parsing
    let mut missing_samples = 0;
    let max_allowed_missing = num_samples * 20 / 100; // Allow up to 20% missing

    for sample_idx in 0..num_samples {
        let h1_idx = crate::filetype::STR_FIXED_COLUMNS + sample_idx * 2;
        let h2_idx = crate::filetype::STR_FIXED_COLUMNS + 1 + sample_idx * 2;

        // Quick missing data check (both alleles missing)
        let h1_missing = fields[h1_idx].eq_ignore_ascii_case("nan") || fields[h1_idx].is_empty();
        let h2_missing = fields[h2_idx].eq_ignore_ascii_case("nan") || fields[h2_idx].is_empty();

        if h1_missing && h2_missing {
            missing_samples += 1;
            // Early termination: if too much missing data, skip entirely
            if missing_samples > max_allowed_missing {
                return None;
            }
        }
    }

    // Skip this locus entirely if too much missing data
    if missing_samples > max_allowed_missing {
        return None;
    }

    // PHASE 2: Full parsing and statistics (only for promising loci)
    for sample_idx in 0..num_samples {
        let h1_idx = crate::filetype::STR_FIXED_COLUMNS + sample_idx * 2;
        let h2_idx = crate::filetype::STR_FIXED_COLUMNS + 1 + sample_idx * 2;

        // Optimized parsing with early NaN check
        let h1_val = parse_numeric_or_nan(fields[h1_idx], h1_idx);
        let h2_val = parse_numeric_or_nan(fields[h2_idx], h2_idx);

        let max_allele = aggregation.aggregate(h1_val, h2_val);

        // Update streaming statistics in single pass
        count += 1;
        sum += max_allele;
        sum_sq += max_allele * max_allele;

        if max_allele > 0.0 {
            non_zero_count += 1;
        }

        if max_allele < min_val {
            min_val = max_allele;
        }
        if max_allele > max_val {
            max_val = max_allele;
        }
    }

    // Skip features with too much missing data
    if count == 0 || non_zero_count < num_samples / 2 {
        return None;
    }

    // Calculate statistics from streaming accumulators
    let mean = sum / count as f64;
    let variance = if count > 1 {
        (sum_sq - sum * sum / count as f64) / (count as f64 - 1.0)
    } else {
        0.0
    };

    let cv = if mean > 0.0 {
        variance.sqrt() / mean
    } else {
        0.0
    };
    let range = max_val - min_val;

    // Combined score
    let coverage_score = non_zero_count as f64 / num_samples as f64;
    let combined_score = variance * coverage_score + cv * 0.5 + range * 0.1;

    Some((region_idx, combined_score))
}

/// Second pass: load only the selected features into memory
fn second_pass_data_loading(
    combined: &std::path::Path,
    num_samples: usize,
    selected_indices: &[usize],
    aggregation: AlleleAggregation,
) -> Array2<f64> {
    // For datasets with many selected features, use parallel processing
    if selected_indices.len() > 1000 {
        return parallel_data_loading(combined, num_samples, selected_indices, aggregation);
    }

    // Sequential version for smaller feature sets
    sequential_data_loading(combined, num_samples, selected_indices, aggregation)
}

/// Sequential data loading for smaller feature sets
fn sequential_data_loading(
    combined: &std::path::Path,
    num_samples: usize,
    selected_indices: &[usize],
    aggregation: AlleleAggregation,
) -> Array2<f64> {
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();

    // Skip metadata lines and header
    let _header = crate::utils::skip_metadata_lines(&mut lines);

    let mut data_matrix = Array2::<f64>::zeros((num_samples, selected_indices.len()));
    let mut selected_idx_map: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();

    // Create mapping from original index to new index
    for (new_idx, &original_idx) in selected_indices.iter().enumerate() {
        selected_idx_map.insert(original_idx, new_idx);
    }

    for (region_idx, line_result) in lines.enumerate() {
        // Only process if this region was selected
        if let Some(&new_idx) = selected_idx_map.get(&region_idx) {
            let line = line_result.expect("IO error reading file");
            let fields: Vec<&str> = line.trim().split('\t').collect();

            // Parse and store data for this selected region
            for sample_idx in 0..num_samples {
                let h1_idx = crate::filetype::STR_FIXED_COLUMNS + sample_idx * 2;
                let h2_idx = crate::filetype::STR_FIXED_COLUMNS + 1 + sample_idx * 2;

                let h1_val = parse_numeric_or_nan(fields[h1_idx], h1_idx);
                let h2_val = parse_numeric_or_nan(fields[h2_idx], h2_idx);

                let max_allele = aggregation.aggregate(h1_val, h2_val);
                data_matrix[[sample_idx, new_idx]] = max_allele;
            }
        }
    }

    data_matrix
}

/// Parallel data loading for large feature sets using streaming approach
fn parallel_data_loading(
    combined: &std::path::Path,
    num_samples: usize,
    selected_indices: &[usize],
    aggregation: AlleleAggregation,
) -> Array2<f64> {
    let num_cores = rayon::current_num_threads();
    println!(
        "Using streaming parallel data loading ({} selected features, {} CPU cores)...",
        selected_indices.len(),
        num_cores
    );

    // Create mapping from original index to new index
    let selected_idx_map: std::collections::HashMap<usize, usize> = selected_indices
        .iter()
        .enumerate()
        .map(|(new_idx, &original_idx)| (original_idx, new_idx))
        .collect();

    // Initialize result matrix
    let mut data_matrix = Array2::<f64>::zeros((num_samples, selected_indices.len()));

    // Stream through file and process relevant lines in chunks
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();

    // Skip metadata lines and header
    let _header = crate::utils::skip_metadata_lines(&mut lines);

    let mut chunk_buffer = Vec::new();
    let chunk_size = 5_000; // Smaller chunks for data loading
    let mut processed_lines = 0;

    println!("Processing data loading in streaming chunks of {} lines...", chunk_size);

    for (region_idx, line_result) in lines.enumerate() {
        let line = line_result.unwrap_or_else(|e| {
            eprintln!("ERROR: Failed to read line {} during data loading: {}", region_idx + 2, e);
            std::process::exit(1);
        });
        // Only collect lines for selected indices
        if selected_idx_map.contains_key(&region_idx) {
            chunk_buffer.push((region_idx, line));
        }

        // Process chunk when full or when we have enough selected lines
        if chunk_buffer.len() >= chunk_size
            || (!chunk_buffer.is_empty() && chunk_buffer.len() >= selected_indices.len().min(1000))
        {
            // Process chunk in parallel
            let parsed_chunk: Vec<(usize, Vec<f64>)> = chunk_buffer
                .par_iter()
                .filter_map(|(region_idx, line)| {
                    if let Some(&new_idx) = selected_idx_map.get(region_idx) {
                        let fields: Vec<&str> = line.trim().split('\t').collect();

                        // Parse sample data for this region
                        let mut sample_data = Vec::with_capacity(num_samples);
                        for sample_idx in 0..num_samples {
                            let h1_idx = crate::filetype::STR_FIXED_COLUMNS + sample_idx * 2;
                            let h2_idx = crate::filetype::STR_FIXED_COLUMNS + 1 + sample_idx * 2;

                            let h1_val = parse_numeric_or_nan(fields[h1_idx], h1_idx);
                            let h2_val = parse_numeric_or_nan(fields[h2_idx], h2_idx);
                            let max_allele = aggregation.aggregate(h1_val, h2_val);
                            sample_data.push(max_allele);
                        }

                        Some((new_idx, sample_data))
                    } else {
                        None
                    }
                })
                .collect();

            // Fill matrix with parsed chunk data
            for (new_idx, sample_data) in parsed_chunk {
                for (sample_idx, &value) in sample_data.iter().enumerate() {
                    data_matrix[[sample_idx, new_idx]] = value;
                }
            }

            processed_lines += chunk_buffer.len();
            chunk_buffer.clear();

            if processed_lines % 10_000 == 0 {
                println!("Processed {} selected lines so far...", processed_lines);
            }
        }
    }

    // Process any remaining lines in the buffer
    if !chunk_buffer.is_empty() {
        let parsed_chunk: Vec<(usize, Vec<f64>)> = chunk_buffer
            .par_iter()
            .filter_map(|(region_idx, line)| {
                if let Some(&new_idx) = selected_idx_map.get(region_idx) {
                    let fields: Vec<&str> = line.trim().split('\t').collect();

                    let mut sample_data = Vec::with_capacity(num_samples);
                    for sample_idx in 0..num_samples {
                        let h1_idx = crate::filetype::STR_FIXED_COLUMNS + sample_idx * 2;
                        let h2_idx = crate::filetype::STR_FIXED_COLUMNS + 1 + sample_idx * 2;

                        let h1_val = parse_numeric_or_nan(fields[h1_idx], h1_idx);
                        let h2_val = parse_numeric_or_nan(fields[h2_idx], h2_idx);
                        let max_allele = aggregation.aggregate(h1_val, h2_val);
                        sample_data.push(max_allele);
                    }

                    Some((new_idx, sample_data))
                } else {
                    None
                }
            })
            .collect();

        for (new_idx, sample_data) in parsed_chunk {
            for (sample_idx, &value) in sample_data.iter().enumerate() {
                data_matrix[[sample_idx, new_idx]] = value;
            }
        }
    }

    data_matrix
}

/// Quick estimate of feature count without loading entire file
fn estimate_feature_count(combined: &std::path::Path) -> usize {
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();

    // Skip metadata lines and header
    let _header = crate::utils::skip_metadata_lines(&mut lines);

    // Count first few hundred lines to estimate total
    let sample_size = 500;
    let mut count = 0;

    for line_result in lines.take(sample_size) {
        let _line = line_result.unwrap_or_else(|e| {
            eprintln!("ERROR: Failed to read line while estimating feature count: {}", e);
            std::process::exit(1);
        });
        count += 1;
    }

    if count < sample_size {
        // File is smaller than sample size, return exact count
        count
    } else {
        // Estimate based on file size (rough approximation)
        // This could be improved with more sophisticated estimation
        let file_size = std::fs::metadata(combined).map(|m| m.len()).unwrap_or(0);
        let estimated_lines = (file_size / 100).max(count as u64) as usize; // Rough estimate
        estimated_lines.min(10_000_000) // Cap at reasonable maximum
    }
}

/// Proper PCA implementation using eigenvalue decomposition
/// Returns principal components and explained variance ratios
fn perform_proper_pca(
    data: &Array2<f64>,
    n_components: usize,
) -> (Array2<f64>, Array1<f64>, Array2<f64>) {
    let n_samples = data.nrows();
    let n_features = data.ncols();

    if n_samples == 0 || n_features == 0 {
        eprintln!(
            "ERROR: No features remained after filtering (0 out of {} samples × {} features).",
            n_samples, n_features
        );
        eprintln!("This typically means all loci were excluded due to high missing data rates.");
        eprintln!("\nPossible causes:");
        eprintln!("  - The input file has too many missing/NaN values");
        eprintln!("  - The data lacks sufficient variation across samples");
        eprintln!("\nSuggestions:");
        eprintln!("  - Check the input file for data quality issues");
        eprintln!("  - Ensure the combined file was generated correctly with 'inquiSTR combine'");
        std::process::exit(1);
    }

    let n_components = n_components.min(n_samples.min(n_features));

    // Step 1: Standardize the data (center and scale each feature to unit variance)
    // This ensures all features contribute equally regardless of their absolute scale,
    // which is critical for STR data where allele lengths vary enormously between loci.
    let mean = data.mean_axis(Axis(0)).unwrap();
    let centered_data = data - &mean.insert_axis(Axis(0));
    let std_dev = centered_data.map_axis(Axis(0), |col| {
        let variance = col.iter().map(|x| x * x).sum::<f64>() / (n_samples as f64 - 1.0);
        variance.sqrt()
    });
    // Avoid division by zero for constant features
    let std_dev = std_dev.mapv(|s| if s == 0.0 { 1.0 } else { s });
    let standardized_data = &centered_data / &std_dev.insert_axis(Axis(0));

    // Step 2: Compute correlation matrix (covariance of standardized data)
    println!("Computing correlation matrix ({} x {} features)...", n_features, n_features);
    let covariance = standardized_data.t().dot(&standardized_data) / ((n_samples - 1) as f64);

    // Total variance = trace of the covariance matrix (sum of ALL eigenvalues)
    let total_variance = covariance.diag().sum();

    // Step 3: Eigenvalue decomposition using nalgebra's SymmetricEigen
    let n = covariance.nrows();
    let nalgebra_matrix = nalgebra::DMatrix::from_iterator(n, n, covariance.iter().copied());
    let eigen = nalgebra::SymmetricEigen::new(nalgebra_matrix);

    // Sort eigenvalues/eigenvectors by descending eigenvalue
    let mut eigen_pairs: Vec<(f64, Vec<f64>)> = eigen
        .eigenvalues
        .iter()
        .enumerate()
        .map(|(i, &val)| {
            let vec: Vec<f64> = eigen.eigenvectors.column(i).iter().copied().collect();
            (val, vec)
        })
        .collect();
    eigen_pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

    let eigenvalues = Array1::from_vec(
        eigen_pairs
            .iter()
            .take(n_components)
            .map(|(v, _)| *v)
            .collect(),
    );
    let mut eigenvectors = Array2::zeros((n, n_components));
    for (i, (_, vec)) in eigen_pairs.iter().take(n_components).enumerate() {
        for (j, &val) in vec.iter().enumerate() {
            eigenvectors[[j, i]] = val;
        }
    }

    // Step 4: Project data onto principal components
    // PC_scores = X_standardized * eigenvectors
    let pca_data = standardized_data.dot(&eigenvectors);

    // Step 5: Calculate explained variance ratios using total variance from trace
    let explained_variance = eigenvalues.mapv(|x| (x / total_variance) * 100.0);

    (pca_data, explained_variance, eigenvectors)
}

/// Write PC scores to a tab-separated file for use as covariates
fn write_pc_scores(
    pca_data: &Array2<f64>,
    sample_names: &[String],
    n_components: usize,
    output_path: &std::path::Path,
) {
    use std::io::Write;

    let mut file = match std::fs::File::create(output_path) {
        Ok(f) => std::io::BufWriter::new(f),
        Err(e) => {
            eprintln!(
                "ERROR: Failed to create scores output file {}: {}",
                output_path.display(),
                e
            );
            std::process::exit(1);
        }
    };

    // Write header
    let mut header = vec!["sample".to_string()];
    for i in 1..=n_components {
        header.push(format!("PC{}", i));
    }
    if let Err(e) = writeln!(file, "{}", header.join("\t")) {
        eprintln!("ERROR: Failed to write to scores output file: {}", e);
        std::process::exit(1);
    }

    // Write PC scores for each sample
    for (sample_idx, sample_name) in sample_names.iter().enumerate() {
        let mut row = vec![sample_name.clone()];
        for pc_idx in 0..n_components {
            let score = pca_data[[sample_idx, pc_idx]];
            row.push(format!("{:.6}", score));
        }
        if let Err(e) = writeln!(file, "{}", row.join("\t")) {
            eprintln!("ERROR: Failed to write to scores output file: {}", e);
            std::process::exit(1);
        }
    }

    println!(
        "PC scores saved to: {} ({} samples × {} PCs)",
        output_path.display(),
        sample_names.len(),
        n_components
    );
}

/// Create a PCA plot using kuva
fn create_pca_plot(
    pca_data: &Array2<f64>,
    sample_names: &[String],
    explained_variance: &Array1<f64>,
    output: &str,
    title: &str,
) {
    // Extract PC1 and PC2 coordinates and zip into (x, y) pairs
    let data: Vec<(f64, f64)> = pca_data
        .column(0)
        .iter()
        .copied()
        .zip(pca_data.column(1).iter().copied())
        .collect();

    // Create scatter plot with sample names shown on hover
    let trace = ScatterPlot::new()
        .with_data(data)
        .with_color("steelblue")
        .with_size(5.0)
        .with_marker_opacity(0.8)
        .with_tooltip_labels(sample_names.iter().cloned())
        .with_tooltips();

    let plots: Vec<Plot> = vec![trace.into()];

    // Set up layout
    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_x_label(format!("PC1 ({:.1}% variance)", explained_variance[0]))
        .with_y_label(format!("PC2 ({:.1}% variance)", explained_variance[1]))
        .with_interactive();

    // Save to SVG file
    let svg = render_to_svg(plots, layout);
    std::fs::write(output, svg).expect("Failed to write PCA plot");
    println!("PCA plot saved to: {}", output);

    // Print summary statistics
    println!("PCA Summary:");
    println!("  PC1 explains {:.1}% of variance", explained_variance[0]);
    println!("  PC2 explains {:.1}% of variance", explained_variance[1]);
    println!(
        "  Total explained by PC1+PC2: {:.1}%",
        explained_variance[0] + explained_variance[1]
    );
}

/// Main PCA function with intelligent feature selection for large datasets
pub fn pca(
    combined: PathBuf,
    output: String,
    _n_components: usize,
    threads: usize,
    aggregation: AlleleAggregation,
    scores_output: Option<PathBuf>,
) {
    if !combined.exists() {
        eprintln!("ERROR: Combined file does not exist: {}", combined.display());
        std::process::exit(1);
    }

    // Validate that input is a combined file, not individual, and detect data type
    let file_type = crate::filetype::read_file_type_metadata(&combined);
    let is_kmer_file = matches!(file_type, Some(crate::filetype::FileType::CombinedKmer));

    if let Some(ref ftype) = file_type
        && !matches!(
            ftype,
            crate::filetype::FileType::CombinedCall | crate::filetype::FileType::CombinedKmer
        )
    {
        eprintln!("ERROR: PCA requires a combined file (combined_call or combined_kmer).");
        eprintln!("The provided file appears to be: {:?}", ftype);
        eprintln!("\nPlease use 'inquiSTR combine' to merge individual sample files first.");
        std::process::exit(1);
    }

    let data_type = if is_kmer_file {
        "kmer frequencies"
    } else {
        "STR genotypes"
    };

    // Configure thread pool based on user input
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("Failed to build thread pool");
        println!("Configured thread pool with {} threads", threads);
    } else {
        println!("Using automatic thread detection ({} threads)", rayon::current_num_threads());
    }

    println!("Reading combined inquiSTR file: {}", combined.display());
    println!("Detected file type: {}", data_type);
    if !is_kmer_file {
        println!("Using '{}' aggregation method for H1/H2 allele lengths", aggregation);
    }

    // Pre-determine if we need feature selection based on a quick file scan
    let estimated_features = estimate_feature_count(combined.as_path());
    let feature_name = if is_kmer_file { "kmers" } else { "STR regions" };
    println!("Estimated {} {} in file", estimated_features, feature_name);

    let max_features = if estimated_features > 100_000 {
        // For massive datasets (>100k loci), be very aggressive
        let suggested_features = 3000;
        println!(
            "Massive dataset detected ({} loci). Using aggressive feature selection to {} most informative loci...",
            estimated_features, suggested_features
        );
        Some(suggested_features)
    } else if estimated_features > 10_000 {
        // For very large datasets (>10k loci), use intelligent selection
        let suggested_features = 5000;
        println!(
            "Large dataset detected ({} loci). Using memory-efficient parsing with feature selection to {} loci...",
            estimated_features, suggested_features
        );
        Some(suggested_features)
    } else if estimated_features > 1000 {
        // For medium datasets, use moderate selection
        let suggested_features = 2000.max(estimated_features / 2);
        println!(
            "Medium dataset detected ({} loci). Using feature selection to {} loci...",
            estimated_features, suggested_features
        );
        Some(suggested_features)
    } else {
        // For small datasets, use all features
        println!("Small dataset detected ({} loci). Loading all features...", estimated_features);
        None
    };

    // Use optimized parsing with or without feature selection
    let (data_matrix, sample_names) =
        parse_combined_file_with_selection(combined.as_path(), max_features, aggregation);

    println!("Loaded data:");
    println!("  {} samples", sample_names.len());
    if is_kmer_file {
        println!("  {} kmers", data_matrix.ncols());
    } else {
        println!("  {} STR regions", data_matrix.ncols());
    }

    if sample_names.len() < 2 {
        eprintln!("ERROR: Need at least 2 samples for PCA, but only {} found.", sample_names.len());
        eprintln!("PCA requires data from at least 2 samples.");
        std::process::exit(1);
    }

    // Show information about feature selection (if it occurred)
    if max_features.is_some() {
        let feature_name = if is_kmer_file { "kmers" } else { "regions" };
        println!("Selected {} include features with:", feature_name);
        println!("  - High variance across samples");
        println!("  - Low missing data rate (<50%)");
        println!("  - Good dynamic range");

        // Feature names not collected to save memory
    }

    let n_components = 2; // For now, stick to 2D visualization

    println!("Performing proper PCA with {} components...", n_components);
    let (pca_data, explained_variance, _eigenvectors) =
        perform_proper_pca(&data_matrix, n_components);

    println!("PCA Results:");
    for i in 0..n_components {
        println!("  PC{}: {:.1}% variance explained", i + 1, explained_variance[i]);
    }
    println!("  Total explained by PC1+PC2: {:.1}%", explained_variance.iter().sum::<f64>());

    // Memory usage summary
    let matrix_size_mb = (data_matrix.len() * std::mem::size_of::<f64>()) as f64 / 1_048_576.0;
    let num_cores = rayon::current_num_threads();
    println!("  Memory usage: {:.1} MB for data matrix", matrix_size_mb);
    println!("  Parallel processing: {} CPU cores utilized", num_cores);

    // Create and save the PCA plot
    let plot_title = if is_kmer_file {
        "PCA of Kmer Frequencies"
    } else {
        "PCA of STR Genotypes"
    };
    create_pca_plot(&pca_data, &sample_names, &explained_variance, &output, plot_title);

    // Write PC scores if requested
    if let Some(scores_path) = scores_output {
        write_pc_scores(&pca_data, &sample_names, n_components, &scores_path);
    }

    println!("PCA analysis complete! Plot saved to: {}", output);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_pca_with_test_data() {
        // This test requires the combined test data file with proper header to exist
        let test_file = PathBuf::from("test_combined_with_header.tsv");
        if test_file.exists() {
            // This should not panic and should create the SVG file
            pca(
                test_file,
                "test_pca_output.svg".to_string(),
                10,
                0,
                AlleleAggregation::Max,
                None,
            );

            // Verify the output file was created
            assert!(PathBuf::from("test_pca_output.svg").exists());

            // Clean up
            std::fs::remove_file("test_pca_output.svg").ok();
        }
    }
}
