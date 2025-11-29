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

use ndarray::{Array1, Array2, Axis};
use plotly::color::{NamedColor, Rgba};
use plotly::common::{Marker, Mode, Title};
use plotly::layout::{Axis as PlotAxis, Layout};
use plotly::{Plot, Scatter};
use rayon::prelude::*;
use std::io::BufRead;
use std::path::PathBuf;
// Removed unused imports for cleaner streaming implementation

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
        panic!(
            "Invalid combined kmer file header. Expected format: kmer\\tsample1\\tsample2\\t..."
        );
    }

    // Extract sample names (all columns after "kmer")
    let sample_names: Vec<String> = header_fields[1..].iter().map(|s| s.to_string()).collect();
    let num_samples = sample_names.len();

    if num_samples < 2 {
        panic!("Need at least 2 samples for PCA, found {}", num_samples);
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
                "Warning: Skipping malformed line {} (expected {} columns, got {})",
                line_num + 2,
                expected_cols,
                fields.len()
            );
            continue;
        }

        // Parse kmer frequencies for this kmer across all samples
        let mut row_data = Vec::with_capacity(num_samples);
        for sample_idx in 0..num_samples {
            let freq_idx = 1 + sample_idx;
            let freq: f64 = fields[freq_idx].parse().unwrap_or(0.0);
            row_data.push(freq);
        }
        data_rows.push(row_data);
    }

    if data_rows.is_empty() {
        panic!("No data lines found after header in kmer file");
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
        let line = match line_result {
            Ok(l) => l,
            Err(_) => continue,
        };

        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 1 + num_samples {
            continue;
        }

        // Calculate variance for this kmer
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        let mut count = 0;

        for sample_idx in 0..num_samples {
            let freq: f64 = fields[1 + sample_idx].parse().unwrap_or(0.0);
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
            let line = match line_result {
                Ok(l) => l,
                Err(_) => continue,
            };

            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() != 1 + num_samples {
                continue;
            }

            for sample_idx in 0..num_samples {
                let freq: f64 = fields[1 + sample_idx].parse().unwrap_or(0.0);
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
                eprintln!("  - STR file: chromosome\\tbegin\\tend\\tsample1_H1\\tsample1_H2...");
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

    if header_fields.len() < 5
        || header_fields[0] != "chromosome"
        || header_fields[1] != "begin"
        || header_fields[2] != "end"
    {
        eprintln!("Error: Invalid STR file header format.");
        eprintln!("Expected: chromosome\\tbegin\\tend\\tsample1_H1\\tsample1_H2\\t...");
        eprintln!("Got: {}", header_line);
        eprintln!("\nThis file does not appear to be a valid combined STR file.");
        eprintln!("If this is a kmer file, it should have been auto-detected.");
        std::process::exit(1);
    }

    // Extract sample names from header (columns 3+ should be sample_H1, sample_H2 pattern)
    let data_cols = &header_fields[3..];
    if !data_cols.len().is_multiple_of(2) {
        panic!("Invalid header: number of sample columns must be even (H1/H2 pairs)");
    }

    let num_samples = data_cols.len() / 2;
    let mut sample_names = Vec::new();

    for i in 0..num_samples {
        let h1_col = &data_cols[i * 2];
        let h2_col = &data_cols[i * 2 + 1];

        // Extract sample name from H1 column (should end with _H1)
        if !h1_col.ends_with("_H1") {
            panic!(
                "Invalid header: expected column {} to end with '_H1', got '{}'",
                3 + i * 2,
                h1_col
            );
        }
        if !h2_col.ends_with("_H2") {
            panic!(
                "Invalid header: expected column {} to end with '_H2', got '{}'",
                3 + i * 2 + 1,
                h2_col
            );
        }

        let sample_name_h1 = h1_col.trim_end_matches("_H1");
        let sample_name_h2 = h2_col.trim_end_matches("_H2");

        if sample_name_h1 != sample_name_h2 {
            panic!(
                "Header error: H1 and H2 columns have different sample names: '{}' vs '{}'",
                sample_name_h1, sample_name_h2
            );
        }

        sample_names.push(sample_name_h1.to_string());
    }

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

        let expected_cols = 3 + num_samples * 2;
        if fields.len() != expected_cols {
            panic!(
                "Malformed line {} (expected {} columns, got {}): {}",
                line_num + 2,
                expected_cols,
                fields.len(),
                line
            );
        }

        // Skip region name creation for memory efficiency

        // Parse STR lengths for this region across all samples
        let mut row_data = Vec::with_capacity(num_samples);
        for sample_idx in 0..num_samples {
            let h1_idx = 3 + sample_idx * 2;
            let h2_idx = 4 + sample_idx * 2;

            let h1_val: f64 = fields[h1_idx]
                .parse()
                .map_err(|e| {
                    format!(
                        "Invalid H1 value '{}' at line {}, column {}: {}",
                        fields[h1_idx],
                        line_num + 2,
                        h1_idx + 1,
                        e
                    )
                })
                .unwrap_or_else(|e| {
                    // Handle NaN values specifically
                    if fields[h1_idx].eq_ignore_ascii_case("nan") {
                        0.0
                    } else {
                        panic!("{}", e);
                    }
                });

            let h2_val: f64 = fields[h2_idx]
                .parse()
                .map_err(|e| {
                    format!(
                        "Invalid H2 value '{}' at line {}, column {}: {}",
                        fields[h2_idx],
                        line_num + 2,
                        h2_idx + 1,
                        e
                    )
                })
                .unwrap_or_else(|e| {
                    // Handle NaN values specifically
                    if fields[h2_idx].eq_ignore_ascii_case("nan") {
                        0.0
                    } else {
                        panic!("{}", e);
                    }
                });

            // Apply selected aggregation method for H1/H2 allele lengths
            let aggregated_value = aggregation.aggregate(h1_val, h2_val);
            row_data.push(aggregated_value);
        }
        data_rows.push(row_data);
    }

    if data_rows.is_empty() {
        panic!("No data lines found after header");
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
        let expected_cols = 3 + num_samples * 2;
        if fields.len() != expected_cols {
            continue; // Skip malformed lines
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
            let h1_idx = 3 + sample_idx * 2;
            let h2_idx = 4 + sample_idx * 2;

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
            let h1_idx = 3 + sample_idx * 2;
            let h2_idx = 4 + sample_idx * 2;

            // Optimized parsing with early NaN check
            let h1_val: f64 =
                if fields[h1_idx].eq_ignore_ascii_case("nan") || fields[h1_idx].is_empty() {
                    0.0
                } else {
                    fields[h1_idx].parse().unwrap_or(0.0)
                };

            let h2_val: f64 =
                if fields[h2_idx].eq_ignore_ascii_case("nan") || fields[h2_idx].is_empty() {
                    0.0
                } else {
                    fields[h2_idx].parse().unwrap_or(0.0)
                };

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
        if let Ok(line) = line_result {
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
    let expected_cols = 3 + num_samples * 2;
    if fields.len() != expected_cols {
        return None; // Skip malformed lines
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
        let h1_idx = 3 + sample_idx * 2;
        let h2_idx = 4 + sample_idx * 2;

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
        let h1_idx = 3 + sample_idx * 2;
        let h2_idx = 4 + sample_idx * 2;

        // Optimized parsing with early NaN check
        let h1_val: f64 = if fields[h1_idx].eq_ignore_ascii_case("nan") || fields[h1_idx].is_empty()
        {
            0.0
        } else {
            fields[h1_idx].parse().unwrap_or(0.0)
        };

        let h2_val: f64 = if fields[h2_idx].eq_ignore_ascii_case("nan") || fields[h2_idx].is_empty()
        {
            0.0
        } else {
            fields[h2_idx].parse().unwrap_or(0.0)
        };

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
                let h1_idx = 3 + sample_idx * 2;
                let h2_idx = 4 + sample_idx * 2;

                let h1_val: f64 = fields[h1_idx].parse().unwrap_or(0.0);
                let h2_val: f64 = fields[h2_idx].parse().unwrap_or(0.0);

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
        if let Ok(line) = line_result {
            // Only collect lines for selected indices
            if selected_idx_map.contains_key(&region_idx) {
                chunk_buffer.push((region_idx, line));
            }

            // Process chunk when full or when we have enough selected lines
            if chunk_buffer.len() >= chunk_size
                || (!chunk_buffer.is_empty()
                    && chunk_buffer.len() >= selected_indices.len().min(1000))
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
                                let h1_idx = 3 + sample_idx * 2;
                                let h2_idx = 4 + sample_idx * 2;

                                let h1_val: f64 = fields[h1_idx].parse().unwrap_or(0.0);
                                let h2_val: f64 = fields[h2_idx].parse().unwrap_or(0.0);
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
                        let h1_idx = 3 + sample_idx * 2;
                        let h2_idx = 4 + sample_idx * 2;

                        let h1_val: f64 = fields[h1_idx].parse().unwrap_or(0.0);
                        let h2_val: f64 = fields[h2_idx].parse().unwrap_or(0.0);
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
        if line_result.is_ok() {
            count += 1;
        }
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
        panic!("Empty data matrix provided to PCA");
    }

    let n_components = n_components.min(n_samples.min(n_features));

    // Step 1: Center the data (subtract mean for each feature)
    let mean = data.mean_axis(Axis(0)).unwrap();
    let centered_data = data - &mean.insert_axis(Axis(0));

    // Step 2: Compute covariance matrix
    // Cov = (X^T * X) / (n - 1)
    println!("Computing covariance matrix ({} x {} features)...", n_features, n_features);
    let covariance = if n_features > 1000 {
        // For large matrices, the matrix multiplication can benefit from parallel BLAS
        // ndarray uses OpenBLAS or similar which should parallelize automatically
        centered_data.t().dot(&centered_data) / ((n_samples - 1) as f64)
    } else {
        centered_data.t().dot(&centered_data) / ((n_samples - 1) as f64)
    };

    // Step 3: Eigenvalue decomposition using power iteration for efficiency
    // For large matrices, we'll use a simplified approach focusing on top components
    let (eigenvalues, eigenvectors) = compute_top_eigenvectors(&covariance, n_components);

    // Step 4: Project data onto principal components
    // PC_scores = X_centered * eigenvectors
    let pca_data = centered_data.dot(&eigenvectors);

    // Step 5: Calculate explained variance ratios
    let total_variance = eigenvalues.sum();
    let explained_variance = eigenvalues.mapv(|x| (x / total_variance) * 100.0);

    (pca_data, explained_variance, eigenvectors)
}

/// Compute top k eigenvectors using power iteration method
/// This is more efficient than full eigendecomposition for large matrices when we only need top components
fn compute_top_eigenvectors(matrix: &Array2<f64>, k: usize) -> (Array1<f64>, Array2<f64>) {
    let n = matrix.nrows();
    let k = k.min(n);

    let mut eigenvalues = Array1::zeros(k);
    let mut eigenvectors = Array2::zeros((n, k));

    // For small matrices, use a simplified approach
    if n <= 1000 {
        // Use nalgebra for proper eigendecomposition if available, otherwise simplified approach
        return simplified_eigen_decomposition(matrix, k);
    }

    // Power iteration for large matrices
    let mut remaining_matrix = matrix.clone();

    for i in 0..k {
        let (eigenval, eigenvec) = power_iteration(&remaining_matrix, 100, 1e-6);

        eigenvalues[i] = eigenval;
        eigenvectors.column_mut(i).assign(&eigenvec);

        // Deflate matrix: remove the found eigenvector's contribution
        let eigenvec_col = eigenvec.clone().insert_axis(Axis(1));
        let eigenvec_row = eigenvec.insert_axis(Axis(0));
        let outer_product = eigenvec_col.dot(&eigenvec_row);
        remaining_matrix = &remaining_matrix - &(outer_product * eigenval);
    }

    (eigenvalues, eigenvectors)
}

/// Simplified eigendecomposition for smaller matrices
fn simplified_eigen_decomposition(matrix: &Array2<f64>, k: usize) -> (Array1<f64>, Array2<f64>) {
    let n = matrix.nrows();
    let k = k.min(n);

    // For demonstration, use power iteration even for small matrices
    // In production, you'd want to use a proper linear algebra library
    let mut eigenvalues = Array1::zeros(k);
    let mut eigenvectors = Array2::zeros((n, k));

    let mut working_matrix = matrix.clone();

    for i in 0..k {
        let (eigenval, eigenvec) = power_iteration(&working_matrix, 50, 1e-8);
        eigenvalues[i] = eigenval;
        eigenvectors.column_mut(i).assign(&eigenvec);

        // Deflation: remove this component
        let eigenvec_col = eigenvec.clone().insert_axis(Axis(1));
        let eigenvec_row = eigenvec.insert_axis(Axis(0));
        let rank_one = eigenvec_col.dot(&eigenvec_row) * eigenval;
        working_matrix = &working_matrix - &rank_one;
    }

    (eigenvalues, eigenvectors)
}

/// Power iteration algorithm to find dominant eigenvector
fn power_iteration(
    matrix: &Array2<f64>,
    max_iterations: usize,
    tolerance: f64,
) -> (f64, Array1<f64>) {
    let n = matrix.nrows();

    // Initialize with random vector
    let mut vector = Array1::from_elem(n, 1.0 / (n as f64).sqrt());
    let mut eigenvalue = 0.0;

    for _iteration in 0..max_iterations {
        // v_new = A * v
        let new_vector = matrix.dot(&vector);

        // Calculate eigenvalue (Rayleigh quotient)
        let new_eigenvalue = vector.dot(&new_vector);

        // Normalize vector
        let norm = new_vector.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 0.0 {
            vector = new_vector / norm;
        }

        // Check convergence
        if (new_eigenvalue - eigenvalue).abs() < tolerance {
            eigenvalue = new_eigenvalue;
            break;
        }

        eigenvalue = new_eigenvalue;
    }

    (eigenvalue, vector)
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

/// Create a PCA plot using Plotly
fn create_pca_plot(
    pca_data: &Array2<f64>,
    sample_names: &[String],
    explained_variance: &Array1<f64>,
    output: &str,
    title: &str,
) {
    // Extract PC1 and PC2 coordinates
    let pc1_coords: Vec<f64> = pca_data.column(0).to_vec();
    let pc2_coords: Vec<f64> = pca_data.column(1).to_vec();

    // Create scatter plot
    let trace = Scatter::new(pc1_coords, pc2_coords)
        .mode(Mode::Markers)
        .name("Samples")
        .text_array(sample_names.to_vec())
        .marker(
            Marker::new()
                .size(10)
                .color(Rgba::new(31, 119, 180, 1.0))
                .line(
                    plotly::common::Line::new()
                        .color(NamedColor::White)
                        .width(1.0),
                ),
        );

    let mut plot = Plot::new();
    plot.add_trace(trace);

    // Set up layout
    let layout = Layout::new()
        .title(Title::with_text(title))
        .x_axis(
            PlotAxis::new()
                .title(Title::with_text(format!("PC1 ({:.1}% variance)", explained_variance[0]))),
        )
        .y_axis(
            PlotAxis::new()
                .title(Title::with_text(format!("PC2 ({:.1}% variance)", explained_variance[1]))),
        )
        .show_legend(false);

    plot.set_layout(layout);

    // Save to HTML file
    plot.write_html(output);
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
        panic!("Combined file does not exist: {}", combined.display());
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
        panic!("Need at least 2 samples for PCA, but only {} found", sample_names.len());
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
            // This should not panic and should create the HTML file
            pca(
                test_file,
                "test_pca_output.html".to_string(),
                10,
                0,
                AlleleAggregation::Max,
                None,
            );

            // Verify the output file was created
            assert!(PathBuf::from("test_pca_output.html").exists());

            // Clean up
            std::fs::remove_file("test_pca_output.html").ok();
        }
    }
}
