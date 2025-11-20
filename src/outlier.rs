use clap::ValueEnum;

use dbscan::Classification::*;
use dbscan::Model;
use log::debug;

use std::cmp::max;
use std::collections::HashSet;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use rayon::prelude::*;

/// Helper function to clean sample names by removing _H1/_H2 suffixes
#[inline]
fn clean_sample_name(sample_name: &str) -> &str {
    if sample_name.ends_with("_H1") || sample_name.ends_with("_H2") {
        &sample_name[..sample_name.len() - 3]
    } else {
        sample_name
    }
}

/// Parse sample input - can be a file path, comma-separated names, or a single name
pub fn parse_sample_input(input: &str) -> Vec<String> {
    let path = Path::new(input);

    // Check if input is a file path
    if path.exists() && path.is_file() {
        eprintln!("Reading sample names from file: {}", input);
        let file = crate::utils::reader(input);
        file.lines()
            .map(|line| line.unwrap().trim().to_string())
            .filter(|line| !line.is_empty())
            .collect()
    } else if input.contains(',') {
        // Comma-separated sample names
        eprintln!("Parsing comma-separated sample names");
        input
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    } else {
        // Single sample name
        eprintln!("Using single sample name: {}", input);
        vec![input.to_string()]
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Method {
    Zscore,
    Dbscan,
}

/// Streaming statistics calculation that avoids intermediate allocations
fn streaming_stats(values: &[f32]) -> (f32, f32) {
    let mut sum = 0.0;
    let mut sum_sq = 0.0;
    let mut count = 0;

    for &value in values {
        if !value.is_nan() && value > 0.0 {
            sum += value;
            sum_sq += value * value;
            count += 1;
        }
    }

    if count == 0 {
        return (0.0, 0.0);
    }

    let count_f = count as f32;
    let mean = sum / count_f;
    let variance = (sum_sq / count_f) - (mean * mean);

    (mean, variance.max(0.0).sqrt())
}

/// Get sample names from the header of the combined file
fn get_sample_names_from_file(file_path: &Path, is_kmer: bool) -> Vec<String> {
    let file_reader = crate::utils::reader(&file_path.to_string_lossy());
    let mut lines = file_reader.lines();

    // Skip metadata lines if present
    let header_line = crate::utils::skip_metadata_lines(&mut lines);
    let skip_columns = if is_kmer { 1 } else { 3 }; // kmer files skip 1, STR files skip 3
    header_line
        .split('\t')
        .skip(skip_columns)
        .map(|s| s.to_string())
        .collect()
}

/// Validate that requested samples exist in the combined file
fn validate_samples(file_path: &Path, requested_samples: &[String], is_kmer: bool) {
    let available_samples = get_sample_names_from_file(file_path, is_kmer);

    // Create a set of available sample names (cleaned) for efficient lookup
    let available_set: HashSet<String> = available_samples
        .iter()
        .map(|s| clean_sample_name(s).to_string())
        .collect();

    // Check each requested sample
    let mut missing_samples = Vec::new();
    for sample in requested_samples {
        if !available_set.contains(sample) {
            missing_samples.push(sample.clone());
        }
    }

    if !missing_samples.is_empty() {
        eprintln!("ERROR: The following requested samples were not found in the combined file:");
        for sample in &missing_samples {
            eprintln!("  - {}", sample);
        }
        eprintln!("\nAvailable samples in file ({} total):", available_samples.len());
        for (i, sample) in available_samples.iter().take(10).enumerate() {
            eprintln!("  - {}", clean_sample_name(sample));
            if i == 9 && available_samples.len() > 10 {
                eprintln!("  ... and {} more", available_samples.len() - 10);
                break;
            }
        }
        std::process::exit(1);
    }

    eprintln!("✓ All {} requested sample(s) found in combined file", requested_samples.len());
}

/// Detect if input file contains kmer frequency data or STR call data
/// Returns true if it's a kmer file, false if it's an STR call file
fn is_kmer_file(file_path: &Path) -> bool {
    let file_reader = crate::utils::reader(&file_path.to_string_lossy());
    let mut lines = file_reader.lines();

    // Skip metadata lines if present
    let first_line = crate::utils::skip_metadata_lines(&mut lines);
    let fields: Vec<&str> = first_line.split('\t').collect();

    // Check if it's a kmer file format
    // Kmer files have "kmer" as first column header
    if fields.len() >= 2 && fields[0] == "kmer" {
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

    eprintln!("ERROR: Unable to determine file format for: {}", file_path.display());
    eprintln!("File must be either STR call data or kmer frequency data from inquiSTR combine.");
    std::process::exit(1);
}

pub fn outlier(
    combined: PathBuf,
    minsize: u32,
    zscore_cutoff: f32,
    method: Method,
    subset: Option<Vec<String>>,
    threads: usize,
) {
    // Configure thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .expect("Failed to build thread pool");

    // Validate that input is a combined file, not individual
    if let Some(file_type) = crate::combine::read_file_type_metadata(&combined) {
        if !matches!(
            file_type,
            crate::combine::FileType::CombinedCall | crate::combine::FileType::CombinedKmer
        ) {
            eprintln!("ERROR: Outlier detection requires a combined file (combined_call or combined_kmer).");
            eprintln!("The provided file appears to be: {:?}", file_type);
            eprintln!("\nPlease use 'inquiSTR combine' to merge individual sample files first.");
            std::process::exit(1);
        }
    }

    // Detect file format
    let is_kmer_format = is_kmer_file(&combined);

    // Validate samples if subset is provided
    if let Some(ref samples) = subset {
        validate_samples(&combined, samples, is_kmer_format);
    }

    if is_kmer_format {
        eprintln!("Detected kmer frequency file format");
        outlier_kmer_analysis(combined, minsize, zscore_cutoff, method, subset);
    } else {
        eprintln!("Detected STR call file format");
        outlier_str_analysis(combined, minsize, zscore_cutoff, method, subset);
    }
}

/// Outlier analysis for STR call data (original functionality)
fn outlier_str_analysis(
    combined: PathBuf,
    minsize: u32,
    zscore_cutoff: f32,
    method: Method,
    subset: Option<Vec<String>>,
) {
    let file = crate::utils::reader(&combined.into_os_string().into_string().unwrap());
    let mut lines = file.lines();
    // Skip metadata lines if present
    let header_line = crate::utils::skip_metadata_lines(&mut lines);
    println!("chrom\tbegin\tend\toutliers");

    // Parse sample names once and store them
    let sample_names: Vec<String> = header_line
        .split('\t')
        .skip(3)
        .map(|s| s.to_string())
        .collect();

    let num_samples = sample_names.len();
    let mincluster = num_samples.ilog2() as usize;

    // Process lines in chunks for better memory efficiency
    let chunk_size = 1000; // Process 1000 lines at a time
    let mut line_buffer = Vec::with_capacity(chunk_size);
    let mut processed_count = 0;

    for line_result in lines {
        let line = line_result.unwrap();
        line_buffer.push(line);

        // When chunk is full, process it in parallel
        if line_buffer.len() == chunk_size {
            process_chunk(
                &line_buffer,
                &sample_names,
                minsize,
                zscore_cutoff,
                method,
                mincluster,
                &subset,
            );

            processed_count += line_buffer.len();
            if processed_count % 10_000 == 0 {
                eprintln!("Processed {} loci...", processed_count);
            }

            line_buffer.clear();
        }
    }

    // Process remaining lines
    if !line_buffer.is_empty() {
        process_chunk(
            &line_buffer,
            &sample_names,
            minsize,
            zscore_cutoff,
            method,
            mincluster,
            &subset,
        );
        processed_count += line_buffer.len();
    }

    eprintln!("Completed processing {} loci", processed_count);
}

/// Outlier analysis for kmer frequency data
fn outlier_kmer_analysis(
    combined: PathBuf,
    minsize: u32,
    zscore_cutoff: f32,
    method: Method,
    subset: Option<Vec<String>>,
) {
    let file = crate::utils::reader(&combined.into_os_string().into_string().unwrap());
    let mut lines = file.lines();
    // Skip metadata lines if present
    let header_line = crate::utils::skip_metadata_lines(&mut lines);
    println!("kmer\toutliers");

    // Parse sample names once and store them
    let sample_names: Vec<String> = header_line
        .split('\t')
        .skip(1) // Skip "kmer" column
        .map(|s| s.to_string())
        .collect();

    let num_samples = sample_names.len();
    let mincluster = num_samples.ilog2() as usize;

    // Process lines in chunks for better memory efficiency
    let chunk_size = 1000; // Process 1000 lines at a time
    let mut line_buffer = Vec::with_capacity(chunk_size);
    let mut processed_count = 0;

    for line_result in lines {
        let line = line_result.unwrap();
        line_buffer.push(line);

        // When chunk is full, process it in parallel
        if line_buffer.len() == chunk_size {
            process_kmer_chunk(
                &line_buffer,
                &sample_names,
                minsize,
                zscore_cutoff,
                method,
                mincluster,
                &subset,
            );

            processed_count += line_buffer.len();
            if processed_count % 10_000 == 0 {
                eprintln!("Processed {} kmers...", processed_count);
            }

            line_buffer.clear();
        }
    }

    // Process remaining lines
    if !line_buffer.is_empty() {
        process_kmer_chunk(
            &line_buffer,
            &sample_names,
            minsize,
            zscore_cutoff,
            method,
            mincluster,
            &subset,
        );
        processed_count += line_buffer.len();
    }

    eprintln!("Completed processing {} kmers", processed_count);
}

/// Process a chunk of lines in parallel
fn process_chunk(
    lines: &[String],
    sample_names: &[String],
    minsize: u32,
    zscore_cutoff: f32,
    method: Method,
    mincluster: usize,
    subset: &Option<Vec<String>>,
) {
    use rayon::prelude::*;

    // Process lines in parallel and collect results
    let results: Vec<_> = lines
        .par_iter()
        .filter_map(|line| {
            let splitline: Vec<&str> = line.split('\t').collect();
            if splitline.len() < 3 {
                return None;
            }

            let (chrom, begin, end) = (splitline[0], splitline[1], splitline[2]);

            if let Some(values) = get_repeat_lengths(&splitline, minsize) {
                let expanded = match method {
                    Method::Zscore => z_score_outliers(&values, sample_names, zscore_cutoff),
                    Method::Dbscan => dbscan_outliers(&values, sample_names, mincluster),
                };

                if !expanded.is_empty() {
                    debug!(
                        "chrom: {}, begin: {}, end: {}, N_expanded: {}, expanded: {:?}",
                        chrom,
                        begin,
                        end,
                        expanded.len(),
                        expanded
                    );

                    // Check subset filtering
                    if let Some(subset) = subset {
                        if expanded.iter().any(|sample| subset.contains(sample)) {
                            return Some((
                                chrom.to_string(),
                                begin.to_string(),
                                end.to_string(),
                                expanded,
                            ));
                        }
                    } else {
                        return Some((
                            chrom.to_string(),
                            begin.to_string(),
                            end.to_string(),
                            expanded,
                        ));
                    }
                }
            }
            None
        })
        .collect();

    // Output results in order (important for deterministic output)
    for (chrom, begin, end, expanded) in results {
        let expanded_str = expanded.join(",");
        println!("{}\t{}\t{}\t{}", chrom, begin, end, expanded_str);
    }
}

/// Process a chunk of kmer lines in parallel
fn process_kmer_chunk(
    lines: &[String],
    sample_names: &[String],
    minsize: u32,
    zscore_cutoff: f32,
    method: Method,
    mincluster: usize,
    subset: &Option<Vec<String>>,
) {
    // Process lines in parallel and collect results
    let results: Vec<_> = lines
        .par_iter()
        .filter_map(|line| {
            let splitline: Vec<&str> = line.split('\t').collect();
            if splitline.len() < 2 {
                return None;
            }

            let kmer = splitline[0];

            if let Some(values) = get_kmer_frequencies(&splitline, minsize) {
                let expanded = match method {
                    Method::Zscore => z_score_outliers(&values, sample_names, zscore_cutoff),
                    Method::Dbscan => dbscan_outliers(&values, sample_names, mincluster),
                };

                if !expanded.is_empty() {
                    debug!(
                        "kmer: {}, N_expanded: {}, expanded: {:?}",
                        kmer,
                        expanded.len(),
                        expanded
                    );

                    // Check subset filtering
                    if let Some(subset) = subset {
                        if expanded.iter().any(|sample| subset.contains(sample)) {
                            return Some((kmer.to_string(), expanded));
                        }
                    } else {
                        return Some((kmer.to_string(), expanded));
                    }
                }
            }
            None
        })
        .collect();

    // Output results in order (important for deterministic output)
    for (kmer, expanded) in results {
        let expanded_str = expanded.join(",");
        println!("{}\t{}", kmer, expanded_str);
    }
}

/// Get repeat lengths for outlier analysis
fn get_repeat_lengths(line: &[&str], minsize: u32) -> Option<Vec<f32>> {
    if line.len() < 4 {
        return None;
    }

    let mut max_value = 0.0f32;
    let mut values = Vec::with_capacity(line.len() - 3);

    // Single pass to parse and find max
    for field in line.iter().skip(3) {
        let value = if field.eq_ignore_ascii_case("nan") || field.is_empty() {
            0.0
        } else {
            field.parse().unwrap_or(0.0)
        };

        if value > max_value {
            max_value = value;
        }
        values.push(value);
    }

    // Early exit if max value is too small
    if max_value < minsize as f32 {
        None
    } else {
        Some(values)
    }
}

/// Get kmer frequencies for outlier analysis
fn get_kmer_frequencies(line: &[&str], minsize: u32) -> Option<Vec<f32>> {
    if line.len() < 2 {
        return None;
    }

    let mut max_value = 0.0f32;
    let mut values = Vec::with_capacity(line.len() - 1);

    // Single pass to parse and find max (skip first column which is kmer name)
    for field in line.iter().skip(1) {
        let value = if field.eq_ignore_ascii_case("nan") || field.is_empty() {
            0.0
        } else {
            field.parse().unwrap_or(0.0)
        };

        if value > max_value {
            max_value = value;
        }
        values.push(value);
    }

    // For kmer frequencies, we use a lower threshold since they're normalized frequencies
    // Convert minsize to a frequency threshold (e.g., minsize/1000000 for per-million normalization)
    let freq_threshold = (minsize as f32) / 1000000.0;

    // Early exit if max value is too small
    if max_value < freq_threshold {
        None
    } else {
        Some(values)
    }
}

/// Z-score outlier detection
fn z_score_outliers(values: &[f32], sample_names: &[String], zscore_cutoff: f32) -> Vec<String> {
    let (values_mean, values_std_dev) = streaming_stats(values);
    debug!("mean: {}, std_dev: {}", values_mean, values_std_dev);

    if values_std_dev == 0.0 {
        return Vec::new();
    }

    let mut outliers = Vec::new();
    for (index, &value) in values.iter().enumerate() {
        if ((value - values_mean) / values_std_dev) >= zscore_cutoff {
            // Optimize string processing
            let sample_name = &sample_names[index];
            let clean_name = clean_sample_name(sample_name);
            outliers.push(clean_name.to_string());
        }
    }
    outliers
}

/// DBSCAN outlier detection
fn dbscan_outliers(values: &[f32], sample_names: &[String], mincluster: usize) -> Vec<String> {
    // the parameters for the dbscan model are as used by the schizophrenia STR outlier paper
    // however, the eps parameter is set as minimally 10
    let eps = max(2 * mode(values), 10) as f64;

    // Convert to format expected by DBSCAN
    let dbscan_values: Vec<Vec<f32>> = values.iter().map(|&value| vec![value]).collect();

    let model = Model::new(eps, mincluster);
    let output = model.run(&dbscan_values);
    debug!("eps: {}, mincluster: {}", eps, mincluster);
    debug!("output: {:?}", output);

    let mut outliers = Vec::new();
    for (index, &classification) in output.iter().enumerate() {
        if matches!(classification, Noise) {
            // Optimize string processing
            let sample_name = &sample_names[index];
            let clean_name = clean_sample_name(sample_name);
            outliers.push(clean_name.to_string());
        }
    }
    outliers
}

/// Calculate mode of values
fn mode(values: &[f32]) -> usize {
    let mut counts = std::collections::HashMap::new();
    for &value in values.iter().filter(|&&value| value > 0.0) {
        *counts.entry(value as usize).or_insert(0) += 1;
    }

    counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(value, _)| value)
        .unwrap_or(10) // Return 10 as fallback instead of panicking
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dbscan_outliers() {
        let values = vec![1.0, 2.0, 2.0, 3.0, 1.0, 5.0, 3.0, 2.0, 2.0, 1.0, 120.0];
        let samples: Vec<String> = vec![
            "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11",
        ]
        .into_iter()
        .map(|s| s.to_string())
        .collect();
        let expected = vec!["s11"];
        let mincluster = values.len().ilog2() as usize;
        assert_eq!(dbscan_outliers(&values, &samples, mincluster), expected);
    }

    #[test]
    fn test_z_score_outliers() {
        let values = vec![1.0, 2.0, 2.0, 3.0, 1.0, 5.0, 3.0, 2.0, 2.0, 1.0, 120.0];
        let samples: Vec<String> = vec![
            "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11",
        ]
        .into_iter()
        .map(|s| s.to_string())
        .collect();
        let expected = vec!["s11"];
        let zscore_cutoff = 2.0;
        assert_eq!(z_score_outliers(&values, &samples, zscore_cutoff), expected);
    }

    #[test]
    fn test_parse_sample_input_single() {
        let input = "Sample123";
        let result = parse_sample_input(input);
        assert_eq!(result, vec!["Sample123"]);
    }

    #[test]
    fn test_parse_sample_input_comma_separated() {
        let input = "Sample1,Sample2,Sample3";
        let result = parse_sample_input(input);
        assert_eq!(result, vec!["Sample1", "Sample2", "Sample3"]);
    }

    #[test]
    fn test_parse_sample_input_comma_separated_with_spaces() {
        let input = "Sample1, Sample2 , Sample3";
        let result = parse_sample_input(input);
        assert_eq!(result, vec!["Sample1", "Sample2", "Sample3"]);
    }

    #[test]
    fn test_clean_sample_name() {
        assert_eq!(clean_sample_name("Sample_H1"), "Sample");
        assert_eq!(clean_sample_name("Sample_H2"), "Sample");
        assert_eq!(clean_sample_name("Sample"), "Sample");
        assert_eq!(clean_sample_name("Sample_H3"), "Sample_H3");
    }
}
