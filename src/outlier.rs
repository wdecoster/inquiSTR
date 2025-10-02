use clap::ValueEnum;

use dbscan::Classification::*;
use dbscan::Model;
use log::debug;

use std::cmp::max;
use std::cmp::Ordering;
use std::io::BufRead;
use std::path::PathBuf;

/// Helper function to clean sample names by removing _H1/_H2 suffixes
#[inline]
fn clean_sample_name(sample_name: &str) -> &str {
    if sample_name.ends_with("_H1") || sample_name.ends_with("_H2") {
        &sample_name[..sample_name.len() - 3]
    } else {
        sample_name
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

pub fn outlier(
    combined: PathBuf,
    minsize: u32,
    zscore_cutoff: f32,
    method: Method,
    subset: Option<Vec<String>>,
    threads: usize,
) {
    // Configure thread pool
    let actual_threads = if threads == 0 {
        rayon::current_num_threads()
    } else {
        // Try to configure the global thread pool, but it might already be configured
        if rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .is_ok()
        {
            threads
        } else {
            // If we can't configure it, use the current number of threads
            eprintln!("Warning: Thread pool already configured, using existing settings");
            rayon::current_num_threads()
        }
    };

    eprintln!("Using {} threads for parallel processing", actual_threads);

    let file = crate::utils::reader(&combined.into_os_string().into_string().unwrap());
    let mut lines = file.lines();
    let header_line = lines.next().unwrap().unwrap();
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

            if let Some(values) = get_repeat_lengths_optimized(&splitline, minsize) {
                let expanded = match method {
                    Method::Zscore => {
                        z_score_outliers_optimized(&values, sample_names, zscore_cutoff)
                    }
                    Method::Dbscan => dbscan_outliers_optimized(&values, sample_names, mincluster),
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

/// Optimized version that reuses buffer and avoids unnecessary allocations
fn get_repeat_lengths_optimized(line: &[&str], minsize: u32) -> Option<Vec<f32>> {
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

#[allow(dead_code)]
fn get_repeat_lengths(line: &[&str], minsize: u32) -> Option<Vec<f32>> {
    let values: Vec<f32> = line
        .iter()
        .skip(3)
        .map(|number| number.parse().expect("Failed to parse number"))
        .collect();
    let values = values
        .iter()
        .map(|&value| if value.is_nan() { 0.0 } else { value })
        .collect::<Vec<f32>>();
    // Check if the maximum value is larger than the minimum size
    // If all values are NaN then the vector will contain only zeroes
    if values
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Less))
        .unwrap()
        < &(minsize as f32)
    {
        None
    } else {
        Some(values)
    }
}

#[allow(dead_code)]
fn z_score_outliers(values: Vec<f32>, samples: &[&str], zscore_cutoff: f32) -> Vec<String> {
    // calculate mean and std deviation of the STR lengths
    let (values_mean, values_std_dev) = streaming_stats(&values);
    debug!("mean: {}, std_dev: {}", values_mean, values_std_dev);
    // calculate the zscore for each haplotype and get the haplotype identifier based on the index if larger zscore > cutoff
    // intentionally this only selects for values that are larger than the mean
    // and therefore only for expansions, not contractions
    values
        .iter()
        .enumerate()
        .filter(|(_, &value)| ((value - values_mean) / values_std_dev) >= zscore_cutoff)
        .map(|(index, _)| samples[index].replace("_H1", "").replace("_H2", ""))
        .collect::<Vec<String>>()
}

/// Optimized z-score outlier detection with reduced allocations
fn z_score_outliers_optimized(
    values: &[f32],
    sample_names: &[String],
    zscore_cutoff: f32,
) -> Vec<String> {
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

#[allow(dead_code)]
fn dbscan_outliers(values: Vec<f32>, samples: &[&str], mincluster: usize) -> Vec<String> {
    // the parameters for the dbscan model are as used by the schizophrenia STR outlier paper (https://doi.org/10.1038/s41380-022-01857-4)
    // however, the eps parameter is set as minimally 10
    let eps = max(2 * mode(&values), 10) as f64;
    let values = values
        .iter()
        .map(|&value| vec![value])
        .collect::<Vec<Vec<f32>>>();
    let model = Model::new(eps, mincluster);
    let output = model.run(&values);
    debug!("eps: {}, mincluster: {}", eps, mincluster);
    debug!("output: {:?}", output);
    output
        .iter()
        .enumerate()
        .filter(|(_, &classification)| matches!(classification, Noise))
        .map(|(index, _)| samples[index].replace("_H1", "").replace("_H2", ""))
        .collect::<Vec<String>>()
}

/// Optimized DBSCAN outlier detection with reduced allocations
fn dbscan_outliers_optimized(
    values: &[f32],
    sample_names: &[String],
    mincluster: usize,
) -> Vec<String> {
    // the parameters for the dbscan model are as used by the schizophrenia STR outlier paper
    // however, the eps parameter is set as minimally 10
    let eps = max(2 * mode_optimized(values), 10) as f64;

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

/// Optimized mode calculation
fn mode_optimized(values: &[f32]) -> usize {
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

#[allow(dead_code)]
fn mode(values: &[f32]) -> usize {
    // calculate the mode of the STR lengths
    // as NaN values are replaced by 0.0, we need to filter out the 0.0 values
    // if not eps will often be 0
    let mut counts = std::collections::HashMap::new();
    for &value in values.iter().filter(|&&value| value > 0.0) {
        *counts.entry(value as usize).or_insert(0) += 1;
    }
    counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(value, _)| value)
        .expect("No mode found for repeat")
}

#[cfg(test)]
#[test]
fn test_dbscan_outliers() {
    let values = vec![1.0, 2.0, 2.0, 3.0, 1.0, 5.0, 3.0, 2.0, 2.0, 1.0, 120.0];
    let samples = vec![
        "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11",
    ];
    let expected = vec!["s11"];
    let mincluster = values.len().ilog2() as usize;
    assert_eq!(dbscan_outliers(values, &samples, mincluster), expected);
}

#[test]
fn test_z_score_outliers() {
    let values = vec![1.0, 2.0, 2.0, 3.0, 1.0, 5.0, 3.0, 2.0, 2.0, 1.0, 120.0];
    let samples = vec![
        "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11",
    ];
    let expected = vec!["s11"];
    let zscore_cutoff = 2.0;
    assert_eq!(z_score_outliers(values, &samples, zscore_cutoff), expected);
}
