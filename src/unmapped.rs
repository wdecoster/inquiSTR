//! Kmer frequency analysis of unmapped reads
//!
//! This module provides functionality to count kmer frequencies in unmapped reads
//! from BAM/CRAM files. It extracts all unmapped reads and counts occurrences of
//! kmers of all sizes from 2 to a specified maximum length.

use crate::bam_utils::setup_index_caching;
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use rust_htslib::bam::{self, Read};
use std::path::Path;
use url::Url;

/// Count kmer frequencies in unmapped reads from a BAM/CRAM file
pub fn count_unmapped_kmers(
    bam_path: String,
    klength: usize,
    sample_name: Option<String>,
    reference: Option<String>,
    threads: usize,
    target_kmer: Option<String>,
) {
    info!("Starting kmer counting for unmapped reads");
    info!("BAM file: {}", bam_path);
    info!("Maximum kmer length: {}", klength);
    info!("Threads: {}", threads);

    if let Some(ref target) = target_kmer {
        info!("Target kmer: {}", target);
        // If target kmer is specified, only count that specific kmer
        count_target_kmer(&bam_path, target, sample_name, &reference, threads);
        return;
    }

    if klength < 2 {
        panic!("klength must be at least 2");
    }

    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .expect("Failed to build thread pool");

    // Determine sample name
    let sample_name = sample_name.unwrap_or_else(|| {
        Path::new(&bam_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("sample")
            .to_string()
    });

    // Try to get unmapped read count from index statistics first
    let counts_from_index = get_counts_from_index(&bam_path);

    // Stream unmapped reads and count kmers in batches
    let (kmer_counts, total_reads) =
        stream_unmapped_reads_and_count_kmers(&bam_path, &reference, counts_from_index, klength);

    info!("Total reads in file: {}", total_reads);
    info!("Processed unmapped reads with streaming approach");

    if kmer_counts
        .iter()
        .all(|counts| counts.iter().all(|&c| c == 0))
    {
        warn!("No valid unmapped reads found in the BAM file");
        output_empty_results(&sample_name, klength);
        return;
    }

    // Canonicalize the counts (sum rotations together)
    let canonical_counts = canonicalize_kmer_counts(&kmer_counts, klength);

    // Output results
    output_results(&canonical_counts, &sample_name, klength, total_reads);
}

/// Get unmapped read count and total read count from BAM index statistics using rust-htslib
/// Returns Some((unmapped_count, total_count)) if successful, None if index is not available or other error
fn get_counts_from_index(bam_path: &str) -> Option<(u64, u64)> {
    // Set up index caching before opening the file
    setup_index_caching(bam_path);

    // Try to create an IndexedReader and get index statistics
    let reader_result = if bam_path.starts_with("s3")
        || bam_path.starts_with("https://")
        || bam_path.starts_with("http://")
        || bam_path.starts_with("ftp://")
    {
        info!("Opening remote BAM file for index statistics...");
        rust_htslib::bam::IndexedReader::from_url(
            &Url::parse(bam_path).expect("Failed to parse URL"),
        )
    } else {
        info!("Opening local BAM file for index statistics...");
        rust_htslib::bam::IndexedReader::from_path(bam_path)
    };

    match reader_result {
        Ok(mut reader) => {
            info!("IndexedReader created successfully, calling index_stats()...");
            match reader.index_stats() {
                Ok(stats) => {
                    // index_stats returns Vec<(i64, u64, u64, u64)> where:
                    // (target_id, length, mapped_reads, unmapped_reads)
                    // The last entry corresponds to unmapped reads for the entire file with tid=-1
                    // Look for the entry with target_id = -1 (which contains total unmapped reads)
                    info!("index_stats() returned {} entries.", stats.len());
                    let mut total_reads = 0u64;
                    let mut total_unmapped = 0u64;

                    for (tid, _, mapped, unmapped) in &stats {
                        total_reads += mapped + unmapped;
                        if *tid == -1 {
                            total_unmapped = *unmapped;
                        }
                    }

                    info!(
                        "Found {} unmapped reads out of {} total reads from index statistics",
                        total_unmapped, total_reads
                    );
                    Some((total_unmapped, total_reads))
                }
                Err(e) => {
                    info!("Failed to get index statistics: {}", e);
                    None
                }
            }
        }
        Err(e) => {
            info!("Failed to open BAM with index: {}.", e);
            None
        }
    }
}

/// Stream unmapped reads and count kmers in batches (memory efficient approach)
/// Uses streaming + rayon parallel iterators (consistent with rest of inquiSTR)
fn stream_unmapped_reads_and_count_kmers(
    bam_path: &str,
    reference: &Option<String>,
    counts_from_index: Option<(u64, u64)>,
    klength: usize,
) -> (Vec<Vec<u64>>, u64) {
    info!("Starting streaming kmer counting for unmapped reads");

    // Initialize global kmer counts
    let mut global_kmer_counts = Vec::new();
    for k in 2..=klength {
        let num_kmers = 4_usize.pow(k as u32);
        global_kmer_counts.push(vec![0u64; num_kmers]);
    }

    // Progress tracking
    let total_unmapped_expected = counts_from_index.map(|(unmapped, _)| unmapped);
    let progress = if let Some(expected) = total_unmapped_expected {
        let pb = ProgressBar::new(expected);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} unmapped reads ({percent}%)")
                .expect("Invalid progress template")
        );
        Some(pb)
    } else {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] Processing reads: {msg}")
                .expect("Invalid progress template"),
        );
        pb.set_message("0 reads processed");
        Some(pb)
    };

    // Collect unmapped sequences in batches
    let batches = collect_unmapped_sequences_in_batches(
        bam_path,
        reference,
        counts_from_index,
        progress.clone(),
    );
    let total_reads = batches.1;
    let sequence_batches = batches.0;

    info!("Collected {} batches of unmapped sequences", sequence_batches.len());

    // Process batches in parallel using rayon (like the rest of inquiSTR)
    let batch_results: Vec<Vec<Vec<u64>>> = sequence_batches
        .into_par_iter()
        .map(|batch| count_kmers_in_batch(&batch, klength))
        .collect();

    // Merge results from all batches
    for batch_counts in batch_results {
        for (k_idx, batch_k_counts) in batch_counts.iter().enumerate() {
            for (kmer_idx, &count) in batch_k_counts.iter().enumerate() {
                global_kmer_counts[k_idx][kmer_idx] += count;
            }
        }
    }

    if let Some(pb) = progress {
        pb.finish_with_message("Completed streaming kmer counting");
    }

    (global_kmer_counts, total_reads)
}

/// Collect unmapped sequences in batches for processing
fn collect_unmapped_sequences_in_batches(
    bam_path: &str,
    reference: &Option<String>,
    counts_from_index: Option<(u64, u64)>,
    progress: Option<ProgressBar>,
) -> (Vec<Vec<Vec<u8>>>, u64) {
    // Try indexed approach first
    if let Some((expected_unmapped, _expected_total)) = counts_from_index {
        if let Some((batches, total_reads)) =
            try_collect_indexed(bam_path, reference, expected_unmapped, progress.clone())
        {
            return (batches, total_reads);
        }
    }

    // Fall back to full file traversal
    collect_fallback(bam_path, reference, progress)
}

/// Count kmers in a batch of sequences
fn count_kmers_in_batch(batch: &[Vec<u8>], klength: usize) -> Vec<Vec<u64>> {
    // Initialize kmer counts for this batch
    let mut batch_kmer_counts = Vec::new();
    for k in 2..=klength {
        let num_kmers = 4_usize.pow(k as u32);
        batch_kmer_counts.push(vec![0u64; num_kmers]);
    }

    // Count kmers in all sequences in this batch
    for sequence in batch {
        for k in 2..=klength {
            if sequence.len() >= k {
                count_kmers_in_read(sequence, k, &mut batch_kmer_counts[k - 2]);
            }
        }
    }

    batch_kmer_counts
}

/// Try to collect sequences using indexed approach, returns batches and total_reads if successful
fn try_collect_indexed(
    bam_path: &str,
    reference: &Option<String>,
    expected_unmapped: u64,
    progress: Option<ProgressBar>,
) -> Option<(Vec<Vec<Vec<u8>>>, u64)> {
    // Set up index caching before opening the file
    setup_index_caching(bam_path);

    // Use the same logic as fetch_unmapped_reads_indexed but stream the results
    let reader_result = if bam_path.starts_with("s3")
        || bam_path.starts_with("https://")
        || bam_path.starts_with("http://")
        || bam_path.starts_with("ftp://")
    {
        info!("Creating IndexedReader for remote file...");
        rust_htslib::bam::IndexedReader::from_url(
            &url::Url::parse(bam_path).expect("Failed to parse URL"),
        )
    } else {
        info!("Creating IndexedReader for local file...");
        rust_htslib::bam::IndexedReader::from_path(bam_path)
    };

    let mut reader = match reader_result {
        Ok(reader) => reader,
        Err(e) => {
            info!("Could not create IndexedReader for indexed approach: {}", e);
            return None;
        }
    };

    // Set reference for CRAM files
    if bam_path.ends_with(".cram") {
        if let Some(ref_path) = reference {
            if let Err(e) = reader.set_reference(ref_path) {
                info!("Failed to set reference for IndexedReader: {}", e);
                return None;
            }
        }
    }

    // Try to fetch unmapped reads specifically
    match reader.fetch(bam::FetchDefinition::Unmapped) {
        Ok(()) => {
            info!("Successfully initiated indexed unmapped read fetch");
            let mut current_batch = Vec::new();
            let mut batches = Vec::new();
            let mut total_unmapped = 0u64;
            const BATCH_SIZE: usize = 1000;

            while let Some(result) = reader.records().next() {
                match result {
                    Ok(rec) => {
                        if rec.is_unmapped() {
                            total_unmapped += 1;
                            let seq = rec.seq().as_bytes();
                            if !seq.is_empty()
                                && seq.len() < 1_000_000
                                && seq.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
                            {
                                current_batch.push(seq);

                                // Complete batch when full
                                if current_batch.len() >= BATCH_SIZE {
                                    batches.push(current_batch.clone());
                                    current_batch.clear();

                                    // Update progress
                                    if let Some(ref pb) = progress {
                                        pb.set_position(total_unmapped);
                                    }
                                }
                            }
                        }

                        // Safety check
                        if total_unmapped > 50_000_000 {
                            info!("Safety limits reached, stopping");
                            break;
                        }
                    }
                    Err(e) => {
                        info!("Error in indexed fetch: {}, falling back", e);
                        return None;
                    }
                }
            }

            // Add final batch if not empty
            if !current_batch.is_empty() {
                batches.push(current_batch);
            }

            // Use expected total from index, as indexed approach may not see all reads
            Some((batches, expected_unmapped))
        }
        Err(e) => {
            info!("Could not fetch unmapped reads with indexed approach: {}", e);
            None
        }
    }
}

/// Collect sequences using fallback full file traversal
fn collect_fallback(
    bam_path: &str,
    reference: &Option<String>,
    progress: Option<ProgressBar>,
) -> (Vec<Vec<Vec<u8>>>, u64) {
    info!("Using fallback file traversal for sequence production");

    // Set up index caching before opening the file
    setup_index_caching(bam_path);

    let mut reader = if bam_path.starts_with("http")
        || bam_path.starts_with("https")
        || bam_path.starts_with("ftp")
        || bam_path.starts_with("s3")
    {
        info!("Opening BAM from URL...");
        bam::Reader::from_url(&bam_path.parse().expect("Invalid URL"))
            .expect("Failed to open BAM from URL")
    } else {
        info!("Opening local BAM file...");
        bam::Reader::from_path(bam_path).expect("Failed to open BAM file")
    };

    // Set reference for CRAM files
    if bam_path.ends_with(".cram") {
        if let Some(ref_path) = reference {
            reader
                .set_reference(ref_path)
                .expect("Failed to set reference");
        }
    }

    let mut current_batch = Vec::new();
    let mut batches = Vec::new();
    let mut total_reads = 0u64;
    let mut unmapped_count = 0u64;
    let mut record = bam::Record::new();
    const BATCH_SIZE: usize = 1000;

    loop {
        match reader.read(&mut record) {
            Some(Ok(())) => {
                total_reads += 1;

                if record.is_unmapped() {
                    unmapped_count += 1;
                    let seq = record.seq().as_bytes();
                    if !seq.is_empty()
                        && seq.len() < 1_000_000
                        && seq.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
                    {
                        current_batch.push(seq);

                        // Complete batch when full
                        if current_batch.len() >= BATCH_SIZE {
                            batches.push(current_batch.clone());
                            current_batch.clear();
                        }
                    }
                }

                // Update progress periodically
                if total_reads.is_multiple_of(10_000) {
                    if let Some(ref pb) = progress {
                        pb.set_message(format!(
                            "{} reads processed, {} unmapped",
                            total_reads, unmapped_count
                        ));
                        pb.tick();
                    }
                }

                // Safety checks
                if total_reads > 1_000_000_000 || unmapped_count > 10_000_000 {
                    warn!("Safety limits reached, stopping");
                    break;
                }
            }
            Some(Err(e)) => {
                warn!("Error reading BAM record: {}", e);
                break;
            }
            None => {
                info!("Reached end of BAM file");
                break;
            }
        }
    }

    // Add final batch if not empty
    if !current_batch.is_empty() {
        batches.push(current_batch);
    }

    (batches, unmapped_count)
}

/// Count kmers of a specific size in a single read
fn count_kmers_in_read(read: &[u8], k: usize, counts: &mut [u64]) {
    for i in 0..=(read.len() - k) {
        let kmer = &read[i..i + k];
        let index = kmer_to_index(kmer);
        counts[index] += 1;
    }
}

/// Convert a kmer to its canonical form (lexicographically smallest rotation)
fn get_canonical_kmer(kmer: &[u8]) -> Vec<u8> {
    let mut rotations = Vec::new();

    // Generate all rotations
    for i in 0..kmer.len() {
        let mut rotation = Vec::with_capacity(kmer.len());
        rotation.extend_from_slice(&kmer[i..]);
        rotation.extend_from_slice(&kmer[..i]);
        rotations.push(rotation);
    }

    // Return the lexicographically smallest
    rotations.into_iter().min().unwrap()
}

/// Count occurrences of a specific target kmer in unmapped reads
fn count_target_kmer(
    bam_path: &str,
    target_kmer: &str,
    sample_name: Option<String>,
    reference: &Option<String>,
    threads: usize,
) {
    // Validate target kmer contains only ACGT
    let target_upper = target_kmer.to_uppercase();
    if !target_upper
        .bytes()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
    {
        panic!("Target kmer must contain only A, C, G, T characters");
    }

    let target_bytes = target_upper.as_bytes();
    let k = target_bytes.len();

    if k < 1 {
        panic!("Target kmer must be at least 1 base long");
    }

    info!("Counting target kmer '{}' (length {})", target_upper, k);

    // Generate all rotations of the target kmer
    let mut rotations = Vec::new();
    for i in 0..k {
        let mut rotation = Vec::with_capacity(k);
        rotation.extend_from_slice(&target_bytes[i..]);
        rotation.extend_from_slice(&target_bytes[..i]);
        rotations.push(rotation);
    }

    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .expect("Failed to build thread pool");

    // Determine sample name
    let sample_name = sample_name.unwrap_or_else(|| {
        Path::new(bam_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("sample")
            .to_string()
    });

    // Try to get unmapped read count from index statistics
    let counts_from_index = get_counts_from_index(bam_path);

    // Stream unmapped reads and count target kmer
    let (total_count, total_reads) =
        stream_and_count_target_kmer(bam_path, reference, counts_from_index, &rotations);

    // Output results
    println!("Sample\tTarget_Kmer\tCanonical_Kmer\tKmer_Length\tCount\tTotal_Reads\tFrequency");
    let canonical = String::from_utf8(get_canonical_kmer(target_bytes)).unwrap();
    let frequency = if total_reads > 0 {
        total_count as f64 / total_reads as f64
    } else {
        0.0
    };
    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
        sample_name, target_upper, canonical, k, total_count, total_reads, frequency
    );
}

/// Stream unmapped reads and count occurrences of target kmer rotations
fn stream_and_count_target_kmer(
    bam_path: &str,
    reference: &Option<String>,
    counts_from_index: Option<(u64, u64)>,
    rotations: &[Vec<u8>],
) -> (u64, u64) {
    const BATCH_SIZE: usize = 10000;

    // Set up progress bar
    let progress_bar = if let Some((unmapped_count, _)) = counts_from_index {
        let pb = ProgressBar::new(unmapped_count);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} reads ({eta})")
                .expect("Failed to set progress bar template")
                .progress_chars("=>-"),
        );
        Some(pb)
    } else {
        None
    };

    // Set up index caching before opening the file
    setup_index_caching(bam_path);

    // Open BAM reader
    let mut reader = if bam_path.starts_with("http")
        || bam_path.starts_with("https")
        || bam_path.starts_with("ftp")
        || bam_path.starts_with("s3")
    {
        info!("Opening BAM from URL...");
        bam::Reader::from_url(&bam_path.parse().expect("Invalid URL"))
            .expect("Failed to open BAM from URL")
    } else {
        info!("Opening local BAM file...");
        bam::Reader::from_path(bam_path).expect("Failed to open BAM file")
    };

    // Set reference for CRAM files
    if bam_path.ends_with(".cram") {
        if let Some(ref_path) = reference {
            reader
                .set_reference(ref_path)
                .expect("Failed to set reference");
        }
    }

    let mut batch = Vec::with_capacity(BATCH_SIZE);
    let mut batch_counts = Vec::new();
    let mut total_reads = 0u64;

    for result in reader.records() {
        match result {
            Ok(record) => {
                if record.is_unmapped() {
                    total_reads += 1;
                    if let Some(ref pb) = progress_bar {
                        pb.inc(1);
                    }
                    batch.push(record.seq().as_bytes());

                    if batch.len() >= BATCH_SIZE {
                        // Process batch and collect count
                        batch_counts.push(batch.clone());
                        batch.clear();
                    }
                }
            }
            Err(e) => {
                warn!("Error reading record: {}", e);
            }
        }
    }

    // Add remaining reads to batches
    if !batch.is_empty() {
        batch_counts.push(batch);
    }

    if let Some(pb) = progress_bar {
        pb.finish_with_message("Done");
    }

    // Count target kmer in parallel batches
    let total_count: u64 = batch_counts
        .into_par_iter()
        .map(|batch| count_target_in_batch(&batch, rotations))
        .sum();

    (total_count, total_reads)
}

/// Count occurrences of target kmer (any rotation) in a batch of reads
fn count_target_in_batch(batch: &[Vec<u8>], rotations: &[Vec<u8>]) -> u64 {
    let mut count = 0u64;
    let k = rotations[0].len();

    for read in batch {
        if read.len() < k {
            continue;
        }

        // Check each position in the read
        for i in 0..=(read.len() - k) {
            let kmer = &read[i..i + k];
            // Check if this kmer matches any rotation
            if rotations.iter().any(|rot| rot.as_slice() == kmer) {
                count += 1;
            }
        }
    }

    count
}

/// Canonicalize kmer counts by summing all rotations together
/// Returns a map from canonical kmer to total count
fn canonicalize_kmer_counts(
    kmer_counts: &[Vec<u64>],
    klength: usize,
) -> Vec<std::collections::HashMap<String, u64>> {
    use std::collections::HashMap;

    let mut canonical_counts = Vec::new();

    for k in 2..=klength {
        let k_idx = k - 2;
        let mut canonical_map = HashMap::new();

        // Iterate through all possible kmers for this k
        let num_kmers = 4_usize.pow(k as u32);
        for kmer_idx in 0..num_kmers {
            let count = kmer_counts[k_idx][kmer_idx];
            if count > 0 {
                let kmer_bytes = index_to_kmer_bytes(kmer_idx, k);
                let canonical_kmer = get_canonical_kmer(&kmer_bytes);
                let canonical_str = String::from_utf8(canonical_kmer).unwrap();

                *canonical_map.entry(canonical_str).or_insert(0) += count;
            }
        }

        canonical_counts.push(canonical_map);
    }

    canonical_counts
}

/// Convert a kmer to its vector index
fn kmer_to_index(kmer: &[u8]) -> usize {
    let mut index = 0;
    for &base in kmer {
        index *= 4;
        index += match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => panic!("Invalid base in kmer: {}", base as char),
        };
    }
    index
}

/// Convert an index back to a kmer string
fn index_to_kmer(index: usize, k: usize) -> String {
    let kmer_bytes = index_to_kmer_bytes(index, k);
    String::from_utf8(kmer_bytes).unwrap()
}

/// Convert an index back to a kmer as bytes (more efficient)
fn index_to_kmer_bytes(mut index: usize, k: usize) -> Vec<u8> {
    let mut kmer = Vec::with_capacity(k);
    for _ in 0..k {
        let base = match index % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
        kmer.push(base);
        index /= 4;
    }
    kmer.reverse();
    kmer
}

/// Output results in TSV format with only canonical kmers
fn output_results(
    canonical_counts: &[std::collections::HashMap<String, u64>],
    sample_name: &str,
    klength: usize,
    total_reads: u64,
) {
    println!("kmer\t{}", sample_name);

    for k in 2..=klength {
        let k_idx = k - 2;
        let mut canonical_kmers: Vec<_> = canonical_counts[k_idx].keys().cloned().collect();
        canonical_kmers.sort(); // Sort alphabetically

        for kmer in canonical_kmers {
            let count = canonical_counts[k_idx][&kmer];
            let normalized_count = count as f64 / total_reads as f64;
            println!("{}\t{:.6}", kmer, normalized_count);
        }
    }
}

/// Output empty results when no unmapped reads are found (only canonical kmers)
fn output_empty_results(sample_name: &str, klength: usize) {
    println!("kmer\t{}", sample_name);

    for k in 2..=klength {
        let all_kmers = generate_all_kmers(k);
        let mut canonical_kmers = std::collections::HashSet::new();

        // Get only canonical forms
        for kmer in all_kmers {
            let canonical = get_canonical_kmer(kmer.as_bytes());
            let canonical_str = String::from_utf8(canonical).unwrap();
            canonical_kmers.insert(canonical_str);
        }

        let mut sorted_canonical: Vec<_> = canonical_kmers.into_iter().collect();
        sorted_canonical.sort();

        for kmer in sorted_canonical {
            println!("{}\t0.000000", kmer);
        }
    }
}

/// Generate all possible kmers of a given length in alphabetical order
fn generate_all_kmers(k: usize) -> Vec<String> {
    let num_kmers = 4_usize.pow(k as u32);
    (0..num_kmers).map(|i| index_to_kmer(i, k)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_to_index() {
        assert_eq!(kmer_to_index(b"AA"), 0);
        assert_eq!(kmer_to_index(b"AC"), 1);
        assert_eq!(kmer_to_index(b"AG"), 2);
        assert_eq!(kmer_to_index(b"AT"), 3);
        assert_eq!(kmer_to_index(b"CA"), 4);
    }

    #[test]
    fn test_index_to_kmer() {
        assert_eq!(index_to_kmer(0, 2), "AA");
        assert_eq!(index_to_kmer(1, 2), "AC");
        assert_eq!(index_to_kmer(2, 2), "AG");
        assert_eq!(index_to_kmer(3, 2), "AT");
        assert_eq!(index_to_kmer(4, 2), "CA");
    }

    #[test]
    fn test_index_to_kmer_bytes() {
        assert_eq!(index_to_kmer_bytes(0, 2), b"AA");
        assert_eq!(index_to_kmer_bytes(1, 2), b"AC");
        assert_eq!(index_to_kmer_bytes(2, 2), b"AG");
        assert_eq!(index_to_kmer_bytes(3, 2), b"AT");
        assert_eq!(index_to_kmer_bytes(4, 2), b"CA");
    }

    #[test]
    fn test_canonical_kmer() {
        assert_eq!(get_canonical_kmer(b"CAG"), b"AGC");
        assert_eq!(get_canonical_kmer(b"AGC"), b"AGC");
        assert_eq!(get_canonical_kmer(b"GCA"), b"AGC");
        assert_eq!(get_canonical_kmer(b"AT"), b"AT");
        assert_eq!(get_canonical_kmer(b"TA"), b"AT");
    }

    #[test]
    fn test_generate_all_kmers() {
        let dimers = generate_all_kmers(2);
        assert_eq!(dimers.len(), 16);
        assert_eq!(dimers[0], "AA");
        assert_eq!(dimers[1], "AC");
        assert_eq!(dimers[15], "TT");
    }

    #[test]
    fn test_kmer_counting_in_read() {
        // Test with a longer, more comprehensive sequence
        let read = b"AACAGTCGTTAGCCCGATTTACG"; // 23 bases
        let mut dimer_counts = vec![0u64; 16]; // For dimers (4^2 = 16)
        let mut trimer_counts = vec![0u64; 64]; // For trimers (4^3 = 64)

        // Count dimers
        count_kmers_in_read(read, 2, &mut dimer_counts);

        // Count trimers
        count_kmers_in_read(read, 3, &mut trimer_counts);

        // Expected dimers from "AACAGTCGTTAGCCCGATTTACG":
        // AA, AC, CA, AG, GT, TC, CG, GT, TT, TA, AG, GC, CC, CC, CG, GA, AT, TT, TT, TA, AC, CG
        // Unique dimers and their expected counts:
        assert_eq!(dimer_counts[kmer_to_index(b"AA")], 1); // position 0
        assert_eq!(dimer_counts[kmer_to_index(b"AC")], 2); // positions 1, 20
        assert_eq!(dimer_counts[kmer_to_index(b"CA")], 1); // position 2
        assert_eq!(dimer_counts[kmer_to_index(b"AG")], 2); // positions 3, 9
        assert_eq!(dimer_counts[kmer_to_index(b"GT")], 2); // positions 4, 7
        assert_eq!(dimer_counts[kmer_to_index(b"TC")], 1); // position 5
        assert_eq!(dimer_counts[kmer_to_index(b"CG")], 3); // positions 6, 14, 21
        assert_eq!(dimer_counts[kmer_to_index(b"TT")], 3); // positions 8, 17, 18
        assert_eq!(dimer_counts[kmer_to_index(b"TA")], 2); // positions 10, 19
        assert_eq!(dimer_counts[kmer_to_index(b"GC")], 1); // position 11
        assert_eq!(dimer_counts[kmer_to_index(b"CC")], 2); // positions 12, 13
        assert_eq!(dimer_counts[kmer_to_index(b"GA")], 1); // position 15
        assert_eq!(dimer_counts[kmer_to_index(b"AT")], 1); // position 16

        // Verify some trimers as well
        assert_eq!(trimer_counts[kmer_to_index(b"AAC")], 1); // position 0
        assert_eq!(trimer_counts[kmer_to_index(b"ACA")], 1); // position 1
        assert_eq!(trimer_counts[kmer_to_index(b"CAG")], 1); // position 2
        assert_eq!(trimer_counts[kmer_to_index(b"AGT")], 1); // position 3
        assert_eq!(trimer_counts[kmer_to_index(b"GCC")], 1); // position 11
        assert_eq!(trimer_counts[kmer_to_index(b"CCC")], 1); // position 12
        assert_eq!(trimer_counts[kmer_to_index(b"TTT")], 1); // position 17
        assert_eq!(trimer_counts[kmer_to_index(b"ACG")], 1); // position 20

        // Verify total counts make sense
        let total_dimers: u64 = dimer_counts.iter().sum();
        assert_eq!(total_dimers, 22); // 23 bases = 22 dimers

        let total_trimers: u64 = trimer_counts.iter().sum();
        assert_eq!(total_trimers, 21); // 23 bases = 21 trimers
    }

    #[test]
    fn test_batch_kmer_counting() {
        // Test the batch kmer counting function used in the streaming approach
        let batch = vec![b"AACC".to_vec(), b"GGTT".to_vec()];

        let batch_counts = count_kmers_in_batch(&batch, 2);

        // Should have counts for dimers only (k=2)
        assert_eq!(batch_counts.len(), 1); // Only k=2
        assert_eq!(batch_counts[0].len(), 16); // 4^2 possible dimers

        // From "AACC": AA, AC, CC
        // From "GGTT": GG, GT, TT
        assert_eq!(batch_counts[0][kmer_to_index(b"AA")], 1);
        assert_eq!(batch_counts[0][kmer_to_index(b"AC")], 1);
        assert_eq!(batch_counts[0][kmer_to_index(b"CC")], 1);
        assert_eq!(batch_counts[0][kmer_to_index(b"GG")], 1);
        assert_eq!(batch_counts[0][kmer_to_index(b"GT")], 1);
        assert_eq!(batch_counts[0][kmer_to_index(b"TT")], 1);

        // Test with multiple k values
        let batch_counts_multi = count_kmers_in_batch(&batch, 3);
        assert_eq!(batch_counts_multi.len(), 2); // k=2 and k=3
        assert_eq!(batch_counts_multi[0].len(), 16); // 4^2 for dimers
        assert_eq!(batch_counts_multi[1].len(), 64); // 4^3 for trimers
    }

    #[test]
    fn test_canonicalization() {
        // Create mock kmer counts with some rotations
        let mut kmer_counts_k2 = vec![0u64; 16]; // 4^2 = 16 dimers

        // Add counts for AC and CA (which should be canonicalized to AC)
        kmer_counts_k2[kmer_to_index(b"AC")] = 3;
        kmer_counts_k2[kmer_to_index(b"CA")] = 2;

        // Add counts for AT and TA (which should be canonicalized to AT)
        kmer_counts_k2[kmer_to_index(b"AT")] = 1;
        kmer_counts_k2[kmer_to_index(b"TA")] = 4;

        let kmer_counts = vec![kmer_counts_k2];
        let canonical_counts = canonicalize_kmer_counts(&kmer_counts, 2);

        // Should have AC with count 3+2=5
        assert_eq!(canonical_counts[0]["AC"], 5);

        // Should have AT with count 1+4=5
        assert_eq!(canonical_counts[0]["AT"], 5);

        // Should not have CA or TA as separate entries
        assert!(!canonical_counts[0].contains_key("CA"));
        assert!(!canonical_counts[0].contains_key("TA"));
    }
}
