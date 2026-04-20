use std::cmp::Ordering;
use std::collections::HashMap;

use crate::repeats::{ChromosomeMapper, RepeatInterval};

/// Compare two strings with natural (human) ordering, so "chr2" < "chr10".
/// Splits each string into runs of digits and non-digits and compares segment by segment.
pub fn human_cmp(a: &str, b: &str) -> Ordering {
    fn segments(s: &str) -> Vec<(bool, &str)> {
        let mut parts = Vec::new();
        let bytes = s.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            let is_digit = bytes[i].is_ascii_digit();
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() == is_digit {
                i += 1;
            }
            parts.push((is_digit, &s[start..i]));
        }
        parts
    }

    let pa = segments(a);
    let pb = segments(b);
    for (sa, sb) in pa.iter().zip(pb.iter()) {
        let ord = if sa.0 && sb.0 {
            // Both numeric: compare by length first (avoids overflow), then lexicographic
            sa.1.len().cmp(&sb.1.len()).then_with(|| sa.1.cmp(sb.1))
        } else {
            sa.1.cmp(sb.1)
        };
        if ord != Ordering::Equal {
            return ord;
        }
    }
    pa.len().cmp(&pb.len())
}

/// Represents a batch of nearby repeats to be processed together
/// Batches are created within chromosome boundaries and within distance threshold
#[derive(Clone)]
pub struct Batch {
    pub chromosome: String,
    pub start: u32,
    pub end: u32,
    pub repeats: Vec<RepeatInterval>,
}

impl Batch {
    /// Create a new batch from a chromosome and its repeats
    ///
    /// The batch region is automatically extended with 10bp padding on both ends
    /// to ensure reads overlapping the boundaries are captured.
    pub fn new(chromosome: String, repeats: Vec<RepeatInterval>) -> Self {
        let start = repeats
            .first()
            .map(|r| r.start.saturating_sub(10))
            .unwrap_or(0);
        let end = repeats
            .last()
            .map(|r| r.end.saturating_add(10))
            .unwrap_or(0);
        Self { chromosome, start, end, repeats }
    }
}

/// Create batches from repeats, grouping nearby repeats within chromosome boundaries
///
/// Groups STR targets that are within `batch_distance_threshold` bases of each other
/// into the same batch to optimize I/O operations. Batches never cross chromosome boundaries.
/// The returned batches are sorted by chromosome (human-sorted) and then by start position
/// for deterministic processing order.
///
/// # Arguments
/// * `repeats` - Vector of STR target intervals to batch
/// * `batch_distance_threshold` - Maximum distance (bp) between repeats to group into same batch
///
/// # Returns
/// Vector of batches, each containing nearby repeats from the same chromosome
pub fn create_batches(
    repeats: Vec<RepeatInterval>,
    batch_distance_threshold: u32,
    chrom_mapper: &ChromosomeMapper,
) -> Vec<Batch> {
    // Group repeats by chromosome first
    let mut by_chromosome: HashMap<u32, Vec<RepeatInterval>> = HashMap::new();
    for repeat in repeats {
        by_chromosome
            .entry(repeat.chrom_id)
            .or_default()
            .push(repeat);
    }

    // Sort repeats within each chromosome
    for chrom_repeats in by_chromosome.values_mut() {
        chrom_repeats.sort_by_key(|a| a.start);
    }

    // Create batches within each chromosome
    // Estimate number of batches (assume ~50 repeats per batch on average)
    let mut all_batches = Vec::new();

    for (chrom_id, chrom_repeats) in by_chromosome {
        let chromosome = chrom_mapper.get_name(chrom_id).to_string();
        let mut current_batch = Vec::new();
        let mut current_end = 0u32;

        for repeat in chrom_repeats {
            // Check if this repeat should be added to current batch or start a new batch
            let should_batch =
                !current_batch.is_empty() && repeat.start <= current_end + batch_distance_threshold;

            if should_batch {
                // Add to current batch and extend the region
                current_end = std::cmp::max(current_end, repeat.end.saturating_add(10));
                current_batch.push(repeat);
            } else {
                // Finish the current batch if it exists, pushing to all_batches
                if !current_batch.is_empty() {
                    all_batches
                        .push(Batch::new(chromosome.clone(), std::mem::take(&mut current_batch)));
                }

                // Start new batch
                current_end = repeat.end.saturating_add(10);
                current_batch = vec![repeat];
            }
        }

        // Add the final batch for this chromosome
        if !current_batch.is_empty() {
            all_batches.push(Batch::new(chromosome, current_batch));
        }
    }

    // Sort batches by chromosome and position for deterministic processing
    all_batches.sort_by(|a, b| human_cmp(&a.chromosome, &b.chromosome).then(a.start.cmp(&b.start)));

    all_batches
}
