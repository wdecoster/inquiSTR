use human_sort::compare as human_compare;
use log::{error, warn};
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;
use std::collections::HashMap;
use std::sync::Mutex;

use crate::bam_utils::{get_bam_reader, get_bam_reader_with_validation, is_accidental_2d};
use crate::call::{Call, Genotype, median_str_length};
use crate::errors::{InquiSTRError, InquiSTRResult};
use crate::repeats::RepeatInterval;
use std::sync::atomic::{AtomicBool, Ordering};

// Global mutex to serialize BAM/CRAM reader creation and fetch operations
// This prevents concurrent htslib index operations that can cause segfaults in cram_index_free
static BAM_READER_LOCK: Mutex<()> = Mutex::new(());

// Track whether phasing validation has been performed (only needed once, on first batch)
static PHASING_VALIDATED: AtomicBool = AtomicBool::new(false);

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
/// Batches are created within configurable distance threshold to optimize I/O
pub fn create_batches(repeats: Vec<RepeatInterval>, batch_distance_threshold: u32) -> Vec<Batch> {
    let total_repeats = repeats.len();

    // Group repeats by chromosome first
    let mut by_chromosome: HashMap<String, Vec<RepeatInterval>> = HashMap::new();

    for repeat in repeats {
        by_chromosome
            .entry(repeat.chrom.clone())
            .or_default()
            .push(repeat);
    }

    // Sort repeats within each chromosome
    for chrom_repeats in by_chromosome.values_mut() {
        chrom_repeats.sort_by(|a, b| a.start.cmp(&b.start));
    }

    // Create batches within each chromosome
    // Estimate number of batches (assume ~50 repeats per batch on average)
    let estimated_batches = total_repeats / 50 + 1;
    let mut all_batches = Vec::with_capacity(estimated_batches);

    for (chromosome, chrom_repeats) in by_chromosome {
        let mut current_batch = Vec::with_capacity(50);
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
                    all_batches.push(Batch::new(chromosome.clone(), current_batch));
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
    all_batches
        .sort_by(|a, b| human_compare(&a.chromosome, &b.chromosome).then(a.start.cmp(&b.start)));

    all_batches
}

/// Unified batch processor that handles a single batch and returns results
pub fn process_batch_worker(
    batch: Batch,
    bamp: &String,
    reference: &Option<String>,
    genotype: crate::call::GenotypeConfig,
) -> InquiSTRResult<Vec<Genotype>> {
    // Serialize BAM reader creation and fetch to prevent concurrent htslib index operations
    // This prevents segfaults in cram_index_free when multiple threads open the same CRAM
    // On the first batch (if not unphased), also validate phasing using the same reader
    let mut bam = {
        let _guard = BAM_READER_LOCK.lock().unwrap();

        // Check if we need to validate phasing (only once, on first batch)
        let need_validation = !genotype.unphased && !PHASING_VALIDATED.load(Ordering::Relaxed);

        if need_validation {
            PHASING_VALIDATED.store(true, Ordering::Relaxed);
            get_bam_reader_with_validation(bamp, reference)?
        } else {
            get_bam_reader(bamp, reference)?
        }
    };
    let mut results = Vec::new();

    // Serialize both header access and fetch for CRAM stability
    // (accessing header and then fetch without lock can cause issues with CRAM internal state)
    let fetch_result = {
        let _guard = BAM_READER_LOCK.lock().unwrap();

        let tid = match bam.header().tid(batch.chromosome.as_bytes()) {
            Some(t) => t,
            None => {
                warn!("Chromosome {} not found in BAM header", batch.chromosome);
                return Ok(results);
            }
        };

        bam.fetch((tid, batch.start, batch.end))
    };

    if let Err(e) = fetch_result {
        warn!(
            "Warning: Failed to fetch region {}:{}-{}: {}",
            batch.chromosome, batch.start, batch.end, e
        );
        return Ok(results);
    }

    // Smart overlap filtering: only store reads that intersect with STR targets
    let target_intervals_with_idx: Vec<(u32, u32, usize)> = batch
        .repeats
        .iter()
        .enumerate()
        .map(|(idx, repeat)| (repeat.start.saturating_sub(100), repeat.end + 100, idx))
        .collect();

    // Pre-allocate with estimated capacity (assume ~20-30 reads per target on average)
    let estimated_reads = batch.repeats.len() * 25;
    let mut batch_reads = Vec::with_capacity(estimated_reads);

    for record_result in bam.rc_records() {
        match record_result {
            Ok(record) => {
                let read_start = record.reference_start() as u32;
                let read_end = record.reference_end() as u32;
                let mapq = record.mapq();

                // Skip low quality reads early
                if mapq <= 10 {
                    continue;
                }

                // Check if read overlaps with any target interval
                // Use a small stack buffer to avoid heap allocation for common case (most reads hit 1-3 targets)
                let mut overlapping_targets_buf: smallvec::SmallVec<[(u32, u32, usize); 8]> =
                    smallvec::SmallVec::new();
                for &interval in &target_intervals_with_idx {
                    let (target_start, target_end, _) = interval;
                    if read_start < target_end && read_end > target_start {
                        overlapping_targets_buf.push(interval);
                    }
                }

                if !overlapping_targets_buf.is_empty() {
                    // Create ReadInfo with pre-computed STR calls for overlapping targets
                    let read_info = ReadInfo::from_record_with_targets(
                        (*record).clone(),
                        genotype.minlen,
                        &overlapping_targets_buf,
                        &batch.repeats,
                    );
                    batch_reads.push(read_info);
                }
            }
            Err(e) => {
                let error_str = e.to_string();
                if error_str.contains("CRC32 failure") || error_str.contains("truncated record") {
                    error!(
                        "CRAM format error in batch {}: {}. This usually indicates that the reference genome doesn't match the CRAM file or CRAM index is corrupted.",
                        batch.chromosome, error_str
                    );
                    return Err(InquiSTRError::new(format!(
                        "CRAM format error in batch {}: {}. This usually indicates that the reference genome doesn't match the CRAM file or CRAM index is corrupted.",
                        batch.chromosome, error_str
                    )));
                } else {
                    error!("Error reading BAM record in batch {}: {}", batch.chromosome, e);
                    return Err(InquiSTRError::new(format!(
                        "Error reading BAM record in batch {}: {}",
                        batch.chromosome, e
                    )));
                }
            }
        }
    }

    // Process each target in the batch using the lightweight read info
    for repeat in &batch.repeats {
        let result = process_target_from_read_info(&batch_reads, repeat, &genotype);

        match result {
            Ok((genotype, _had_hp_tags)) => {
                results.push(genotype);
            }
            Err(_e) => {
                // Output NaN genotype for failed targets (matches original behavior)
                let failed_genotype =
                    Genotype { repeat: repeat.clone(), phase1: f64::NAN, phase2: f64::NAN };
                results.push(failed_genotype);
            }
        }
    }

    Ok(results)
}

/// Memory-efficient read information - stores pre-computed STR calls instead of full BAM record
#[derive(Clone)]
struct ReadInfo {
    start: u32,
    end: u32,
    hp_tag: Option<u8>,
    // Instead of storing the full record, store pre-computed STR calls for overlapping targets
    str_calls: Vec<(u32, u32, Call)>, // (target_start, target_end, str_call)
}

impl ReadInfo {
    fn from_record_with_targets(
        record: rust_htslib::bam::Record,
        minlen: u32,
        target_intervals: &[(u32, u32, usize)], // (start, end, batch_index)
        batch: &[RepeatInterval],               // Access to original repeat coordinates
    ) -> Self {
        // Extract HP tag if present
        let hp_tag = record.aux(b"HP").ok().and_then(|aux| match aux {
            rust_htslib::bam::record::Aux::U8(val) => Some(val),
            rust_htslib::bam::record::Aux::I32(val) if (0..=255).contains(&val) => Some(val as u8),
            _ => None,
        });

        let read_start = record.reference_start() as u32;
        let read_end = record.reference_end() as u32;

        // Pre-compute STR calls for all overlapping targets
        let mut str_calls = Vec::new();
        for &(target_start, target_end, batch_idx) in target_intervals {
            if read_start < target_end && read_end > target_start {
                // Get the original repeat coordinates (without padding)
                let original_repeat = &batch[batch_idx];
                let original_start = original_repeat.start;
                let original_end = original_repeat.end;

                // Calculate STR call using the original coordinates + analysis padding
                let analysis_start = original_start.saturating_sub(10);
                let analysis_end = original_end + 10;
                let str_call = Self::calculate_str_call_for_region_static(
                    &record,
                    minlen,
                    analysis_start,
                    analysis_end,
                );

                // Store with original coordinates for later matching
                str_calls.push((analysis_start, analysis_end, str_call));
            }
        }

        ReadInfo { start: read_start, end: read_end, hp_tag, str_calls }
    }

    /// Get pre-computed STR call for a specific target region
    #[inline]
    fn get_str_call_for_region(&self, target_start: u32, target_end: u32) -> Option<Call> {
        self.str_calls
            .iter()
            .find(|(start, end, _)| *start == target_start && *end == target_end)
            .map(|(_, _, call)| *call)
    }

    /// Static method to calculate STR call from BAM record (moved from instance method)
    fn calculate_str_call_for_region_static(
        record: &rust_htslib::bam::Record,
        minlen: u32,
        start: u32,
        end: u32,
    ) -> Call {
        let mut call: i64 = 0;
        let mut reference_position = record.reference_start() as u32;
        let mut clipped = false;

        for entry in record.cigar().iter() {
            match entry {
                rust_htslib::bam::record::Cigar::Match(len)
                | rust_htslib::bam::record::Cigar::Equal(len)
                | rust_htslib::bam::record::Cigar::Diff(len) => {
                    reference_position += *len;
                }
                rust_htslib::bam::record::Cigar::Del(len) => {
                    if *len > minlen && start < reference_position && reference_position < end {
                        call -= i64::from(*len);
                    }
                    reference_position += *len;
                }
                rust_htslib::bam::record::Cigar::SoftClip(len) => {
                    if !is_accidental_2d(record)
                        && *len > minlen
                        && start < reference_position
                        && reference_position < end
                    {
                        call += i64::from(*len);
                        clipped = true;
                    }
                }
                rust_htslib::bam::record::Cigar::Ins(len) => {
                    if *len > minlen && start < reference_position && reference_position < end {
                        call += i64::from(*len);
                    }
                }
                rust_htslib::bam::record::Cigar::RefSkip(len) => reference_position += *len,
                _ => (),
            }
        }

        if clipped {
            Call::Clip(call)
        } else {
            Call::Span(call)
        }
    }
}

/// Process a single target using lightweight read information
fn process_target_from_read_info(
    read_infos: &[ReadInfo],
    repeat: &RepeatInterval,
    genotype: &crate::call::GenotypeConfig,
) -> Result<(Genotype, bool), String> {
    let start_ext = repeat.start.saturating_sub(10);
    let end_ext = repeat.end + 10;

    if genotype.unphased {
        let mut calls = Vec::with_capacity(50);

        for read_info in read_infos {
            // Check if this read overlaps with our target region
            if read_info.start > end_ext || read_info.end < start_ext {
                continue;
            }

            // Get pre-computed STR call for this target region
            if let Some(str_call) = read_info.get_str_call_for_region(start_ext, end_ext) {
                calls.push(str_call);
            }
        }

        if calls.len() < genotype.support {
            return Err(format!("Insufficient support: {} < {}", calls.len(), genotype.support));
        }

        calls.sort_unstable_by_key(|call| call.value());
        let (hap1, hap2) = calls.split_at(calls.len() / 2);

        Ok((
            Genotype {
                repeat: repeat.clone(),
                phase1: median_str_length(hap1, genotype.support),
                phase2: median_str_length(hap2, genotype.support),
            },
            false,
        )) // unphased mode never has HP tags
    } else {
        // Phased processing
        let mut phase1_calls = Vec::with_capacity(25);
        let mut phase2_calls = Vec::with_capacity(25);
        let mut unphased_calls = Vec::with_capacity(10);
        let mut found_hp_tags = false;

        for read_info in read_infos {
            // Check if this read overlaps with our target region
            if read_info.start > end_ext || read_info.end < start_ext {
                continue;
            }

            // Check for phasing information
            if let Some(hp_tag) = read_info.hp_tag {
                found_hp_tags = true;

                // Get pre-computed STR call for this target region
                if let Some(str_call) = read_info.get_str_call_for_region(start_ext, end_ext) {
                    match hp_tag {
                        1 => phase1_calls.push(str_call),
                        2 => phase2_calls.push(str_call),
                        _ => unphased_calls.push(str_call),
                    }
                }
            } else {
                // Get pre-computed STR call for this target region
                if let Some(str_call) = read_info.get_str_call_for_region(start_ext, end_ext) {
                    unphased_calls.push(str_call);
                }
            }
        }

        let total_reads = phase1_calls.len() + phase2_calls.len() + unphased_calls.len();
        if total_reads < genotype.support {
            return Err(format!("Insufficient support: {total_reads} < {}", genotype.support));
        }

        phase1_calls.sort_unstable_by_key(|call| call.value());
        phase2_calls.sort_unstable_by_key(|call| call.value());
        unphased_calls.sort_unstable_by_key(|call| call.value());

        let is_phased = !phase1_calls.is_empty() && !phase2_calls.is_empty();

        if is_phased {
            Ok((
                Genotype {
                    repeat: repeat.clone(),
                    phase1: median_str_length(&phase1_calls, genotype.support),
                    phase2: median_str_length(&phase2_calls, genotype.support),
                },
                found_hp_tags,
            ))
        } else {
            // If no phased reads found, should return NaN (not split unphased reads)
            if phase1_calls.is_empty() && phase2_calls.is_empty() {
                Ok((
                    Genotype { repeat: repeat.clone(), phase1: f64::NAN, phase2: f64::NAN },
                    found_hp_tags,
                ))
            } else {
                // Fall back to unphased processing only if we have some phased reads but insufficient in one phase
                let mut all_calls = phase1_calls;
                all_calls.extend(phase2_calls);
                all_calls.extend(unphased_calls);
                all_calls.sort_unstable_by_key(|call| call.value());

                let (hap1, hap2) = all_calls.split_at(all_calls.len() / 2);
                Ok((
                    Genotype {
                        repeat: repeat.clone(),
                        phase1: median_str_length(hap1, genotype.support),
                        phase2: median_str_length(hap2, genotype.support),
                    },
                    found_hp_tags,
                ))
            }
        }
    }
}
