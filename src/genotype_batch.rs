use log::{error, warn};
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;
use std::collections::HashMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::bam_utils::{get_bam_reader, get_bam_reader_with_validation, is_accidental_2d};
use crate::batching::Batch;
use crate::call::Genotype;
use crate::errors::{InquiSTRError, InquiSTRResult};
use crate::repeats::RepeatInterval;

// Global mutex to serialize BAM/CRAM reader creation and fetch operations
// This prevents concurrent htslib index operations that can cause segfaults in cram_index_free
static BAM_READER_LOCK: Mutex<()> = Mutex::new(());

// Track whether phasing validation has been performed (only needed once, on first batch)
static PHASING_VALIDATED: AtomicBool = AtomicBool::new(false);

/// Represents a STR call from a single read
#[derive(Clone, Copy)] // Make it Copy for better performance
pub enum CallType {
    Span(i64),
    Clip(i64),
}

pub enum Phase {
    Phase1,
    Phase2,
}

impl Phase {
    pub fn from_u8(value: u8) -> Option<Phase> {
        match value {
            1 => Some(Phase::Phase1),
            2 => Some(Phase::Phase2),
            _ => None,
        }
    }
}

pub struct Call {
    value: CallType,
    hp_tag: Option<Phase>,
}

impl Call {
    /// Extract the numeric STR length difference regardless of call type
    ///
    /// Returns the signed length difference from reference for both
    /// spanning reads (Call::Span) and soft-clipped reads (Call::Clip).
    #[inline]
    pub fn value(&self) -> i64 {
        match self.value {
            CallType::Span(v) | CallType::Clip(v) => v,
        }
    }
}

/// Calculate the median STR length difference from reference genome
///
/// Spanning reads (Call::Span) are preferred over soft-clipped reads (Call::Clip)
/// for more accurate genotyping. The algorithm:
///
/// 1. If fewer than `support` calls, return NaN (insufficient data)
/// 2. Separate calls into spanning and clipped groups
/// 3. If ≥ `support` spanning reads exist, use only those for median
/// 4. If `require_spanning` is set and spanning reads < `support`, return NaN
/// 5. Otherwise, use all spanning reads + largest soft-clips to reach `support` threshold
/// 6. Return median of the selected reads
///
/// # Arguments
/// * `array` - Array of STR calls (sorted or unsorted)
/// * `support` - Minimum number of reads required
/// * `require_spanning` - If true, return NaN when fewer than `support` spanning reads are available
///
/// # Returns
/// Median STR length difference from reference, or NaN if insufficient support
pub fn median_str_length(array: &[&Call], support: usize, require_spanning: bool) -> f64 {
    if array.len() < support {
        return f64::NAN;
    }

    // Separate spanning and clipped reads
    // Pre-size to actual array length since we'll likely use most/all of them
    let mut spanning = Vec::with_capacity(array.len());
    let mut clipped = Vec::with_capacity(array.len());

    for call in array {
        match call.value {
            CallType::Span(v) => spanning.push(v),
            CallType::Clip(v) => clipped.push(v),
        }
    }

    // Use spanning reads if we have enough, otherwise supplement with largest clips
    let mut values = if spanning.len() >= support {
        spanning
    } else if require_spanning {
        return f64::NAN;
    } else {
        // Sort clipped from large to small to get the largest clips
        clipped.sort_unstable_by(|a, b| b.cmp(a));
        spanning.extend_from_slice(&clipped[0..(support - spanning.len()).min(clipped.len())]);
        spanning
    };

    values.sort_unstable();

    // Calculate median
    if (values.len() % 2) == 0 {
        let ind_left = values.len() / 2 - 1;
        let ind_right = values.len() / 2;
        (values[ind_left] + values[ind_right]) as f64 / 2.0
    } else {
        values[values.len() / 2] as f64
    }
}

fn from_record_with_target(
    record: &rust_htslib::bam::Record,
    minlen: u32,
    target_interval: &(u32, u32), // (start, end)
    no_extend: bool,
) -> Call {
    // Pre-compute STR calls for all overlapping targets
    let &(target_start, target_end) = target_interval;

    // Calculate STR call using the original coordinates, optionally with 10bp analysis padding
    let (analysis_start, analysis_end) = if no_extend {
        (target_start, target_end)
    } else {
        (target_start.saturating_sub(10), target_end + 10)
    };
    calculate_str_call_for_region(record, minlen, analysis_start, analysis_end)
}

/// Retrieve pre-computed STR call for a specific target region
///
/// Looks up the STR call that was computed during `from_record_with_targets()`
/// for the given target coordinates. Returns `None` if this read doesn't overlap
/// the specified target (overlap was checked during ReadInfo creation).
///
/// # Arguments
/// * `target_start` - Start coordinate of the target region (with 10bp padding)
/// * `target_end` - End coordinate of the target region (with 10bp padding)
///
/// # Returns
/// The pre-computed STR call if this read overlaps the target, None otherwise
/// Calculate STR call from BAM record by analyzing CIGAR string
///
/// Walks through the CIGAR string to identify insertions, deletions, and soft-clips
/// within the specified target region that exceed the minimum length threshold.
///
/// The algorithm:
/// - Insertions (I) and soft-clips (S) add to the call value (expansions)
/// - Deletions (D) subtract from the call value (contractions)
/// - Only events occurring within the target region and exceeding minlen are counted
/// - Returns `Call::Clip` if soft-clips contributed, `Call::Span` otherwise
///
/// # Arguments
/// * `record` - BAM record containing CIGAR string
/// * `minlen` - Minimum length for indels/clips to be counted
/// * `start` - Start of target region
/// * `end` - End of target region
///
/// # Returns
/// Call::Span or Call::Clip with the computed STR length difference
fn calculate_str_call_for_region(
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

    let hp_tag = record.aux(b"HP").ok().and_then(|aux| match aux {
        rust_htslib::bam::record::Aux::U8(val) => Some(Phase::from_u8(val)?),
        rust_htslib::bam::record::Aux::I32(val) if (0..=255).contains(&val) => {
            Some(Phase::from_u8(val as u8)?)
        }
        _ => None,
    });

    if clipped {
        Call { value: CallType::Clip(call), hp_tag }
    } else {
        Call { value: CallType::Span(call), hp_tag }
    }
}

fn get_str_call_for_region(
    calls: &HashMap<(u32, u32), Vec<Call>>,
    target_start: u32,
    target_end: u32,
) -> Option<&Vec<Call>> {
    calls.get(&(target_start, target_end))
}

/// Process a single STR target using pre-computed read information
///
/// Uses the lightweight ReadInfo structures (which contain pre-computed STR calls)
/// to genotype a single STR target. This avoids re-parsing BAM records and CIGAR strings.
///
/// Depending on the genotype configuration:
/// - **Unphased mode**: Splits all reads into two groups and returns median of each
/// - **Phased mode**: Separates reads by HP tag (1 vs 2) and returns median per haplotype
///   - Falls back to unphased if insufficient reads in either phase
///   - Returns NaN if no phased reads found at all
///
/// # Arguments
/// * `read_infos` - Pre-computed read information for all reads in the batch
/// * `repeat` - The STR target to genotype
/// * `genotype` - Configuration (minlen, support threshold, phasing)
///
/// # Returns
/// Tuple of (Genotype, found_hp_tags) or error if insufficient read support
fn process_target_from_read_info(
    batch_calls: &mut HashMap<(u32, u32), Vec<Call>>,
    repeat: &RepeatInterval,
    genotype: &crate::call::GenotypeConfig,
    chromosome: &str,
) -> Result<(Genotype, bool), String> {
    // Haploid path: pool all reads regardless of HP tags, produce single allele + NaN
    if genotype
        .haploid
        .as_ref()
        .is_some_and(|chroms| chroms.iter().any(|c| c == chromosome))
    {
        let repeat_calls = match batch_calls.get(&(repeat.start, repeat.end)) {
            Some(calls) => calls,
            None => return Err(format!("Insufficient support: {} < {}", 0, genotype.support)),
        };

        if repeat_calls.len() < genotype.support {
            return Err(format!(
                "Insufficient support: {} < {}",
                repeat_calls.len(),
                genotype.support
            ));
        }

        let all_refs: Vec<&Call> = repeat_calls.iter().collect();
        let median = median_str_length(&all_refs, genotype.support, genotype.require_spanning);
        return Ok((Genotype { repeat: repeat.clone(), phase1: median, phase2: f64::NAN }, false));
    }

    if genotype.unphased {
        let repeat_calls: &mut Vec<Call> = match batch_calls.get_mut(&(repeat.start, repeat.end)) {
            Some(calls) => calls,
            None => return Err(format!("Insufficient support: {} < {}", 0, genotype.support)),
        };

        if repeat_calls.len() < genotype.support {
            return Err(format!(
                "Insufficient support: {} < {}",
                repeat_calls.len(),
                genotype.support
            ));
        }

        repeat_calls.sort_unstable_by_key(|call| call.value());
        let top_fraction = genotype.imbalance.unwrap_or(0.5);
        let cut: usize = ((repeat_calls.len() as f64) * (1.0 - top_fraction)) as usize;
        // Clamp to valid range to avoid empty slices
        let cut = cut.clamp(1, repeat_calls.len() - 1);
        let (hap1, hap2) = repeat_calls.split_at(cut);

        Ok((
            Genotype {
                repeat: repeat.clone(),
                // TODO: code below looks awful
                phase1: median_str_length(
                    &hap1.iter().collect::<Vec<_>>(),
                    genotype.support,
                    genotype.require_spanning,
                ),
                phase2: median_str_length(
                    &hap2.iter().collect::<Vec<_>>(),
                    genotype.support,
                    genotype.require_spanning,
                ),
            },
            false,
        )) // unphased mode never has HP tags
    } else {
        // Phased processing
        let mut phase1_calls = Vec::with_capacity(25);
        let mut phase2_calls = Vec::with_capacity(25);
        let mut unphased_calls = Vec::with_capacity(10);
        let mut found_hp_tags = false;

        let repeat_calls = match get_str_call_for_region(batch_calls, repeat.start, repeat.end) {
            Some(calls) => calls,
            None => return Err(format!("Insufficient support: {} < {}", 0, genotype.support)),
        };

        for c in repeat_calls.iter() {
            // Check for phasing information
            match c.hp_tag {
                Some(Phase::Phase1) => phase1_calls.push(c),
                Some(Phase::Phase2) => phase2_calls.push(c),
                None => unphased_calls.push(c),
            }
            if c.hp_tag.is_some() {
                found_hp_tags = true;
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
                    phase1: median_str_length(
                        &phase1_calls,
                        genotype.support,
                        genotype.require_spanning,
                    ),
                    phase2: median_str_length(
                        &phase2_calls,
                        genotype.support,
                        genotype.require_spanning,
                    ),
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
                        phase1: median_str_length(
                            hap1,
                            genotype.support,
                            genotype.require_spanning,
                        ),
                        phase2: median_str_length(
                            hap2,
                            genotype.support,
                            genotype.require_spanning,
                        ),
                    },
                    found_hp_tags,
                ))
            }
        }
    }
}

/// Process a batch with an existing BAM reader (more efficient, avoids reopening file)
///
/// This is the preferred method for processing multiple batches from the same file,
/// as it reuses the BAM reader and avoids file descriptor exhaustion.
///
/// The function:
/// 1. Fetches all reads in the batch region from the BAM file
/// 2. Filters reads by quality (MAPQ > 10) and overlap with STR targets
/// 3. Pre-computes STR calls for all overlapping targets during read processing
/// 4. Processes each target using the pre-computed read information
///
/// # Arguments
/// * `batch` - The batch of STR targets to process
/// * `bam` - Mutable reference to an indexed BAM/CRAM reader
/// * `genotype` - Configuration for genotyping (minlen, support, phasing)
///
/// # Returns
/// Vector of genotypes for each target in the batch (NaN for failed targets)
pub fn process_batch_with_reader(
    batch: &Batch,
    bam: &mut rust_htslib::bam::IndexedReader,
    genotype: &crate::call::GenotypeConfig,
) -> InquiSTRResult<Vec<Genotype>> {
    let mut results = Vec::new();

    // Serialize both header access and fetch for CRAM stability
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

    // Create a hashmap with keys as (start, end) for quick lookup, and value is a vector of Calls
    let mut calls_map = std::collections::HashMap::new();

    for record_result in bam.rc_records() {
        match record_result {
            Ok(record) => {
                // Skip low quality reads early
                if record.mapq() <= 10 {
                    continue;
                }
                let read_start = record.reference_start() as u32;
                let read_end = record.reference_end() as u32;

                // Check if read overlaps with a target interval
                for repeat in &batch.repeats {
                    // Repeats are sorted by start position, so we can break early
                    if repeat.start >= read_end {
                        break;
                    }
                    if read_start < repeat.end && read_end > repeat.start {
                        let call = from_record_with_target(
                            &record,
                            genotype.minlen,
                            &(repeat.start, repeat.end),
                            genotype.no_extend,
                        );

                        calls_map
                            .entry((repeat.start, repeat.end))
                            .or_insert_with(Vec::new)
                            .push(call);
                    }
                }
            }
            Err(e) => {
                let error_str = e.to_string();
                if error_str.contains("CRC32 failure") {
                    error!(
                        "CRAM decoding error (CRC32 mismatch) in batch {}: {}. The reference genome likely doesn't match the CRAM file.",
                        batch.chromosome, error_str
                    );
                    return Err(InquiSTRError::new(format!(
                        "CRAM decoding error (CRC32 mismatch) in batch {}: {}. The reference genome likely doesn't match the CRAM file.",
                        batch.chromosome, error_str
                    )));
                } else if error_str.contains("truncated record") {
                    error!(
                        "Truncated read in batch {}: {}. This is usually a transient network interruption for remote files.",
                        batch.chromosome, error_str
                    );
                    return Err(InquiSTRError::new(format!(
                        "truncated record in batch {}: {}",
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
        let result =
            process_target_from_read_info(&mut calls_map, repeat, genotype, &batch.chromosome);

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

/// Process a batch using a thread-local BAM reader (no global lock)
///
/// This version is designed for use with per-thread BAM readers where no
/// synchronization is needed because each thread has exclusive access to its reader.
/// This enables true parallel processing without lock contention.
///
/// Each reader is independent with its own file handle, so concurrent fetch()
/// operations are safe. Reader creation is still serialized in bam_pool.rs.
///
/// **IMPORTANT**: Only use this when you can guarantee the reader is not shared
/// between threads (e.g., thread-local storage or indexed reader pool).
///
/// # Arguments
/// * `batch` - The batch of STR targets to process
/// * `bam` - Mutable reference to an indexed BAM/CRAM reader (thread-local)
/// * `genotype` - Configuration for genotyping (minlen, support, phasing)
///
/// # Returns
/// Vector of genotypes for each target in the batch (NaN for failed targets)
pub fn process_batch_with_dedicated_reader(
    batch: &Batch,
    bam: &mut rust_htslib::bam::IndexedReader,
    genotype: &crate::call::GenotypeConfig,
) -> InquiSTRResult<Vec<Genotype>> {
    let mut results = Vec::new();

    // No lock - each thread has its own independent reader
    // Concurrent fetch() calls on different readers are safe
    let tid = match bam.header().tid(batch.chromosome.as_bytes()) {
        Some(t) => t,
        None => {
            warn!("Chromosome {} not found in BAM header", batch.chromosome);
            return Ok(results);
        }
    };

    if let Err(e) = bam.fetch((tid, batch.start, batch.end)) {
        warn!(
            "Warning: Failed to fetch region {}:{}-{}: {}",
            batch.chromosome, batch.start, batch.end, e
        );
        return Ok(results);
    }

    // Create a hashmap with keys as (start, end) for quick lookup, and value is a vector of Calls
    let mut calls_map = std::collections::HashMap::new();

    for record_result in bam.rc_records() {
        match record_result {
            Ok(record) => {
                // Skip low quality reads early
                if record.mapq() <= 10 {
                    continue;
                }
                let read_start = record.reference_start() as u32;
                let read_end = record.reference_end() as u32;

                // Check if read overlaps with a target interval
                for repeat in &batch.repeats {
                    // Repeats are sorted by start position, so we can break early
                    if repeat.start >= read_end {
                        break;
                    }
                    if read_start < repeat.end && read_end > repeat.start {
                        let call = from_record_with_target(
                            &record,
                            genotype.minlen,
                            &(repeat.start, repeat.end),
                            genotype.no_extend,
                        );

                        calls_map
                            .entry((repeat.start, repeat.end))
                            .or_insert_with(Vec::new)
                            .push(call);
                    }
                }
            }
            Err(e) => {
                let error_str = e.to_string();
                if error_str.contains("CRC32 failure") {
                    error!(
                        "CRAM decoding error (CRC32 mismatch) in batch {}: {}. The reference genome likely doesn't match the CRAM file.",
                        batch.chromosome, error_str
                    );
                    return Err(InquiSTRError::new(format!(
                        "CRAM decoding error (CRC32 mismatch) in batch {}: {}. The reference genome likely doesn't match the CRAM file.",
                        batch.chromosome, error_str
                    )));
                } else if error_str.contains("truncated record") {
                    error!(
                        "Truncated read in batch {}: {}. This is usually a transient network interruption for remote files.",
                        batch.chromosome, error_str
                    );
                    return Err(InquiSTRError::new(format!(
                        "truncated record in batch {}: {}",
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
        let result =
            process_target_from_read_info(&mut calls_map, repeat, genotype, &batch.chromosome);

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

/// Unified batch processor that handles a single batch and returns results
///
/// **NOTE**: This creates a new BAM reader for each batch. For processing many batches,
/// use `process_batch_with_reader()` instead to reuse a single reader and avoid
/// file descriptor exhaustion.
///
/// This function is thread-safe and serializes BAM reader creation to prevent
/// concurrent htslib index operations that can cause segfaults. On the first batch
/// (if phasing is enabled), it also validates that HP tags are present in the BAM.
///
/// # Arguments
/// * `batch` - The batch of STR targets to process
/// * `bamp` - Path to the BAM/CRAM file
/// * `reference` - Optional path to reference genome (required for CRAM)
/// * `genotype` - Configuration for genotyping (minlen, support, phasing)
///
/// # Returns
/// Vector of genotypes for each target in the batch (NaN for failed targets)
pub fn process_batch_worker(
    batch: &Batch,
    bamp: &String,
    reference: &Option<String>,
    genotype: &crate::call::GenotypeConfig,
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

    // Use the new efficient function
    let results = process_batch_with_reader(batch, &mut bam, genotype)?;

    // Explicitly drop the BAM reader to free file descriptors immediately
    drop(bam);

    Ok(results)
}
