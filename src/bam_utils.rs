use hts_sys;
use log::{debug, error, warn};

use crate::errors::{InquiSTRError, InquiSTRResult};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Aux;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::env;
use std::path::PathBuf;
use url::Url;

/// Get chromosome lengths; for CRAM prefer .fai from the provided reference to avoid early CRAM opens
pub fn get_chrom_lengths_from_bam_header(
    bam: String,
    reference: &Option<String>,
) -> InquiSTRResult<HashMap<String, u64>> {
    // Prefer FAI for CRAM to avoid opening CRAM before main processing
    if bam.ends_with(".cram") {
        if let Some(ref_path) = reference {
            if let Some(map) = read_fai_lengths(ref_path)
                && !map.is_empty()
            {
                debug!(
                    "get_chrom_lengths_from_bam_header: using reference FAI for chrom lengths ({} entries)",
                    map.len()
                );
                return Ok(map);
            }
            warn!(
                "Reference FAI not found or empty for '{}', falling back to CRAM header read",
                ref_path
            );
        } else {
            warn!("CRAM input without --reference provided; falling back to CRAM header read");
        }
    }

    debug!("get_chrom_lengths_from_bam_header: opening file to read header");
    let reader = get_bam_reader(&bam, reference)?;
    let header = bam::Header::from_template(reader.header());

    let mut chrom_lengts = HashMap::new();
    for (key, records) in header.to_hashmap() {
        if key != "SQ" {
            continue;
        }
        for record in records {
            chrom_lengts.insert(
                record["SN"].clone(),
                record["LN"]
                    .parse()
                    .expect("Failed to parse length of chromosome"),
            );
        }
    }

    debug!(
        "get_chrom_lengths_from_bam_header: completed successfully, found {} chromosomes",
        chrom_lengts.len()
    );
    Ok(chrom_lengts)
}

/// Read chromosome lengths from a .fai next to the reference FASTA
fn read_fai_lengths(reference_fasta: &str) -> Option<HashMap<String, u64>> {
    let fai_path = if reference_fasta.ends_with(".fai") {
        PathBuf::from(reference_fasta)
    } else {
        PathBuf::from(format!("{}.fai", reference_fasta))
    };

    if !fai_path.exists() {
        debug!("FAI file not found: {}", fai_path.display());
        return None;
    }

    let mut map = HashMap::new();
    let content = match std::fs::read_to_string(&fai_path) {
        Ok(s) => s,
        Err(e) => {
            warn!("Failed to read FAI file {}: {}", fai_path.display(), e);
            return None;
        }
    };
    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        if let (Some(name), Some(len_str)) = (fields.next(), fields.next())
            && let Ok(len) = len_str.parse::<u64>()
        {
            map.insert(name.to_string(), len);
        }
    }
    Some(map)
}

/// Sets up CRAM reference caching for remote CRAM files
/// Configures HTS_REF_CACHE to cache reference sequences in ~/.cache/inquistr/
/// Only sets up caching when accessing remote CRAM files since local files don't need it
pub fn setup_reference_caching(bam_path: &str) {
    // Only set up caching for remote CRAM files
    // Local CRAM files with local references don't benefit from caching
    let is_remote = bam_path.starts_with("http://")
        || bam_path.starts_with("https://")
        || bam_path.starts_with("ftp://")
        || bam_path.starts_with("s3://");

    let is_cram = bam_path.ends_with(".cram");

    if !is_remote || !is_cram {
        return;
    }

    // Set up cache directory for CRAM reference sequences
    let cache_dir = if let Ok(home) = env::var("HOME") {
        PathBuf::from(home).join(".cache").join("inquistr")
    } else {
        PathBuf::from(".inquistr_cache")
    };

    // Create the cache directory if it doesn't exist
    if let Err(e) = std::fs::create_dir_all(&cache_dir) {
        warn!("Failed to create cache directory: {}", e);
        return;
    }

    let cache_path = cache_dir.to_string_lossy().to_string();

    // Set htslib caching environment variables for CRAM reference sequences
    // SAFETY: These set_var calls are safe because:
    // 1. They occur during initialization before any parallel processing
    // 2. These variables are only read by htslib, not by other Rust threads
    // 3. We check if they're already set to avoid unnecessary modifications
    if env::var("HTS_REF_CACHE").is_err() {
        unsafe {
            env::set_var("HTS_REF_CACHE", &cache_path);
        }
        debug!("Set HTS_REF_CACHE to: {} (for remote CRAM reference caching)", cache_path);
    }

    if env::var("REF_CACHE").is_err() {
        unsafe {
            env::set_var("REF_CACHE", &cache_path);
        }
    }
}

/// Clean up downloaded index files from the current working directory
/// htslib downloads .bai/.crai files when accessing remote BAM/CRAM files
/// This function removes those files to avoid cluttering the working directory
pub fn cleanup_index_files(bam_path: &str) {
    // Only clean up if the path is a remote URL
    if let Ok(url) = Url::parse(bam_path)
        && ["http", "https", "ftp", "s3"].contains(&url.scheme())
    {
        // Try common index extensions that htslib might download
        let extensions = [".bai", ".crai"];

        for ext in &extensions {
            // Extract filename from URL and construct index filename
            if let Some(filename) = url.path().split('/').next_back() {
                let index_name = format!("{}{}", filename, ext);
                let index_path = PathBuf::from(&index_name);

                if index_path.exists() {
                    if let Err(e) = std::fs::remove_file(&index_path) {
                        debug!("Failed to remove index file {}: {}", index_path.display(), e);
                    } else {
                        debug!("Cleaned up downloaded index file: {}", index_path.display());
                    }
                }
            }
        }
    }
}

/// Sets up the CURL_CA_BUNDLE environment variable for HTTPS/S3 access
/// Tries to use a CA bundle from standard locations, with appropriate fallbacks
fn setup_ssl_certificates() {
    // Only configure if not already set by the user
    if env::var("CURL_CA_BUNDLE").is_ok() {
        return;
    }

    // Common CA bundle locations across different systems
    let possible_paths = vec![
        "/etc/ssl/certs/ca-certificates.crt",     // Debian/Ubuntu
        "/etc/pki/tls/certs/ca-bundle.crt",       // RHEL/CentOS/Amazon Linux
        "/etc/ssl/ca-bundle.pem",                 // SUSE
        "/usr/local/share/certs/ca-root-nss.crt", // FreeBSD
        "/usr/local/etc/openssl/cert.pem",        // macOS Homebrew
        "/etc/ssl/cert.pem",                      // macOS/OpenBSD
    ];

    // Try each path in order
    for path in possible_paths {
        if std::path::Path::new(path).exists() {
            // SAFETY: This set_var call is safe because:
            // 1. It occurs during initialization, before parallel processing
            // 2. The variable is only used by libcurl for SSL certificate validation
            // 3. It's set once and not modified afterward
            unsafe {
                env::set_var("CURL_CA_BUNDLE", path);
            }
            return;
        }
    }

    // None of the paths exist, warn the user
    warn!(
        "Could not find a valid CA certificate bundle for HTTPS/S3 access. \
        HTTPS/S3 connections may fail. Set the CURL_CA_BUNDLE environment \
        variable to the path of your system's CA certificate bundle."
    );
}

/// Create a BAM reader from path or URL with appropriate configuration
/// Note: This returns IndexedReader for compatibility with existing code, but we should
/// consider using Reader for sequential access
pub fn get_bam_reader(
    bamp: &String,
    reference: &Option<String>,
) -> InquiSTRResult<bam::IndexedReader> {
    get_bam_reader_internal(bamp, reference, false)
}

/// Create a BAM reader with phasing validation (for first batch in non-unphased mode)
pub fn get_bam_reader_with_validation(
    bamp: &String,
    reference: &Option<String>,
) -> InquiSTRResult<bam::IndexedReader> {
    get_bam_reader_internal(bamp, reference, true)
}

/// Internal function to create a BAM reader with optional phasing validation
fn get_bam_reader_internal(
    bamp: &String,
    reference: &Option<String>,
    validate_phasing: bool,
) -> InquiSTRResult<bam::IndexedReader> {
    debug!("Opening BAM/CRAM file: {}", bamp);

    // Set up reference caching for remote CRAM files
    setup_reference_caching(bamp);

    let mut bam = if bamp.starts_with("s3")
        || bamp.starts_with("https://")
        || bamp.starts_with("ftp://")
    {
        setup_ssl_certificates();
        bam::IndexedReader::from_url(&Url::parse(bamp.as_str()).expect("Failed to parse URL"))
            .map_err(|err| {
                error!("Error opening remote BAM file: {err}");
                InquiSTRError::new(format!("Error opening remote BAM file: {err}"))
            })?
    } else {
        bam::IndexedReader::from_path(bamp).map_err(|err| {
                let err_msg = err.to_string();

                // Check if this is an index-related error
                if err_msg.contains("index") {
                    // Check if the index file actually exists
                    let possible_indices = [
                        format!("{}.bai", bamp),
                        format!("{}.crai", bamp),
                        bamp.trim_end_matches(".bam").to_string() + ".bai",
                        bamp.trim_end_matches(".cram").to_string() + ".crai",
                    ];

                    let index_exists = possible_indices.iter().any(|idx| PathBuf::from(idx).exists());

                    if index_exists {
                        error!("Error opening local BAM/CRAM file: {err}");
                        error!("Note: An index file was found but appears to be corrupted or incompatible.");
                        error!("Try regenerating the index with: samtools index {}", bamp);
                        InquiSTRError::new(format!(
                            "Error opening local BAM/CRAM file: {err}\n\
                            Note: An index file was found but appears to be corrupted or incompatible.\n\
                            Try regenerating the index with: samtools index {}", 
                            bamp
                        ))
                    } else {
                        error!("Error opening local BAM/CRAM file: {err}");
                        error!("Note: No index file found. Create one with: samtools index {}", bamp);
                        InquiSTRError::new(format!(
                            "Error opening local BAM/CRAM file: {err}\n\
                            Note: No index file found. Create one with: samtools index {}", 
                            bamp
                        ))
                    }
                } else if err_msg.contains("No file descriptors available") 
                    || err_msg.contains("Too many open files") {
                    error!("Error opening local BAM/CRAM file: {err}");
                    error!("Note: System has run out of available file descriptors.");
                    error!("This typically happens with CRAM files when processing too many samples in parallel.");
                    error!("Each CRAM file opens 10-20 file descriptors (file + index + reference + buffers).");
                    error!("Solutions (in order of effectiveness):");
                    error!("  1. Reduce --parallel-samples to 1 (process samples sequentially)");
                    error!("  2. Increase --batch-size to 100-500 (reduces number of BAM reopens per sample)");
                    error!("  3. Increase system limit: ulimit -n 4096");
                    error!("  4. Convert CRAM to BAM (uses fewer file descriptors)");
                    InquiSTRError::new(format!(
                        "Error opening local BAM/CRAM file: {err}\n\
                        Note: System has run out of available file descriptors.\n\
                        This typically happens with CRAM files when processing too many samples in parallel.\n\
                        Each CRAM file opens 10-20 file descriptors (file + index + reference + buffers).\n\
                        Solutions (in order of effectiveness):\n\
                          1. Reduce --parallel-samples to 1 (process samples sequentially)\n\
                          2. Increase --batch-size to 100-500 (reduces number of BAM reopens per sample)\n\
                          3. Increase system limit: ulimit -n 4096\n\
                          4. Convert CRAM to BAM (uses fewer file descriptors)"
                    ))
                } else {
                    error!("Error opening local BAM/CRAM file: {err}");
                    InquiSTRError::new(format!("Error opening local BAM/CRAM file: {err}"))
                }
            })?
    };
    if bamp.ends_with(".cram") {
        debug!("Detected CRAM file, setting CRAM options...");
        bam.set_cram_options(
            hts_sys::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
            hts_sys::sam_fields_SAM_AUX
                | hts_sys::sam_fields_SAM_MAPQ
                | hts_sys::sam_fields_SAM_CIGAR
                | hts_sys::sam_fields_SAM_POS
                | hts_sys::sam_fields_SAM_TLEN,
        )
        .expect("Failed setting cram options");
        debug!("CRAM options set successfully");

        if reference.is_some() {
            let ref_path = reference.as_ref().unwrap();
            debug!("Setting CRAM reference to: {}", ref_path);
            bam.set_reference(ref_path)
                .map_err(|err| {
                    error!("Failed to set CRAM reference '{}': {}. This usually means the reference genome doesn't match the CRAM file.", ref_path, err);
                    InquiSTRError::new(format!("Failed to set CRAM reference '{}': {}. This usually means the reference genome doesn't match the CRAM file.", ref_path, err))
                })?;
            debug!("CRAM reference set successfully");
        } else {
            error!(
                "No reference provided for CRAM file. Use --reference option to specify the reference genome."
            );
            return Err(InquiSTRError::new("No reference provided for CRAM file. Use --reference option to specify the reference genome.".to_string()));
        }
    }

    debug!("Finished setting up CRAM/BAM reader, validate_phasing={}", validate_phasing);

    // If validation is requested, check for phasing in first reads
    if validate_phasing {
        debug!("Performing phasing validation on IndexedReader...");
        debug!("About to call validate_phasing_on_reader");
        if let Err(e) = validate_phasing_on_reader(&mut bam, 10000) {
            error!("ERROR: {}", e);
            return Err(InquiSTRError::new(e));
        }
        debug!("validate_phasing_on_reader completed successfully");
    }

    debug!("Returning BAM reader from get_bam_reader_internal");
    Ok(bam)
}

/// Validate phasing using an existing IndexedReader (avoids opening file twice)
fn validate_phasing_on_reader(
    bam: &mut bam::IndexedReader,
    max_reads: usize,
) -> Result<(), String> {
    debug!(
        "Starting HP tag validation on existing reader, scanning up to {} reads...",
        max_reads
    );

    let mut reads_checked = 0;

    debug!("About to get header from reader");
    // Get the header to find the first chromosome/contig
    let header = bam.header().clone();
    debug!("Header cloned successfully");

    let target_names = header.target_names();
    debug!("Got {} target names from header", target_names.len());

    if target_names.is_empty() {
        return Err("BAM/CRAM file has no reference sequences in header".to_string());
    }

    // Try to fetch from all chromosomes until we find HP tags or reach max_reads
    // No limit on chromosome count - we return immediately when HP tag is found
    debug!("Starting loop through chromosomes");
    for (tid, target_name) in target_names.iter().enumerate() {
        let chrom = String::from_utf8_lossy(target_name);
        debug!("Attempting to fetch reads from chromosome: {} (tid={})", chrom, tid);

        // Fetch reads from this chromosome (entire chromosome - we'll stop after max_reads anyway)
        debug!("Calling bam.fetch for tid={}", tid);
        match bam.fetch((tid as i32, 0, i32::MAX)) {
            Ok(_) => {
                debug!("Fetch successful for {}, iterating records...", chrom);
                let mut chrom_reads = 0;

                debug!("About to call bam.records() iterator");
                // Iterate through fetched records
                for record_result in bam.records() {
                    if reads_checked >= max_reads {
                        debug!("Reached max_reads limit ({}), stopping validation", max_reads);
                        break;
                    }

                    match record_result {
                        Ok(record) => {
                            chrom_reads += 1;

                            // Check for HP tag
                            if record.aux(b"HP").is_ok() {
                                debug!(
                                    "Found HP tag in read {} on {} (total reads checked: {}), validation successful",
                                    chrom_reads,
                                    chrom,
                                    reads_checked + 1
                                );
                                return Ok(());
                            }

                            reads_checked += 1;
                        }
                        Err(e) => {
                            let error_str = e.to_string();
                            if error_str.contains("CRC32 failure")
                                || error_str.contains("truncated record")
                            {
                                return Err(format!(
                                    "CRAM format error: {}. This usually indicates that the reference genome doesn't match the CRAM file.",
                                    error_str
                                ));
                            }
                            warn!("Error reading record during phasing validation: {}", e);
                            reads_checked += 1;
                        }
                    }
                }

                debug!(
                    "Checked {} reads on {}, total so far: {}",
                    chrom_reads, chrom, reads_checked
                );

                // If we've checked enough reads across chromosomes, stop
                if reads_checked >= max_reads {
                    break;
                }
            }
            Err(e) => {
                debug!("Could not fetch from {}: {}", chrom, e);
                continue;
            }
        }
    }

    Err(format!(
        "No phasing information (HP tags) found in the first {} reads. \
        This suggests the BAM/CRAM file lacks phasing information. \
        Use the --unphased option if you want to proceed without phasing.",
        reads_checked
    ))
}

/// Create a sequential BAM reader for reading all records
fn get_sequential_bam_reader(bamp: &String, reference: &Option<String>) -> bam::Reader {
    debug!("Opening BAM/CRAM file for sequential reading: {}", bamp);

    // Set up reference caching for remote CRAM files
    setup_reference_caching(bamp);

    let mut bam =
        if bamp.starts_with("s3") || bamp.starts_with("https://") || bamp.starts_with("ftp://") {
            setup_ssl_certificates();
            bam::Reader::from_url(&Url::parse(bamp.as_str()).expect("Failed to parse s3 URL"))
                .unwrap_or_else(|err| {
                    error!("Error opening remote BAM file: {err}");
                    std::process::exit(1);
                })
        } else {
            bam::Reader::from_path(bamp).unwrap_or_else(|err| {
                error!("Error opening local BAM file: {err}");
                std::process::exit(1);
            })
        };

    if bamp.ends_with(".cram") {
        debug!("Detected CRAM file for sequential reading, setting CRAM options...");
        bam.set_cram_options(
            hts_sys::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
            hts_sys::sam_fields_SAM_AUX
                | hts_sys::sam_fields_SAM_MAPQ
                | hts_sys::sam_fields_SAM_CIGAR
                | hts_sys::sam_fields_SAM_POS
                | hts_sys::sam_fields_SAM_TLEN,
        )
        .expect("Failed setting cram options");
        debug!("CRAM options set successfully for sequential reader");

        if reference.is_some() {
            let ref_path = reference.as_ref().unwrap();
            debug!("Setting CRAM reference for sequential reader to: {}", ref_path);
            bam.set_reference(ref_path)
                .unwrap_or_else(|err| {
                    error!("Failed to set CRAM reference '{}' for sequential reader: {}. This usually means the reference genome doesn't match the CRAM file.", ref_path, err);
                    std::process::exit(1);
                });
            debug!("CRAM reference set successfully for sequential reader");
        } else {
            error!(
                "No reference provided for CRAM file. Use --reference option to specify the reference genome."
            );
            std::process::exit(1);
        }
    }

    bam
}

/// Validates that the BAM/CRAM file contains phasing information (HP tags) by scanning early reads
/// Returns immediately upon finding the first HP tag, or errors after scanning max_reads without finding any
/// This provides fast failure detection before expensive processing begins
///
/// Note: This function detects HP tags regardless of their Aux type (U8, I32, etc.) by simply
/// checking for the presence of the tag. The actual value parsing and validation happens elsewhere.
pub fn validate_phasing_early(
    bam_path: &str,
    reference: &Option<String>,
    max_reads: usize,
) -> Result<(), String> {
    let bam_path_string = bam_path.to_string();

    // Check if file exists (skip check for URLs)
    if !bam_path_string.starts_with("http")
        && !bam_path_string.starts_with("ftp")
        && !bam_path_string.starts_with("s3")
        && !std::path::Path::new(&bam_path_string).exists()
    {
        return Err(format!("BAM/CRAM file does not exist: {}", bam_path_string));
    }

    debug!("Starting phasing validation for file: {}", bam_path);
    debug!("Reference path: {:?}", reference);
    debug!("File exists, attempting to open BAM/CRAM reader...");

    // Wrap reader creation and usage in a scope to ensure it's dropped before
    // the main processing opens an IndexedReader to the same file

    {
        let mut seq_bam = get_sequential_bam_reader(&bam_path_string, reference);
        debug!("Sequential BAM/CRAM reader opened successfully");

        // Try to read the header to verify the file is accessible
        let header = seq_bam.header();
        let header_text = std::str::from_utf8(header.as_bytes()).unwrap_or("Invalid UTF-8");
        debug!("Header read successfully, length: {} bytes", header_text.len());

        let mut reads_checked = 0;
        let records_iterator = seq_bam.records();
        debug!("Starting HP tag validation, scanning up to {} reads...", max_reads);

        let mut result = Err(format!(
            "No phasing information (HP tags) found in the first {} reads. \
            This suggests the BAM/CRAM file lacks phasing information. \
            Use the --unphased option if you want to proceed without phasing.",
            max_reads
        ));

        for record_result in records_iterator {
            if reads_checked >= max_reads {
                break;
            }

            match record_result {
                Ok(record) => {
                    // Check for HP tag - simple and direct
                    if record.aux(b"HP").is_ok() {
                        debug!("Found HP tag in read {}, validation successful", reads_checked + 1);
                        result = Ok(());
                        break;
                    }

                    reads_checked += 1;
                }
                Err(e) => {
                    // Check if this is a CRAM format error indicating incompatible reference
                    let error_str = e.to_string();
                    if error_str.contains("CRC32 failure") || error_str.contains("truncated record")
                    {
                        error!(
                            "CRAM format error detected: {}. This usually indicates that the reference genome doesn't match the CRAM file. Please verify that you're using the correct reference genome.",
                            error_str
                        );
                        std::process::exit(1);
                    }
                    warn!("Error reading record during phasing validation: {}", e);
                    reads_checked += 1;
                }
            }
        }

        result
        // seq_bam is explicitly dropped here when this scope ends
    }
}

/// Checks if read alignment might be an accidental 2D read (ONT artefact)
///
/// This function determines if a read is an accidental 2D read, which means that right after
/// the template strand, the complement strand was also sequenced. This is a common artifact
/// in ONT data. The read will align in two pieces of similar length to the reference genome,
/// with the second piece on the opposite strand. In that case, softclipped fragments should
/// not be considered as they represent the overlap between template and complement.
pub fn is_accidental_2d(record: &bam::Record) -> bool {
    // An entry in the SA tag consists of: rname, POS, strand, CIGAR, mapQ, NM
    let read_strand = if record.is_reverse() { '-' } else { '+' };
    let sa = record.aux(b"SA");

    // If the SA tag is not present, the read has no supplementary alignments
    if sa.is_err() {
        return false;
    }

    let sa_tag = sa.unwrap();
    let sa_tag = match sa_tag {
        Aux::String(s) => s,
        _ => {
            warn!("Unexpected type of SA auxiliary tag: {sa_tag:?}");
            return false;
        }
    };

    // Split the SA tag into its entries, separated by ';', but remove any empty entries
    let sa_entries = sa_tag
        .split(';')
        .filter(|x| !x.is_empty())
        .collect::<Vec<&str>>();

    // While not conclusive, if there are multiple entries in the SA tag,
    // it is likely that the read is not just a 2D read
    if sa_entries.len() > 1 {
        return false;
    }

    // Ensure we have at least one SA entry
    if sa_entries.is_empty() {
        return false;
    }

    let sa_entry = sa_entries[0].split(',').collect::<Vec<&str>>();

    // SA entry format: rname,pos,strand,CIGAR,mapQ,NM - need at least 4 fields
    if sa_entry.len() < 4 {
        warn!("Malformed SA tag entry: insufficient fields");
        return false;
    }

    // Check if the read is on the opposite strand. If it is on the same strand,
    // it is not an accidental 2D read
    if read_strand == sa_entry[2].chars().next().unwrap() {
        return false;
    }

    // Check if the supplementary alignment overlaps with the original alignment
    // If it does overlap, the read could be an accidental 2D read
    // (alternatively, it could indicate an inverted duplication, but that is not of interest to inquiSTR)
    let start = record.reference_start();
    let end = record.reference_end();
    let sa_start = sa_entry[1].parse::<i64>().unwrap();
    let sa_end = sa_start + cigar_to_rlen(sa_entry[3]);

    // Check if the max of the start values is smaller than the min of the end values
    // If that is the case, the two alignments overlap
    if std::cmp::max(start, sa_start) < std::cmp::min(end, sa_end) {
        debug!("Identified read as accidental 2D read");
        return true;
    }

    false
}

/// Convert CIGAR string to reference length consumed
pub fn cigar_to_rlen(cigar: &str) -> i64 {
    let mut rlen: i64 = 0;
    let mut num = String::new();
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            let n = match num.parse::<i64>() {
                Ok(val) => val,
                Err(_) => {
                    error!("Failed to parse CIGAR number: {}", num);
                    std::process::exit(1);
                }
            };
            match c {
                'M' | '=' | 'X' | 'D' | 'N' => rlen = rlen.saturating_add(n),
                _ => (),
            }
            num.clear();
        }
    }
    rlen
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::{bam, bam::Read};

    #[test]
    #[ignore]
    // the test data contains a 2D-candidates_test_set.bam file
    // this one should have reads that are identified as 2D reads
    // this test is ignored because the test data is not included in the repository
    fn test_is_accidental_2d() {
        let test_file = "test-data/2D-candidates_test_set.bam";

        // Skip test if file doesn't exist (not included in repository)
        if !std::path::Path::new(test_file).exists() {
            println!(
                "Skipping test: {} not found (test data not included in repository)",
                test_file
            );
            return;
        }

        let mut bam = bam::Reader::from_path(test_file).unwrap();
        let mut count = 0;
        let mut all_reads = 0;
        for r in bam.records() {
            let r = r.unwrap();
            all_reads += 1;
            if is_accidental_2d(&r) {
                count += 1;
            }
            if all_reads > 100 {
                break;
            }
        }
        println!("Found {count} 2D reads out of {all_reads} reads");
        assert_eq!(count, all_reads);
    }

    #[test]
    fn test_get_chrom_lengths_from_bam_header() {
        let bam = String::from("test-data/small-test.bam");
        let chrom_lengths = get_chrom_lengths_from_bam_header(bam, &None)
            .expect("Failed to get chromosome lengths");
        assert_eq!(chrom_lengths.get("chr7").unwrap(), &159345973);
    }
}
