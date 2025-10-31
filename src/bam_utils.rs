use hts_sys;
use log::{debug, error, warn};

use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Aux;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::env;
use std::path::PathBuf;
use url::Url;

/// Get chromosome lengths from BAM header
pub fn get_chrom_lengths_from_bam_header(
    bam: String,
    reference: &Option<String>,
) -> HashMap<String, u64> {
    let bam = get_bam_reader(&bam, reference);
    let header = bam::Header::from_template(bam.header());
    let mut chrom_lengts = HashMap::new();
    for (key, records) in header.to_hashmap() {
        for record in records {
            if key != "SQ" {
                continue;
            }
            chrom_lengts.insert(
                record["SN"].clone(),
                record["LN"]
                    .parse()
                    .expect("Failed to parse length of chromosome"),
            );
        }
    }

    chrom_lengts
}

/// Check if a local index file exists for the given BAM/CRAM file
/// Returns the path to the local index if it exists
fn find_local_index(bam_path: &str) -> Option<PathBuf> {
    // Try common index extensions
    let extensions = [".bai", ".crai", ".bam.bai", ".cram.crai"];

    for ext in &extensions {
        let index_path = PathBuf::from(format!("{}{}", bam_path, ext));
        if index_path.exists() {
            debug!("Found local index file: {}", index_path.display());
            return Some(index_path);
        }
    }

    None
}

/// Clean up old cached index files to prevent cache bloat
/// Removes files older than the specified number of days
fn cleanup_old_cache_files(cache_dir: &PathBuf, max_age_days: u64) {
    use std::time::{SystemTime, UNIX_EPOCH};

    let Ok(entries) = std::fs::read_dir(cache_dir) else {
        return;
    };

    let current_time = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs();

    let max_age_seconds = max_age_days * 24 * 60 * 60;
    let mut removed_count = 0;

    for entry in entries.flatten() {
        if let Ok(metadata) = entry.metadata() {
            if metadata.is_file() {
                if let Ok(modified) = metadata.modified() {
                    if let Ok(modified_duration) = modified.duration_since(UNIX_EPOCH) {
                        let age_seconds = current_time.saturating_sub(modified_duration.as_secs());

                        if age_seconds > max_age_seconds
                            && std::fs::remove_file(entry.path()).is_ok()
                        {
                            removed_count += 1;
                            debug!("Removed old cache file: {:?}", entry.path());
                        }
                    }
                }
            }
        }
    }

    if removed_count > 0 {
        debug!("Cleaned up {} old cache files", removed_count);
    }
}

/// Sets up caching for remote files (both index files and CRAM reference sequences)
/// Checks if local index exists first, then configures cache directory for reference sequences
/// Also performs periodic cleanup of old cache files (>30 days)
///
/// Note: htslib automatically caches downloaded index files (.bai/.crai) in the current
/// working directory and reuses them on subsequent runs. CRAM reference sequences are
/// cached in ~/.cache/inquistr/ via HTS_REF_CACHE.
pub fn setup_index_caching(bam_path: &str) {
    // Check if caching is explicitly disabled
    if env::var("INQUISTR_NO_CACHE").is_ok() {
        debug!("Caching disabled by INQUISTR_NO_CACHE environment variable");
        return;
    }

    // First check if there's a local index file in the current directory
    if let Some(local_index) = find_local_index(bam_path) {
        debug!("Using existing local index file: {}", local_index.display());
        return;
    }

    // Only set up remote caching if the path is a URL
    if let Ok(url) = Url::parse(bam_path) {
        if ["http", "https", "ftp", "s3"].contains(&url.scheme()) {
            debug!(
                "Remote file detected - index will be downloaded and cached in current directory"
            );

            // Set up cache directory for reference sequences (CRAM files)
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

            // Clean up old cache files (older than 30 days by default)
            let max_age_days = env::var("INQUISTR_CACHE_MAX_AGE_DAYS")
                .ok()
                .and_then(|s| s.parse::<u64>().ok())
                .unwrap_or(30);

            cleanup_old_cache_files(&cache_dir, max_age_days);

            let cache_path = cache_dir.to_string_lossy().to_string();

            // Set htslib caching environment variables for CRAM reference sequences
            if env::var("HTS_REF_CACHE").is_err() {
                env::set_var("HTS_REF_CACHE", &cache_path);
                debug!("Set HTS_REF_CACHE to: {} (for CRAM reference caching)", cache_path);
            }

            if env::var("REF_CACHE").is_err() {
                env::set_var("REF_CACHE", &cache_path);
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
            env::set_var("CURL_CA_BUNDLE", path);
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
pub fn get_bam_reader(bamp: &String, reference: &Option<String>) -> bam::IndexedReader {
    debug!("Opening BAM/CRAM file: {}", bamp);

    // Set up index caching before opening the file
    setup_index_caching(bamp);

    let mut bam =
        if bamp.starts_with("s3") || bamp.starts_with("https://") || bamp.starts_with("ftp://") {
            setup_ssl_certificates();
            bam::IndexedReader::from_url(&Url::parse(bamp.as_str()).expect("Failed to parse URL"))
                .unwrap_or_else(|err| {
                    error!("Error opening remote BAM file: {err}");
                    std::process::exit(1);
                })
        } else {
            bam::IndexedReader::from_path(bamp).unwrap_or_else(|err| {
                error!("Error opening local BAM file: {err}");
                std::process::exit(1);
            })
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
                .unwrap_or_else(|err| {
                    error!("Failed to set CRAM reference '{}': {}. This usually means the reference genome doesn't match the CRAM file.", ref_path, err);
                    std::process::exit(1);
                });
            debug!("CRAM reference set successfully");
        } else {
            error!("No reference provided for CRAM file. Use --reference option to specify the reference genome.");
            std::process::exit(1);
        }
    }

    bam
}

/// Create a sequential BAM reader for reading all records
fn get_sequential_bam_reader(bamp: &String, reference: &Option<String>) -> bam::Reader {
    debug!("Opening BAM/CRAM file for sequential reading: {}", bamp);

    // Set up index caching before opening the file (some readers may still use index)
    setup_index_caching(bamp);

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
            error!("No reference provided for CRAM file. Use --reference option to specify the reference genome.");
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

    let mut seq_bam = get_sequential_bam_reader(&bam_path_string, reference);
    debug!("Sequential BAM/CRAM reader opened successfully");

    // Try to read the header to verify the file is accessible
    let header = seq_bam.header();
    let header_text = std::str::from_utf8(header.as_bytes()).unwrap_or("Invalid UTF-8");
    debug!("Header read successfully, length: {} bytes", header_text.len());

    let mut reads_checked = 0;
    let records_iterator = seq_bam.records();
    debug!("Starting HP tag validation, scanning up to {} reads...", max_reads);

    for record_result in records_iterator {
        if reads_checked >= max_reads {
            break;
        }

        match record_result {
            Ok(record) => {
                // Check for HP tag - simple and direct
                if record.aux(b"HP").is_ok() {
                    debug!("Found HP tag in read {}, validation successful", reads_checked + 1);
                    return Ok(());
                }

                reads_checked += 1;
            }
            Err(e) => {
                // Check if this is a CRAM format error indicating incompatible reference
                let error_str = e.to_string();
                if error_str.contains("CRC32 failure") || error_str.contains("truncated record") {
                    error!("CRAM format error detected: {}. This usually indicates that the reference genome doesn't match the CRAM file. Please verify that you're using the correct reference genome.", error_str);
                    std::process::exit(1);
                }
                warn!("Error reading record during phasing validation: {}", e);
                reads_checked += 1;
            }
        }
    }

    // No HP tags found in the first max_reads reads
    Err(format!(
        "No phasing information (HP tags) found in the first {} reads. \
        This suggests the BAM/CRAM file lacks phasing information. \
        Use the --unphased option if you want to proceed without phasing.",
        max_reads
    ))
}

/// Checks if read alignment might be an accidental 2D read (ONT artefact)
pub fn is_accidental_2d(record: &bam::Record) -> bool {
    // this function will determine if a read is an accidental 2D read
    // this means that right after the template strand also the complement strand was sequenced
    // this is a common artifact in ONT data
    // the read will then align in two pieces of similar length to the reference genome, with the second piece on the opposite strand
    // in that case, softclipped fragments are not to be considered
    // An entry in the SA tag consist of rname, POS, strand, CIGAR, mapQ, NM
    let read_strand = if record.is_reverse() { '-' } else { '+' };
    let sa = record.aux(b"SA");
    // if the SA tag is not present, the read has no supplementary alignments and is thus not an accidental 2D read
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
    // split the SA tag into its entries, separated by ';', but remove any empty entries
    let sa_entries = sa_tag
        .split(';')
        .filter(|x| !x.is_empty())
        .collect::<Vec<&str>>();
    // while not conclusive, if there are multiple entries in the SA tag, it is likely that the read is not just a 2D read
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
    // check if the read is on the opposite strand. If it is on the same strand, it is not an accidental 2D read
    if read_strand == sa_entry[2].chars().next().unwrap() {
        return false;
    }
    // check if the supplementary alignment overlaps with the original alignment
    // if it does overlap the read could be an accidental 2D read
    // alternatively, it could indicate an inverted duplication
    // but that is not of interst to inquiSTR
    let start = record.reference_start();
    let end = record.reference_end();
    let sa_start = sa_entry[1].parse::<i64>().unwrap();
    let sa_end = sa_start + cigar_to_rlen(sa_entry[3]);
    // check if the max of the start values is smaller than the min of the end values
    // if that is the case, the two alignments overlap
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

/// Clean the index cache directory
pub fn clean_cache(dry_run: bool, all: bool, max_age_days: u64) {
    use std::time::{SystemTime, UNIX_EPOCH};

    // Determine cache directory location
    let cache_dir = if let Ok(home) = env::var("HOME") {
        PathBuf::from(home).join(".cache").join("inquistr")
    } else {
        PathBuf::from(".inquistr_cache")
    };

    if !cache_dir.exists() {
        println!("Cache directory does not exist: {}", cache_dir.display());
        println!("Nothing to clean.");
        return;
    }

    println!("Cleaning cache directory: {}", cache_dir.display());

    let Ok(entries) = std::fs::read_dir(&cache_dir) else {
        eprintln!("Failed to read cache directory");
        return;
    };

    let current_time = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs();

    let max_age_seconds = max_age_days * 24 * 60 * 60;
    let mut total_size = 0u64;
    let mut total_files = 0;
    let mut removed_count = 0;
    let mut removed_size = 0u64;

    for entry in entries.flatten() {
        if let Ok(metadata) = entry.metadata() {
            if metadata.is_file() {
                let file_size = metadata.len();
                total_size += file_size;
                total_files += 1;

                let should_delete = if all {
                    true
                } else if let Ok(modified) = metadata.modified() {
                    if let Ok(modified_duration) = modified.duration_since(UNIX_EPOCH) {
                        let age_seconds = current_time.saturating_sub(modified_duration.as_secs());
                        age_seconds > max_age_seconds
                    } else {
                        false
                    }
                } else {
                    false
                };

                if should_delete {
                    if dry_run {
                        println!(
                            "  [DRY RUN] Would delete: {} ({} bytes)",
                            entry.path().display(),
                            file_size
                        );
                    } else if std::fs::remove_file(entry.path()).is_ok() {
                        println!("  Deleted: {} ({} bytes)", entry.path().display(), file_size);
                        removed_count += 1;
                        removed_size += file_size;
                    }
                }
            }
        }
    }

    println!("\n=== Cache Summary ===");
    println!("Total files: {}", total_files);
    println!("Total size: {} bytes ({:.2} MB)", total_size, total_size as f64 / 1_048_576.0);

    if dry_run {
        println!(
            "\nWould remove {} files ({} bytes, {:.2} MB)",
            removed_count,
            removed_size,
            removed_size as f64 / 1_048_576.0
        );
        println!("\nRun without --dry-run to actually delete files.");
    } else if removed_count > 0 {
        println!(
            "\nRemoved {} files ({} bytes, {:.2} MB)",
            removed_count,
            removed_size,
            removed_size as f64 / 1_048_576.0
        );
        println!(
            "Remaining: {} files ({} bytes, {:.2} MB)",
            total_files - removed_count,
            total_size - removed_size,
            (total_size - removed_size) as f64 / 1_048_576.0
        );
    } else {
        println!("\nNo files removed.");
    }
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
        let mut bam = bam::Reader::from_path("test-data/2D-candidates_test_set.bam").unwrap();
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
        let chrom_lengths = get_chrom_lengths_from_bam_header(bam, &None);
        assert_eq!(chrom_lengths.get("chr7").unwrap(), &159345973);
    }
}
