use crate::call::{GenotypeConfig, ProcessingConfig};
use crate::repeats::TargetConfig;
use indicatif::{ProgressBar, ProgressStyle};
use log::debug;
use rayon::prelude::*;
use std::fs;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

/// Helper function to catch panics and convert them to string errors
/// This unifies panic handling across different sample processing functions
fn catch_sample_panic<F>(f: F) -> Result<(), String>
where
    F: FnOnce() -> crate::errors::InquiSTRResult<()> + std::panic::UnwindSafe,
{
    match std::panic::catch_unwind(f) {
        Ok(Ok(())) => Ok(()),
        Ok(Err(e)) => Err(e.message),
        Err(e) => {
            // Try to extract error message from panic
            if let Some(s) = e.downcast_ref::<&str>() {
                Err(s.to_string())
            } else if let Some(s) = e.downcast_ref::<String>() {
                Err(s.clone())
            } else {
                Err("Unknown panic occurred during sample processing".to_string())
            }
        }
    }
}

/// Configuration for unmapped kmer counting
#[derive(Debug, Clone)]
pub struct UnmappedConfig {
    pub klength: usize,
    pub target_kmer: Option<String>,
    pub combine_revcomp: bool,
}

/// Batch processing configuration
#[derive(Debug, Clone)]
pub struct BatchConfig {
    pub manifest: PathBuf,
    pub output: PathBuf,
    pub save_individual: Option<PathBuf>,
    pub tmpdir: Option<PathBuf>,
    pub resume: bool,
    pub dry_run: bool,
    pub keep_going: bool,
    pub reference: Option<String>,
    pub parallel_samples: usize,
    pub profile: bool,
    pub skip_validation: bool,
}

/// Performance profiling data for a sample
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct SampleProfile {
    sample_name: String,
    duration_secs: f64,
    success: bool,
}

/// Mode-specific configuration for batch processing
#[derive(Debug, Clone)]
pub enum BatchMode {
    StrGenotyping {
        target_config: TargetConfig,
        genotype_config: GenotypeConfig,
        processing_config: ProcessingConfig,
    },
    UnmappedKmer {
        unmapped_config: UnmappedConfig,
        processing_config: ProcessingConfig,
    },
}

/// Sample information from the manifest file
#[derive(Debug, Clone)]
struct SampleInfo {
    bam_path: String,
    sample_name: String,
}

/// Parse the manifest TSV file
/// Format: bam_path\tsample_name (header required, sample_name optional)
fn parse_manifest(manifest_path: &Path) -> Result<Vec<SampleInfo>, String> {
    let file = std::fs::File::open(manifest_path)
        .map_err(|e| format!("Failed to open manifest file: {}", e))?;
    let reader = io::BufReader::new(file);
    let mut lines = reader.lines();

    // Read and validate header
    let header = lines
        .next()
        .ok_or("Manifest file is empty")?
        .map_err(|e| format!("Failed to read header: {}", e))?;

    let header_fields: Vec<&str> = header.split('\t').collect();
    if header_fields.is_empty() || header_fields[0] != "bam_path" {
        return Err(format!(
            "Invalid header. Expected 'bam_path\tsample_name' (tab-separated), got: '{}'",
            header
        ));
    }

    let has_sample_name_column = header_fields.len() > 1 && header_fields[1] == "sample_name";

    // Parse sample lines
    let mut samples = Vec::new();
    for (line_num, line_result) in lines.enumerate() {
        let line =
            line_result.map_err(|e| format!("Failed to read line {}: {}", line_num + 2, e))?;

        if line.trim().is_empty() {
            continue; // Skip empty lines
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }

        let bam_path = fields[0].to_string();

        let sample_name = if has_sample_name_column && fields.len() > 1 && !fields[1].is_empty() {
            fields[1].to_string()
        } else {
            // Extract sample name from BAM path
            crate::utils::extract_sample_name_from_path(&bam_path)
        };

        samples.push(SampleInfo { bam_path, sample_name });
    }

    if samples.is_empty() {
        return Err("No samples found in manifest file".to_string());
    }

    Ok(samples)
}

/// Check if a URL is accessible by performing a HEAD request
fn validate_url(url: &str) -> Result<(), String> {
    let client = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(10))
        .build()
        .map_err(|e| format!("Failed to create HTTP client: {}", e))?;

    let response = client
        .head(url)
        .send()
        .map_err(|e| format!("Failed to connect: {}", e))?;

    if response.status().is_success() {
        Ok(())
    } else {
        Err(format!("HTTP error: {}", response.status()))
    }
}

/// Check if a file path is accessible (local file or URL)
fn validate_file_path(path: &str) -> Result<(), String> {
    if path.starts_with("http://")
        || path.starts_with("https://")
        || path.starts_with("ftp://")
        || path.starts_with("ftps://")
        || path.starts_with("s3://")
    {
        validate_url(path)
    } else {
        // Local file
        if Path::new(path).exists() {
            Ok(())
        } else {
            Err("File does not exist".to_string())
        }
    }
}

/// Process a single sample for STR genotyping and write to output file
/// Writes directly to file to avoid stdout redirection conflicts in parallel execution
fn process_sample(
    sample: &SampleInfo,
    target_config: &TargetConfig,
    genotype_config: &GenotypeConfig,
    processing_config: &ProcessingConfig,
    reference: &Option<String>,
    output_path: &Path,
) -> Result<(), String> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    // Catch panics and convert to errors using helper
    catch_sample_panic(|| {
        // Get all genotypes
        let all_repeats = target_config.get_targets(&sample.bam_path, reference)?;
        let batches =
            crate::batching::create_batches(all_repeats, processing_config.batch_size_kb * 1000);

        // CRITICAL FIX: Create a SINGLE BAM reader and reuse it for all batches
        // This prevents file descriptor exhaustion with CRAM files
        let mut bam = if !genotype_config.unphased && !batches.is_empty() {
            // First batch needs phasing validation
            crate::bam_utils::get_bam_reader_with_validation(&sample.bam_path, reference)?
        } else {
            crate::bam_utils::get_bam_reader(&sample.bam_path, reference)?
        };

        let results: Result<Vec<Vec<crate::call::Genotype>>, crate::errors::InquiSTRError> =
            if processing_config.threads > 1 {
                // For parallel processing, we CANNOT share the BAM reader across threads
                // So we must fall back to the old approach (each batch opens its own reader)
                // This is a known limitation with CRAM files + parallel processing
                use rayon::prelude::*;
                let thread_pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(processing_config.threads)
                    .build()
                    .expect("Failed to build thread pool");

                thread_pool.install(|| {
                    batches
                        .into_par_iter()
                        .map(|batch| {
                            crate::genotype_batch::process_batch_worker(
                                batch,
                                &sample.bam_path,
                                reference,
                                *genotype_config,
                            )
                        })
                        .collect()
                })
            } else {
                // Sequential processing: Reuse the single BAM reader for all batches
                // This dramatically reduces file descriptor usage
                batches
                    .iter()
                    .map(|batch| {
                        crate::genotype_batch::process_batch_with_reader(
                            batch,
                            &mut bam,
                            genotype_config,
                        )
                    })
                    .collect()
            };

        // Explicitly drop the BAM reader now that we're done with all batches
        drop(bam);

        // Propagate the error if any batch failed
        let results = results?;

        let mut all_genotypes: Vec<crate::call::Genotype> = results.into_iter().flatten().collect();
        all_genotypes.sort_unstable();

        // Write to output file
        let file = File::create(output_path)
            .map_err(|e| format!("Failed to create output file {}: {}", output_path.display(), e))
            .unwrap();
        let mut writer = BufWriter::new(file);

        debug!("Writing individual file: {}", output_path.display());

        let sample_name = sample.sample_name.as_str();

        // Write metadata header
        writeln!(writer, "# file_type=individual_call").unwrap();
        writeln!(writer, "# version={}", crate::VERSION).unwrap();
        writeln!(writer, "# command=batch").unwrap();
        writeln!(writer, "# sample={}", sample_name).unwrap();
        writeln!(writer, "# minlen={}", genotype_config.minlen).unwrap();
        writeln!(writer, "# support={}", genotype_config.support).unwrap();
        writeln!(writer, "# unphased={}", genotype_config.unphased).unwrap();

        // Write data header
        writeln!(writer, "chromosome\tbegin\tend\tinfo\t{}_H1\t{}_H2", sample_name, sample_name)
            .unwrap();

        // Write genotypes
        for genotype in &all_genotypes {
            writeln!(writer, "{}", genotype).unwrap();
        }

        // Ensure all data is flushed to disk
        writer
            .flush()
            .map_err(|e| format!("Failed to flush output file: {}", e))
            .unwrap();
        debug!("Successfully wrote individual file: {}", output_path.display());

        Ok(())
    })
}

/// Process a single sample for unmapped kmer counting and write to output file
fn process_sample_unmapped(
    sample: &SampleInfo,
    unmapped_config: &UnmappedConfig,
    processing_config: &ProcessingConfig,
    reference: &Option<String>,
    output_path: &Path,
) -> Result<(), String> {
    use std::fs::OpenOptions;

    // Open output file
    let output_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_path)
        .map_err(|e| format!("Failed to create output file {}: {}", output_path.display(), e))?;

    // Redirect stdout to the file
    let _redirect = gag::Redirect::stdout(output_file)
        .map_err(|e| format!("Failed to redirect stdout: {}", e))?;

    // Catch panics and convert to errors using helper
    catch_sample_panic(|| {
        crate::unmapped::count_unmapped_kmers(
            sample.bam_path.clone(),
            unmapped_config.klength,
            Some(sample.sample_name.clone()),
            reference.clone(),
            processing_config.threads,
            unmapped_config.target_kmer.clone(),
            unmapped_config.combine_revcomp,
            false, // Don't show progress - batch mode has its own progress bar
        )
    })
    // redirect is dropped here, restoring stdout
}

/// Unified function to extract sample names from any combined output file
/// Automatically detects file format (STR, kmer regular, or kmer target) and extracts samples
/// Returns None if file doesn't exist or isn't a valid combined file
fn get_samples_from_combined_file(file_path: &Path) -> Option<Vec<String>> {
    if !file_path.exists() {
        return None;
    }

    let file = match fs::File::open(file_path) {
        Ok(f) => f,
        Err(_) => return None,
    };

    let reader = io::BufReader::new(file);
    let mut lines = reader.lines();

    // Read first line (might be metadata or header)
    let first_line = match lines.next() {
        Some(Ok(line)) => line,
        _ => return None,
    };

    // Skip metadata line if present
    let header = if first_line.starts_with("# file_type=") {
        match lines.next() {
            Some(Ok(line)) => line,
            _ => return None,
        }
    } else {
        first_line
    };

    let fields: Vec<&str> = header.split('\t').collect();

    // Detect format and extract samples accordingly
    match fields.first() {
        Some(&"chromosome") => extract_str_samples(&fields),
        Some(&"kmer") => extract_kmer_regular_samples(&fields),
        Some(&"Sample") => extract_kmer_target_samples(lines),
        _ => None,
    }
}

/// Extract sample names from STR file header columns
fn extract_str_samples(fields: &[&str]) -> Option<Vec<String>> {
    // Check if it has more than 5 columns (indicating it's a combined file)
    if fields.len() <= 5 {
        return None; // Individual file, not combined
    }

    // Format: chromosome, begin, end, info, then pairs of (sample_H1, sample_H2)
    let mut sample_names = Vec::new();
    let mut seen_samples = std::collections::HashSet::new();
    let mut i = 4; // Start after chromosome, begin, end, info

    while i < fields.len() {
        if i + 1 < fields.len() {
            // Extract base sample name by removing _H1 or _H2 suffix
            let col_name = fields[i];
            let base_name = col_name
                .strip_suffix("_H1")
                .or_else(|| col_name.strip_suffix("_H2"))
                .unwrap_or(col_name);

            // Only add each unique sample name once
            if !seen_samples.contains(base_name) {
                sample_names.push(base_name.to_string());
                seen_samples.insert(base_name.to_string());
            }

            i += 2; // Move to next pair
        } else {
            break;
        }
    }

    Some(sample_names)
}

/// Extract sample names from regular kmer file header (columns after "kmer")
fn extract_kmer_regular_samples(fields: &[&str]) -> Option<Vec<String>> {
    // Sample names are all columns after the first (kmer column)
    Some(fields[1..].iter().map(|s| s.to_string()).collect())
}

/// Extract sample names from target kmer file (samples listed in rows)
fn extract_kmer_target_samples<I>(lines: I) -> Option<Vec<String>>
where
    I: Iterator<Item = Result<String, std::io::Error>>,
{
    // Read all rows to get sample names from first column
    let sample_names: Vec<String> = lines
        .map_while(Result::ok)
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            fields.first().map(|s| s.to_string())
        })
        .collect();

    if sample_names.is_empty() {
        None
    } else {
        Some(sample_names)
    }
}

/// Legacy wrapper for backward compatibility - delegates to unified function
#[allow(dead_code)]
fn get_samples_from_combined_str_file(file_path: &Path) -> Option<Vec<String>> {
    get_samples_from_combined_file(file_path)
}

/// Legacy wrapper for backward compatibility - delegates to unified function
#[allow(dead_code)]
fn get_samples_from_combined_kmer_file(file_path: &Path) -> Option<Vec<String>> {
    get_samples_from_combined_file(file_path)
}

/// Main batch processing function
pub fn batch_process(config: BatchConfig, mode: BatchMode) {
    use std::time::Instant;

    let start_time = Instant::now();
    eprintln!("Starting batch processing...");

    // Get thread count and calculate threads per sample
    let total_threads = match &mode {
        BatchMode::UnmappedKmer { processing_config, .. } => processing_config.threads,
        BatchMode::StrGenotyping { processing_config, .. } => processing_config.threads,
    };

    let parallel_samples = config.parallel_samples;
    let threads_per_sample = (total_threads / parallel_samples).max(1);

    if parallel_samples > 1 {
        eprintln!("Parallel processing: {} samples at a time", parallel_samples);
        eprintln!(
            "Thread allocation: {} threads total, {} threads per sample",
            total_threads, threads_per_sample
        );

        if threads_per_sample == 1 {
            eprintln!(
                "  Mode: Samples in parallel, batches sequential per sample (optimal for CRAM)"
            );
        } else {
            eprintln!(
                "  Mode: Samples in parallel, {} threads per sample for batch parallelism",
                threads_per_sample
            );
        }
    } else {
        eprintln!("Sequential processing: 1 sample at a time with {} thread(s)", total_threads);
    }

    match &mode {
        BatchMode::UnmappedKmer { unmapped_config, .. } => {
            eprintln!("Mode: Unmapped kmer counting");
            eprintln!("  Max kmer length: {}", unmapped_config.klength);
            if let Some(ref target) = unmapped_config.target_kmer {
                eprintln!("  Target kmer: {}", target);
            }
            if unmapped_config.combine_revcomp {
                eprintln!("  Combining with reverse complements");
            }
        }
        BatchMode::StrGenotyping { target_config, .. } => {
            eprintln!("Mode: STR genotyping");
            if target_config.region.is_none()
                && target_config.region_file.is_none()
                && target_config.preset.is_none()
            {
                eprintln!(
                    "ERROR: For STR genotyping mode, you must provide --region, --region-file, or --preset"
                );
                std::process::exit(1);
            }
        }
    }
    eprintln!("Reading manifest: {}", config.manifest.display());

    // Parse manifest
    let mut samples = match parse_manifest(&config.manifest) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("ERROR: Failed to parse manifest: {}", e);
            std::process::exit(1);
        }
    };

    // Check for duplicate sample names
    let mut seen_names = std::collections::HashSet::new();
    let mut duplicates = std::collections::HashSet::new();
    for sample in &samples {
        if !seen_names.insert(&sample.sample_name) {
            duplicates.insert(sample.sample_name.clone());
        }
    }
    if !duplicates.is_empty() {
        eprintln!("ERROR: Duplicate sample names found in manifest:");
        let mut dup_list: Vec<_> = duplicates.into_iter().collect();
        dup_list.sort();
        for dup in &dup_list {
            eprintln!("  - {}", dup);
        }
        eprintln!(
            "\nSample names must be unique. Please fix the manifest or ensure unique sample_name column values."
        );
        std::process::exit(1);
    }

    // Validate all file paths (especially URLs) before starting processing
    if !config.skip_validation {
        eprintln!("Validating file paths...");
        let mut invalid_paths = Vec::new();
        for sample in &samples {
            if let Err(e) = validate_file_path(&sample.bam_path) {
                invalid_paths.push((sample.sample_name.clone(), sample.bam_path.clone(), e));
            }
        }
        if !invalid_paths.is_empty() {
            eprintln!("ERROR: {} invalid file path(s) found:", invalid_paths.len());
            for (sample_name, path, error) in &invalid_paths {
                eprintln!("  - {}: {} ({})", sample_name, path, error);
            }
            eprintln!("\nPlease fix the file paths in the manifest and try again.");
            std::process::exit(1);
        }
        eprintln!("  All {} file paths validated successfully", samples.len());
    } else {
        eprintln!("Skipping file path validation (--skip-validation enabled)");
    }

    // Determine output directory for individual files (need this early for resume logic)
    let individual_dir = if let Some(dir) = config.save_individual.as_ref() {
        // Create directory if it doesn't exist
        if !dir.exists() {
            if config.dry_run {
                eprintln!("Would create output directory: {}", dir.display());
            } else if let Err(e) = fs::create_dir_all(dir) {
                eprintln!("ERROR: Failed to create output directory {}: {}", dir.display(), e);
                std::process::exit(1);
            }
        }
        dir.clone()
    } else {
        // Use temporary directory
        config
            .tmpdir
            .clone()
            .or_else(|| std::env::var("TMPDIR").ok().map(PathBuf::from))
            .unwrap_or_else(|| PathBuf::from("."))
    };

    // If resume mode is enabled, check for completed samples
    let mut completed_samples = Vec::new();
    let mut samples_in_output = 0;
    let mut samples_with_individual_files = 0;

    if config.resume {
        eprintln!("Resume mode: Checking for existing output in: {}", individual_dir.display());

        // First check the combined output file
        let completed_from_output = get_samples_from_combined_file(&config.output);

        if let Some(completed_list) = completed_from_output {
            samples_in_output = completed_list.len();
            eprintln!(
                "Resume mode: Found {} samples in existing combined output file",
                samples_in_output
            );
            completed_samples.extend(completed_list);
        }

        // Also check for existing individual files
        for sample in &samples {
            let individual_file = individual_dir.join(format!("{}.inq", sample.sample_name));
            if individual_file.exists() && !completed_samples.contains(&sample.sample_name) {
                completed_samples.push(sample.sample_name.clone());
                samples_with_individual_files += 1;
                eprintln!("  Found existing file: {}", individual_file.display());
            }
        }

        if samples_with_individual_files > 0 {
            eprintln!(
                "Resume mode: Found {} existing individual .inq files (not yet in combined output)",
                samples_with_individual_files
            );
        }

        if !completed_samples.is_empty() {
            // Filter out already completed samples
            let original_count = samples.len();
            samples.retain(|s| !completed_samples.contains(&s.sample_name));
            let skipped = original_count - samples.len();

            eprintln!(
                "  Total skipped: {} (in output: {}, individual files: {})",
                skipped, samples_in_output, samples_with_individual_files
            );

            if samples.is_empty() {
                // All samples already processed
                if config.output.exists() {
                    // Combined output already exists, nothing to do
                    eprintln!("\nAll samples already processed and combined output exists.");
                    eprintln!("Output file: {}", config.output.display());
                    return;
                } else if samples_with_individual_files > 0 {
                    // Individual files exist but combined output doesn't - need to combine
                    eprintln!("\nAll samples already processed as individual files.");
                    eprintln!("Creating combined output file...");
                    // Will skip processing and jump directly to combine step
                } else {
                    // All samples in output but output file doesn't exist (shouldn't happen)
                    eprintln!(
                        "\nERROR: Inconsistent state - samples marked as complete but no files found."
                    );
                    std::process::exit(1);
                }
            } else {
                eprintln!("  Will process {} remaining sample(s)", samples.len());
            }
        } else {
            eprintln!(
                "Resume mode: No existing output or individual files found, will process all samples"
            );
        }
    } else {
        eprintln!("Found {} samples in manifest", samples.len());
    }

    // Dry-run: validate manifest and report what would be processed
    if config.dry_run {
        eprintln!("\nDRY RUN MODE - No processing will occur");
        eprintln!("\nValidating samples:");
        let mut valid_samples = 0;
        let mut invalid_samples = Vec::new();
        let mut existing_samples = 0;

        for sample in &samples {
            let individual_file = individual_dir.join(format!("{}.inq", sample.sample_name));

            // Check if individual file already exists (not yet in combined output)
            if individual_file.exists() && !completed_samples.contains(&sample.sample_name) {
                eprintln!("  ○ {}: Individual file exists, would combine", sample.sample_name);
                existing_samples += 1;
            // Check if BAM file exists (for local files)
            } else if !sample.bam_path.starts_with("http://")
                && !sample.bam_path.starts_with("https://")
                && !sample.bam_path.starts_with("ftp://")
                && !sample.bam_path.starts_with("s3://")
                && !std::path::Path::new(&sample.bam_path).exists()
            {
                eprintln!("  ✗ {}: BAM file not found: {}", sample.sample_name, sample.bam_path);
                invalid_samples.push(sample.sample_name.clone());
            } else {
                eprintln!("  ✓ {}: Would process", sample.sample_name);
                valid_samples += 1;
            }
        }

        eprintln!("\nSummary:");
        eprintln!("  Total samples in manifest: {}", samples.len() + completed_samples.len());
        if config.resume && samples_in_output > 0 {
            eprintln!("  Already in combined output: {}", samples_in_output);
        }
        if config.resume && samples_with_individual_files > 0 {
            eprintln!("  Individual files to combine: {}", samples_with_individual_files);
        }
        if existing_samples > 0 {
            eprintln!("  Individual files discovered in dry-run: {}", existing_samples);
        }
        eprintln!("  Would process from scratch: {}", valid_samples);
        if !invalid_samples.is_empty() {
            eprintln!("  Invalid: {}", invalid_samples.len());
        }

        if config.resume && config.output.exists() {
            eprintln!("\nOutput: Would append to existing file: {}", config.output.display());
        } else if config.resume && samples_with_individual_files > 0 {
            eprintln!(
                "\nOutput: Would combine {} individual files and create: {}",
                samples_with_individual_files,
                config.output.display()
            );
        } else {
            eprintln!("\nOutput: Would create new file: {}", config.output.display());
        }
        return;
    }

    // Create adjusted mode with threads per sample
    let adjusted_mode = match mode {
        BatchMode::UnmappedKmer { unmapped_config, processing_config } => BatchMode::UnmappedKmer {
            unmapped_config,
            processing_config: ProcessingConfig {
                threads: threads_per_sample,
                ..processing_config
            },
        },
        BatchMode::StrGenotyping { target_config, genotype_config, processing_config } => {
            BatchMode::StrGenotyping {
                target_config,
                genotype_config,
                processing_config: ProcessingConfig {
                    threads: threads_per_sample,
                    ..processing_config
                },
            }
        }
    };

    // Create progress tracking
    let individual_files = Arc::new(Mutex::new(Vec::new()));
    let failed_samples = Arc::new(Mutex::new(Vec::new()));
    let sample_profiles = if config.profile {
        Some(Arc::new(Mutex::new(Vec::new())))
    } else {
        None
    };

    // IMPORTANT WARNING: Check for CRAM files with parallel processing
    if parallel_samples > 1 {
        let cram_count = samples
            .iter()
            .filter(|s| s.bam_path.ends_with(".cram"))
            .count();
        if cram_count > 0 && threads_per_sample > 1 {
            // Warn only if using parallel batches per sample (threads_per_sample > 1)
            // This is the problematic configuration that opens many readers
            eprintln!(
                "\n⚠️  WARNING: Processing {} CRAM file(s) with {} threads per sample",
                cram_count, threads_per_sample
            );
            eprintln!("   Each sample will process batches in parallel, opening multiple readers.");
            eprintln!("   This may cause file descriptor exhaustion.");
            eprintln!("   RECOMMENDED: Set --parallel-samples equal to --threads (e.g., both 4)");
            eprintln!(
                "   This processes samples in parallel but batches sequentially (FD-efficient)."
            );
            eprintln!();
        } else if cram_count > 0 && threads_per_sample == 1 {
            // Good configuration! Acknowledge it
            eprintln!(
                "\n✓ Optimal CRAM configuration: {} samples in parallel, sequential batches per sample",
                parallel_samples
            );
            eprintln!(
                "  Each sample uses 1 thread (batches processed sequentially with single reader)"
            );
            eprintln!("  Estimated file descriptors: ~{} (very manageable)", parallel_samples * 15);
            eprintln!();
        }
    }

    // If all samples are already processed, skip processing and collect existing files
    if samples.is_empty() && samples_with_individual_files > 0 {
        // Collect all existing individual files
        for sample_name in &completed_samples {
            let individual_file = individual_dir.join(format!("{}.inq", sample_name));
            if individual_file.exists() {
                individual_files.lock().unwrap().push(individual_file);
            }
        }
        // Jump directly to combine step (skip the processing loop)
    } else if parallel_samples == 1 {
        // Sequential processing with simple progress bar
        let pb = ProgressBar::new(samples.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} (ETA: {eta}) {msg}")
                .unwrap()
                .progress_chars("=>-"),
        );

        for sample in &samples {
            pb.set_message(format!("Processing {}", sample.sample_name));

            let individual_file = individual_dir.join(format!("{}.inq", sample.sample_name));

            // Pre-validate BAM file exists (for local files)
            if !sample.bam_path.starts_with("http://")
                && !sample.bam_path.starts_with("https://")
                && !sample.bam_path.starts_with("ftp://")
                && !sample.bam_path.starts_with("s3://")
                && !std::path::Path::new(&sample.bam_path).exists()
            {
                eprintln!(
                    "\nERROR: Failed to process sample '{}': BAM file does not exist",
                    sample.sample_name
                );
                eprintln!("       BAM file: {}", sample.bam_path);
                failed_samples
                    .lock()
                    .unwrap()
                    .push(sample.sample_name.clone());
                pb.inc(1);
                continue;
            }

            let sample_start = if config.profile {
                Some(Instant::now())
            } else {
                None
            };
            let result = match &adjusted_mode {
                BatchMode::UnmappedKmer { unmapped_config, processing_config } => {
                    process_sample_unmapped(
                        sample,
                        unmapped_config,
                        processing_config,
                        &config.reference,
                        &individual_file,
                    )
                }
                BatchMode::StrGenotyping { target_config, genotype_config, processing_config } => {
                    process_sample(
                        sample,
                        target_config,
                        genotype_config,
                        processing_config,
                        &config.reference,
                        &individual_file,
                    )
                }
            };

            match result {
                Ok(()) => {
                    individual_files.lock().unwrap().push(individual_file);
                    if let (Some(profiles), Some(start)) = (&sample_profiles, sample_start) {
                        profiles.lock().unwrap().push(SampleProfile {
                            sample_name: sample.sample_name.clone(),
                            duration_secs: start.elapsed().as_secs_f64(),
                            success: true,
                        });
                    }
                }
                Err(e) => {
                    eprintln!("\nERROR: Failed to process sample '{}': {}", sample.sample_name, e);
                    eprintln!("       BAM file: {}", sample.bam_path);
                    failed_samples
                        .lock()
                        .unwrap()
                        .push(sample.sample_name.clone());
                    if let (Some(profiles), Some(start)) = (&sample_profiles, sample_start) {
                        profiles.lock().unwrap().push(SampleProfile {
                            sample_name: sample.sample_name.clone(),
                            duration_secs: start.elapsed().as_secs_f64(),
                            success: false,
                        });
                    }
                }
            }

            pb.inc(1);
        }

        pb.finish_with_message("Sample processing complete");
    } else {
        // Parallel processing with overall progress bar
        let overall_pb = ProgressBar::new(samples.len() as u64);
        overall_pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} (ETA: {eta}) {msg}")
                .unwrap()
                .progress_chars("=>-"),
        );
        overall_pb.set_message(format!("Processing {} samples in parallel", parallel_samples));

        // Set up rayon thread pool
        rayon::ThreadPoolBuilder::new()
            .num_threads(parallel_samples)
            .build()
            .unwrap()
            .install(|| {
                samples.par_iter().for_each(|sample| {

                    let individual_file =
                        individual_dir.join(format!("{}.inq", sample.sample_name));

                    // Pre-validate BAM file exists (for local files)
                    if !sample.bam_path.starts_with("http://")
                        && !sample.bam_path.starts_with("https://")
                        && !sample.bam_path.starts_with("ftp://")
                        && !sample.bam_path.starts_with("s3://")
                        && !std::path::Path::new(&sample.bam_path).exists()
                    {
                        eprintln!(
                            "\nERROR: Failed to process sample '{}': BAM file does not exist",
                            sample.sample_name
                        );
                        eprintln!("       BAM file: {}", sample.bam_path);
                        failed_samples
                            .lock()
                            .unwrap()
                            .push(sample.sample_name.clone());
                        overall_pb.inc(1);
                        return;
                    }

                    let sample_start = if config.profile { Some(Instant::now()) } else { None };
                    let result = match &adjusted_mode {
                        BatchMode::UnmappedKmer { unmapped_config, processing_config } => {
                            process_sample_unmapped(
                                sample,
                                unmapped_config,
                                processing_config,
                                &config.reference,
                                &individual_file,
                            )
                        }
                        BatchMode::StrGenotyping {
                            target_config,
                            genotype_config,
                            processing_config,
                        } => process_sample(
                            sample,
                            target_config,
                            genotype_config,
                            processing_config,
                            &config.reference,
                            &individual_file,
                        ),
                    };

                    match result {
                        Ok(()) => {
                            individual_files.lock().unwrap().push(individual_file);
                            if let (Some(profiles), Some(start)) = (&sample_profiles, sample_start) {
                                profiles.lock().unwrap().push(SampleProfile {
                                    sample_name: sample.sample_name.clone(),
                                    duration_secs: start.elapsed().as_secs_f64(),
                                    success: true,
                                });
                            }
                        }
                        Err(e) => {
                            eprintln!(
                                "\nERROR: Failed to process sample '{}': {}",
                                sample.sample_name, e
                            );
                            eprintln!("       BAM file: {}", sample.bam_path);
                            failed_samples
                                .lock()
                                .unwrap()
                                .push(sample.sample_name.clone());
                            if let (Some(profiles), Some(start)) = (&sample_profiles, sample_start) {
                                profiles.lock().unwrap().push(SampleProfile {
                                    sample_name: sample.sample_name.clone(),
                                    duration_secs: start.elapsed().as_secs_f64(),
                                    success: false,
                                });
                            }

                            // If keep_going is not set, stop processing immediately
                            if !config.keep_going {
                                overall_pb.finish_with_message("Stopping due to error (use --keep-going to continue)");
                                eprintln!("\nStopping batch processing due to error.");
                                eprintln!("Use --keep-going to continue processing remaining samples despite failures.");
                                std::process::exit(1);
                            }
                        }
                    }

                    overall_pb.inc(1);
                });
            });

        overall_pb.finish_with_message("Sample processing complete");

        // Give runtime a chance to clean up file descriptors from parallel processing
        // This is especially important with CRAM files which open many FDs per sample
        std::thread::sleep(std::time::Duration::from_millis(100));
    }

    // Extract results from Arc<Mutex<>>
    let individual_files = Arc::try_unwrap(individual_files)
        .unwrap()
        .into_inner()
        .unwrap();
    let failed_samples = Arc::try_unwrap(failed_samples)
        .unwrap()
        .into_inner()
        .unwrap();

    // Report failures
    if !failed_samples.is_empty() {
        eprintln!("\n{} sample(s) failed:", failed_samples.len());
        for failed in &failed_samples {
            eprintln!("  - {}", failed);
        }
    }

    // Check if we should proceed with combining files
    if individual_files.is_empty() {
        eprintln!("\nERROR: No samples were successfully processed. Cannot create combined file.");
        std::process::exit(1);
    }

    // If some samples failed, warn but continue to create combined file with successful samples
    if !failed_samples.is_empty() {
        eprintln!(
            "\nWARNING: {} sample(s) failed, but {} succeeded.",
            failed_samples.len(),
            individual_files.len()
        );
        eprintln!("Creating combined file with successful samples only.");
    }

    // Get thread count from adjusted_mode
    let threads = match &adjusted_mode {
        BatchMode::UnmappedKmer { processing_config, .. } => processing_config.threads,
        BatchMode::StrGenotyping { processing_config, .. } => processing_config.threads,
    };

    // Combine files - either create new or append to existing
    if config.resume && config.output.exists() {
        // Resume mode with existing output: append new samples to existing combined file
        eprintln!(
            "\nAppending {} new sample(s) to existing combined file: {}",
            individual_files.len(),
            config.output.display()
        );

        // combine supports adding individual files to a combined file
        // We pass the existing combined file + new individual files
        let mut all_files = vec![config.output.clone()];
        all_files.extend(individual_files.clone());

        // Redirect stdout to temporary file first
        use std::fs::OpenOptions;
        let temp_output = config.output.with_extension("tmp");
        let output_file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&temp_output)
            .unwrap_or_else(|e| {
                eprintln!("ERROR: Failed to create temporary output file: {}", e);
                std::process::exit(1);
            });

        let _redirect = gag::Redirect::stdout(output_file).unwrap_or_else(|e| {
            eprintln!("ERROR: Failed to redirect stdout: {}", e);
            std::process::exit(1);
        });

        crate::combine::combine(all_files, threads);
        drop(_redirect); // Explicitly drop to restore stdout before file operations

        // Replace original with updated file
        if let Err(e) = fs::rename(&temp_output, &config.output) {
            eprintln!("ERROR: Failed to update output file: {}", e);
            std::process::exit(1);
        }
    } else {
        // Normal mode: create new combined file
        eprintln!(
            "\nCombining {} sample(s) into {}",
            individual_files.len(),
            config.output.display()
        );

        use std::fs::OpenOptions;
        let output_file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&config.output)
            .unwrap_or_else(|e| {
                eprintln!("ERROR: Failed to create output file {}: {}", config.output.display(), e);
                std::process::exit(1);
            });

        let _redirect = gag::Redirect::stdout(output_file).unwrap_or_else(|e| {
            eprintln!("ERROR: Failed to redirect stdout: {}", e);
            std::process::exit(1);
        });

        crate::combine::combine(individual_files.clone(), threads);
    } // redirect is dropped here, restoring stdout

    // Clean up temporary files if not saving individual files
    if config.save_individual.is_none() {
        eprintln!("Cleaning up temporary files...");
        for file in &individual_files {
            if let Err(e) = fs::remove_file(file) {
                eprintln!("Warning: Failed to remove temporary file {}: {}", file.display(), e);
            }
        }
    } else {
        eprintln!(
            "Individual sample files saved to: {}",
            config.save_individual.as_ref().unwrap().display()
        );
    }

    let total_duration = start_time.elapsed();
    eprintln!("\nBatch processing complete!");
    eprintln!("  Successful: {}", individual_files.len());
    eprintln!("  Failed: {}", failed_samples.len());
    eprintln!("  Output: {}", config.output.display());
    eprintln!("  Total time: {:.1}s", total_duration.as_secs_f64());

    // Generate profiling report if requested
    if config.profile
        && let Some(profiles) = sample_profiles
    {
        let profiles = profiles.lock().unwrap();
        generate_profile_report(
            &profiles,
            total_duration.as_secs_f64(),
            total_threads,
            parallel_samples,
            &adjusted_mode,
        );
    }

    // If we reach here with failures, it means --keep-going was set
    // (otherwise we would have exited immediately on first failure)
    if !failed_samples.is_empty() {
        eprintln!(
            "\nWARNING: {} sample(s) failed but --keep-going was specified.",
            failed_samples.len()
        );
        eprintln!(
            "Combined file contains only the {} successful sample(s).",
            individual_files.len()
        );
        std::process::exit(1); // Exit with error code to indicate failures occurred
    }
}

/// Generate a performance profile report with recommendations
fn generate_profile_report(
    profiles: &[SampleProfile],
    total_duration: f64,
    total_threads: usize,
    parallel_samples: usize,
    mode: &BatchMode,
) {
    if profiles.is_empty() {
        return;
    }

    eprintln!("\n╔════════════════════════════════════════════════════════════╗");
    eprintln!("║              PERFORMANCE PROFILE REPORT                    ║");
    eprintln!("╚════════════════════════════════════════════════════════════╝");

    // Calculate statistics
    let successful: Vec<_> = profiles.iter().filter(|p| p.success).collect();
    let n = successful.len();

    if n == 0 {
        eprintln!("\nNo successful samples to profile.");
        return;
    }

    let durations: Vec<f64> = successful.iter().map(|p| p.duration_secs).collect();
    let mean = durations.iter().sum::<f64>() / n as f64;
    let min = durations.iter().copied().fold(f64::INFINITY, f64::min);
    let max = durations.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    // Calculate standard deviation
    let variance = durations.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / n as f64;
    let stddev = variance.sqrt();
    let cv = if mean > 0.0 {
        (stddev / mean) * 100.0
    } else {
        0.0
    };

    eprintln!("\nSample Processing Times:");
    eprintln!("  Samples:    {}", n);
    eprintln!("  Mean:       {:.1}s", mean);
    eprintln!("  Std Dev:    {:.1}s ({:.1}% CV)", stddev, cv);
    eprintln!("  Min:        {:.1}s", min);
    eprintln!("  Max:        {:.1}s", max);
    eprintln!("  Range:      {:.1}s ({:.1}x slower)", max - min, max / min);

    // Configuration info
    eprintln!("\nConfiguration:");
    eprintln!("  Total threads:      {}", total_threads);
    eprintln!("  Parallel samples:   {}", parallel_samples);
    eprintln!("  Threads per sample: {}", total_threads / parallel_samples);

    if let BatchMode::StrGenotyping { processing_config, .. } = mode {
        eprintln!("  Batch size:         {}kb", processing_config.batch_size_kb);
    }

    // Efficiency metrics
    let theoretical_min = mean * n as f64 / parallel_samples as f64;
    let parallel_efficiency = if total_duration > 0.0 {
        (theoretical_min / total_duration) * 100.0
    } else {
        0.0
    };

    eprintln!("\nEfficiency:");
    eprintln!(
        "  Theoretical min:    {:.1}s ({} samples / {} parallel)",
        theoretical_min, n, parallel_samples
    );
    eprintln!("  Actual:             {:.1}s", total_duration);
    eprintln!("  Parallel efficiency: {:.1}%", parallel_efficiency);

    // Generate recommendations
    eprintln!("\n╔════════════════════════════════════════════════════════════╗");
    eprintln!("║                    RECOMMENDATIONS                         ║");
    eprintln!("╚════════════════════════════════════════════════════════════╝");

    let mut recommendations = Vec::new();

    // High variability suggests I/O contention
    if cv > 30.0 && parallel_samples > 1 {
        recommendations.push(format!(
            "⚠️  High variability (CV={:.1}%) with parallel samples suggests I/O contention.\n\
             →  Try: Reduce --parallel-samples to {} or increase --batch-size to {} (from {}kb)",
            cv,
            (parallel_samples / 2).max(1),
            if let BatchMode::StrGenotyping { processing_config, .. } = mode {
                processing_config.batch_size_kb * 2
            } else {
                100
            },
            if let BatchMode::StrGenotyping { processing_config, .. } = mode {
                processing_config.batch_size_kb
            } else {
                50
            }
        ));
    }

    // Low parallel efficiency
    if parallel_efficiency < 70.0 && parallel_samples > 1 {
        recommendations.push(format!(
            "⚠️  Low parallel efficiency ({:.1}%) indicates overhead or I/O bottleneck.\n\
             →  Try: Process samples sequentially (--parallel-samples 1) with more threads per sample",
            parallel_efficiency
        ));
    }

    // Good efficiency but single sample processing
    if parallel_samples == 1 && parallel_efficiency > 90.0 && n > 5 {
        recommendations.push(format!(
            "✓  Excellent single-sample performance! For large cohorts, consider:\n\
             →  Try: --parallel-samples {} --threads {} to process multiple samples simultaneously",
            (total_threads / 2).max(2),
            total_threads
        ));
    }

    // Very fast samples - might benefit from more parallelization
    if mean < 30.0 && parallel_samples == 1 && n > 10 {
        recommendations.push(format!(
            "✓  Fast sample processing ({:.1}s avg). You could process more samples in parallel:\n\
             →  Try: --parallel-samples {} (will still use {} threads total)",
            mean,
            (total_threads / 2).max(2),
            total_threads
        ));
    }

    // Slow samples with low thread count
    if mean > 300.0 && total_threads < 8 {
        recommendations.push(format!(
            "⚠️  Slow sample processing ({:.1}s avg) with only {} threads.\n\
             →  Try: Increase --threads to use more CPU cores if available",
            mean, total_threads
        ));
    }

    // High variability in sequential mode
    if cv > 20.0
        && parallel_samples == 1
        && let BatchMode::StrGenotyping { processing_config, .. } = mode
    {
        recommendations.push(format!(
            "⚠️  High variability (CV={:.1}%) suggests varying sample complexity or I/O issues.\n\
                 →  Try: Increase --batch-size to {} (currently {}kb) to reduce I/O overhead",
            cv,
            processing_config.batch_size_kb * 2,
            processing_config.batch_size_kb
        ));
    }

    if recommendations.is_empty() {
        eprintln!("\n✓  Performance looks good! Configuration appears well-tuned for your system.");
    } else {
        eprintln!();
        for (i, rec) in recommendations.iter().enumerate() {
            if i > 0 {
                eprintln!();
            }
            eprintln!("{}", rec);
        }
    }

    eprintln!("\n╔════════════════════════════════════════════════════════════╗");
    eprintln!("║  Note: Run with --profile again after adjusting settings  ║");
    eprintln!("╚════════════════════════════════════════════════════════════╝\n");
}
