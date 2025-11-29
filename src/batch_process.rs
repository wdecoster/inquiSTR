use crate::call::{GenotypeConfig, ProcessingConfig, TargetConfig};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use log::debug;
use rayon::prelude::*;
use std::fs;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

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
    use std::panic;

    // Catch panics and convert to errors
    let result = panic::catch_unwind(panic::AssertUnwindSafe(
        || -> Result<(), crate::errors::InquiSTRError> {
            // Get all genotypes
            let repeats = target_config.get_targets(&sample.bam_path, reference)?;
            let all_repeats: Vec<crate::repeats::RepeatInterval> = repeats.collect();
            let batches = crate::locus_batching::create_batches(
                all_repeats,
                processing_config.batch_size_kb * 1000,
            );

            let results: Result<Vec<Vec<crate::call::Genotype>>, crate::errors::InquiSTRError> =
                if processing_config.threads > 1 {
                    use rayon::prelude::*;
                    let thread_pool = rayon::ThreadPoolBuilder::new()
                        .num_threads(processing_config.threads)
                        .build()
                        .expect("Failed to build thread pool");

                    thread_pool.install(|| {
                        batches
                            .into_par_iter()
                            .map(|batch| {
                                crate::locus_batching::process_batch_worker(
                                    batch,
                                    &sample.bam_path,
                                    reference,
                                    *genotype_config,
                                )
                            })
                            .collect()
                    })
                } else {
                    batches
                        .into_iter()
                        .map(|batch| {
                            crate::locus_batching::process_batch_worker(
                                batch,
                                &sample.bam_path,
                                reference,
                                *genotype_config,
                            )
                        })
                        .collect()
                };

            // Propagate the error if any batch failed
            let results = results?;

            let mut all_genotypes: Vec<crate::call::Genotype> =
                results.into_iter().flatten().collect();
            all_genotypes.sort_unstable();

            // Write to output file
            let file = File::create(output_path)
                .map_err(|e| {
                    format!("Failed to create output file {}: {}", output_path.display(), e)
                })
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
            writeln!(
                writer,
                "chromosome\tbegin\tend\tinfo\t{}_H1\t{}_H2",
                sample_name, sample_name
            )
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
        },
    ));

    match result {
        Ok(Ok(())) => Ok(()),
        Ok(Err(e)) => Err(e.message), // Propagate the InquiSTRError
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

/// Process a single sample for unmapped kmer counting and write to output file
fn process_sample_unmapped(
    sample: &SampleInfo,
    unmapped_config: &UnmappedConfig,
    processing_config: &ProcessingConfig,
    reference: &Option<String>,
    output_path: &Path,
) -> Result<(), String> {
    use std::fs::OpenOptions;
    use std::panic;

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

    // Catch panics and convert to errors
    let result = panic::catch_unwind(panic::AssertUnwindSafe(
        || -> Result<(), crate::errors::InquiSTRError> {
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
        },
    ));

    // redirect is dropped here, restoring stdout

    match result {
        Ok(Ok(())) => Ok(()),
        Ok(Err(e)) => Err(e.message), // Propagate the InquiSTRError
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

/// Extract sample names from a combined STR file header
/// Returns None if file doesn't exist or isn't a valid combined file
fn get_samples_from_combined_str_file(file_path: &Path) -> Option<Vec<String>> {
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

    // Check if it's a valid STR file header
    let fields: Vec<&str> = header.split('\t').collect();
    if fields.first() != Some(&"chromosome") {
        return None;
    }

    // Check if it has more than 5 columns (indicating it's a combined file)
    if fields.len() <= 5 {
        return None; // Individual file, not combined
    }

    // Extract sample names from column headers
    // Format: chromosome, begin, end, then pairs of (sample_H1, sample_H2) for each sample
    let mut sample_names = Vec::new();
    let mut seen_samples = std::collections::HashSet::new();
    let mut i = 3; // Start after chromosome, begin, end

    while i < fields.len() {
        if i + 1 < fields.len() {
            // Extract base sample name by removing _H1 or _H2 suffix
            let col_name = fields[i + 1];
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

/// Extract sample names from a combined kmer file header
/// Returns None if file doesn't exist or isn't a valid combined kmer file
fn get_samples_from_combined_kmer_file(file_path: &Path) -> Option<Vec<String>> {
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

    // Check if it's a valid kmer file header
    let fields: Vec<&str> = header.split('\t').collect();

    // Check for regular kmer file format (kmer, sample1, sample2, ...)
    if fields.first() == Some(&"kmer") {
        // Sample names are all columns after the first
        return Some(fields[1..].iter().map(|s| s.to_string()).collect());
    }

    // Check for target kmer file format (Sample, count, ...)
    if fields.first() == Some(&"Sample") {
        // This format lists samples in rows, not columns
        // Read all rows to get sample names
        let mut sample_names = Vec::new();
        for line in lines.map_while(Result::ok) {
            let line_fields: Vec<&str> = line.split('\t').collect();
            if !line_fields.is_empty() {
                sample_names.push(line_fields[0].to_string());
            }
        }
        return Some(sample_names);
    }

    None
}

/// Main batch processing function
pub fn batch_process(config: BatchConfig, mode: BatchMode) {
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
        let completed_from_output = match &mode {
            BatchMode::StrGenotyping { .. } => get_samples_from_combined_str_file(&config.output),
            BatchMode::UnmappedKmer { .. } => get_samples_from_combined_kmer_file(&config.output),
        };

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
            if individual_file.exists()
                && !completed_samples.contains(&sample.sample_name)
            {
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
                }
                Err(e) => {
                    eprintln!("\nERROR: Failed to process sample '{}': {}", sample.sample_name, e);
                    eprintln!("       BAM file: {}", sample.bam_path);
                    failed_samples
                        .lock()
                        .unwrap()
                        .push(sample.sample_name.clone());
                }
            }

            pb.inc(1);
        }

        pb.finish_with_message("Sample processing complete");
    } else {
        // Parallel processing with multi-progress
        let multi = MultiProgress::new();
        let overall_pb = multi.add(ProgressBar::new(samples.len() as u64));
        overall_pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} (ETA: {eta}) Overall progress")
                .unwrap()
                .progress_chars("=>-"),
        );

        // Set up rayon thread pool
        rayon::ThreadPoolBuilder::new()
            .num_threads(parallel_samples)
            .build()
            .unwrap()
            .install(|| {
                samples.par_iter().for_each(|sample| {
                    let pb = multi.add(ProgressBar::new_spinner());
                    pb.set_style(
                        ProgressStyle::default_spinner()
                            .template("{spinner:.green} {msg}")
                            .unwrap(),
                    );
                    pb.set_message(format!("Processing {}", sample.sample_name));

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
                        pb.finish_with_message(format!("✗ {}", sample.sample_name));
                        overall_pb.inc(1);
                        return;
                    }

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
                            pb.finish_with_message(format!("✓ {}", sample.sample_name));
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
                            pb.finish_with_message(format!("✗ {}", sample.sample_name));

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
    }

    eprintln!("\nBatch processing complete!");
    eprintln!("  Successful: {}", individual_files.len());
    eprintln!("  Failed: {}", failed_samples.len());
    eprintln!("  Output: {}", config.output.display());

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
