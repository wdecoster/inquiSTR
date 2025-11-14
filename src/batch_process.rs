use crate::call::{GenotypeConfig, ProcessingConfig, TargetConfig};
use indicatif::{ProgressBar, ProgressStyle};
use std::fs;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};

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
            "Invalid header. Expected 'bam_path' as first column, got: '{}'",
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
            extract_sample_name_from_path(&bam_path)
        };

        samples.push(SampleInfo { bam_path, sample_name });
    }

    if samples.is_empty() {
        return Err("No samples found in manifest file".to_string());
    }

    Ok(samples)
}

/// Extract sample name from file path (same logic as in call.rs)
fn extract_sample_name_from_path(path: &str) -> String {
    // Handle URLs
    let path_for_name = if path.starts_with("http://")
        || path.starts_with("https://")
        || path.starts_with("ftp://")
        || path.starts_with("ftps://")
        || path.starts_with("s3://")
    {
        // Extract filename from URL
        path.split('/').next_back().unwrap_or(path)
    } else {
        path
    };

    Path::new(path_for_name)
        .file_stem()
        .and_then(|s| s.to_str())
        .map(|s| {
            // Remove common extensions (.bam, .cram, .bai, .crai)
            s.trim_end_matches(".bam")
                .trim_end_matches(".cram")
                .trim_end_matches(".bai")
                .trim_end_matches(".crai")
                .to_string()
        })
        .unwrap_or_else(|| "sample".to_string())
}

/// Process a single sample and write to output file
/// This uses gag crate to redirect stdout to a file
fn process_sample(
    sample: &SampleInfo,
    target_config: &TargetConfig,
    genotype_config: &GenotypeConfig,
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
    let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
        crate::call::genotype_repeats(
            sample.bam_path.clone(),
            target_config.clone(),
            *genotype_config,
            *processing_config,
            Some(sample.sample_name.clone()),
            reference.clone(),
        );
    }));

    // redirect is dropped here, restoring stdout

    match result {
        Ok(_) => Ok(()),
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

/// Main batch processing function
pub fn batch_process(
    manifest: PathBuf,
    output: PathBuf,
    save_individual: Option<PathBuf>,
    tmpdir: Option<PathBuf>,
    target_config: TargetConfig,
    genotype_config: GenotypeConfig,
    processing_config: ProcessingConfig,
    reference: Option<String>,
) {
    eprintln!("Starting batch processing...");
    eprintln!("Reading manifest: {}", manifest.display());

    // Parse manifest
    let samples = match parse_manifest(&manifest) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("ERROR: Failed to parse manifest: {}", e);
            std::process::exit(1);
        }
    };

    eprintln!("Found {} samples in manifest", samples.len());

    // Determine output directory for individual files
    let individual_dir = if let Some(dir) = save_individual.as_ref() {
        // Create directory if it doesn't exist
        if !dir.exists() {
            if let Err(e) = fs::create_dir_all(dir) {
                eprintln!("ERROR: Failed to create output directory {}: {}", dir.display(), e);
                std::process::exit(1);
            }
        }
        dir.clone()
    } else {
        // Use temporary directory
        tmpdir
            .or_else(|| std::env::var("TMPDIR").ok().map(PathBuf::from))
            .unwrap_or_else(|| PathBuf::from("."))
    };

    // Create progress bar
    let pb = ProgressBar::new(samples.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("=>-"),
    );

    let mut individual_files = Vec::new();
    let mut failed_samples = Vec::new();

    // Process each sample
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
            failed_samples.push(sample.sample_name.clone());
            pb.inc(1);
            continue;
        }

        match process_sample(
            sample,
            &target_config,
            &genotype_config,
            &processing_config,
            &reference,
            &individual_file,
        ) {
            Ok(()) => {
                individual_files.push(individual_file);
            }
            Err(e) => {
                eprintln!("\nERROR: Failed to process sample '{}': {}", sample.sample_name, e);
                eprintln!("       BAM file: {}", sample.bam_path);
                failed_samples.push(sample.sample_name.clone());
            }
        }

        pb.inc(1);
    }

    pb.finish_with_message("Sample processing complete");

    // Report failures
    if !failed_samples.is_empty() {
        eprintln!("\n{} sample(s) failed:", failed_samples.len());
        for failed in &failed_samples {
            eprintln!("  - {}", failed);
        }
    }

    // Combine individual files
    if individual_files.is_empty() {
        eprintln!("\nERROR: No samples were successfully processed. Cannot create combined file.");
        std::process::exit(1);
    }

    eprintln!(
        "\nCombining {} successful sample(s) into {}",
        individual_files.len(),
        output.display()
    );

    // Redirect stdout to output file for combine
    {
        use std::fs::OpenOptions;
        let output_file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&output)
            .unwrap_or_else(|e| {
                eprintln!("ERROR: Failed to create output file {}: {}", output.display(), e);
                std::process::exit(1);
            });

        let _redirect = gag::Redirect::stdout(output_file).unwrap_or_else(|e| {
            eprintln!("ERROR: Failed to redirect stdout: {}", e);
            std::process::exit(1);
        });

        crate::combine::combine(individual_files.clone(), processing_config.threads);
    } // redirect is dropped here, restoring stdout

    // Clean up temporary files if not saving individual files
    if save_individual.is_none() {
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
    eprintln!("  Output: {}", output.display());
}
