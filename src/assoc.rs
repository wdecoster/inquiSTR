//! Association testing module
//!
//! This module provides an interface to perform STR association testing by wrapping
//! the R script STR_regression.R. It handles R environment detection, dependency checking,
//! and execution of the statistical analysis.

use log::{error, info, warn};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

// Embed the R script directly in the binary
const STR_REGRESSION_SCRIPT: &str = include_str!("../scripts/STR_regression.R");

// Required R packages
const REQUIRED_R_PACKAGES: &[&str] = &["data.table", "argparser", "parallel"];

/// Main function to run STR association testing
#[allow(clippy::too_many_arguments)]
pub fn run_association(
    input: PathBuf,
    phenocovar: PathBuf,
    phenotype: String,
    out: PathBuf,
    str_mode: String,
    outcometype: String,
    covnames: Option<String>,
    missing_cutoff: f64,
    minimal_length: Option<f64>,
    threads: usize,
    chunk_size: usize,
    binary_order: Option<String>,
    quiet: bool,
) {
    if !quiet {
        info!("Starting STR association testing");
    }

    // Validate input arguments
    if let Err(e) = validate_arguments(&input, &phenocovar, &str_mode, &outcometype, &binary_order)
    {
        error!("Argument validation failed: {}", e);
        std::process::exit(1);
    }

    // Check R environment
    if let Err(e) = check_r_environment(quiet) {
        error!("R environment check failed: {}", e);
        print_r_setup_instructions();
        std::process::exit(1);
    }

    // Create temporary R script file
    let temp_script = match create_temp_r_script() {
        Ok(path) => path,
        Err(e) => {
            error!("Failed to create temporary R script: {}", e);
            std::process::exit(1);
        }
    };

    // Execute R script
    let result = execute_r_script(
        &temp_script,
        input,
        phenocovar,
        phenotype,
        out,
        str_mode,
        outcometype,
        covnames,
        missing_cutoff,
        minimal_length,
        threads,
        chunk_size,
        binary_order,
        quiet,
    );

    // Clean up temporary file
    if let Err(e) = fs::remove_file(&temp_script) {
        warn!("Failed to remove temporary R script: {}", e);
    }

    match result {
        Ok(()) => {
            if !quiet {
                info!("Association testing completed successfully");
            }
        }
        Err(e) => {
            error!("Association testing failed: {}", e);
            std::process::exit(1);
        }
    }
}

/// Validate command-line arguments
fn validate_arguments(
    input: &Path,
    phenocovar: &Path,
    str_mode: &str,
    outcometype: &str,
    binary_order: &Option<String>,
) -> Result<(), String> {
    // Check input files exist
    if !input.exists() {
        return Err(format!("Input file does not exist: {}", input.display()));
    }

    if !phenocovar.exists() {
        return Err(format!("Phenotype file does not exist: {}", phenocovar.display()));
    }

    // Validate STR mode
    if !["MEAN", "MAX", "MIN"].contains(&str_mode) {
        return Err(format!("STR mode must be MEAN, MAX, or MIN, got: {}", str_mode));
    }

    // Validate outcome type
    if !["binary", "continuous"].contains(&outcometype) {
        return Err(format!("Outcome type must be binary or continuous, got: {}", outcometype));
    }

    // Check binary order for binary outcomes
    if outcometype == "binary" && binary_order.is_none() {
        return Err("Binary order is required for binary outcomes. Use --binary-order".to_string());
    }

    Ok(())
}

/// Check if R is installed and required packages are available
fn check_r_environment(quiet: bool) -> Result<(), String> {
    if !quiet {
        info!("Checking R environment...");
    }

    // Check if R is installed
    let r_version = Command::new("R").args(["--version"]).output();

    match r_version {
        Ok(output) => {
            if !output.status.success() {
                return Err("R is installed but not functioning properly".to_string());
            }
            if !quiet {
                info!("R installation found");
            }
        }
        Err(_) => {
            return Err("R is not installed or not in PATH".to_string());
        }
    }

    // Check for required packages
    let missing_packages = check_r_packages(quiet)?;

    if !missing_packages.is_empty() {
        return Err(format!("Missing required R packages: {}", missing_packages.join(", ")));
    }

    if !quiet {
        info!("All required R packages are available");
    }

    Ok(())
}

/// Check which required R packages are missing
fn check_r_packages(quiet: bool) -> Result<Vec<String>, String> {
    let mut missing_packages = Vec::new();

    for package in REQUIRED_R_PACKAGES {
        if !quiet {
            info!("Checking R package: {}", package);
        }

        let check_cmd = format!(
            "if (!require('{}', quietly = TRUE)) quit(status = 1) else quit(status = 0)",
            package
        );

        let result = Command::new("R")
            .args(["--slave", "-e", &check_cmd])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status();

        match result {
            Ok(status) => {
                if !status.success() {
                    missing_packages.push(package.to_string());
                }
            }
            Err(_) => {
                return Err("Failed to check R packages".to_string());
            }
        }
    }

    Ok(missing_packages)
}

/// Create a temporary R script file
fn create_temp_r_script() -> Result<PathBuf, std::io::Error> {
    let temp_dir = std::env::temp_dir();
    let temp_file = temp_dir.join("str_regression_temp.R");

    let mut file = fs::File::create(&temp_file)?;
    file.write_all(STR_REGRESSION_SCRIPT.as_bytes())?;

    Ok(temp_file)
}

/// Execute the R script with provided arguments
#[allow(clippy::too_many_arguments)]
fn execute_r_script(
    script_path: &Path,
    input: PathBuf,
    phenocovar: PathBuf,
    phenotype: String,
    out: PathBuf,
    str_mode: String,
    outcometype: String,
    covnames: Option<String>,
    missing_cutoff: f64,
    minimal_length: Option<f64>,
    threads: usize,
    chunk_size: usize,
    binary_order: Option<String>,
    quiet: bool,
) -> Result<(), String> {
    let script_path_str = script_path.to_string_lossy();
    let input_str = input.to_string_lossy();
    let phenocovar_str = phenocovar.to_string_lossy();
    let out_str = out.to_string_lossy();
    let missing_cutoff_str = missing_cutoff.to_string();
    let threads_str = threads.to_string();
    let chunk_size_str = chunk_size.to_string();

    let mut args = vec![
        script_path_str.as_ref(),
        "--input",
        input_str.as_ref(),
        "--phenocovar",
        phenocovar_str.as_ref(),
        "--phenotype",
        &phenotype,
        "--out",
        out_str.as_ref(),
        "--STRmode",
        &str_mode,
        "--outcometype",
        &outcometype,
        "--missing_cutoff",
        &missing_cutoff_str,
        "--threads",
        &threads_str,
        "--chunk_size",
        &chunk_size_str,
    ];

    // Add optional arguments
    if let Some(ref covnames_val) = covnames {
        args.extend(&["--covnames", covnames_val]);
    }

    let minimal_length_str = minimal_length.map(|val| val.to_string());
    if let Some(ref minimal_length_val) = minimal_length_str {
        args.extend(&["--minimal_length", minimal_length_val]);
    }

    if let Some(ref binary_order_val) = binary_order {
        args.extend(&["--binaryOrder", binary_order_val]);
    }

    if quiet {
        args.push("--quiet");
    }

    info!("Executing R script with arguments: {:?}", args);

    let output = Command::new("Rscript")
        .args(&args)
        .output()
        .map_err(|e| format!("Failed to execute R script: {}", e))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);

        error!("R script execution failed");
        error!("STDOUT: {}", stdout);
        error!("STDERR: {}", stderr);

        return Err(format!("R script failed with exit code: {:?}", output.status.code()));
    }

    if !quiet {
        let stdout = String::from_utf8_lossy(&output.stdout);
        if !stdout.is_empty() {
            println!("{}", stdout);
        }
    }

    Ok(())
}

/// Print instructions for setting up R environment
fn print_r_setup_instructions() {
    eprintln!("\n{}", "=".repeat(60));
    eprintln!("  R ENVIRONMENT SETUP INSTRUCTIONS");
    eprintln!("{}", "=".repeat(60));
    eprintln!("\nTo use the association testing feature, you need:");
    eprintln!("\n1. R installed and available in your PATH");
    eprintln!("   - Ubuntu/Debian: sudo apt-get install r-base");
    eprintln!("   - CentOS/RHEL: sudo yum install R");
    eprintln!("   - macOS with Homebrew: brew install r");
    eprintln!("   - Windows: Download from https://cran.r-project.org/");

    eprintln!("\n2. Required R packages installed:");
    eprintln!("   Start R and run:");
    eprintln!("   install.packages(c('data.table', 'argparser'))");
    eprintln!("   # Note: 'parallel' is part of base R");

    eprintln!("\n3. Alternative: Install using command line:");
    eprintln!("   Rscript -e \"install.packages(c('data.table', 'argparser'), repos='https://cran.rstudio.com/')\"");

    eprintln!("\nOnce R is properly set up, you can run:");
    eprintln!("   inquiSTR association --help");
    eprintln!("\nFor more details, see: https://github.com/wdecoster/inquiSTR#association-testing");
    eprintln!("{}", "=".repeat(60));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_arguments_valid() {
        // Create temporary test files
        let temp_dir = std::env::temp_dir();
        let input_file = temp_dir.join("test_input.txt");
        let phenocovar_file = temp_dir.join("test_pheno.txt");

        // Create the files
        std::fs::write(&input_file, "test").unwrap();
        std::fs::write(&phenocovar_file, "test").unwrap();

        let result = validate_arguments(&input_file, &phenocovar_file, "MAX", "continuous", &None);

        assert!(result.is_ok());

        // Clean up
        std::fs::remove_file(&input_file).unwrap();
        std::fs::remove_file(&phenocovar_file).unwrap();
    }

    #[test]
    fn test_validate_arguments_invalid_str_mode() {
        let temp_dir = std::env::temp_dir();
        let input_file = temp_dir.join("test_input2.txt");
        let phenocovar_file = temp_dir.join("test_pheno2.txt");

        std::fs::write(&input_file, "test").unwrap();
        std::fs::write(&phenocovar_file, "test").unwrap();

        let result =
            validate_arguments(&input_file, &phenocovar_file, "INVALID", "continuous", &None);

        assert!(result.is_err());
        assert!(result.unwrap_err().contains("STR mode must be"));

        std::fs::remove_file(&input_file).unwrap();
        std::fs::remove_file(&phenocovar_file).unwrap();
    }

    #[test]
    fn test_validate_arguments_binary_missing_order() {
        let temp_dir = std::env::temp_dir();
        let input_file = temp_dir.join("test_input3.txt");
        let phenocovar_file = temp_dir.join("test_pheno3.txt");

        std::fs::write(&input_file, "test").unwrap();
        std::fs::write(&phenocovar_file, "test").unwrap();

        let result = validate_arguments(&input_file, &phenocovar_file, "MAX", "binary", &None);

        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Binary order is required"));

        std::fs::remove_file(&input_file).unwrap();
        std::fs::remove_file(&phenocovar_file).unwrap();
    }

    #[test]
    fn test_create_temp_r_script() {
        let result = create_temp_r_script();
        assert!(result.is_ok());

        let temp_file = result.unwrap();
        assert!(temp_file.exists());

        // Check content
        let content = std::fs::read_to_string(&temp_file).unwrap();
        assert!(content.contains("#!/usr/bin/env Rscript"));
        assert!(content.contains("STR_regression"));

        // Clean up
        std::fs::remove_file(&temp_file).unwrap();
    }
}
