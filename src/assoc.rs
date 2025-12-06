//! Association testing module
//!
//! This module provides an interface to perform STR association testing by wrapping
//! the R script STR_regression.R. It handles R environment detection, dependency checking,
//! and execution of the statistical analysis.

use log::{error, info, warn};
use std::env;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

// Embed the R script directly in the binary
const STR_REGRESSION_SCRIPT: &str = include_str!("../scripts/STR_regression.R");

// Required R packages (qqman only needed for --plot)
const REQUIRED_R_PACKAGES: &[&str] = &["data.table", "argparser", "parallel"];

// Portable R installation constants
const MINIFORGE_URL: &str = "https://github.com/conda-forge/miniforge/releases/download/24.11.0-0/Miniforge3-24.11.0-0-Linux-x86_64.sh";
const R_VERSION: &str = "4.4.1";

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
    plot: Option<String>,
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
        print_error_help_message();
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
        plot.as_deref(),
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
            print_error_help_message();
            std::process::exit(1);
        }
    }
}

/// Print helpful error message directing users to report issues
fn print_error_help_message() {
    eprintln!("\n{}", "=".repeat(60));
    eprintln!("  ASSOCIATION TESTING ERROR");
    eprintln!("{}", "=".repeat(60));
    eprintln!("\nIf you're experiencing issues with association testing:");
    eprintln!("\n1. Check the error message above for specific details");
    eprintln!("2. Verify your phenotype file format matches the documentation");
    eprintln!("3. Ensure your R environment is properly configured");
    eprintln!("\nIf the problem persists, please open an issue on GitHub:");
    eprintln!("  https://github.com/wdecoster/inquiSTR/issues/new");
    eprintln!("\nInclude in your report:");
    eprintln!("  - The full error message above");
    eprintln!("  - Your inquiSTR version: {}", env!("CARGO_PKG_VERSION"));
    eprintln!("  - Your R version (run: R --version)");
    eprintln!("  - Your operating system");
    eprintln!("  - The command you ran (you can redact sensitive file paths)");
    eprintln!("{}", "=".repeat(60));
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
/// Uses a fallback cascade: system R -> conda R -> cached portable R -> offer installation
fn check_r_environment(quiet: bool) -> Result<(), String> {
    if !quiet {
        info!("Checking R environment...");
    }

    // Step 1: Check for system R
    if let Ok(version) = check_system_r() {
        if !quiet {
            info!("System R found: {}", version);
        }
        // Check for required packages and offer to install if missing
        ensure_r_packages(quiet)?;
        return Ok(());
    }

    // Step 2: Check for system conda with inquiSTR environment
    match check_system_conda_r() {
        Ok((r_path, version, _conda_cmd)) => {
            if !quiet {
                info!("System conda R found in inquiSTR environment: {}", version);
                info!("R path: {}", r_path);
            }
            return Ok(());
        }
        Err(e) if e.contains("found but no 'inquiSTR' environment") => {
            // Conda exists but no inquiSTR env - offer to create it
            let conda_cmd = e.split_whitespace().next().unwrap_or("conda");
            return offer_conda_env_creation(conda_cmd, quiet);
        }
        Err(_) => {
            // No conda found, continue to next fallback
        }
    }

    // Step 3: Check for cached portable R
    let cache_dir = get_cache_dir();
    let conda_prefix = cache_dir.join("miniforge");
    let r_bin = conda_prefix.join("bin/R");

    if r_bin.exists() {
        if !quiet {
            info!("Cached portable R found at: {}", conda_prefix.display());
        }
        return Ok(());
    }

    // Step 4: No R found - offer to install portable R
    offer_portable_r_installation(&conda_prefix, quiet)
}

/// Ensure R packages are installed, offering to install missing ones
fn ensure_r_packages(quiet: bool) -> Result<(), String> {
    let missing_packages = check_r_packages(quiet)?;

    if missing_packages.is_empty() {
        return Ok(());
    }

    // Packages are missing - offer to install
    eprintln!("\n{}", "=".repeat(60));
    eprintln!("  Missing R packages: {}", missing_packages.join(", "));
    eprintln!("{}", "=".repeat(60));
    eprintln!("\nThe association testing feature requires the following R packages:");
    for pkg in &missing_packages {
        eprintln!("  - {}", pkg);
    }
    eprintln!("\nWould you like to install them now? This may take 1-2 minutes.");
    eprintln!("Type 'y' to install, or any other key to cancel: ");

    let mut response = String::new();
    io::stdin()
        .read_line(&mut response)
        .map_err(|e| format!("Failed to read user input: {}", e))?;

    let response = response.trim().to_lowercase();

    if response != "y" {
        eprintln!("\nInstallation cancelled. You can install the packages manually:");
        eprintln!(
            "  Rscript -e \"install.packages(c('{}'), repos='https://cran.rstudio.com/')\"",
            missing_packages.join("', '")
        );
        return Err("R package installation cancelled by user".to_string());
    }

    // User confirmed - proceed with installation
    install_r_packages(&missing_packages)?;

    // Verify installation succeeded
    let still_missing = check_r_packages(true)?;
    if !still_missing.is_empty() {
        return Err(format!(
            "Installation completed but packages still missing: {}. Please try installing manually.",
            still_missing.join(", ")
        ));
    }

    eprintln!("\n✓ R packages installed successfully!\n");
    Ok(())
}

/// Install R packages automatically
fn install_r_packages(packages: &[String]) -> Result<(), String> {
    eprintln!("\nInstalling R packages...");
    eprintln!("This may take 1-2 minutes depending on your system.");

    let pkg_list = packages.join("', '");
    let install_cmd = format!(
        "install.packages(c('{}'), repos='https://cran.rstudio.com/', quiet=FALSE)",
        pkg_list
    );

    let output = Command::new("Rscript")
        .args(["-e", &install_cmd])
        .output()
        .map_err(|e| format!("Failed to run R package installation: {}", e))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);

        eprintln!("\nR package installation failed:");
        if !stdout.is_empty() {
            eprintln!("STDOUT: {}", stdout);
        }
        if !stderr.is_empty() {
            eprintln!("STDERR: {}", stderr);
        }

        return Err(format!(
            "Failed to install R packages. Exit code: {:?}. \
             You may need to install them manually or check your R installation.",
            output.status.code()
        ));
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
    plot: Option<&str>,
) -> Result<(), String> {
    let script_path_str = script_path.to_string_lossy();
    let input_str = input.to_string_lossy();
    let phenocovar_str = phenocovar.to_string_lossy();
    let out_str = out.to_string_lossy();
    let missing_cutoff_str = missing_cutoff.to_string();
    let threads_str = threads.to_string();
    let chunk_size_str = chunk_size.to_string();

    let mut script_args = vec![
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
        script_args.extend(&["--covnames", covnames_val]);
    }

    let minimal_length_str = minimal_length.map(|val| val.to_string());
    if let Some(ref minimal_length_val) = minimal_length_str {
        script_args.extend(&["--minimal_length", minimal_length_val]);
    }

    if let Some(ref binary_order_val) = binary_order {
        script_args.extend(&["--binaryOrder", binary_order_val]);
    }

    if quiet {
        script_args.push("--quiet");
    }

    let plot_prefix_str;
    if let Some(plot_prefix) = plot {
        script_args.push("--plot");
        plot_prefix_str = plot_prefix.to_string();
        script_args.push(&plot_prefix_str);
    }

    // Determine which R to use based on detection cascade
    let (r_command, r_args) = get_r_command(script_path_str.as_ref(), &script_args)?;

    info!("Executing R script: {} {:?}", r_command, r_args);

    let output = Command::new(&r_command)
        .args(&r_args)
        .output()
        .map_err(|e| format!("Failed to execute R script: {}", e))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);

        eprintln!("\n{}", "=".repeat(60));
        eprintln!("R SCRIPT EXECUTION FAILED");
        eprintln!("{}", "=".repeat(60));

        if !stderr.is_empty() {
            eprintln!("\nError details:");
            eprintln!("{}", stderr);
        }

        if !stdout.is_empty() {
            eprintln!("\nOutput before failure:");
            eprintln!("{}", stdout);
        }

        eprintln!("\n{}", "=".repeat(60));

        return Err(format!(
            "R script failed with exit code: {}",
            output
                .status
                .code()
                .map_or("unknown".to_string(), |c| c.to_string())
        ));
    }

    // Always show stderr (warnings, plot messages, etc.) even on success
    let stderr = String::from_utf8_lossy(&output.stderr);
    if !stderr.is_empty() && !quiet {
        eprintln!("{}", stderr);
    }

    if !quiet {
        let stdout = String::from_utf8_lossy(&output.stdout);
        if !stdout.is_empty() {
            println!("{}", stdout);
        }
    }

    Ok(())
}

/// Get the appropriate R command and arguments based on which R installation is available
fn get_r_command(script_path: &str, script_args: &[&str]) -> Result<(String, Vec<String>), String> {
    // Try system R first
    if Command::new("Rscript").arg("--version").output().is_ok() {
        let mut args = vec![script_path.to_string()];
        args.extend(script_args.iter().map(|s| s.to_string()));
        return Ok(("Rscript".to_string(), args));
    }

    // Try conda R in inquiSTR environment
    for conda_cmd in &["conda", "mamba", "micromamba"] {
        let check = Command::new(conda_cmd)
            .args(["run", "-n", "inquiSTR", "Rscript", "--version"])
            .output();

        if let Ok(output) = check
            && output.status.success()
        {
            let mut args = vec![
                "run".to_string(),
                "-n".to_string(),
                "inquiSTR".to_string(),
                "Rscript".to_string(),
                script_path.to_string(),
            ];
            args.extend(script_args.iter().map(|s| s.to_string()));
            return Ok((conda_cmd.to_string(), args));
        }
    }

    // Try cached portable R
    let cache_dir = get_cache_dir();
    let r_script = cache_dir.join("miniforge/bin/Rscript");
    if r_script.exists() {
        let mut args = vec![script_path.to_string()];
        args.extend(script_args.iter().map(|s| s.to_string()));
        return Ok((r_script.to_string_lossy().to_string(), args));
    }

    Err("No R installation found. This should not happen after check_r_environment".to_string())
}

/// Check for system R installation
fn check_system_r() -> Result<String, String> {
    let output = Command::new("R")
        .args(["--version"])
        .output()
        .map_err(|_| "R not found in PATH".to_string())?;

    if !output.status.success() {
        return Err("R command failed".to_string());
    }

    let version_output = String::from_utf8_lossy(&output.stdout);
    let version = version_output
        .lines()
        .next()
        .unwrap_or("unknown version")
        .to_string();

    Ok(version)
}

/// Check for system conda/mamba with inquiSTR environment
fn check_system_conda_r() -> Result<(String, String, String), String> {
    let conda_commands = ["conda", "mamba", "micromamba"];

    for conda_cmd in &conda_commands {
        // Check if conda command exists
        let conda_check = Command::new(conda_cmd).arg("--version").output();

        if conda_check.is_err() {
            continue;
        }

        // Check if inquiSTR environment exists
        let env_check = Command::new(conda_cmd).args(["env", "list"]).output();

        let has_inquistr_env = if let Ok(output) = env_check {
            String::from_utf8_lossy(&output.stdout)
                .lines()
                .any(|line| line.split_whitespace().next() == Some("inquiSTR"))
        } else {
            false
        };

        if !has_inquistr_env {
            return Err(format!("{} found but no 'inquiSTR' environment", conda_cmd));
        }

        // Try to run R from the inquiSTR environment
        let r_check = Command::new(conda_cmd)
            .args(["run", "-n", "inquiSTR", "R", "--version"])
            .output();

        if let Ok(output) = r_check
            && output.status.success()
        {
            let version_info = String::from_utf8_lossy(&output.stdout);
            let version = version_info.lines().next().unwrap_or("unknown").to_string();

            // Check for required packages
            let pkg_check = Command::new(conda_cmd)
                .args(["run", "-n", "inquiSTR", "R", "-e", 
                       "if (!require('data.table', quietly=TRUE)) quit(status=1); if (!require('argparser', quietly=TRUE)) quit(status=1)"])
                .output();

            if let Ok(pkg_output) = pkg_check
                && pkg_output.status.success()
            {
                let r_path = format!("{} run -n inquiSTR R", conda_cmd);
                return Ok((
                    r_path,
                    format!("{} (conda env: inquiSTR)", version),
                    conda_cmd.to_string(),
                ));
            }
        }
    }

    Err("No conda/mamba/micromamba found in PATH".to_string())
}

/// Get cache directory for portable R installation
fn get_cache_dir() -> PathBuf {
    if let Ok(xdg_cache) = env::var("XDG_CACHE_HOME") {
        PathBuf::from(xdg_cache).join("inquiSTR")
    } else if let Some(home) = env::var_os("HOME") {
        PathBuf::from(home).join(".cache/inquiSTR")
    } else {
        PathBuf::from("/tmp/inquiSTR")
    }
}

/// Offer to create inquiSTR conda environment in existing conda installation
fn offer_conda_env_creation(conda_cmd: &str, quiet: bool) -> Result<(), String> {
    eprintln!("\n{}", "=".repeat(60));
    eprintln!("  CONDA ENVIRONMENT SETUP");
    eprintln!("{}", "=".repeat(60));
    eprintln!("\nFound {} but no 'inquiSTR' environment.", conda_cmd);
    eprintln!("\nWould you like to create an inquiSTR environment with R and required packages?");
    eprintln!(
        "This will install R {} with data.table, argparser, and qqman (for plotting) (~250 MB)",
        R_VERSION
    );
    eprintln!("\nType 'y' to create environment, or any other key to cancel: ");

    let mut response = String::new();
    io::stdin()
        .read_line(&mut response)
        .map_err(|e| format!("Failed to read user input: {}", e))?;

    if response.trim().to_lowercase() != "y" {
        eprintln!("\nCancelled. You can create the environment later with:");
        eprintln!(
            "  {} create -n inquiSTR -c conda-forge r-base={} r-data.table r-argparser r-qqman",
            conda_cmd, R_VERSION
        );
        return Err("User cancelled conda environment creation".to_string());
    }

    // Create the environment
    if !quiet {
        info!("Creating inquiSTR conda environment...");
    }
    eprintln!("\nCreating inquiSTR conda environment...");
    eprintln!("This may take 60-90 seconds...");

    let create_status = Command::new(conda_cmd)
        .args([
            "create",
            "-n",
            "inquiSTR",
            "-y",
            "-q",
            "-c",
            "conda-forge",
            &format!("r-base={}", R_VERSION),
            "r-data.table",
            "r-argparser",
            "r-qqman",
        ])
        .status()
        .map_err(|e| format!("Failed to create environment: {}", e))?;

    if !create_status.success() {
        return Err("Environment creation failed".to_string());
    }

    eprintln!("\n✓ inquiSTR conda environment created successfully!");
    if !quiet {
        info!("inquiSTR conda environment created");
    }

    Ok(())
}

/// Offer to install portable R via miniforge
fn offer_portable_r_installation(conda_prefix: &Path, quiet: bool) -> Result<(), String> {
    eprintln!("\n{}", "=".repeat(60));
    eprintln!("  PORTABLE R INSTALLATION");
    eprintln!("{}", "=".repeat(60));
    eprintln!("\nNo R installation found.");
    eprintln!("\ninquiSTR can download a self-contained R environment using conda.");
    eprintln!("This works across all Linux distributions without system dependencies.");
    eprintln!("\nDetails:");
    eprintln!("  - Download size: ~150 MB (miniforge) + ~200 MB (R + packages)");
    eprintln!("  - Total disk space: ~350 MB");
    eprintln!("  - Install location: {}", conda_prefix.display());
    eprintln!("  - No system changes required");
    eprintln!(
        "  - Includes R {} with data.table, argparser, and qqman (for plotting)",
        R_VERSION
    );
    eprintln!("  - One-time download, ~2-3 minutes");
    eprintln!("\nWould you like to download and install portable R now?");
    eprintln!("Type 'y' to proceed, or any other key to cancel: ");

    let mut response = String::new();
    io::stdin()
        .read_line(&mut response)
        .map_err(|e| format!("Failed to read user input: {}", e))?;

    if response.trim().to_lowercase() != "y" {
        eprintln!("\nInstallation cancelled.");
        print_r_setup_instructions();
        return Err("User cancelled portable R installation".to_string());
    }

    // Install miniforge and R
    install_portable_r(conda_prefix, quiet)
}

/// Install portable R via miniforge
fn install_portable_r(conda_prefix: &Path, quiet: bool) -> Result<(), String> {
    if !quiet {
        info!("Installing portable R environment...");
    }

    let temp_dir = env::temp_dir();
    let installer_path = temp_dir.join("miniforge-installer.sh");

    // Step 1: Download miniforge
    eprintln!("\n[1/4] Downloading miniforge installer...");
    eprintln!("      URL: {}", MINIFORGE_URL);

    let download_status = Command::new("curl")
        .args([
            "-L",
            "-f",
            "-o",
            installer_path.to_str().unwrap(),
            "--progress-bar",
            MINIFORGE_URL,
        ])
        .status();

    if download_status.is_err() || !download_status.unwrap().success() {
        eprintln!("      curl failed, trying wget...");
        let wget_status = Command::new("wget")
            .args([
                "-O",
                installer_path.to_str().unwrap(),
                "--show-progress",
                MINIFORGE_URL,
            ])
            .status()
            .map_err(|_| "Failed to download: curl and wget both unavailable".to_string())?;

        if !wget_status.success() {
            return Err("Download failed".to_string());
        }
    }

    eprintln!(
        "      ✓ Downloaded ({:.1} MB)",
        fs::metadata(&installer_path).unwrap().len() as f64 / 1_000_000.0
    );

    // Step 2: Install miniforge
    eprintln!("\n[2/4] Installing miniforge to {}...", conda_prefix.display());
    eprintln!("      This may take 30-60 seconds...");

    let install_status = Command::new("bash")
        .args([
            installer_path.to_str().unwrap(),
            "-b",
            "-p",
            conda_prefix.to_str().unwrap(),
        ])
        .status()
        .map_err(|e| format!("Failed to run miniforge installer: {}", e))?;

    if !install_status.success() {
        return Err("Miniforge installation failed".to_string());
    }

    eprintln!("      ✓ Miniforge installed");

    // Step 3: Install R and packages
    eprintln!("\n[3/4] Installing R {} and packages...", R_VERSION);
    eprintln!("      This may take 60-90 seconds...");

    let conda_bin = conda_prefix.join("bin/conda");
    let install_r_status = Command::new(conda_bin.to_str().unwrap())
        .args([
            "install",
            "-y",
            "-q",
            "-c",
            "conda-forge",
            &format!("r-base={}", R_VERSION),
            "r-data.table",
            "r-argparser",
            "r-qqman",
        ])
        .env("CONDA_PREFIX", conda_prefix)
        .status()
        .map_err(|e| format!("Failed to install R packages: {}", e))?;

    if !install_r_status.success() {
        return Err("R package installation failed".to_string());
    }

    eprintln!("      ✓ R and packages installed");

    // Step 4: Verify installation
    eprintln!("\n[4/4] Verifying installation...");

    let r_bin = conda_prefix.join("bin/R");
    if !r_bin.exists() {
        return Err(format!("R binary not found at: {}", r_bin.display()));
    }

    let test_status = Command::new(r_bin.to_str().unwrap())
        .args(["--slave", "-e", "cat('OK')"])
        .output()
        .map_err(|e| format!("Failed to test R: {}", e))?;

    if !test_status.status.success() {
        return Err("R execution test failed".to_string());
    }

    eprintln!("      ✓ R is working correctly");

    // Clean up installer
    let _ = fs::remove_file(&installer_path);

    eprintln!("\n✓ Portable R installation complete!");
    eprintln!("  Location: {}", conda_prefix.display());
    eprintln!("  R binary: {}", r_bin.display());

    if !quiet {
        info!("Portable R installation complete");
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
    eprintln!("   # Optional: install.packages('qqman') for --plot functionality");

    eprintln!("\n3. Alternative: Install using command line:");
    eprintln!(
        "   Rscript -e \"install.packages(c('data.table', 'argparser'), repos='https://cran.rstudio.com/')\""
    );
    eprintln!(
        "   # Optional for plotting: Rscript -e \"install.packages('qqman', repos='https://cran.rstudio.com/')\""
    );

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
