use rayon::prelude::*;
use std::io::{BufRead, BufWriter, Write};
use std::path::PathBuf;

/// Run native Rust association testing for STR or kmer frequency data.
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
    binary_order: Option<String>,
    quiet: bool,
    sort: bool,
    plot: Option<PathBuf>,
) {
    // Validate arguments
    if !input.exists() {
        eprintln!("ERROR: Input file does not exist: {}", input.display());
        std::process::exit(1);
    }
    if !phenocovar.exists() {
        eprintln!("ERROR: Phenotype file does not exist: {}", phenocovar.display());
        std::process::exit(1);
    }
    if !(0.0..=1.0).contains(&missing_cutoff) {
        eprintln!("ERROR: missing_cutoff must be between 0 and 1");
        std::process::exit(1);
    }
    if outcometype != "binary" && outcometype != "continuous" {
        eprintln!("ERROR: outcometype must be 'binary' or 'continuous'");
        std::process::exit(1);
    }
    if outcometype == "binary" && binary_order.is_none() {
        eprintln!("ERROR: --binary-order is required for binary outcomes");
        std::process::exit(1);
    }

    // Detect input type
    let input_type = detect_input_type(&input);

    if input_type == InputType::Str && !["MEAN", "MAX", "MIN"].contains(&str_mode.as_str()) {
        eprintln!("ERROR: --str-mode must be MEAN, MAX, or MIN for STR data");
        std::process::exit(1);
    }

    // Parse binary levels
    let binary_levels: Option<Vec<String>> = binary_order.as_ref().map(|bo| {
        bo.split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    });

    // Prepare phenotype data
    let pheno =
        prepare_phenotype_data(&phenocovar, &phenotype, &covnames, binary_levels.as_deref(), quiet);

    if !quiet {
        eprintln!("Loaded {} samples with phenotype data", pheno.sample_ids.len());
    }

    // Parse header and extract sample names
    let file = crate::utils::reader(&input.to_string_lossy());
    let mut lines = file.lines();
    let header_line = crate::utils::skip_metadata_lines(&mut lines, &input.to_string_lossy());
    let (sample_names, input_type_confirmed) = parse_header(&header_line, &input_type);

    if !quiet {
        eprintln!("Found {} samples in input file", sample_names.len());
        eprintln!("Input type: {}", input_type_confirmed);
    }

    // Build sample index mapping: input column position -> pheno row index
    let sample_mapping: Vec<Option<usize>> = sample_names
        .iter()
        .map(|name| pheno.sample_ids.iter().position(|s| s == name))
        .collect();

    // Count overlapping samples
    let overlap_count = sample_mapping.iter().filter(|m| m.is_some()).count();
    if overlap_count == 0 {
        eprintln!("ERROR: No overlapping samples between input file and phenotype file.");
        if sample_names.len() <= 10 {
            eprintln!("  Sample IDs in input file: {:?}", sample_names);
        } else {
            eprintln!("  First 10 sample IDs in input file: {:?}", &sample_names[..10]);
            eprintln!("  ({} total)", sample_names.len());
        }
        if pheno.sample_ids.len() <= 10 {
            eprintln!("  Sample IDs in phenotype file: {:?}", pheno.sample_ids);
        } else {
            eprintln!("  First 10 sample IDs in phenotype file: {:?}", &pheno.sample_ids[..10]);
            eprintln!("  ({} total)", pheno.sample_ids.len());
        }
        std::process::exit(1);
    }
    if !quiet {
        eprintln!("{} samples overlap between input and phenotype files", overlap_count);
    }

    // Pre-compute total phenotyped count (constant across variants)
    let total_phenotyped = overlap_count;

    // Read all variant lines
    let variant_lines: Vec<String> = lines.map_while(Result::ok).collect();
    let total_variants = variant_lines.len();

    if !quiet {
        eprintln!("Processing {} variants...", total_variants);
    }

    // Set up thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Failed to build thread pool");

    // Process variants in parallel
    let results: Vec<Option<VariantResult>> = pool.install(|| {
        variant_lines
            .par_iter()
            .map(|line| {
                process_variant_line(
                    line,
                    &sample_names,
                    &sample_mapping,
                    &pheno,
                    &str_mode,
                    &outcometype,
                    missing_cutoff,
                    minimal_length,
                    &input_type_confirmed,
                    total_phenotyped,
                )
            })
            .collect()
    });

    // Collect non-None results
    let mut valid_results: Vec<VariantResult> = results.into_iter().flatten().collect();
    let variants_written = valid_results.len();

    if !quiet {
        let pass_rate = if total_variants > 0 {
            (variants_written as f64 / total_variants as f64 * 100.0) as u32
        } else {
            0
        };
        eprintln!(
            "Processed {} variants, {} passed filters ({}% pass rate)",
            total_variants, variants_written, pass_rate
        );
    }

    // Add Bonferroni correction
    let n_tests = variants_written;
    for r in &mut valid_results {
        r.pvalue_bonf = Some((r.pvalue * n_tests as f64).min(1.0));
    }

    // Sort if requested
    if sort {
        valid_results.sort_by(|a, b| {
            a.pvalue
                .partial_cmp(&b.pvalue)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    // Write output
    write_results(&out, &valid_results, &outcometype, binary_levels.as_deref(), quiet);

    if !quiet {
        let bonf_sig = valid_results
            .iter()
            .filter(|r| r.pvalue_bonf.unwrap_or(1.0) < 0.05)
            .count();
        eprintln!("Bonferroni significant (p < 0.05): {}", bonf_sig);
        eprintln!("Results written to: {}", out.display());
    }

    // Generate Manhattan plot if requested
    if let Some(plot_path) = plot {
        write_manhattan_plot(&valid_results, &plot_path, quiet);
    }
}

// ---------------------------------------------------------------------------
// Input type detection
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
enum InputType {
    Str,
    Kmer,
}

impl std::fmt::Display for InputType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InputType::Str => write!(f, "STR"),
            InputType::Kmer => write!(f, "kmer"),
        }
    }
}

fn detect_input_type(path: &std::path::Path) -> InputType {
    let file = crate::utils::reader(&path.to_string_lossy());
    for line_result in file.lines() {
        let line = line_result.expect("Error reading input file");
        if line.starts_with('#') {
            if line.contains("file_type=kmer") || line.contains("file_type = kmer") {
                return InputType::Kmer;
            }
            if line.contains("file_type=call") || line.contains("file_type = call") {
                return InputType::Str;
            }
            continue;
        }
        // First non-metadata line is the header
        let first_field = line.split('\t').next().unwrap_or("");
        return match first_field {
            "kmer" => InputType::Kmer,
            "chromosome" => InputType::Str,
            _ => {
                eprintln!(
                    "ERROR: Unrecognized file format. Expected 'kmer' or 'chromosome' as first column, got: {}",
                    first_field
                );
                std::process::exit(1);
            }
        };
    }
    eprintln!("ERROR: Empty input file");
    std::process::exit(1);
}

// ---------------------------------------------------------------------------
// Header parsing
// ---------------------------------------------------------------------------

fn parse_header(header_line: &str, _expected_type: &InputType) -> (Vec<String>, InputType) {
    let fields: Vec<&str> = header_line.trim().split('\t').collect();

    match fields.first().copied() {
        Some("kmer") => {
            let sample_names: Vec<String> = fields[1..].iter().map(|s| s.to_string()).collect();
            (sample_names, InputType::Kmer)
        }
        Some("chromosome") => {
            if fields.len() < 5 {
                eprintln!("ERROR: STR file must have at least 5 columns");
                std::process::exit(1);
            }
            if fields[1] != "begin" || fields[2] != "end" {
                eprintln!("ERROR: Expected 'begin' and 'end' as columns 2 and 3");
                std::process::exit(1);
            }
            let sample_cols = &fields[4..];
            if !sample_cols.len().is_multiple_of(2) {
                eprintln!("ERROR: Must have even number of sample columns (H1/H2 pairs)");
                std::process::exit(1);
            }
            let n_samples = sample_cols.len() / 2;
            let mut sample_names = Vec::with_capacity(n_samples);
            for i in 0..n_samples {
                let h1_col = sample_cols[i * 2];
                let h2_col = sample_cols[i * 2 + 1];
                if !h1_col.ends_with("_H1") || !h2_col.ends_with("_H2") {
                    eprintln!(
                        "ERROR: Expected _H1/_H2 suffixed columns, got: {} / {}",
                        h1_col, h2_col
                    );
                    std::process::exit(1);
                }
                let name_h1 = h1_col.strip_suffix("_H1").unwrap();
                let name_h2 = h2_col.strip_suffix("_H2").unwrap();
                if name_h1 != name_h2 {
                    eprintln!("ERROR: H1/H2 pair mismatch: {} vs {}", name_h1, name_h2);
                    std::process::exit(1);
                }
                sample_names.push(name_h1.to_string());
            }
            (sample_names, InputType::Str)
        }
        _ => {
            eprintln!(
                "ERROR: Expected 'kmer' or 'chromosome' as first column, got: {}",
                fields.first().unwrap_or(&"<empty>")
            );
            std::process::exit(1);
        }
    }
}

// ---------------------------------------------------------------------------
// Phenotype data
// ---------------------------------------------------------------------------

struct PhenotypeData {
    /// Sample IDs in order
    sample_ids: Vec<String>,
    /// Phenotype values (0/1 for binary, continuous for continuous)
    phenotype_values: Vec<f64>,
    /// Covariate matrix: rows = samples, cols = expanded covariates
    /// Categorical covariates are dummy-encoded (N-1 columns per N levels)
    covariate_matrix: Vec<Vec<f64>>,
    /// Expanded covariate column names (e.g. "sex_M", "pop_EUR", "pop_AFR")
    covariate_names: Vec<String>,
    /// For binary: the two level names
    binary_levels: Option<Vec<String>>,
}

fn prepare_phenotype_data(
    phenocovar: &std::path::Path,
    phenotype_col: &str,
    covnames: &Option<String>,
    binary_levels: Option<&[String]>,
    quiet: bool,
) -> PhenotypeData {
    let file = crate::utils::reader(&phenocovar.to_string_lossy());
    let mut lines_iter = file.lines();

    // Read header
    let header = lines_iter
        .next()
        .expect("Phenotype file is empty")
        .expect("Error reading phenotype file");
    let header_fields: Vec<String> = header
        .trim()
        .split('\t')
        .map(|s| s.trim().to_string())
        .collect();

    if header_fields.len() < 2 {
        eprintln!("ERROR: Phenotype file must have at least 2 columns (sample ID + phenotype)");
        std::process::exit(1);
    }

    let _sample_id_col_name = &header_fields[0];
    let pheno_col_idx = header_fields
        .iter()
        .position(|h| h == phenotype_col)
        .unwrap_or_else(|| {
            eprintln!(
                "ERROR: Phenotype column '{}' not found in phenotype file. Available columns: {}",
                phenotype_col,
                header_fields.join(", ")
            );
            std::process::exit(1);
        });

    // Parse covariate names
    let cov_names: Vec<String> = covnames
        .as_ref()
        .map(|c| {
            c.split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect()
        })
        .unwrap_or_default();

    let cov_indices: Vec<usize> = cov_names
        .iter()
        .map(|name| {
            header_fields.iter().position(|h| h == name).unwrap_or_else(|| {
                eprintln!(
                    "ERROR: Covariate column '{}' not found in phenotype file. Available columns: {}",
                    name,
                    header_fields.join(", ")
                );
                std::process::exit(1);
            })
        })
        .collect();

    let mut sample_ids = Vec::new();
    let mut pheno_values = Vec::new();
    let mut raw_covariates: Vec<Vec<String>> = Vec::new(); // rows of raw string values

    for (line_num, line_result) in lines_iter.enumerate() {
        let line = line_result.expect("Error reading phenotype file");
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() <= pheno_col_idx {
            eprintln!(
                "ERROR: Line {} in phenotype file has {} columns, but phenotype column '{}' is at index {}",
                line_num + 2,
                fields.len(),
                phenotype_col,
                pheno_col_idx + 1
            );
            std::process::exit(1);
        }

        let sample_id = fields[0].trim().to_string();
        let pheno_str = fields[pheno_col_idx].trim();

        // Skip samples with missing phenotype (NA/NaN/empty) — this is expected
        if pheno_str.is_empty()
            || pheno_str.eq_ignore_ascii_case("na")
            || pheno_str.eq_ignore_ascii_case("nan")
        {
            continue;
        }

        // For binary: filter to only samples in binary_levels, encode as 0/1
        if let Some(levels) = binary_levels {
            if !levels.iter().any(|l| l == pheno_str) {
                eprintln!(
                    "ERROR: Sample '{}' has phenotype value '{}' which is not in --binary-order '{}'. \
                     All non-missing phenotype values must match one of the specified binary levels.",
                    sample_id,
                    pheno_str,
                    levels.join(",")
                );
                std::process::exit(1);
            }
            let pheno_val = if pheno_str == levels[0] { 0.0 } else { 1.0 };
            pheno_values.push(pheno_val);
        } else {
            // Continuous phenotype
            let pheno_val: f64 = match pheno_str.parse() {
                Ok(v) => v,
                Err(_) => {
                    eprintln!(
                        "ERROR: Sample '{}' has non-numeric phenotype value '{}' for continuous outcome",
                        sample_id, pheno_str
                    );
                    std::process::exit(1);
                }
            };
            pheno_values.push(pheno_val);
        }

        // Collect raw covariate strings
        let mut cov_row = Vec::with_capacity(cov_indices.len());
        for (cov_i, &ci) in cov_indices.iter().enumerate() {
            if ci >= fields.len() {
                eprintln!(
                    "ERROR: Sample '{}' is missing covariate column '{}' (line has {} columns, need column {})",
                    sample_id,
                    cov_names[cov_i],
                    fields.len(),
                    ci + 1
                );
                std::process::exit(1);
            }
            let val_str = fields[ci].trim();
            if val_str.is_empty()
                || val_str.eq_ignore_ascii_case("na")
                || val_str.eq_ignore_ascii_case("nan")
            {
                eprintln!(
                    "ERROR: Sample '{}' has missing value for covariate '{}'",
                    sample_id, cov_names[cov_i]
                );
                std::process::exit(1);
            }
            cov_row.push(val_str.to_string());
        }

        sample_ids.push(sample_id);
        raw_covariates.push(cov_row);
    }

    if sample_ids.is_empty() {
        // Re-read the file to collect all sample IDs (first column) for diagnostics
        let diag_file = crate::utils::reader(&phenocovar.to_string_lossy());
        let mut diag_lines = diag_file.lines();
        diag_lines.next(); // skip header
        let all_ids: Vec<String> = diag_lines
            .map_while(Result::ok)
            .filter_map(|l| l.split('\t').next().map(|s| s.trim().to_string()))
            .filter(|s| !s.is_empty())
            .collect();
        eprintln!("ERROR: No valid samples found in phenotype file after filtering.");
        eprintln!("  Possible causes:");
        eprintln!("    - All phenotype values are missing (NA/NaN/empty)");
        eprintln!("    - The phenotype file may not be tab-separated");
        if all_ids.len() <= 10 {
            eprintln!("  Sample IDs in phenotype file: {:?}", all_ids);
        } else {
            eprintln!("  First 10 sample IDs in phenotype file: {:?}", &all_ids[..10]);
            eprintln!("  ({} total)", all_ids.len());
        }
        std::process::exit(1);
    }

    // Validate binary groups
    if let Some(levels) = binary_levels {
        let n_group0 = pheno_values.iter().filter(|&&v| v == 0.0).count();
        let n_group1 = pheno_values.iter().filter(|&&v| v == 1.0).count();
        if n_group0 < 2 || n_group1 < 2 {
            eprintln!(
                "ERROR: Insufficient samples in binary groups. {} has {} samples, {} has {} samples. Need at least 2 per group.",
                levels[0], n_group0, levels[1], n_group1
            );
            std::process::exit(1);
        }
        if !quiet {
            eprintln!(
                "Sample counts per group: {} = {}, {} = {}",
                levels[0], n_group0, levels[1], n_group1
            );
        }
    }

    // Encode covariates: numeric covariates stay as-is, categorical get dummy-encoded
    let (cov_matrix, expanded_cov_names) = encode_covariates(&raw_covariates, &cov_names, quiet);

    PhenotypeData {
        sample_ids,
        phenotype_values: pheno_values,
        covariate_matrix: cov_matrix,
        covariate_names: expanded_cov_names,
        binary_levels: binary_levels.map(|l| l.to_vec()),
    }
}

/// Encode covariates into a numeric matrix.
/// Numeric covariates are kept as-is. Categorical covariates (non-numeric) are
/// dummy-encoded using N-1 indicator columns (the first level is the reference).
fn encode_covariates(
    raw: &[Vec<String>],
    cov_names: &[String],
    quiet: bool,
) -> (Vec<Vec<f64>>, Vec<String>) {
    if cov_names.is_empty() || raw.is_empty() {
        return (vec![Vec::new(); raw.len()], Vec::new());
    }

    let n_cov = cov_names.len();
    let n_samples = raw.len();

    // Determine which covariates are numeric vs categorical
    let is_numeric: Vec<bool> = (0..n_cov)
        .map(|j| raw.iter().all(|row| row[j].parse::<f64>().is_ok()))
        .collect();

    // For categorical covariates, collect sorted unique levels
    let cat_levels: Vec<Option<Vec<String>>> = (0..n_cov)
        .map(|j| {
            if is_numeric[j] {
                None
            } else {
                let mut levels: Vec<String> = raw.iter().map(|row| row[j].clone()).collect();
                levels.sort();
                levels.dedup();
                Some(levels)
            }
        })
        .collect();

    // Build expanded column names
    let mut expanded_names = Vec::new();
    for j in 0..n_cov {
        if is_numeric[j] {
            expanded_names.push(cov_names[j].clone());
        } else {
            let levels = cat_levels[j].as_ref().unwrap();
            if !quiet {
                eprintln!(
                    "Covariate '{}' is categorical with {} levels (reference: '{}')",
                    cov_names[j],
                    levels.len(),
                    levels[0]
                );
            }
            // Skip the first level (reference), create a column for each remaining level
            for level in &levels[1..] {
                expanded_names.push(format!("{}_{}", cov_names[j], level));
            }
        }
    }

    // Build numeric matrix
    let mut matrix = Vec::with_capacity(n_samples);
    for row in raw {
        let mut numeric_row = Vec::with_capacity(expanded_names.len());
        for j in 0..n_cov {
            if is_numeric[j] {
                numeric_row.push(row[j].parse::<f64>().unwrap());
            } else {
                let levels = cat_levels[j].as_ref().unwrap();
                // One-hot for levels[1..], reference level gets all zeros
                for level in &levels[1..] {
                    numeric_row.push(if row[j] == *level { 1.0 } else { 0.0 });
                }
            }
        }
        matrix.push(numeric_row);
    }

    (matrix, expanded_names)
}

// ---------------------------------------------------------------------------
// Variant result
// ---------------------------------------------------------------------------

struct VariantResult {
    variant_id: String,
    effect: f64, // OR for binary, Beta for continuous
    effect_l95: f64,
    effect_u95: f64,
    std_err: f64,
    pvalue: f64,
    n: usize,
    avg_size: f64,
    // Binary-specific
    group_stats: Option<Vec<(usize, f64, f64)>>, // (N, AvgSize, MaxSize) per group
    // Continuous-specific
    min_size: Option<f64>,
    max_size: Option<f64>,
    // Bonferroni
    pvalue_bonf: Option<f64>,
}

// ---------------------------------------------------------------------------
// Per-variant processing
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn process_variant_line(
    line: &str,
    sample_names: &[String],
    sample_mapping: &[Option<usize>],
    pheno: &PhenotypeData,
    str_mode: &str,
    outcometype: &str,
    missing_cutoff: f64,
    minimal_length: Option<f64>,
    input_type: &InputType,
    total_phenotyped: usize,
) -> Option<VariantResult> {
    let fields: Vec<&str> = line.trim().split('\t').collect();

    let (variant_id, values) = match input_type {
        InputType::Kmer => {
            if fields.len() < 1 + sample_names.len() {
                return None;
            }
            let kmer = fields[0].to_string();
            let vals: Vec<Option<f64>> = fields[1..1 + sample_names.len()]
                .iter()
                .map(|f| parse_field(f))
                .collect();
            (kmer, vals)
        }
        InputType::Str => {
            if fields.len() < 4 + 2 * sample_names.len() {
                return None;
            }
            let vid = format!("{}:{}-{}", fields[0], fields[1], fields[2]);
            let mut vals = Vec::with_capacity(sample_names.len());
            for i in 0..sample_names.len() {
                let h1 = parse_field(fields[4 + i * 2]);
                let h2 = parse_field(fields[4 + i * 2 + 1]);
                let combined = match (h1, h2) {
                    (Some(a), Some(b)) => Some(combine_alleles(a, b, str_mode)),
                    (Some(a), None) | (None, Some(a)) => Some(a),
                    (None, None) => None,
                };
                vals.push(combined);
            }
            (vid, vals)
        }
    };

    // Map variant values to phenotyped samples
    let mut variant_vals = Vec::new();
    let mut pheno_vals = Vec::new();
    let mut cov_rows: Vec<&Vec<f64>> = Vec::new();

    for (sample_idx, mapping) in sample_mapping.iter().enumerate() {
        if let Some(pheno_idx) = mapping
            && let Some(val) = values[sample_idx]
        {
            variant_vals.push(val);
            pheno_vals.push(pheno.phenotype_values[*pheno_idx]);
            if !pheno.covariate_matrix.is_empty() {
                cov_rows.push(&pheno.covariate_matrix[*pheno_idx]);
            }
        }
    }

    let n = variant_vals.len();

    // Check call rate
    if total_phenotyped == 0 || (n as f64 / total_phenotyped as f64) < missing_cutoff {
        return None;
    }

    // Check variance
    let first = variant_vals[0];
    if variant_vals
        .iter()
        .all(|&v| (v - first).abs() < f64::EPSILON)
    {
        return None;
    }

    // Check minimal length filter
    if let Some(min_len) = minimal_length {
        let max_val = variant_vals
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        if max_val <= min_len {
            return None;
        }
    }

    // Build design matrix: intercept + variant_value + covariates
    let n_cov = pheno.covariate_names.len();
    let p = 2 + n_cov; // intercept + variant + covariates
    let mut x = vec![0.0; n * p];
    for i in 0..n {
        x[i * p] = 1.0; // intercept
        x[i * p + 1] = variant_vals[i]; // variant value
        for j in 0..n_cov {
            x[i * p + 2 + j] = cov_rows[i][j];
        }
    }

    // Fit model
    let fit_result = if outcometype == "binary" {
        fit_logistic(&x, &pheno_vals, n, p)
    } else {
        fit_gaussian(&x, &pheno_vals, n, p)
    };

    let (beta, se, pval) = fit_result?;

    // Extract the variant coefficient (index 1)
    let variant_beta = beta[1];
    let variant_se = se[1];
    let variant_pval = pval[1];

    // Calculate effect and CI
    let avg_size = variant_vals.iter().sum::<f64>() / n as f64;

    if outcometype == "binary" {
        let or = variant_beta.exp();
        let or_l95 = (variant_beta - 1.96 * variant_se).exp();
        let or_u95 = (variant_beta + 1.96 * variant_se).exp();

        // Group statistics
        let group_stats = pheno.binary_levels.as_ref().map(|_levels| {
            let mut stats = Vec::new();
            for level_val in [0.0, 1.0] {
                let group_variant_vals: Vec<f64> = pheno_vals
                    .iter()
                    .zip(variant_vals.iter())
                    .filter(|&(&p, _)| (p - level_val).abs() < f64::EPSILON)
                    .map(|(_, &v)| v)
                    .collect();
                let group_n = group_variant_vals.len();
                let (group_avg, group_max) = if group_n > 0 {
                    (
                        group_variant_vals.iter().sum::<f64>() / group_n as f64,
                        group_variant_vals
                            .iter()
                            .cloned()
                            .fold(f64::NEG_INFINITY, f64::max),
                    )
                } else {
                    (f64::NAN, f64::NAN)
                };
                stats.push((group_n, group_avg, group_max));
            }
            stats
        });

        Some(VariantResult {
            variant_id,
            effect: or,
            effect_l95: or_l95,
            effect_u95: or_u95,
            std_err: variant_se,
            pvalue: variant_pval,
            n,
            avg_size,
            group_stats,
            min_size: None,
            max_size: None,
            pvalue_bonf: None,
        })
    } else {
        let beta_l95 = variant_beta - 1.96 * variant_se;
        let beta_u95 = variant_beta + 1.96 * variant_se;
        let min_size = variant_vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_size = variant_vals
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);

        Some(VariantResult {
            variant_id,
            effect: variant_beta,
            effect_l95: beta_l95,
            effect_u95: beta_u95,
            std_err: variant_se,
            pvalue: variant_pval,
            n,
            avg_size,
            group_stats: None,
            min_size: Some(min_size),
            max_size: Some(max_size),
            pvalue_bonf: None,
        })
    }
}

fn parse_field(s: &str) -> Option<f64> {
    let s = s.trim();
    if s.is_empty() || s.eq_ignore_ascii_case("nan") || s.eq_ignore_ascii_case("na") {
        None
    } else {
        s.parse().ok()
    }
}

fn combine_alleles(h1: f64, h2: f64, mode: &str) -> f64 {
    match mode {
        "MEAN" => (h1 + h2) / 2.0,
        "MAX" => h1.max(h2),
        "MIN" => h1.min(h2),
        _ => h1.max(h2),
    }
}

// ---------------------------------------------------------------------------
// Linear algebra helpers (using nalgebra)
// ---------------------------------------------------------------------------

use nalgebra::{DMatrix, DVector};

/// Solve ordinary least squares: β = (X'X)^{-1} X'y
/// Returns (coefficients, standard errors, p-values) or None on failure.
/// Uses a flat row-major slice to avoid per-variant DMatrix allocation for small systems.
fn fit_gaussian(
    x_flat: &[f64],
    y: &[f64],
    n: usize,
    p: usize,
) -> Option<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    if n <= p {
        return None;
    }

    // For small p (typical: 2-10), use direct fixed-size math to avoid DMatrix overhead
    // X'X is p×p, X'y is p×1
    let mut xtx = vec![0.0; p * p];
    let mut xty = vec![0.0; p];

    for i in 0..n {
        let row = &x_flat[i * p..(i + 1) * p];
        for j in 0..p {
            xty[j] += row[j] * y[i];
            for k in j..p {
                let val = row[j] * row[k];
                xtx[j * p + k] += val;
                if k != j {
                    xtx[k * p + j] += val;
                }
            }
        }
    }

    // Use nalgebra only for the p×p inversion (p is small)
    let xtx_mat = DMatrix::from_row_slice(p, p, &xtx);
    let xtx_inv = xtx_mat.try_inverse()?;
    let xty_vec = DVector::from_column_slice(&xty);
    let beta = &xtx_inv * &xty_vec;

    // Compute RSS without forming full n-vectors
    let mut rss = 0.0;
    for i in 0..n {
        let row = &x_flat[i * p..(i + 1) * p];
        let mut y_hat_i = 0.0;
        for j in 0..p {
            y_hat_i += row[j] * beta[j];
        }
        let r = y[i] - y_hat_i;
        rss += r * r;
    }

    let df = (n - p) as f64;
    let sigma2 = rss / df;
    let cov = &xtx_inv * sigma2;

    let mut betas = Vec::with_capacity(p);
    let mut ses = Vec::with_capacity(p);
    let mut pvals = Vec::with_capacity(p);

    for j in 0..p {
        let b = beta[j];
        let se = cov[(j, j)].sqrt();
        let t_stat = if se > 0.0 { b / se } else { 0.0 };
        let pval = if se > 0.0 {
            2.0 * t_cdf_upper(t_stat.abs(), df)
        } else {
            1.0
        };
        betas.push(b);
        ses.push(se);
        pvals.push(pval);
    }

    Some((betas, ses, pvals))
}

/// Fit logistic regression via IRLS (iteratively reweighted least squares).
/// Returns (coefficients, standard errors, p-values) or None on failure.
/// Avoids constructing n×n diagonal weight matrices by using element-wise operations.
fn fit_logistic(
    x_flat: &[f64],
    y: &[f64],
    n: usize,
    p: usize,
) -> Option<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    if n <= p {
        return None;
    }

    // Initialize beta to zeros
    let mut beta = vec![0.0; p];

    let max_iter = 25;
    let tol = 1e-8;

    // Scratch buffers reused across iterations
    let mut eta = vec![0.0; n];
    let mut mu = vec![0.0; n];
    let mut w = vec![0.0; n];
    let mut z = vec![0.0; n];
    let mut xtwx = vec![0.0; p * p];
    let mut xtwz = vec![0.0; p];

    for _iter in 0..max_iter {
        // η = X β, μ = sigmoid(η), w = μ(1-μ)
        for i in 0..n {
            let row = &x_flat[i * p..(i + 1) * p];
            let mut e = 0.0;
            for j in 0..p {
                e += row[j] * beta[j];
            }
            eta[i] = e;
            mu[i] = sigmoid(e);
            w[i] = mu[i] * (1.0 - mu[i]);
        }

        // Clamp near-zero weights to avoid division by zero (quasi-separation)
        // Instead of dropping the locus, this produces large coefficients (Inf OR),
        // matching R's glm() behavior.
        for wi in w.iter_mut() {
            if *wi < 1e-10 {
                *wi = 1e-10;
            }
        }

        // Working response: z = η + (y - μ) / w
        for i in 0..n {
            z[i] = eta[i] + (y[i] - mu[i]) / w[i];
        }

        // X'WX and X'Wz — no n×n matrix needed
        xtwx.fill(0.0);
        xtwz.fill(0.0);
        for i in 0..n {
            let row = &x_flat[i * p..(i + 1) * p];
            let wi = w[i];
            for j in 0..p {
                let wj = row[j] * wi;
                xtwz[j] += wj * z[i];
                for k in j..p {
                    let val = wj * row[k];
                    xtwx[j * p + k] += val;
                    if k != j {
                        xtwx[k * p + j] += val;
                    }
                }
            }
        }

        let xtwx_mat = DMatrix::from_row_slice(p, p, &xtwx);
        let xtwx_inv = xtwx_mat.try_inverse()?;
        let xtwz_vec = DVector::from_column_slice(&xtwz);
        let beta_new = &xtwx_inv * &xtwz_vec;

        // Check convergence
        let mut diff_sq = 0.0;
        for j in 0..p {
            let d = beta_new[j] - beta[j];
            diff_sq += d * d;
        }
        for j in 0..p {
            beta[j] = beta_new[j];
        }

        if diff_sq.sqrt() < tol {
            break;
        }
    }

    // Final variance: recompute weights at converged beta, invert X'WX
    for i in 0..n {
        let row = &x_flat[i * p..(i + 1) * p];
        let mut e = 0.0;
        for j in 0..p {
            e += row[j] * beta[j];
        }
        let m = sigmoid(e);
        w[i] = m * (1.0 - m);
    }

    for wi in w.iter_mut() {
        if *wi < 1e-10 {
            *wi = 1e-10;
        }
    }

    xtwx.fill(0.0);
    for i in 0..n {
        let row = &x_flat[i * p..(i + 1) * p];
        let wi = w[i];
        for j in 0..p {
            let wj = row[j] * wi;
            for k in j..p {
                let val = wj * row[k];
                xtwx[j * p + k] += val;
                if k != j {
                    xtwx[k * p + j] += val;
                }
            }
        }
    }

    let xtwx_mat = DMatrix::from_row_slice(p, p, &xtwx);
    let cov = xtwx_mat.try_inverse()?;

    let mut betas = Vec::with_capacity(p);
    let mut ses = Vec::with_capacity(p);
    let mut pvals = Vec::with_capacity(p);

    for j in 0..p {
        let b = beta[j];
        let se = cov[(j, j)].sqrt();
        let z_stat = if se > 0.0 { b / se } else { 0.0 };
        let pval = if se > 0.0 {
            2.0 * normal_cdf_upper(z_stat.abs())
        } else {
            1.0
        };
        betas.push(b);
        ses.push(se);
        pvals.push(pval);
    }

    Some((betas, ses, pvals))
}

fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

// ---------------------------------------------------------------------------
// Statistical distribution functions
// ---------------------------------------------------------------------------

use statrs::distribution::{ContinuousCDF, Normal, StudentsT};

/// Upper-tail probability of the standard normal: P(Z > z)
fn normal_cdf_upper(z: f64) -> f64 {
    let n = Normal::new(0.0, 1.0).unwrap();
    1.0 - n.cdf(z)
}

/// Upper-tail probability of Student's t distribution: P(T > t)
fn t_cdf_upper(t: f64, df: f64) -> f64 {
    let dist = StudentsT::new(0.0, 1.0, df).unwrap();
    1.0 - dist.cdf(t)
}

// ---------------------------------------------------------------------------
// Output writing
// ---------------------------------------------------------------------------

fn write_results(
    out: &std::path::Path,
    results: &[VariantResult],
    outcometype: &str,
    binary_levels: Option<&[String]>,
    _quiet: bool,
) {
    let file = std::fs::File::create(out).unwrap_or_else(|e| {
        eprintln!("ERROR: Could not create output file {}: {}", out.display(), e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(file);

    // Write header
    if outcometype == "binary" {
        let mut header = "VariantID\tOR\tOR_L95\tOR_U95\tOR_stdErr\tPvalue\tN\tAvgSize".to_string();
        if let Some(levels) = binary_levels {
            for level in levels {
                header.push_str(&format!("\t{}_N\t{}_AvgSize\t{}_MaxSize", level, level, level));
            }
        }
        header.push_str("\tPvalue_bonf");
        writeln!(writer, "{}", header).unwrap();
    } else {
        writeln!(
            writer,
            "VariantID\tBeta\tBeta_L95\tBeta_U95\tBeta_stdErr\tPvalue\tN\tAvgSize\tMinSize\tMaxSize\tPvalue_bonf"
        )
        .unwrap();
    }

    for r in results {
        if outcometype == "binary" {
            let mut line = format!(
                "{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{:.3}",
                r.variant_id,
                r.effect,
                r.effect_l95,
                r.effect_u95,
                r.std_err,
                r.pvalue,
                r.n,
                r.avg_size,
            );
            if let Some(ref gs) = r.group_stats {
                for (gn, gavg, gmax) in gs {
                    line.push_str(&format!("\t{}\t{:.3}\t{:.3}", gn, gavg, gmax));
                }
            }
            line.push_str(&format!("\t{}", r.pvalue_bonf.unwrap_or(f64::NAN)));
            writeln!(writer, "{}", line).unwrap();
        } else {
            writeln!(
                writer,
                "{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}",
                r.variant_id,
                r.effect,
                r.effect_l95,
                r.effect_u95,
                r.std_err,
                r.pvalue,
                r.n,
                r.avg_size,
                r.min_size.unwrap_or(f64::NAN),
                r.max_size.unwrap_or(f64::NAN),
                r.pvalue_bonf.unwrap_or(f64::NAN),
            )
            .unwrap();
        }
    }

    writer.flush().unwrap();
}

// ---------------------------------------------------------------------------
// Manhattan plot
// ---------------------------------------------------------------------------

fn parse_str_coords(variant_id: &str) -> Option<(String, f64)> {
    let colon = variant_id.find(':')?;
    let chrom = variant_id[..colon].to_string();
    let rest = &variant_id[colon + 1..];
    let dash = rest.find('-')?;
    let begin: f64 = rest[..dash].parse().ok()?;
    Some((chrom, begin))
}

fn write_manhattan_plot(results: &[VariantResult], plot_path: &std::path::Path, quiet: bool) {
    use kuva::prelude::*;

    let data: Vec<(String, f64, f64)> = results
        .iter()
        .filter_map(|r| {
            let (chrom, pos) = parse_str_coords(&r.variant_id)?;
            Some((chrom, pos, r.pvalue))
        })
        .collect();

    if data.is_empty() {
        eprintln!(
            "WARNING: Cannot generate Manhattan plot: no genomic coordinates found in variant IDs \
             (Manhattan plots require STR data, not kmer data)."
        );
        return;
    }

    let n_tests = results.len();
    let bonf_threshold = -(0.05_f64 / n_tests as f64).log10();

    let mp = ManhattanPlot::new()
        .with_data_bp(data, GenomeBuild::Hg38)
        .with_genome_wide(bonf_threshold)
        .with_legend("Association thresholds");

    let plots = vec![Plot::Manhattan(mp)];
    let layout = Layout::auto_from_plots(&plots)
        .with_title("inquiSTR STR Association")
        .with_x_label("Chromosome")
        .with_y_label("−log₁₀(p-value)");

    let svg = render_to_svg(plots, layout);
    std::fs::write(plot_path, svg).unwrap_or_else(|e| {
        eprintln!("ERROR: Could not write Manhattan plot to {}: {}", plot_path.display(), e);
    });

    if !quiet {
        eprintln!("Manhattan plot written to: {}", plot_path.display());
    }
}
