use plotly::common::{Mode, Title};
use plotly::layout::{Axis, Layout};
use plotly::{Plot, Scatter};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
struct InquiSTRRecord {
    #[allow(dead_code)]
    chromosome: String,
    #[allow(dead_code)]
    begin: u32,
    #[allow(dead_code)]
    end: u32,
    h1: f64,
    h2: f64,
}

#[derive(Debug, Clone)]
struct TruthRecord {
    #[allow(dead_code)]
    chromosome: String,
    #[allow(dead_code)]
    pos: u32,
    h1: f64,
    h2: f64,
}

/// Parse BED file with 9 columns (last 2 are haplotype lengths)
fn parse_bed_file(
    file_path: &Path,
) -> Result<HashMap<String, TruthRecord>, Box<dyn std::error::Error>> {
    let reader = crate::utils::reader(file_path.to_string_lossy().as_ref());
    let mut records = HashMap::new();
    let mut line_count = 0;

    for line in reader.lines() {
        let line = line?;
        line_count += 1;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 9 {
            eprintln!(
                "Warning: Skipping malformed BED line {} (expected 9 fields, got {}): {}",
                line_count,
                fields.len(),
                line
            );
            continue;
        }

        let chromosome = fields[0].to_string();
        let begin: u32 = fields[1].parse()?;

        // Parse the last 2 columns as haplotype lengths and flip their signs
        let h1: f64 = -fields[7].parse::<f64>()?;
        let h2: f64 = -fields[8].parse::<f64>()?;

        let record = TruthRecord {
            chromosome: chromosome.clone(),
            pos: begin, // Using begin position for matching
            h1,
            h2,
        };

        // Use chromosome:begin as key for matching
        let key = format!("{}:{}", chromosome, begin);
        records.insert(key, record);
    }

    println!("BED file processed: {} lines, {} records", line_count, records.len());

    Ok(records)
}

/// Parse inquiSTR output file
fn parse_inquistr_file(
    file_path: &Path,
) -> Result<HashMap<String, InquiSTRRecord>, Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut records = HashMap::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;

        // Skip header line
        if line_num == 0 && line.starts_with("chromosome") {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 5 {
            eprintln!("Warning: Skipping malformed line {}: {}", line_num + 1, line);
            continue;
        }

        let chromosome = fields[0].to_string();
        let begin: u32 = fields[1].parse()?;
        let end: u32 = fields[2].parse()?;

        // Parse H1 and H2, handling NaN values
        let h1 = match fields[3].parse::<f64>() {
            Ok(val) => val,
            Err(_) => f64::NAN,
        };
        let h2 = match fields[4].parse::<f64>() {
            Ok(val) => val,
            Err(_) => f64::NAN,
        };

        let record = InquiSTRRecord { chromosome: chromosome.clone(), begin, end, h1, h2 };

        // Use chromosome:begin as key for matching
        let key = format!("{}:{}", chromosome, begin);
        records.insert(key, record);
    }

    Ok(records)
}

/// Parse VCF file (supports compressed files)
fn parse_vcf_file(
    file_path: &Path,
) -> Result<HashMap<String, TruthRecord>, Box<dyn std::error::Error>> {
    let reader = crate::utils::reader(&file_path.to_string_lossy());
    let mut records = HashMap::new();
    let mut total_variants = 0;
    let mut snp_variants = 0;

    for line in reader.lines() {
        let line = line?;

        // Skip header lines
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let chromosome = fields[0].to_string();
        let pos: u32 = fields[1].parse()?;
        let ref_seq = fields[3].to_string();
        let alt_field = fields[4];

        total_variants += 1;

        // Parse ALT field (can contain multiple alleles separated by commas)
        let alt_seqs: Vec<String> = alt_field.split(',').map(|s| s.to_string()).collect();

        // Filter out SNPs: variants where REF and all ALT alleles are single nucleotides
        let is_snp = ref_seq.len() == 1 && alt_seqs.iter().all(|alt| alt.len() == 1);

        if is_snp {
            snp_variants += 1;
            continue; // Skip SNPs as they're not STR variants
        }

        // Calculate allele lengths as ALT - REF for each ALT
        let ref_len = ref_seq.len() as f64;
        let alt_lengths: Vec<f64> = alt_seqs
            .iter()
            .map(|alt| alt.len() as f64 - ref_len)
            .collect();

        // For VCF, we need to determine H1 and H2 from the available ALT alleles
        // Since we don't have phasing info, we'll use the available ALT lengths
        let (h1, h2) = if alt_lengths.len() >= 2 {
            (alt_lengths[0], alt_lengths[1])
        } else if alt_lengths.len() == 1 {
            (alt_lengths[0], alt_lengths[0]) // Duplicate if only one ALT
        } else {
            (0.0, 0.0) // No ALT alleles
        };

        let record = TruthRecord { chromosome: chromosome.clone(), pos, h1, h2 };

        // Use chromosome:pos as key for matching
        let key = format!("{}:{}", chromosome, pos);
        records.insert(key, record);
    }

    println!(
        "VCF variants processed: {} total, {} SNPs filtered out, {} STR variants retained",
        total_variants,
        snp_variants,
        records.len()
    );

    Ok(records)
}

/// Select allele based on mode (MAX or MIN)
fn select_allele(values: &[f64], mode: &str) -> Option<f64> {
    let finite_values: Vec<f64> = values.iter().filter(|&&x| x.is_finite()).copied().collect();

    if finite_values.is_empty() {
        return None;
    }

    match mode.to_uppercase().as_str() {
        "MAX" => finite_values
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .copied(),
        "MIN" => finite_values
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .copied(),
        _ => {
            eprintln!("ERROR: Invalid mode '{}'. Must be 'MAX' or 'MIN'", mode);
            std::process::exit(1);
        }
    }
}

/// Calculate Pearson correlation coefficient
fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() || x.is_empty() {
        return f64::NAN;
    }

    let n = x.len() as f64;
    let sum_x: f64 = x.iter().sum();
    let sum_y: f64 = y.iter().sum();
    let sum_x2: f64 = x.iter().map(|&xi| xi * xi).sum();
    let sum_y2: f64 = y.iter().map(|&yi| yi * yi).sum();
    let sum_xy: f64 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi * yi).sum();

    let numerator = n * sum_xy - sum_x * sum_y;
    let denominator = ((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y)).sqrt();

    if denominator == 0.0 {
        f64::NAN
    } else {
        numerator / denominator
    }
}

/// Main benchmark function
pub fn benchmark(
    inquistr_file: PathBuf,
    vcf_file: Option<PathBuf>,
    bed_file: Option<PathBuf>,
    mode: String,
    plot_file: PathBuf,
) {
    // Validate that exactly one truth file is provided
    match (&vcf_file, &bed_file) {
        (Some(_), Some(_)) => {
            eprintln!("Error: Cannot specify both --vcf and --bed options. Choose one.");
            std::process::exit(1);
        }
        (None, None) => {
            eprintln!("Error: Must specify either --vcf or --bed option for truth data.");
            std::process::exit(1);
        }
        _ => {} // Exactly one is specified, which is correct
    }

    println!("Loading inquiSTR file: {}", inquistr_file.display());
    let inquistr_records = match parse_inquistr_file(&inquistr_file) {
        Ok(records) => records,
        Err(e) => {
            eprintln!("Error parsing inquiSTR file: {}", e);
            std::process::exit(1);
        }
    };

    let truth_records = if let Some(vcf_path) = vcf_file {
        println!("Loading VCF file: {}", vcf_path.display());
        match parse_vcf_file(&vcf_path) {
            Ok(records) => records,
            Err(e) => {
                eprintln!("Error parsing VCF file: {}", e);
                std::process::exit(1);
            }
        }
    } else if let Some(bed_path) = bed_file {
        println!("Loading BED file: {}", bed_path.display());
        match parse_bed_file(&bed_path) {
            Ok(records) => records,
            Err(e) => {
                eprintln!("Error parsing BED file: {}", e);
                std::process::exit(1);
            }
        }
    } else {
        unreachable!("Should have been caught by validation above");
    };

    println!("inquiSTR records loaded: {}", inquistr_records.len());
    println!("Truth records loaded: {}", truth_records.len());

    // Find matching loci and collect data for correlation
    let mut inquistr_values = Vec::new();
    let mut truth_values = Vec::new();
    let mut nan_count = 0;
    let mut matched_loci = Vec::new();

    for (key, inquistr_record) in &inquistr_records {
        if let Some(truth_record) = truth_records.get(key) {
            // Calculate inquiSTR value (MAX or MIN of H1 and H2)
            let inquistr_alleles = vec![inquistr_record.h1, inquistr_record.h2];

            if let Some(inquistr_value) = select_allele(&inquistr_alleles, &mode) {
                // Calculate truth value (MAX or MIN of H1 and H2)
                let truth_alleles = vec![truth_record.h1, truth_record.h2];

                if let Some(truth_value) = select_allele(&truth_alleles, &mode) {
                    inquistr_values.push(inquistr_value);
                    truth_values.push(truth_value);
                    matched_loci.push(key.clone());
                }
            } else {
                // inquiSTR locus was targeted but has NaN values
                nan_count += 1;
            }
        }
    }

    let matched_count = inquistr_values.len();
    let inquistr_only = inquistr_records.len() - matched_count - nan_count;
    let truth_only = truth_records.len() - matched_count;

    println!("Loci found in inquiSTR only: {}", inquistr_only);
    println!("Loci found in truth data only: {}", truth_only);
    println!("Loci found in both: {}", matched_count);
    println!("Loci with NaN in inquiSTR (but found in truth data): {}", nan_count);

    if matched_count == 0 {
        eprintln!("Error: No matching loci found between inquiSTR and truth files");
        std::process::exit(1);
    }

    // Calculate correlation
    let correlation = pearson_correlation(&inquistr_values, &truth_values);
    let r_squared = correlation * correlation;

    println!("Pearson correlation coefficient: {:.4}", correlation);
    println!("R² value: {:.4}", r_squared);

    // Create scatter plot
    let trace = Scatter::new(inquistr_values.clone(), truth_values.clone())
        .mode(Mode::Markers)
        .name("Data points");

    let layout = Layout::new()
        .title(Title::with_text(format!(
            "inquiSTR vs Truth Genotypes (Mode: {}, R² = {:.4})",
            mode, r_squared
        )))
        .x_axis(Axis::new().title(Title::with_text("inquiSTR genotypes")))
        .y_axis(Axis::new().title(Title::with_text("Truth genotypes")));

    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.set_layout(layout);

    // Save plot
    let html_content = plot.to_html();
    match std::fs::write(&plot_file, html_content) {
        Ok(_) => println!("Plot saved to: {}", plot_file.display()),
        Err(e) => {
            eprintln!("Error saving plot: {}", e);
            std::process::exit(1);
        }
    }
}
