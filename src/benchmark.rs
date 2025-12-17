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
    tier: String,
}

/// Parse BED file with 9 columns (last 2 are haplotype lengths)
fn parse_bed_file(
    file_path: &Path,
    max_locus: Option<u32>,
) -> Result<HashMap<String, TruthRecord>, Box<dyn std::error::Error>> {
    let reader = crate::utils::reader(file_path.to_string_lossy().as_ref());
    let mut records = HashMap::new();
    let mut line_count = 0;
    let mut filtered_count = 0;

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
        let end: u32 = fields[2].parse()?;
        let tier = fields[3].to_string();

        // Filter by max_locus size if specified
        if let Some(max_size) = max_locus {
            let locus_size = end - begin;
            if locus_size > max_size {
                filtered_count += 1;
                continue;
            }
        }

        // Parse the last 2 columns as haplotype lengths and flip their signs
        // Handle zero specially to avoid -0.0
        let h1_raw: f64 = fields[7].parse()?;
        let h1 = if h1_raw == 0.0 { 0.0 } else { -h1_raw };

        let h2_raw: f64 = fields[8].parse()?;
        let h2 = if h2_raw == 0.0 { 0.0 } else { -h2_raw };

        let record = TruthRecord {
            chromosome: chromosome.clone(),
            pos: begin, // Using begin position for matching
            h1,
            h2,
            tier,
        };

        // Use chromosome:begin as key for matching
        let key = format!("{}:{}", chromosome, begin);
        records.insert(key, record);
    }

    if filtered_count > 0 {
        println!(
            "BED file: Filtered out {} intervals larger than {} bp (max-locus limit)",
            filtered_count,
            max_locus.unwrap()
        );
    }
    println!("BED file processed: {} lines, {} records retained", line_count, records.len());

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
        if fields.len() != 6 {
            eprintln!("Warning: Skipping malformed line {}: {}", line_num + 1, line);
            continue;
        }

        let chromosome = fields[0].to_string();
        let begin: u32 = fields[1].parse()?;
        let end: u32 = fields[2].parse()?;
        // Skip info field at index 3

        // Parse H1 and H2, handling NaN values
        let h1 = fields[4].parse::<f64>().unwrap_or(f64::NAN);
        let h2 = fields[5].parse::<f64>().unwrap_or(f64::NAN);

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
    max_locus: Option<u32>,
) -> Result<HashMap<String, TruthRecord>, Box<dyn std::error::Error>> {
    let reader = crate::utils::reader(&file_path.to_string_lossy());
    let mut records = HashMap::new();
    let mut total_variants = 0;
    let mut snp_variants = 0;
    let mut filtered_count = 0;

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

        // Filter by max_locus size if specified
        // For VCF, we use the REF allele length as the locus size
        if let Some(max_size) = max_locus
            && ref_seq.len() as u32 > max_size
        {
            filtered_count += 1;
            continue;
        }

        // For VCF, we need to determine H1 and H2 from the available ALT alleles
        // Since we don't have phasing info, we'll use the available ALT lengths
        let (h1, h2) = if alt_lengths.len() >= 2 {
            (alt_lengths[0], alt_lengths[1])
        } else if alt_lengths.len() == 1 {
            (alt_lengths[0], alt_lengths[0]) // Duplicate if only one ALT
        } else {
            (0.0, 0.0) // No ALT alleles
        };

        let record =
            TruthRecord { chromosome: chromosome.clone(), pos, h1, h2, tier: "Tier1".to_string() };

        // Use chromosome:pos as key for matching
        let key = format!("{}:{}", chromosome, pos);
        records.insert(key, record);
    }

    if filtered_count > 0 {
        println!(
            "VCF file: Filtered out {} variants with REF length > {} bp (max-locus limit)",
            filtered_count,
            max_locus.unwrap()
        );
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
#[allow(clippy::too_many_arguments)]
pub fn benchmark(
    inquistr_file: PathBuf,
    vcf_file: Option<PathBuf>,
    bed_file: Option<PathBuf>,
    mode: String,
    plot_file: PathBuf,
    max_plot_length: f64,
    tier1_only: bool,
    diff_out: Option<PathBuf>,
    max_locus: Option<u32>,
    nonzero: bool,
    tolerance: u32,
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
        match parse_vcf_file(&vcf_path, max_locus) {
            Ok(records) => records,
            Err(e) => {
                eprintln!("Error parsing VCF file: {}", e);
                std::process::exit(1);
            }
        }
    } else if let Some(bed_path) = bed_file {
        println!("Loading BED file: {}", bed_path.display());
        match parse_bed_file(&bed_path, max_locus) {
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

    // Filter truth records by tier if requested, but keep the full set for counting
    let (truth_records_filtered, all_truth_records) = if tier1_only {
        let original_count = truth_records.len();
        let all_truth = truth_records.clone(); // Keep full set for later
        let filtered: HashMap<String, TruthRecord> = truth_records
            .into_iter()
            .filter(|(_, record)| record.tier == "Tier1")
            .collect();
        let tier2_count = original_count - filtered.len();
        println!(
            "Filtered to Tier1 only: {} Tier1 variants, {} Tier2 variants ignored",
            filtered.len(),
            tier2_count
        );
        (filtered, all_truth)
    } else {
        let cloned = truth_records.clone();
        (truth_records, cloned)
    };

    // Find matching loci and collect data for correlation
    let mut inquistr_values = Vec::new();
    let mut truth_values = Vec::new();
    let mut loci_info = Vec::new(); // Store (chromosome, begin, end) for each matched locus
    let mut nan_count = 0;
    let mut matched_loci = Vec::new();
    let mut matched_tier2_only = 0; // Track inquiSTR loci that only match Tier2 variants
    let mut zero_zero_filtered = 0; // Track zero-zero pairs filtered out by --nonzero flag

    for (key, inquistr_record) in &inquistr_records {
        if let Some(truth_record) = truth_records_filtered.get(key) {
            // Calculate inquiSTR value (MAX or MIN of H1 and H2)
            let inquistr_alleles = vec![inquistr_record.h1, inquistr_record.h2];

            if let Some(inquistr_value) = select_allele(&inquistr_alleles, &mode) {
                // Calculate truth value (MAX or MIN of H1 and H2)
                let truth_alleles = vec![truth_record.h1, truth_record.h2];

                if let Some(truth_value) = select_allele(&truth_alleles, &mode) {
                    // Skip zero-zero pairs if --nonzero flag is set
                    if nonzero && inquistr_value == 0.0 && truth_value == 0.0 {
                        zero_zero_filtered += 1;
                        continue;
                    }

                    inquistr_values.push(inquistr_value);
                    truth_values.push(truth_value);
                    loci_info.push((
                        inquistr_record.chromosome.clone(),
                        inquistr_record.begin,
                        inquistr_record.end,
                    ));
                    matched_loci.push(key.clone());
                }
            } else {
                // inquiSTR locus was targeted but has NaN values
                nan_count += 1;
            }
        } else if tier1_only && all_truth_records.contains_key(key) {
            // This inquiSTR locus matches a Tier2 variant (not in filtered set, but in full set)
            matched_tier2_only += 1;
        }
    }

    let matched_count = inquistr_values.len();
    let inquistr_only = inquistr_records.len() - matched_count - nan_count - matched_tier2_only;
    let truth_only = truth_records_filtered.len() - matched_count;

    println!("\n=== Loci Assessment Summary ===");
    println!("Total inquiSTR loci: {}", inquistr_records.len());
    println!("Total truth loci (after max-locus filtering): {}", truth_records_filtered.len());
    println!("\nBreakdown:");
    println!("  Loci found in inquiSTR only: {}", inquistr_only);
    if matched_tier2_only > 0 {
        println!("  Loci in inquiSTR matching Tier2 variants only: {}", matched_tier2_only);
    }
    println!("  Loci found in truth data only: {}", truth_only);
    println!("  Loci with NaN in inquiSTR (excluded): {}", nan_count);
    if nonzero && zero_zero_filtered > 0 {
        println!("  Zero-zero pairs filtered out (--nonzero): {}", zero_zero_filtered);
    }
    println!("\n  Loci successfully matched and assessed: {} ✓", matched_count);
    println!("===============================\n");

    if matched_count == 0 {
        eprintln!("Error: No matching loci found between inquiSTR and truth files");
        std::process::exit(1);
    }

    // Calculate correlation - behavior depends on --nonzero flag
    let (correlation_all, r_squared_all, zero_zero_count) = if nonzero {
        // When --nonzero is used, data is already filtered, so just calculate correlation
        let corr = pearson_correlation(&inquistr_values, &truth_values);
        (corr, corr * corr, zero_zero_filtered)
    } else {
        // Calculate both all and nonzero correlations for comparison
        let corr_all = pearson_correlation(&inquistr_values, &truth_values);
        let r2_all = corr_all * corr_all;

        // Count zero-zero pairs in the data
        let zero_count = inquistr_values
            .iter()
            .zip(truth_values.iter())
            .filter(|&(inq, truth)| *inq == 0.0 && *truth == 0.0)
            .count();

        (corr_all, r2_all, zero_count)
    };

    // For nonzero mode, also calculate what correlation would have been with zeros
    let (correlation_nonzero, r_squared_nonzero, nonzero_count) = if !nonzero {
        let mut inquistr_nonzero = Vec::new();
        let mut truth_nonzero = Vec::new();

        for (&inq, &truth) in inquistr_values.iter().zip(truth_values.iter()) {
            if !(inq == 0.0 && truth == 0.0) {
                inquistr_nonzero.push(inq);
                truth_nonzero.push(truth);
            }
        }

        let corr = if !inquistr_nonzero.is_empty() {
            pearson_correlation(&inquistr_nonzero, &truth_nonzero)
        } else {
            f64::NAN
        };
        (corr, corr * corr, inquistr_nonzero.len())
    } else {
        // In nonzero mode, the "nonzero" values are just the main values
        (correlation_all, r_squared_all, matched_count)
    };

    println!("\n=== Correlation Analysis ===");
    if nonzero {
        println!(
            "Nonzero loci only (n={}, {} zero-zero pairs excluded): R={:.4}, R²={:.4}",
            matched_count, zero_zero_count, correlation_all, r_squared_all
        );
    } else {
        println!(
            "All loci (n={}): R={:.4}, R²={:.4}",
            matched_count, correlation_all, r_squared_all
        );
        if zero_zero_count > 0 {
            println!("  - Including {} zero-zero pairs (unchanged alleles)", zero_zero_count);
            println!(
                "Excluding zeros (n={}): R={:.4}, R²={:.4}",
                nonzero_count, correlation_nonzero, r_squared_nonzero
            );
        }
    }
    println!("============================");

    // Calculate exact matches and matches within tolerance
    let mut exact_matches = 0;
    let mut within_tolerance = 0;

    for (&inq_val, &truth_val) in inquistr_values.iter().zip(truth_values.iter()) {
        let diff = (inq_val - truth_val).abs();
        if diff == 0.0 {
            exact_matches += 1;
            within_tolerance += 1; // Exact matches are also within tolerance
        } else if diff <= tolerance as f64 {
            within_tolerance += 1;
        }
    }

    let exact_percent = (exact_matches as f64 / matched_count as f64) * 100.0;
    let within_tolerance_percent = (within_tolerance as f64 / matched_count as f64) * 100.0;

    println!("\n=== Accuracy Analysis ===");
    println!("Exact matches: {} ({:.2}%)", exact_matches, exact_percent);
    println!(
        "Within {} bp tolerance: {} ({:.2}%)",
        tolerance, within_tolerance, within_tolerance_percent
    );
    println!("=========================");

    // Output top 100 loci with largest differences if requested
    if let Some(diff_out_path) = diff_out {
        // Reduce type complexity by using a local type alias for readability
        type LocusInfo<'a> = &'a (String, u32, u32);
        type DiffEntry<'a> = (f64, LocusInfo<'a>, f64, f64);

        let mut differences: Vec<DiffEntry> = inquistr_values
            .iter()
            .zip(truth_values.iter())
            .zip(loci_info.iter())
            .map(|((&inq_val, &truth_val), locus)| {
                let diff = (inq_val - truth_val).abs();
                (diff, locus, inq_val, truth_val)
            })
            .collect();

        // Sort by absolute difference (largest first)
        differences.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

        // Take top 100 and write to file
        let top_100 = differences.iter().take(100);

        match std::fs::File::create(&diff_out_path) {
            Ok(file) => {
                use std::io::Write;
                let mut writer = std::io::BufWriter::new(file);

                // Write header
                if writeln!(
                    writer,
                    "chromosome\tbegin\tend\tinquistr_genotype\ttruth_genotype\tabsolute_difference"
                )
                .is_err()
                {
                    eprintln!("Error writing header to diff output file");
                } else {
                    // Write data
                    let mut count = 0;
                    for (diff, locus, inq_val, truth_val) in top_100 {
                        if writeln!(
                            writer,
                            "{}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}",
                            locus.0, locus.1, locus.2, inq_val, truth_val, diff
                        )
                        .is_err()
                        {
                            eprintln!("Error writing data to diff output file");
                            break;
                        }
                        count += 1;
                    }
                    println!(
                        "Top {} loci with largest differences written to: {}",
                        count,
                        diff_out_path.display()
                    );
                }
            }
            Err(e) => {
                eprintln!("Error creating diff output file: {}", e);
            }
        }
    }

    // Filter points for plotting based on max_plot_length
    let mut plot_inquistr = Vec::new();
    let mut plot_truth = Vec::new();
    let mut plot_hover_text = Vec::new();
    let mut hidden_count = 0;

    for ((inq_val, truth_val), locus) in inquistr_values
        .iter()
        .zip(truth_values.iter())
        .zip(loci_info.iter())
    {
        if inq_val.abs() <= max_plot_length && truth_val.abs() <= max_plot_length {
            plot_inquistr.push(*inq_val);
            plot_truth.push(*truth_val);
            plot_hover_text.push(format!(
                "{}:{}-{}<br>inquiSTR: {:.1}<br>Truth: {:.1}",
                locus.0, locus.1, locus.2, inq_val, truth_val
            ));
        } else {
            hidden_count += 1;
        }
    }

    if hidden_count > 0 {
        println!(
            "{} points that are larger than {} are not shown on the plot",
            hidden_count, max_plot_length
        );
    }

    // Determine axis range for square plot
    let axis_limit = max_plot_length;

    // Create scatter plot with custom hover text
    let trace = Scatter::new(plot_inquistr.clone(), plot_truth.clone())
        .mode(Mode::Markers)
        .name("Data points")
        .text_array(plot_hover_text)
        .hover_template("%{text}<extra></extra>");

    let title_text = if nonzero {
        format!(
            "inquiSTR vs Truth Genotypes (Mode: {}, Nonzero only, R² = {:.4})",
            mode, r_squared_all
        )
    } else {
        format!("inquiSTR vs Truth Genotypes (Mode: {}, R² = {:.4})", mode, r_squared_all)
    };

    let layout = Layout::new()
        .title(Title::with_text(title_text))
        .x_axis(
            Axis::new()
                .title(Title::with_text("inquiSTR genotypes"))
                .range(vec![-axis_limit, axis_limit])
                .constrain(plotly::layout::AxisConstrain::Domain)
                .scale_anchor("y"),
        )
        .y_axis(
            Axis::new()
                .title(Title::with_text("Truth genotypes"))
                .range(vec![-axis_limit, axis_limit])
                .scale_anchor("x"),
        )
        .width(800)
        .height(800);

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

    // Output summary in parseable format for scripting
    println!("\n=== BENCHMARK SUMMARY ===");
    if nonzero {
        println!("NONZERO_MODE: true");
        println!("LOCI_ASSESSED: {}", matched_count);
        println!("ZERO_ZERO_PAIRS_EXCLUDED: {}", zero_zero_count);
        println!("PEARSON_R: {:.6}", correlation_all);
        println!("R_SQUARED: {:.6}", r_squared_all);
    } else {
        println!("LOCI_ASSESSED: {}", matched_count);
        println!("ZERO_ZERO_PAIRS: {}", zero_zero_count);
        println!("NONZERO_LOCI: {}", nonzero_count);
        println!("PEARSON_R_ALL: {:.6}", correlation_all);
        println!("R_SQUARED_ALL: {:.6}", r_squared_all);
        println!("PEARSON_R_NONZERO: {:.6}", correlation_nonzero);
        println!("R_SQUARED_NONZERO: {:.6}", r_squared_nonzero);
    }
    println!("EXACT_MATCHES: {}", exact_matches);
    println!("EXACT_MATCHES_PERCENT: {:.2}", exact_percent);
    println!("WITHIN_TOLERANCE: {}", within_tolerance);
    println!("WITHIN_TOLERANCE_PERCENT: {:.2}", within_tolerance_percent);
    println!("TOLERANCE_BP: {}", tolerance);
    println!("=========================");
}
