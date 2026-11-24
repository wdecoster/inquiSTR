use kuva::plot::scatter::TrendLine;
use kuva::prelude::*;
use std::collections::HashMap;
use std::io::BufRead;
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
    chromosome: String,
    /// Interval start, in the same 0-based coordinate space as the test file.
    pos: u32,
    /// Interval end, used only by `--match-overlap`.
    end: u32,
    h1: f64,
    h2: f64,
    tier: String,
}

/// Intervals of one chromosome, sorted by start, supporting overlap queries.
///
/// `max_end_prefix[i]` is the largest end among intervals `0..=i`, which bounds how far back
/// a backwards scan has to go and keeps the query correct even if intervals overlap each other.
#[derive(Default)]
struct ChromIntervals {
    begins: Vec<u32>,
    ends: Vec<u32>,
    keys: Vec<String>,
    max_end_prefix: Vec<u32>,
}

/// Interval index over truth loci, for matching by overlap instead of by exact start.
#[derive(Default)]
struct IntervalIndex {
    by_chromosome: HashMap<String, ChromIntervals>,
}

impl IntervalIndex {
    /// Index every truth record for which `keep` returns true.
    fn build<F>(records: &HashMap<String, TruthRecord>, keep: F) -> Self
    where
        F: Fn(&TruthRecord) -> bool,
    {
        let mut by_chromosome: HashMap<String, ChromIntervals> = HashMap::new();
        for (key, record) in records {
            if !keep(record) {
                continue;
            }
            let entry = by_chromosome.entry(record.chromosome.clone()).or_default();
            entry.begins.push(record.pos);
            entry.ends.push(record.end);
            entry.keys.push(key.clone());
        }

        for intervals in by_chromosome.values_mut() {
            let mut order: Vec<usize> = (0..intervals.begins.len()).collect();
            order.sort_by_key(|&i| (intervals.begins[i], intervals.ends[i]));
            intervals.begins = order.iter().map(|&i| intervals.begins[i]).collect();
            intervals.ends = order.iter().map(|&i| intervals.ends[i]).collect();
            intervals.keys = order.iter().map(|&i| intervals.keys[i].clone()).collect();

            let mut running = 0;
            intervals.max_end_prefix = intervals
                .ends
                .iter()
                .map(|&end| {
                    running = running.max(end);
                    running
                })
                .collect();
        }

        IntervalIndex { by_chromosome }
    }

    /// The indexed interval sharing the most bases with `[begin, end)`, and the reciprocal
    /// overlap of that pair: shared bases divided by the length of the *longer* of the two
    /// intervals, so 1.0 means identical and a small value means one interval dwarfs the other.
    ///
    /// Ties are broken on the interval start, so the result does not depend on hash order.
    fn best_overlap(&self, chromosome: &str, begin: u32, end: u32) -> Option<(&str, f64)> {
        let intervals = self.by_chromosome.get(chromosome)?;
        if end <= begin {
            return None;
        }

        // Intervals starting at or after `end` cannot overlap.
        let upper = intervals.begins.partition_point(|&b| b < end);
        let mut best: Option<(u32, usize)> = None;
        for i in (0..upper).rev() {
            // Nothing at or before `i` can reach `begin` any more.
            if intervals.max_end_prefix[i] <= begin {
                break;
            }
            if intervals.ends[i] <= begin {
                continue;
            }
            let shared = intervals.ends[i].min(end) - intervals.begins[i].max(begin);
            if best.is_none_or(|(best_shared, _)| shared > best_shared) {
                best = Some((shared, i));
            }
        }

        let (shared, i) = best?;
        let longer = (end - begin).max(intervals.ends[i] - intervals.begins[i]);
        Some((intervals.keys[i].as_str(), f64::from(shared) / f64::from(longer)))
    }
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
            return Err(format!(
                "Error: Malformed BED line {} (expected 9 fields, got {}): {}",
                line_count,
                fields.len(),
                line
            )
            .into());
        }

        let chromosome = fields[0].to_string();
        let begin: u32 = fields[1].parse().map_err(|e| {
            format!("Error parsing begin position on BED line {}: {} - {}", line_count, line, e)
        })?;
        let end: u32 = fields[2].parse().map_err(|e| {
            format!("Error parsing end position on BED line {}: {} - {}", line_count, line, e)
        })?;
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
            end,
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
    // Validate file type (should be individual call, not combined)
    if let Some(file_type) = crate::filetype::read_file_type_metadata(file_path)
        && file_type != crate::filetype::FileType::IndividualCall
    {
        return Err(format!(
            "Error: Expected individual call file, but got {:?}. Benchmark requires individual call files with 6 columns (chromosome, begin, end, info, H1, H2).",
            file_type
        )
        .into());
    }

    let file_reader = crate::utils::reader(&file_path.to_string_lossy());
    let mut lines = file_reader.lines();
    let mut records = HashMap::new();
    let mut line_num = 0;

    // Skip metadata lines and read header
    let header_line = crate::utils::skip_metadata_lines(&mut lines, &file_path.to_string_lossy());
    line_num += 1;

    // Validate header format
    let header_fields: Vec<&str> = header_line.split('\t').collect();
    match crate::filetype::validate_str_header(&header_fields) {
        Ok(n_samples) => {
            if n_samples != 1 {
                return Err(format!(
                    "Error: Expected individual call file with 1 sample, but found {} samples. Benchmark requires individual call files.",
                    n_samples
                )
                .into());
            }
        }
        Err(e) => {
            return Err(format!("Error: Invalid inquiSTR file header: {}", e).into());
        }
    }

    // Process data lines
    for line_result in lines {
        let line = line_result?;
        line_num += 1;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 6 {
            return Err(format!(
                "Error: Malformed line {} in inquiSTR file (expected 6 fields, got {}): {}",
                line_num,
                fields.len(),
                line
            )
            .into());
        }

        let chromosome = fields[0].to_string();
        let begin: u32 = fields[1].parse().map_err(|e| {
            format!("Error parsing begin position on line {}: {} - {}", line_num, line, e)
        })?;
        let end: u32 = fields[2].parse().map_err(|e| {
            format!("Error parsing end position on line {}: {} - {}", line_num, line, e)
        })?;
        // Skip info field at index 3

        // Parse H1 and H2, handling NaN values
        let h1 = fields[crate::filetype::STR_FIXED_COLUMNS]
            .parse::<f64>()
            .unwrap_or(f64::NAN);
        let h2 = fields[crate::filetype::STR_FIXED_COLUMNS + 1]
            .parse::<f64>()
            .unwrap_or(f64::NAN);

        let record = InquiSTRRecord { chromosome: chromosome.clone(), begin, end, h1, h2 };

        // Use chromosome:begin as key for matching
        let key = format!("{}:{}", chromosome, begin);
        records.insert(key, record);
    }

    println!("inquiSTR file processed: {} records loaded", records.len());

    Ok(records)
}

/// Parse another inquiSTR individual call file as truth data.
/// Caller is responsible for verifying the file type via metadata before calling this function.
fn parse_inquistr_as_truth(
    file_path: &Path,
    max_locus: Option<u32>,
) -> Result<HashMap<String, TruthRecord>, Box<dyn std::error::Error>> {
    let file_reader = crate::utils::reader(&file_path.to_string_lossy());
    let mut lines = file_reader.lines();
    let mut records = HashMap::new();
    let mut filtered_count = 0;

    // Skip metadata lines and read header
    let header_line = crate::utils::skip_metadata_lines(&mut lines, &file_path.to_string_lossy());

    let header_fields: Vec<&str> = header_line.split('\t').collect();
    match crate::filetype::validate_str_header(&header_fields) {
        Ok(n_samples) => {
            if n_samples != 1 {
                return Err(format!(
                    "Error: Expected individual call file with 1 sample as truth, but found {} samples.",
                    n_samples
                )
                .into());
            }
        }
        Err(e) => {
            return Err(format!("Error: Invalid inquiSTR truth file header: {}", e).into());
        }
    }

    for line_result in lines {
        let line = line_result?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 6 {
            return Err(format!(
                "Malformed line in inquiSTR truth file (expected 6 fields, got {}): {}",
                fields.len(),
                line
            )
            .into());
        }

        let chromosome = fields[0].to_string();
        let begin: u32 = fields[1].parse()?;
        let end: u32 = fields[2].parse()?;

        // Filter by locus size if requested
        if let Some(max_size) = max_locus
            && end - begin > max_size
        {
            filtered_count += 1;
            continue;
        }

        let h1 = fields[crate::filetype::STR_FIXED_COLUMNS]
            .parse::<f64>()
            .unwrap_or(f64::NAN);
        let h2 = fields[crate::filetype::STR_FIXED_COLUMNS + 1]
            .parse::<f64>()
            .unwrap_or(f64::NAN);

        let record = TruthRecord {
            chromosome: chromosome.clone(),
            pos: begin,
            end,
            h1,
            h2,
            tier: "Tier1".to_string(),
        };
        let key = format!("{}:{}", chromosome, begin);
        records.insert(key, record);
    }

    if filtered_count > 0 {
        println!(
            "inquiSTR truth file: Filtered out {} intervals larger than {} bp (max-locus limit)",
            filtered_count,
            max_locus.unwrap()
        );
    }
    println!("inquiSTR truth file processed: {} records loaded", records.len());

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

        // VCF POS is 1-based and inclusive; `pos`/`end` are the same interval in the 0-based
        // half-open space the test file uses, so `--match-overlap` compares like with like.
        // The exact-match key below deliberately keeps the 1-based POS it has always used.
        let record = TruthRecord {
            chromosome: chromosome.clone(),
            pos: pos.saturating_sub(1),
            end: pos.saturating_sub(1) + ref_seq.len() as u32,
            h1,
            h2,
            tier: "Tier1".to_string(),
        };

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

/// How test loci were paired with truth loci, before genotypes are looked at.
#[derive(Default)]
struct Pairing<'a> {
    pairs: Vec<(&'a InquiSTRRecord, &'a TruthRecord)>,
    /// Test loci with no truth locus at that position at all.
    test_only: usize,
    /// Test loci whose only truth partner is a Tier2 variant (`--tier1` only).
    tier2_only: usize,
    /// Test loci overlapping a truth locus too weakly to accept (`--match-overlap` only).
    below_min_overlap: usize,
    /// Test loci dropped because a better-overlapping test locus claimed the same truth
    /// locus (`--match-overlap` only), so one truth locus is never counted twice.
    contested: usize,
}

/// Pair test and truth loci on an identical `chromosome:begin` key.
fn pair_on_exact_start<'a>(
    test_records: &'a HashMap<String, InquiSTRRecord>,
    truth_records: &'a HashMap<String, TruthRecord>,
    all_truth_records: &'a HashMap<String, TruthRecord>,
    tier1_only: bool,
) -> Pairing<'a> {
    let mut pairing = Pairing::default();
    for (key, test_record) in test_records {
        if let Some(truth_record) = truth_records.get(key) {
            pairing.pairs.push((test_record, truth_record));
        } else if tier1_only && all_truth_records.contains_key(key) {
            pairing.tier2_only += 1;
        } else {
            pairing.test_only += 1;
        }
    }
    pairing
}

/// Pair each test locus with the truth locus it shares the most bases with.
///
/// Catalogs and truth sets that describe the same repeats rarely agree on exact interval
/// boundaries, which makes matching on the start coordinate alone discard most loci. A pair is
/// kept only if its reciprocal overlap reaches `min_reciprocal_overlap`, and each truth locus is
/// awarded to its single best test locus so that no truth genotype is scored more than once.
fn pair_on_overlap<'a>(
    test_records: &'a HashMap<String, InquiSTRRecord>,
    truth_records: &'a HashMap<String, TruthRecord>,
    all_truth_records: &'a HashMap<String, TruthRecord>,
    tier1_only: bool,
    min_reciprocal_overlap: f64,
) -> Pairing<'a> {
    let index = IntervalIndex::build(truth_records, |_| true);
    let tier2_index = tier1_only
        .then(|| IntervalIndex::build(all_truth_records, |record| record.tier != "Tier1"));

    let mut pairing = Pairing::default();
    // Best test locus per truth locus, keyed by truth key: (reciprocal overlap, test key).
    let mut claims: HashMap<&str, (f64, &str)> = HashMap::new();

    for (test_key, test_record) in test_records {
        let best = index.best_overlap(&test_record.chromosome, test_record.begin, test_record.end);
        let Some((truth_key, reciprocal)) = best else {
            let overlaps_tier2 = tier2_index.as_ref().is_some_and(|idx| {
                idx.best_overlap(&test_record.chromosome, test_record.begin, test_record.end)
                    .is_some()
            });
            if overlaps_tier2 {
                pairing.tier2_only += 1;
            } else {
                pairing.test_only += 1;
            }
            continue;
        };

        if reciprocal < min_reciprocal_overlap {
            pairing.below_min_overlap += 1;
            continue;
        }

        // Keep the strongest claim; fall back to the key so hash order cannot change the result.
        match claims.get(truth_key) {
            Some(&(best_reciprocal, best_key))
                if (best_reciprocal, best_key) >= (reciprocal, test_key.as_str()) =>
            {
                pairing.contested += 1;
            }
            Some(_) => {
                pairing.contested += 1;
                claims.insert(truth_key, (reciprocal, test_key.as_str()));
            }
            None => {
                claims.insert(truth_key, (reciprocal, test_key.as_str()));
            }
        }
    }

    pairing.pairs = claims
        .into_iter()
        .map(|(truth_key, (_, test_key))| (&test_records[test_key], &truth_records[truth_key]))
        .collect();
    pairing
}

/// Main benchmark function
#[allow(clippy::too_many_arguments)]
pub fn benchmark(
    inquistr_file: PathBuf,
    truth_file: PathBuf,
    mode: String,
    plot_file: Option<PathBuf>,
    max_plot_length: f64,
    tier1_only: bool,
    diff_out: Option<PathBuf>,
    max_locus: Option<u32>,
    nonzero: bool,
    tolerance: u32,
    match_overlap: bool,
    min_reciprocal_overlap: f64,
) {
    println!("Loading inquiSTR file: {}", inquistr_file.display());
    let inquistr_records = match parse_inquistr_file(&inquistr_file) {
        Ok(records) => records,
        Err(e) => {
            eprintln!("Error parsing inquiSTR file: {}", e);
            std::process::exit(1);
        }
    };

    // Auto-detect truth file format:
    //   1. VCF if extension is .vcf or ends with .vcf.gz
    //   2. inquiSTR individual call file if it carries inquiSTR metadata
    //   3. Adotto-style BED file otherwise
    let truth_name = truth_file.to_string_lossy();
    let is_vcf = truth_name.ends_with(".vcf") || truth_name.ends_with(".vcf.gz");

    let truth_records = if is_vcf {
        println!("Detected VCF file as truth: {}", truth_file.display());
        match parse_vcf_file(&truth_file, max_locus) {
            Ok(records) => records,
            Err(e) => {
                eprintln!("Error parsing VCF file: {}", e);
                std::process::exit(1);
            }
        }
    } else if let Some(file_type) = crate::filetype::read_file_type_metadata(&truth_file)
        && file_type == crate::filetype::FileType::IndividualCall
    {
        println!("Detected inquiSTR individual call file as truth: {}", truth_file.display());
        match parse_inquistr_as_truth(&truth_file, max_locus) {
            Ok(records) => records,
            Err(e) => {
                eprintln!("Error parsing inquiSTR truth file: {}", e);
                std::process::exit(1);
            }
        }
    } else {
        println!("Detected BED file as truth: {}", truth_file.display());
        match parse_bed_file(&truth_file, max_locus) {
            Ok(records) => records,
            Err(e) => {
                eprintln!("Error parsing BED file: {}", e);
                std::process::exit(1);
            }
        }
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

    // Pair loci first, then look at the genotypes of each pair. Keeping the two steps apart is
    // what lets every test locus be accounted for exactly once below.
    let pairing = if match_overlap {
        println!(
            "Matching loci by interval overlap (minimum reciprocal overlap: {:.2})",
            min_reciprocal_overlap
        );
        pair_on_overlap(
            &inquistr_records,
            &truth_records_filtered,
            &all_truth_records,
            tier1_only,
            min_reciprocal_overlap,
        )
    } else {
        println!("Matching loci on an exact chromosome:begin coordinate");
        pair_on_exact_start(
            &inquistr_records,
            &truth_records_filtered,
            &all_truth_records,
            tier1_only,
        )
    };

    // Collect data for correlation
    let mut inquistr_values = Vec::new();
    let mut truth_values = Vec::new();
    let mut loci_info = Vec::new(); // Store (chromosome, begin, end) for each matched locus
    let mut nan_count = 0; // Paired, but the test locus has no callable allele
    let mut truth_nan_count = 0; // Paired, but the truth locus has no callable allele
    let mut zero_zero_filtered = 0; // Track zero-zero pairs filtered out by --nonzero flag

    for (inquistr_record, truth_record) in &pairing.pairs {
        // Calculate inquiSTR value (MAX or MIN of H1 and H2)
        let inquistr_alleles = vec![inquistr_record.h1, inquistr_record.h2];
        let Some(inquistr_value) = select_allele(&inquistr_alleles, &mode) else {
            nan_count += 1;
            continue;
        };

        // Calculate truth value (MAX or MIN of H1 and H2)
        let truth_alleles = vec![truth_record.h1, truth_record.h2];
        let Some(truth_value) = select_allele(&truth_alleles, &mode) else {
            truth_nan_count += 1;
            continue;
        };

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
    }

    let matched_count = inquistr_values.len();
    // Both pairing strategies award each truth locus to at most one test locus, so the number
    // of pairs is also the number of truth loci that found a partner.
    let truth_only = truth_records_filtered.len() - pairing.pairs.len();

    println!("\n=== Loci Assessment Summary ===");
    println!("Total inquiSTR loci: {}", inquistr_records.len());
    println!("Total truth loci (after max-locus filtering): {}", truth_records_filtered.len());
    println!("\nBreakdown of the {} inquiSTR loci:", inquistr_records.len());
    println!("  Paired with a truth locus: {}", pairing.pairs.len());
    println!("    - assessed: {}", matched_count);
    println!("    - uncalled in inquiSTR (NaN): {}", nan_count);
    println!("    - uncalled in truth data (NaN): {}", truth_nan_count);
    if nonzero {
        println!("    - zero-zero pairs excluded (--nonzero): {}", zero_zero_filtered);
    }
    if pairing.tier2_only > 0 {
        println!("  Matching Tier2 variants only: {}", pairing.tier2_only);
    }
    if match_overlap {
        println!(
            "  Overlapping a truth locus below the minimum reciprocal overlap: {}",
            pairing.below_min_overlap
        );
        println!(
            "  Dropped, truth locus claimed by a better-overlapping locus: {}",
            pairing.contested
        );
    }
    println!("  Without a truth locus: {}", pairing.test_only);
    println!("\nBreakdown of the {} truth loci:", truth_records_filtered.len());
    println!("  Paired with an inquiSTR locus: {}", pairing.pairs.len());
    println!("  Without an inquiSTR locus: {}", truth_only);
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
    let mut off_by_one = 0;
    let mut within_tolerance = 0;

    for (&inq_val, &truth_val) in inquistr_values.iter().zip(truth_values.iter()) {
        let diff = (inq_val - truth_val).abs();
        if diff == 0.0 {
            exact_matches += 1;
        }
        // A genotype is a median, so it can be a half-integer and a difference can be 0.5.
        // These are thresholds rather than equality tests, or such a pair counts as neither
        // exact nor off-by-one. --tolerance is at least 1, so each count nests in the next.
        if diff <= 1.0 {
            off_by_one += 1;
        }
        if diff <= f64::from(tolerance) {
            within_tolerance += 1;
        }
    }

    let exact_percent = (exact_matches as f64 / matched_count as f64) * 100.0;
    let off_by_one_percent = (off_by_one as f64 / matched_count as f64) * 100.0;
    let within_tolerance_percent = (within_tolerance as f64 / matched_count as f64) * 100.0;

    println!("\n=== Accuracy Analysis ===");
    println!("Exact matches: {} ({:.2}%)", exact_matches, exact_percent);
    println!("Maximally off by one: {} ({:.2}%)", off_by_one, off_by_one_percent);
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
                    "chromosome\tbegin\tend\ttest_genotype\ttruth_genotype\tabsolute_difference"
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
    let mut plot_data: Vec<(f64, f64)> = Vec::new();
    let mut hidden_count = 0;

    for (inq_val, truth_val) in inquistr_values.iter().zip(truth_values.iter()) {
        if inq_val.abs() <= max_plot_length && truth_val.abs() <= max_plot_length {
            plot_data.push((*inq_val, *truth_val));
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

    // Create scatter plot with trend line and R²
    let trace = ScatterPlot::new()
        .with_data(plot_data)
        .with_color("steelblue")
        .with_size(3.0)
        .with_marker_opacity(0.5)
        .with_trend(TrendLine::Linear)
        .with_correlation()
        .with_legend("Data points");

    let title_text = if nonzero {
        format!(
            "inquiSTR vs Truth Genotypes (Mode: {}, Nonzero only, R² = {:.4})",
            mode, r_squared_all
        )
    } else {
        format!("inquiSTR vs Truth Genotypes (Mode: {}, R² = {:.4})", mode, r_squared_all)
    };

    let layout = Layout::new((-axis_limit, axis_limit), (-axis_limit, axis_limit))
        .with_title(title_text)
        .with_x_label("inquiSTR genotypes")
        .with_y_label("Truth genotypes")
        .with_width(800.0)
        .with_height(800.0)
        .with_interactive();

    // Save plot if output file is specified
    if let Some(plot_path) = plot_file {
        let plots: Vec<Plot> = vec![trace.into()];
        let svg = render_to_svg(plots, layout);
        match std::fs::write(&plot_path, svg) {
            Ok(_) => println!("Plot saved to: {}", plot_path.display()),
            Err(e) => {
                eprintln!("Error saving plot: {}", e);
                std::process::exit(1);
            }
        }
    } else {
        println!("Skipping plot generation (no --plot file specified)");
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
    println!(
        "MATCH_MODE: {}",
        if match_overlap {
            "overlap"
        } else {
            "exact_start"
        }
    );
    if match_overlap {
        println!("MIN_RECIPROCAL_OVERLAP: {:.2}", min_reciprocal_overlap);
    }
    println!("LOCI_PAIRED: {}", pairing.pairs.len());
    println!("TEST_NOCALL: {}", nan_count);
    println!("TRUTH_NOCALL: {}", truth_nan_count);
    println!("TEST_ONLY: {}", pairing.test_only);
    println!("TRUTH_ONLY: {}", truth_only);
    println!("=========================");
}

#[cfg(test)]
mod tests {
    use super::*;

    fn truth(chromosome: &str, pos: u32, end: u32, tier: &str) -> TruthRecord {
        TruthRecord {
            chromosome: chromosome.to_string(),
            pos,
            end,
            h1: 0.0,
            h2: 0.0,
            tier: tier.to_string(),
        }
    }

    fn test_locus(chromosome: &str, begin: u32, end: u32) -> InquiSTRRecord {
        InquiSTRRecord { chromosome: chromosome.to_string(), begin, end, h1: 0.0, h2: 0.0 }
    }

    fn truth_map(records: Vec<(&str, TruthRecord)>) -> HashMap<String, TruthRecord> {
        records
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect()
    }

    fn test_map(records: Vec<(&str, InquiSTRRecord)>) -> HashMap<String, InquiSTRRecord> {
        records
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect()
    }

    #[test]
    fn best_overlap_prefers_the_largest_shared_span() {
        let records = truth_map(vec![
            ("chr1:100", truth("chr1", 100, 160, "Tier1")),
            ("chr1:150", truth("chr1", 150, 400, "Tier1")),
        ]);
        let index = IntervalIndex::build(&records, |_| true);

        // [140, 400) shares 20 bp with the first and 250 bp with the second; the reciprocal
        // overlap is then 250 shared bases over the 260 bp query, the longer of the two.
        let (key, reciprocal) = index.best_overlap("chr1", 140, 400).unwrap();
        assert_eq!(key, "chr1:150");
        assert!((reciprocal - 250.0 / 260.0).abs() < 1e-9);
    }

    #[test]
    fn best_overlap_reports_the_fraction_of_the_longer_interval() {
        let records = truth_map(vec![("chr1:100", truth("chr1", 100, 110, "Tier1"))]);
        let index = IntervalIndex::build(&records, |_| true);

        // A 10 bp truth locus fully inside a 100 bp test locus: 10 / 100.
        let (_, reciprocal) = index.best_overlap("chr1", 50, 150).unwrap();
        assert!((reciprocal - 0.1).abs() < 1e-9);
    }

    #[test]
    fn best_overlap_ignores_touching_and_distant_intervals() {
        let records = truth_map(vec![("chr1:100", truth("chr1", 100, 200, "Tier1"))]);
        let index = IntervalIndex::build(&records, |_| true);

        assert!(index.best_overlap("chr1", 200, 300).is_none(), "half-open ends must not touch");
        assert!(index.best_overlap("chr1", 0, 100).is_none(), "half-open starts must not touch");
        assert!(index.best_overlap("chr2", 100, 200).is_none(), "other chromosome");
    }

    #[test]
    fn exact_pairing_separates_tier2_and_unmatched_loci() {
        let tests = test_map(vec![
            ("chr1:100", test_locus("chr1", 100, 200)),
            ("chr1:300", test_locus("chr1", 300, 400)),
            ("chr1:500", test_locus("chr1", 500, 600)),
        ]);
        let all_truth = truth_map(vec![
            ("chr1:100", truth("chr1", 100, 200, "Tier1")),
            ("chr1:300", truth("chr1", 300, 400, "Tier2")),
        ]);
        let tier1: HashMap<String, TruthRecord> = all_truth
            .iter()
            .filter(|(_, r)| r.tier == "Tier1")
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();

        let pairing = pair_on_exact_start(&tests, &tier1, &all_truth, true);
        assert_eq!(pairing.pairs.len(), 1);
        assert_eq!(pairing.tier2_only, 1);
        assert_eq!(pairing.test_only, 1);
    }

    #[test]
    fn overlap_pairing_rejects_weak_overlaps() {
        let tests = test_map(vec![("chr1:0", test_locus("chr1", 0, 1000))]);
        // 10 bp of truth inside a 1000 bp test locus is a reciprocal overlap of 0.01.
        let truth_records = truth_map(vec![("chr1:100", truth("chr1", 100, 110, "Tier1"))]);
        let all_truth = truth_records.clone();

        let accepted = pair_on_overlap(&tests, &truth_records, &all_truth, false, 0.005);
        assert_eq!(accepted.pairs.len(), 1);
        assert_eq!(accepted.below_min_overlap, 0);

        let rejected = pair_on_overlap(&tests, &truth_records, &all_truth, false, 0.5);
        assert_eq!(rejected.pairs.len(), 0);
        assert_eq!(rejected.below_min_overlap, 1);
        assert_eq!(rejected.test_only, 0, "a weak overlap is not the same as no overlap");
    }

    #[test]
    fn overlap_pairing_awards_a_truth_locus_to_one_test_locus_only() {
        // Both test loci overlap the same truth locus; the second one matches it exactly.
        let tests = test_map(vec![
            ("chr1:100", test_locus("chr1", 100, 260)),
            ("chr1:200", test_locus("chr1", 200, 300)),
        ]);
        let truth_records = truth_map(vec![("chr1:200", truth("chr1", 200, 300, "Tier1"))]);
        let all_truth = truth_records.clone();

        let pairing = pair_on_overlap(&tests, &truth_records, &all_truth, false, 0.1);
        assert_eq!(pairing.pairs.len(), 1, "the truth locus must not be scored twice");
        assert_eq!(pairing.pairs[0].0.begin, 200, "the better overlap wins");
        assert_eq!(pairing.contested, 1);
    }

    #[test]
    fn overlap_pairing_accounts_for_every_test_locus_exactly_once() {
        let tests = test_map(vec![
            ("chr1:100", test_locus("chr1", 100, 200)), // exact overlap
            ("chr1:150", test_locus("chr1", 150, 250)), // contests the same truth locus
            ("chr1:400", test_locus("chr1", 400, 4000)), // overlaps too weakly
            ("chr1:9000", test_locus("chr1", 9000, 9100)), // no truth locus at all
        ]);
        let truth_records = truth_map(vec![
            ("chr1:100", truth("chr1", 100, 200, "Tier1")),
            ("chr1:410", truth("chr1", 410, 420, "Tier1")),
        ]);
        let all_truth = truth_records.clone();

        let pairing = pair_on_overlap(&tests, &truth_records, &all_truth, false, 0.5);
        let total = pairing.pairs.len()
            + pairing.tier2_only
            + pairing.below_min_overlap
            + pairing.contested
            + pairing.test_only;
        assert_eq!(total, tests.len(), "the breakdown must partition the test loci");
        assert_eq!(pairing.pairs.len(), 1);
        assert_eq!(pairing.contested, 1);
        assert_eq!(pairing.below_min_overlap, 1);
        assert_eq!(pairing.test_only, 1);
    }
}
