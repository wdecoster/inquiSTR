use bio::io::bed;
use std::collections::HashMap;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

/// Simple BED interval for filtering
#[derive(Debug)]
struct BedInterval {
    start: u32,
    end: u32,
}

/// Parse a BED file and return intervals grouped by chromosome (sorted within each chromosome)
fn parse_bed_file(bed_path: &Path) -> io::Result<HashMap<String, Vec<BedInterval>>> {
    let file_reader = crate::utils::reader(&bed_path.to_string_lossy());
    let mut reader = bed::Reader::new(file_reader);
    let mut intervals: HashMap<String, Vec<BedInterval>> = HashMap::new();

    for record in reader.records() {
        let rec = record.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let chrom = rec.chrom().to_string();
        let start: u32 = rec.start().try_into().unwrap();
        let end: u32 = rec.end().try_into().unwrap();

        intervals
            .entry(chrom.clone())
            .or_default()
            .push(BedInterval { start, end });
    }

    Ok(intervals)
}

/// Check if a locus overlaps with any BED interval using binary search for efficiency
/// Both the locus file and BED file are sorted, so we can skip intervals that are before the locus
fn overlaps_bed(
    chrom: &str,
    start: u32,
    end: u32,
    bed_intervals: &HashMap<String, Vec<BedInterval>>,
    bed_index: &mut HashMap<String, usize>,
) -> bool {
    if let Some(intervals) = bed_intervals.get(chrom) {
        let idx = bed_index.entry(chrom.to_string()).or_insert(0);

        // Start from the last position we checked for this chromosome
        // Skip intervals that end before our locus starts
        while *idx < intervals.len() && intervals[*idx].end <= start {
            *idx += 1;
        }

        // Check remaining intervals until we pass the locus
        let mut check_idx = *idx;
        while check_idx < intervals.len() {
            let interval = &intervals[check_idx];

            // If this interval starts after our locus ends, no more overlaps possible
            if interval.start >= end {
                break;
            }

            // Check for overlap: locus.start < interval.end && locus.end > interval.start
            if start < interval.end && end > interval.start {
                return true;
            }

            check_idx += 1;
        }
    }
    false
}

/// Verify that a file is sorted by chromosome and position
fn verify_sorted<R: BufRead>(reader: &mut R, file_name: &str, skip_header: bool) -> io::Result<()> {
    let mut prev_chrom: Option<String> = None;
    let mut prev_start: Option<u32> = None;
    let mut line_num = 0;

    for line in reader.lines() {
        let line = line?;
        line_num += 1;

        // Skip header line if requested
        if line_num == 1 && skip_header {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let chrom = fields[0].to_string();
        let start: u32 = fields[1].parse().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid start position at line {} in {}", line_num, file_name),
            )
        })?;

        // Check if sorted
        if let Some(ref prev_chr) = prev_chrom
            && chrom == *prev_chr
            && let Some(prev_pos) = prev_start
            && start < prev_pos
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "File {} is not sorted: found {}:{} after {}:{}",
                    file_name, chrom, start, prev_chr, prev_pos
                ),
            ));
        }

        prev_chrom = Some(chrom);
        prev_start = Some(start);
    }

    Ok(())
}

/// Calculate coefficient of variation for a set of values
fn calculate_cv(values: &[f64]) -> f64 {
    let valid_values: Vec<f64> = values
        .iter()
        .filter(|v| v.is_finite() && **v > 0.0)
        .copied()
        .collect();

    if valid_values.len() < 2 {
        return 0.0;
    }

    let mean: f64 = valid_values.iter().sum::<f64>() / valid_values.len() as f64;

    if mean == 0.0 {
        return 0.0;
    }

    let variance: f64 =
        valid_values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / valid_values.len() as f64;

    let stdev = variance.sqrt();
    stdev / mean
}

/// Count the total number of lines in a file (for progress indicator)
fn count_lines(file_path: &Path) -> io::Result<usize> {
    let reader = crate::utils::reader(&file_path.to_string_lossy());
    Ok(reader.lines().count())
}

/// Verify file has a valid header (must start with "chromosome")
fn verify_header(file_path: &Path) -> io::Result<()> {
    let reader = crate::utils::reader(&file_path.to_string_lossy());
    let mut lines = reader.lines();

    // Skip metadata lines starting with #
    let first_line = crate::utils::skip_metadata_lines(&mut lines);

    if !first_line.starts_with("chromosome") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Invalid inquiSTR file: missing header. First line should start with 'chromosome', got: {}",
                first_line.split('\t').next().unwrap_or("<empty>")
            ),
        ));
    }
    Ok(())
}

/// Determine if input is a single-sample or combined file based on number of columns
/// First tries to read file_type metadata, then falls back to heuristics for backward compatibility
fn is_combined_file(file_path: &Path) -> io::Result<bool> {
    // Try reading metadata first
    if let Some(file_type) = crate::filetype::read_file_type_metadata(file_path) {
        return Ok(file_type == crate::filetype::FileType::CombinedCall);
    }

    // Fall back to heuristic detection for files without metadata
    let reader = crate::utils::reader(&file_path.to_string_lossy());
    let mut lines = reader.lines();

    // Skip metadata lines starting with #
    let first_line = crate::utils::skip_metadata_lines(&mut lines);

    // Count number of columns
    // Single sample: chr, start, end, info, H1, H2 = 6 columns
    // Combined: chr, start, end, info, sample1_H1, sample1_H2, sample2_H1, sample2_H2, ... = > 6 columns
    let num_cols = first_line.trim_end().split('\t').count();
    Ok(num_cols > 6)
}

pub fn filter(
    input: PathBuf,
    minlen: Option<u32>,
    minchange: Option<u32>,
    bed: Option<PathBuf>,
    call_rate: Option<f64>,
    min_cv: Option<f64>,
) {
    // Validate call_rate range
    if let Some(rate) = call_rate
        && !(0.0..=1.0).contains(&rate)
    {
        eprintln!("ERROR: call-rate must be between 0.0 and 1.0");
        std::process::exit(1);
    }

    // Check if input exists
    if !input.exists() {
        eprintln!("ERROR: Input file does not exist: {}", input.display());
        std::process::exit(1);
    }

    // Verify input file has a valid header
    verify_header(&input).expect("Invalid input file format");

    // Determine if single-sample or combined file
    let is_combined = is_combined_file(&input).expect("Failed to read input file");

    // Validate min_cv is only used with combined files
    if min_cv.is_some() && !is_combined {
        eprintln!("ERROR: --min-cv filter can only be used with combined files (multi-sample)");
        std::process::exit(1);
    }

    // Warn about call_rate with single-sample files
    if call_rate.is_some() && !is_combined {
        eprintln!("WARNING: --call-rate filter is ignored for single-sample files");
    }

    // Verify input file is sorted if bed filter is specified
    if bed.is_some() {
        eprintln!("Verifying input file is sorted...");
        let mut reader = crate::utils::reader(&input.to_string_lossy());
        verify_sorted(&mut reader, &input.to_string_lossy(), true)
            .expect("Input file is not sorted by chromosome and position");
        eprintln!("✓ Input file is sorted");
    }

    // Load and verify BED file if specified
    let bed_intervals = if let Some(ref bed_path) = bed {
        if !bed_path.exists() {
            eprintln!("ERROR: BED file does not exist: {}", bed_path.display());
            std::process::exit(1);
        }

        eprintln!("Verifying BED file is sorted...");
        let mut bed_reader = crate::utils::reader(&bed_path.to_string_lossy());
        verify_sorted(&mut bed_reader, &bed_path.to_string_lossy(), false)
            .expect("BED file is not sorted by chromosome and position");
        eprintln!("✓ BED file is sorted");

        eprintln!("Loading BED intervals...");
        let intervals = parse_bed_file(bed_path).expect("Failed to parse BED file");
        let total_intervals: usize = intervals.values().map(|v| v.len()).sum();
        eprintln!("✓ Loaded {} intervals from {} chromosomes", total_intervals, intervals.len());
        Some(intervals)
    } else {
        None
    };

    // Count total lines for progress indicator
    eprintln!("Counting lines...");
    let total_lines = count_lines(&input).expect("Failed to count lines");
    eprintln!("✓ Processing {} lines", total_lines);

    // Initialize counters for statistics
    let mut total_loci = 0;
    let mut filtered_by_bed = 0;
    let mut filtered_by_minlen = 0;
    let mut filtered_by_minchange = 0;
    let mut filtered_by_call_rate = 0;
    let mut filtered_by_cv = 0;
    let mut passed_loci = 0;

    // Initialize BED interval index for efficient sorted traversal
    let mut bed_index: HashMap<String, usize> = HashMap::new();

    // Process file line by line
    let reader = crate::utils::reader(&input.to_string_lossy());
    let stdout = io::stdout();
    let mut writer = io::BufWriter::new(stdout.lock());

    let mut line_num = 0;
    let progress_interval = (total_lines / 100).max(1000); // Update every 1% or 1000 lines
    let mut header_written = false;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        line_num += 1;

        // Print progress
        if line_num % progress_interval == 0 {
            eprint!(
                "\rProcessing: {}/{} ({:.1}%)",
                line_num,
                total_lines,
                (line_num as f64 / total_lines as f64) * 100.0
            );
        }

        // Skip and write metadata lines
        if line.starts_with('#') {
            writeln!(writer, "{}", line).expect("Failed to write metadata");
            continue;
        }

        // Write header (first non-metadata line)
        if !header_written {
            // Validate and store expected field count from header
            let header_fields: Vec<&str> = line.split('\t').collect();
            if header_fields.len() < 6 {
                eprintln!("ERROR: Invalid header format");
                eprintln!("Expected at least 6 columns (chromosome, begin, end, info, sample_H1, sample_H2)");
                eprintln!("Got {} columns", header_fields.len());
                std::process::exit(1);
            }
            
            // Use validation function if this is a STR file with proper header
            if header_fields[0] == "chromosome"
                && let Err(e) = crate::filetype::validate_str_header(&header_fields) {
                    eprintln!("ERROR: Invalid STR file header: {}", e);
                    std::process::exit(1);
                }
            writeln!(writer, "{}", line).expect("Failed to write header");
            header_written = true;
            continue;
        }

        // Parse line
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            eprintln!("ERROR: Invalid line {} in input file", line_num);
            eprintln!(
                "Expected at least 6 columns (chr, start, end, info, H1, H2), got {}",
                fields.len()
            );
            eprintln!("Line content: '{}'", line);
            eprintln!("\nThe file may be corrupted or in an incorrect format.");
            std::process::exit(1);
        }

        total_loci += 1;

        let chrom = fields[0];
        let start: u32 = fields[1].parse().expect("Invalid start position");
        let end: u32 = fields[2].parse().expect("Invalid end position");
        // Skip info field at index 3

        // Apply filters in order

        // 1. BED overlap filter (using efficient sorted traversal)
        if let Some(ref intervals) = bed_intervals
            && !overlaps_bed(chrom, start, end, intervals, &mut bed_index)
        {
            filtered_by_bed += 1;
            continue;
        }

        // Parse allele values (columns 4 onwards, skipping info at column 3)
        let alleles: Vec<&str> = fields[4..].to_vec();

        // 2. minlen filter (only positive values >= minlen)
        if let Some(min) = minlen {
            let has_expansion = alleles.iter().any(|&allele| {
                if let Ok(val) = allele.parse::<f64>() {
                    val.is_finite() && val >= min as f64
                } else {
                    false
                }
            });

            if !has_expansion {
                filtered_by_minlen += 1;
                continue;
            }
        }

        // 3. minchange filter (absolute value >= minchange)
        if let Some(min) = minchange {
            let has_change = alleles.iter().any(|&allele| {
                if let Ok(val) = allele.parse::<f64>() {
                    val.is_finite() && val.abs() >= min as f64
                } else {
                    false
                }
            });

            if !has_change {
                filtered_by_minchange += 1;
                continue;
            }
        }

        // 4. Call rate filter (only for combined files)
        if let Some(rate) = call_rate
            && is_combined
        {
            let total_samples = alleles.len();
            let genotyped_samples = alleles
                .iter()
                .filter(|&&a| {
                    if let Ok(val) = a.parse::<f64>() {
                        val.is_finite()
                    } else {
                        false
                    }
                })
                .count();

            let actual_rate = genotyped_samples as f64 / total_samples as f64;
            if actual_rate < rate {
                filtered_by_call_rate += 1;
                continue;
            }
        }

        // 5. Coefficient of variation filter (only for combined files)
        if let Some(min_cv_threshold) = min_cv
            && is_combined
        {
            let values: Vec<f64> = alleles
                .iter()
                .filter_map(|&a| a.parse::<f64>().ok())
                .collect();

            let cv = calculate_cv(&values);
            if cv < min_cv_threshold {
                filtered_by_cv += 1;
                continue;
            }
        }

        // Locus passed all filters
        passed_loci += 1;
        writeln!(writer, "{}", line).expect("Failed to write output");
    }

    eprint!("\r"); // Clear progress line

    // Print statistics to stderr
    eprintln!("\n=== Filtering Statistics ===");
    eprintln!("Total loci processed: {}", total_loci);
    if bed.is_some() {
        eprintln!(
            "  Filtered by BED overlap: {} ({:.1}%)",
            filtered_by_bed,
            (filtered_by_bed as f64 / total_loci as f64) * 100.0
        );
    }
    if minlen.is_some() {
        eprintln!(
            "  Filtered by minlen: {} ({:.1}%)",
            filtered_by_minlen,
            (filtered_by_minlen as f64 / total_loci as f64) * 100.0
        );
    }
    if minchange.is_some() {
        eprintln!(
            "  Filtered by minchange: {} ({:.1}%)",
            filtered_by_minchange,
            (filtered_by_minchange as f64 / total_loci as f64) * 100.0
        );
    }
    if call_rate.is_some() && is_combined {
        eprintln!(
            "  Filtered by call-rate: {} ({:.1}%)",
            filtered_by_call_rate,
            (filtered_by_call_rate as f64 / total_loci as f64) * 100.0
        );
    }
    if min_cv.is_some() && is_combined {
        eprintln!(
            "  Filtered by min-cv: {} ({:.1}%)",
            filtered_by_cv,
            (filtered_by_cv as f64 / total_loci as f64) * 100.0
        );
    }
    eprintln!(
        "Loci passing filters: {} ({:.1}%)",
        passed_loci,
        (passed_loci as f64 / total_loci as f64) * 100.0
    );
}
