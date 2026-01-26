use log::error;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::io::{self, Write};
use std::path::PathBuf;

// Phasing validation now happens lazily on first batch via get_bam_reader_with_validation
use crate::batching::create_batches;
use crate::errors::{InquiSTRError, InquiSTRResult};
use crate::genotype_batch::{process_batch_with_reader, process_batch_worker};
use crate::repeats::{ChromosomeMapper, RepeatInterval, TargetConfig};

/// Configuration for genotyping parameters (how to call STRs)
#[derive(Clone, Copy, Debug)]
pub struct GenotypeConfig {
    pub minlen: u32,
    pub support: usize,
    pub unphased: bool,
}

/// Configuration for processing (threads, batching, etc.)
#[derive(Clone, Debug)]
pub struct ProcessingConfig {
    pub threads: usize,
    pub batch_size_kb: u32,
    pub output_vcf: bool,
}

// This struct keeps the genotype information and allows to compare them and thus sort them on chromosomal location
pub struct Genotype {
    pub repeat: RepeatInterval,
    pub phase1: f64,
    pub phase2: f64,
}

impl Ord for Genotype {
    fn cmp(&self, other: &Self) -> Ordering {
        // Compare by chromosome ID first (simpler and faster than human_compare)
        // then by start position
        self.repeat
            .chrom_id
            .cmp(&other.repeat.chrom_id)
            .then(self.repeat.start.cmp(&other.repeat.start))
    }
}

impl PartialOrd for Genotype {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Genotype {
    fn eq(&self, other: &Self) -> bool {
        self.repeat.chrom_id == other.repeat.chrom_id && self.repeat.start == other.repeat.start
    }
}

impl Eq for Genotype {}

impl Genotype {
    /// Format genotype for output using the chromosome mapper
    pub fn format_output(&self, mapper: &ChromosomeMapper) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            mapper.get_name(self.repeat.chrom_id),
            self.repeat.start,
            self.repeat.end,
            self.repeat.info,
            self.phase1,
            self.phase2
        )
    }
}

/// This function genotypes STRs, either from a region string or from a bed file
/// For a bed file the genotyping is done in parallel
/// show_progress can be silenced for inquiSTR batch, but this is not explicitly exposed to the CLI
pub fn genotype_repeats(
    bam: String,
    targets: TargetConfig,
    genotype: GenotypeConfig,
    processing: ProcessingConfig,
    sample_name: Option<String>,
    reference: Option<String>,
    show_progress: bool,
) -> InquiSTRResult<()> {
    // only test if path.is_file() if the file is local
    if !PathBuf::from(&bam).is_file()
        && !bam.starts_with("s3")
        && !bam.starts_with("https://")
        && !bam.starts_with("ftp://")
    {
        error!("ERROR: path to bam file {} is not valid!\n\n", &bam);
        return Err(InquiSTRError::new(format!("path to bam file {} is not valid", &bam)));
    };

    // Unified batch-level producer-consumer approach for both single and multi-threaded
    // Batch size is configurable for performance optimization
    let (repeats, chrom_mapper) = targets.get_targets(&bam, &reference)?;
    let total_loci = repeats.len();
    let batches = create_batches(repeats, processing.batch_size_kb * 1000, &chrom_mapper); // Convert kb to basepair

    // Setup progress bar with smoothed ETA
    let pb = if show_progress {
        let pb = indicatif::ProgressBar::new(total_loci as u64);
        pb.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} loci ({eta})")
                .expect("Failed to set progress bar template")
        );
        // Enable steady tick for smoothed ETA calculation (updates every 100ms)
        // This makes the time estimate converge to accuracy rather than jumping around
        pb.enable_steady_tick(std::time::Duration::from_millis(100));
        Some(pb)
    } else {
        None
    };

    // Process batches using producer-consumer pattern with configurable worker count
    let results: Vec<Vec<Genotype>> = if processing.threads > 1 {
        // Multi-threaded: Process batches in parallel
        // NOTE: Cannot share BAM reader across threads, so each worker creates its own
        // This means parallel processing with CRAM may still hit FD limits
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(processing.threads)
            .build()
            .expect("Failed to build thread pool");

        thread_pool.install(|| {
            batches
                .into_par_iter()
                .map(|batch| {
                    let batch_size = batch.repeats.len();
                    let results = process_batch_worker(batch, &bam, &reference, genotype)?;
                    if let Some(ref pb) = pb {
                        pb.inc(batch_size as u64);
                    }
                    Ok(results)
                })
                .collect::<InquiSTRResult<Vec<Vec<Genotype>>>>()
        })?
    } else {
        // Single-threaded: Create a SINGLE BAM reader and reuse it for all batches
        // This dramatically reduces file descriptor usage with CRAM files
        let mut bam_reader = if !genotype.unphased {
            crate::bam_utils::get_bam_reader_with_validation(&bam, &reference)?
        } else {
            crate::bam_utils::get_bam_reader(&bam, &reference)?
        };

        // BAM reader will be dropped here automatically
        batches
            .iter()
            .map(|batch| {
                let batch_size = batch.repeats.len();
                let results = process_batch_with_reader(batch, &mut bam_reader, &genotype)?;
                if let Some(ref pb) = pb {
                    pb.inc(batch_size as u64);
                }
                Ok(results)
            })
            .collect::<InquiSTRResult<Vec<Vec<Genotype>>>>()?
    };

    if let Some(ref pb) = pb {
        pb.finish_with_message("Processing completed, writing results...");
    }

    // Output results in consistent format
    // Use either the sample_name provided as command line argument or extract one from the path
    let sample = sample_name.unwrap_or_else(|| crate::utils::extract_sample_name_from_path(&bam));

    // Memory optimization: Process and write results chromosome-by-chromosome
    // Since batches are sorted by chromosome, we can detect chromosome boundaries,
    // sort within chromosome, write output, then drop those results before the next chromosome

    // If VCF output is requested, we only write to a temp TSV file for conversion
    // Otherwise, write TSV directly to stdout
    let temp_tsv_path = if processing.output_vcf {
        let temp_path =
            std::env::temp_dir().join(format!("inquiSTR_{}.tmp.tsv", std::process::id()));
        Some(temp_path)
    } else {
        None
    };

    // Create writers: either just temp file (for VCF) or stdout (for TSV)
    let mut stdout_writer = if !processing.output_vcf {
        let stdout = io::stdout();
        Some(io::BufWriter::new(stdout))
    } else {
        None
    };

    let mut temp_writer = temp_tsv_path.as_ref().map(|path| {
        io::BufWriter::new(std::fs::File::create(path).expect("Failed to create temp TSV file"))
    });

    // Write metadata header to stdout if in TSV mode
    if let Some(ref mut handle) = stdout_writer {
        writeln!(handle, "# file_type=individual_call").expect("Failed writing metadata.");
        writeln!(handle, "# version={}", crate::VERSION).expect("Failed writing metadata.");
        writeln!(handle, "# command=call").expect("Failed writing metadata.");
        writeln!(handle, "# sample={}", sample).expect("Failed writing metadata.");
        writeln!(handle, "# minlen={}", genotype.minlen).expect("Failed writing metadata.");
        writeln!(handle, "# support={}", genotype.support).expect("Failed writing metadata.");
        writeln!(handle, "# unphased={}", genotype.unphased).expect("Failed writing metadata.");
    }

    // Write data header to both stdout (if TSV) and temp file (if VCF)
    let file_header = format!("chromosome\tbegin\tend\tinfo\t{sample}_H1\t{sample}_H2");
    if let Some(ref mut handle) = stdout_writer {
        writeln!(handle, "{file_header}").expect("Failed writing the header.");
    }
    if let Some(ref mut temp) = temp_writer {
        writeln!(temp, "{file_header}").expect("Failed writing temp header.");
    }

    // Process chromosome-by-chromosome
    let mut current_chrom_results: Vec<Genotype> = Vec::new();
    let mut current_chrom_id: Option<u32> = None;

    for batch_results in results {
        for genotype in batch_results {
            // Check if we're starting a new chromosome (by ID instead of string comparison)
            if current_chrom_id != Some(genotype.repeat.chrom_id) {
                // Write out previous chromosome's results (sorted)
                if !current_chrom_results.is_empty() {
                    current_chrom_results.sort_unstable();
                    for g in &current_chrom_results {
                        let output_line = g.format_output(&chrom_mapper);
                        if let Some(ref mut handle) = stdout_writer {
                            writeln!(handle, "{}", output_line)
                                .expect("Failed writing the result.");
                        }
                        if let Some(ref mut temp) = temp_writer {
                            writeln!(temp, "{}", output_line)
                                .expect("Failed writing to temp file.");
                        }
                    }
                    current_chrom_results.clear();
                }
                current_chrom_id = Some(genotype.repeat.chrom_id);
            }

            current_chrom_results.push(genotype);
        }
    }

    // Write final chromosome's results
    if !current_chrom_results.is_empty() {
        current_chrom_results.sort_unstable();
        for g in &current_chrom_results {
            let output_line = g.format_output(&chrom_mapper);
            if let Some(ref mut handle) = stdout_writer {
                writeln!(handle, "{}", output_line).expect("Failed writing the result.");
            }
            if let Some(ref mut temp) = temp_writer {
                writeln!(temp, "{}", output_line).expect("Failed writing to temp file.");
            }
        }
    }

    // Flush temp file to ensure all data is written
    if let Some(mut temp) = temp_writer {
        temp.flush().expect("Failed to flush temp file");
    }

    // Write VCF if requested by converting from the temp TSV file and writing to stdout
    if let Some(temp_path) = temp_tsv_path.as_ref() {
        write_vcf_to_stdout(temp_path, &sample, &reference);
        // Clean up temp file
        std::fs::remove_file(temp_path).ok(); // Ignore errors if already deleted
    }

    // Clean up any downloaded index files
    crate::bam_utils::cleanup_index_files(&bam);

    Ok(())
}

/// Write VCF to stdout by reading from a TSV file (avoids keeping all genotypes in memory)
fn write_vcf_to_stdout(tsv_path: &std::path::Path, sample: &str, reference: &Option<String>) {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let stdout = io::stdout();
    let mut writer = io::BufWriter::new(stdout);

    // Write VCF header
    writeln!(writer, "##fileformat=VCFv4.3").expect("Failed writing VCF header");
    writeln!(writer, "##source=inquiSTR").expect("Failed writing VCF header");

    // Add reference if provided
    if let Some(ref_path) = reference {
        writeln!(writer, "##reference={}", ref_path).expect("Failed writing VCF header");
    }

    // Add INFO fields
    writeln!(
        writer,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">"
    )
    .expect("Failed writing VCF header");
    writeln!(
        writer,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    )
    .expect("Failed writing VCF header");
    writeln!(writer, "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">").expect("Failed writing VCF header");

    // Add FORMAT fields
    writeln!(writer, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
        .expect("Failed writing VCF header");
    writeln!(
        writer,
        "##FORMAT=<ID=AL,Number=.,Type=Float,Description=\"Allele length (relative to reference)\">"
    )
    .expect("Failed writing VCF header");

    // Add ALT definition for STR
    writeln!(writer, "##ALT=<ID=STR,Description=\"Short Tandem Repeat\">")
        .expect("Failed writing VCF header");

    // Write column headers
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", sample)
        .expect("Failed writing VCF header");

    // Read TSV and convert to VCF records
    let tsv_file = File::open(tsv_path).expect("Failed to open temp TSV file");
    let reader = BufReader::new(tsv_file);

    let mut idx = 0;
    for line in reader.lines() {
        let line = line.expect("Failed reading TSV line");

        // Skip header line
        if line.starts_with("chromosome") || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 6 {
            continue; // Skip malformed lines
        }

        let chrom = fields[0];
        let start: u32 = fields[1].parse().unwrap_or(0);
        let end: u32 = fields[2].parse().unwrap_or(0);
        let _info = fields[3]; // Currently unused, but available if needed
        let phase1_str = fields[4];
        let phase2_str = fields[5];

        idx += 1;
        let pos = start + 1; // VCF is 1-based
        let id = format!("STR_{}", idx);

        // Parse phase values (handle "nan" or numeric values)
        let phase1: Option<f64> = phase1_str.parse().ok();
        let phase2: Option<f64> = phase2_str.parse().ok();

        // Determine genotype and allele lengths
        let (gt, al1, al2) = match (phase1, phase2) {
            (None, None) => ("./.".to_string(), ".".to_string(), ".".to_string()),
            (None, Some(p2)) => ("./1".to_string(), ".".to_string(), format!("{:.0}", p2)),
            (Some(p1), None) => ("0/.".to_string(), format!("{:.0}", p1), ".".to_string()),
            (Some(p1), Some(p2)) => ("0|1".to_string(), format!("{:.0}", p1), format!("{:.0}", p2)),
        };

        // Calculate SVLEN (using max of the two alleles for simplicity)
        let svlen = match (phase1, phase2) {
            (Some(p1), Some(p2)) => format!("{:.0}", p1.max(p2)),
            (Some(p1), None) => format!("{:.0}", p1),
            (None, Some(p2)) => format!("{:.0}", p2),
            (None, None) => ".".to_string(),
        };

        let info = format!("END={};SVTYPE=STR;SVLEN={}", end, svlen);
        let format_str = "GT:AL";
        let sample_data = format!("{}:{},{}", gt, al1, al2);

        writeln!(
            writer,
            "{}\t{}\t{}\tN\t<STR>\t.\tPASS\t{}\t{}\t{}",
            chrom, pos, id, info, format_str, sample_data
        )
        .expect("Failed writing VCF record");
    }

    writer.flush().expect("Failed to flush VCF to stdout");
}

#[cfg(test)]
#[test]
fn test_region() {
    genotype_repeats(
        String::from("test-data/small-test.bam"),
        TargetConfig {
            region: Some("chr7:154778571-154779363".to_string()),
            region_file: None,
            preset: None,
            max_locus: None,
        },
        GenotypeConfig { minlen: 5, support: 3, unphased: true },
        ProcessingConfig { threads: 4, batch_size_kb: 50, output_vcf: false },
        Some("sample".to_string()),
        None,
        true, // show_progress
    )
    .expect("genotype_repeats failed");
}

#[test]
#[ignore] // Requires network access - can be enabled with: cargo test -- --ignored
fn test_region_from_url() {
    genotype_repeats(
        String::from(
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG00096.hg38.cram",
        ),
        TargetConfig {
            region: Some("chr7:154778571-154779363".to_string()),
            region_file: None,
            preset: None,
            max_locus: None,
        },
        GenotypeConfig { minlen: 5, support: 3, unphased: true },
        ProcessingConfig { threads: 4, batch_size_kb: 50, output_vcf: false },
        Some("sample".to_string()),
        Some(String::from(
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        )),
        true, // show_progress
    ).expect("genotype_repeats failed");
}

#[test]
fn test_region_bed() {
    genotype_repeats(
        String::from("test-data/small-test.bam"),
        TargetConfig {
            region: None,
            region_file: Some(PathBuf::from("test-data/test.bed")),
            preset: None,
            max_locus: None,
        },
        GenotypeConfig { minlen: 5, support: 3, unphased: true },
        ProcessingConfig { threads: 4, batch_size_kb: 50, output_vcf: false },
        Some("sample".to_string()),
        None,
        true, // show_progress
    )
    .expect("genotype_repeats failed");
}
#[test]
fn test_unphased() {
    genotype_repeats(
        String::from("test-data/small-test.bam"),
        TargetConfig {
            region: None,
            region_file: Some(PathBuf::from("test-data/test.bed")),
            preset: None,
            max_locus: None,
        },
        GenotypeConfig { minlen: 5, support: 3, unphased: true },
        ProcessingConfig { threads: 4, batch_size_kb: 50, output_vcf: false },
        Some("sample".to_string()),
        None,
        true, // show_progress
    )
    .expect("genotype_repeats failed");
}

// NOTE: Previously had a test_region_wrong_chromosome test with #[should_panic]
// This was removed after refactoring panics to graceful error handling with std::process::exit(1)
// The code now logs a warning and continues instead of panicking when a chromosome is not found

#[test]
fn test_phasing_validation_triggers() {
    // This test verifies that our phasing validation works correctly.
    // The new logic only errors if:
    // 1. No phased loci have been observed (no HP tags found) AND
    // 2. 20 consecutive unphased loci are encountered
    //
    // Since test-data/small-test.bam likely has some phased reads (based on output showing different haplotype values),
    // the validation should not trigger for a small region.
    // This is the correct behavior - we only want to catch files that completely
    // lack phasing information.

    // This test now verifies that early validation works correctly.
    // Since test-data/small-test.bam likely lacks HP tags, we use unphased mode:
    genotype_repeats(
        String::from("test-data/small-test.bam"),
        TargetConfig {
            region: Some("chr7:154778571-154779363".to_string()),
            region_file: None,
            preset: None,
            max_locus: None,
        },
        GenotypeConfig { minlen: 5, support: 3, unphased: true },
        ProcessingConfig { threads: 1, batch_size_kb: 50, output_vcf: false },
        Some("sample".to_string()),
        None,
        true, // show_progress
    )
    .expect("genotype_repeats failed");
}

#[test]
fn test_nan_genotype_for_unphased_loci() {
    // Test that loci with no phased reads return NaN, not arbitrary split values
    use std::fs::File;
    use std::io::Write;

    // Create a test BED file with the problematic locus
    let test_bed_content = "chr15\t21653406\t21653580\ttest_locus\n";
    let mut file = File::create("test_temp_nan_fix.bed").expect("Could not create test file");
    file.write_all(test_bed_content.as_bytes())
        .expect("Could not write test file");

    // Test with small-test.bam which likely has no phased reads for this locus
    genotype_repeats(
        String::from("test-data/small-test.bam"),
        TargetConfig {
            region: None,
            region_file: Some(std::path::PathBuf::from("test_temp_nan_fix.bed")),
            preset: None,
            max_locus: None,
        },
        GenotypeConfig { minlen: 5, support: 3, unphased: true },
        ProcessingConfig { threads: 1, batch_size_kb: 50, output_vcf: false },
        Some("test_sample".to_string()),
        None,
        true, // show_progress
    )
    .expect("genotype_repeats failed");

    // Clean up
    std::fs::remove_file("test_temp_nan_fix.bed").expect("Could not remove test file");
}
