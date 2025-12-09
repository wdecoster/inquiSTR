//! # inquiSTR
//!
//! A high-performance toolkit for genotyping and analyzing Short Tandem Repeats (STRs)
//! from long-read sequencing data. Optimized for Oxford Nanopore Technologies (ONT) data
//! with support for phased and unphased analysis.
//!
//! ## Features
//!
//! - **Fast STR genotyping** with memory-efficient batch processing
//! - **Phased analysis** using HP tags from phased BAM files  
//! - **Multi-sample analysis** with outlier detection
//! - **Multiple output formats** including TSV and plots
//! - **Quality control** with comprehensive validation
//!
//! ## Performance
//!
//! inquiSTR uses several optimizations for processing millions of STR targets:
//! - Region batching to minimize I/O operations
//! - Memory-efficient read processing
//! - Parallel processing for multi-sample workflows
//! - Optimized data structures for STR call storage

#![allow(clippy::uninlined_format_args)]

use clap::{Args, Parser, Subcommand};
use log::info;
use std::path::PathBuf;

/// inquiSTR version from Cargo.toml
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Common arguments for batch processing
#[derive(Args, Debug)]
struct BatchCommonArgs {
    /// TSV manifest file with bam_path (required) and sample_name (optional) columns
    #[clap(value_parser, required = true)]
    manifest: PathBuf,

    /// Output file for combined results
    #[clap(short, long, value_parser, required = true)]
    output: PathBuf,

    /// Number of parallel threads to use. By default, this parallelizes across samples
    /// (processes N samples with 1 thread each - optimal for CRAM files).
    /// Use --parallel-samples to override this behavior.
    #[clap(short, long, value_parser, default_value_t = 1)]
    threads: usize,

    /// Number of samples to process in parallel (threads will be divided among them).
    /// If not specified, defaults to --threads value (optimal for CRAM: N samples with 1 thread each).
    /// Set explicitly for advanced use cases (e.g., --threads 16 --parallel-samples 4 = 4 threads per sample).
    #[clap(long, value_parser)]
    parallel_samples: Option<usize>,

    /// Reference fasta for CRAM decoding (applies to all samples)
    #[clap(long, value_parser)]
    reference: Option<String>,

    /// Save individual sample files to this directory (optional)
    #[clap(long, value_parser)]
    save_individual: Option<PathBuf>,

    /// Temporary directory for intermediate files (default: $TMPDIR or current directory)
    #[clap(long, value_parser)]
    tmpdir: Option<PathBuf>,

    /// Resume: skip already processed samples by checking the output file
    #[clap(long, value_parser)]
    resume: bool,

    /// Show what would be processed without actually running
    #[clap(long, value_parser)]
    dry_run: bool,

    /// Continue processing even if some samples fail (by default, exits with error if any sample fails)
    #[clap(long, value_parser)]
    keep_going: bool,

    /// Enable performance profiling and recommendations (minimal overhead)
    #[clap(long, value_parser)]
    profile: bool,

    /// Skip file path validation before processing (faster startup, but errors will occur during processing if paths are invalid)
    #[clap(long, value_parser)]
    skip_validation: bool,
}

/// Mode selection for batch processing
#[derive(Args, Debug)]
#[group(id = "mode")]
#[command(next_help_heading = "Mode Selection")]
struct BatchModeArgs {
    /// Process unmapped reads instead of genotyping STRs
    #[clap(long, value_parser)]
    unmapped: bool,
}

/// STR genotyping mode arguments
#[derive(Args, Debug)]
#[group(id = "str_mode")]
#[command(next_help_heading = "STR Genotyping Options")]
struct BatchStrArgs {
    /// Region string to genotype expansion in
    #[clap(short, long, value_parser)]
    region: Option<String>,

    /// Bed file with region(s) to genotype expansion(s) in
    #[clap(short = 'R', long, value_parser)]
    region_file: Option<PathBuf>,

    /// Use a predefined TR catalog (pathogenic, adotto, trexplorer, or codis)
    #[clap(long, value_parser)]
    preset: Option<repeats::TRPreset>,

    /// Minimal length of insertion/deletion operation
    #[clap(short, long, value_parser, default_value_t = 1)]
    minlen: u32,

    /// Minimal number of supporting reads
    #[clap(short, long, value_parser, default_value_t = 3)]
    support: usize,

    /// If reads have to be considered unphased
    #[clap(short, long, value_parser)]
    unphased: bool,

    /// Maximum locus size to consider (intervals larger than this will be filtered out)
    #[clap(long, value_parser)]
    max_locus: Option<u32>,

    /// Batch size in KB for grouping nearby STR targets
    #[clap(long, value_parser, default_value_t = 50)]
    batch_size: u32,
}

/// Unmapped kmer mode arguments
#[derive(Args, Debug)]
#[group(id = "unmapped_mode")]
#[command(next_help_heading = "Unmapped Kmer Options")]
struct BatchUnmappedArgs {
    /// Maximum kmer length to count (counts all sizes from 2 to klength)
    #[clap(short = 'k', long, value_parser, default_value_t = 6)]
    klength: usize,

    /// Target kmer to specifically quantify (optional, can be any length)
    #[clap(long, value_parser)]
    target_kmer: Option<String>,

    /// Combine kmers with their reverse complements
    #[clap(long, value_parser, default_value_t = false)]
    combine_revcomp: bool,
}

pub mod assoc;
pub mod bam_utils;
pub mod batch_process;
pub mod batching;
pub mod benchmark;
pub mod call;
pub mod combine;
pub mod convert;
pub mod errors;
pub mod filetype;
pub mod filter;
pub mod genotype_batch;
pub mod histogram;
pub mod locus_search;
pub mod outlier;
pub mod pca;
pub mod plot;
pub mod query;
pub mod relate;
pub mod repeats;
pub mod sample_info;
pub mod unmapped;
pub mod utils;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[clap(author, version, about="Tool to genotype STRs from long reads", long_about = None)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}
// Every subcommand is a variation of the Commands Enum, and has its arguments defined below
#[derive(Debug, Subcommand)]
enum Commands {
    /// Call lengths
    #[clap(arg_required_else_help = true)]
    Call {
        /// bam file to call STRs in
        #[clap()]
        bam: String,

        /// region string to genotype expansion in
        #[clap(short, long, value_parser)]
        region: Option<String>,

        /// Bed file with region(s) to genotype expansion(s) in
        #[clap(short = 'R', long, value_parser)]
        region_file: Option<PathBuf>,

        /// Use a predefined TR catalog (pathogenic, adotto, trexplorer, or codis)
        #[clap(long, value_parser)]
        preset: Option<repeats::TRPreset>,

        /// minimal length of insertion/deletion operation
        #[clap(short, long, value_parser, default_value_t = 1)]
        minlen: u32,

        /// minimal number of supporting reads
        #[clap(short, long, value_parser, default_value_t = 3)]
        support: usize,

        /// Number of parallel threads to use
        #[clap(short, long, value_parser, default_value_t = 1)]
        threads: usize,

        /// If reads have to be considered unphased
        #[clap(short, long, value_parser)]
        unphased: bool,

        /// sample name to use in output
        #[clap(long, value_parser)]
        sample_name: Option<String>,

        /// reference fasta for cram decoding
        #[clap(long, value_parser)]
        reference: Option<String>,

        /// maximum locus size to consider (intervals larger than this will be filtered out)
        #[clap(long, value_parser)]
        max_locus: Option<u32>,

        /// Batch size in KB for grouping nearby STR targets (default: 50). Larger values use more memory but reduce I/O operations. Use 20-35 for memory-constrained systems, 80-100 for high-performance setups.
        #[clap(long, value_parser, default_value_t = 50)]
        batch_size: u32,

        /// Output VCF file path (optional, TSV still written to stdout)
        #[clap(long, value_parser)]
        vcf: Option<PathBuf>,
    },
    /// Process multiple samples in batch and combine results
    #[clap(arg_required_else_help = true)]
    Batch {
        #[clap(flatten)]
        common: BatchCommonArgs,

        #[clap(flatten)]
        mode: BatchModeArgs,

        #[clap(flatten)]
        str_args: BatchStrArgs,

        #[clap(flatten)]
        unmapped_args: BatchUnmappedArgs,
    },
    /// Combine STR lengths or kmer frequency files from multiple samples
    Combine {
        /// files from inquiSTR call or inquiSTR unmapped
        // this validator gets applied to each element from the Vec separately
        #[clap(value_parser, required = true)]
        calls: Vec<PathBuf>,

        /// Number of threads to use for parallel processing
        #[clap(short = 't', long, value_parser, default_value_t = 1)]
        threads: usize,
    },
    /// Filter inquiSTR output by various criteria
    #[clap(arg_required_else_help = true)]
    Filter {
        /// Input file from inquiSTR call or inquiSTR combine
        #[clap(value_parser, required = true)]
        input: PathBuf,

        /// Minimum allele length (expansion/insertion only, ignores negative values)
        #[clap(long, value_parser)]
        minlen: Option<u32>,

        /// Minimum absolute allele length change (uses absolute value)
        #[clap(long, value_parser)]
        minchange: Option<u32>,

        /// BED file to filter by overlap (requires both files to be sorted)
        #[clap(long, value_parser)]
        bed: Option<PathBuf>,

        /// Minimum call rate (fraction 0.0-1.0) for combined files
        #[clap(long, value_parser)]
        call_rate: Option<f64>,

        /// Minimum coefficient of variation (only for combined files)
        #[clap(long, value_parser)]
        min_cv: Option<f64>,
    },
    /// Find outliers from TSV
    Outlier {
        /// combined file of calls
        #[clap(value_parser, required = true)]
        combined: PathBuf,

        /// minimal length of expansion to be present in cohort
        #[clap(long, value_parser, default_value_t = 10)]
        minsize: u32,

        /// zscore cutoff to decide if a value is an outlier
        #[clap(short, long, value_parser, default_value_t = 3.0)]
        zscore: f32,

        /// method to test for outliers
        #[clap(long, value_enum, value_parser, default_value_t = outlier::Method::Zscore)]
        method: outlier::Method,

        /// sample(s) to consider: can be a single sample name, comma-separated sample names, or a file path containing sample names (one per line)
        #[clap(short = 's', long, value_parser)]
        sample: Option<String>,

        /// Number of threads to use for parallel processing
        #[clap(short = 't', long, value_parser, default_value_t = 1)]
        threads: usize,

        /// Output file to write outlier counts per sample (tab-separated: sample\tcount)
        #[clap(short = 'c', long, value_parser)]
        count: Option<PathBuf>,
    },
    /// Lookup genotypes and display
    Query {
        /// combined file of calls
        #[clap(value_parser, required = true)]
        combined: PathBuf,

        /// region to query or file with regions to query
        #[clap(required = true)]
        region: String,
    },
    /// Generate text-based histogram for a specific repeat
    Histogram {
        /// combined file of calls
        #[clap(value_parser, required = true)]
        combined: PathBuf,

        /// region to query
        #[clap(required = true)]
        region: String,
    },
    /// Perform association testing for STRs or kmer frequencies using embedded R script
    #[clap(arg_required_else_help = true)]
    Association {
        /// Combined STR or kmer frequency file from inquiSTR combine command
        #[clap(short, long, value_parser)]
        input: PathBuf,

        /// Phenotype and covariate file with header, first column is individual ID
        #[clap(short, long, value_parser)]
        phenocovar: PathBuf,

        /// Column name of phenotype in phenocovar file
        #[clap(long, value_parser)]
        phenotype: String,

        /// Output file name for association results
        #[clap(short, long, value_parser)]
        out: PathBuf,

        /// STR mode: MEAN, MAX, or MIN for H1/H2 combination (not used for kmer data)
        #[clap(long, value_parser, default_value = "MAX")]
        str_mode: String,

        /// Outcome type: binary or continuous
        #[clap(long, value_parser)]
        outcometype: String,

        /// Covariate names, comma separated (optional)
        #[clap(long, value_parser)]
        covnames: Option<String>,

        /// Call rate cutoff for variants (default 0.80)
        #[clap(long, value_parser, default_value_t = 0.80)]
        missing_cutoff: f64,

        /// Minimum maximum STR length across samples for variant to be included
        #[clap(long, value_parser)]
        minimal_length: Option<f64>,

        /// Number of threads for parallel processing
        #[clap(short, long, value_parser, default_value_t = 1)]
        threads: usize,

        /// Number of variants to process in each chunk (default 1000)
        #[clap(long, value_parser, default_value_t = 1000)]
        chunk_size: usize,

        /// Binary phenotype order, comma separated (e.g., Control,Patient) - required for binary outcomes
        #[clap(long, value_parser)]
        binary_order: Option<String>,

        /// Do not print progress messages
        #[clap(long)]
        quiet: bool,

        /// Generate QQ plot and Manhattan plot with custom prefix (e.g., --plot myresults). If no prefix given, uses output filename stem.
        #[clap(long, value_parser)]
        plot: Option<String>,
    },
    /// Show a histogram with multiple groups for a specific repeat
    Plot {
        /// combined file of calls
        #[clap(value_parser, required = true)]
        combined: PathBuf,

        /// file with sample IDs, phenotypes, and covariates (TSV format)
        #[clap(value_parser, required = true)]
        sample_metadata: PathBuf,

        /// test column and groups to plot e.g. group:PAT,CON with <group> the name of the column containing <PAT> and <CON>
        #[clap(short, long, value_parser)]
        condition: String,

        /// region to query
        #[clap(required = true)]
        region: String,

        /// HTML output file name
        #[clap(short, long, value_parser, default_value_t=String::from("groupplot.html"))]
        output: String,
    },
    /// Perform PCA analysis on combined STR or kmer data
    Pca {
        /// Combined file from inquiSTR combine command (supports both STR calls and kmer frequencies)
        #[clap(value_parser, required = true)]
        combined: PathBuf,

        /// HTML output file name for interactive PCA plot
        #[clap(short, long, value_parser, default_value_t=String::from("pca_plot.html"))]
        output: String,

        /// Number of principal components to compute (currently only first 2 are plotted)
        #[clap(short, long, value_parser, default_value_t = 10)]
        components: usize,

        /// Number of threads to use for parallel processing
        #[clap(short, long, value_parser, default_value_t = 1)]
        threads: usize,

        /// Method for aggregating H1/H2 allele lengths: max (default), min, or sum (only applies to STR files, ignored for kmer files)
        #[clap(short, long, value_parser, default_value_t = pca::AlleleAggregation::Max)]
        aggregation: pca::AlleleAggregation,

        /// Output file for PC scores (tab-separated: sample, PC1, PC2, ...) to use as covariates
        #[clap(long, value_parser)]
        scores: Option<PathBuf>,
    },
    /// Count kmer frequencies in unmapped reads
    #[clap(arg_required_else_help = true)]
    Unmapped {
        /// BAM/CRAM file to analyze unmapped reads from
        #[clap()]
        bam: String,

        /// Maximum kmer length to count (counts all sizes from 2 to klength)
        #[clap(short = 'k', long, value_parser, default_value_t = 6)]
        klength: usize,

        /// Sample name to use in output header
        #[clap(long, value_parser)]
        sample_name: Option<String>,

        /// Reference fasta for CRAM decoding
        #[clap(long, value_parser)]
        reference: Option<String>,

        /// Number of parallel threads to use
        #[clap(short, long, value_parser, default_value_t = 1)]
        threads: usize,

        /// Target kmer to specifically quantify (optional, can be any length). Supports shorthand notation: (CT)4 = CTCTCTCT
        #[clap(long, value_parser)]
        target_kmer: Option<String>,

        /// Combine kmers with their reverse complements (e.g., CTCTCT and AGAGAG counted together)
        #[clap(long, value_parser, default_value_t = false)]
        combine_revcomp: bool,
    },
    /// Benchmark inquiSTR calls against a truth VCF or BED file
    Benchmark {
        /// inquiSTR call output file (.inq format)
        #[clap(value_parser, required = true)]
        inquistr: PathBuf,

        /// VCF file with truth genotypes (can be compressed)
        #[clap(long, value_parser)]
        vcf: Option<PathBuf>,

        /// BED file with truth genotypes (can be compressed, 9 columns with last 2 being haplotype lengths)
        #[clap(long, value_parser)]
        bed: Option<PathBuf>,

        /// Mode for selecting alleles: MAX (default) or MIN
        #[clap(short, long, value_parser, default_value_t = String::from("MAX"))]
        mode: String,

        /// Output file for correlation plot
        #[clap(short, long, value_parser, required = true)]
        plot: PathBuf,

        /// Maximum length to display on plot (points beyond this are counted but not shown)
        #[clap(long, value_parser, default_value_t = 5000.0)]
        max_plot_length: f64,

        /// Only use Tier1 variants from BED file (ignores Tier2 variants)
        #[clap(long)]
        tier1: bool,

        /// Create an output file for largest differences between truth and calls
        #[clap(short, long, value_parser)]
        diff_out: Option<PathBuf>,

        /// Maximum locus size in bp to include in benchmark (filters out large intervals from truth data)
        #[clap(long, value_parser)]
        max_locus: Option<u32>,

        /// Exclude zero-zero pairs (unchanged alleles) from correlation and plot
        #[clap(long)]
        nonzero: bool,

        /// Tolerance in bp for considering calls as matching (default: 5)
        #[clap(long, value_parser, default_value_t = 5)]
        tolerance: u32,
    },
    /// Compute relatedness between samples
    #[clap(arg_required_else_help = true)]
    Relate {
        /// Combined file of STR calls from inquiSTR combine command
        #[clap(value_parser, required = true)]
        combined: PathBuf,

        /// Output file for relatedness matrix
        #[clap(short, long, value_parser, required = true)]
        output: PathBuf,

        /// Number of threads to use for parallel processing
        #[clap(short, long, value_parser, default_value_t = 1)]
        threads: usize,
    },
    /// Convert VCF files to inquiSTR format
    #[clap(arg_required_else_help = true)]
    Convert {
        /// VCF file(s) to convert (can be compressed). Single sample VCF produces individual_call file, multiple samples or multiple VCFs produce combined_call file.
        #[clap(value_parser, required = true)]
        vcf: Vec<PathBuf>,
    },
}

fn main() {
    env_logger::init();
    let args = Cli::parse();
    info!("Collected arguments");
    match args.command {
        Commands::Call {
            bam,
            region,
            region_file,
            preset,
            minlen,
            support,
            threads,
            unphased,
            sample_name,
            reference,
            max_locus,
            batch_size,
            vcf,
        } => {
            if let Err(e) = call::genotype_repeats(
                bam,
                repeats::TargetConfig { region, region_file, preset, max_locus },
                call::GenotypeConfig { minlen, support, unphased },
                call::ProcessingConfig { threads, batch_size_kb: batch_size, output_vcf: vcf },
                sample_name,
                reference,
                true, // show_progress
            ) {
                eprintln!("ERROR: {}", e.message);
                std::process::exit(1);
            }
        }
        Commands::Batch { common, mode, str_args, unmapped_args } => {
            // Smart default: if parallel_samples not specified, set it equal to threads
            // This gives optimal CRAM performance (N samples in parallel, 1 thread each)
            let parallel_samples = common.parallel_samples.unwrap_or(common.threads);

            let batch_config = batch_process::BatchConfig {
                manifest: common.manifest,
                output: common.output,
                save_individual: common.save_individual,
                tmpdir: common.tmpdir,
                resume: common.resume,
                dry_run: common.dry_run,
                keep_going: common.keep_going,
                reference: common.reference,
                parallel_samples,
                profile: common.profile,
                skip_validation: common.skip_validation,
            };

            let batch_mode = if mode.unmapped {
                batch_process::BatchMode::UnmappedKmer {
                    unmapped_config: batch_process::UnmappedConfig {
                        klength: unmapped_args.klength,
                        target_kmer: unmapped_args.target_kmer,
                        combine_revcomp: unmapped_args.combine_revcomp,
                    },
                    processing_config: call::ProcessingConfig {
                        threads: common.threads,
                        batch_size_kb: str_args.batch_size,
                        output_vcf: None,
                    },
                }
            } else {
                batch_process::BatchMode::StrGenotyping {
                    target_config: repeats::TargetConfig {
                        region: str_args.region,
                        region_file: str_args.region_file,
                        preset: str_args.preset,
                        max_locus: str_args.max_locus,
                    },
                    genotype_config: call::GenotypeConfig {
                        minlen: str_args.minlen,
                        support: str_args.support,
                        unphased: str_args.unphased,
                    },
                    processing_config: call::ProcessingConfig {
                        threads: common.threads,
                        batch_size_kb: str_args.batch_size,
                        output_vcf: None,
                    },
                }
            };

            batch_process::batch_process(batch_config, batch_mode);
        }
        Commands::Combine { calls, threads } => {
            combine::combine(calls, threads);
        }
        Commands::Filter { input, minlen, minchange, bed, call_rate, min_cv } => {
            filter::filter(input, minlen, minchange, bed, call_rate, min_cv);
        }
        Commands::Outlier { combined, minsize, zscore, method, sample, threads, count } => {
            if !combined.exists() {
                eprintln!("ERROR: Combined file does not exist: {}", combined.display());
                std::process::exit(1);
            }
            let subset = sample.map(|s| outlier::parse_sample_input(&s));
            outlier::outlier(combined, minsize, zscore, method, subset, threads, count);
        }
        Commands::Query { combined, region } => {
            query::query(combined, region);
        }
        Commands::Histogram { combined, region } => {
            histogram::histogram(combined, region);
        }
        Commands::Association {
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
            plot,
        } => {
            assoc::run_association(
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
                plot,
            );
        }
        Commands::Plot { combined, sample_metadata, condition, region, output } => {
            plot::plot(combined, sample_metadata, condition, region, output)
        }
        Commands::Pca { combined, output, components, threads, aggregation, scores } => {
            if !combined.exists() {
                eprintln!("ERROR: Combined file does not exist: {}", combined.display());
                std::process::exit(1);
            }
            pca::pca(combined, output, components, threads, aggregation, scores);
        }
        Commands::Unmapped {
            bam,
            klength,
            sample_name,
            reference,
            threads,
            target_kmer,
            combine_revcomp,
        } => {
            if let Err(e) = unmapped::count_unmapped_kmers(
                bam,
                klength,
                sample_name,
                reference,
                threads,
                target_kmer,
                combine_revcomp,
                true, // show_progress
            ) {
                eprintln!("ERROR: {}", e.message);
                std::process::exit(1);
            }
        }
        Commands::Benchmark {
            inquistr,
            vcf,
            bed,
            mode,
            plot,
            max_plot_length,
            tier1,
            diff_out,
            max_locus,
            nonzero,
            tolerance,
        } => {
            benchmark::benchmark(
                inquistr,
                vcf,
                bed,
                mode,
                plot,
                max_plot_length,
                tier1,
                diff_out,
                max_locus,
                nonzero,
                tolerance,
            );
        }
        Commands::Relate { combined, output, threads } => {
            relate::compute_relatedness(combined, output, threads);
        }
        Commands::Convert { vcf } => {
            if let Err(e) = convert::convert_vcf(vcf) {
                eprintln!("ERROR: {}", e.message);
                std::process::exit(1);
            }
        }
    }
}

#[cfg(test)]
#[ctor::ctor]
fn init() {
    env_logger::init();
}

#[test]
fn verify_app() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
