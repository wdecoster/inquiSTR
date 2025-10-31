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

use clap::{Parser, Subcommand};
use log::info;
use std::path::PathBuf;

pub mod assoc;
pub mod bam_utils;
pub mod batch;
pub mod benchmark;
pub mod call;
pub mod combine;
pub mod histogram;
pub mod locus_search;
pub mod metadata;
pub mod outlier;
pub mod pca;
pub mod plot;
pub mod query;
pub mod repeats;
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

        /// minimal length of insertion/deletion operation
        #[clap(short, long, value_parser, default_value_t = 5)]
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
    },
    /// Combine lengths from multiple bams to a TSV, or combine kmer frequency files from unmapped
    Combine {
        /// files from inquiSTR call or inquiSTR unmapped
        // this validator gets applied to each element from the Vec separately
        #[clap(value_parser, required = true)]
        calls: Vec<PathBuf>,

        /// Number of threads to use for parallel processing
        #[clap(short = 't', long, value_parser, default_value_t = 1)]
        threads: usize,
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
    Histogram {
        /// combined file of calls
        #[clap(value_parser, required = true)]
        combined: PathBuf,

        /// region to query
        #[clap(required = true)]
        region: String,
    },
    /// Perform association testing for STRs using embedded R script
    #[clap(arg_required_else_help = true)]
    Association {
        /// Combined STR file from inquiSTR combine command
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

        /// STR mode: MEAN, MAX, or MIN for H1/H2 combination
        #[clap(long, value_parser)]
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
    },
    /// Show a histogram with multiple groups for a specific repeat
    Plot {
        /// combined file of calls
        #[clap(value_parser, required = true)]
        combined: PathBuf,

        /// file with sample_id, phenotype and covariates
        #[clap(value_parser, required = true)]
        metadata: PathBuf,

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
    /// Perform PCA analysis on combined STR data
    Pca {
        /// Combined file of STR calls from inquiSTR combine command
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

        /// Method for aggregating H1/H2 allele lengths: max (default), min, or sum
        #[clap(short, long, value_parser, default_value_t = pca::AlleleAggregation::Max)]
        aggregation: pca::AlleleAggregation,
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
    },
    /// Clean the index cache directory (hidden command)
    #[clap(hide = true)]
    CleanCache {
        /// Show what would be deleted without actually deleting
        #[clap(long)]
        dry_run: bool,

        /// Delete all cached files regardless of age
        #[clap(long)]
        all: bool,

        /// Maximum age in days (files older than this will be deleted)
        #[clap(long, default_value_t = 30)]
        max_age_days: u64,
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
            minlen,
            support,
            threads,
            unphased,
            sample_name,
            reference,
            max_locus,
            batch_size,
        } => call::genotype_repeats(
            bam,
            region,
            region_file,
            minlen,
            support,
            threads,
            unphased,
            sample_name,
            reference,
            max_locus,
            batch_size,
        ),
        Commands::Combine { calls, threads } => {
            combine::combine(calls, threads);
        }
        Commands::Outlier { combined, minsize, zscore, method, sample, threads } => {
            if !combined.exists() {
                eprintln!("ERROR: Combined file does not exist: {}", combined.display());
                std::process::exit(1);
            }
            let subset = sample.map(|s| outlier::parse_sample_input(&s));
            outlier::outlier(combined, minsize, zscore, method, subset, threads);
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
            );
        }
        Commands::Plot { combined, metadata, condition, region, output } => {
            plot::plot(combined, metadata, condition, region, output)
        }
        Commands::Pca { combined, output, components, threads, aggregation } => {
            if !combined.exists() {
                eprintln!("ERROR: Combined file does not exist: {}", combined.display());
                std::process::exit(1);
            }
            pca::pca(combined, output, components, threads, aggregation);
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
            unmapped::count_unmapped_kmers(
                bam,
                klength,
                sample_name,
                reference,
                threads,
                target_kmer,
                combine_revcomp,
            );
        }
        Commands::Benchmark { inquistr, vcf, bed, mode, plot, max_plot_length, tier1, diff_out } => {
            benchmark::benchmark(inquistr, vcf, bed, mode, plot, max_plot_length, tier1, diff_out);
        }
        Commands::CleanCache { dry_run, all, max_age_days } => {
            bam_utils::clean_cache(dry_run, all, max_age_days);
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
