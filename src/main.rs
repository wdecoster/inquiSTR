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
use std::{io::BufRead, path::PathBuf};

pub mod bam_utils;
pub mod batch;
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
    /// Combine lengths from multiple bams to a TSV
    Combine {
        /// files from inquiSTR call
        // this validator gets applied to each element from the Vec separately
        #[clap(value_parser, required = true)]
        calls: Vec<PathBuf>,
    },
    /// Search for regions potentially containing a polymorphic repeat
    Scan {},
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

        /// sample to consider
        #[clap(short = 's', long, value_parser)]
        sample: Option<String>,

        /// file with subset of samples to consider
        #[clap(short = 'S', long, value_parser)]
        subset: Option<PathBuf>,
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
    // /// Test for association of repeat length by comparing two cohorts
    // Association {
    //     /// combined file of calls
    //     #[clap(parse(from_os_str), required = true, validator=is_file)]
    //     combined: PathBuf,

    //     /// file with sample_id, phenotype and covariates
    //     #[clap(parse(from_os_str), required = true, validator=is_file)]
    //     metadata: PathBuf,

    //     /// missing genotypes cutoff
    //     #[clap(long, value_parser, default_value_t = 0.8)]
    //     missing_cutoff: f32,

    //     /// association mode
    //     #[clap(short, long, value_enum, value_parser, default_value_t = assoc::Mode::Max)]
    //     mode: assoc::Mode,

    //     /// test column and 2 groups to test e.g. group:PAT,CON with <group> the name of the column containing <PAT> and <CON>
    //     #[clap(short, long, value_parser)]
    //     condition: String,

    //     /// covariates, comma separated
    //     #[clap(long, value_parser)]
    //     covariates: Option<String>,
    //     // p <- add_argument(p, "--outcometype", help = "Select a outcome variable type: binary or continuous", nargs = 1)
    // },
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

        /// Number of threads to use for parallel processing (0 = auto-detect)
        #[clap(short, long, value_parser, default_value_t = 0)]
        threads: usize,
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
        Commands::Combine { calls } => {
            combine::combine(calls);
        }
        Commands::Scan {} => {
            unimplemented!();
        }
        Commands::Outlier { combined, minsize, zscore, method, sample, subset } => {
            if !combined.exists() {
                panic!("Combined file does not exist!");
            }
            let subset = match (sample, subset) {
                (Some(_), Some(_)) => {
                    panic!("Cannot use both -s and -S arguments");
                }
                (Some(sample), None) => Some(vec![sample]),
                (None, Some(subset)) => {
                    let file =
                        crate::utils::reader(&subset.into_os_string().into_string().unwrap());
                    Some(file.lines().map(|line| line.unwrap()).collect())
                }
                (None, None) => None,
            };
            outlier::outlier(combined, minsize, zscore, method, subset);
        }
        Commands::Query { combined, region } => {
            query::query(combined, region);
        }
        Commands::Histogram { combined, region } => {
            histogram::histogram(combined, region);
        }
        // Commands::Association {
        //     combined,
        //     metadata,
        //     missing_cutoff,
        //     mode,
        //     condition,
        //     covariates,
        // } => {
        //     assoc::assocation(
        //         combined,
        //         metadata,
        //         missing_cutoff,
        //         mode,
        //         condition,
        //         covariates,
        //     );
        // }
        Commands::Plot { combined, metadata, condition, region, output } => {
            plot::plot(combined, metadata, condition, region, output)
        }
        Commands::Pca { combined, output, components, threads } => {
            if !combined.exists() {
                panic!("Combined file does not exist!");
            }
            pca::pca(combined, output, components, threads);
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
