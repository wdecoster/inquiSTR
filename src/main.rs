use clap::AppSettings::DeriveDisplayOrder;
use clap::{Parser, Subcommand};
use log::info;
use std::path::PathBuf;

pub mod assoc;
pub mod call;
pub mod combine;
pub mod histogram;
pub mod outlier;
pub mod query;
pub mod utils;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[structopt(global_settings=&[DeriveDisplayOrder])]
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
        #[clap(parse(from_os_str), validator=is_file)]
        bam: PathBuf,

        /// region string to genotype expansion in
        #[clap(short, long, value_parser)]
        region: Option<String>,

        /// Bed file with region(s) to genotype expansion(s) in
        #[clap(short = 'R', long, value_parser, validator=is_file)]
        region_file: Option<PathBuf>,

        /// minimal length of insertion/deletion operation
        #[clap(short, long, value_parser, default_value_t = 5)]
        minlen: u32,

        /// minimal number of supporting reads
        #[clap(short, long, value_parser, default_value_t = 3)]
        support: usize,

        /// Number of parallel threads to use
        #[clap(short, long, value_parser, default_value_t = 8)]
        threads: usize,

        /// If reads have to be considered unphased
        #[clap(short, long, value_parser)]
        unphased: bool,
    },
    /// Combine lengths from multiple bams to a TSV
    Combine {
        /// files from inquiSTR call
        // this validator gets applied to each element from the Vec separately
        #[clap(parse(from_os_str), multiple_values = true, required = true, validator=is_file)]
        calls: Vec<PathBuf>,

        /// If reads were unphased for inquiSTR call
        #[clap(short, long, value_parser)]
        unphased: bool,
    },
    /// Search for regions potentially containing a polymorphic repeat
    Scan {},
    /// Find outliers from TSV
    Outlier {
        /// combined file of calls
        #[clap(parse(from_os_str), required = true, validator=is_file)]
        combined: PathBuf,

        /// minimal length of expansion to be present in cohort
        #[clap(short, long, value_parser, default_value_t = 10)]
        minsize: u32,

        /// zscore cutoff to decide if a value is an outlier
        #[clap(short, long, value_parser, default_value_t = 3.0)]
        zscore: f32,
    },
    /// Lookup genotypes and display
    Query {
        /// combined file of calls
        #[clap(parse(from_os_str), required = true, validator=is_file)]
        combined: PathBuf,

        /// region to query
        #[clap(required = true)]
        region: String,
    },
    Histogram {
        /// combined file of calls
        #[clap(parse(from_os_str), required = true, validator=is_file)]
        combined: PathBuf,

        /// region to query
        #[clap(required = true)]
        region: String,
    },
    /// Test for association of repeat length by comparing two cohorts
    Association {
        /// combined file of calls
        #[clap(parse(from_os_str), required = true, validator=is_file)]
        combined: PathBuf,

        /// file with sample_id, phenotype and covariates
        #[clap(parse(from_os_str), required = true, validator=is_file)]
        metadata: PathBuf,

        /// missing genotypes cutoff
        #[clap(long, value_parser, default_value_t = 0.8)]
        missing_cutoff: f32,

        /// association mode
        #[clap(short, long, value_enum, value_parser, default_value_t = assoc::Mode::Max)]
        mode: assoc::Mode,

        /// test column and groups e.g. group:PAT,CON with <group> the name of the column containing <PAT> and <CON>
        #[clap(short, long, value_parser)]
        condition: String,

        /// covariates, comma separated
        #[clap(long, value_parser)]
        covariates: Option<String>,
        // p <- add_argument(p, "--outcometype", help = "Select a outcome variable type: binary or continuous", nargs = 1)
    },
}

fn is_file(pathname: &str) -> Result<(), String> {
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
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
        } => call::genotype_repeats(bam, region, region_file, minlen, support, threads, unphased),
        Commands::Combine { calls, unphased } => {
            combine::combine(calls, unphased);
        }
        Commands::Scan {} => {
            unimplemented!();
        }
        Commands::Outlier {
            combined,
            minsize,
            zscore,
        } => {
            outlier::outlier(combined, minsize, zscore);
        }
        Commands::Query { combined, region } => {
            query::query(combined, region);
        }
        Commands::Histogram { combined, region } => {
            histogram::histogram(combined, region);
        }
        Commands::Association {
            combined,
            metadata,
            missing_cutoff,
            mode,
            condition,
            covariates,
        } => {
            assoc::assocation(
                combined,
                metadata,
                missing_cutoff,
                mode,
                condition,
                covariates,
            );
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
