use clap::AppSettings::DeriveDisplayOrder;
use clap::{Parser, Subcommand};
use log::info;
use std::path::PathBuf;

pub mod call;
pub mod combine;
pub mod outlier;
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
    },
    /// Test for association of repeat length by comparing two cohorts
    Association {},
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
            threads,
            unphased,
        } => call::genotype_repeats(bam, region, region_file, minlen, threads, unphased),
        Commands::Combine { calls, unphased } => {
            combine::combine(calls, unphased);
        }
        Commands::Scan {} => {
            unimplemented!();
        }
        Commands::Outlier { combined } => {
            outlier::outlier(combined);
        }
        Commands::Association {} => {
            unimplemented!();
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
