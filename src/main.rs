use clap::AppSettings::DeriveDisplayOrder;
use clap::{Parser, Subcommand};
use log::info;
use std::path::PathBuf;

pub mod call;
pub mod combine;

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
        #[clap(parse(from_os_str))]
        bam: PathBuf,

        /// region string to genotype expansion in
        #[clap(short, long, value_parser)]
        region: Option<String>,

        /// Bed file with region(s) to genotype expansion(s) in
        #[clap(short = 'R', long, value_parser)]
        region_file: Option<PathBuf>,

        /// minimal length of insertion/deletion operation
        #[clap(short, long, value_parser, default_value_t = 5)]
        minlen: u32,

        /// Number of parallel threads to use
        #[clap(short, long, value_parser, default_value_t = 8)]
        threads: usize,
    },
    /// Combine lengths from multiple bams to a TSV
    Combine {
        /// bam file to call STRs in
        #[clap(parse(from_os_str), multiple_values = true)]
        calls: Vec<PathBuf>,
    },
    /// Search for regions potentially containing a polymorphic repeat
    Scan {},
    /// Find outliers from TSV
    Outlier {},
    /// Test for association of repeat length by comparing two cohorts
    Association {},
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
        } => call::genotype_repeats(bam, region, region_file, minlen, threads),
        Commands::Combine { calls } => {
            combine::combine(calls);
        }
        Commands::Scan {} => {
            unimplemented!();
        }
        Commands::Outlier {} => {
            unimplemented!();
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
