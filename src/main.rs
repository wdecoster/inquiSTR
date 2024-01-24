use clap::{Parser, Subcommand};
use log::info;
use std::path::PathBuf;

// pub mod assoc;
pub mod call;
pub mod combine;
pub mod histogram;
pub mod metadata;
pub mod outlier;
pub mod plot;
pub mod query;
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
        #[clap(short, long, value_parser, default_value_t = 8)]
        threads: usize,

        /// If reads have to be considered unphased
        #[clap(short, long, value_parser)]
        unphased: bool,

        /// sample name to use in output
        #[clap(long, value_parser)]
        sample_name: Option<String>,
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
        } => call::genotype_repeats(
            bam,
            region,
            region_file,
            minlen,
            support,
            threads,
            unphased,
            sample_name,
        ),
        Commands::Combine { calls } => {
            combine::combine(calls);
        }
        Commands::Scan {} => {
            unimplemented!();
        }
        Commands::Outlier {
            combined,
            minsize,
            zscore,
            method,
        } => {
            outlier::outlier(combined, minsize, zscore, method);
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
        Commands::Plot {
            combined,
            metadata,
            condition,
            region,
            output,
        } => plot::plot(combined, metadata, condition, region, output),
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
