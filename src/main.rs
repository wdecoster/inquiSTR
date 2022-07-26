use clap::{Parser, Subcommand};
use log::info;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(author, version, about="Tool to genotype STRs from long reads", long_about = None)]
struct Cli {
    /// Bam file to genotype
    #[clap(subcommand)]
    command: Commands,
}

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

        /// fraction to extend the region intervals
        #[clap(short, long, value_parser, default_value_t = 0.5)]
        wobble: f64,
    },
    /// Combine lengths from multiple bams to a TSV
    Combine {},
    /// Search for regions potentially containing a polymorphic repeat
    Scan {},
    /// Find outliers from TSV
    Outlier {},
    /// Test for association of repeat length by comparing two cohorts
    Association {},
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let args = Cli::parse();
    info!("Collected arguments");
    match args.command {
        Commands::Call {
            bam,
            region,
            region_file,
            minlen,
            wobble,
        } => match genotype_repeats(bam, region, region_file, wobble, minlen) {
            Ok(out) => {
                println!("{}", out)
            }
            Err(error) => {
                return Err(error);
            }
        },
        Commands::Combine {} => {
            println!("Not implemented!")
        }
        Commands::Scan {} => {
            println!("Not implemented!")
        }
        Commands::Outlier {} => {
            println!("Not implemented!")
        }
        Commands::Association {} => {
            println!("Not implemented!")
        }
    }

    Ok(())
}

fn genotype_repeats(
    bamp: PathBuf,
    region: Option<String>,
    region_file: Option<PathBuf>,
    wobble: f64,
    minlen: u32,
) -> Result<String, Box<dyn std::error::Error>> {
    assert!(bamp.is_file(), "ERROR: No such file {}!", bamp.display());
    match (region, region_file) {
        (Some(_region), Some(_region_file)) => {
            panic!("ERROR: Specify either a region (-r) or region_file (-R), not both!")
        }
        (None, None) => {
            panic!("ERROR: Specify one of region (-r) or region_file (-R)!")
        }
        (Some(region), None) => {
            let (chrom, start, end) = process_region(region, wobble).unwrap();
            let bamf = bamp.into_os_string().into_string().unwrap();
            genotype_repeat(bamf, chrom, start, end, minlen)
        }
        (None, Some(region_file)) => panic!("Not implemented to use {}!", region_file.display()),
    }
}

fn genotype_repeat(
    bamf: String,
    chrom: String,
    start: i64,
    end: i64,
    minlen: u32,
) -> Result<String, Box<dyn std::error::Error>> {
    let mut bam = bam::IndexedReader::from_path(&bamf).expect("Error opening indexed BAM.");

    if let Some(tid) = bam.header().tid(chrom.as_bytes()) {
        bam.fetch((tid, start, end)).unwrap();
    } else {
        panic!("Chromosome {chrom} not found in the bam file")
    }
    let mut calls: HashMap<u8, Vec<i64>> = HashMap::new();
    calls.insert(1, Vec::new());
    calls.insert(2, Vec::new());
    calls.insert(0, Vec::new());
    info!("Checks passed, genotyping repeat");

    for r in bam.records() {
        let r = r.expect("Error reading BAM file in region.");
        if start < r.reference_start() || r.reference_end() < end || r.mapq() <= 10 {
            continue;
        }

        let mut reference_position = r.reference_start() + 1;
        let phase = get_phase(&r);
        let mut call: i64 = 0;
        for entry in r.cigar().iter() {
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    reference_position += *len as i64;
                }
                Cigar::Del(len) => {
                    if *len > minlen && start < reference_position && reference_position < end {
                        call -= *len as i64;
                    }
                    reference_position += *len as i64;
                }
                Cigar::SoftClip(len) => {
                    if *len > minlen && start < reference_position && reference_position < end {
                        call += *len as i64;
                    }
                }
                Cigar::Ins(len) => {
                    if *len > minlen && start < reference_position && reference_position < end {
                        call += *len as i64;
                    }
                }
                Cigar::RefSkip(len) => reference_position += *len as i64,
                _ => (),
            }
        }
        calls.get_mut(&phase).unwrap().push(call);
    }
    info!(
        "Used {}+{} reads for genotyping",
        calls[&1].len(),
        calls[&2].len()
    );
    Ok(format!(
        "{chrom}:{start}-{end}\t{:?}\t{:?}",
        median(&calls[&1]),
        median(&calls[&2])
    ))
}

/// parse a region string and extend the start and begin by a wobble space
/// defined by a wobble fraction relative to the interval length
fn process_region(
    reg: String,
    wobble: f64,
) -> Result<(String, i64, i64), Box<dyn std::error::Error>> {
    assert!(
        wobble >= 0.0,
        "Wobble argument should be a positive float, received {}",
        wobble
    );

    let chrom = reg.split(':').collect::<Vec<&str>>()[0];
    let interval = reg.split(':').collect::<Vec<&str>>()[1];
    let start: f64 = interval.split('-').collect::<Vec<&str>>()[0]
        .parse()
        .unwrap();
    let end: f64 = interval.split('-').collect::<Vec<&str>>()[1]
        .parse()
        .unwrap();
    assert!(
        end - start > 0.0,
        r#"Invalid region: begin has to be smaller than end."#
    );
    let wobble_length = ((end - start) * wobble) / 2.0;
    Ok((
        chrom.to_string(),
        (start - wobble_length) as i64,
        (end + wobble_length) as i64,
    ))
}

fn get_phase(record: &bam::Record) -> u8 {
    match record.aux(b"HP") {
        Ok(value) => {
            if let Aux::U8(v) = value {
                v
            } else {
                panic!("Unexpected type of Aux {:?}", value)
            }
        }
        Err(_e) => 0,
    }
}

fn median(array: &Vec<i64>) -> f64 {
    if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) as f64 / 2.0
    } else {
        array[(array.len() / 2)] as f64
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

#[test]
fn test_region() -> Result<(), Box<dyn std::error::Error>> {
    genotype_repeats(
        PathBuf::from("/home/wdecoster/test-data/test.bam"),
        Some("chr7:154778571-154779363".to_string()),
        None,
        0.5,
        5,
    )?;
    Ok(())
}

#[test]
#[should_panic]
fn test_wrong_bam_path() {
    match genotype_repeats(
        PathBuf::from("/home/wdecoster/wrong_path_to_test-data/test.bam"),
        Some("chr7:154778571-154779363".to_string()),
        None,
        0.5,
        5,
    ) {
        Ok(it) => it,
        Err(_) => "yiek".to_string(),
    };
}

#[test]
#[should_panic]
fn test_wrong_interval() {
    match process_region("chr7:154779363-154778571".to_string(), 0.5) {
        Ok(it) => it,
        Err(_) => ("yiek".to_string(), 12, 12),
    };
}
#[test]
#[should_panic]
fn test_negative_wobble() {
    match process_region("chr7:154778571-154779363".to_string(), -0.5) {
        Ok(it) => it,
        Err(_) => ("yiek".to_string(), 12, 12),
    };
}

#[test]
#[should_panic]
fn test_region_wrong_chromosome() {
    match genotype_repeats(
        PathBuf::from("/home/wdecoster/test-data/test.bam"),
        Some("7:154778571-154779363".to_string()),
        None,
        0.5,
        5,
    ) {
        Ok(it) => it,
        Err(_) => "yiek".to_string(),
    };
}
