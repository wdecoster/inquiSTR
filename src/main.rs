use clap::Parser;
use log::info;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Bam file to genotype
    #[clap(parse(from_os_str))]
    bam: PathBuf,

    /// region string to genotype expansion in
    #[clap(short, long, value_parser)]
    region: String,

    // /// name of the locus in --region
    // #[clap(short, long, value_parser)]
    // locus: String,

    // /// reference genome
    // #[clap(long, value_parser)]
    // reffas: String,
    /// minimal length of insertion/deletion operation
    #[clap(short, long, value_parser, default_value_t = 5)]
    minlen: u32,

    /// fraction to extend the region intervals
    #[clap(short, long, value_parser, default_value_t = 0.5)]
    wobble: f64,
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let args = Args::parse();
    info!("Collected arguments");
    match genotype_repeat(args.bam, args.region, args.wobble, args.minlen) {
        Ok(out) => {
            println!("{}", out)
        }
        Err(error) => {
            return Err(error);
        }
    }
    Ok(())
}

fn genotype_repeat(
    bamp: PathBuf,
    region: String,
    wobble: f64,
    minlen: u32,
) -> Result<String, Box<dyn std::error::Error>> {
    let (chrom, start, end) = process_region(region, wobble).unwrap();
    assert!(bamp.is_file());

    let bamf = bamp.into_os_string().into_string().unwrap();

    let mut bam = bam::IndexedReader::from_path(&bamf).expect("Error opening indexed BAM.");

    if let Some(tid) = bam.header().tid(chrom.as_bytes()) {
        bam.fetch((tid, start, end)).unwrap();
    } else {
        panic!("Chromosome {} not found in the bam file", chrom)
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

        // let mut read_position = 0;
        let mut reference_position = r.reference_start() + 1;
        let phase = get_phase(&r);
        let mut call: i64 = 0;
        for entry in r.cigar().iter() {
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    // read_position += *len as i64;
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
                    // read_position += *len as i64;
                }
                Cigar::Ins(len) => {
                    if *len > minlen && start < reference_position && reference_position < end {
                        call += *len as i64;
                    }
                    // read_position += *len as i64;
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
        "collected H1: {:?}, H2: {:?}",
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
    Args::command().debug_assert()
}

#[test]
fn test_region() -> Result<(), Box<dyn std::error::Error>> {
    genotype_repeat(
        PathBuf::from("/home/wdecoster/test-data/test.bam"),
        "chr7:154778571-154779363".to_string(),
        0.5,
        5,
    )?;
    Ok(())
}

#[test]
#[should_panic]
fn test_wrong_bam_path() {
    match genotype_repeat(
        PathBuf::from("/home/wdecoster/wrong_path_to_test-data/test.bam"),
        "chr7:154778571-154779363".to_string(),
        0.5,
        5,
    ) {
        Ok(it) => it,
        Err(_) => return (),
    };
}

#[test]
#[should_panic]
fn test_wrong_interval() {
    match process_region("chr7:154779363-154778571".to_string(), 0.5) {
        Ok(it) => it,
        Err(_) => return (),
    };
}
#[test]
#[should_panic]
fn test_negative_wobble() {
    match process_region("chr7:154778571-154779363".to_string(), -0.5) {
        Ok(it) => it,
        Err(_) => return (),
    };
}

#[test]
#[should_panic]
fn test_region_wrong_chromosome() {
    match genotype_repeat(
        PathBuf::from("/home/wdecoster/test-data/test.bam"),
        "7:154778571-154779363".to_string(),
        0.5,
        5,
    ) {
        Ok(it) => it,
        Err(_) => return (),
    };
}
