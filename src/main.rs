use bio::io::bed;
use clap::{Parser, Subcommand};
use log::{error, info};
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::f64::NAN;
use std::path::PathBuf;
use std::sync::Mutex;

// #[derive(Eq)]
// struct Genotype {
//     chrom: String,
//     start: u64,
//     end: u64,
//     phase1: f64,
//     phase2: f64
// }

// impl Ord for Genotype {
//     fn cmp(&self, other: &Self) -> Ordering {
//         self.height.cmp(&other.height)
//     }
// }

// impl PartialOrd for Person {
//     fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
//         Some(self.cmp(other))
//     }
// }

// impl PartialEq for Person {
//     fn eq(&self, other: &Self) -> bool {
//         self.height == other.height
//     }
// }

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

        /// Number of parallel threads to use
        #[clap(short, long, value_parser, default_value_t = 8)]
        threads: usize,
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
            wobble,
            threads,
        } => genotype_repeats(bam, region, region_file, wobble, minlen, threads),
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
}

fn genotype_repeats(
    bamp: PathBuf,
    region: Option<String>,
    region_file: Option<PathBuf>,
    wobble: f64,
    minlen: u32,
    threads: usize,
) {
    if !bamp.is_file() {
        error!(
            "ERROR: path to bam file {} is not valid!\n\n",
            bamp.display()
        );
        panic!();
    };
    match (region, region_file) {
        (Some(_region), Some(_region_file)) => {
            error!("ERROR: Specify either a region (-r) or region_file (-R), not both!\n\n");
            panic!();
        }
        (None, None) => {
            error!("ERROR: Specify one of region (-r) or region_file (-R)!\n\n");
            panic!();
        }
        (Some(region), None) => {
            let (chrom, start, end) = process_region(region, wobble).unwrap();
            let bamf = bamp.into_os_string().into_string().unwrap();
            match genotype_repeat(&bamf, chrom, start, end, minlen) {
                Ok(output) => write_genotype(output),
                Err(chrom) => error!("Contig {chrom} not found in bam file"),
            };
        }
        (None, Some(region_file)) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();
            let mut reader =
                bed::Reader::from_file(region_file.into_os_string().into_string().unwrap())
                    .unwrap();
            let bamf = bamp.into_os_string().into_string().unwrap();
            let chrom_reported = Mutex::new(Vec::new());
            let genotypes = Mutex::new(Vec::new());
            reader.records().par_bridge().for_each(|record| {
                let rec = record.expect("Error reading bed record.");
                match genotype_repeat(
                    &bamf,
                    rec.chrom().to_string(),
                    rec.start().try_into().unwrap(),
                    rec.end().try_into().unwrap(),
                    minlen,
                ) {
                    Ok(output) => {
                        let mut geno = genotypes.lock().unwrap();
                        geno.push(output);
                    }
                    Err(chrom) => {
                        let mut chroms_reported = chrom_reported.lock().unwrap();
                        if !chroms_reported.contains(&chrom) {
                            error!("Contig {chrom} not found in bam file");
                            chroms_reported.push(chrom);
                        }
                    }
                };
            });
            let mut genotypes = genotypes.lock().unwrap();
            genotypes.sort_unstable_by(|chrom, start, end, _, _| (k.0, k.1, k.2));
        }
    }
}

fn genotype_repeat(
    bamf: &String,
    chrom: String,
    start: i64,
    end: i64,
    minlen: u32,
) -> Result<(String, i64, i64, f64, f64), String> {
    let mut bam = match bam::IndexedReader::from_path(&bamf) {
        Ok(handle) => handle,
        Err(e) => {
            error!("Error opening BAM {}.\n{}", bamf, e);
            panic!();
        }
    };

    if let Some(tid) = bam.header().tid(chrom.as_bytes()) {
        bam.fetch((tid, start, end)).unwrap();

        let mut calls: HashMap<u8, Vec<i64>> =
            HashMap::from([(1, Vec::new()), (2, Vec::new()), (0, Vec::new())]);

        info!("Checks passed, genotyping repeat");

        for r in bam.records() {
            let r = r.expect("Error reading BAM file in region {chrom}:{start}-{end}.");
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
        Ok((chrom, start, end, median(&calls[&1]), median(&calls[&2])))
    } else {
        Err(chrom)
    }
}

fn write_genotype(genotype: (String, i64, i64, f64, f64)) {
    let (chrom, start, end, phase1, phase2) = genotype;
    println!("{chrom}\t{start}\t{end}\t{phase1}\t{phase2}");
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
    if array.is_empty() {
        return NAN;
    }
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
fn test_region() {
    genotype_repeats(
        PathBuf::from("/home/wdecoster/test-data/test.bam"),
        Some("chr7:154778571-154779363".to_string()),
        None,
        0.5,
        5,
        4,
    );
}

#[test]
fn test_region_bed() {
    genotype_repeats(
        PathBuf::from("/home/wdecoster/test-data/test.bam"),
        None,
        Some(PathBuf::from("/home/wdecoster/test-data/test.bed")),
        0.5,
        5,
        4,
    );
}

#[test]
#[should_panic]
fn test_no_region() {
    genotype_repeats(
        PathBuf::from("/home/wdecoster/test-data/test.bam"),
        None,
        None,
        0.5,
        5,
        4,
    );
}

#[test]
#[should_panic]
fn test_wrong_bam_path() {
    genotype_repeats(
        PathBuf::from("/home/wdecoster/wrong_path_to_test-data/test.bam"),
        Some("chr7:154778571-154779363".to_string()),
        None,
        0.5,
        5,
        4,
    );
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
fn test_region_wrong_chromosome() {
    genotype_repeats(
        PathBuf::from("/home/wdecoster/test-data/test.bam"),
        Some("7:154778571-154779363".to_string()),
        None,
        0.5,
        5,
        4,
    );
}
