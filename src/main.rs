use bio::io::bed;
use clap::AppSettings::DeriveDisplayOrder;
use clap::{Parser, Subcommand};
use human_sort::compare as human_compare;
use log::{error, info};
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::f64::NAN;
use std::fmt;
use std::path::PathBuf;
use std::sync::Mutex;

// This struct keeps the genotype information and allows to compare them and thus sort them on chromosomal location
struct Genotype {
    chrom: String,
    start: u32,
    end: u32,
    phase1: f64,
    phase2: f64,
}

impl Ord for Genotype {
    fn cmp(&self, other: &Self) -> Ordering {
        human_compare(&self.chrom, &other.chrom).then(self.start.cmp(&other.start))
    }
}

impl PartialOrd for Genotype {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Genotype {
    fn eq(&self, other: &Self) -> bool {
        (self.chrom.clone(), &self.start) == (other.chrom.clone(), &other.start)
    }
}

impl Eq for Genotype {}

// How to print the struct, in bed-like format
impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.chrom, self.start, self.end, self.phase1, self.phase2
        )
    }
}
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
            unimplemented!();
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

/// This function genotypes STRs, either from a region string or from a bed file
/// For a bed file the genotyping is done in parallel
/// The minlen argument indicates the smallest CIGAR operation that is considered
/// Wobble extents the interval to make sure INDEL operations at the borders aren't missed
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
                Ok(output) => println!("{}", output),
                Err(chrom) => error!("Contig {chrom} not found in bam file"),
            };
        }
        (None, Some(region_file)) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();
            // TODO: check if bed file is okay
            let mut reader =
                bed::Reader::from_file(region_file.into_os_string().into_string().unwrap())
                    .unwrap();
            let bamf = bamp.into_os_string().into_string().unwrap();
            // chrom_reported and genotypes are vectors that are used by multiple threads to add findings, therefore as a Mutex
            // chrom_reported contains those chromosomes for which an error (absence in the bam) was already reported
            // to avoid reporthing the same error multiple times
            let chrom_reported = Mutex::new(Vec::new());
            // genotypes contains the output of the genotyping, a struct instance
            let genotypes = Mutex::new(Vec::new());
            // par_bridge does not guarantee that results are returned in order
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
                        // For now the Err is only used for when a chromosome from the bed file does not appear in the bam file
                        // this error is reported once per chromosome
                        let mut chroms_reported = chrom_reported.lock().unwrap();
                        if !chroms_reported.contains(&chrom) {
                            error!("Contig {chrom} not found in bam file");
                            chroms_reported.push(chrom);
                        }
                    }
                };
            });
            let mut genotypes_vec = genotypes.lock().unwrap();
            // The final output is sorted by chrom, start and end
            genotypes_vec.sort_unstable();
            for g in &mut *genotypes_vec {
                println!("{}", g);
            }
        }
    }
}

/// This function genotypes a particular repeat defined by chrom, start and end in the specified bam file
/// All indel cigar operations longer than minlen are considered
/// The bam file is expected to be phased using the HP tag
fn genotype_repeat(
    bamf: &String,
    chrom: String,
    start: u32,
    end: u32,
    minlen: u32,
) -> Result<Genotype, String> {
    let mut bam = match bam::IndexedReader::from_path(&bamf) {
        Ok(handle) => handle,
        Err(e) => {
            error!("Error opening BAM {}.\n{}", bamf, e);
            panic!();
        }
    };

    info!("Checks passed, genotyping repeat");
    if let Some(tid) = bam.header().tid(chrom.as_bytes()) {
        bam.fetch((tid, start, end)).unwrap();
        // Per haplotype the difference with the reference genome is kept in a dictionary
        let mut calls: HashMap<u8, Vec<i64>> =
            HashMap::from([(1, Vec::new()), (2, Vec::new()), (0, Vec::new())]);

        // CIGAR operations are assessed per read
        for r in bam.records() {
            let r = r.expect("Error reading BAM file in region {chrom}:{start}-{end}.");
            // reads with either end inside the window are ignored or if mapping quality is low
            if start < (r.reference_start() as u32)
                || (r.reference_end() as u32) < end
                || r.mapq() <= 10
            {
                continue;
            }
            // move the cursor for the reference position for all cigar operations that consume the reference
            let mut reference_position = (r.reference_start() + 1) as u32;
            let phase = get_phase(&r);
            let mut call: i64 = 0;
            for entry in r.cigar().iter() {
                match entry {
                    Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                        reference_position += *len;
                    }
                    Cigar::Del(len) => {
                        if *len > minlen && start < reference_position && reference_position < end {
                            call -= *len as i64;
                        }
                        reference_position += *len;
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
                    Cigar::RefSkip(len) => reference_position += *len,
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
        // For now, the unphased calls are ignored
        let output = Genotype {
            chrom,
            start,
            end,
            phase1: median_str_length(&calls[&1]),
            phase2: median_str_length(&calls[&2]),
        };
        Ok(output)
    } else {
        Err(chrom)
    }
}

/// parse a region string and extend the start and begin by a wobble space
/// defined by a wobble fraction relative to the interval length
fn process_region(
    reg: String,
    wobble: f64,
) -> Result<(String, u32, u32), Box<dyn std::error::Error>> {
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
        (start - wobble_length) as u32,
        (end + wobble_length) as u32,
    ))
}

/// Get the phase of a read by parsing the HP tag
/// The outcome should always be a u8
/// If the tag is absent '0' is returned, indicating unphased
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

/// Take the median of the lengths of the STRs, relative to the reference genome
/// If the vector is empty then return NAN
fn median_str_length(array: &Vec<i64>) -> f64 {
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
