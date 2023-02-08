use bio::io::bed;
use human_sort::compare as human_compare;
use indicatif::ParallelProgressIterator;
use indicatif::ProgressIterator;
use linecount::count_lines;
use log::{error, info, warn};
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read};
use std::cmp::max;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::f64::NAN;
use std::fmt;
use std::path::PathBuf;
use std::rc::Rc;
use std::sync::Mutex;

// This struct keeps the genotype information and allows to compare them and thus sort them on chromosomal location
struct Genotype {
    chrom: String,
    start: u32,
    end: u32,
    phase1: f64,
    phase2: f64,
    unphased: f64,
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
// self.unphased should be 0, except if explicitly opted in to use unphased calls
impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.unphased > 0.0 {
            write!(
                f,
                "{}\t{}\t{}\t{}",
                self.chrom, self.start, self.end, self.unphased
            )
        } else {
            write!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                self.chrom, self.start, self.end, self.phase1, self.phase2
            )
        }
    }
}

#[derive(Clone)]
enum Call {
    Span(i64),
    Clip(i64),
}

/// This function genotypes STRs, either from a region string or from a bed file
/// For a bed file the genotyping is done in parallel
/// The minlen argument indicates the smallest CIGAR operation that is considered
pub fn genotype_repeats(
    bamp: PathBuf,
    region: Option<String>,
    region_file: Option<PathBuf>,
    minlen: u32,
    support: usize,
    threads: usize,
    unphased: bool,
    sample_name: Option<String>,
) {
    if !bamp.is_file() {
        error!(
            "ERROR: path to bam file {} is not valid!\n\n",
            bamp.display()
        );
        panic!();
    };
    let sample = sample_name.unwrap_or_else(|| {
        bamp.clone()
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .replace(".bam", "")
    });
    let header = if unphased {
        format!("chromosome\tbegin\tend\t{sample}")
    } else {
        format!("chromosome\tbegin\tend\t{sample}_H1\t{sample}_H2")
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
            let (chrom, start, end) = crate::utils::process_region(region).unwrap();
            let bamf = bamp.into_os_string().into_string().unwrap();

            match genotype_repeat(&bamf, chrom, start, end, minlen, support, unphased) {
                Ok(output) => {
                    println!("{header}");
                    println!("{output}")
                }
                Err(chrom) => error!("Contig {chrom} not found in bam file"),
            };
        }
        (None, Some(region_file)) => {
            if threads > 1 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                // TODO: check if bed file is okay
                let mut reader = bed::Reader::from_file(region_file.clone()).unwrap();
                let bamf = bamp.into_os_string().into_string().unwrap();
                // chrom_reported and genotypes are vectors that are used by multiple threads to add findings, therefore as a Mutex
                // chrom_reported contains those chromosomes for which an error (absence in the bam) was already reported
                // to avoid reporting the same error multiple times
                let chrom_reported = Mutex::new(Vec::new());
                // genotypes contains the output of the genotyping, a struct instance
                let genotypes = Mutex::new(Vec::new());
                // par_bridge does not guarantee that results are returned in order
                let lines: usize = count_lines(std::fs::File::open(region_file).unwrap()).unwrap();
                println!("{header}");
                reader
                    .records()
                    .par_bridge()
                    .progress_count(lines as u64)
                    .for_each(|record| {
                        let rec = record.expect("Error reading bed record.");
                        match genotype_repeat(
                            &bamf,
                            rec.chrom().to_string(),
                            rec.start().try_into().unwrap(),
                            rec.end().try_into().unwrap(),
                            minlen,
                            support,
                            unphased,
                        ) {
                            Ok(output) => {
                                let mut geno = genotypes.lock().unwrap();
                                geno.push(output);
                            }
                            Err(locus) => {
                                // For now the Err is only used for when a chromosome or (extended) interval from the bed file does not appear in the bam file
                                // this error is reported once per locus
                                let mut chroms_reported = chrom_reported.lock().unwrap();
                                if !chroms_reported.contains(&locus) {
                                    warn!("{locus} not found in bam file");
                                    chroms_reported.push(locus);
                                }
                            }
                        };
                    });
                let mut genotypes_vec = genotypes.lock().unwrap();
                // The final output is sorted by chrom, start and end
                genotypes_vec.sort_unstable();
                for g in &mut *genotypes_vec {
                    println!("{g}");
                }
            } else {
                // When running single threaded things become easier and the tool will require less memory
                // Output is returned in the same order as the bed, and therefore not sorted before writing to stdout
                // TODO: check if bed file is okay
                let mut reader = bed::Reader::from_file(region_file.clone()).unwrap();
                let bamf = bamp.into_os_string().into_string().unwrap();
                // chrom_reported contains those chromosomes for which an error (absence in the bam) was already reported
                // to avoid reporting the same error multiple times
                let mut chrom_reported = Vec::new();
                // genotypes contains the output of the genotyping, a struct instance
                let lines: usize = count_lines(std::fs::File::open(region_file).unwrap()).unwrap();
                for record in reader.records().progress_count(lines as u64) {
                    let rec = record.expect("Error reading bed record.");
                    match genotype_repeat(
                        &bamf,
                        rec.chrom().to_string(),
                        rec.start().try_into().unwrap(),
                        rec.end().try_into().unwrap(),
                        minlen,
                        support,
                        unphased,
                    ) {
                        Ok(output) => {
                            println!("{output}");
                        }
                        Err(locus) => {
                            // For now the Err is only used for when a chromosome or (extended) interval from the bed file does not appear in the bam file
                            // this error is reported once per locus
                            if !chrom_reported.contains(&locus) {
                                warn!("{locus} not found in bam file");
                                chrom_reported.push(locus);
                            }
                        }
                    };
                }
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
    support: usize,
    unphased: bool,
) -> Result<Genotype, String> {
    let bam = match bam::IndexedReader::from_path(bamf) {
        Ok(handle) => handle,
        Err(e) => {
            error!("Error opening BAM {}.\n{}", bamf, e);
            panic!();
        }
    };
    let start = max(start - 10, 0);
    let end = end + 10;
    info!("Checks passed, genotyping repeat");
    if unphased {
        genotype_repeat_unphased(bam, chrom, start, end, minlen, support)
    } else {
        genotype_repeat_phased(bam, chrom, start, end, minlen, support)
    }
}

fn genotype_repeat_unphased(
    mut bam: bam::IndexedReader,
    chrom: String,
    start: u32,
    end: u32,
    minlen: u32,
    support: usize,
) -> Result<Genotype, String> {
    if let Some(tid) = bam.header().tid(chrom.as_bytes()) {
        bam.fetch((tid, start, end)).unwrap();
        // Per haplotype the difference with the reference genome is kept in a dictionary
        let mut calls = vec![];

        // CIGAR operations are assessed per read
        for r in bam.rc_records() {
            let r = r.expect("Error reading BAM file in region {chrom}:{start}-{end}.");
            // reads with either end inside the window are ignored or if mapping quality is low
            if start < (r.reference_start() as u32)
                || (r.reference_end() as u32) < end
                || r.mapq() <= 10
            {
                continue;
            }
            let call = call_from_cigar(r, minlen, start, end);
            calls.push(call);
        }
        info!("Found {} reads for genotyping", calls.len(),);
        // unphased is set to 0 if those are to be ignored and vice versa
        // just taking the median of unphased reads is not optimal
        let output = Genotype {
            chrom,
            start,
            end,
            phase1: 0.0,
            phase2: 0.0,
            unphased: median_str_length(&calls.clone(), support),
        };
        Ok(output)
    } else {
        Err(chrom)
    }
}

fn genotype_repeat_phased(
    mut bam: bam::IndexedReader,
    chrom: String,
    start: u32,
    end: u32,
    minlen: u32,
    support: usize,
) -> Result<Genotype, String> {
    if let Some(tid) = bam.header().tid(chrom.as_bytes()) {
        match bam.fetch((tid, start, end)) {
            Ok(_b) => (),
            Err(e) => return Err(e.to_string()),
        };

        // Per haplotype the difference with the reference genome is kept in a dictionary
        let mut calls: HashMap<u8, Vec<Call>> =
            HashMap::from([(1, Vec::new()), (2, Vec::new()), (0, Vec::new())]);

        // CIGAR operations are assessed per read
        for r in bam.rc_records() {
            let r = r.expect("Error reading BAM file in region {chrom}:{start}-{end}.");
            // reads with either end inside the window are ignored or if mapping quality is low
            // if the bam is supposed to be phased, ignore all unphased reads
            let phase = get_phase(&r);
            if start < (r.reference_start() as u32)
                || (r.reference_end() as u32) < end
                || r.mapq() <= 10
                || phase == 0
            {
                continue;
            }

            let call = call_from_cigar(r, minlen, start, end);
            calls.get_mut(&phase).unwrap().push(call);
        }
        info!(
            "Found {}[H1]+{}[H2] reads for genotyping",
            calls[&1].len(),
            calls[&2].len()
        );
        // unphased is set to 0 if those are to be ignored and vice versa
        let output = Genotype {
            chrom,
            start,
            end,
            phase1: median_str_length(&calls[&1].clone(), support),
            phase2: median_str_length(&calls[&2].clone(), support),
            unphased: 0.0,
        };
        Ok(output)
    } else {
        Err(chrom)
    }
}

fn call_from_cigar(r: Rc<bam::Record>, minlen: u32, start: u32, end: u32) -> Call {
    let mut call: i64 = 0;
    // move the cursor for the reference position for all cigar operations that consume the reference
    let mut reference_position = (r.reference_start() + 1) as u32;
    let mut clipped = false;
    for entry in r.cigar().iter() {
        match entry {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                reference_position += *len;
            }
            Cigar::Del(len) => {
                if *len > minlen && start < reference_position && reference_position < end {
                    call -= i64::from(*len);
                }
                reference_position += *len;
            }
            Cigar::SoftClip(len) => {
                if *len > minlen && start < reference_position && reference_position < end {
                    call += i64::from(*len);
                    clipped = true
                }
            }
            Cigar::Ins(len) => {
                if *len > minlen && start < reference_position && reference_position < end {
                    call += i64::from(*len);
                }
            }
            Cigar::RefSkip(len) => reference_position += *len,
            _ => (),
        }
    }
    if clipped {
        Call::Clip(call)
    } else {
        Call::Span(call)
    }
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
                panic!("Unexpected type of Aux {value:?}")
            }
        }
        Err(_e) => 0,
    }
}

/// Take the median of the lengths of the STRs, relative to the reference genome
/// If the vector has fewer than <support> calls then return NAN
/// Spanning reads have the preference, so if more than <support> spanning reads are present the median is calculated for those
/// Otherwise, the longest softclipped reads are added up to <support> reads
fn median_str_length(array: &Vec<Call>, support: usize) -> f64 {
    if array.len() <= support {
        return NAN;
    }
    let mut spanning = vec![];
    let mut clipped = vec![];
    for a in array {
        match a {
            Call::Span(v) => spanning.push(*v),
            Call::Clip(v) => clipped.push(*v),
        };
    }
    if spanning.len() <= support {
        // Sort clipped from large to small to the largest clips
        clipped.sort_unstable_by_key(|k| -k);
        spanning.extend(&clipped[0..support - spanning.len()]);
    }
    spanning.sort_unstable();
    if (spanning.len() % 2) == 0 {
        let ind_left = spanning.len() / 2 - 1;
        let ind_right = spanning.len() / 2;
        (spanning[ind_left] + spanning[ind_right]) as f64 / 2.0
    } else {
        spanning[(spanning.len() / 2)] as f64
    }
}

#[cfg(test)]
#[test]
fn test_region() {
    genotype_repeats(
        PathBuf::from("test-data/small-test.bam"),
        Some("chr7:154778571-154779363".to_string()),
        None,
        5,
        3,
        4,
        false,
        Some("sample".to_string()),
    );
}

#[test]
fn test_region_bed() {
    genotype_repeats(
        PathBuf::from("test-data/small-test.bam"),
        None,
        Some(PathBuf::from("test-data/test.bed")),
        5,
        3,
        4,
        false,
        Some("sample".to_string()),
    );
}
#[test]
fn test_unphased() {
    genotype_repeats(
        PathBuf::from("test-data/small-test.bam"),
        None,
        Some(PathBuf::from("test-data/test.bed")),
        5,
        3,
        4,
        true,
        Some("sample".to_string()),
    );
}

#[test]
#[should_panic]
fn test_no_region() {
    genotype_repeats(
        PathBuf::from("test-data/small-test.bam"),
        None,
        None,
        5,
        3,
        4,
        false,
        Some("sample".to_string()),
    );
}

#[test]
#[should_panic]
fn test_wrong_bam_path() {
    genotype_repeats(
        PathBuf::from("test-data/wrong-test.bam"),
        Some("chr7:154778571-154779363".to_string()),
        None,
        5,
        3,
        4,
        false,
        Some("sample".to_string()),
    );
}

#[test]
fn test_region_wrong_chromosome() {
    genotype_repeats(
        PathBuf::from("test-data/small-test.bam"),
        Some("7:154778571-154779363".to_string()),
        None,
        5,
        3,
        4,
        false,
        Some("sample".to_string()),
    );
}
