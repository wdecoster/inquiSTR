use hts_sys;
use human_sort::compare as human_compare;
use indicatif::ParallelProgressIterator;
use indicatif::ProgressIterator;
use log::{error, info, warn};
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read};
use std::cmp::max;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::env;
use std::f64::NAN;
use std::fmt;
use std::io::{self, Write};
use std::path::PathBuf;
use std::rc::Rc;
use std::sync::Mutex;
use url::Url;

use crate::repeats::RepeatInterval;
use crate::repeats::RepeatIntervalIterator;

// This struct keeps the genotype information and allows to compare them and thus sort them on chromosomal location
struct Genotype {
    repeat: RepeatInterval,
    phase1: f64,
    phase2: f64,
}

impl Ord for Genotype {
    fn cmp(&self, other: &Self) -> Ordering {
        human_compare(&self.repeat.chrom, &other.repeat.chrom)
            .then(self.repeat.start.cmp(&other.repeat.start))
    }
}

impl PartialOrd for Genotype {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Genotype {
    fn eq(&self, other: &Self) -> bool {
        (self.repeat.chrom.clone(), &self.repeat.start)
            == (other.repeat.chrom.clone(), &other.repeat.start)
    }
}

impl Eq for Genotype {}

// How to print the struct, in bed-like format
// self.unphased should be 0, except if explicitly opted in to use unphased calls
impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.repeat.chrom, self.repeat.start, self.repeat.end, self.phase1, self.phase2
        )
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
    bamp: String,
    region: Option<String>,
    region_file: Option<PathBuf>,
    minlen: u32,
    support: usize,
    threads: usize,
    unphased: bool,
    sample_name: Option<String>,
) {
    if !PathBuf::from(&bamp).is_file() && !bamp.starts_with("s3") && !bamp.starts_with("https://") {
        error!("ERROR: path to bam file {} is not valid!\n\n", &bamp);
        std::process::exit(1);
    };
    let sample = sample_name.unwrap_or_else(|| {
        PathBuf::from(&bamp)
            .clone()
            .file_stem()
            .expect("Failed to get file stem")
            .to_str()
            .expect("Failed to convert to string")
            .replace(".bam", "")
            .replace(".cram", "")
    });
    let file_header = format!("chromosome\tbegin\tend\t{sample}_H1\t{sample}_H2");
    let repeats = get_targets(region, region_file, &bamp);
    if threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to build thread pool");
        // chrom_reported and genotypes are vectors that are used by multiple threads to add findings, therefore as a Mutex
        // chrom_reported contains those chromosomes for which an error (absence in the bam) was already reported
        // to avoid reporting the same error multiple times
        let chrom_reported = Mutex::new(Vec::new());
        // genotypes contains the output of the genotyping, a struct instance
        let genotypes = Mutex::new(Vec::new());
        let num_intervals = repeats.num_intervals;
        repeats
            .par_bridge()
            .progress_count(num_intervals as u64)
            .for_each(|repeat| {
                match genotype_repeat_multithreaded(&bamp, repeat, minlen, support, unphased) {
                    Ok(output) => {
                        let mut geno = genotypes.lock().expect("Failed to lock genotypes");
                        geno.push(output);
                    }
                    Err(locus) => {
                        // For now the Err is only used for when a chromosome or (extended) interval from the bed file does not appear in the bam file
                        // this error is reported once per locus
                        let mut chroms_reported = chrom_reported
                            .lock()
                            .expect("Failed to lock chrom_reported");
                        if !chroms_reported.contains(&locus) {
                            warn!("{locus} not found in bam file");
                            chroms_reported.push(locus);
                        }
                    }
                };
            });
        let mut genotypes_vec = genotypes.lock().expect("Failed to lock genotypes");
        // The final output is sorted by chrom, start and end
        let stdout = io::stdout(); // get the global stdout entity
        let mut handle = io::BufWriter::new(stdout); // optional: wrap that handle in a buffer
        genotypes_vec.sort_unstable();
        writeln!(handle, "{file_header}").expect("Failed writing the header.");
        for g in &mut *genotypes_vec {
            writeln!(handle, "{g}").expect("Failed writing the result.");
        }
    } else {
        let mut bam = get_bam_reader(&bamp);
        let num_intervals = repeats.num_intervals;
        println!("{file_header}");
        for repeat in repeats.progress_count(num_intervals as u64) {
            match genotype_repeat(&mut bam, repeat, minlen, support, unphased) {
                Ok(output) => {
                    println!("{output}")
                }
                Err(chrom) => error!("Contig {chrom} not found in bam file"),
            };
        }
    }
}

fn get_chrom_lengths_from_bam_header(bam: String) -> HashMap<String, u64> {
    let bam = get_bam_reader(&bam);
    let header = bam::Header::from_template(bam.header());
    let mut chrom_lengts = HashMap::new();
    for (key, records) in header.to_hashmap() {
        for record in records {
            if key != "SQ" {
                continue;
            }
            chrom_lengts.insert(
                record["SN"].clone(),
                record["LN"]
                    .parse()
                    .expect("Failed to parse length of chromosome"),
            );
        }
    }

    chrom_lengts
}

fn get_targets(
    region: Option<String>,
    region_file: Option<PathBuf>,
    bam: &str,
) -> RepeatIntervalIterator {
    let chrom_lengths = get_chrom_lengths_from_bam_header(bam.to_string());
    match (&region, &region_file) {
        // a region string
        (Some(region), None) => RepeatIntervalIterator::from_string(region, chrom_lengths),
        // a region file
        (None, Some(region_file)) => RepeatIntervalIterator::from_bed(
            &region_file.to_string_lossy().to_string(),
            chrom_lengths,
        ),
        // invalid input
        _ => {
            eprintln!("ERROR: Specify a region string (-r) or a region_file (-R)!\n");
            std::process::exit(1);
        }
    }
}

/// This function genotypes a particular repeat defined by chrom, start and end in the specified bam file
/// All indel cigar operations longer than minlen are considered
/// The bam file is expected to be phased using the HP tag
/// This function is specific to multithreaded use, as it takes a String for the bam rather than the Reader
/// The function below is the single threaded version
fn genotype_repeat_multithreaded(
    bamf: &String,
    repeat: RepeatInterval,
    minlen: u32,
    support: usize,
    unphased: bool,
) -> Result<Genotype, String> {
    let mut bam = get_bam_reader(bamf);
    info!("Checks passed, genotyping repeat");
    if unphased {
        genotype_repeat_unphased(&mut bam, repeat, minlen, support)
    } else {
        genotype_repeat_phased(&mut bam, repeat, minlen, support)
    }
}

fn get_bam_reader(bamp: &String) -> bam::IndexedReader {
    let mut bam = if bamp.starts_with("s3") || bamp.starts_with("https://") {
        if env::var("CURL_CA_BUNDLE").is_err() {
            env::set_var("CURL_CA_BUNDLE", "/etc/ssl/certs/ca-certificates.crt");
        }
        bam::IndexedReader::from_url(&Url::parse(bamp.as_str()).expect("Failed to parse s3 URL"))
            .unwrap_or_else(|err| panic!("Error opening remote BAM: {err}"))
    } else {
        bam::IndexedReader::from_path(bamp)
            .unwrap_or_else(|err| panic!("Error opening local BAM: {err}"))
    };
    if bamp.ends_with(".cram") {
        bam.set_cram_options(
            hts_sys::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
            hts_sys::sam_fields_SAM_AUX
                | hts_sys::sam_fields_SAM_MAPQ
                | hts_sys::sam_fields_SAM_CIGAR,
        )
        .expect("Failed setting cram options");
    }
    bam
}

fn genotype_repeat(
    bam: &mut bam::IndexedReader,
    repeat: RepeatInterval,
    minlen: u32,
    support: usize,
    unphased: bool,
) -> Result<Genotype, String> {
    info!("Checks passed, genotyping repeat");
    if unphased {
        genotype_repeat_unphased(bam, repeat, minlen, support)
    } else {
        genotype_repeat_phased(bam, repeat, minlen, support)
    }
}

fn genotype_repeat_unphased(
    bam: &mut bam::IndexedReader,
    repeat: RepeatInterval,
    minlen: u32,
    support: usize,
) -> Result<Genotype, String> {
    let start_ext = max(repeat.start - 10, 0);
    let end_ext = repeat.end + 10;
    if let Some(tid) = bam.header().tid(repeat.chrom.as_bytes()) {
        bam.fetch((tid, start_ext, end_ext)).unwrap();
        // Per haplotype the difference with the reference genome is kept in a dictionary
        let mut calls = vec![];

        // CIGAR operations are assessed per read
        for r in bam.rc_records() {
            let r = r.expect("Error reading BAM file in region {chrom}:{start}-{end}.");
            // reads with either end inside the window are ignored or if mapping quality is low
            if start_ext < (r.reference_start() as u32)
                || (r.reference_end() as u32) < end_ext
                || r.mapq() <= 10
            {
                continue;
            }
            let call = call_from_cigar(r, minlen, start_ext, end_ext);
            calls.push(call);
        }
        info!("Found {} reads for genotyping", calls.len(),);
        // sort the vec of calls based on the value
        let f = |c: &Call| match c {
            Call::Span(v) => *v,
            Call::Clip(v) => *v,
        };
        calls.sort_unstable_by_key(f);
        // split both haplotypes with median split, split_at divides one slice into two at an index.
        let (h1, h2) = calls.split_at(calls.len() / 2);

        // unphased is set to 0 if those are to be ignored and vice versa
        // just taking the median of unphased reads is not optimal
        let output = Genotype {
            repeat,
            phase1: median_str_length(&h1.to_vec(), support),
            phase2: median_str_length(&h2.to_vec(), support),
        };
        Ok(output)
    } else {
        Err(repeat.chrom)
    }
}

fn genotype_repeat_phased(
    bam: &mut bam::IndexedReader,
    repeat: RepeatInterval,
    minlen: u32,
    support: usize,
) -> Result<Genotype, String> {
    let start_ext = max(repeat.start - 10, 0);
    let end_ext = repeat.end + 10;
    if let Some(tid) = bam.header().tid(repeat.chrom.as_bytes()) {
        match bam.fetch((tid, start_ext, end_ext)) {
            Ok(_b) => (),
            Err(e) => return Err(e.to_string()),
        };

        // Per haplotype the difference with the reference genome is kept in a dictionary
        let mut calls: HashMap<u8, Vec<Call>> =
            HashMap::from([(1, Vec::new()), (2, Vec::new()), (0, Vec::new())]);

        // CIGAR operations are assessed per read
        for r in bam.rc_records() {
            let r = r.expect("Error reading BAM file in region {chrom}:{start}-{end}.");
            // reads with both ends inside the window are ignored or if mapping quality is low
            // since the bam is supposed to be phased, ignore all unphased reads
            let phase = get_phase(&r);
            if phase.is_none() 
                || start_ext < (r.reference_start() as u32) && (r.reference_end() as u32) < end_ext
                || r.mapq() <= 10
            {
                continue;
            }

            let call = call_from_cigar(r, minlen, start_ext, end_ext);
            calls.get_mut(&phase.expect("Couldn't get phase - this shouldn't happen")).unwrap().push(call);
        }
        info!(
            "Found {}[H1]+{}[H2] reads for genotyping",
            calls[&1].len(),
            calls[&2].len()
        );
        let output = Genotype {
            repeat,
            phase1: median_str_length(&calls[&1].clone(), support),
            phase2: median_str_length(&calls[&2].clone(), support),
        };
        Ok(output)
    } else {
        Err(repeat.chrom)
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
fn get_phase(record: &bam::Record) -> Option<u8> {
    match record.aux(b"HP") {
        Ok(value) => match value {
            Aux::U8(v) => Some(v),
            Aux::I32(v) => Some(v as u8),
            _ => panic!("Unexpected type of Aux {value:?}"),
        },
        Err(_e) => None,
    }
}

/// Take the median of the lengths of the STRs, relative to the reference genome
/// If the vector has fewer than <support> calls then return NAN
/// Spanning reads have the preference, so if more than <support> spanning reads are present the median is calculated for those
/// Otherwise, the longest softclipped reads are added up to <support> reads
fn median_str_length(array: &Vec<Call>, support: usize) -> f64 {
    if array.len() < support {
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
        spanning[spanning.len() / 2] as f64
    }
}

#[cfg(test)]
#[test]
fn test_region() {
    genotype_repeats(
        String::from("test-data/small-test.bam"),
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
fn test_region_from_url() {
    genotype_repeats(
        String::from("https://s3.amazonaws.com/1000g-ont/FIRST_100_FREEZE/minimap2_2.24_alignment_data/GM18501/GM18501.LSK110.R9.guppy646.sup.with5mC.pass.phased.bam"),
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
        String::from("test-data/small-test.bam"),
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
        String::from("test-data/small-test.bam"),
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
fn test_region_wrong_chromosome() {
    genotype_repeats(
        String::from("test-data/small-test.bam"),
        Some("7:154778571-154779363".to_string()),
        None,
        5,
        3,
        4,
        false,
        Some("sample".to_string()),
    );
}

#[test]
fn test_get_chrom_lengths_from_bam_header() {
    let bam = String::from("test-data/small-test.bam");
    let chrom_lengths = get_chrom_lengths_from_bam_header(bam);
    assert_eq!(chrom_lengths.get("chr7").unwrap(), &159345973);
}
