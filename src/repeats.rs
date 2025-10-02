use bio::io::bed;
use std::{collections::HashMap, fmt, path::PathBuf};

use crate::bam_utils::get_chrom_lengths_from_bam_header;

#[derive(Debug)]
pub struct RepeatIntervalIterator {
    current_index: usize,
    data: Vec<RepeatInterval>,
    pub num_intervals: usize,
}

impl RepeatIntervalIterator {
    // parse a region string
    pub fn from_string(reg: &str, chrom_lengths: HashMap<String, u64>) -> Self {
        let chrom = reg.split(':').collect::<Vec<&str>>()[0].to_string();
        let interval = reg.split(':').collect::<Vec<&str>>()[1];
        let start: u32 = interval.split('-').collect::<Vec<&str>>()[0]
            .parse()
            .unwrap();
        let end: u32 = interval.split('-').collect::<Vec<&str>>()[1]
            .parse()
            .unwrap();
        let repeat = RepeatInterval::new_interval(chrom, start, end, &chrom_lengths)
            .expect("Failed to create repeat interval");
        RepeatIntervalIterator { current_index: 0, data: vec![repeat], num_intervals: 1 }
    }
    pub fn from_bed(
        region_file: &String,
        chrom_lengths: HashMap<String, u64>,
        max_locus: Option<u32>,
    ) -> Self {
        let mut reader = bed::Reader::from_file(region_file).expect("Problem reading bed file!");
        let mut data = Vec::new();
        let mut filtered_count = 0;
        for record in reader.records() {
            let rec =
                record.expect("Error reading bed record. Is the file valid and tab-delimited?");
            let repeat = RepeatInterval::from_bed(&rec, &chrom_lengths);
            if let Some(repeat) = repeat {
                // Filter by max_locus size if specified
                let locus_size = repeat.end - repeat.start;
                if let Some(max_size) = max_locus {
                    if locus_size > max_size {
                        filtered_count += 1;
                        continue;
                    }
                }
                data.push(repeat);
            }
        }
        if filtered_count > 0 {
            eprintln!(
                "INFO: Filtered out {} intervals larger than {} bp (max-locus limit)",
                filtered_count,
                max_locus.unwrap()
            );
        }
        RepeatIntervalIterator { current_index: 0, data: data.clone(), num_intervals: data.len() }
    }
}

impl Clone for RepeatInterval {
    fn clone(&self) -> Self {
        RepeatInterval { chrom: self.chrom.clone(), start: self.start, end: self.end }
    }
}

impl Iterator for RepeatIntervalIterator {
    type Item = RepeatInterval;

    fn next(&mut self) -> Option<Self::Item> {
        // Implement the logic to get the next RepeatInterval here.
        // This is a simple example that gets the next item from a vector.
        if self.current_index < self.data.len() {
            let result = self.data[self.current_index].clone();
            self.current_index += 1;
            Some(result)
        } else {
            None
        }
    }
}

#[derive(Debug)]
pub struct RepeatInterval {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

impl fmt::Display for RepeatInterval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

impl RepeatInterval {
    // parse a bed record
    pub fn from_bed(rec: &bed::Record, chrom_lengths: &HashMap<String, u64>) -> Option<Self> {
        let chrom = rec.chrom().to_string();
        let start = rec.start().try_into().unwrap();
        let end = rec.end().try_into().unwrap();
        RepeatInterval::new_interval(chrom, start, end, chrom_lengths)
    }

    fn new_interval(
        chrom: String,
        start: u32,
        end: u32,
        chrom_lengths: &HashMap<String, u64>,
    ) -> Option<Self> {
        if end < start {
            panic!("End coordinate is smaller than start coordinate for {chrom}:{start}-{end}")
        }

        // check if the chromosome exists in the chrom lengths hashmap
        // and if the end coordinate is within the chromosome length
        if chrom_lengths.contains_key(&chrom) && (end as u64) < chrom_lengths[&chrom] {
            return Some(Self { chrom, start, end });
        }
        // if the chromosome is not in the fai file or the end does not fit the interval, return None
        panic!(
            "Chromosome {chrom} is not in the fasta file or the end coordinate is out of bounds"
        );
    }
    pub fn new(chrom: &str, start: u32, end: u32) -> Self {
        Self { chrom: chrom.to_string(), start, end }
    }
}

/// Get targets from region string or region file
pub fn get_targets(
    region: Option<String>,
    region_file: Option<PathBuf>,
    bam: &str,
    max_locus: Option<u32>,
) -> RepeatIntervalIterator {
    let chrom_lengths = get_chrom_lengths_from_bam_header(bam.to_string());
    match (&region, &region_file) {
        // a region string
        (Some(region), None) => RepeatIntervalIterator::from_string(region, chrom_lengths),
        // a region file
        (None, Some(region_file)) => RepeatIntervalIterator::from_bed(
            &region_file.to_string_lossy().to_string(),
            chrom_lengths,
            max_locus,
        ),
        // invalid input
        _ => {
            eprintln!("ERROR: Specify a region string (-r) or a region_file (-R)!\n");
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam_utils::get_chrom_lengths_from_bam_header;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn test_max_locus_filter() {
        // Test filtering functionality by creating a test BED file and checking results

        // Create a temporary test BED file
        let test_bed_content =
            "chr7\t154778571\t154779363\tsmall_interval\nchr7\t154780000\t154900000\thuge_interval\n";
        let mut file = File::create("test_temp_max_locus.bed").expect("Could not create test file");
        file.write_all(test_bed_content.as_bytes())
            .expect("Could not write test file");

        let chrom_lengths =
            get_chrom_lengths_from_bam_header(String::from("test-data/small-test.bam"));

        // Test without max_locus - should include both intervals
        let repeats_no_filter = RepeatIntervalIterator::from_bed(
            &String::from("test_temp_max_locus.bed"),
            chrom_lengths.clone(),
            None,
        );
        assert_eq!(repeats_no_filter.num_intervals, 2);

        // Test with max_locus 1000 - should filter out the huge interval (120000 bp)
        let repeats_with_filter = RepeatIntervalIterator::from_bed(
            &String::from("test_temp_max_locus.bed"),
            chrom_lengths,
            Some(1000),
        );
        assert_eq!(repeats_with_filter.num_intervals, 1);

        // Clean up
        std::fs::remove_file("test_temp_max_locus.bed").expect("Could not remove test file");
    }
}
