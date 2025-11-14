use bio::io::bed;
use std::{collections::HashMap, fmt};

use crate::bam_utils::get_chrom_lengths_from_bam_header;

#[derive(Debug)]
pub struct RepeatIntervalIterator {
    current_index: usize,
    data: Vec<RepeatInterval>,
    pub num_intervals: usize,
}

impl RepeatIntervalIterator {
    // parse a region string in format "chr:start-end"
    pub fn from_string(reg: &str, chrom_lengths: HashMap<String, u64>) -> Self {
        let reg = reg.trim();

        // Validate format
        let parts: Vec<&str> = reg.split(':').collect();
        if parts.len() != 2 {
            panic!("Invalid region format '{}'. Expected format: chr:start-end", reg);
        }

        let chrom = parts[0].to_string();
        let interval = parts[1];

        let interval_parts: Vec<&str> = interval.split('-').collect();
        if interval_parts.len() != 2 {
            panic!("Invalid region format '{}'. Expected format: chr:start-end", reg);
        }

        let start: u32 = interval_parts[0]
            .parse()
            .unwrap_or_else(|_| panic!("Invalid start position in region '{}'", reg));
        let end: u32 = interval_parts[1]
            .parse()
            .unwrap_or_else(|_| panic!("Invalid end position in region '{}'", reg));

        let repeat = RepeatInterval::new_interval(chrom, start, end, &chrom_lengths)
            .expect("Failed to create repeat interval");
        RepeatIntervalIterator { current_index: 0, data: vec![repeat], num_intervals: 1 }
    }
    pub fn from_bed(
        region_file: &str,
        chrom_lengths: HashMap<String, u64>,
        max_locus: Option<u32>,
    ) -> Self {
        // Use utils::reader to handle gzipped files
        let file_reader = crate::utils::reader(region_file);
        let mut reader = bed::Reader::new(file_reader);
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

    /// Download and cache a predefined TR catalog preset
    ///
    /// This method handles downloading TR catalogs from remote URLs, caching them
    /// locally for 7 days, and handling network failures gracefully by falling back
    /// to cached versions when available.
    pub fn from_preset(
        preset: crate::call::TRPreset,
        chrom_lengths: HashMap<String, u64>,
        max_locus: Option<u32>,
    ) -> Self {
        let (url, cache_filename) = preset.metadata();
        let preset_name = preset.display_name();

        let cache_dir = dirs::cache_dir()
            .unwrap_or_else(std::env::temp_dir)
            .join("inquistr");

        // Create cache directory if it doesn't exist
        if let Err(e) = std::fs::create_dir_all(&cache_dir) {
            eprintln!("ERROR: Failed to create cache directory {}: {}", cache_dir.display(), e);
            std::process::exit(1);
        }

        let cache_file = cache_dir.join(cache_filename);

        // Check if cached file exists and is recent (less than 7 days old)
        let needs_download = if cache_file.exists() {
            match std::fs::metadata(&cache_file) {
                Ok(metadata) => {
                    if let Ok(modified) = metadata.modified() {
                        let age = std::time::SystemTime::now()
                            .duration_since(modified)
                            .unwrap_or(std::time::Duration::from_secs(0));
                        age > std::time::Duration::from_secs(7 * 24 * 60 * 60) // 7 days
                    } else {
                        true // Can't get modification time, re-download
                    }
                }
                Err(_) => true, // Can't get metadata, re-download
            }
        } else {
            true // File doesn't exist, download
        };

        // Download if needed
        if needs_download {
            eprintln!("Downloading {} catalog...", preset_name);
            match reqwest::blocking::get(url) {
                Ok(resp) => {
                    if !resp.status().is_success() {
                        let status = resp.status();
                        // Fall back to cache if download fails
                        if cache_file.exists() {
                            eprintln!(
                                "Warning: Download returned status {}, using cached version",
                                status
                            );
                        } else {
                            eprintln!(
                                "ERROR: Failed to download {} catalog: HTTP status {}",
                                preset_name, status
                            );
                            eprintln!(
                                "The catalog URL may have changed. Please check for updates or report this issue."
                            );
                            std::process::exit(1);
                        }
                    } else {
                        match resp.bytes() {
                            Ok(body) => {
                                if let Err(e) = std::fs::write(&cache_file, &body) {
                                    eprintln!(
                                        "Warning: Failed to cache {} catalog: {}",
                                        preset_name, e
                                    );
                                    // Continue with in-memory data - decompress if gzipped
                                    if cache_filename.ends_with(".gz") {
                                        match Self::decompress_gzip(&body) {
                                            Ok(decompressed) => {
                                                return Self::from_string_data(
                                                    &decompressed,
                                                    chrom_lengths,
                                                    max_locus,
                                                );
                                            }
                                            Err(e) => {
                                                eprintln!("ERROR: Failed to decompress gzipped catalog: {}", e);
                                                std::process::exit(1);
                                            }
                                        }
                                    } else {
                                        let text = String::from_utf8_lossy(&body).to_string();
                                        return Self::from_string_data(
                                            &text,
                                            chrom_lengths,
                                            max_locus,
                                        );
                                    }
                                }
                                eprintln!(
                                    "Cached {} catalog to: {}",
                                    preset_name,
                                    cache_file.display()
                                );
                            }
                            Err(e) => {
                                // Fall back to cache if reading response fails
                                if cache_file.exists() {
                                    eprintln!("Warning: Failed to read response body ({}), using cached version", e);
                                } else {
                                    eprintln!(
                                        "ERROR: Failed to download {} catalog: {}",
                                        preset_name, e
                                    );
                                    eprintln!(
                                        "Please check your network connection and try again."
                                    );
                                    std::process::exit(1);
                                }
                            }
                        }
                    }
                }
                Err(e) => {
                    // If download fails but we have a cached version, use it even if old
                    if cache_file.exists() {
                        eprintln!("Warning: Download failed ({}), using cached version", e);
                    } else {
                        eprintln!("ERROR: Failed to download {} catalog: {}", preset_name, e);
                        eprintln!("Please check your network connection and try again.");
                        eprintln!("URL: {}", url);
                        std::process::exit(1);
                    }
                }
            }
        }

        // Read from cached file (possibly gzipped)
        if cache_filename.ends_with(".gz") {
            Self::from_gzipped_bed(&cache_file.to_string_lossy(), chrom_lengths, max_locus)
        } else {
            Self::from_bed(&cache_file.to_string_lossy(), chrom_lengths, max_locus)
        }
    }

    /// Decompress gzipped data
    fn decompress_gzip(data: &[u8]) -> Result<String, std::io::Error> {
        use flate2::read::GzDecoder;
        use std::io::Read;

        let mut decoder = GzDecoder::new(data);
        let mut decompressed = String::new();
        decoder.read_to_string(&mut decompressed)?;
        Ok(decompressed)
    }

    /// Read from a gzipped BED file
    fn from_gzipped_bed(
        path: &str,
        chrom_lengths: HashMap<String, u64>,
        max_locus: Option<u32>,
    ) -> Self {
        use flate2::read::GzDecoder;
        use std::io::Read;

        let file = std::fs::File::open(path).unwrap_or_else(|e| {
            eprintln!("ERROR: Failed to open gzipped BED file {}: {}", path, e);
            std::process::exit(1);
        });

        let mut decoder = GzDecoder::new(file);
        let mut contents = String::new();
        decoder.read_to_string(&mut contents).unwrap_or_else(|e| {
            eprintln!("ERROR: Failed to decompress gzipped BED file {}: {}", path, e);
            std::process::exit(1);
        });

        Self::from_string_data(&contents, chrom_lengths, max_locus)
    }

    // Helper function to process data from string (for fallback)
    fn from_string_data(
        data: &str,
        chrom_lengths: HashMap<String, u64>,
        max_locus: Option<u32>,
    ) -> Self {
        let reader = std::io::BufReader::new(data.as_bytes());
        let mut bed_reader = bed::Reader::new(reader);
        let mut data_vec = Vec::new();
        let mut filtered_count = 0;

        for record in bed_reader.records() {
            let rec = record.expect("Error reading bed record from downloaded data");
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
                data_vec.push(repeat);
            }
        }

        if filtered_count > 0 {
            eprintln!(
                "INFO: Filtered out {} intervals larger than {} bp (max-locus limit)",
                filtered_count,
                max_locus.unwrap()
            );
        }

        RepeatIntervalIterator {
            current_index: 0,
            data: data_vec.clone(),
            num_intervals: data_vec.len(),
        }
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
            eprintln!("ERROR: Invalid coordinates. End position ({}) is less than start position ({}) for {}:{}-{}", 
                end, start, chrom, start, end);
            std::process::exit(1);
        }

        // check if the chromosome exists in the chrom lengths hashmap
        // and if the end coordinate is within the chromosome length
        if chrom_lengths.contains_key(&chrom) && (end as u64) < chrom_lengths[&chrom] {
            return Some(Self { chrom, start, end });
        }
        // if the chromosome is not in the fai file or the end does not fit the interval, return None
        eprintln!(
            "ERROR: Chromosome '{}' not found in reference or coordinate {} is out of bounds",
            chrom, end
        );
        eprintln!(
            "Available chromosomes: {}",
            chrom_lengths
                .keys()
                .map(|s| s.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        );
        std::process::exit(1);
    }
    pub fn new(chrom: &str, start: u32, end: u32) -> Self {
        Self { chrom: chrom.to_string(), start, end }
    }
}

/// Get targets from region string, region file, or preset catalog
pub fn get_targets(
    targets: crate::call::TargetConfig,
    bam: &str,
    reference: &Option<String>,
) -> RepeatIntervalIterator {
    let chrom_lengths = get_chrom_lengths_from_bam_header(bam.to_string(), reference);
    match (&targets.region, &targets.region_file, &targets.preset) {
        // a region string
        (Some(region), None, None) => RepeatIntervalIterator::from_string(region, chrom_lengths),
        // a region file
        (None, Some(region_file), None) => RepeatIntervalIterator::from_bed(
            &region_file.to_string_lossy(),
            chrom_lengths,
            targets.max_locus,
        ),
        // preset catalog
        (None, None, Some(preset)) => {
            RepeatIntervalIterator::from_preset(*preset, chrom_lengths, targets.max_locus)
        }
        // invalid input
        _ => {
            eprintln!(
                "ERROR: Specify a region string (-r), a region_file (-R), or --preset <pathogenic|adotto|trexplorer>!\n"
            );
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
            get_chrom_lengths_from_bam_header(String::from("test-data/small-test.bam"), &None);

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
