use bio::io::bed;
use clap::ValueEnum;
use std::{collections::HashMap, path::PathBuf};

use crate::bam_utils::get_chrom_lengths_from_bam_header;
use crate::errors::{InquiSTRError, InquiSTRResult};

/// Maps chromosome names to numeric IDs to reduce memory usage
/// Instead of storing chromosome names as Strings in every RepeatInterval,
/// we store a small u32 ID and use this mapper for lookups during I/O
#[derive(Debug, Clone, Default)]
pub struct ChromosomeMapper {
    // Name to ID lookup
    name_to_id: HashMap<String, u32>,
    // ID to name lookup (stored as Vec for fast indexed access)
    id_to_name: Vec<String>,
}

impl ChromosomeMapper {
    pub fn new() -> Self {
        Self::default()
    }

    /// Get or create an ID for a chromosome name
    pub fn get_or_insert(&mut self, chrom: &str) -> u32 {
        if let Some(&id) = self.name_to_id.get(chrom) {
            id
        } else {
            let id = self.id_to_name.len() as u32;
            self.id_to_name.push(chrom.to_string());
            self.name_to_id.insert(chrom.to_string(), id);
            id
        }
    }

    /// Get chromosome name from ID
    pub fn get_name(&self, id: u32) -> &str {
        &self.id_to_name[id as usize]
    }

    /// Get chromosome ID from name (returns None if not found)
    pub fn get_id(&self, chrom: &str) -> Option<u32> {
        self.name_to_id.get(chrom).copied()
    }
}

/// Predefined tandem repeat (TR) catalogs for genotyping
///
/// These presets allow quick access to well-known TR catalogs without manually
/// downloading and specifying BED files. Downloaded catalogs are cached locally
/// for 7 days to avoid repeated downloads.
///
/// # Adding new presets
/// To add a new TR catalog preset:
/// 1. Add a new variant to this enum with a descriptive doc comment
/// 2. Add its metadata (URL, cache filename) to the `TRPreset::metadata()` method
/// 3. Update tests if needed
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum TRPreset {
    /// STRchive pathogenic disease-associated STRs
    Pathogenic,
    /// ADOTTO TR regions catalog v1.2.1
    Adotto,
    /// Broad Institute TR Explorer catalog (1-1000bp motifs)
    Trexplorer,
    /// CODIS forensic STR markers (USAT catalog)
    Codis,
}

impl TRPreset {
    /// Get metadata for this preset (URL and cache filename)
    pub fn metadata(&self) -> (&'static str, &'static str) {
        match self {
            TRPreset::Pathogenic => (
                "https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/catalogs/STRchive-disease-loci.hg38.longTR.bed",
                "STRchive-disease-loci.hg38.TRGT.bed",
            ),
            TRPreset::Adotto => (
                "https://zenodo.org/records/13987414/files/adotto_TRregions_v1.2.1.bed.gz",
                "adotto_TRregions_v1.2.1.bed.gz",
            ),
            TRPreset::Trexplorer => (
                "https://github.com/broadinstitute/tandem-repeat-catalog/releases/download/v1.0/repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz",
                "repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz",
            ),
            TRPreset::Codis => (
                "https://raw.githubusercontent.com/XuewenWangUGA/USAT/refs/heads/main/settings/STRRegionsV5xwlinuxBest.bed",
                "USAT-CODIS-STRRegionsV5.bed",
            ),
        }
    }

    /// Get a human-readable name for this preset
    pub fn display_name(&self) -> &'static str {
        match self {
            TRPreset::Pathogenic => "STRchive pathogenic STRs",
            TRPreset::Adotto => "ADOTTO TR regions v1.2.1",
            TRPreset::Trexplorer => "TR Explorer catalog (1-1000bp motifs)",
            TRPreset::Codis => "CODIS forensic markers (USAT)",
        }
    }
}

/// Configuration for target selection (what STRs to genotype)
#[derive(Clone, Debug)]
pub struct TargetConfig {
    pub region: Option<String>,
    pub region_file: Option<PathBuf>,
    pub preset: Option<TRPreset>,
    pub max_locus: Option<u32>,
}

impl TargetConfig {
    /// Get target intervals based on the configuration
    /// Returns both the repeat intervals and the chromosome mapper for ID->name lookups
    pub fn get_targets(
        &self,
        bam: &str,
        reference: &Option<String>,
    ) -> InquiSTRResult<(Vec<RepeatInterval>, ChromosomeMapper)> {
        get_targets(self.clone(), bam, reference)
    }
}

#[derive(Debug)]
pub struct RepeatIntervalIterator {
    current_index: usize,
    data: Vec<RepeatInterval>,
    pub num_intervals: usize,
    pub chrom_mapper: ChromosomeMapper,
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

        let chrom = parts[0];
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

        let mut chrom_mapper = ChromosomeMapper::new();
        let chrom_id = chrom_mapper.get_or_insert(chrom);
        let repeat = RepeatInterval::new_interval(
            chrom_id,
            chrom,
            start,
            end,
            ".".to_string(),
            &chrom_lengths,
        )
        .expect("Failed to create repeat interval");
        RepeatIntervalIterator {
            current_index: 0,
            data: vec![repeat],
            num_intervals: 1,
            chrom_mapper,
        }
    }
    pub fn from_bed(
        region_file: &str,
        chrom_lengths: HashMap<String, u64>,
        max_locus: Option<u32>,
    ) -> Self {
        use std::io::{BufRead, BufReader};

        // Use utils::reader to handle gzipped files
        let file_reader = crate::utils::reader(region_file);
        let buf_reader = BufReader::new(file_reader);

        // Check only the first non-empty, non-comment line to see if it's a header
        // This avoids overhead on large BED files
        let lines: Vec<String> = buf_reader
            .lines()
            .map(|l| l.expect("Error reading line from BED file"))
            .collect();

        let mut skipped_headers = 0;
        let mut start_idx = 0;

        // Find first non-empty, non-comment line
        for (idx, line) in lines.iter().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            // Check if this looks like a header (case-insensitive)
            let first_field = trimmed.split('\t').next().unwrap_or("").to_lowercase();
            let is_header = first_field.contains("chrom")
                || first_field.contains("chr") && !first_field.starts_with("chr")
                || first_field == "name"
                || first_field == "id";

            if is_header {
                skipped_headers = 1;
                start_idx = idx + 1;
            }
            break;
        }

        // Remove skipped lines
        let filtered_lines: Vec<String> = lines
            .into_iter()
            .skip(start_idx)
            .filter(|line| {
                let trimmed = line.trim();
                !trimmed.is_empty() && !trimmed.starts_with('#')
            })
            .collect();

        // Join lines and parse as BED
        let filtered_content = filtered_lines.join("\n");
        let cursor = std::io::Cursor::new(filtered_content.as_bytes());
        let mut reader = bed::Reader::new(cursor);

        let mut data = Vec::new();
        let mut filtered_count = 0;
        let mut chrom_mapper = ChromosomeMapper::new();

        for record in reader.records() {
            let rec =
                record.expect("Error reading bed record. Is the file valid and tab-delimited?");

            let repeat = RepeatInterval::from_bed(&rec, &chrom_lengths, &mut chrom_mapper);
            if let Some(repeat) = repeat {
                // Filter by max_locus size if specified
                let locus_size = repeat.end - repeat.start;
                if let Some(max_size) = max_locus
                    && locus_size > max_size
                {
                    filtered_count += 1;
                    continue;
                }
                data.push(repeat);
            }
        }
        if skipped_headers > 0 {
            eprintln!("INFO: Skipped {} header line(s) in BED file", skipped_headers);
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
            data: data.clone(),
            num_intervals: data.len(),
            chrom_mapper,
        }
    }

    /// Download and cache a predefined TR catalog preset
    ///
    /// This method handles downloading TR catalogs from remote URLs, caching them
    /// locally for 7 days, and handling network failures gracefully by falling back
    /// to cached versions when available.
    pub fn from_preset(
        preset: TRPreset,
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
                                                eprintln!(
                                                    "ERROR: Failed to decompress gzipped catalog: {}",
                                                    e
                                                );
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
                                    eprintln!(
                                        "Warning: Failed to read response body ({}), using cached version",
                                        e
                                    );
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
        // Check only the first non-empty, non-comment line to see if it's a header
        let lines: Vec<&str> = data.lines().collect();
        let mut skipped_headers = 0;
        let mut start_idx = 0;

        // Find first non-empty, non-comment line
        for (idx, line) in lines.iter().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            // Check if this looks like a header (case-insensitive)
            let first_field = trimmed.split('\t').next().unwrap_or("").to_lowercase();
            let is_header = first_field.contains("chrom")
                || first_field.contains("chr") && !first_field.starts_with("chr")
                || first_field == "name"
                || first_field == "id";

            if is_header {
                skipped_headers = 1;
                start_idx = idx + 1;
            }
            break;
        }

        // Filter lines: skip header and empty/comment lines
        let filtered_lines: Vec<&str> = lines
            .into_iter()
            .skip(start_idx)
            .filter(|line| {
                let trimmed = line.trim();
                !trimmed.is_empty() && !trimmed.starts_with('#')
            })
            .collect();

        // Join lines and parse as BED
        let filtered_content = filtered_lines.join("\n");
        let cursor = std::io::Cursor::new(filtered_content.as_bytes());
        let mut bed_reader = bed::Reader::new(cursor);

        let mut data_vec = Vec::new();
        let mut filtered_count = 0;
        let mut chrom_mapper = ChromosomeMapper::new();

        for record in bed_reader.records() {
            let rec = record.expect("Error reading bed record from downloaded data");

            let repeat = RepeatInterval::from_bed(&rec, &chrom_lengths, &mut chrom_mapper);
            if let Some(repeat) = repeat {
                // Filter by max_locus size if specified
                let locus_size = repeat.end - repeat.start;
                if let Some(max_size) = max_locus
                    && locus_size > max_size
                {
                    filtered_count += 1;
                    continue;
                }
                data_vec.push(repeat);
            }
        }

        if skipped_headers > 0 {
            eprintln!("INFO: Skipped {} header line(s) in downloaded catalog", skipped_headers);
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
            chrom_mapper,
        }
    }
}

impl Clone for RepeatInterval {
    fn clone(&self) -> Self {
        RepeatInterval {
            chrom_id: self.chrom_id,
            start: self.start,
            end: self.end,
            info: self.info.clone(),
        }
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
    pub chrom_id: u32,
    pub start: u32,
    pub end: u32,
    pub info: String,
}

impl RepeatInterval {
    /// Get chromosome name using the mapper
    pub fn chrom_name<'a>(&self, mapper: &'a ChromosomeMapper) -> &'a str {
        mapper.get_name(self.chrom_id)
    }
}

impl RepeatInterval {
    // parse a bed record
    pub fn from_bed(
        rec: &bed::Record,
        chrom_lengths: &HashMap<String, u64>,
        chrom_mapper: &mut ChromosomeMapper,
    ) -> Option<Self> {
        let chrom = rec.chrom();
        let chrom_id = chrom_mapper.get_or_insert(chrom);
        let start = rec.start().try_into().unwrap();
        let end = rec.end().try_into().unwrap();
        // Extract 4th column (name field) or use "." if not present
        let info = rec
            .name()
            .map(|s| s.to_string())
            .unwrap_or_else(|| ".".to_string());
        RepeatInterval::new_interval(chrom_id, chrom, start, end, info, chrom_lengths)
    }

    fn new_interval(
        chrom_id: u32,
        chrom: &str,
        start: u32,
        end: u32,
        info: String,
        chrom_lengths: &HashMap<String, u64>,
    ) -> Option<Self> {
        if end < start {
            eprintln!(
                "ERROR: Invalid coordinates. End position ({}) is less than start position ({}) for {}:{}-{}",
                end, start, chrom, start, end
            );
            std::process::exit(1);
        }

        // check if the chromosome exists in the chrom lengths hashmap
        // and if the end coordinate is within the chromosome length
        if chrom_lengths.contains_key(chrom) && (end as u64) < chrom_lengths[chrom] {
            return Some(Self { chrom_id, start, end, info });
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
    pub fn new(chrom_id: u32, start: u32, end: u32) -> Self {
        Self { chrom_id, start, end, info: ".".to_string() }
    }
}

/// Get targets from region string, region file, or preset catalog
/// Returns both the repeat intervals and the chromosome mapper for ID->name lookups
pub fn get_targets(
    targets: TargetConfig,
    bam: &str,
    reference: &Option<String>,
) -> InquiSTRResult<(Vec<RepeatInterval>, ChromosomeMapper)> {
    let chrom_lengths = get_chrom_lengths_from_bam_header(bam.to_string(), reference)?;
    let iterator = match (&targets.region, &targets.region_file, &targets.preset) {
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
            return Err(InquiSTRError::new(
                "Specify a region string (-r), a region_file (-R), or --preset <pathogenic|adotto|trexplorer>".to_string()
            ))
        }
    };
    let chrom_mapper = iterator.chrom_mapper.clone();
    Ok((iterator.collect(), chrom_mapper))
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
        let test_bed_content = "chr7\t154778571\t154779363\tsmall_interval\nchr7\t154780000\t154900000\thuge_interval\n";
        let mut file = File::create("test_temp_max_locus.bed").expect("Could not create test file");
        file.write_all(test_bed_content.as_bytes())
            .expect("Could not write test file");

        let chrom_lengths =
            get_chrom_lengths_from_bam_header(String::from("test-data/small-test.bam"), &None)
                .expect("Failed to get chromosome lengths");

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

    /// Test that all preset catalog URLs are accessible
    /// This test requires network access and is ignored by default
    /// Run with: cargo test test_preset_urls -- --ignored
    #[test]
    #[ignore]
    fn test_preset_urls_accessible() {
        use super::TRPreset;

        // Request only the first 1KB to minimize data transfer while still validating URL accessibility
        const TEST_RANGE_BYTES: usize = 1023;

        let presets = vec![TRPreset::Pathogenic, TRPreset::Adotto, TRPreset::Trexplorer];

        for preset in presets {
            let (url, _cache_filename) = preset.metadata();
            let preset_name = preset.display_name();

            eprintln!("Testing URL for {}: {}", preset_name, url);

            // Try to make a HEAD request first (faster and more polite)
            // Add user agent following RFC 7231 format to help identify requests and reduce blocking
            let client = reqwest::blocking::Client::builder()
                .timeout(std::time::Duration::from_secs(30))
                .user_agent("inquiSTR-test/1.0 (+https://github.com/wdecoster/inquiSTR)")
                .build()
                .expect("Failed to build HTTP client");

            let response = client.head(url).send().unwrap_or_else(|e| {
                panic!("Failed to connect to {} ({}): {}", preset_name, url, e)
            });

            let status = response.status();

            // Some servers (like Zenodo) may block HEAD requests with 403/405
            // but still allow GET requests. Try GET if HEAD fails with these codes.
            let is_accessible = if status.is_success() {
                eprintln!("✓ {} URL is accessible via HEAD (status: {})", preset_name, status);
                true
            } else if status == reqwest::StatusCode::FORBIDDEN
                || status == reqwest::StatusCode::METHOD_NOT_ALLOWED
            {
                eprintln!("  HEAD request returned {}, trying GET request...", status);
                // Try a GET request with a range header to minimize data transfer
                match client
                    .get(url)
                    .header("Range", format!("bytes=0-{}", TEST_RANGE_BYTES))
                    .send()
                {
                    Ok(get_response) => {
                        let get_status = get_response.status();
                        // Accept 200 (full response), 206 (Partial Content), or 416 (Range Not Satisfiable)
                        // Status 416 occurs when the file is smaller than our requested range, which still
                        // indicates the URL is accessible and the file exists
                        if get_status.is_success()
                            || get_status == reqwest::StatusCode::RANGE_NOT_SATISFIABLE
                        {
                            eprintln!(
                                "✓ {} URL is accessible via GET (status: {})",
                                preset_name, get_status
                            );
                            true
                        } else {
                            eprintln!(
                                "✗ {} URL GET request returned status: {}",
                                preset_name, get_status
                            );
                            false
                        }
                    }
                    Err(e) => {
                        eprintln!("✗ {} URL GET request failed: {}", preset_name, e);
                        false
                    }
                }
            } else {
                eprintln!("✗ {} URL returned status: {}", preset_name, status);
                false
            };

            assert!(is_accessible, "{} URL ({}) is not accessible", preset_name, url);
        }
    }

    /// Test that preset catalog URLs return valid content
    /// This is a more thorough test that actually downloads the full files
    /// Run with: cargo test test_preset_urls_content -- --ignored
    /// Note: ADOTTO and TRexplorer are large files (90MB+ and 38MB+) and may take time to download
    #[test]
    #[ignore]
    fn test_preset_urls_content() {
        use super::TRPreset;

        let presets = vec![TRPreset::Pathogenic, TRPreset::Adotto, TRPreset::Trexplorer];
        let mut failed_presets = Vec::new();

        for preset in presets {
            let (url, cache_filename) = preset.metadata();
            let preset_name = preset.display_name();
            let is_gzipped = cache_filename.ends_with(".gz");

            eprintln!("Testing content for {}: {}", preset_name, url);

            let client = reqwest::blocking::Client::builder()
                .timeout(std::time::Duration::from_secs(180)) // Increased to 3 minutes for large files
                .connect_timeout(std::time::Duration::from_secs(30))
                .build()
                .expect("Failed to build HTTP client");

            let response = match client.get(url).send() {
                Ok(r) => r,
                Err(e) => {
                    let error_msg = format!("Failed to download {} ({}): {}", preset_name, url, e);
                    eprintln!("✗ {}", error_msg);
                    failed_presets.push((preset_name.to_string(), error_msg));
                    continue;
                }
            };

            if !response.status().is_success() {
                let error_msg =
                    format!("{} URL ({}) returned status: {}", preset_name, url, response.status());
                eprintln!("✗ {}", error_msg);
                failed_presets.push((preset_name.to_string(), error_msg));
                continue;
            }

            let bytes = match response.bytes() {
                Ok(b) => b,
                Err(e) => {
                    let error_msg =
                        format!("Failed to read response for {} ({}): {}", preset_name, url, e);
                    eprintln!("✗ {}", error_msg);
                    failed_presets.push((preset_name.to_string(), error_msg));
                    continue;
                }
            };

            if bytes.is_empty() {
                let error_msg = format!("{} returned empty content", preset_name);
                eprintln!("✗ {}", error_msg);
                failed_presets.push((preset_name.to_string(), error_msg));
                continue;
            }

            // Verify content type based on file format
            if is_gzipped {
                // Check for gzip magic bytes (0x1f 0x8b)
                if bytes.len() < 2 || bytes[0] != 0x1f || bytes[1] != 0x8b {
                    let error_msg =
                        format!("{} does not appear to be a valid gzip file", preset_name);
                    eprintln!("✗ {}", error_msg);
                    failed_presets.push((preset_name.to_string(), error_msg));
                    continue;
                }
                eprintln!("✓ {} content is valid gzip ({} bytes)", preset_name, bytes.len());
            } else {
                // For plain BED files, check if it looks like valid BED format
                let content = String::from_utf8_lossy(&bytes[..std::cmp::min(1000, bytes.len())]);
                // BED files should have tab-separated columns with chromosome names
                let is_valid_bed = content.lines().any(|line| {
                    let fields: Vec<&str> = line.split('\t').collect();
                    fields.len() >= 3 && !line.starts_with('#')
                });
                if !is_valid_bed {
                    let error_msg =
                        format!("{} does not appear to be a valid BED file", preset_name);
                    eprintln!("✗ {}", error_msg);
                    failed_presets.push((preset_name.to_string(), error_msg));
                    continue;
                }
                eprintln!("✓ {} content is valid BED format ({} bytes)", preset_name, bytes.len());
            }
        }

        // Report all failures at once with detailed information
        if !failed_presets.is_empty() {
            let failure_summary: Vec<String> = failed_presets
                .iter()
                .map(|(name, msg)| format!("  - {}: {}", name, msg))
                .collect();
            panic!(
                "\n\n{} preset catalog(s) failed validation:\n{}\n\nFailed presets: {}\n",
                failed_presets.len(),
                failure_summary.join("\n"),
                failed_presets
                    .iter()
                    .map(|(name, _)| name.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
    }
}
