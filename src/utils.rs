use flate2::read::MultiGzDecoder;
use noodles_bgzf as bgzf;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

/// Extract a sample name from a file path by removing path, extension, and common BAM/CRAM suffixes
pub fn extract_sample_name_from_path(path: &str) -> String {
    let path_buf = PathBuf::from(path);

    // Handle URLs by extracting just the filename part
    let filename = if path.starts_with("http") || path.starts_with("ftp") || path.starts_with("s3")
    {
        path.rsplit('/').next().unwrap_or(path)
    } else {
        path_buf
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or(path)
    };

    // Remove common extensions and suffixes
    let mut result = filename.to_string();

    // Remove BAM/CRAM extensions and index files (check longer suffixes first)
    if result.ends_with(".bam.bai") {
        result = result.strip_suffix(".bam.bai").unwrap().to_string();
    } else if result.ends_with(".cram.crai") {
        result = result.strip_suffix(".cram.crai").unwrap().to_string();
    } else if result.ends_with(".bam") {
        result = result.strip_suffix(".bam").unwrap().to_string();
    } else if result.ends_with(".cram") {
        result = result.strip_suffix(".cram").unwrap().to_string();
    }

    result
}

/// Detect if a file starts with gzip magic bytes (0x1f 0x8b)
fn is_gzip_file(filename: &str) -> bool {
    if let Ok(mut file) = File::open(filename) {
        let mut magic = [0u8; 2];
        if file.read_exact(&mut magic).is_ok() {
            return magic[0] == 0x1f && magic[1] == 0x8b;
        }
    }
    false
}

/// Detect if a file is bgzip compressed by checking magic bytes
fn is_bgzip_file(filename: &str) -> bool {
    if !filename.ends_with(".gz") {
        return false;
    }

    if let Ok(mut file) = File::open(filename) {
        let mut magic = [0u8; 16];
        if file.read_exact(&mut magic).is_ok() {
            // Bgzip has specific magic bytes: 1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00
            // Check for bgzip signature (simplified check)
            return magic[0] == 0x1f && magic[1] == 0x8b && magic[2] == 0x08 && magic[3] == 0x04;
        }
    }
    false
}

/// Read normal, gzip, or bgzip compressed files seamlessly
/// Uses file extension and magic bytes to decide decompression method
pub fn reader(filename: &str) -> BufReader<Box<dyn Read>> {
    let path = Path::new(filename);

    if filename.ends_with(".gz") && is_bgzip_file(filename) {
        // Handle bgzip files using noodles-bgzf
        let file = File::open(path).unwrap_or_else(|e| {
            eprintln!("Error: Failed to open bgzip file {}", path.display());
            eprintln!("Reason: {}", e);
            eprintln!("\nPlease check that the file exists and you have permission to read it.");
            std::process::exit(1);
        });
        let bgzf_reader = bgzf::io::Reader::new(file);
        BufReader::new(Box::new(bgzf_reader) as Box<dyn Read>)
    } else if filename.ends_with(".gz") {
        // Handle regular gzip files
        let file = File::open(path).unwrap_or_else(|e| {
            eprintln!("Error: Failed to open gzip file {}", path.display());
            eprintln!("Reason: {}", e);
            eprintln!("\nPlease check that the file exists and you have permission to read it.");
            std::process::exit(1);
        });
        BufReader::new(Box::new(MultiGzDecoder::new(file)) as Box<dyn Read>)
    } else if is_gzip_file(filename) {
        // Handle gzip files without .gz extension
        let file = File::open(path).unwrap_or_else(|e| {
            eprintln!("Error: Failed to open file {}", path.display());
            eprintln!("Reason: {}", e);
            eprintln!("\nPlease check that the file exists and you have permission to read it.");
            std::process::exit(1);
        });
        BufReader::new(Box::new(MultiGzDecoder::new(file)) as Box<dyn Read>)
    } else {
        // Handle uncompressed files
        let file = File::open(path).unwrap_or_else(|e| {
            eprintln!("Error: Failed to open file {}", path.display());
            eprintln!("Reason: {}", e);
            eprintln!("\nPlease check that the file exists and you have permission to read it.");
            std::process::exit(1);
        });
        BufReader::new(Box::new(file) as Box<dyn Read>)
    }
}

/// Parse sample input - can be a file path, comma-separated names, or a single name
pub fn parse_sample_input(input: &str) -> Vec<String> {
    let path = Path::new(input);

    // Check if input is a file path
    if path.exists() && path.is_file() {
        eprintln!("Reading sample names from file: {}", input);
        let file = reader(input);
        file.lines()
            .map(|line| line.unwrap().trim().to_string())
            .filter(|line| !line.is_empty())
            .collect()
    } else if input.contains(',') {
        // Comma-separated sample names
        eprintln!("Parsing comma-separated sample names");
        input
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    } else {
        // Single sample name
        eprintln!("Using single sample name: {}", input);
        vec![input.to_string()]
    }
}

/// Skip all metadata lines (lines starting with #) and return the actual data/header line
/// This allows for future expansion of metadata types beyond just file_type
pub fn skip_metadata_lines(
    lines: &mut dyn Iterator<Item = Result<String, std::io::Error>>,
    filename: &str,
) -> String {
    loop {
        match lines.next() {
            Some(Ok(line)) => {
                if !line.starts_with('#') {
                    return line;
                }
                // Otherwise keep looping to skip metadata lines
            }
            Some(Err(e)) => {
                eprintln!("Error: Failed to read file '{}': {}", filename, e);
                eprintln!("The file may be truncated or corrupt.");
                std::process::exit(1);
            }
            None => {
                eprintln!("Error: File '{}' is empty or contains only metadata lines.", filename);
                std::process::exit(1);
            }
        }
    }
}

/// parse a region string
pub fn process_region(reg: String) -> Result<(String, u32, u32), Box<dyn std::error::Error>> {
    let reg = reg.replace(',', "");
    if reg.matches(':').count() != 1 {
        eprintln!(
            "ERROR: Invalid region format. Expected format: 'chr:start-end' (e.g., 'chr1:1000-2000')"
        );
        eprintln!("Got: {}", reg);
        std::process::exit(1);
    }
    if reg.matches('-').count() != 1 {
        eprintln!(
            "ERROR: Invalid region format. Could not find exactly one '-' character separating start and end"
        );
        eprintln!("Expected format: 'chr:start-end' (e.g., 'chr1:1000-2000')");
        eprintln!("Got: {}", reg);
        std::process::exit(1);
    }
    let chrom = reg.split(':').collect::<Vec<&str>>()[0];
    let interval = reg.split(':').collect::<Vec<&str>>()[1];
    let start: u32 = interval.split('-').collect::<Vec<&str>>()[0]
        .parse()
        .expect("\n\nError while parsing interval start coordinate!\n\n");
    let end: u32 = interval.split('-').collect::<Vec<&str>>()[1]
        .parse()
        .expect("\n\nError while parsing interval end coordinate!\n\n");
    assert!(
        start < end,
        r#"\n\nInvalid region: start coordinate has to be smaller than end.\n\n"#
    );
    Ok((chrom.to_string(), start, end))
}

#[cfg(test)]
#[test]
fn test_extract_sample_name_from_path() {
    // Test local file paths
    assert_eq!(extract_sample_name_from_path("test-data/sample.bam"), "sample");
    assert_eq!(extract_sample_name_from_path("test-data/sample.cram"), "sample");
    assert_eq!(extract_sample_name_from_path("/path/to/HG00096.hg38.cram"), "HG00096.hg38");

    // Test URL paths
    assert_eq!(extract_sample_name_from_path("https://example.com/data/sample.bam"), "sample");
    assert_eq!(
        extract_sample_name_from_path("s3://bucket/data/HG00096.hg38.cram"),
        "HG00096.hg38"
    );

    // Test multiple extensions
    assert_eq!(extract_sample_name_from_path("sample.sorted.dedup.bam"), "sample.sorted.dedup");

    // Test edge cases
    assert_eq!(extract_sample_name_from_path("sample"), "sample");
    assert_eq!(extract_sample_name_from_path("sample.bam.bai"), "sample");
}
