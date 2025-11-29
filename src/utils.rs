use flate2::read::GzDecoder;
use noodles_bgzf as bgzf;
use std::fs::File;
use std::io::{BufReader, Read};
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
        BufReader::new(Box::new(GzDecoder::new(file)) as Box<dyn Read>)
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

/// Skip all metadata lines (lines starting with #) and return the actual data/header line
/// This allows for future expansion of metadata types beyond just file_type
pub fn skip_metadata_lines(
    lines: &mut dyn Iterator<Item = Result<String, std::io::Error>>,
) -> String {
    loop {
        let line = lines
            .next()
            .expect("File is empty or contains only metadata")
            .expect("Error reading line");

        if !line.starts_with('#') {
            return line;
        }
        // Otherwise keep looping to skip metadata lines
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
