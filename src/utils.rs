use flate2::read::GzDecoder;
use noodles_bgzf as bgzf;
use std::fs::File;
use std::io::{BufReader, Read};

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
    if filename.ends_with(".gz") && is_bgzip_file(filename) {
        // Handle bgzip files using noodles-bgzf
        let file = File::open(filename).expect("Problem opening bgzip file");
        let bgzf_reader = bgzf::Reader::new(file);
        BufReader::new(Box::new(bgzf_reader) as Box<dyn Read>)
    } else if filename.ends_with(".gz") {
        // Handle regular gzip files
        let file = File::open(filename).expect("Problem opening gzip file");
        BufReader::new(Box::new(GzDecoder::new(file)) as Box<dyn Read>)
    } else {
        // Handle uncompressed files
        let file = File::open(filename).expect("Problem opening file");
        BufReader::new(Box::new(file) as Box<dyn Read>)
    }
}

/// parse a region string
pub fn process_region(reg: String) -> Result<(String, u32, u32), Box<dyn std::error::Error>> {
    let reg = reg.replace(',', "");
    if reg.matches(':').count() != 1 {
        eprintln!("ERROR: Invalid region format. Expected format: 'chr:start-end' (e.g., 'chr1:1000-2000')");
        eprintln!("Got: {}", reg);
        std::process::exit(1);
    }
    if reg.matches('-').count() != 1 {
        eprintln!("ERROR: Invalid region format. Could not find exactly one '-' character separating start and end");
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
