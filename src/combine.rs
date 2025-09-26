use flate2::read;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::path::PathBuf;

/// Read normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(128 * 1024, read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}

pub fn combine(calls: Vec<PathBuf>) {
    // check if all files exist
    for file in &calls {
        if !file.exists() {
            panic!("File {} does not exist!", file.display());
        }
    }
    
    // Open all files
    let mut file_readers: Vec<_> = calls
        .iter()
        .map(|file| reader(&file.clone().into_os_string().into_string().unwrap()).lines())
        .collect();
    
    // Read headers from all files
    let mut headers: Vec<String> = Vec::new();
    for (i, file_reader) in file_readers.iter_mut().enumerate() {
        let header = file_reader.next()
            .unwrap_or_else(|| panic!("File {} is empty", calls[i].display()))
            .unwrap_or_else(|e| panic!("Error reading header from {}: {}", calls[i].display(), e));
        headers.push(header);
    }
    
    // Check if files have headers (look for "chromosome" in first field)
    let has_headers = headers[0].split('\t').next() == Some("chromosome");
    
    if has_headers {
        // Construct combined header
        let first_header_fields: Vec<&str> = headers[0].split('\t').collect();
        if first_header_fields.len() < 5 {
            panic!("Invalid header format in first file: {}", headers[0]);
        }
        
        // Start with chr, begin, end
        let mut combined_header = vec![first_header_fields[0], first_header_fields[1], first_header_fields[2]];
        
        // Add sample columns from first file
        combined_header.extend(&first_header_fields[3..]);
        
        // Add sample columns from other files (skip chr, begin, end)
        for header in &headers[1..] {
            let fields: Vec<&str> = header.split('\t').collect();
            if fields.len() < 5 {
                panic!("Invalid header format: {}", header);
            }
            combined_header.extend(&fields[3..]);
        }
        
        println!("{}", combined_header.join("\t"));
    }
    
    // Process data lines
    loop {
        let mut data_lines: Vec<String> = Vec::new();
        let mut all_done = true;
        
        // Read one line from each file
        for file_reader in &mut file_readers {
            match file_reader.next() {
                Some(Ok(line)) => {
                    data_lines.push(line);
                    all_done = false;
                }
                Some(Err(e)) => panic!("Error reading file: {}", e),
                None => {
                    // File is finished - but we need all files to have same number of lines
                    if !all_done {
                        panic!("Files have different number of data lines");
                    }
                }
            }
        }
        
        if all_done {
            break;
        }
        
        // Construct combined line
        let first_line_fields: Vec<&str> = data_lines[0].split('\t').collect();
        if first_line_fields.len() < 3 {
            panic!("Invalid data line format: {}", data_lines[0]);
        }
        
        // Start with chr, begin, end from first file
        let mut combined_line = vec![first_line_fields[0], first_line_fields[1], first_line_fields[2]];
        
        // Add all sample data from first file
        combined_line.extend(&first_line_fields[3..]);
        
        // Add sample data from other files (skip chr, begin, end)
        for line in &data_lines[1..] {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                panic!("Invalid data line format: {}", line);
            }
            combined_line.extend(&fields[3..]);
        }
        
        println!("{}", combined_line.join("\t"));
    }
}

#[cfg(test)]
#[test]
fn test_combine() {
    combine(vec![
        PathBuf::from("test-data/file1.inq"),
        PathBuf::from("test-data/file2.inq"),
        PathBuf::from("test-data/file3.inq"),
    ]);
}

#[test]
fn test_combine_gzipped() {
    combine(vec![
        PathBuf::from("test-data/file1.inq.gz"),
        PathBuf::from("test-data/file2.inq.gz"),
        PathBuf::from("test-data/file3.inq.gz"),
    ]);
}
