use flate2::read;
use rayon::prelude::*;
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

pub fn combine(calls: Vec<PathBuf>, threads: usize) {
    // Configure thread pool
    let actual_threads = if threads == 0 {
        rayon::current_num_threads()
    } else if rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .is_ok()
    {
        threads
    } else {
        eprintln!("Warning: Thread pool already configured, using existing settings");
        rayon::current_num_threads()
    };

    eprintln!(
        "Using {} threads for parallel processing of {} files",
        actual_threads,
        calls.len()
    );

    // Check if all files exist
    for file in &calls {
        if !file.exists() {
            panic!("File {} does not exist!", file.display());
        }
    }

    // Read and validate headers first
    let headers = read_and_validate_headers(&calls);
    let has_headers = headers[0].split('\t').next() == Some("chromosome");

    // Output combined header if present
    if has_headers {
        output_combined_header(&headers);
    } else {
        eprintln!("Warning: No headers detected in input files");
    }

    // Process data with parallel chunked reading
    process_data_parallel(&calls, has_headers);

    eprintln!("Completed combining {} files", calls.len());
}

/// Read and validate headers from all files
fn read_and_validate_headers(calls: &[PathBuf]) -> Vec<String> {
    let headers: Vec<String> = calls
        .par_iter()
        .enumerate()
        .map(|(_i, file)| {
            let mut file_reader =
                reader(&file.clone().into_os_string().into_string().unwrap()).lines();
            file_reader
                .next()
                .unwrap_or_else(|| panic!("File {} is empty", file.display()))
                .unwrap_or_else(|e| panic!("Error reading header from {}: {}", file.display(), e))
        })
        .collect();

    // Validate that all files have headers or none do
    let first_has_header = headers[0].split('\t').next() == Some("chromosome");
    for (i, header) in headers.iter().enumerate() {
        let has_header = header.split('\t').next() == Some("chromosome");
        if has_header != first_has_header {
            panic!(
                "Inconsistent header presence: file {} differs from first file",
                calls[i].display()
            );
        }
    }

    headers
}

/// Output the combined header line
fn output_combined_header(headers: &[String]) {
    let first_header_fields: Vec<&str> = headers[0].split('\t').collect();
    if first_header_fields.len() < 5 {
        panic!("Invalid header format in first file: {}", headers[0]);
    }

    // Start with chr, begin, end
    let mut combined_header = vec![
        first_header_fields[0],
        first_header_fields[1],
        first_header_fields[2],
    ];

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

/// Process data lines with optimized sequential reading and parallel line processing
fn process_data_parallel(calls: &[PathBuf], has_headers: bool) {
    // Create file readers
    let mut file_readers: Vec<_> = calls
        .iter()
        .map(|file| reader(&file.clone().into_os_string().into_string().unwrap()).lines())
        .collect();

    // Skip headers if present
    if has_headers {
        for file_reader in &mut file_readers {
            file_reader.next(); // Skip header line
        }
    }

    let mut line_count = 0;
    let chunk_size = 1000; // Process 1000 lines at a time for progress reporting

    loop {
        let mut data_lines: Vec<String> = Vec::with_capacity(calls.len());
        let mut all_done = true;

        // Read one line from each file
        for (file_idx, file_reader) in file_readers.iter_mut().enumerate() {
            match file_reader.next() {
                Some(Ok(line)) => {
                    data_lines.push(line);
                    all_done = false;
                }
                Some(Err(e)) => panic!("Error reading file {}: {}", calls[file_idx].display(), e),
                None => {
                    // File is finished - check if all files are done
                    if !all_done {
                        panic!("Files have different number of data lines at line {}", line_count);
                    }
                }
            }
        }

        if all_done {
            break;
        }

        // Validate variant consistency and combine lines in parallel if we have enough lines
        let combined_line = combine_data_lines(&data_lines, line_count);
        println!("{}", combined_line);

        line_count += 1;
        if line_count % chunk_size == 0 {
            eprintln!("Processed {} lines...", line_count);
        }
    }

    eprintln!("Processed {} total lines", line_count);
}

/// Combine data lines from all files with variant consistency checking
fn combine_data_lines(data_lines: &[String], line_number: usize) -> String {
    if data_lines.is_empty() {
        panic!("No data lines to combine at line {}", line_number);
    }

    // Parse first line to get coordinates
    let first_line_fields: Vec<&str> = data_lines[0].split('\t').collect();
    if first_line_fields.len() < 3 {
        panic!("Invalid data line format at line {}: {}", line_number, data_lines[0]);
    }

    let (chr, start, end) = (first_line_fields[0], first_line_fields[1], first_line_fields[2]);

    // Validate that all files have the same variant coordinates
    for (file_idx, line) in data_lines.iter().enumerate().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            panic!(
                "Invalid data line format in file {} at line {}: {}",
                file_idx, line_number, line
            );
        }

        if fields[0] != chr || fields[1] != start || fields[2] != end {
            panic!(
                "Variant mismatch at line {}: file 0 has {}:{}-{}, file {} has {}:{}-{}",
                line_number, chr, start, end, file_idx, fields[0], fields[1], fields[2]
            );
        }
    }

    // Pre-allocate combined line capacity to reduce allocations
    let estimated_capacity =
        data_lines.iter().map(|line| line.len()).sum::<usize>() + data_lines.len() * 10;
    let mut combined_fields = Vec::with_capacity(estimated_capacity / 8); // Rough estimate of field count

    // Add coordinates from first file
    combined_fields.push(chr);
    combined_fields.push(start);
    combined_fields.push(end);

    // Add all sample data from first file
    combined_fields.extend(&first_line_fields[3..]);

    // Add sample data from other files (skip coordinates)
    for line in &data_lines[1..] {
        let fields: Vec<&str> = line.split('\t').collect();
        combined_fields.extend(&fields[3..]);
    }

    combined_fields.join("\t")
}

#[cfg(test)]
#[test]
fn test_combine() {
    combine(
        vec![
            PathBuf::from("test-data/file1.inq"),
            PathBuf::from("test-data/file2.inq"),
            PathBuf::from("test-data/file3.inq"),
        ],
        1,
    );
}

#[test]
fn test_combine_gzipped() {
    combine(
        vec![
            PathBuf::from("test-data/file1.inq.gz"),
            PathBuf::from("test-data/file2.inq.gz"),
            PathBuf::from("test-data/file3.inq.gz"),
        ],
        1,
    );
}
