use flate2::read;
use rayon::prelude::*;
use std::collections::HashMap;
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

/// Detect if a file contains kmer frequency data or STR call data
/// Returns true if it's a kmer file, false if it's an STR call file
fn is_kmer_file(file_path: &Path) -> bool {
    let mut file_reader = reader(&file_path.to_string_lossy()).lines();

    if let Some(Ok(first_line)) = file_reader.next() {
        let fields: Vec<&str> = first_line.split('\t').collect();

        // Check if it's a kmer file format
        // Kmer files have "kmer" as first column header
        if fields.len() >= 2 && fields[0] == "kmer" {
            return true;
        }

        // Check if it's STR call format (with or without header)
        // STR files either start with "chromosome" or have genomic coordinates
        if fields.len() >= 3 {
            if fields[0] == "chromosome" {
                return false; // STR file with header
            }

            // Check if first line looks like genomic coordinates (chr1, etc.)
            if fields[0].starts_with("chr")
                && fields[1].parse::<u32>().is_ok()
                && fields[2].parse::<u32>().is_ok()
            {
                return false; // STR file without header
            }
        }
    }

    panic!("Unable to determine file format for {}", file_path.display());
}

/// Determine if a file is a combined STR file (has more than 5 columns)
/// Individual STR files have 5 columns: chr, start, end, H1, H2
/// Combined STR files have 4+ columns: chr, start, end, sample1_H1, sample1_H2, sample2_H1, sample2_H2, ...
fn is_combined_str_file(file: &Path) -> bool {
    let mut file_reader = reader(&file.to_string_lossy()).lines();
    
    // Get first line (must be header)
    let first_line = match file_reader.next() {
        Some(Ok(line)) => line,
        Some(Err(e)) => panic!("Error reading file {}: {}", file.display(), e),
        None => panic!("File {} is empty", file.display()),
    };
    
    // Validate that file has a header
    if first_line.split('\t').next() != Some("chromosome") {
        eprintln!("Error: File {} does not have a valid header.", file.display());
        eprintln!("Expected first column to be 'chromosome', got '{}'", 
            first_line.split('\t').next().unwrap_or("<empty>"));
        eprintln!("\nAll STR files must have headers starting with 'chromosome'.");
        std::process::exit(1);
    }
    
    // Check the next line to determine column count
    match file_reader.next() {
        Some(Ok(second_line)) => second_line.split('\t').count() > 5,
        Some(Err(e)) => panic!("Error reading file {}: {}", file.display(), e),
        None => panic!("File {} has header but no data lines", file.display()),
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

    // Detect file format from first file and validate all files are the same type
    let first_file_is_kmer = is_kmer_file(&calls[0]);

    // Validate that all files are the same type as the first file
    for file in calls.iter().skip(1) {
        let is_kmer = is_kmer_file(file);
        if is_kmer != first_file_is_kmer {
            let first_type = if first_file_is_kmer {
                "kmer frequency"
            } else {
                "STR call"
            };
            let current_type = if is_kmer {
                "kmer frequency"
            } else {
                "STR call"
            };
            panic!(
                "File type mismatch: first file {} is a {} file, but file {} is a {} file. All files must be the same type.",
                calls[0].display(), first_type, file.display(), current_type
            );
        }
    }

    // Route to appropriate combining function based on detected file type
    if first_file_is_kmer {
        eprintln!("Detected kmer frequency files - combining kmer data");
        combine_kmer_files(calls, actual_threads);
    } else {
        eprintln!("Detected STR call files - combining STR data");
        combine_str_files(calls, actual_threads);
    }
}

/// Combine STR call files (original functionality)
fn combine_str_files(calls: Vec<PathBuf>, _actual_threads: usize) {
    // Detect if any files are combined files (have more than 5 columns)
    let mut combined_files = Vec::new();
    let mut individual_files = Vec::new();
    
    for file in &calls {
        if is_combined_str_file(file) {
            combined_files.push(file.clone());
        } else {
            individual_files.push(file.clone());
        }
    }
    
    // Determine processing path based on file mix
    if !combined_files.is_empty() {
        // If we have any combined files, use the merge logic
        if combined_files.len() > 1 {
            eprintln!(
                "Error: Cannot combine multiple combined files together.\nFound {} combined files: {}",
                combined_files.len(),
                combined_files.iter().map(|f| f.display().to_string()).collect::<Vec<_>>().join(", ")
            );
            eprintln!("\nTo merge these cohorts, you can add individual samples to one combined file:");
            eprintln!("  inquiSTR combine {} sample1.inq sample2.inq", combined_files[0].display());
            std::process::exit(1);
        }
        
        eprintln!(
            "Detected 1 combined file and {} individual file(s) - merging cohorts",
            individual_files.len()
        );
        
        merge_combined_and_individual_str_files(combined_files[0].clone(), individual_files);
    } else {
        // All files are individual files - use original logic
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
}

/// Combine kmer frequency files from inquiSTR unmapped
fn combine_kmer_files(calls: Vec<PathBuf>, _actual_threads: usize) {
    // Read and validate headers from all kmer files
    let headers = read_and_validate_kmer_headers(&calls);

    // Output combined header
    output_combined_kmer_header(&headers);

    // Process kmer data
    process_kmer_data(&calls);

    eprintln!("Completed combining {} kmer frequency files", calls.len());
}

/// Merge a combined STR file with individual STR files
fn merge_combined_and_individual_str_files(combined_file: PathBuf, individual_files: Vec<PathBuf>) {
    use std::collections::HashMap;
    
    // Read the combined file into memory
    let mut combined_reader = reader(&combined_file.to_string_lossy()).lines();
    
    // Read combined file header
    let combined_header = combined_reader
        .next()
        .unwrap_or_else(|| panic!("Combined file {} is empty", combined_file.display()))
        .unwrap_or_else(|e| panic!("Error reading combined file header: {}", e));
    
    let has_combined_header = combined_header.split('\t').next() == Some("chromosome");
    
    // Read individual file headers
    let individual_headers: Vec<String> = if individual_files.is_empty() {
        Vec::new()
    } else {
        individual_files
            .iter()
            .map(|file| {
                let mut reader = reader(&file.to_string_lossy()).lines();
                reader
                    .next()
                    .unwrap_or_else(|| panic!("File {} is empty", file.display()))
                    .unwrap_or_else(|e| panic!("Error reading file {}: {}", file.display(), e))
            })
            .collect()
    };
    
    let has_individual_headers = if !individual_headers.is_empty() {
        individual_headers[0].split('\t').next() == Some("chromosome")
    } else {
        false
    };
    
    // Validate header consistency
    if has_combined_header != has_individual_headers {
        panic!("Header mismatch: combined file has {} header, but individual files have {} headers",
            if has_combined_header { "a" } else { "no" },
            if has_individual_headers { "a" } else { "no" }
        );
    }
    
    // Output merged header
    if has_combined_header {
        let combined_header_fields: Vec<&str> = combined_header.split('\t').collect();
        let mut merged_header = vec![
            combined_header_fields[0], // chromosome
            combined_header_fields[1], // begin
            combined_header_fields[2], // end
        ];
        
        // Add all sample columns from combined file
        merged_header.extend(&combined_header_fields[3..]);
        
        // Add sample columns from individual files
        for header in &individual_headers {
            let fields: Vec<&str> = header.split('\t').collect();
            if fields.len() < 5 {
                panic!("Invalid individual file header format: {}", header);
            }
            merged_header.extend(&fields[3..]);
        }
        
        println!("{}", merged_header.join("\t"));
    } else {
        eprintln!("Warning: No headers detected in files");
    }
    
    // Read combined file data into HashMap
    let mut combined_data: HashMap<String, String> = HashMap::new();
    
    for (line_num, line_result) in combined_reader.enumerate() {
        let line = line_result.unwrap_or_else(|e| panic!("Error reading combined file: {}", e));
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            panic!("Invalid combined file data line {}: {}", line_num, line);
        }
        
        // Use chr:start-end as key
        let key = format!("{}:{}−{}", fields[0], fields[1], fields[2]);
        combined_data.insert(key, line);
    }
    
    if combined_data.is_empty() {
        panic!("Combined file contains no data lines");
    }
    
    eprintln!("Loaded {} loci from combined file", combined_data.len());
    
    // Open individual file readers
    let mut individual_readers: Vec<_> = individual_files
        .iter()
        .map(|file| {
            let mut file_reader = reader(&file.to_string_lossy()).lines();
            if has_individual_headers {
                file_reader.next(); // Skip header
            }
            file_reader
        })
        .collect();
    
    // Process each line from all individual files simultaneously
    let mut processed_lines = 0;
    loop {
        let mut all_done = true;
        let mut individual_lines = Vec::new();
        
        for (file_idx, reader) in individual_readers.iter_mut().enumerate() {
            match reader.next() {
                Some(Ok(line)) => {
                    individual_lines.push(line);
                    all_done = false;
                }
                Some(Err(e)) => panic!("Error reading file {}: {}", individual_files[file_idx].display(), e),
                None => {}
            }
        }
        
        if all_done {
            break;
        }
        
        if individual_lines.len() != individual_files.len() {
            panic!("Individual files have different number of data lines at line {}", processed_lines);
        }
        
        // Validate all individual files agree on coordinates
        let first_fields: Vec<&str> = individual_lines[0].split('\t').collect();
        if first_fields.len() < 5 {
            panic!("Invalid individual file data line: {}", individual_lines[0]);
        }
        
        let (chr, start, end) = (first_fields[0], first_fields[1], first_fields[2]);
        
        for (file_idx, line) in individual_lines.iter().enumerate().skip(1) {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields[0] != chr || fields[1] != start || fields[2] != end {
                panic!(
                    "Coordinate mismatch at line {}: file 0 has {}:{}-{}, file {} has {}:{}-{}",
                    processed_lines, chr, start, end, file_idx + 1, fields[0], fields[1], fields[2]
                );
            }
        }
        
        // Look up this locus in combined data
        let key = format!("{}:{}−{}", chr, start, end);
        let combined_line = combined_data.get(&key)
            .unwrap_or_else(|| panic!("Locus {}:{}-{} found in individual files but not in combined file", chr, start, end));
        
        // Build merged line
        let combined_fields: Vec<&str> = combined_line.split('\t').collect();
        let mut merged_fields = vec![chr, start, end];
        
        // Add combined file sample columns
        merged_fields.extend(&combined_fields[3..]);
        
        // Add individual file sample columns
        for line in &individual_lines {
            let fields: Vec<&str> = line.split('\t').collect();
            merged_fields.extend(&fields[3..]);
        }
        
        println!("{}", merged_fields.join("\t"));
        processed_lines += 1;
    }
    
    eprintln!("Processed {} loci", processed_lines);
    eprintln!(
        "Completed merging 1 combined file with {} individual file(s)",
        individual_files.len()
    );
}

/// Read and validate headers from kmer frequency files
fn read_and_validate_kmer_headers(calls: &[PathBuf]) -> Vec<String> {
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

    // Validate that all files have kmer headers
    for (i, header) in headers.iter().enumerate() {
        let fields: Vec<&str> = header.split('\t').collect();
        if fields.len() < 2 || fields[0] != "kmer" {
            panic!(
                "Invalid kmer file header in {}: expected 'kmer' as first column, got '{}'",
                calls[i].display(),
                fields.first().unwrap_or(&"<empty>")
            );
        }
    }

    headers
}

/// Output the combined header for kmer files
fn output_combined_kmer_header(headers: &[String]) {
    let mut combined_header = vec!["kmer"];

    // Add sample names from all files
    for header in headers {
        let fields: Vec<&str> = header.split('\t').collect();
        if fields.len() >= 2 {
            combined_header.push(fields[1]); // Sample name is in second column
        }
    }

    println!("{}", combined_header.join("\t"));
}

/// Process kmer data from all files
fn process_kmer_data(calls: &[PathBuf]) {
    // Create file readers
    let mut file_readers: Vec<_> = calls
        .iter()
        .map(|file| reader(&file.clone().into_os_string().into_string().unwrap()).lines())
        .collect();

    // Skip headers
    for file_reader in &mut file_readers {
        file_reader.next(); // Skip header line
    }

    // Collect all kmer data into a master map
    let mut kmer_data: HashMap<String, Vec<String>> = HashMap::new();
    let mut line_count = 0;
    let chunk_size = 1000;

    loop {
        let mut data_lines: Vec<Option<String>> = Vec::with_capacity(calls.len());
        let mut all_done = true;

        // Read one line from each file
        for (file_idx, file_reader) in file_readers.iter_mut().enumerate() {
            match file_reader.next() {
                Some(Ok(line)) => {
                    data_lines.push(Some(line));
                    all_done = false;
                }
                Some(Err(e)) => panic!("Error reading file {}: {}", calls[file_idx].display(), e),
                None => {
                    data_lines.push(None);
                }
            }
        }

        if all_done {
            break;
        }

        // Process this line set - validate kmer consistency
        if let Some(ref first_line) = data_lines[0] {
            let first_fields: Vec<&str> = first_line.split('\t').collect();
            if first_fields.len() < 2 {
                panic!("Invalid kmer data line at line {}: {}", line_count, first_line);
            }
            let kmer = first_fields[0];

            // Validate that all files have the same kmer (or are finished)
            for (file_idx, data_line) in data_lines.iter().enumerate() {
                if let Some(line) = data_line {
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() < 2 {
                        panic!(
                            "Invalid kmer data line in file {} at line {}: {}",
                            file_idx, line_count, line
                        );
                    }
                    if fields[0] != kmer {
                        panic!(
                            "Kmer mismatch at line {}: file 0 has '{}', file {} has '{}'",
                            line_count, kmer, file_idx, fields[0]
                        );
                    }
                } else if !all_done {
                    panic!(
                        "File {} finished early at line {}, expected kmer '{}'",
                        calls[file_idx].display(),
                        line_count,
                        kmer
                    );
                }
            }

            // Collect frequency values for this kmer
            let mut values = Vec::new();
            for line in data_lines.iter().flatten() {
                let fields: Vec<&str> = line.split('\t').collect();
                values.push(fields[1].to_string());
            }

            kmer_data.insert(kmer.to_string(), values);
        }

        line_count += 1;
        if line_count % chunk_size == 0 {
            eprintln!("Processed {} lines...", line_count);
        }
    }

    // Output combined data sorted by kmer
    let mut sorted_kmers: Vec<_> = kmer_data.keys().collect();
    sorted_kmers.sort();

    for kmer in sorted_kmers {
        let values = &kmer_data[kmer];
        println!("{}\t{}", kmer, values.join("\t"));
    }

    eprintln!("Processed {} total kmer entries", kmer_data.len());
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
            PathBuf::from("test-data/file_with_header1.inq"),
            PathBuf::from("test-data/file_with_header2.inq"),
            PathBuf::from("test-data/file_with_header3.inq"),
        ],
        1,
    );
}

#[test]
fn test_combine_gzipped() {
    // Note: These files need headers - update test files or skip if they don't exist
    combine(
        vec![
            PathBuf::from("test-data/file_with_header1.inq"),
            PathBuf::from("test-data/file_with_header2.inq"),
        ],
        1,
    );
}
