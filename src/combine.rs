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
        Err(why) => {
            eprintln!("Error: Failed to open file {}", path.display());
            eprintln!("Reason: {}", why);
            eprintln!("\nPlease check that the file exists and you have permission to read it.");
            std::process::exit(1);
        }
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
            eprintln!("Error: File does not exist: {}", file.display());
            eprintln!("\nPlease check the file path and try again.");
            std::process::exit(1);
        }
    }

    // Detect file format from first file and validate all files are the same type
    let first_file_is_kmer = crate::filetype::is_kmer_file(&calls[0]);

    // Validate that all files are the same type as the first file
    for file in calls.iter().skip(1) {
        let is_kmer = crate::filetype::is_kmer_file(file);
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
            eprintln!("Error: File type mismatch detected.");
            eprintln!("  First file: {} ({})", calls[0].display(), first_type);
            eprintln!("  This file: {} ({})", file.display(), current_type);
            eprintln!(
                "\nAll files must be the same type. Cannot combine kmer and STR files together."
            );
            std::process::exit(1);
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
        if crate::filetype::is_combined_str_file(file) {
            combined_files.push(file.clone());
        } else {
            individual_files.push(file.clone());
        }
    }

    // Determine processing path based on file mix
    if !combined_files.is_empty() {
        // If we have any combined files, use the merge logic
        eprintln!(
            "Detected {} combined file(s) and {} individual file(s) - merging cohorts",
            combined_files.len(),
            individual_files.len()
        );

        merge_combined_files(combined_files, individual_files);
    } else {
        // All files are individual files - use original logic
        // Read and validate headers first
        let headers = read_and_validate_headers(&calls);
        let has_headers = headers[0].split('\t').next() == Some("chromosome");

        // Validate that headers are present
        if !has_headers {
            eprintln!("Error: No headers detected in input files.");
            eprintln!("All inquiSTR output files should have headers starting with 'chromosome'.");
            eprintln!("\nPlease ensure you are combining valid inquiSTR output files.");
            std::process::exit(1);
        }

        // Output combined header
        output_combined_header(&headers);

        // Process data with parallel chunked reading
        process_data_parallel(&calls);

        eprintln!("Completed combining {} files", calls.len());
    }
}

/// Combine kmer frequency files from inquiSTR unmapped
fn combine_kmer_files(calls: Vec<PathBuf>, _actual_threads: usize) {
    // Read and validate headers from all kmer files
    let headers = read_and_validate_kmer_headers(&calls);

    // Check if these are target kmer files or regular kmer files
    let first_header_fields: Vec<&str> = headers[0].split('\t').collect();
    let is_target_kmer = first_header_fields[0] == "Sample";

    if is_target_kmer {
        // Combine target kmer files (each file has one row per sample)
        combine_target_kmer_files(&calls, &headers);
    } else {
        // Output combined header for regular kmer files
        output_combined_kmer_header(&headers);

        // Process kmer data for regular kmer files
        process_kmer_data(&calls);
    }

    eprintln!("Completed combining {} kmer frequency files", calls.len());
}

/// Combine target kmer files (each file typically has one data row per sample)
fn combine_target_kmer_files(calls: &[PathBuf], headers: &[String]) {
    // Parse all files and collect data
    let mut all_data: Vec<HashMap<String, String>> = Vec::new();

    for file in calls {
        let mut file_reader = reader(&file.to_string_lossy()).lines();

        // Skip metadata lines and header
        while let Some(Ok(line)) = file_reader.next() {
            if line.starts_with('#') {
                continue; // Skip metadata lines
            } else {
                // This is the header line, consume it and stop
                break;
            }
        }

        let mut file_data = HashMap::new();
        for line in file_reader.map_while(Result::ok) {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 2 {
                let sample_name = fields[0].to_string();
                // Store entire row for this sample
                file_data.insert(sample_name, line);
            }
        }
        all_data.push(file_data);
    }

    // Output header (use the first file's header as template)
    println!("{}", headers[0]);

    // Get all unique sample names across all files
    let mut all_samples: std::collections::HashSet<String> = std::collections::HashSet::new();
    for data in &all_data {
        for sample in data.keys() {
            all_samples.insert(sample.clone());
        }
    }

    // Output data rows
    for sample in all_samples.iter() {
        // Find the data for this sample from the first file that has it
        for data in &all_data {
            if let Some(row) = data.get(sample) {
                println!("{}", row);
                break; // Only output each sample once
            }
        }
    }
}

/// Merge multiple combined STR files and/or individual STR files
fn merge_combined_files(combined_files: Vec<PathBuf>, individual_files: Vec<PathBuf>) {
    use std::collections::HashMap;

    // Read the first combined file into memory
    let mut combined_reader = reader(&combined_files[0].to_string_lossy()).lines();

    // Skip metadata lines and read first combined file header
    let mut combined_header = String::new();
    for line_result in &mut combined_reader {
        let line = line_result.unwrap_or_else(|e| {
            eprintln!("Error: Failed to read combined file: {}", combined_files[0].display());
            eprintln!("Reason: {}", e);
            std::process::exit(1);
        });
        if !line.starts_with('#') {
            combined_header = line;
            break;
        }
    }
    
    if combined_header.is_empty() {
        eprintln!("Error: Combined file is empty or contains only metadata: {}", combined_files[0].display());
        std::process::exit(1);
    }

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
                    .unwrap_or_else(|| {
                        eprintln!("Error: File is empty: {}", file.display());
                        std::process::exit(1);
                    })
                    .unwrap_or_else(|e| {
                        eprintln!("Error: Failed to read file: {}", file.display());
                        eprintln!("Reason: {}", e);
                        std::process::exit(1);
                    })
            })
            .collect()
    };

    let has_individual_headers = if !individual_headers.is_empty() {
        individual_headers[0].split('\t').next() == Some("chromosome")
    } else {
        // If there are no individual files, assume they match the combined file header style
        has_combined_header
    };

    // Validate header consistency (only if we have individual files to check)
    if !individual_files.is_empty() && has_combined_header != has_individual_headers {
        eprintln!("Error: Header mismatch between files.");
        eprintln!("  Combined file has {} header", if has_combined_header { "a" } else { "no" });
        eprintln!(
            "  Individual files have {} headers",
            if has_individual_headers { "a" } else { "no" }
        );
        eprintln!("\nAll files must have consistent header format.");
        std::process::exit(1);
    }

    // Output merged header
    if has_combined_header {
        // Read headers from all combined files first
        let mut all_combined_headers = vec![combined_header.clone()];
        
        for combined_file in combined_files.iter().skip(1) {
            let mut file_reader = reader(&combined_file.to_string_lossy()).lines();
            
            // Skip metadata lines and read header
            let mut header = String::new();
            for line_result in &mut file_reader {
                let line = line_result.unwrap_or_else(|e| {
                    eprintln!("Error: Failed to read combined file: {}", combined_file.display());
                    eprintln!("Reason: {}", e);
                    std::process::exit(1);
                });
                if !line.starts_with('#') {
                    header = line;
                    break;
                }
            }
            
            if header.is_empty() {
                eprintln!("Error: Combined file contains only metadata: {}", combined_file.display());
                std::process::exit(1);
            }
            
            all_combined_headers.push(header);
        }
        
        let combined_header_fields: Vec<&str> = all_combined_headers[0].split('\t').collect();

        // Validate combined file header
        if combined_header_fields.len() < 5 {
            eprintln!(
                "ERROR: Combined file header must have at least 5 columns (chr, begin, end, info, sample_H1). Got {} columns.",
                combined_header_fields.len()
            );
            std::process::exit(1);
        }

        let mut merged_header = vec![
            combined_header_fields[0], // chromosome
            combined_header_fields[1], // begin
            combined_header_fields[2], // end
            combined_header_fields[3], // info - IMPORTANT!
        ];

        // Add all sample columns from combined file (skip chr, begin, end, info)
        merged_header.extend(&combined_header_fields[4..]);

        // Add sample columns from additional combined files
        for header in all_combined_headers.iter().skip(1) {
            let header_fields: Vec<&str> = header.split('\t').collect();
            if header_fields.len() < 5 {
                eprintln!(
                    "ERROR: Combined file header must have at least 5 columns. Got {} columns",
                    header_fields.len()
                );
                std::process::exit(1);
            }
            // Add sample columns (skip chr, begin, end, info)
            merged_header.extend(&header_fields[4..]);
        }

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

    // Read first combined file data into HashMap
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

    eprintln!("Loaded {} loci from first combined file", combined_data.len());

    // Merge additional combined files (if any)
    for (file_idx, combined_file) in combined_files.iter().enumerate().skip(1) {
        eprintln!("Merging combined file {} of {}: {}", file_idx + 1, combined_files.len(), combined_file.display());
        
        let mut file_reader = reader(&combined_file.to_string_lossy()).lines();
        
        // Skip metadata lines and header
        for line_result in &mut file_reader {
            let line = line_result.unwrap_or_else(|e| {
                eprintln!("Error reading combined file: {}", e);
                std::process::exit(1);
            });
            // Skip metadata lines
            if line.starts_with('#') {
                continue;
            }
            // First non-metadata line is the header, skip it and break
            break;
        }
        
        // Read and merge data lines
        let mut merged_loci = 0;
        let mut new_loci = 0;
        
        for line_result in file_reader {
            let line = line_result.unwrap_or_else(|e| {
                eprintln!("Error reading combined file: {}", e);
                std::process::exit(1);
            });
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                eprintln!("Warning: Skipping invalid line in {}", combined_file.display());
                continue;
            }
            
            let key = format!("{}:{}−{}", fields[0], fields[1], fields[2]);
            
            if let Some(existing) = combined_data.get_mut(&key) {
                // Locus exists - merge sample columns (skip first 4 columns: chr, start, end, info)
                let new_samples = &fields[4..];
                existing.push('\t');
                existing.push_str(&new_samples.join("\t"));
                merged_loci += 1;
            } else {
                // New locus - add it
                combined_data.insert(key, line);
                new_loci += 1;
            }
        }
        
        eprintln!("  Merged {} loci, added {} new loci", merged_loci, new_loci);
    }

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
                Some(Err(e)) => {
                    panic!("Error reading file {}: {}", individual_files[file_idx].display(), e)
                }
                None => {}
            }
        }

        if all_done {
            break;
        }

        if individual_lines.len() != individual_files.len() {
            panic!(
                "Individual files have different number of data lines at line {}",
                processed_lines
            );
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
                    processed_lines,
                    chr,
                    start,
                    end,
                    file_idx + 1,
                    fields[0],
                    fields[1],
                    fields[2]
                );
            }
        }

        // Look up this locus in combined data
        let key = format!("{}:{}−{}", chr, start, end);
        let combined_line = combined_data.get(&key).unwrap_or_else(|| {
            panic!(
                "Locus {}:{}-{} found in individual files but not in combined file",
                chr, start, end
            )
        });

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
        "Completed merging {} combined file(s) with {} individual file(s)",
        combined_files.len(),
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

            // Skip metadata lines (lines starting with #) and return first non-metadata line
            loop {
                match file_reader.next() {
                    Some(Ok(line)) => {
                        if line.starts_with('#') {
                            continue; // Skip metadata line
                        } else {
                            return line; // Return first non-metadata line (header)
                        }
                    }
                    Some(Err(e)) => panic!("Error reading from {}: {}", file.display(), e),
                    None => panic!("File {} is empty or contains only metadata", file.display()),
                }
            }
        })
        .collect();

    // Validate that all files have kmer headers
    for (i, header) in headers.iter().enumerate() {
        let fields: Vec<&str> = header.split('\t').collect();
        // Accept both "kmer" (regular kmer files) and "Sample" (target kmer files)
        if fields.len() < 2 || (fields[0] != "kmer" && fields[0] != "Sample") {
            panic!(
                "Invalid kmer file header in {}: expected 'kmer' or 'Sample' as first column, got '{}'",
                calls[i].display(),
                fields.first().unwrap_or(&"<empty>")
            );
        }
    }

    headers
}

/// Output the combined header for kmer files
fn output_combined_kmer_header(headers: &[String]) {
    // Output file type metadata
    println!("# file_type=combined_kmer");
    println!("# version={}", crate::VERSION);
    println!("# command=combine");

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

    // Skip metadata lines (starting with #) and the header line in all files
    for file_reader in &mut file_readers {
        while let Some(Ok(line)) = file_reader.next() {
            if line.starts_with('#') {
                continue; // Skip metadata lines
            } else {
                // This is the header line (first non-metadata line), consume it and stop
                break;
            }
        }
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

            // Skip metadata lines (lines starting with #) and return first non-metadata line
            loop {
                match file_reader.next() {
                    Some(Ok(line)) => {
                        if line.starts_with('#') {
                            continue; // Skip metadata line
                        } else {
                            return line; // Return first non-metadata line (header or data)
                        }
                    }
                    Some(Err(e)) => panic!("Error reading from {}: {}", file.display(), e),
                    None => panic!("File {} is empty or contains only metadata", file.display()),
                }
            }
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
    if first_header_fields.len() < 6 {
        eprintln!("Error: Invalid header format in first file.");
        eprintln!(
            "Expected at least 6 columns: chromosome, begin, end, info, sample_H1, sample_H2"
        );
        eprintln!("Got {} columns: {}", first_header_fields.len(), headers[0]);
        std::process::exit(1);
    }

    // Output file type metadata
    println!("# file_type=combined_call");
    println!("# version={}", crate::VERSION);
    println!("# command=combine");

    // Start with chr, begin, end, info
    let mut combined_header = vec![
        first_header_fields[0],
        first_header_fields[1],
        first_header_fields[2],
        first_header_fields[3],
    ];

    // Add sample columns from first file (skip chr, begin, end, info)
    combined_header.extend(&first_header_fields[4..]);

    // Add sample columns from other files (skip chr, begin, end, info)
    for header in &headers[1..] {
        let fields: Vec<&str> = header.split('\t').collect();
        if fields.len() < 6 {
            eprintln!("Error: Invalid header format in one of the input files.");
            eprintln!(
                "Expected at least 6 columns: chromosome, begin, end, info, sample_H1, sample_H2"
            );
            eprintln!("Got {} columns: {}", fields.len(), header);
            std::process::exit(1);
        }
        combined_header.extend(&fields[4..]);
    }

    println!("{}", combined_header.join("\t"));
}

/// Process data lines with optimized sequential reading and parallel line processing
/// Get the maximum number of open files this process can have
fn get_max_open_files() -> usize {
    #[cfg(unix)]
    {
        use std::process::Command;

        // Try to get file descriptor limit using ulimit
        // This works on Linux, macOS, and most Unix systems
        let result = Command::new("sh")
            .arg("-c")
            .arg("ulimit -n")
            .output()
            .ok()
            .and_then(|output| String::from_utf8(output.stdout).ok())
            .and_then(|limit_str| limit_str.trim().parse::<usize>().ok());

        if let Some(limit) = result {
            if limit < 100 {
                eprintln!("  WARNING: Very low file descriptor limit detected: {}", limit);
                eprintln!(
                    "  This may cause issues with large cohorts. Consider increasing with 'ulimit -n 1024'"
                );
            }

            // Use 60% of limit to leave room for other file operations
            // Clamp to maximum 1000 (no minimum - respect system limits)
            let safe_limit = ((limit as f64 * 0.6) as usize).min(1000);
            eprintln!("  File descriptor limit: {} (using {})", limit, safe_limit);
            return safe_limit;
        }
    }

    // Fallback: conservative default for non-Unix systems or if ulimit fails
    eprintln!("  Could not detect file descriptor limit, using default: 200");
    200
}

/// Process STR data from a large number of files in batches to avoid exhausting file descriptors
/// Each batch writes to its own temp file (batch 0 has chr/start/end/info, others only genotypes).
/// Final output merges all temp files line-by-line. Memory usage is minimal.
fn process_data_batched(calls: &[PathBuf], batch_size: usize) {
    use std::fs::{File, remove_file};
    use std::io::{BufRead, BufReader, Write};

    let num_batches = calls.len().div_ceil(batch_size);
    eprintln!(
        "  Processing {} files in {} batches of up to {} files each",
        calls.len(),
        num_batches,
        batch_size
    );

    let mut temp_files: Vec<PathBuf> = Vec::new();

    // Process each batch and write to separate temp file
    for (batch_idx, file_chunk) in calls.chunks(batch_size).enumerate() {
        eprintln!(
            "  Processing batch {}/{} ({} files)...",
            batch_idx + 1,
            num_batches,
            file_chunk.len()
        );

        let temp_file_path = std::env::temp_dir().join(format!(
            "inquistr_combine_{}_{}.tmp",
            std::process::id(),
            batch_idx
        ));
        temp_files.push(temp_file_path.clone());

        let mut temp_file = File::create(&temp_file_path).expect("Failed to create temporary file");

        let mut file_readers: Vec<_> = file_chunk
            .iter()
            .map(|file| reader(&file.clone().into_os_string().into_string().unwrap()).lines())
            .collect();

        // Skip metadata and header
        for file_reader in &mut file_readers {
            while let Some(Ok(line)) = file_reader.next() {
                if line.starts_with('#') {
                    continue;
                } else {
                    break;
                }
            }
        }

        let mut line_count = 0;
        let chunk_size = 10000;

        loop {
            let mut data_lines: Vec<String> = Vec::with_capacity(file_chunk.len());
            let mut all_done = true;

            for file_reader in file_readers.iter_mut() {
                match file_reader.next() {
                    Some(Ok(line)) => {
                        data_lines.push(line);
                        all_done = false;
                    }
                    None => {}
                    Some(Err(e)) => panic!("Error reading file: {}", e),
                }
            }

            if all_done {
                break;
            }

            if batch_idx == 0 {
                // First batch: write chr/start/end/info + genotypes
                let combined_line = combine_data_lines(&data_lines, line_count);
                writeln!(temp_file, "{}", combined_line).expect("Failed to write to temp file");
            } else {
                // Subsequent batches: only write genotypes (no chr/start/end/info)
                for (idx, line) in data_lines.iter().enumerate() {
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() >= 2 {
                        if idx > 0 {
                            write!(temp_file, "\t").expect("Failed to write tab");
                        }
                        write!(
                            temp_file,
                            "{}\t{}",
                            fields[fields.len() - 2],
                            fields[fields.len() - 1]
                        )
                        .expect("Failed to write genotypes");
                    }
                }
                writeln!(temp_file).expect("Failed to write newline");
            }

            line_count += 1;
            if line_count % chunk_size == 0 {
                eprintln!("    Processed {} variants...", line_count);
            }
        }

        eprintln!(
            "  Batch {} complete: {} variants written to {}",
            batch_idx + 1,
            line_count,
            temp_file_path.display()
        );

        drop(file_readers);
        drop(temp_file);
    }

    // Merge all temp files line-by-line and output
    eprintln!("  Merging {} temporary files...", temp_files.len());

    let mut temp_readers: Vec<_> = temp_files
        .iter()
        .map(|path| BufReader::new(File::open(path).expect("Failed to open temp file")).lines())
        .collect();

    let mut line_count = 0;
    loop {
        let mut merged_line = String::new();
        let mut all_done = true;

        for (idx, temp_reader) in temp_readers.iter_mut().enumerate() {
            match temp_reader.next() {
                Some(Ok(line)) => {
                    if idx == 0 {
                        merged_line = line;
                    } else {
                        merged_line.push('\t');
                        merged_line.push_str(&line);
                    }
                    all_done = false;
                }
                None => {}
                Some(Err(e)) => panic!("Error reading temp file {}: {}", idx, e),
            }
        }

        if all_done {
            break;
        }

        println!("{}", merged_line);
        line_count += 1;
    }

    // Clean up temp files
    for temp_file in temp_files {
        remove_file(&temp_file).ok();
    }

    eprintln!("  Done! Wrote {} variants", line_count);
}

fn process_data_parallel(calls: &[PathBuf]) {
    // Dynamically determine safe number of open files
    let max_open_files = get_max_open_files();

    if calls.len() > max_open_files {
        eprintln!(
            "  Large cohort ({} files), processing in batches of {} to avoid file descriptor exhaustion",
            calls.len(),
            max_open_files
        );
        process_data_batched(calls, max_open_files);
        return;
    }

    // Create file readers for small cohorts (original approach)
    let mut file_readers: Vec<_> = calls
        .iter()
        .map(|file| reader(&file.clone().into_os_string().into_string().unwrap()).lines())
        .collect();

    // Skip metadata lines (starting with #) and the header line in all files
    // All files must have headers (validated earlier), so we skip until we've consumed the header
    for file_reader in &mut file_readers {
        while let Some(Ok(line)) = file_reader.next() {
            if line.starts_with('#') {
                continue; // Skip metadata lines
            } else {
                // This is the header line (first non-metadata line), consume it and stop
                break;
            }
        }
    }

    let mut line_count = 0;
    let chunk_size = 100000; // Report progress every 10k lines

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
            eprintln!("  Processed {} variants...", line_count);
        }
    }

    eprintln!("  Completed processing {} variants", line_count);
}

/// Combine data lines from all files with variant consistency checking
fn combine_data_lines(data_lines: &[String], line_number: usize) -> String {
    if data_lines.is_empty() {
        panic!("No data lines to combine at line {}", line_number);
    }

    // Parse first line to get coordinates
    let first_line_fields: Vec<&str> = data_lines[0].split('\t').collect();
    if first_line_fields.len() < 6 {
        panic!(
            "Invalid data line format at line {} (expected at least 6 columns): {}",
            line_number, data_lines[0]
        );
    }

    let (chr, start, end) = (first_line_fields[0], first_line_fields[1], first_line_fields[2]);
    let expected_field_count = first_line_fields.len();

    // Validate that all files have the same variant coordinates and field count
    for (file_idx, line) in data_lines.iter().enumerate().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();

        // Validate field count matches
        if fields.len() != expected_field_count {
            panic!(
                "Field count mismatch at line {}: file 0 has {} fields, file {} has {} fields",
                line_number,
                expected_field_count,
                file_idx + 1,
                fields.len()
            );
        }

        // Validate coordinates match
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

    // Add coordinates and info from first file
    combined_fields.push(chr);
    combined_fields.push(start);
    combined_fields.push(end);
    combined_fields.push(first_line_fields[3]); // info column

    // Add all sample data from first file (skip chr, start, end, info)
    combined_fields.extend(&first_line_fields[4..]);

    // Add sample data from other files (skip chr, start, end, info)
    for line in &data_lines[1..] {
        let fields: Vec<&str> = line.split('\t').collect();
        combined_fields.extend(&fields[4..]);
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
