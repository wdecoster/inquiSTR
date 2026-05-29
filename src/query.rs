use crate::locus_search::{OverlapStrategy, extract_sample_names, find_multiple_loci};
use log::debug;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::PathBuf;

pub fn query(combined: PathBuf, region: String) {
    if !combined.exists() {
        eprintln!("ERROR: Combined file does not exist: {}", combined.display());
        std::process::exit(1);
    }

    // Validate that input is a combined file, not individual
    if let Some(file_type) = crate::filetype::read_file_type_metadata(&combined)
        && !matches!(
            file_type,
            crate::filetype::FileType::CombinedCall | crate::filetype::FileType::CombinedKmer
        )
    {
        eprintln!("ERROR: Query requires a combined file (combined_call or combined_kmer).");
        eprintln!("The provided file appears to be: {:?}", file_type);
        eprintln!("\nPlease use 'inquiSTR combine' to merge individual sample files first.");
        std::process::exit(1);
    }

    // Read header to get sample names
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();
    // Skip metadata lines if present
    let header_line = crate::utils::skip_metadata_lines(&mut lines, &combined.to_string_lossy());
    let samples = extract_sample_names(&header_line);
    debug!("Samples: {:?}", samples);

    // Use the new locus search utility with overlap strategy (original behavior)
    let matches = find_multiple_loci(combined, region, OverlapStrategy::Overlap);

    if matches.is_empty() {
        eprintln!("No matching intervals found in file");
        return;
    }

    // Organize data by sample for display
    let mut lengths: HashMap<String, Vec<f64>> = HashMap::new();
    let mut matching_intervals = Vec::new();

    for locus_match in matches {
        // Create interval name
        let interval_name =
            format!("{}:{}-{}", locus_match.chromosome, locus_match.start, locus_match.end);
        matching_intervals.push(interval_name);

        // Collect values by sample
        for (sample_name, &value) in samples.iter().zip(locus_match.values.iter()) {
            lengths.entry(sample_name.clone()).or_default().push(value);
        }
    }

    // Display results based on number of matching intervals
    match matching_intervals.len() {
        0 => eprintln!("No matching intervals found in file"), // Already handled above
        1 => {
            // Single interval: print sorted by value
            println!("name\t{}", matching_intervals[0]);

            // Sort samples by their STR length (descending, NaN values last)
            let mut zipped: Vec<(String, Vec<f64>)> = lengths.into_iter().collect();
            zipped.sort_by(|(_, a), (_, b)| {
                let (a, b) = (a[0], b[0]);
                match (a.is_nan(), b.is_nan()) {
                    (true, true) => std::cmp::Ordering::Equal,
                    (true, false) => std::cmp::Ordering::Greater, // NaN sorts last
                    (false, true) => std::cmp::Ordering::Less,
                    (false, false) => b.total_cmp(&a), // descending by value
                }
            });

            for (name, val) in zipped {
                println!("{name}\t{length}", length = val[0]);
            }
        }
        _ => {
            // Multiple intervals: print as table
            println!("name\t{}", matching_intervals.join("\t"));

            for (name, val) in lengths {
                let length_strings: Vec<String> = val.iter().map(|x| x.to_string()).collect();
                println!("{name}\t{length}", length = length_strings.join("\t"));
            }
        }
    }
}
