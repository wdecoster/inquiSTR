use crate::locus_search::{find_locus, LocusSearchConfig, OverlapStrategy};
use histo_fp::Histogram;
use std::path::PathBuf;

pub fn histogram(combined: PathBuf, region: String) {
    if !combined.exists() {
        eprintln!("ERROR: Combined file does not exist: {}", combined.display());
        std::process::exit(1);
    }

    // Validate that input is a combined file, not individual
    if let Some(file_type) = crate::filetype::read_file_type_metadata(&combined) {
        if !matches!(
            file_type,
            crate::filetype::FileType::CombinedCall | crate::filetype::FileType::CombinedKmer
        ) {
            eprintln!(
                "ERROR: Histogram requires a combined file (combined_call or combined_kmer)."
            );
            eprintln!("The provided file appears to be: {:?}", file_type);
            eprintln!("\nPlease use 'inquiSTR combine' to merge individual sample files first.");
            std::process::exit(1);
        }
    }

    // Use the new locus search utility with containment strategy (original behavior)
    let config = LocusSearchConfig {
        combined_file: combined,
        target_region: region,
        overlap_strategy: OverlapStrategy::Containment,
    };

    if let Some(locus_match) = find_locus(config) {
        let mut histogram = Histogram::with_buckets(100, Some(2));

        for value in locus_match.values {
            if !value.is_nan() {
                histogram.add(value);
            }
        }

        println!("{histogram}");
    } else {
        eprintln!("No matching interval found");
    }
}
