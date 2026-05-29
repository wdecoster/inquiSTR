use crate::locus_search::{LocusSearchConfig, OverlapStrategy, find_locus};
use histo_fp::Histogram;
use std::path::PathBuf;

pub fn histogram(combined: PathBuf, region: String) {
    if !combined.exists() {
        eprintln!("ERROR: Combined file does not exist: {}", combined.display());
        std::process::exit(1);
    }

    // Validate that input is a combined file, not individual; capture type for axis label
    let is_kmer = match crate::filetype::read_file_type_metadata(&combined) {
        Some(file_type) => {
            if !matches!(
                file_type,
                crate::filetype::FileType::CombinedCall | crate::filetype::FileType::CombinedKmer
            ) {
                eprintln!(
                    "ERROR: Histogram requires a combined file (combined_call or combined_kmer)."
                );
                eprintln!("The provided file appears to be: {:?}", file_type);
                eprintln!(
                    "\nPlease use 'inquiSTR combine' to merge individual sample files first."
                );
                std::process::exit(1);
            }
            matches!(file_type, crate::filetype::FileType::CombinedKmer)
        }
        None => false,
    };

    let x_label = if is_kmer {
        "kmer frequency"
    } else {
        "STR length (bp)"
    };

    // Match any locus overlapping the requested region (the first one found).
    let combined_path = combined.clone();
    let config = LocusSearchConfig {
        combined_file: combined,
        target_region: region.clone(),
        overlap_strategy: OverlapStrategy::Overlap,
    };

    if let Some(locus_match) = find_locus(config) {
        let mut histogram = Histogram::with_buckets(100, Some(2));

        for value in locus_match.values {
            if !value.is_nan() {
                histogram.add(value);
            }
        }

        println!("# Locus: {}:{}-{}", locus_match.chromosome, locus_match.start, locus_match.end);
        println!("# x-axis: {x_label}");
        println!("# y-axis: count");
        println!("{histogram}");
    } else {
        eprintln!(
            "ERROR: No locus overlapping region '{}' was found in {}.\n\
             Check the coordinates (expected chrom:begin-end) and that a locus in this interval \
             is present in the combined file.",
            region,
            combined_path.display()
        );
        std::process::exit(1);
    }
}
