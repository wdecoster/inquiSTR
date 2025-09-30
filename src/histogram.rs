use crate::locus_search::{find_locus, LocusSearchConfig, OverlapStrategy};
use histo_fp::Histogram;
use std::path::PathBuf;

pub fn histogram(combined: PathBuf, region: String) {
    if !combined.exists() {
        panic!("Combined file does not exist!");
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
