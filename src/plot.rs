use crate::locus_search::{
    extract_clean_sample_names, find_locus, LocusSearchConfig, OverlapStrategy,
};
use plotly::{Histogram, Plot};
use std::collections::HashMap;
use std::path::PathBuf;

pub fn plot(
    combined: PathBuf,
    sample_metadata: PathBuf,
    condition: String,
    region: String,
    output: String,
) {
    if !combined.exists() {
        eprintln!("ERROR: Combined file does not exist: {}", combined.display());
        std::process::exit(1);
    }
    if !sample_metadata.exists() {
        eprintln!("ERROR: Sample metadata file does not exist: {}", sample_metadata.display());
        std::process::exit(1);
    }

    // Validate that input is a combined file, not individual
    if let Some(file_type) = crate::filetype::read_file_type_metadata(&combined) {
        if !matches!(
            file_type,
            crate::filetype::FileType::CombinedCall | crate::filetype::FileType::CombinedKmer
        ) {
            eprintln!("ERROR: Plot requires a combined file (combined_call or combined_kmer).");
            eprintln!("The provided file appears to be: {:?}", file_type);
            eprintln!("\nPlease use 'inquiSTR combine' to merge individual sample files first.");
            std::process::exit(1);
        }
    }

    // Read header to get sample names
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = std::io::BufRead::lines(file);
    // Skip metadata lines if present
    let header_line = crate::utils::skip_metadata_lines(&mut lines);
    let samples: Vec<String> = extract_clean_sample_names(&header_line);

    let samples_of_interest = crate::sample_info::parse_phenotypes(&sample_metadata, &condition)
        .expect("Problem parsing sample metadata file");
    let mut samples_map = HashMap::with_capacity(samples_of_interest.len());
    for s in samples_of_interest {
        samples_map.insert(s.identifier, s.group);
    }

    // Use the new locus search utility with containment strategy (original behavior)
    let config = LocusSearchConfig {
        combined_file: combined,
        target_region: region,
        overlap_strategy: OverlapStrategy::Containment,
    };

    let locus_match = find_locus(config).expect("Specified interval not found!");

    let mut lengths_for_plot: HashMap<String, Vec<f64>> = HashMap::new();
    let mut ids_for_plot: HashMap<String, Vec<&String>> = HashMap::new();

    for (sample, length) in samples.iter().zip(locus_match.values) {
        if samples_map.contains_key(sample) {
            lengths_for_plot
                .entry(samples_map[sample].clone())
                .or_default()
                .push(length);
            ids_for_plot
                .entry(samples_map[sample].clone())
                .or_default()
                .push(sample);
        }
    }
    plot_hist(lengths_for_plot, ids_for_plot, output);
}

fn plot_hist(
    lengths_map: HashMap<String, Vec<f64>>,
    ids_map: HashMap<String, Vec<&String>>,
    output: String,
) {
    let mut plot = Plot::new();
    for (group, lengths) in lengths_map {
        let trace = Histogram::new(lengths)
            .name(&group)
            .opacity(0.5)
            .text_array(ids_map[&group].clone());
        plot.add_trace(trace);
    }

    // plot.set_layout(Layout::new().bar_mode(BarMode::Overlay));
    plot.write_html(output);
}
