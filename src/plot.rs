use plotly::{Histogram, Plot};
use std::collections::HashMap;
use std::path::PathBuf;
use crate::locus_search::{LocusSearchConfig, find_locus, extract_clean_sample_names, OverlapStrategy};

pub fn plot(
    combined: PathBuf,
    metadata: PathBuf,
    condition: String,
    region: String,
    output: String,
) {
    if !combined.exists() {
        panic!("Combined file does not exist!");
    }
    if !metadata.exists() {
        panic!("Metadata file does not exist!");
    }

    // Read header to get sample names
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = std::io::BufRead::lines(file);
    let header_line = lines.next().unwrap().unwrap();
    let samples: Vec<String> = extract_clean_sample_names(&header_line);

    let samples_of_interest = crate::metadata::parse_phenotypes(&metadata, &condition)
        .expect("Problem parsing metadata file");
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
