use crate::locus_search::{
    LocusSearchConfig, OverlapStrategy, extract_clean_sample_names, find_locus,
};
use kuva::prelude::*;
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
    if let Some(file_type) = crate::filetype::read_file_type_metadata(&combined)
        && !matches!(
            file_type,
            crate::filetype::FileType::CombinedCall | crate::filetype::FileType::CombinedKmer
        )
    {
        eprintln!("ERROR: Plot requires a combined file (combined_call or combined_kmer).");
        eprintln!("The provided file appears to be: {:?}", file_type);
        eprintln!("\nPlease use 'inquiSTR combine' to merge individual sample files first.");
        std::process::exit(1);
    }

    // Read header to get sample names
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = std::io::BufRead::lines(file);
    // Skip metadata lines if present
    let header_line = crate::utils::skip_metadata_lines(&mut lines, &combined.to_string_lossy());
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

    for (sample, length) in samples.iter().zip(locus_match.values) {
        if samples_map.contains_key(sample) {
            lengths_for_plot
                .entry(samples_map[sample].clone())
                .or_default()
                .push(length);
        }
    }
    plot_hist(lengths_for_plot, output);
}

fn plot_hist(lengths_map: HashMap<String, Vec<f64>>, output: String) {
    // Compute shared range across all groups
    let mut global_min = f64::INFINITY;
    let mut global_max = f64::NEG_INFINITY;
    for lengths in lengths_map.values() {
        for &v in lengths {
            if v < global_min {
                global_min = v;
            }
            if v > global_max {
                global_max = v;
            }
        }
    }
    let range = (global_min, global_max);

    // Semi-transparent colors for overlapping histograms (~50% alpha)
    let colors = [
        "#4682b480",
        "#dc143c80",
        "#2e8b5780",
        "#ff8c0080",
        "#8a2be280",
        "#00ced180",
    ];

    let mut plots: Vec<Plot> = Vec::new();
    for (i, (group, lengths)) in lengths_map.into_iter().enumerate() {
        let color = colors[i % colors.len()];
        let hist = Histogram::new()
            .with_data(lengths)
            .with_bins(100)
            .with_range(range)
            .with_color(color)
            .with_legend(&group);
        plots.push(hist.into());
    }

    let layout = Layout::auto_from_plots(&plots).with_interactive();

    let svg = render_to_svg(plots, layout);
    std::fs::write(&output, svg).expect("Failed to write SVG plot");
}
