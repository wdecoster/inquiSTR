use crate::locus_search::{
    LocusSearchConfig, OverlapStrategy, extract_clean_sample_names, select_overlapping_locus,
};
use kuva::plot::legend::{LegendEntry, LegendPosition, LegendShape};
use kuva::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;

pub fn plot(
    combined: PathBuf,
    sample_metadata: Option<PathBuf>,
    condition: Option<String>,
    region: String,
    min: Option<f64>,
    max: Option<f64>,
    output: String,
) {
    if let (Some(min), Some(max)) = (min, max)
        && min > max
    {
        eprintln!("ERROR: --min must be less than or equal to --max.");
        std::process::exit(1);
    }

    if !combined.exists() {
        eprintln!("ERROR: Combined file does not exist: {}", combined.display());
        std::process::exit(1);
    }
    if let Some(sample_metadata) = sample_metadata.as_ref()
        && !sample_metadata.exists()
    {
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

    // Match any locus overlapping the requested region.
    let combined_for_header = combined.clone();
    let config = LocusSearchConfig {
        combined_file: combined,
        target_region: region.clone(),
        overlap_strategy: OverlapStrategy::Overlap,
    };

    let locus_match = select_overlapping_locus(config);

    let mut lengths_for_plot: HashMap<String, Vec<f64>> = HashMap::new();

    if let (Some(sample_metadata), Some(condition)) = (sample_metadata, condition) {
        // Read header to map sample columns to metadata groups.
        let file = crate::utils::reader(&combined_for_header.to_string_lossy());
        let mut lines = std::io::BufRead::lines(file);
        let header_line =
            crate::utils::skip_metadata_lines(&mut lines, &combined_for_header.to_string_lossy());
        let samples: Vec<String> = extract_clean_sample_names(&header_line);

        let samples_of_interest =
            crate::sample_info::parse_phenotypes(&sample_metadata, &condition)
                .expect("Problem parsing sample metadata file");
        let mut samples_map = HashMap::with_capacity(samples_of_interest.len());
        for s in samples_of_interest {
            samples_map.insert(s.identifier, s.group);
        }

        for (sample, length) in samples.iter().zip(locus_match.values) {
            if min.is_some_and(|v| length < v) || max.is_some_and(|v| length > v) {
                continue;
            }
            if let Some(group) = samples_map.get(sample) {
                lengths_for_plot
                    .entry(group.clone())
                    .or_default()
                    .push(length);
            }
        }
    } else {
        let filtered: Vec<f64> = locus_match
            .values
            .into_iter()
            .filter(|length| !min.is_some_and(|v| *length < v) && !max.is_some_and(|v| *length > v))
            .collect();
        lengths_for_plot.insert(String::from("all"), filtered);
    }

    plot_hist(lengths_for_plot, output);
}

fn plot_hist(lengths_map: HashMap<String, Vec<f64>>, output: String) {
    if lengths_map.is_empty() || lengths_map.values().all(|lengths| lengths.is_empty()) {
        eprintln!("ERROR: No values available to plot for the requested inputs.");
        std::process::exit(1);
    }

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
    let nice_range = compute_nice_range(range.0, range.1);

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
    let mut legend_entries: Vec<LegendEntry> = Vec::new();
    for (i, (group, lengths)) in lengths_map.into_iter().enumerate() {
        let color = colors[i % colors.len()];
        let legend_color = color.get(..7).unwrap_or(color).to_string();
        let hist = Histogram::new()
            .with_data(lengths)
            .with_bins(100)
            .with_range(nice_range)
            .with_color(color)
            .with_legend(&group);
        legend_entries.push(LegendEntry {
            label: group.clone(),
            color: legend_color,
            shape: LegendShape::Rect,
            dasharray: None,
        });
        plots.push(hist.into());
    }

    let layout = Layout::auto_from_plots(&plots)
        .with_legend_position(LegendPosition::InsideTopRight)
        .with_legend_title("Groups")
        .with_legend_entries(legend_entries)
        .with_interactive();

    let svg = render_to_svg(plots, layout);
    std::fs::write(&output, svg).expect("Failed to write SVG plot");
}

fn compute_nice_range(min: f64, max: f64) -> (f64, f64) {
    if !min.is_finite() || !max.is_finite() {
        return (0.0, 1.0);
    }

    if (max - min).abs() < f64::EPSILON {
        let pad = if min.abs() < 1.0 {
            1.0
        } else {
            10f64.powf(min.abs().log10().floor()) * 0.1
        };
        return (min - pad, max + pad);
    }

    let span = max - min;
    let target_step = span / 5.0;
    let magnitude = 10f64.powf(target_step.log10().floor());
    let normalized = target_step / magnitude;

    let nice_multiplier = if normalized <= 1.0 {
        1.0
    } else if normalized <= 2.0 {
        2.0
    } else if normalized <= 5.0 {
        5.0
    } else {
        10.0
    };

    let step = nice_multiplier * magnitude;
    let nice_min = (min / step).floor() * step;
    let nice_max = (max / step).ceil() * step;
    (nice_min, nice_max)
}
