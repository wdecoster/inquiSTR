use plotly::{Histogram, Plot};
use std::collections::HashMap;
use std::io::BufRead;
use std::io::{BufReader, Lines, Read};
use std::path::PathBuf;

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
    let file = crate::utils::reader(&combined.into_os_string().into_string().unwrap());
    let mut lines = file.lines();
    let header_line = lines.next().unwrap().unwrap();
    let samples: Vec<String> = header_line
        .split('\t')
        .skip(3)
        .map(|s| s.replace("_H1", "").replace("_H2", ""))
        .collect();

    let samples_of_interest = crate::metadata::parse_phenotypes(&metadata, &condition)
        .expect("Problem parsing metadata file");
    let mut samples_map = HashMap::with_capacity(samples_of_interest.len());
    for s in samples_of_interest {
        samples_map.insert(s.identifier, s.group);
    }
    let lengths = get_str_lengths(region, lines).expect("Specified interval not found!");
    let mut lengths_for_plot: HashMap<String, Vec<f64>> = HashMap::new();
    let mut ids_for_plot: HashMap<String, Vec<&String>> = HashMap::new();
    for (sample, length) in samples.iter().zip(lengths) {
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

fn get_str_lengths(region: String, lines: Lines<BufReader<Box<dyn Read>>>) -> Option<Vec<f64>> {
    let (chrom, reg_start, reg_end) = crate::utils::process_region(region).unwrap();
    // Add a tab character to the chromosome so we can search for this with starts_with below (to make sure chr1 does not match chr15)
    let reg_chrom = format!("{chrom}\t");
    for line in lines {
        let line = line.unwrap();
        if line.starts_with(&reg_chrom) {
            let splitline = line.split('\t').collect::<Vec<&str>>();
            let begin: u32 = splitline[1].parse().expect("Failed parsing interval");
            let end: u32 = splitline[2].parse().expect("Failed parsing interval");
            if reg_start <= begin && end <= reg_end {
                let values = splitline
                    .iter()
                    .skip(3)
                    .map(|number| number.parse::<f64>().expect("Failed parsing lengths"))
                    .collect::<Vec<f64>>();
                return Some(values);
            }
        }
    }
    None
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
