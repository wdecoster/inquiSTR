use clap::ValueEnum;
use itertools::Itertools;
use log::error;
use std::cmp::Ordering;
use std::io::BufRead;
use std::path::{Path, PathBuf};

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Mode {
    /// Use the sum of H1 and H2
    Sum,
    /// Use the max of H1 and H2
    Max,
    /// Use the min of H1 and H2
    Min,
}

pub fn assocation(
    combined: PathBuf,
    metadata: PathBuf,
    missing_cutoff: f32,
    mode: Mode,
    condition: String,
    covariates: Option<String>,
) {
    // Taking out the samples
    let samples_of_interest = crate::metadata::parse_phenotypes(&metadata, &condition)
        .expect("Problem parsing metadata file");
    if samples_of_interest.len() < 2 {
        error!(
            "{}",
            format!(
                "Not enough samples in {} matching {}",
                metadata.display(),
                condition
            )
        );
        panic!()
    }

    let file = crate::utils::reader(&combined.into_os_string().into_string().unwrap());
    let mut lines = file.lines();
    let line = lines.next().unwrap().unwrap();
    // TODO: may need to be changed if the sample name is used differently by call and combine
    // The vector below should contain each sample twice, for each haplotype...
    let samples: Vec<String> = line
        .split('\t')
        .skip(3)
        .map(|s| s.replace(".inq_H1", "").replace(".inq_H2", ""))
        .collect();
    // iterate over the Individual instances to get a vector of samples of interest
    let sample_names_of_interest = samples_of_interest
        .iter()
        .map(|s| &s.identifier)
        .collect::<Vec<&String>>();
    // Using the vector of samples and the vector of samples of interest, prepare a vector of indices of interest
    let sample_indices_of_interest: Vec<usize> = samples
        .iter()
        .enumerate()
        .filter(|(_, s)| sample_names_of_interest.contains(s))
        .map(|(index, _)| index)
        .collect();
    assert!(sample_indices_of_interest.len() == 2 * sample_names_of_interest.len());
    // For each variant, take out the values that correspond with the sample_indices_of_interest
    for line in lines {
        let line = line.unwrap();
        let splitline = line.split('\t').collect::<Vec<&str>>();
        let chrom = &splitline[0];
        let begin = &splitline[1];
        let end = &splitline[2];
        let values: Vec<f32> = splitline
            .iter()
            .skip(3)
            .enumerate()
            .filter(|&(i, _)| sample_indices_of_interest.contains(&i))
            .map(|(_, number)| number.parse().unwrap())
            .collect();
        let summarized_values = summarize_values(values, mode);
    }
}

fn summarize_values(values: Vec<f32>, mode: Mode) -> Vec<f32> {
    match mode {
        Mode::Max => values
            .iter()
            .chunks(2)
            .into_iter()
            .map(|chunk| {
                chunk
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Less))
                    .unwrap()
                    .to_owned()
            })
            .collect::<Vec<f32>>(),
        Mode::Min => values
            .iter()
            .chunks(2)
            .into_iter()
            .map(|chunk| {
                chunk
                    .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Less))
                    .unwrap()
                    .to_owned()
            })
            .collect::<Vec<f32>>(),
        Mode::Sum => values
            .iter()
            .chunks(2)
            .into_iter()
            .map(|chunk| chunk.sum())
            .collect::<Vec<f32>>(),
    }
}
