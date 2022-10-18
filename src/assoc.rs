use clap::ValueEnum;
use log::error;
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

struct Individual {
    identifier: String,
    group: String,
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
    let samples_of_interest =
        parse_phenotypes(&metadata, &condition).expect("Problem parsing metadata file");
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
    unimplemented!();
    match mode {
        Mode::Max => vec![1.0],
        Mode::Min => vec![1.0],
        Mode::Sum => vec![1.0],
    }
}

fn parse_phenotypes(
    metadata: &Path,
    condition: &str,
) -> Result<Vec<Individual>, Box<dyn std::error::Error>> {
    let pheno_column = condition
        .split(':')
        .next()
        .expect("Issue parsing condition string");
    let pheno_values = condition
        .split(':')
        .nth(1)
        .unwrap()
        .split(',')
        .collect::<Vec<&str>>();
    let meta_file = crate::utils::reader(metadata.to_str().unwrap());
    let mut lines = meta_file.lines();
    let header = lines.next().unwrap().unwrap();
    let pheno_column_index = header
        .split('\t')
        .enumerate()
        .filter(|(_, col)| col == &pheno_column)
        .map(|(index, _)| index)
        .next()
        .unwrap_or_else(|| {
            panic!(
                "Could not find column {} in {}",
                pheno_column,
                metadata.display()
            )
        });
    let mut samples_of_interest: Vec<Individual> = vec![];
    for line in lines {
        let line = line.unwrap();
        let splitline = line.split('\t').collect::<Vec<&str>>();
        let pheno_value = splitline.get(pheno_column_index).unwrap();
        if pheno_values.contains(pheno_value) {
            samples_of_interest.push(Individual {
                identifier: splitline.first().unwrap().to_string(),
                group: pheno_value.to_string(),
            })
        }
    }
    Ok(samples_of_interest)
}
