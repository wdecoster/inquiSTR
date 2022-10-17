use clap::ValueEnum;
use std::io::BufRead;
use std::path::PathBuf;

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
    let samples_of_interest = parse_phenotypes(metadata, condition);

    let file = crate::utils::reader(&combined.into_os_string().into_string().unwrap());
    let mut lines = file.lines();
    let line = lines.next().unwrap().unwrap();
    let samples: Vec<&str> = line.split('\t').skip(3).collect();
    let sample_indices_of_interest: Vec<usize> = vec![1];
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

fn parse_phenotypes(metadata: PathBuf, condition: String) -> Vec<String> {
    unimplemented!()
}
