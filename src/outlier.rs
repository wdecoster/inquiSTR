use clap::ValueEnum;

use dbscan::Classification::*;
use dbscan::Model;
use log::debug;

use std::cmp::max;
use std::cmp::Ordering;
use std::io::BufRead;
use std::path::PathBuf;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Method {
    Zscore,
    Dbscan,
}

fn std_deviation_and_mean(data: &Vec<f32>) -> (f32, f32) {
    let sum = data.iter().sum::<f32>();
    let count = data.len() as f32;
    let data_mean = sum / count;
    let variance = data
        .iter()
        .map(|value| {
            let diff = data_mean - *value;
            diff * diff
        })
        .sum::<f32>()
        / count;
    (data_mean, variance.sqrt())
}

pub fn outlier(combined: PathBuf, minsize: u32, zscore_cutoff: f32, method: Method) {
    if !combined.exists() {
        panic!("Combined file does not exist!");
    }
    let file = crate::utils::reader(&combined.into_os_string().into_string().unwrap());
    let mut lines = file.lines();
    let line = lines.next().unwrap().unwrap();
    println!("chrom\tbegin\tend\toutliers");
    let samples: Vec<&str> = line.split('\t').skip(3).collect();
    let mincluster = samples.len().ilog2() as usize;
    for line in lines {
        let line = line.unwrap();
        let splitline = line.split('\t').collect::<Vec<&str>>();
        let (chrom, begin, end) = (&splitline[0], &splitline[1], &splitline[2]);
        if let Some(values) = get_repeat_lengths(&splitline, minsize) {
            let expanded = match method {
                Method::Zscore => z_score_outliers(values, &samples, zscore_cutoff),
                Method::Dbscan => dbscan_outliers(values, &samples, mincluster),
            };
            if !expanded.is_empty() {
                debug!(
                    "chrom: {}, begin: {}, end: {}, N_expanded: {}, expanded: {:?}",
                    chrom,
                    begin,
                    end,
                    expanded.len(),
                    expanded
                );
                let expanded = expanded.join(",");
                println!("{chrom}\t{begin}\t{end}\t{expanded}")
            }
        }
    }
}

fn get_repeat_lengths(line: &[&str], minsize: u32) -> Option<Vec<f32>> {
    let values: Vec<f32> = line
        .iter()
        .skip(3)
        .map(|number| number.parse().expect("Failed to parse number"))
        .collect();
    let values = values
        .iter()
        .map(|&value| if value.is_nan() { 0.0 } else { value })
        .collect::<Vec<f32>>();
    // Check if the maximum value is larger than the minimum size
    // If all values are NaN then the vector will contain only zeroes
    if values
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Less))
        .unwrap()
        < &(minsize as f32)
    {
        None
    } else {
        Some(values)
    }
}

fn z_score_outliers<'a>(values: Vec<f32>, samples: &[&'a str], zscore_cutoff: f32) -> Vec<&'a str> {
    // calculate mean and std deviation of the STR lengths
    let (values_mean, values_std_dev) = std_deviation_and_mean(&values);
    debug!("mean: {}, std_dev: {}", values_mean, values_std_dev);
    // calculate the zscore for each haplotype and get the haplotype identifier based on the index if larger zscore > cutoff
    // intentionally this only selects for values that are larger than the mean
    // and therefore only for expansions, not contractions
    values
        .iter()
        .enumerate()
        .filter(|(_, &value)| ((value - values_mean) / values_std_dev) >= zscore_cutoff)
        .map(|(index, _)| samples[index])
        .collect::<Vec<&str>>()
}

fn dbscan_outliers<'a>(values: Vec<f32>, samples: &[&'a str], mincluster: usize) -> Vec<&'a str> {
    // the parameters for the dbscan model are as used by the schizophrenia STR outlier paper (https://doi.org/10.1038/s41380-022-01857-4)
    // however, the eps parameter is set as minimally 10
    let eps = max(2 * mode(&values), 10) as f64;
    let values = values
        .iter()
        .map(|&value| vec![value])
        .collect::<Vec<Vec<f32>>>();
    let model = Model::new(eps, mincluster);
    let output = model.run(&values);
    debug!("eps: {}, mincluster: {}", eps, mincluster);
    debug!("output: {:?}", output);
    output
        .iter()
        .enumerate()
        .filter(|(_, &classification)| matches!(classification, Noise))
        .map(|(index, _)| samples[index])
        .collect::<Vec<&str>>()
}

fn mode(values: &[f32]) -> usize {
    // calculate the mode of the STR lengths
    // as NaN values are replaced by 0.0, we need to filter out the 0.0 values
    // if not eps will often be 0
    let mut counts = std::collections::HashMap::new();
    for &value in values.iter().filter(|&&value| value > 0.0) {
        *counts.entry(value as usize).or_insert(0) += 1;
    }
    counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(value, _)| value)
        .expect("No mode found for repeat")
}

#[cfg(test)]
#[test]
fn test_dbscan_outliers() {
    let values = vec![1.0, 2.0, 2.0, 3.0, 1.0, 5.0, 3.0, 2.0, 2.0, 1.0, 120.0];
    let samples = vec![
        "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11",
    ];
    let expected = vec!["s11"];
    let mincluster = values.len().ilog2() as usize;
    assert_eq!(dbscan_outliers(values, &samples, mincluster), expected);
}

#[test]
fn test_z_score_outliers() {
    let values = vec![1.0, 2.0, 2.0, 3.0, 1.0, 5.0, 3.0, 2.0, 2.0, 1.0, 120.0];
    let samples = vec![
        "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11",
    ];
    let expected = vec!["s11"];
    let zscore_cutoff = 2.0;
    assert_eq!(z_score_outliers(values, &samples, zscore_cutoff), expected);
}
