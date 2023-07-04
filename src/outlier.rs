use clap::ValueEnum;
use linfa::traits::Transformer;
use linfa_clustering::Dbscan;
use ndarray::prelude::*;

use std::cmp::Ordering;
use std::io::BufRead;
use std::path::PathBuf;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Method {
    Zscore,
    // Dbscan,
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
    let file = crate::utils::reader(&combined.into_os_string().into_string().unwrap());
    let mut lines = file.lines();
    let line = lines.next().unwrap().unwrap();
    println!("chrom\tbegin\tend\toutliers");
    let samples: Vec<&str> = line.split('\t').skip(3).collect();
    for line in lines {
        let line = line.unwrap();
        let splitline = line.split('\t').collect::<Vec<&str>>();
        let chrom = &splitline[0];
        let begin = &splitline[1];
        let end = &splitline[2];
        let values: Vec<f32> = splitline
            .iter()
            .skip(3)
            .map(|number| number.parse().unwrap())
            .filter(|v: &f32| !v.is_nan())
            .collect();
        // If all values are NaN the vector is empty
        if values.is_empty() {
            continue;
        }
        if values
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Less))
            .unwrap()
            < &(minsize as f32)
        {
            continue;
        }
        match method {
            Method::Zscore => {
                // calculate mean and std deviation of the STR lengths
                let (values_mean, values_std_dev) = std_deviation_and_mean(&values);
                // calculate the zscore for each haplotype and get the haplotype identifier based on the index if larger zscore > 3.0
                let expanded = values
                    .iter()
                    .enumerate()
                    .filter(|(_, &value)| ((value - values_mean) / values_std_dev) > zscore_cutoff)
                    .map(|(index, _)| samples[index])
                    .collect::<Vec<&str>>();
                if !expanded.is_empty() {
                    let expanded = expanded.join(",");
                    println!("{chrom}\t{begin}\t{end}\t{expanded}")
                }
            }
            // Method::Dbscan => {
            //     let arr = array![values];
            //     let clusters = Dbscan::params(3).tolerance(1e-2).transform(&arr);
            // }
        }
    }
}
