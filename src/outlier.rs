use flate2::read;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::path::PathBuf;

/// Read normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(
            128 * 1024,
            read::GzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}

fn std_deviation_and_mean(data: &Vec<f32>) -> (f32, f32) {
    let sum = data.iter().sum::<f32>();
    let count = data.len() as f32;
    let data_mean = sum / count;
    let variance = data
        .iter()
        .map(|value| {
            let diff = data_mean - (*value as f32);

            diff * diff
        })
        .sum::<f32>()
        / count;

    (data_mean, variance.sqrt())
}

pub fn outlier(combined: PathBuf) {
    let file = reader(&combined.into_os_string().into_string().unwrap());
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
            .collect();
        // calculate mean and std deviation of the STR lengths
        let (values_mean, values_std_dev) = std_deviation_and_mean(&values);
        // calculate the zscore for each haplotype and store the index if larger than 3.0
        let expanded: Vec<&str> = values
            .iter()
            .enumerate()
            .filter(|(_, &value)| ((value - values_mean) / values_std_dev) > 3.0)
            .map(|(index, _)| samples[index])
            .collect();
        if !expanded.is_empty() {
            println!("{}\t{}\t{}\t{}", chrom, begin, end, expanded.join(","))
        }
    }
}
