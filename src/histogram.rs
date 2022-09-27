use histo_fp::Histogram;
use std::io::BufRead;
use std::path::PathBuf;

pub fn histogram(combined: PathBuf, region: String) {
    let file = crate::utils::reader(&combined.into_os_string().into_string().unwrap());
    let lines = file.lines();

    let (chrom, reg_start, reg_end) = crate::utils::process_region(region).unwrap();
    // Add a tab character to the chromosome so we can search for this with starts_with below (to make sure chr1 does not match chr15)
    let reg_chrom = format!("{}\t", chrom);

    for line in lines {
        let line = line.unwrap();
        if line.starts_with(&reg_chrom) {
            let splitline = line.split('\t').collect::<Vec<&str>>();
            let begin: u32 = splitline[1].parse().expect("Failed parsing interval");
            let end: u32 = splitline[2].parse().expect("Failed parsing interval");
            if reg_start <= begin && end <= reg_end {
                let mut histogram = Histogram::with_buckets(100, Some(2));
                for value in splitline
                    .iter()
                    .skip(3)
                    .map(|number| number.parse::<f64>().expect("Failed parsing lengths"))
                {
                    histogram.add(value);
                }
                println!("{}", histogram);
                break;
            }
        }
    }
}
