use std::io::BufRead;
use std::path::PathBuf;

pub fn query(combined: PathBuf, region: String) {
    let file = crate::utils::reader(&combined.into_os_string().into_string().unwrap());
    let mut lines = file.lines();
    let header_line = lines.next().unwrap().unwrap();
    let samples = header_line.split('\t').skip(3);

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
                println!("{}:{}-{}", chrom, begin, end);
                let values = splitline
                    .iter()
                    .skip(3)
                    .map(|number| number.parse::<f64>().expect("Failed parsing lengths"));
                let mut zipped = samples.zip(values).into_iter().collect::<Vec<_>>();
                zipped.sort_by_key(
                    |&(_, val)| {
                        if !val.is_nan() {
                            -val as i64
                        } else {
                            i64::MAX
                        }
                    },
                );
                for (name, val) in zipped {
                    println!("{name}\t{val}");
                }
                break;
            }
        }
    }
}
