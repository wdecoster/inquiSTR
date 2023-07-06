use log::debug;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;
use std::path::PathBuf;

pub fn query(combined: PathBuf, region: String) {
    let filename = combined
        .file_name()
        .expect("Problem getting filename")
        .to_str()
        .expect("Problem converting filename to string");
    let file = crate::utils::reader(filename);
    let mut lines = file.lines();
    let header_line = lines
        .next()
        .expect("Problem parsing file to get header")
        .expect("Problem parsing file to get header");
    let samples = header_line.split('\t').skip(3);
    debug!("Samples: {:?}", samples.clone().collect::<Vec<_>>());
    // check if the region variable is a file
    let intervals = match Path::new(&region).exists() {
        true => {
            let mut intervals = vec![];
            for line in crate::utils::reader(&region).lines() {
                intervals.push(crate::utils::process_region(line.unwrap()).unwrap());
            }
            intervals
        }
        false => vec![crate::utils::process_region(region).unwrap()],
    };
    // a vector in which the matching intervals will be stored
    let mut matching_intervals = vec![];
    let mut lengths = HashMap::new();
    for (chrom, reg_start, reg_end) in intervals {
        // note that we cannot make any assumptions about either the calls or the intervals to be sorted
        // so we have to search through the entire file for every interval
        // if this ever were to change things could be sped up
        debug!("Searching for region: {}:{}-{}", chrom, reg_start, reg_end);
        // Add a tab character to the chromosome so we can search for this with starts_with below (to make sure chr1 does not match chr15)
        let reg_chrom = format!("{chrom}\t");

        for line in crate::utils::reader(filename).lines() {
            let line = line.unwrap();
            debug!("Found line: {}", line);
            if line.starts_with(&reg_chrom) {
                debug!("Found right chromosome: {}", line);
                let splitline = line.split('\t').collect::<Vec<&str>>();
                let begin: u32 = splitline[1].parse().expect("Failed parsing interval");
                let end: u32 = splitline[2].parse().expect("Failed parsing interval");
                // test if there is an overlap between the region of interest and the region in the file
                if std::cmp::max(reg_start, begin) < std::cmp::min(reg_end, end) {
                    // store the matching interval to be used as header later
                    matching_intervals.push(format!("{chrom}:{begin}-{end}"));
                    let values = splitline
                        .iter()
                        .skip(3)
                        .map(|number| number.parse::<f64>().expect("Failed parsing lengths"));
                    for (sample, value) in samples.clone().zip(values) {
                        let entry = lengths.entry(sample).or_insert(vec![]);
                        entry.push(value);
                    }
                }
            }
        }
    }
    // if there was only a single matching interval, the output is printed sorted
    match matching_intervals.len() {
        0 => eprintln!("No matching intervals found in file"),
        1 => {
            // print the header
            println!("name\t{}", matching_intervals[0]);
            // sort the hashmap by value
            let mut zipped = lengths.clone().into_iter().collect::<Vec<_>>();

            zipped.sort_by_key(|(_, val)| {
                if !val[0].is_nan() {
                    // values that are not NaN are sorted by their value, descending
                    -val[0] as i64
                } else {
                    // NaN values are sorted last
                    i64::MAX
                }
            });
            for (name, val) in zipped {
                println!("{name}\t{length}", length = val[0]);
            }
        }
        _ => {
            // if there were multiple matching intervals, the output is printed in a tab separated table
            // print the header
            println!("name\t{}", matching_intervals.join("\t"));
            // for every sample in the hashmap, print all the values
            for (name, val) in lengths {
                println!(
                    "{name}\t{length}",
                    name = name,
                    length = val
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join("\t")
                );
            }
        }
    }
}
