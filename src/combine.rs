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
    let file = match File::open(path) {
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

pub fn combine(calls: Vec<PathBuf>, unphased: bool) {
    // open the first file, regardless if it is gzipped or not
    let file1 = reader(&calls[0].clone().into_os_string().into_string().unwrap());

    // make a vector of readers over all other files
    let mut files: Vec<_> = calls[1..]
        .iter()
        .map(|file| reader(&file.clone().into_os_string().into_string().unwrap()).lines())
        .collect();
    // write out a header based on the filenames
    write_header(calls, unphased);
    for line in file1.lines() {
        // Get the full line for the first file, create a vector to collect all data
        let line = line.unwrap();
        let mut line_out = vec![line.as_str()];
        // Get one line from every other file
        let rest_of_files: Vec<String> = files
            .iter_mut()
            .map(|file2| file2.next().unwrap().unwrap())
            .collect();
        // Only get the scores, removing cols 1-3 for every other file
        let mut scores: Vec<&str> = rest_of_files
            .iter()
            .flat_map(|rec| rec.split('\t').skip(3))
            .collect();
        line_out.append(&mut scores);
        println!("{}", line_out.join("\t"));
    }
}

fn write_header(filenames: Vec<PathBuf>, unphased: bool) {
    // make the header vector
    let mut header = vec!["chrom".to_string(), "begin".to_string(), "end".to_string()];
    // take all filenames without the last extension (file_stem())
    if unphased {
        for name in filenames
            .iter()
            .map(|f| f.file_stem().unwrap().to_str().unwrap().to_owned())
        {
            header.push(name);
        }
    } else {
        // push twice to the header for every sample if phased
        for name in filenames
            .iter()
            .map(|f| f.file_stem().unwrap().to_str().unwrap().to_owned())
        {
            header.push(format!("{name}_H1"));
            header.push(format!("{name}_H2"));
        }
    }
    println!("{}", header.join("\t"));
}

#[cfg(test)]
#[test]
fn test_combine() {
    combine(
        vec![
            PathBuf::from("test-data/file1.inq"),
            PathBuf::from("test-data/file2.inq"),
            PathBuf::from("test-data/file3.inq"),
        ],
        false,
    );
}

#[test]
fn test_combine_gzipped() {
    combine(
        vec![
            PathBuf::from("test-data/file1.inq.gz"),
            PathBuf::from("test-data/file2.inq.gz"),
            PathBuf::from("test-data/file3.inq.gz"),
        ],
        false,
    );
}
