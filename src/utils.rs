use niffler;
use std::fs::File;
use std::io::{BufReader, Read};

/// Read normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
pub fn reader(filename: &str) -> BufReader<Box<dyn Read>> {
    let (reader, _) = niffler::get_reader(Box::new(
        File::open(filename).expect("Problem opening file"),
    ))
    .expect("Problem reading file");
    BufReader::new(reader)
}

/// parse a region string
pub fn process_region(reg: String) -> Result<(String, u32, u32), Box<dyn std::error::Error>> {
    let reg = reg.replace(',', "");
    if reg.matches(':').count() != 1 {
        panic!(
            "\n\nError while parsing interval, could not find exactly one `:` character separating chromosome and start\nGot {reg}"
        );
    }
    if reg.matches('-').count() != 1 {
        panic!(
            "\n\nError while parsing interval, could not find exactly one `-` character separating start and end"
        );
    }
    let chrom = reg.split(':').collect::<Vec<&str>>()[0];
    let interval = reg.split(':').collect::<Vec<&str>>()[1];
    let start: u32 = interval.split('-').collect::<Vec<&str>>()[0]
        .parse()
        .expect("\n\nError while parsing interval start coordinate!\n\n");
    let end: u32 = interval.split('-').collect::<Vec<&str>>()[1]
        .parse()
        .expect("\n\nError while parsing interval end coordinate!\n\n");
    assert!(
        start < end,
        r#"\n\nInvalid region: start coordinate has to be smaller than end.\n\n"#
    );
    Ok((chrom.to_string(), start, end))
}
