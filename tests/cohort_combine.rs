//! Incremental cohort building, end to end through the binary.
//!
//! Both forms of this workflow were silently broken: adding a sample to an existing cohort
//! failed with a header mismatch because the `#` metadata block was mistaken for the column
//! header, and merging two cohort files assembled the result in memory and then discarded it,
//! emitting a header-only file with a zero exit status.
//!
//! Driven through the binary rather than the library because `combine` writes to stdout, and
//! the test harness captures `println!` even from spawned threads.

use std::path::{Path, PathBuf};
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_inquiSTR");

/// A minimal individual-call file, metadata block included
fn individual(dir: &Path, sample: &str) -> PathBuf {
    let path = dir.join(format!("{}.inq", sample));
    std::fs::write(
        &path,
        format!(
            "# file_type=individual_call\n# version=test\n# command=call\n\
             chromosome\tbegin\tend\tinfo\t{s}_H1\t{s}_H2\n\
             chr1\t100\t130\tCAG\t10\t20\n\
             chr1\t200\t230\tCAG\t0\t5\n",
            s = sample
        ),
    )
    .unwrap();
    path
}

fn combine(inputs: &[&Path], out: &Path) {
    let output = Command::new(BIN)
        .arg("combine")
        .args(inputs)
        .output()
        .expect("failed to run inquiSTR");
    assert!(
        output.status.success(),
        "combine failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::write(out, &output.stdout).unwrap();
}

fn rows(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with("chr") && !l.starts_with("chromosome"))
        .map(str::to_string)
        .collect()
}

fn header(path: &Path) -> String {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .find(|l| l.starts_with("chromosome"))
        .expect("no column header")
        .to_string()
}

#[test]
fn cohort_expansion_and_merge_keep_their_data() {
    let dir = tempfile::tempdir().unwrap();
    let (a, b, c) = (
        individual(dir.path(), "sA"),
        individual(dir.path(), "sB"),
        individual(dir.path(), "sC"),
    );

    // Baseline: two individual files
    let cohort = dir.path().join("cohort.inq");
    combine(&[&a, &b], &cohort);
    assert_eq!(rows(&cohort).len(), 2, "combining individual files must emit every locus");

    // Adding a sample to an existing cohort
    let expanded = dir.path().join("expanded.inq");
    combine(&[&cohort, &c], &expanded);
    let h = header(&expanded);
    assert_eq!(h.matches("info").count(), 1, "`info` must appear exactly once");
    assert!(h.contains("sC_H1"), "the added sample must reach the output");
    let data = rows(&expanded);
    assert_eq!(data.len(), 2, "expansion must keep every locus");
    for row in &data {
        assert_eq!(
            row.split('\t').count(),
            h.split('\t').count(),
            "header and data must agree on field count"
        );
    }

    // Merging two cohorts, with no individual file to drive the writer
    let cohort2 = dir.path().join("cohort2.inq");
    combine(&[&c, &a], &cohort2);
    let merged = dir.path().join("merged.inq");
    combine(&[&cohort, &cohort2], &merged);
    assert_eq!(rows(&merged).len(), 2, "merging cohorts must not discard the data");
}
