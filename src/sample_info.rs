//! Sample/phenotype information parsing
//!
//! This module provides functionality to parse sample metadata files
//! containing phenotype/grouping information for samples.
//! This is distinct from file header metadata (lines starting with #).

use std::io::BufRead;
use std::path::Path;

/// Individual sample with identifier and group assignment
pub struct Individual {
    pub identifier: String,
    pub group: String,
}

/// Parse sample phenotypes from a sample metadata file
///
/// Reads a TSV file with sample identifiers and phenotype information,
/// filtering for samples matching the specified condition.
pub fn parse_phenotypes(
    sample_metadata: &Path,
    condition: &str,
) -> Result<Vec<Individual>, Box<dyn std::error::Error>> {
    let pheno_column = condition
        .split(':')
        .next()
        .expect("Issue parsing condition string");
    let pheno_values = condition
        .split(':')
        .nth(1)
        .unwrap()
        .split(',')
        .collect::<Vec<&str>>();
    let meta_file = crate::utils::reader(sample_metadata.to_str().unwrap());
    let mut lines = meta_file.lines();
    let header = lines.next().unwrap().unwrap();
    let columns: Vec<&str> = header.split('\t').collect();
    let pheno_column_index = columns
        .iter()
        .position(|col| *col == pheno_column)
        .ok_or_else(|| {
            // A common slip is mixing '-' and '_' in the column name; suggest a
            // column that matches once both are normalised.
            let normalize = |s: &str| s.replace('-', "_");
            let suggestion = columns
                .iter()
                .find(|col| normalize(col) == normalize(pheno_column))
                .map(|col| format!(" Did you mean '{col}'?"))
                .unwrap_or_default();
            format!(
                "Could not find column '{}' in sample metadata file: {}.{}\nAvailable columns: {}",
                pheno_column,
                sample_metadata.display(),
                suggestion,
                columns.join(", ")
            )
        })?;
    let mut samples_of_interest: Vec<Individual> = vec![];
    for line in lines {
        let line = line.unwrap();
        let splitline = line.split('\t').collect::<Vec<&str>>();
        let pheno_value = splitline.get(pheno_column_index).unwrap();
        if pheno_values.contains(pheno_value) {
            samples_of_interest.push(Individual {
                identifier: splitline.first().unwrap().to_string(),
                group: pheno_value.to_string(),
            })
        }
    }
    Ok(samples_of_interest)
}
