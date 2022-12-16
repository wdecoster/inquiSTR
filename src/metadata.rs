use std::io::BufRead;
use std::path::Path;

pub struct Individual {
    pub identifier: String,
    pub group: String,
}

pub fn parse_phenotypes(
    metadata: &Path,
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
    let meta_file = crate::utils::reader(metadata.to_str().unwrap());
    let mut lines = meta_file.lines();
    let header = lines.next().unwrap().unwrap();
    let pheno_column_index = header
        .split('\t')
        .enumerate()
        .filter(|(_, col)| col == &pheno_column)
        .map(|(index, _)| index)
        .next()
        .unwrap_or_else(|| {
            panic!(
                "Could not find column {} in {}",
                pheno_column,
                metadata.display()
            )
        });
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
