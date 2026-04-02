//! Convert VCF files to inquiSTR format
//!
//! This module provides functionality to convert VCF files (typically from TRGT or similar
//! STR genotyping tools) to inquiSTR's TSV format. It handles:
//! - Single sample VCF → individual call file
//! - Multiple sample VCFs or multiple single-sample VCFs → combined file
//! - Allele length calculation relative to reference
//! - Proper handling of missing genotypes

use crate::errors::InquiSTRError;
use crate::filetype::FileType;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// Represents a single STR variant from a VCF file
#[derive(Debug, Clone)]
struct VcfVariant {
    chromosome: String,
    pos: u32,
    end: u32,
    id: String,
    ref_len: usize,
    alt_seqs: Vec<String>,
    samples: Vec<Genotype>,
}

/// Represents a genotype for a single sample
#[derive(Debug, Clone)]
struct Genotype {
    sample_name: String,
    allele1_idx: Option<i32>,
    allele2_idx: Option<i32>,
}

impl Genotype {
    /// Calculate allele lengths relative to reference
    fn get_allele_lengths(&self, variant: &VcfVariant) -> (String, String) {
        let ref_len = variant.ref_len as i32;

        let allele1_len = match self.allele1_idx {
            None | Some(-1) => "NA".to_string(),
            Some(0) => "0".to_string(), // Reference allele
            Some(idx) => {
                let alt_idx = (idx - 1) as usize;
                if alt_idx < variant.alt_seqs.len() {
                    let alt_len = variant.alt_seqs[alt_idx].len() as i32;
                    (alt_len - ref_len).to_string()
                } else {
                    "NA".to_string()
                }
            }
        };

        let allele2_len = match self.allele2_idx {
            None | Some(-1) => "NA".to_string(),
            Some(0) => "0".to_string(), // Reference allele
            Some(idx) => {
                let alt_idx = (idx - 1) as usize;
                if alt_idx < variant.alt_seqs.len() {
                    let alt_len = variant.alt_seqs[alt_idx].len() as i32;
                    (alt_len - ref_len).to_string()
                } else {
                    "NA".to_string()
                }
            }
        };

        (allele1_len, allele2_len)
    }
}

/// Parse a VCF file and extract STR variants
fn parse_vcf(
    vcf_path: &Path,
    off_by_one: bool,
) -> Result<(Vec<String>, Vec<VcfVariant>), InquiSTRError> {
    let reader = crate::utils::reader(&vcf_path.to_string_lossy());
    let mut sample_names: Vec<String> = Vec::new();
    let mut variants: Vec<VcfVariant> = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| InquiSTRError {
            message: format!("Failed to read VCF file {}: {}", vcf_path.display(), e),
        })?;

        // Parse header to get sample names
        if line.starts_with("#CHROM") {
            let fields: Vec<&str> = line.split('\t').collect();
            // Sample names start from column 9 (0-indexed)
            if fields.len() > 9 {
                sample_names = fields[9..].iter().map(|s| s.to_string()).collect();
            }
            continue;
        }

        // Skip other header lines
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue; // Skip lines without genotype data
        }

        let chromosome = fields[0].to_string();
        let pos: u32 = fields[1].parse().map_err(|_| InquiSTRError {
            message: format!("Invalid position in VCF: {}", fields[1]),
        })?;
        let id = if fields[2] == "." {
            format!("{}:{}", chromosome, pos)
        } else {
            fields[2].to_string()
        };
        let ref_seq = fields[3].to_string();
        let alt_field = fields[4];

        // Parse ALT alleles
        let alt_seqs: Vec<String> = alt_field.split(',').map(|s| s.to_string()).collect();

        // Get END from INFO field if present
        let info = fields[7];
        let end = if let Some(end_str) = info.split(';').find(|s| s.starts_with("END=")) {
            end_str
                .trim_start_matches("END=")
                .parse()
                .unwrap_or(pos + ref_seq.len() as u32)
        } else {
            pos + ref_seq.len() as u32
        };

        // Parse FORMAT field to find GT index
        let format = fields[8];
        let format_fields: Vec<&str> = format.split(':').collect();
        let gt_idx = format_fields.iter().position(|&s| s == "GT");

        if gt_idx.is_none() {
            continue; // Skip if no genotype information
        }
        let gt_idx = gt_idx.unwrap();

        // Parse genotypes for each sample
        let mut genotypes = Vec::new();
        for (i, sample_field) in fields[9..].iter().enumerate() {
            let sample_name = if i < sample_names.len() {
                sample_names[i].clone()
            } else {
                format!("Sample{}", i + 1)
            };

            let sample_fields: Vec<&str> = sample_field.split(':').collect();
            if gt_idx >= sample_fields.len() {
                continue;
            }

            let gt = sample_fields[gt_idx];
            let (allele1_idx, allele2_idx) = parse_genotype(gt);

            genotypes.push(Genotype { sample_name, allele1_idx, allele2_idx });
        }

        variants.push(VcfVariant {
            chromosome,
            pos: if off_by_one { pos } else { pos - 1 },
            end,
            id,
            ref_len: ref_seq.len(),
            alt_seqs,
            samples: genotypes,
        });
    }

    if sample_names.is_empty() {
        return Err(InquiSTRError {
            message: format!("No samples found in VCF file: {}", vcf_path.display()),
        });
    }

    Ok((sample_names, variants))
}

/// Parse a genotype string (e.g., "0/1", "1|2", "./.")
fn parse_genotype(gt: &str) -> (Option<i32>, Option<i32>) {
    let separator = if gt.contains('|') { '|' } else { '/' };
    let alleles: Vec<&str> = gt.split(separator).collect();

    if alleles.len() != 2 {
        return (None, None);
    }

    let allele1 = if alleles[0] == "." {
        None
    } else {
        alleles[0].parse::<i32>().ok()
    };
    let allele2 = if alleles[1] == "." {
        None
    } else {
        alleles[1].parse::<i32>().ok()
    };

    (allele1, allele2)
}

/// Convert VCF files to inquiSTR format
pub fn convert_vcf(vcf_files: Vec<PathBuf>, off_by_one: bool) -> Result<(), InquiSTRError> {
    if vcf_files.is_empty() {
        return Err(InquiSTRError { message: "No VCF files provided".to_string() });
    }

    // Parse all VCF files
    let mut all_samples: Vec<String> = Vec::new();
    let mut all_variants_by_locus: HashMap<String, Vec<VcfVariant>> = HashMap::new();
    let mut sample_counts: HashMap<String, usize> = HashMap::new();

    for vcf_path in &vcf_files {
        eprintln!("Processing VCF file: {}", vcf_path.display());
        let (sample_names, mut variants) = parse_vcf(vcf_path, off_by_one)?;

        // Handle duplicate sample names by appending a suffix
        for sample in &sample_names {
            let count = sample_counts.entry(sample.clone()).or_insert(0);
            *count += 1;

            let final_sample_name = if *count > 1 {
                let suffix_name = format!("{}_{}", sample, count);
                eprintln!(
                    "Warning: Sample '{}' appears in multiple VCF files. Renaming to '{}' for file: {}",
                    sample,
                    suffix_name,
                    vcf_path.display()
                );
                suffix_name
            } else {
                sample.clone()
            };

            // Update sample names in variants
            for variant in &mut variants {
                for genotype in &mut variant.samples {
                    if &genotype.sample_name == sample {
                        genotype.sample_name = final_sample_name.clone();
                    }
                }
            }

            all_samples.push(final_sample_name);
        }

        // Group variants by locus
        for variant in variants {
            let locus_key = format!("{}:{}:{}", variant.chromosome, variant.pos, variant.end);
            all_variants_by_locus
                .entry(locus_key)
                .or_default()
                .push(variant);
        }
    }

    let is_single_sample = all_samples.len() == 1;
    let file_type = if is_single_sample {
        FileType::IndividualCall
    } else {
        FileType::CombinedCall
    };

    // Write metadata header
    println!(
        "# file_type={}",
        match file_type {
            FileType::IndividualCall => "individual_call",
            FileType::CombinedCall => "combined_call",
            _ => unreachable!(),
        }
    );
    println!("# source=inquiSTR convert");
    for (i, vcf_path) in vcf_files.iter().enumerate() {
        println!("# input_vcf_{}={}", i + 1, vcf_path.display());
    }
    if !all_samples.is_empty() {
        println!("# samples={}", all_samples.join(","));
    }

    // Write column headers
    print!("chromosome\tbegin\tend\tinfo");
    for sample in &all_samples {
        print!("\t{}_H1\t{}_H2", sample, sample);
    }
    println!();

    // Sort loci for consistent output
    let mut locus_keys: Vec<String> = all_variants_by_locus.keys().cloned().collect();
    locus_keys.sort_by(|a, b| {
        let parse_locus = |s: &str| -> (String, u32, u32) {
            let parts: Vec<&str> = s.split(':').collect();
            (
                parts[0].to_string(),
                parts[1].parse().unwrap_or(0),
                parts[2].parse().unwrap_or(0),
            )
        };
        let (chr_a, pos_a, end_a) = parse_locus(a);
        let (chr_b, pos_b, end_b) = parse_locus(b);
        crate::batching::human_cmp(&chr_a, &chr_b)
            .then(pos_a.cmp(&pos_b))
            .then(end_a.cmp(&end_b))
    });

    // Write variant data
    for locus_key in locus_keys {
        let variants = &all_variants_by_locus[&locus_key];
        if variants.is_empty() {
            continue;
        }

        // Use the first variant for locus information
        let first_variant = &variants[0];
        print!(
            "{}\t{}\t{}\t{}",
            first_variant.chromosome, first_variant.pos, first_variant.end, first_variant.id
        );

        // Merge genotypes from all variants at this locus
        let mut genotypes_by_sample: HashMap<String, (String, String)> = HashMap::new();
        for variant in variants {
            for genotype in &variant.samples {
                let (h1, h2) = genotype.get_allele_lengths(variant);
                genotypes_by_sample.insert(genotype.sample_name.clone(), (h1, h2));
            }
        }

        // Write genotypes in sample order
        // Missing variants (locus not present in a sample's VCF) are filled with NA
        for sample in &all_samples {
            if let Some((h1, h2)) = genotypes_by_sample.get(sample) {
                print!("\t{}\t{}", h1, h2);
            } else {
                // Sample doesn't have this locus (variant not in their VCF)
                print!("\tNA\tNA");
            }
        }
        println!();
    }

    eprintln!(
        "Conversion complete: {} loci, {} sample(s)",
        all_variants_by_locus.len(),
        all_samples.len()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_genotype() {
        assert_eq!(parse_genotype("0/0"), (Some(0), Some(0)));
        assert_eq!(parse_genotype("0/1"), (Some(0), Some(1)));
        assert_eq!(parse_genotype("1|2"), (Some(1), Some(2)));
        assert_eq!(parse_genotype("./."), (None, None));
        assert_eq!(parse_genotype(".|1"), (None, Some(1)));
    }

    #[test]
    fn test_get_allele_lengths() {
        let variant = VcfVariant {
            chromosome: "chr1".to_string(),
            pos: 1000,
            end: 1050,
            id: "test".to_string(),
            ref_len: 50,
            alt_seqs: vec!["ACGTACGT".to_string(), "ACGTACGTACGT".to_string()],
            samples: vec![],
        };

        let gt = Genotype {
            sample_name: "Sample1".to_string(),
            allele1_idx: Some(0), // Reference
            allele2_idx: Some(1), // First ALT
        };

        let (h1, h2) = gt.get_allele_lengths(&variant);
        assert_eq!(h1, "0");
        assert_eq!(h2, "-42"); // 8 - 50 = -42

        let gt_missing = Genotype {
            sample_name: "Sample2".to_string(),
            allele1_idx: None,
            allele2_idx: Some(2), // Second ALT
        };

        let (h1, h2) = gt_missing.get_allele_lengths(&variant);
        assert_eq!(h1, "NA");
        assert_eq!(h2, "-38"); // 12 - 50 = -38
    }
}
