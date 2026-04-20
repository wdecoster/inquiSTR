//! # Relatedness Computation for STR Genotypes
//!
//! This module implements kinship computation between samples using STR genotype data.
//! It calculates the kinship coefficient between all pairs of samples, using an
//! IBS-based estimator that is robust to population structure.
//!
//! ## Method
//!
//! The kinship coefficient is computed using frequency-corrected IBS (Identity-by-State)
//! with per-pair baseline estimation for robustness to population structure:
//!
//! - For each pair (i,j), the expected IBS/2 for unrelated individuals is estimated as
//!   E0 = (hom_rate_i + hom_rate_j) / 2, where hom_rate is the fraction of loci at which
//!   a sample is homozygous. This adapts to each pair's population background.
//! - Per-locus homozygosity h_l (fraction of homozygous samples) provides locus-specific
//!   scaling in the denominator.
//! - kinship = Σ(IBS/2 - E0) / Σ(1 + h_l - 2*E0) / 2
//!
//! Expected values:
//! - Identical twins / duplicates: ~0.5
//! - Parent-child / full siblings: ~0.25
//! - Half-siblings / grandparent-grandchild: ~0.125
//! - First cousins: ~0.0625
//! - Unrelated: ~0.0 (or slightly negative)
//!
//! ## Output Format
//!
//! TSV file with columns:
//! - sample1: First sample name
//! - sample2: Second sample name
//! - kinship: Kinship coefficient (~0 for unrelated, ~0.25 for parent-child)
//! - n_loci: Number of informative loci used
//! - ibs0: Number of loci with 0 shared alleles
//! - ibs1: Number of loci with 1 shared allele
//! - ibs2: Number of loci with 2 shared alleles

use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::PathBuf;

/// Type alias for genotype matrix: Vec<Vec<(H1, H2)>>
type GenotypeMatrix = Vec<Vec<(f64, f64)>>;

/// Compute relatedness between all pairs of samples
pub fn compute_relatedness(
    combined: PathBuf,
    output: PathBuf,
    threads: usize,
    min_spacing: u32,
    tolerance: u32,
) {
    // Set number of threads
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap_or_else(|_| {
            eprintln!("Warning: Failed to set thread count, using default");
        });

    eprintln!("Reading combined STR file...");
    let (genotype_matrix, sample_names, locus_names) = parse_combined_file(&combined);
    let (genotype_matrix, locus_names) =
        apply_min_spacing_filter(&genotype_matrix, &locus_names, min_spacing);
    let locus_hom = compute_locus_homozygosity(&genotype_matrix);
    let sample_hom_rates = compute_sample_hom_rates(&genotype_matrix);

    let n_samples = sample_names.len();
    let n_loci = locus_names.len();

    eprintln!("Computing relatedness for {} samples across {} loci...", n_samples, n_loci);

    // Generate all sample pairs (excluding self-comparisons)
    let mut pairs = Vec::new();
    for i in 0..n_samples {
        for j in (i + 1)..n_samples {
            pairs.push((i, j));
        }
    }

    eprintln!("Computing {} pairwise comparisons...", pairs.len());

    // Compute relatedness for all pairs in parallel
    let results: Vec<_> = pairs
        .par_iter()
        .map(|(i, j)| {
            let (kinship, n_loci, ibs0, ibs1, ibs2) = compute_pairwise_relatedness(
                &genotype_matrix,
                &locus_hom,
                sample_hom_rates[*i],
                sample_hom_rates[*j],
                *i,
                *j,
                tolerance,
            );
            (
                sample_names[*i].clone(),
                sample_names[*j].clone(),
                kinship,
                n_loci,
                ibs0,
                ibs1,
                ibs2,
            )
        })
        .collect();

    // Write output
    eprintln!("Writing results to {}...", output.display());
    let mut out_file = File::create(&output).expect("Failed to create output file");

    // Write header
    writeln!(out_file, "sample1\tsample2\tkinship\tn_loci\tibs0\tibs1\tibs2")
        .expect("Failed to write header");

    // Write results sorted by relatedness (descending)
    let mut sorted_results = results;
    sorted_results.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

    for (sample1, sample2, relatedness, n_loci, ibs0, ibs1, ibs2) in sorted_results {
        writeln!(
            out_file,
            "{}\t{}\t{:.6}\t{}\t{}\t{}\t{}",
            sample1, sample2, relatedness, n_loci, ibs0, ibs1, ibs2
        )
        .expect("Failed to write result");
    }

    eprintln!("Done!");
}

/// Parse combined file and extract genotype matrix
/// Returns: (genotype_matrix[locus][sample][haplotype], sample_names, locus_names)
fn parse_combined_file(combined: &std::path::Path) -> (GenotypeMatrix, Vec<String>, Vec<String>) {
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();

    // Skip metadata lines and get header
    let header_line = crate::utils::skip_metadata_lines(&mut lines, &combined.to_string_lossy());
    let header_fields: Vec<&str> = header_line.trim().split('\t').collect();

    // Use centralized validation function
    let n_samples = match crate::filetype::validate_str_header(&header_fields) {
        Ok(n) => n,
        Err(e) => {
            eprintln!("ERROR: Invalid combined file header format");
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };

    if n_samples < 2 {
        eprintln!("ERROR: Need at least 2 samples for relatedness computation");
        std::process::exit(1);
    }

    // Extract sample names from validated header
    let sample_names: Vec<String> = (0..n_samples)
        .map(|i| {
            let h1_col = header_fields[crate::filetype::STR_FIXED_COLUMNS + i * 2];
            h1_col.trim_end_matches("_H1").to_string()
        })
        .collect();

    eprintln!("Found {} samples", n_samples);

    // Parse data lines
    let mut locus_names = Vec::new();
    let mut genotype_matrix = Vec::new();

    for (line_num, line_result) in lines.enumerate() {
        let line = line_result.expect("Error reading line");
        let fields: Vec<&str> = line.trim().split('\t').collect();

        let expected_fields = crate::filetype::STR_FIXED_COLUMNS + n_samples * 2; // chr, start, end, info + H1/H2 pairs
        if fields.len() != expected_fields {
            eprintln!(
                "ERROR: Malformed line {} in combined file: expected {} fields, got {}",
                line_num + 2, // +2 because we skipped header line and lines are 0-indexed
                expected_fields,
                fields.len()
            );
            eprintln!(
                "Line content: {}:{}-{}",
                fields.first().unwrap_or(&"?"),
                fields.get(1).unwrap_or(&"?"),
                fields.get(2).unwrap_or(&"?")
            );
            std::process::exit(1);
        }

        // Store locus identifier
        let locus_id = format!("{}:{}-{}", fields[0], fields[1], fields[2]);
        locus_names.push(locus_id);

        // Parse genotypes for this locus (all samples)
        let mut locus_genotypes = Vec::new();
        for i in 0..n_samples {
            let h1_idx = crate::filetype::STR_FIXED_COLUMNS + i * 2;
            let h2_idx = crate::filetype::STR_FIXED_COLUMNS + i * 2 + 1;

            let h1 = fields
                .get(h1_idx)
                .and_then(|s| s.parse::<f64>().ok())
                .unwrap_or(f64::NAN);
            let h2 = fields
                .get(h2_idx)
                .and_then(|s| s.parse::<f64>().ok())
                .unwrap_or(f64::NAN);

            locus_genotypes.push((h1, h2));
        }
        genotype_matrix.push(locus_genotypes);
    }

    eprintln!("Parsed {} loci", genotype_matrix.len());

    (genotype_matrix, sample_names, locus_names)
}

/// Thin loci based on minimum spacing (bp) along each chromosome.
///
/// This is a simple greedy spacing filter intended to reduce LD inflation in relatedness.
fn apply_min_spacing_filter(
    genotype_matrix: &GenotypeMatrix,
    locus_names: &[String],
    min_spacing: u32,
) -> (GenotypeMatrix, Vec<String>) {
    if locus_names.is_empty() || min_spacing == 0 {
        return (genotype_matrix.clone(), locus_names.to_vec());
    }

    let mut filtered_genotypes = Vec::new();
    let mut filtered_loci = Vec::new();

    let mut last_chr = String::new();
    let mut last_pos = 0u32;
    let mut first = true;

    for (locus_idx, locus_name) in locus_names.iter().enumerate() {
        // locus_name format: chr:start-end
        let mut parts = locus_name.split(':');
        let chr = match parts.next() {
            Some(c) => c,
            None => continue,
        };
        let range = match parts.next() {
            Some(r) => r,
            None => continue,
        };
        let start = match range.split('-').next().and_then(|s| s.parse::<u32>().ok()) {
            Some(s) => s,
            None => continue,
        };

        let keep = if first || chr != last_chr {
            true
        } else {
            start.saturating_sub(last_pos) >= min_spacing
        };

        if keep {
            filtered_loci.push(locus_name.clone());
            filtered_genotypes.push(genotype_matrix[locus_idx].clone());
            last_chr = chr.to_string();
            last_pos = start;
            first = false;
        }
    }

    let dropped = locus_names.len() - filtered_loci.len();
    if dropped > 0 {
        eprintln!("Thinned loci: dropped {} loci with min_spacing={} bp", dropped, min_spacing);
    }

    (filtered_genotypes, filtered_loci)
}

/// Compute per-locus homozygosity fraction (fraction of samples that are homozygous).
///
/// Uses exact allele comparison (no tolerance) since this measures locus-level
/// polymorphism, not cross-sample matching.
fn compute_locus_homozygosity(genotype_matrix: &GenotypeMatrix) -> Vec<f64> {
    genotype_matrix
        .iter()
        .map(|locus_genotypes| {
            let mut n_hom = 0usize;
            let mut n_valid = 0usize;
            for &(h1, h2) in locus_genotypes {
                if h1.is_nan() || h2.is_nan() {
                    continue;
                }
                n_valid += 1;
                if h1.round() as i64 == h2.round() as i64 {
                    n_hom += 1;
                }
            }
            if n_valid == 0 {
                1.0
            } else {
                n_hom as f64 / n_valid as f64
            }
        })
        .collect()
}

/// Compute per-sample genome-wide homozygosity rate.
///
/// For each sample, returns the fraction of non-missing loci where the sample is
/// homozygous. Uses exact allele comparison (no tolerance).
fn compute_sample_hom_rates(genotype_matrix: &GenotypeMatrix) -> Vec<f64> {
    if genotype_matrix.is_empty() {
        return Vec::new();
    }
    let n_samples = genotype_matrix[0].len();
    let mut n_hom = vec![0usize; n_samples];
    let mut n_valid = vec![0usize; n_samples];

    for locus_genotypes in genotype_matrix {
        for (idx, &(h1, h2)) in locus_genotypes.iter().enumerate() {
            if h1.is_nan() || h2.is_nan() {
                continue;
            }
            n_valid[idx] += 1;
            if h1.round() as i64 == h2.round() as i64 {
                n_hom[idx] += 1;
            }
        }
    }

    n_hom
        .iter()
        .zip(n_valid.iter())
        .map(|(&h, &v)| if v == 0 { 1.0 } else { h as f64 / v as f64 })
        .collect()
}

/// Compute pairwise kinship between two samples.
///
/// Uses frequency-corrected IBS with per-pair baseline estimation:
///   E0 = (hom_rate_i + hom_rate_j) / 2
///   kinship = Σ(IBS/2 - E0) / Σ(1 + h_l - 2*E0) / 2
///
/// Returns: (kinship, n_informative_loci, ibs0, ibs1, ibs2)
fn compute_pairwise_relatedness(
    genotype_matrix: &[Vec<(f64, f64)>],
    locus_hom: &[f64],
    hom_rate_i: f64,
    hom_rate_j: f64,
    sample_i: usize,
    sample_j: usize,
    tolerance: u32,
) -> (f64, usize, usize, usize, usize) {
    let mut ibs0 = 0usize;
    let mut ibs1 = 0usize;
    let mut ibs2 = 0usize;
    let mut n_informative = 0usize;
    let mut numerator = 0.0;
    let mut denominator = 0.0;

    let e0 = (hom_rate_i + hom_rate_j) / 2.0;
    let within_tolerance = |a: i64, b: i64| (a - b).abs() <= tolerance as i64;

    for (locus_idx, locus_genotypes) in genotype_matrix.iter().enumerate() {
        let (h1_i, h2_i) = locus_genotypes[sample_i];
        let (h1_j, h2_j) = locus_genotypes[sample_j];

        // Skip loci where either sample has missing data
        if h1_i.is_nan() || h2_i.is_nan() || h1_j.is_nan() || h2_j.is_nan() {
            continue;
        }

        n_informative += 1;

        let h1_i_r = h1_i.round() as i64;
        let h2_i_r = h2_i.round() as i64;
        let h1_j_r = h1_j.round() as i64;
        let h2_j_r = h2_j.round() as i64;

        // Count matching alleles (IBS)
        let matches = if (within_tolerance(h1_i_r, h1_j_r) && within_tolerance(h2_i_r, h2_j_r))
            || (within_tolerance(h1_i_r, h2_j_r) && within_tolerance(h2_i_r, h1_j_r))
        {
            2
        } else if within_tolerance(h1_i_r, h1_j_r)
            || within_tolerance(h1_i_r, h2_j_r)
            || within_tolerance(h2_i_r, h1_j_r)
            || within_tolerance(h2_i_r, h2_j_r)
        {
            1
        } else {
            0
        };

        match matches {
            0 => ibs0 += 1,
            1 => ibs1 += 1,
            2 => ibs2 += 1,
            _ => unreachable!(),
        }

        let ibs_half = matches as f64 / 2.0;
        let h_l = locus_hom[locus_idx];
        numerator += ibs_half - e0;
        denominator += 1.0 + h_l - 2.0 * e0;
    }

    if n_informative == 0 || denominator == 0.0 {
        return (f64::NAN, 0, 0, 0, 0);
    }

    // Kinship = relatedness / 2
    // relatedness = Σ(IBS/2 - E0) / Σ(1 + h - 2*E0) gives ~0.5 for parent-child
    // dividing by 2 gives kinship ~0.25 for parent-child
    let kinship = numerator / denominator / 2.0;

    (kinship, n_informative, ibs0, ibs1, ibs2)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: compute kinship from test data with tolerance=0
    fn test_relatedness(
        genotype_matrix: &GenotypeMatrix,
        i: usize,
        j: usize,
    ) -> (f64, usize, usize, usize, usize) {
        test_relatedness_with_tolerance(genotype_matrix, i, j, 0)
    }

    /// Helper: compute kinship with explicit tolerance.
    fn test_relatedness_with_tolerance(
        genotype_matrix: &GenotypeMatrix,
        i: usize,
        j: usize,
        tolerance: u32,
    ) -> (f64, usize, usize, usize, usize) {
        let locus_hom = compute_locus_homozygosity(genotype_matrix);
        let sample_hom_rates = compute_sample_hom_rates(genotype_matrix);
        compute_pairwise_relatedness(
            genotype_matrix,
            &locus_hom,
            sample_hom_rates[i],
            sample_hom_rates[j],
            i,
            j,
            tolerance,
        )
    }

    #[test]
    fn test_identical_samples() {
        // Two identical heterozygous samples should have kinship ~0.5
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (10.0, 12.0)],
            vec![(15.0, 17.0), (15.0, 17.0)],
            vec![(8.0, 9.0), (8.0, 9.0)],
        ];

        let (kinship, n_loci, _ibs0, _ibs1, ibs2) = test_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 3);
        assert_eq!(ibs2, 3);
        assert!(kinship > 0.4, "Identical samples should have kinship > 0.4, got {}", kinship);
    }

    #[test]
    fn test_unrelated_samples() {
        // Two completely different samples with no shared alleles
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (20.0, 22.0)],
            vec![(15.0, 17.0), (25.0, 27.0)],
            vec![(8.0, 9.0), (18.0, 19.0)],
        ];

        let (kinship, n_loci, ibs0, _ibs1, _ibs2) = test_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 3);
        assert_eq!(ibs0, 3);
        assert!(
            kinship < 0.1,
            "Completely unrelated samples should have low kinship, got {}",
            kinship
        );
    }

    #[test]
    fn test_ibs_counts_partial_sharing() {
        // Samples sharing one allele per locus
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (10.0, 22.0)], // Share first allele
            vec![(15.0, 17.0), (25.0, 17.0)], // Share second allele
        ];

        let (_kinship, n_loci, _ibs0, ibs1, _ibs2) = test_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 2);
        assert_eq!(ibs1, 2);
    }

    #[test]
    fn test_ibs_counts_parent_child() {
        // Parent-child: child inherits one allele from parent at each locus
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (10.0, 14.0)],
            vec![(15.0, 17.0), (17.0, 19.0)],
            vec![(8.0, 8.0), (8.0, 10.0)],
            vec![(20.0, 22.0), (20.0, 24.0)],
        ];

        let (_kinship, n_loci, _ibs0, ibs1, _ibs2) = test_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 4);
        assert_eq!(ibs1, 4);
    }

    #[test]
    fn test_ibs_counts_full_siblings() {
        // Full siblings: mix of IBS0, IBS1, IBS2
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (10.0, 12.0)], // IBS2
            vec![(15.0, 17.0), (15.0, 19.0)], // IBS1
            vec![(8.0, 9.0), (10.0, 11.0)],   // IBS0
            vec![(20.0, 22.0), (20.0, 24.0)], // IBS1
        ];

        let (_kinship, n_loci, ibs0, ibs1, ibs2) = test_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 4);
        assert_eq!(ibs2, 1);
        assert_eq!(ibs1, 2);
        assert_eq!(ibs0, 1);
    }

    #[test]
    fn test_ibs_counts_half_siblings() {
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (10.0, 14.0)], // IBS1
            vec![(15.0, 17.0), (19.0, 21.0)], // IBS0
            vec![(8.0, 9.0), (8.0, 11.0)],    // IBS1
            vec![(20.0, 22.0), (24.0, 26.0)], // IBS0
        ];

        let (_kinship, n_loci, ibs0, ibs1, _ibs2) = test_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 4);
        assert_eq!(ibs1, 2);
        assert_eq!(ibs0, 2);
    }

    #[test]
    fn test_ibs_counts_first_cousins() {
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (10.0, 14.0)], // IBS1
            vec![(15.0, 17.0), (19.0, 21.0)], // IBS0
            vec![(8.0, 9.0), (11.0, 13.0)],   // IBS0
            vec![(20.0, 22.0), (24.0, 26.0)], // IBS0
            vec![(30.0, 32.0), (34.0, 36.0)], // IBS0
            vec![(40.0, 42.0), (44.0, 46.0)], // IBS0
            vec![(50.0, 52.0), (54.0, 56.0)], // IBS0
            vec![(60.0, 62.0), (60.0, 64.0)], // IBS1
        ];

        let (_kinship, n_loci, ibs0, ibs1, _ibs2) = test_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 8);
        assert_eq!(ibs1, 2);
        assert_eq!(ibs0, 6);
    }

    #[test]
    fn test_missing_data() {
        // Test handling of missing data (NaN values)
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (10.0, 12.0)],       // Valid, IBS2
            vec![(f64::NAN, 15.0), (15.0, 17.0)],   // Missing in sample 1, skip
            vec![(8.0, 9.0), (f64::NAN, f64::NAN)], // Missing in sample 2, skip
            vec![(20.0, 22.0), (20.0, 24.0)],       // Valid, IBS1
        ];

        let (_kinship, n_loci, _ibs0, ibs1, ibs2) = test_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 2); // Only 2 loci are informative
        assert_eq!(ibs2, 1);
        assert_eq!(ibs1, 1);
    }

    #[test]
    fn test_homozygous_loci() {
        // Test with homozygous genotypes
        let genotype_matrix = vec![
            vec![(10.0, 10.0), (10.0, 10.0)], // IBS2
            vec![(15.0, 15.0), (15.0, 17.0)], // IBS1
            vec![(20.0, 20.0), (22.0, 22.0)], // IBS0
        ];

        let (_kinship, n_loci, ibs0, ibs1, ibs2) = test_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 3);
        assert_eq!(ibs2, 1);
        assert_eq!(ibs1, 1);
        assert_eq!(ibs0, 1);
    }

    #[test]
    fn test_phased_vs_unphased() {
        // Test that phase doesn't matter for IBS calculation
        // (10, 12) should match (12, 10)
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (12.0, 10.0)], // Same genotype, different phase
            vec![(15.0, 17.0), (17.0, 15.0)], // Same genotype, different phase
        ];

        let (kinship, n_loci, _ibs0, _ibs1, ibs2) = test_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 2);
        assert_eq!(ibs2, 2); // Both should be IBS2
        assert!(kinship > 0.4, "Identical het samples should have high kinship, got {}", kinship);
    }

    #[test]
    fn test_tolerance_1bp() {
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (11.0, 12.0)], // differ by 1 on one allele
            vec![(15.0, 17.0), (16.0, 17.0)],
        ];

        let (kinship, n_loci, _ibs0, _ibs1, ibs2) =
            test_relatedness_with_tolerance(&genotype_matrix, 0, 1, 1);

        assert_eq!(n_loci, 2);
        assert_eq!(ibs2, 2);
        assert!(kinship > 0.0);
    }

    #[test]
    fn test_locus_homozygosity() {
        // Monomorphic: all hom
        let genotype_matrix = vec![vec![(10.0, 10.0), (10.0, 10.0), (10.0, 10.0)]];
        let hom = compute_locus_homozygosity(&genotype_matrix);
        assert!((hom[0] - 1.0).abs() < 1e-6);

        // All het
        let genotype_matrix = vec![vec![(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]];
        let hom = compute_locus_homozygosity(&genotype_matrix);
        assert!(hom[0].abs() < 1e-6);

        // Mixed: 1 hom out of 3
        let genotype_matrix = vec![vec![(1.0, 1.0), (3.0, 4.0), (5.0, 6.0)]];
        let hom = compute_locus_homozygosity(&genotype_matrix);
        assert!((hom[0] - 1.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_sample_hom_rates() {
        // Sample 0: 2/3 hom, Sample 1: 1/3 hom
        let genotype_matrix = vec![
            vec![(10.0, 10.0), (1.0, 2.0)],
            vec![(15.0, 15.0), (3.0, 3.0)],
            vec![(8.0, 9.0), (5.0, 6.0)],
        ];
        let rates = compute_sample_hom_rates(&genotype_matrix);
        assert!((rates[0] - 2.0 / 3.0).abs() < 1e-6);
        assert!((rates[1] - 1.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_min_spacing_filter() {
        let genotype_matrix = vec![vec![(1.0, 1.0)], vec![(1.0, 1.0)], vec![(1.0, 1.0)]];
        let locus_names = vec![
            "chr1:100-100".to_string(),
            "chr1:150-150".to_string(),
            "chr1:260-260".to_string(),
        ];

        let (filtered_genotypes, filtered_loci) =
            apply_min_spacing_filter(&genotype_matrix, &locus_names, 100);
        assert_eq!(filtered_loci, vec!["chr1:100-100".to_string(), "chr1:260-260".to_string()]);
        assert_eq!(filtered_genotypes.len(), 2);
    }

    #[test]
    fn test_relatedness_ordering() {
        // Verify that the estimator gives higher kinship for related pairs.
        // 6 samples, 10 loci. Each locus uses unique alleles per sample except:
        // Samples 4 and 5 always share exactly 1 allele (parent-child pattern).
        let mut genotype_matrix = Vec::new();
        for l in 0..10u64 {
            let base = l * 20;
            genotype_matrix.push(vec![
                ((base + 1) as f64, (base + 2) as f64),  // sample 0
                ((base + 3) as f64, (base + 4) as f64),  // sample 1
                ((base + 5) as f64, (base + 6) as f64),  // sample 2
                ((base + 7) as f64, (base + 8) as f64),  // sample 3
                ((base + 9) as f64, (base + 10) as f64), // sample 4 (parent)
                ((base + 9) as f64, (base + 11) as f64), // sample 5 (child)
            ]);
        }

        let locus_hom = compute_locus_homozygosity(&genotype_matrix);
        let sample_hom_rates = compute_sample_hom_rates(&genotype_matrix);

        let (kin_parent_child, _, _, ibs1_pc, _) = compute_pairwise_relatedness(
            &genotype_matrix,
            &locus_hom,
            sample_hom_rates[4],
            sample_hom_rates[5],
            4,
            5,
            0,
        );
        let (kin_unrelated, _, _, _, _) = compute_pairwise_relatedness(
            &genotype_matrix,
            &locus_hom,
            sample_hom_rates[0],
            sample_hom_rates[2],
            0,
            2,
            0,
        );

        assert_eq!(ibs1_pc, 10, "Parent-child should share 1 allele at all loci");
        assert!(
            kin_parent_child > kin_unrelated,
            "Parent-child ({:.4}) should exceed unrelated ({:.4})",
            kin_parent_child,
            kin_unrelated
        );
    }
}
