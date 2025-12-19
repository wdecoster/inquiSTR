//! # Relatedness Computation for STR Genotypes
//!
//! This module implements relatedness computation between samples using STR genotype data.
//! It calculates the coefficient of relationship (kinship coefficient) between all pairs
//! of samples, similar to tools like somalier.
//!
//! ## Method
//!
//! The relatedness coefficient is computed using the IBS (Identity-by-State) method:
//! - For each locus, count the number of matching alleles between two samples
//! - Average across all informative loci (where both samples have valid calls)
//! - The coefficient ranges from 0 (unrelated) to 1 (identical)
//!
//! Expected values:
//! - Identical twins / duplicates: ~1.0
//! - Parent-child / full siblings: ~0.5
//! - Half-siblings / grandparent-grandchild: ~0.25
//! - First cousins: ~0.125
//! - Unrelated: ~0.0
//!
//! ## Output Format
//!
//! TSV file with columns:
//! - sample1: First sample name
//! - sample2: Second sample name
//! - relatedness: Coefficient of relationship (0-1)
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
pub fn compute_relatedness(combined: PathBuf, output: PathBuf, threads: usize) {
    // Set number of threads
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap_or_else(|_| {
            eprintln!("Warning: Failed to set thread count, using default");
        });

    eprintln!("Reading combined STR file...");
    let (genotype_matrix, sample_names, locus_names) = parse_combined_file(&combined);

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
            let (relatedness, n_loci, ibs0, ibs1, ibs2) =
                compute_pairwise_relatedness(&genotype_matrix, *i, *j);
            (
                sample_names[*i].clone(),
                sample_names[*j].clone(),
                relatedness,
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
    writeln!(out_file, "sample1\tsample2\trelatedness\tn_loci\tibs0\tibs1\tibs2")
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
    let header_line = crate::utils::skip_metadata_lines(&mut lines);
    let header_fields: Vec<&str> = header_line.trim().split('\t').collect();

    if header_fields.len() < 5 {
        eprintln!("ERROR: Invalid combined file format - insufficient columns");
        std::process::exit(1);
    }

    // Extract sample names from header (skip chromosome, begin, end, info columns)
    // Sample names appear as sample_H1, sample_H2 pairs
    let mut sample_names = Vec::new();
    let mut seen_samples = std::collections::HashSet::new();

    for field in header_fields.iter().skip(4) {
        if let Some(sample) = field.strip_suffix("_H1")
            && !seen_samples.contains(sample)
        {
            sample_names.push(sample.to_string());
            seen_samples.insert(sample.to_string());
        }
    }

    let n_samples = sample_names.len();
    if n_samples < 2 {
        eprintln!("ERROR: Need at least 2 samples for relatedness computation");
        std::process::exit(1);
    }

    eprintln!("Found {} samples", n_samples);

    // Parse data lines
    let mut locus_names = Vec::new();
    let mut genotype_matrix = Vec::new();

    for line_result in lines {
        let line = line_result.expect("Error reading line");
        let fields: Vec<&str> = line.trim().split('\t').collect();

        if fields.len() < 3 + n_samples * 2 {
            continue; // Skip malformed lines
        }

        // Store locus identifier
        let locus_id = format!("{}:{}-{}", fields[0], fields[1], fields[2]);
        locus_names.push(locus_id);

        // Parse genotypes for this locus (all samples)
        let mut locus_genotypes = Vec::new();
        for i in 0..n_samples {
            let h1_idx = 3 + i * 2;
            let h2_idx = 3 + i * 2 + 1;

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

/// Compute pairwise relatedness between two samples
/// Returns: (relatedness, n_informative_loci, ibs0, ibs1, ibs2)
fn compute_pairwise_relatedness(
    genotype_matrix: &[Vec<(f64, f64)>],
    sample_i: usize,
    sample_j: usize,
) -> (f64, usize, usize, usize, usize) {
    let mut ibs0 = 0; // No shared alleles
    let mut ibs1 = 0; // One shared allele
    let mut ibs2 = 0; // Two shared alleles
    let mut n_informative = 0;

    for locus_genotypes in genotype_matrix.iter() {
        let (h1_i, h2_i) = locus_genotypes[sample_i];
        let (h1_j, h2_j) = locus_genotypes[sample_j];

        // Skip loci where either sample has missing data
        if h1_i.is_nan() || h2_i.is_nan() || h1_j.is_nan() || h2_j.is_nan() {
            continue;
        }

        n_informative += 1;

        // Count shared alleles using IBS (Identity By State)
        // Compare all allele pairs and count matches
        let mut matches = 0;

        // Round to nearest integer for comparison (STR alleles are integer repeat counts)
        let h1_i_round = h1_i.round();
        let h2_i_round = h2_i.round();
        let h1_j_round = h1_j.round();
        let h2_j_round = h2_j.round();

        // Count matching alleles using IBS (Identity By State)
        // For diploid genotypes, we can have 0, 1, or 2 matching alleles
        // This follows the same approach as somalier (see: https://github.com/brentp/somalier)

        // Check if genotypes are identical (both alleles match)
        if (h1_i_round == h1_j_round && h2_i_round == h2_j_round)
            || (h1_i_round == h2_j_round && h2_i_round == h1_j_round)
        {
            matches = 2;
        } else {
            // Check for one matching allele
            if h1_i_round == h1_j_round
                || h1_i_round == h2_j_round
                || h2_i_round == h1_j_round
                || h2_i_round == h2_j_round
            {
                matches = 1;
            }
            // Otherwise matches stays at 0
        }

        match matches {
            0 => ibs0 += 1,
            1 => ibs1 += 1,
            2 => ibs2 += 1,
            _ => unreachable!(),
        }
    }

    if n_informative == 0 {
        return (f64::NAN, 0, 0, 0, 0);
    }

    // Calculate relatedness coefficient
    // IBS2 = 2 shared alleles, IBS1 = 1 shared allele, IBS0 = 0 shared alleles
    // Relatedness = (IBS2 + 0.5 * IBS1) / n_informative
    let relatedness = (ibs2 as f64 + 0.5 * ibs1 as f64) / n_informative as f64;

    (relatedness, n_informative, ibs0, ibs1, ibs2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_samples() {
        // Two identical samples should have relatedness = 1.0
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (10.0, 12.0)],
            vec![(15.0, 15.0), (15.0, 15.0)],
            vec![(8.0, 9.0), (8.0, 9.0)],
        ];

        let (relatedness, n_loci, _ibs0, _ibs1, ibs2) =
            compute_pairwise_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 3);
        assert_eq!(ibs2, 3);
        assert!((relatedness - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_unrelated_samples() {
        // Two completely different samples should have low relatedness
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (20.0, 22.0)],
            vec![(15.0, 17.0), (25.0, 27.0)],
            vec![(8.0, 9.0), (18.0, 19.0)],
        ];

        let (relatedness, n_loci, ibs0, _ibs1, _ibs2) =
            compute_pairwise_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 3);
        assert_eq!(ibs0, 3);
        assert!((relatedness - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_partial_sharing() {
        // Samples sharing one allele per locus
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (10.0, 22.0)], // Share first allele
            vec![(15.0, 17.0), (25.0, 17.0)], // Share second allele
        ];

        let (relatedness, n_loci, _ibs0, ibs1, _ibs2) =
            compute_pairwise_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 2);
        assert_eq!(ibs1, 2);
        assert_eq!(relatedness, 0.5);
    }

    #[test]
    fn test_parent_child_relationship() {
        // Parent-child: child inherits one allele from parent at each locus
        // Expected relatedness: 0.5 (all loci should be IBS1)
        let genotype_matrix = vec![
            // Parent has (10, 12), child has one from parent (10) and one new (14)
            vec![(10.0, 12.0), (10.0, 14.0)],
            // Parent has (15, 17), child has one from parent (17) and one new (19)
            vec![(15.0, 17.0), (17.0, 19.0)],
            // Parent has (8, 8), child has one from parent (8) and one new (10)
            vec![(8.0, 8.0), (8.0, 10.0)],
            // Parent has (20, 22), child has one from parent (20) and one new (24)
            vec![(20.0, 22.0), (20.0, 24.0)],
        ];

        let (relatedness, n_loci, _ibs0, ibs1, _ibs2) =
            compute_pairwise_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 4);
        assert_eq!(ibs1, 4); // All loci share exactly 1 allele
        assert!((relatedness - 0.5).abs() < 1e-6, "Parent-child relatedness should be ~0.5");
    }

    #[test]
    fn test_full_siblings() {
        // Full siblings: share both parents, so on average share 1 allele per locus
        // But can share 0, 1, or 2 alleles at different loci
        // Expected relatedness: 0.5 overall
        let genotype_matrix = vec![
            // Both inherit same alleles from parents
            vec![(10.0, 12.0), (10.0, 12.0)], // IBS2
            // Share one allele
            vec![(15.0, 17.0), (15.0, 19.0)], // IBS1
            // Share no alleles
            vec![(8.0, 9.0), (10.0, 11.0)], // IBS0
            // Share one allele
            vec![(20.0, 22.0), (20.0, 24.0)], // IBS1
        ];

        let (relatedness, n_loci, ibs0, ibs1, ibs2) =
            compute_pairwise_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 4);
        assert_eq!(ibs2, 1);
        assert_eq!(ibs1, 2);
        assert_eq!(ibs0, 1);
        // Relatedness = (1*2 + 2*1 + 1*0) / 4 = 4/4 = 1.0... wait, that's wrong
        // Relatedness = (IBS2 + 0.5*IBS1) / n = (1 + 0.5*2) / 4 = 2/4 = 0.5
        assert!((relatedness - 0.5).abs() < 1e-6, "Full sibling relatedness should be ~0.5");
    }

    #[test]
    fn test_half_siblings() {
        // Half-siblings: share one parent, so on average share 0.5 alleles per locus
        // Expected relatedness: 0.25
        let genotype_matrix = vec![
            // Share one allele (from shared parent)
            vec![(10.0, 12.0), (10.0, 14.0)], // IBS1
            // Share no alleles
            vec![(15.0, 17.0), (19.0, 21.0)], // IBS0
            // Share one allele
            vec![(8.0, 9.0), (8.0, 11.0)], // IBS1
            // Share no alleles
            vec![(20.0, 22.0), (24.0, 26.0)], // IBS0
        ];

        let (relatedness, n_loci, ibs0, ibs1, _ibs2) =
            compute_pairwise_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 4);
        assert_eq!(ibs1, 2);
        assert_eq!(ibs0, 2);
        // Relatedness = (0 + 0.5*2) / 4 = 1/4 = 0.25
        assert!((relatedness - 0.25).abs() < 1e-6, "Half-sibling relatedness should be ~0.25");
    }

    #[test]
    fn test_first_cousins() {
        // First cousins: share grandparents, expected relatedness: 0.125
        let genotype_matrix = vec![
            // Occasionally share alleles from shared grandparents
            vec![(10.0, 12.0), (10.0, 14.0)], // IBS1
            vec![(15.0, 17.0), (19.0, 21.0)], // IBS0
            vec![(8.0, 9.0), (11.0, 13.0)],   // IBS0
            vec![(20.0, 22.0), (24.0, 26.0)], // IBS0
            vec![(30.0, 32.0), (34.0, 36.0)], // IBS0
            vec![(40.0, 42.0), (44.0, 46.0)], // IBS0
            vec![(50.0, 52.0), (54.0, 56.0)], // IBS0
            vec![(60.0, 62.0), (60.0, 64.0)], // IBS1
        ];

        let (relatedness, n_loci, ibs0, ibs1, _ibs2) =
            compute_pairwise_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 8);
        assert_eq!(ibs1, 2);
        assert_eq!(ibs0, 6);
        // Relatedness = (0 + 0.5*2) / 8 = 1/8 = 0.125
        assert!((relatedness - 0.125).abs() < 1e-6, "First cousin relatedness should be ~0.125");
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

        let (relatedness, n_loci, _ibs0, ibs1, ibs2) =
            compute_pairwise_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 2); // Only 2 loci are informative
        assert_eq!(ibs2, 1);
        assert_eq!(ibs1, 1);
        // Relatedness = (1 + 0.5*1) / 2 = 1.5/2 = 0.75
        assert!((relatedness - 0.75).abs() < 1e-6);
    }

    #[test]
    fn test_homozygous_loci() {
        // Test with homozygous genotypes
        let genotype_matrix = vec![
            // Both homozygous for same allele
            vec![(10.0, 10.0), (10.0, 10.0)], // IBS2
            // One homozygous, one heterozygous with matching allele
            vec![(15.0, 15.0), (15.0, 17.0)], // IBS1
            // Both homozygous for different alleles
            vec![(20.0, 20.0), (22.0, 22.0)], // IBS0
        ];

        let (relatedness, n_loci, ibs0, ibs1, ibs2) =
            compute_pairwise_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 3);
        assert_eq!(ibs2, 1);
        assert_eq!(ibs1, 1);
        assert_eq!(ibs0, 1);
        // Relatedness = (1 + 0.5*1) / 3 = 1.5/3 = 0.5
        assert!((relatedness - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_phased_vs_unphased() {
        // Test that phase doesn't matter for IBS calculation
        // (10, 12) should match (12, 10)
        let genotype_matrix = vec![
            vec![(10.0, 12.0), (12.0, 10.0)], // Same genotype, different phase
            vec![(15.0, 17.0), (17.0, 15.0)], // Same genotype, different phase
        ];

        let (relatedness, n_loci, _ibs0, _ibs1, ibs2) =
            compute_pairwise_relatedness(&genotype_matrix, 0, 1);

        assert_eq!(n_loci, 2);
        assert_eq!(ibs2, 2); // Both should be IBS2
        assert!((relatedness - 1.0).abs() < 1e-6);
    }
}
