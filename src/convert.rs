//! Convert VCF files to inquiSTR format
//!
//! This module converts STR-genotyper VCFs (TRGT or LongTR) to inquiSTR's TSV format.
//! It is **format-aware**: rather than reconstructing allele lengths from REF/ALT string
//! arithmetic (which is fragile across anchoring/coordinate conventions), it uses each
//! tool's native length field and repeat coordinates:
//!
//! - **TRGT**: locus coordinates are `POS`..`END` (TRGT writes the 0-based catalog start
//!   directly into `POS`, matching inquiSTR `call`), and per-haplotype lengths come from the
//!   `AL` FORMAT field (absolute allele length), reported relative to the interval width
//!   `END - POS` to match inquiSTR's reference-relative output.
//! - **LongTR**: locus coordinates are `INFO/START`..`INFO/END` (the repetitive portion of
//!   the reference allele). LongTR's `POS` is **not** the repeat start — it is shifted by
//!   flanking context — so `START` must be used for the match coordinate. Per-haplotype
//!   lengths come straight from the `GB` FORMAT field (bp difference from reference), which
//!   already matches inquiSTR's convention.
//!
//! It handles single-sample VCFs (→ individual call file), multiple samples / multiple VCFs
//! (→ combined file), and missing genotypes (→ NA).

use crate::errors::InquiSTRError;
use crate::filetype::FileType;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// Which genotyper produced the VCF. Determines coordinate and length handling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VcfSource {
    Trgt,
    LongTr,
}

impl VcfSource {
    fn label(self) -> &'static str {
        match self {
            VcfSource::Trgt => "TRGT",
            VcfSource::LongTr => "LongTR",
        }
    }
}

/// A single STR locus with per-sample genotype lengths already expressed as
/// signed bp differences from the reference (or "NA").
#[derive(Debug, Clone)]
struct VcfVariant {
    chromosome: String,
    begin: u32,
    end: u32,
    id: String,
    samples: Vec<SampleCall>,
}

/// A genotype for a single sample at a single locus, as two inquiSTR length strings.
#[derive(Debug, Clone)]
struct SampleCall {
    sample_name: String,
    h1: String,
    h2: String,
}

/// Look up a `KEY=value` entry in a VCF INFO field, returning `value`.
fn info_value<'a>(info: &'a str, key_eq: &str) -> Option<&'a str> {
    info.split(';').find_map(|kv| kv.strip_prefix(key_eq))
}

/// Detect the genotyper from the VCF meta-header (`##` lines).
fn detect_source(header_lines: &[String]) -> Option<VcfSource> {
    let has = |needle: &str| header_lines.iter().any(|l| l.contains(needle));
    // TRGT: AL FORMAT field plus MOTIFS/STRUC INFO fields.
    if has("ID=AL,") && (has("ID=MOTIFS,") || has("ID=STRUC,")) {
        return Some(VcfSource::Trgt);
    }
    // LongTR (HipSTR-derived): GB FORMAT field plus START INFO field.
    if has("ID=GB,") && has("ID=START,") {
        return Some(VcfSource::LongTr);
    }
    None
}

/// Compute the two inquiSTR haplotype length strings for one sample field.
///
/// `width` is the locus interval width (`end - begin`), used only for TRGT to turn the
/// absolute `AL` length into a reference-relative one. Returns ("NA", "NA") for missing data.
fn haplotype_lengths(
    source: VcfSource,
    sample_fields: &[&str],
    gt_idx: usize,
    len_idx: Option<usize>,
    width: u32,
) -> (String, String) {
    let gt = sample_fields.get(gt_idx).copied().unwrap_or(".");
    let sep = if gt.contains('|') { '|' } else { '/' };
    let gt_alleles: Vec<&str> = gt.split(sep).collect();

    // Length values are in genotype (haplotype) order. TRGT separates AL with ','; LongTR
    // separates GB with the phase character. Split on any of them to handle both.
    let len_field = len_idx
        .and_then(|i| sample_fields.get(i))
        .copied()
        .unwrap_or(".");
    let len_vals: Vec<&str> = len_field.split([',', '|', '/']).collect();

    let value_at = |k: usize| -> String {
        // A missing GT allele means a missing haplotype.
        if gt_alleles.get(k).map(|a| *a == ".").unwrap_or(true) {
            return "NA".to_string();
        }
        let raw = match len_vals.get(k) {
            Some(v) if *v != "." && !v.is_empty() => *v,
            _ => return "NA".to_string(),
        };
        let parsed: i64 = match raw.parse() {
            Ok(v) => v,
            Err(_) => return "NA".to_string(),
        };
        match source {
            // AL is the absolute allele length; report relative to the interval width.
            VcfSource::Trgt => (parsed - width as i64).to_string(),
            // GB is already the bp difference from reference.
            VcfSource::LongTr => parsed.to_string(),
        }
    };

    (value_at(0), value_at(1))
}

/// Parse a VCF file and extract STR variants with per-sample inquiSTR lengths.
/// Returns the sample names, the variants, and the detected genotyper.
fn parse_vcf(vcf_path: &Path) -> Result<(Vec<String>, Vec<VcfVariant>, VcfSource), InquiSTRError> {
    let reader = crate::utils::reader(&vcf_path.to_string_lossy());
    let mut sample_names: Vec<String> = Vec::new();
    let mut variants: Vec<VcfVariant> = Vec::new();
    let mut header_lines: Vec<String> = Vec::new();
    let mut source: Option<VcfSource> = None;

    for line in reader.lines() {
        let line = line.map_err(|e| InquiSTRError {
            message: format!("Failed to read VCF file {}: {}", vcf_path.display(), e),
        })?;

        // Collect meta-header for format detection.
        if line.starts_with("##") {
            header_lines.push(line);
            continue;
        }

        // Column header: detect the genotyper and read sample names.
        if line.starts_with("#CHROM") {
            source = Some(detect_source(&header_lines).ok_or_else(|| InquiSTRError {
                message: format!(
                    "Unrecognized VCF format in {}: expected TRGT (AL + MOTIFS/STRUC) or \
                     LongTR (GB + START) header fields.",
                    vcf_path.display()
                ),
            })?);
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() > 9 {
                sample_names = fields[9..].iter().map(|s| s.to_string()).collect();
            }
            continue;
        }

        if line.starts_with('#') {
            continue;
        }

        let source = source.ok_or_else(|| InquiSTRError {
            message: format!("VCF {} has data lines before the #CHROM header", vcf_path.display()),
        })?;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue; // Skip lines without genotype data
        }

        let chromosome = fields[0].to_string();
        let pos: u32 = fields[1].parse().map_err(|_| InquiSTRError {
            message: format!("Invalid position in VCF: {}", fields[1]),
        })?;
        let info = fields[7];
        let ref_len = fields[3].len() as u32;

        // Locus coordinates: use the field that matches the catalog / inquiSTR `call`.
        let (begin, end) = match source {
            VcfSource::Trgt => {
                // TRGT writes the 0-based catalog start directly into POS.
                let end = info_value(info, "END=")
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(pos + ref_len);
                (pos, end)
            }
            VcfSource::LongTr => {
                // LongTR POS is shifted by flanking context; START is the repeat start.
                let begin = info_value(info, "START=")
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(pos);
                let end = info_value(info, "END=")
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(begin + ref_len);
                (begin, end)
            }
        };

        let id = if fields[2] == "." {
            format!("{}:{}", chromosome, begin)
        } else {
            fields[2].to_string()
        };

        // Locate GT and the source-specific length field in FORMAT.
        let format_fields: Vec<&str> = fields[8].split(':').collect();
        let gt_idx = match format_fields.iter().position(|&s| s == "GT") {
            Some(i) => i,
            None => continue, // No genotype information
        };
        let len_tag = match source {
            VcfSource::Trgt => "AL",
            VcfSource::LongTr => "GB",
        };
        let len_idx = format_fields.iter().position(|&s| s == len_tag);

        let width = end.saturating_sub(begin);
        let mut samples = Vec::with_capacity(fields.len().saturating_sub(9));
        for (i, sample_field) in fields[9..].iter().enumerate() {
            let sample_name = if i < sample_names.len() {
                sample_names[i].clone()
            } else {
                format!("Sample{}", i + 1)
            };
            let sf: Vec<&str> = sample_field.split(':').collect();
            let (h1, h2) = haplotype_lengths(source, &sf, gt_idx, len_idx, width);
            samples.push(SampleCall { sample_name, h1, h2 });
        }

        variants.push(VcfVariant { chromosome, begin, end, id, samples });
    }

    if sample_names.is_empty() {
        return Err(InquiSTRError {
            message: format!("No samples found in VCF file: {}", vcf_path.display()),
        });
    }

    let source = source.ok_or_else(|| InquiSTRError {
        message: format!("VCF {} has no #CHROM header line", vcf_path.display()),
    })?;

    Ok((sample_names, variants, source))
}

/// Convert VCF files to inquiSTR format
pub fn convert_vcf(vcf_files: Vec<PathBuf>) -> Result<(), InquiSTRError> {
    if vcf_files.is_empty() {
        return Err(InquiSTRError { message: "No VCF files provided".to_string() });
    }

    // Parse all VCF files
    let mut all_samples: Vec<String> = Vec::new();
    let mut all_variants_by_locus: HashMap<String, Vec<VcfVariant>> = HashMap::new();
    let mut sample_counts: HashMap<String, usize> = HashMap::new();
    let mut sources: Vec<VcfSource> = Vec::new();

    for vcf_path in &vcf_files {
        eprintln!("Processing VCF file: {}", vcf_path.display());
        let (sample_names, mut variants, source) = parse_vcf(vcf_path)?;
        eprintln!("  detected genotyper: {}", source.label());
        sources.push(source);

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
                    if genotype.sample_name == *sample {
                        genotype.sample_name = final_sample_name.clone();
                    }
                }
            }

            all_samples.push(final_sample_name);
        }

        // Group variants by locus
        for variant in variants {
            let locus_key = format!("{}:{}:{}", variant.chromosome, variant.begin, variant.end);
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
        let from = sources.get(i).map(|s| s.label()).unwrap_or("unknown");
        println!("# input_vcf_{}={} (converted_from={})", i + 1, vcf_path.display(), from);
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
            let parts: Vec<&str> = s.rsplitn(3, ':').collect();
            // rsplitn yields parts reversed: [end, begin, chrom]
            (
                parts.get(2).copied().unwrap_or("").to_string(),
                parts.get(1).and_then(|v| v.parse().ok()).unwrap_or(0),
                parts.first().and_then(|v| v.parse().ok()).unwrap_or(0),
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
            first_variant.chromosome, first_variant.begin, first_variant.end, first_variant.id
        );

        // Merge genotypes from all variants at this locus
        let mut genotypes_by_sample: HashMap<String, (String, String)> = HashMap::new();
        for variant in variants {
            for genotype in &variant.samples {
                genotypes_by_sample.insert(
                    genotype.sample_name.clone(),
                    (genotype.h1.clone(), genotype.h2.clone()),
                );
            }
        }

        // Write genotypes in sample order; loci absent from a sample's VCF are NA.
        for sample in &all_samples {
            if let Some((h1, h2)) = genotypes_by_sample.get(sample) {
                print!("\t{}\t{}", h1, h2);
            } else {
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
    fn test_info_value() {
        let info = "START=16682;END=16774;MOTIF=TGGTGGGGG;PERIOD=9";
        assert_eq!(info_value(info, "START="), Some("16682"));
        assert_eq!(info_value(info, "END="), Some("16774"));
        assert_eq!(info_value(info, "AC="), None);
    }

    #[test]
    fn test_detect_source() {
        let trgt = vec![
            "##INFO=<ID=MOTIFS,Number=.,Type=String,Description=\"..\">".to_string(),
            "##FORMAT=<ID=AL,Number=.,Type=Integer,Description=\"Length of each allele\">"
                .to_string(),
        ];
        assert_eq!(detect_source(&trgt), Some(VcfSource::Trgt));

        let longtr = vec![
            "##INFO=<ID=START,Number=1,Type=Integer,Description=\"..\">".to_string(),
            "##FORMAT=<ID=GB,Number=1,Type=String,Description=\"..\">".to_string(),
        ];
        assert_eq!(detect_source(&longtr), Some(VcfSource::LongTr));

        let plain = vec!["##FORMAT=<ID=GT,Number=1,Type=String>".to_string()];
        assert_eq!(detect_source(&plain), None);
    }

    #[test]
    fn test_trgt_lengths_relative_to_width() {
        // chr1:20798-20893 -> width 95. AL=95,95 (ref) => 0,0. AL=143,165 => -... relative.
        // FORMAT GT:AL -> indices gt=0, al=1
        let sf = ["0/0", "95,95"];
        assert_eq!(
            haplotype_lengths(VcfSource::Trgt, &sf, 0, Some(1), 95),
            ("0".to_string(), "0".to_string())
        );
        let sf = ["1/2", "143,165"];
        assert_eq!(
            haplotype_lengths(VcfSource::Trgt, &sf, 0, Some(1), 165),
            ("-22".to_string(), "0".to_string())
        );
        // Fully missing
        let sf = ["."];
        assert_eq!(
            haplotype_lengths(VcfSource::Trgt, &sf, 0, Some(1), 95),
            ("NA".to_string(), "NA".to_string())
        );
    }

    #[test]
    fn test_longtr_lengths_from_gb() {
        // GB is the bp diff directly; width is irrelevant for LongTR.
        // FORMAT GT:GB -> indices gt=0, gb=1
        let sf = ["0|1", "0|-86"];
        assert_eq!(
            haplotype_lengths(VcfSource::LongTr, &sf, 0, Some(1), 86),
            ("0".to_string(), "-86".to_string())
        );
        let sf = ["1|2", "-77|-10"];
        assert_eq!(
            haplotype_lengths(VcfSource::LongTr, &sf, 0, Some(1), 77),
            ("-77".to_string(), "-10".to_string())
        );
        // Half missing
        let sf = [".|1", ".|5"];
        assert_eq!(
            haplotype_lengths(VcfSource::LongTr, &sf, 0, Some(1), 50),
            ("NA".to_string(), "5".to_string())
        );
    }
}
