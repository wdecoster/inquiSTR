#!/usr/bin/env python3
"""
Convert TRGT VCF output to inquiSTR-style TSV format.

This script reads a VCF file (typically from TRGT) and converts it to the
inquiSTR output format with columns: chromosome, begin, end, info, sample_H1, sample_H2

The allele lengths are calculated as the difference between ALT and REF lengths,
representing the length relative to the reference genome.
"""

import sys
import argparse
from cyvcf2 import VCF


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert TRGT VCF to inquiSTR TSV format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "vcf",
        help="Input VCF file (can be gzipped)"
    )
    return parser.parse_args()


def get_allele_lengths(variant, sample_idx=0):
    """
    Extract allele lengths relative to reference for a sample.
    
    Returns tuple of (allele1_length, allele2_length) where lengths are
    relative to the reference (ALT length - REF length).
    Returns NA for missing alleles.
    """
    ref_len = len(variant.REF)
    
    # Get genotype for the sample
    gt = variant.genotypes[sample_idx]
    allele1_idx = gt[0]
    allele2_idx = gt[1]
    
    # Calculate lengths relative to reference
    if allele1_idx == -1 or allele1_idx is None:
        allele1_len = "NA"
    elif allele1_idx == 0:
        allele1_len = 0  # Reference allele
    else:
        alt_seq = variant.ALT[allele1_idx - 1]
        allele1_len = len(alt_seq) - ref_len
    
    if allele2_idx == -1 or allele2_idx is None:
        allele2_len = "NA"
    elif allele2_idx == 0:
        allele2_len = 0  # Reference allele
    else:
        alt_seq = variant.ALT[allele2_idx - 1]
        allele2_len = len(alt_seq) - ref_len
    
    return str(allele1_len), str(allele2_len)


def vcf_to_inquistr(vcf_path):
    """
    Convert VCF to inquiSTR format.
    
    Args:
        vcf_path: Path to input VCF file
    """
    vcf = VCF(vcf_path)

    # Get sample name from VCF
    if len(vcf.samples) == 0:
        print("Error: No samples found in VCF file", file=sys.stderr)
        sys.exit(1)

    sample_name = vcf.samples[0]

    # Write header
    print("# file_type=individual_call")
    print(f"# source=vcf_to_inquistr.py")
    print(f"# input_vcf={vcf_path}")
    print(f"# sample={sample_name}")
    print(f"chromosome\tbegin\tend\tinfo\t{sample_name}_H1\t{sample_name}_H2")
    # Process variants
    for variant in vcf:
        chrom = variant.CHROM
        start = variant.POS - 1  # Convert to 0-based

        # Get END from INFO field, or use POS + REF length
        end = variant.INFO['END'] if 'END' in variant.INFO else start + len(variant.REF)

        # Get variant ID as info, or use '.'
        info = variant.ID if variant.ID else '.'

        # Get allele lengths
        allele1_len, allele2_len = get_allele_lengths(variant)

        print(f"{chrom}\t{start}\t{end}\t{info}\t{allele1_len}\t{allele2_len}")

    vcf.close()


def main():
    args = parse_args()
    vcf_to_inquistr(args.vcf)


if __name__ == "__main__":
    main()
