# inquiSTR

[![CI](https://github.com/wdecoster/inquiSTR/actions/workflows/test.yml/badge.svg)](https://github.com/wdecoster/inquiSTR/actions/workflows/test.yml)
[![Security Audit](https://github.com/wdecoster/inquiSTR/actions/workflows/security.yml/badge.svg)](https://github.com/wdecoster/inquiSTR/actions/workflows/security.yml)
[![Crates.io](https://img.shields.io/crates/v/inquiSTR.svg)](https://crates.io/crates/inquiSTR)
[![Documentation](https://docs.rs/inquiSTR/badge.svg)](https://docs.rs/inquiSTR)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A high-performance toolset to genotype and analyze Short Tandem Repeats (STRs) from long-read sequencing data. Optimized for Oxford Nanopore Technologies (ONT) data with support for both phased and unphased analysis.

## 🚀 Features

- **Fast STR length determination**: Memory-efficient batch processing for millions of targets
- **Phased Analysis**: Support for HP-tagged BAM files with comprehensive validation
- **Downstream Analysis**: Association testing and outlier detection across cohorts  
- **Cross-platform**: Pre-built binaries for Linux (glibc/musl) and macOS

## 📊 Performance

inquiSTR is optimized for large-scale STR analysis:
- **Region batching**: 10-100x reduction in I/O operations
- **Memory efficiency**: ~99% reduction in memory usage vs. naive approaches
- **Optimized data structures**: Custom STR call storage and processing

## 📦 Installation

### Pre-built Binaries (Recommended)

Download the latest binary for your system from the [releases page](https://github.com/wdecoster/inquiSTR/releases). Make sure the binary is in your $PATH or use the full path to execute.

```bash
# Linux (glibc)
curl -L https://github.com/wdecoster/inquiSTR/releases/latest/download/inquiSTR-linux -o inquiSTR
chmod +x inquiSTR

# Linux (musl - static binary)  
curl -L https://github.com/wdecoster/inquiSTR/releases/latest/download/inquiSTR-linux-musl -o inquiSTR
chmod +x inquiSTR

# macOS
curl -L https://github.com/wdecoster/inquiSTR/releases/latest/download/inquiSTR-macos -o inquiSTR
chmod +x inquiSTR
```

### From Source

```bash
git clone https://github.com/wdecoster/inquiSTR.git
cd inquiSTR  
cargo build --release
# Binary will be in target/release/inquiSTR
```

### Using Cargo

```bash
cargo install inquiSTR
```

This repository contains Rust code for inquiSTR, a toolset to genotype and analyze STRs from long read sequencing data, and has been tested with ONT data.

## Usage

The inquiSTR tool has several subcommands, as detailed below. All commands write to stdout.

```text
Usage: inquiSTR <COMMAND>

Commands:
  call       Call lengths
  combine    Combine lengths from multiple bams to a TSV
  outlier    Find outliers from TSV
  query      Lookup genotypes and display
  histogram
  plot       Show a histogram with multiple groups for a specific repeat
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

### `inquiSTR call` - STR Genotyping

Variants are genotyped with `inquiSTR call`, which will determine the length of each haplotype of each STR locus.

```text
Usage: inquiSTR call [OPTIONS] <BAM>

Arguments:
  <BAM>  bam file to call STRs in

Options:
  -r, --region <REGION>            region string to genotype expansion in
  -R, --region-file <REGION_FILE>  Bed file with region(s) to genotype expansion(s) in
  -m, --minlen <MINLEN>            minimal length of insertion/deletion operation [default: 5]
  -s, --support <SUPPORT>          minimal number of supporting reads [default: 3]
  -t, --threads <THREADS>          Number of parallel threads to use [default: 1]
  -u, --unphased                   If reads have to be considered unphased
      --sample-name <SAMPLE_NAME>  sample name to use in output
      --reference <REFERENCE>      reference fasta for cram decoding
  -h, --help                       Print help
```

**Examples:**

```bash
# Single region genotyping
inquiSTR call sample.bam -r chr1:1000-1100

# Multiple regions from BED file
inquiSTR call sample.bam -R regions.bed

# Multithreaded processing with custom parameters
inquiSTR call sample.bam -R regions.bed --threads 8 --minlen 10 --support 5

# CRAM file with reference
inquiSTR call sample.cram --reference genome.fa -R regions.bed

# Unphased analysis with custom sample name
inquiSTR call sample.bam -R regions.bed --unphased --sample-name "Sample123"
```

### `inquiSTR combine` - Multi-sample Analysis

Variants from multiple samples can be combined with `inquiSTR combine`.

```text
Usage: inquiSTR combine <CALLS>...

Arguments:
  <CALLS>...  files from inquiSTR call

Options:
  -h, --help  Print help
```

**Examples:**

```bash
# Combine STR calls from multiple samples
inquiSTR combine sample1.inq sample2.inq sample3.inq > combined.tsv

# Combine all .inq files in current directory
inquiSTR combine *.inq > cohort_combined.tsv
```

### `inquiSTR query` - Genotype Lookup

Querying genotypes from a combined file can be done with `inquiSTR query`, taking a region string or a file with regions to query. When querying a single locus, output will be sorted by repeat length (descending, with highest expansions first).

```text
Usage: inquiSTR query <COMBINED> <REGION>

Arguments:
  <COMBINED>  combined file of calls
  <REGION>    region to query or file with regions to query

Options:
  -h, --help  Print help
```

**Examples:**

```bash
# Query a specific region (single locus - output sorted by length)
inquiSTR query combined.tsv chr1:1000-2000

# Query multiple regions from a BED file (multi-locus table format)
inquiSTR query combined.tsv regions.bed

# Query using a text file with multiple coordinate strings
inquiSTR query combined.tsv regions.txt
```

### `inquiSTR outlier` - Outlier Detection

Identifying outliers from a combined file can be done with `inquiSTR outlier`, using either z-scores or DBSCAN.

```text
Usage: inquiSTR outlier [OPTIONS] <COMBINED>

Arguments:
  <COMBINED>  combined file of calls

Options:
      --minsize <MINSIZE>  minimal length of expansion to be present in cohort [default: 10]
  -z, --zscore <ZSCORE>    zscore cutoff to decide if a value is an outlier [default: 3]
      --method <METHOD>    method to test for outliers [default: zscore] [possible values: zscore, dbscan]
  -s, --sample <SAMPLE>    sample to consider
  -S, --subset <SUBSET>    file with subset of samples to consider
  -h, --help               Print help
```

**Examples:**

```bash
# Find outliers using default z-score method
inquiSTR outlier combined.tsv

# Find outliers with stricter z-score threshold
inquiSTR outlier combined.tsv --zscore 2.5

# Find outliers using DBSCAN method
inquiSTR outlier combined.tsv --method dbscan

# Find outliers for specific sample only
inquiSTR outlier combined.tsv --sample "Sample123"

# Find outliers for subset of samples
inquiSTR outlier combined.tsv --subset sample_list.txt

# Custom minimum expansion size
inquiSTR outlier combined.tsv --minsize 15
```

### `inquiSTR histogram` - Data Visualization

Generate histograms for STR length distributions at specific loci.

```text
Usage: inquiSTR histogram <COMBINED> <REGION>

Arguments:
  <COMBINED>  combined file of calls
  <REGION>    region to query

Options:
  -h, --help  Print help
```

**Examples:**

```bash
# Generate histogram for specific region
inquiSTR histogram combined.tsv chr1:1000-2000

# Generate histogram for trinucleotide repeat
inquiSTR histogram combined.tsv chr4:3074877-3074933
```

### `inquiSTR plot` - Group Comparison Plots  

Show a histogram with multiple groups for a specific repeat, useful for comparing cohorts.

```text
Usage: inquiSTR plot [OPTIONS] <COMBINED> <METADATA> <REGION>

Arguments:
  <COMBINED>  combined file of calls
  <METADATA>  file with sample_id, phenotype and covariates
  <REGION>    region to query

Options:
  -c, --condition <CONDITION>  test column and groups to plot e.g. group:PAT,CON
  -o, --output <OUTPUT>        HTML output file name [default: groupplot.html]
  -h, --help                   Print help
```

**Examples:**

```bash
# Compare patient vs control groups
inquiSTR plot combined.tsv metadata.tsv chr1:1000-2000 --condition "group:PAT,CON"

# Compare multiple groups with custom output
inquiSTR plot combined.tsv metadata.tsv chr4:3074877-3074933 \
  --condition "diagnosis:AD,PD,Control" --output "disease_comparison.html"

# Simple binary comparison  
inquiSTR plot combined.tsv metadata.tsv chr15:34419414-34419461 \
  --condition "status:case,control" --output "case_control_plot.html"
```

## Usage for Association Testing

This repository contains `STR_regression.R`, code to perform association testing of STRs with built-in parallelization. The code is written in R and can be found in the scripts folder. This version provides improved performance and memory efficiency for large datasets.

Below are some worked usage examples for "MAX" STRmode for binary phenotypes, without covariates. Please open an issue if the usage of this script is unclear after reading the examples below.

### Full (genome-wide) run  

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run full --out full_genome_wide_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```

### Chromosome-wide run  

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run chromosome --chr chr15 --out chr15_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```

### Chromosome interval run

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run chr_interval --chr chr15 --chr_begin 34419410 --chr_end 34419465 --out chr15_34419410_34419465_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```

### Bed interval run

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run bed_interval --bed chr15_roi.bed --out bed_chr15_roi_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```

### Single variant (Expanded Allele) run

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run single_variant --single_variant chr15_34419414_34419461 --expandedAllele 201 --out singleVariant_chr15_34419414_34419461_expandedAllele201_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```

or

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run single_variant --single_variant chr15:34419414-34419461 --expandedAllele 201 --out singleVariant_chr15_34419414_34419461_expandedAllele201_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```
