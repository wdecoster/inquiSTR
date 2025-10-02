# inquiSTR

[![CI](https://github.com/wdecoster/inquiSTR/actions/workflows/test.yml/badge.svg)](https://github.com/wdecoster/inquiSTR/actions/workflows/test.yml)
[![Security Audit](https://github.com/wdecoster/inquiSTR/actions/workflows/security.yml/badge.svg)](https://github.com/wdecoster/inquiSTR/actions/workflows/security.yml)
[![Crates.io](https://img.shields.io/crates/v/inquiSTR.svg)](https://crates.io/crates/inquiSTR)
[![Documentation](https://docs.rs/inquiSTR/badge.svg)](https://docs.rs/inquiSTR)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A toolkit for genotyping and analyzing Short Tandem Repeats (STRs) from long-read sequencing data. Works with Oxford Nanopore Technologies (ONT) BAM/CRAM files and supports both phased and unphased analysis.

## Features

- **STR length genotyping**: Determine repeat lengths at specific genomic loci
- **Phased analysis**: Use HP tags from phased BAM files to analyze haplotypes separately
- **Multi-sample analysis**: Combine results across samples for cohort studies
- **Quality control**: Built-in filtering and validation
- **Association testing**: R scripts for statistical analysis of STR variations
- **Cross-platform**: Pre-built binaries for Linux and macOS

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
  pca        Perform Principal Component Analysis on combined STR data
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
      --max-locus <MAX_LOCUS>      maximum locus size to consider (intervals larger than this will be filtered out)
      --batch-size <BATCH_SIZE>    batch size in KB for grouping nearby STR targets [default: 50]
  -h, --help                       Print help
```

**Examples:**

```bash
# Single region genotyping
inquiSTR call sample.bam -r chr1:1000-1100

# Multiple regions from BED file
inquiSTR call sample.bam -R regions.bed

# Filter out large intervals (>10kb) that may span problematic regions
inquiSTR call sample.bam -R regions.bed --max-locus 10000

# Multithreaded processing with custom parameters
inquiSTR call sample.bam -R regions.bed --threads 8 --minlen 10 --support 5

# CRAM file with reference
inquiSTR call sample.cram --reference genome.fa -R regions.bed

# Unphased analysis with custom sample name
inquiSTR call sample.bam -R regions.bed --unphased --sample-name "Sample123"

# Memory-constrained system (use smaller batch size)
inquiSTR call sample.bam -R regions.bed --batch-size 30 --threads 4

# High-performance system with ample RAM (use larger batch size)
inquiSTR call sample.bam -R regions.bed --batch-size 100 --threads 16
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
  -t, --threads <THREADS>  Number of threads to use for parallel processing (0 = auto-detect) [default: 0]
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

# Use 8 threads for parallel processing (for large datasets)
inquiSTR outlier combined.tsv --threads 8

# Auto-detect optimal number of threads (default behavior)
inquiSTR outlier combined.tsv --threads 0

# Combine thread control with other options
inquiSTR outlier combined.tsv --method dbscan --threads 4 --zscore 2.5
```

**Performance Optimization:**

For large datasets (>100K loci), inquiSTR outlier uses parallel processing to improve performance:

- **Automatic threading**: By default (`--threads 0`), uses all available CPU cores
- **Thread control**: Use `--threads N` to limit CPU usage when running multiple analyses
- **Memory efficiency**: Processes data in chunks to minimize memory usage
- **Progress reporting**: Shows processing progress for large datasets

**Threading Guidelines:**

- Use fewer threads (`--threads 2-4`) when running multiple analyses simultaneously
- Use more threads (`--threads 8+`) for single large-scale analyses on high-performance systems
- For very large datasets (>1M loci), consider limiting threads to avoid memory pressure

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

### `inquiSTR pca` - Principal Component Analysis

Perform Principal Component Analysis (PCA) on combined STR data to identify population structure and sample relationships. This is particularly useful for large-scale genomic analyses and quality control.

```text
Usage: inquiSTR pca [OPTIONS] <COMBINED>

Arguments:
  <COMBINED>  Combined file of STR calls from inquiSTR combine command

Options:
  -o, --output <OUTPUT>            HTML output file name for interactive PCA plot [default: pca_plot.html]
  -c, --components <COMPONENTS>    Number of principal components to compute (currently only first 2 are plotted) [default: 10]
  -t, --threads <THREADS>          Number of threads to use for parallel processing (0 = auto-detect) [default: 0]
  -a, --aggregation <AGGREGATION>  Method for aggregating H1/H2 allele lengths: max (default), min, or sum [default: max]
  -h, --help                       Print help
```

**Performance and Memory Optimization:**

For large datasets (>100K loci), inquiSTR uses intelligent feature selection and parallel processing:

- **Automatic feature selection**: Selects the most informative STR loci based on variance and completeness
- **Streaming processing**: Processes data in chunks to minimize memory usage
- **Parallel computing**: Utilizes multiple CPU cores for feature scoring and data loading
- **Thread control**: Use `--threads` to limit CPU usage (0 = use all available cores)

**Allele Aggregation Methods:**

The `--aggregation` parameter controls how the two allele lengths (H1/H2) are combined for each sample:

- **max**: Uses the longer allele length (default, good for detecting expansions)
- **min**: Uses the shorter allele length (useful for identifying deletions)
- **sum**: Uses the sum of both alleles (captures total repeat burden)

**Examples:**

```bash
# Basic PCA analysis
inquiSTR pca combined.tsv

# Use maximum aggregation with 8 threads and custom output
inquiSTR pca combined.tsv --aggregation max --threads 8 --output population_pca.html

# Use minimum aggregation to focus on deletions
inquiSTR pca combined.tsv --aggregation min --output deletion_pca.html

# Use sum aggregation to capture total repeat burden
inquiSTR pca combined.tsv --aggregation sum --threads 4 --output total_burden_pca.html

# Compute more principal components (only first 2 are plotted)
inquiSTR pca combined.tsv --components 20 --output detailed_pca.html
```

**Performance Tips:**

- For very large datasets (>1M loci), consider using fewer threads to avoid memory pressure
- The sum aggregation method may be more sensitive to technical artifacts
- Max aggregation is recommended for most population genetics applications

## Performance Tuning

The `--batch-size` parameter controls how inquiSTR groups nearby STR targets for processing, affecting both memory usage and I/O efficiency:

- **Default (50KB)**: Works well for most systems and datasets
- **Smaller values (20-35KB)**: Use on memory-constrained systems or when processing many samples simultaneously
- **Larger values (80-100KB)**: Optimal for high-performance systems with ample RAM, especially beneficial under heavy system load

**Guidelines:**

- Increase batch size if you have >32GB RAM and fast storage (SSD)
- Decrease batch size if experiencing memory pressure or running many parallel jobs
- Test different values with your specific data to find the optimal setting

## Usage for Association Testing

This repository contains `STR_regression.R`, code to perform association testing of STRs with built-in parallelization. The code is written in R and can be found in the scripts folder. This version provides improved performance and memory efficiency for large datasets. In [STR_regression_examples.md](STR_regression_examples.md) are some examples. Please open an issue if the usage of this script is unclear after reading the examples below.
