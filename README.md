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
  combine    Combine STR calls or kmer frequencies from multiple samples to a TSV
  outlier    Find outliers from combined STR or kmer data
  query      Lookup genotypes and display
  histogram
  plot       Show a histogram with multiple groups for a specific repeat
  pca        Perform Principal Component Analysis on combined STR data
  unmapped   Count kmer frequencies in unmapped reads
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

# Use multiple threads for faster processing
inquiSTR call sample.bam -R regions.bed --threads 8 --minlen 10 --support 5

# CRAM file with reference
inquiSTR call sample.cram --reference genome.fa -R regions.bed

# Unphased analysis with custom sample name
inquiSTR call sample.bam -R regions.bed --unphased --sample-name "Sample123"

# Custom batch size (in kb) and threads
inquiSTR call sample.bam -R regions.bed --batch-size 30 --threads 4
```

### `inquiSTR combine` - Multi-sample Analysis

Combine data from multiple samples with `inquiSTR combine`. This command supports both STR call files (from `inquiSTR call`) and kmer frequency files (from `inquiSTR unmapped`), automatically detecting the input format.

```text
Usage: inquiSTR combine [OPTIONS] <CALLS>...

Arguments:
  <CALLS>...  files from inquiSTR call or inquiSTR unmapped

Options:
  -t, --threads <THREADS>  Number of threads to use for parallel processing [default: 1]
  -h, --help               Print help
```

**Examples:**

```bash
# Combine STR calls from multiple samples
inquiSTR combine sample1.inq sample2.inq sample3.inq > str_combined.tsv

# Combine kmer frequency files from unmapped analysis
inquiSTR combine sample1_kmers.tsv sample2_kmers.tsv sample3_kmers.tsv > kmer_combined.tsv

# Combine all .inq files in current directory (STR data)
inquiSTR combine *.inq > cohort_combined.tsv

# Use multiple threads 
inquiSTR combine *.inq --threads 8 > combined.tsv
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

Identify outliers from combined data using `inquiSTR outlier`. This command works with both STR call data and kmer frequency data (automatically detecting the format), using either z-scores or DBSCAN algorithms.

```text
Usage: inquiSTR outlier [OPTIONS] <COMBINED>

Arguments:
  <COMBINED>  combined file from inquiSTR combine (STR calls or kmer frequencies)

Options:
      --minsize <MINSIZE>  minimal length of expansion to be present in cohort [default: 10]
  -z, --zscore <ZSCORE>    zscore cutoff to decide if a value is an outlier [default: 3]
      --method <METHOD>    method to test for outliers [default: zscore] [possible values: zscore, dbscan]
  -s, --sample <SAMPLE>    sample to consider
  -S, --subset <SUBSET>    file with subset of samples to consider
  -t, --threads <THREADS>  Number of threads to use for parallel processing [default: 1]
  -h, --help               Print help
```

**Examples:**

```bash
# Find outliers in STR call data using default z-score method
inquiSTR outlier str_combined.tsv

# Find outliers in kmer frequency data
inquiSTR outlier kmer_combined.tsv --zscore 2.5

# Find outliers using DBSCAN method (works with both data types)
inquiSTR outlier combined.tsv --method dbscan

# Find outliers for specific sample only
inquiSTR outlier combined.tsv --sample "Sample123"

# Find outliers for subset of samples
inquiSTR outlier combined.tsv --subset sample_list.txt

# Custom minimum expansion size
inquiSTR outlier combined.tsv --minsize 15

# Use multiple threads
inquiSTR outlier combined.tsv --threads 8
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

### `inquiSTR unmapped` - Kmer Frequency Analysis

Count kmer frequencies in unmapped reads from BAM/CRAM files. This analysis extracts all unmapped reads and counts occurrences of kmers of all sizes from 2 to a specified maximum length. All kmer rotations (e.g., CAG, AGC, GCA) are represented by their lexicographically smallest form (AGC) for consistent analysis.

```text
Usage: inquiSTR unmapped [OPTIONS] <BAM>

Arguments:
  <BAM>  BAM/CRAM file to analyze unmapped reads from

Options:
  -k, --klength <KLENGTH>          Maximum kmer length to count (counts all sizes from 2 to klength) [default: 6]
      --sample-name <SAMPLE_NAME>  Sample name to use in output header
      --reference <REFERENCE>      Reference fasta for CRAM decoding
  -t, --threads <THREADS>          Number of parallel threads to use [default: 1]
  -h, --help                       Print help
```

**Key Features:**

- **Comprehensive kmer analysis**: Counts all kmer sizes from 2 to `--klength` (default 6)
- **Canonical representation**: All rotations of the same kmer are represented by the lexicographically smallest form
- **Complete output**: Includes all possible kmers (including zeros) sorted alphabetically by k-size
- **Normalized counts**: Results are normalized by total read count in the entire BAM/CRAM file
- **Parallel processing**: Batch-based parallel processing for efficient analysis of large files
- **Quality filtering**: Automatically excludes reads containing N bases and empty sequences

**Output Format:**

The output is a TSV file with two columns:

- **kmer**: All possible kmers of each length (2-mers, then 3-mers, etc.)
- **sample_name**: Normalized frequency counts (kmer_count / total_reads)

**Examples:**

```bash
# Basic kmer analysis with default settings (k=6, 1 thread)
inquiSTR unmapped sample.bam > kmers.tsv

# Custom kmer length and sample name
inquiSTR unmapped sample.bam --klength 4 --sample-name MySample > kmers_k4.tsv

# Use multiple threads
inquiSTR unmapped sample.bam --klength 5 --threads 8 --sample-name MySample

# CRAM file with reference
inquiSTR unmapped sample.cram --reference genome.fasta --klength 3 --sample-name CramSample

# Remote file analysis
inquiSTR unmapped https://example.com/sample.bam --threads 4 --sample-name RemoteSample

# Test with provided example data
inquiSTR unmapped test-data/unmapped.bam --klength 3 --sample-name test_sample
```

**Example Output:**

```tsv
kmer    demo_sample
AA      1080.441548
AC      1586.379696
AG      1649.190091
AT      1545.704501
CC      734.356865
CG      765.888925
CT      1648.491418
GG      737.300032
GT      1598.415641
TT      1079.760363
AAA     ...
```

Note that only canonical forms are shown. Rotations like CA, GA, GC, TA, TC, TG are represented by their canonical forms (AC, AG, CG, AT, CT, GT respectively) and their counts are combined.

**Use Cases:**

- **Contamination detection**: Identify unexpected sequences in unmapped reads
- **Species identification**: Compare kmer profiles against known organisms
- **Quality control**: Assess the composition of unmapped reads
- **Metagenomic analysis**: Analyze microbial content in unmapped fractions
- **Technical artifact detection**: Identify adapter sequences or other technical contamination

**Integration with Other Commands:**

The kmer frequency output from `inquiSTR unmapped` can be processed with other inquiSTR commands:

```bash
# Complete kmer analysis workflow
inquiSTR unmapped sample1.bam > sample1_kmers.tsv
inquiSTR unmapped sample2.bam > sample2_kmers.tsv
inquiSTR unmapped sample3.bam > sample3_kmers.tsv

# Combine kmer data from multiple samples
inquiSTR combine sample*_kmers.tsv > combined_kmers.tsv

# Find outlier samples based on kmer frequencies
inquiSTR outlier combined_kmers.tsv --method dbscan
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
  -t, --threads <THREADS>          Number of threads to use for parallel processing [default: 1]
  -a, --aggregation <AGGREGATION>  Method for aggregating H1/H2 allele lengths: max (default), min, or sum [default: max]
  -h, --help                       Print help
```

**Allele Aggregation Methods:**

The `--aggregation` parameter controls how the two allele lengths (H1/H2) are combined for each sample:

- **max**: Uses the longer allele length (default, good for detecting expansions)
- **min**: Uses the shorter allele length (useful for identifying deletions)
- **sum**: Uses the sum of both alleles (captures total repeat burden)

**Examples:**

```bash
# Basic PCA analysis
inquiSTR pca combined.tsv

# Use maximum aggregation (default) with 8 threads and custom output
inquiSTR pca combined.tsv --aggregation max --threads 8 --output population_pca.html

# Use minimum aggregation to focus on deletions
inquiSTR pca combined.tsv --aggregation min --output deletion_pca.html

# Use sum aggregation to capture total repeat burden
inquiSTR pca combined.tsv --aggregation sum --threads 4 --output total_burden_pca.html

# Compute more principal components (only first 2 are plotted)
inquiSTR pca combined.tsv --components 20 --output detailed_pca.html
```

## Usage for Association Testing

This repository contains `STR_regression.R`, an R script for association testing of STRs. The script can be found in the scripts folder. Examples are provided in [STR_regression_examples.md](STR_regression_examples.md). Please open an issue if the usage is unclear.
