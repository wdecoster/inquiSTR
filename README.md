# inquiSTR

[![CI](https://github.com/wdecoster/inquiSTR/actions/workflows/test.yml/badge.svg)](https://github.com/wdecoster/inquiSTR/actions/workflows/test.yml)
[![Security Audit](https://github.com/wdecoster/inquiSTR/actions/workflows/security.yml/badge.svg)](https://github.com/wdecoster/inquiSTR/actions/workflows/security.yml)
[![Crates.io](https://img.shields.io/crates/v/inquiSTR.svg)](https://crates.io/crates/inquiSTR)
[![Documentation](https://docs.rs/inquiSTR/badge.svg)](https://docs.rs/inquiSTR)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A toolkit for lightning-fast Short Tandem Repeat (STR) length genotyping and downstream analysis from long-read sequencing data. inquiSTR works with Oxford Nanopore Technologies and PacBio BAM/CRAM files and supports both phased and unphased data.

## Table of Contents

- [Features](#features)
  - [Remote File Access](#remote-file-access)
- [Installation](#-installation)
  - [Pre-built Binaries](#pre-built-binaries-recommended)
  - [From Source](#from-source)
  - [Using Cargo](#using-cargo)
- [Usage](#usage)
  - [inquiSTR call - STR Genotyping](#inquistr-call---str-genotyping)
  - [inquiSTR combine - Multi-sample Analysis](#inquistr-combine---multi-sample-analysis)
  - [inquiSTR query - Genotype Lookup](#inquistr-query---genotype-lookup)
  - [inquiSTR outlier - Outlier Detection](#inquistr-outlier---outlier-detection)
  - [inquiSTR histogram - Data Visualization](#inquistr-histogram---data-visualization)
  - [inquiSTR plot - Group Comparison Plots](#inquistr-plot---group-comparison-plots)
  - [inquiSTR unmapped - Kmer Frequency Analysis](#inquistr-unmapped---kmer-frequency-analysis)
  - [inquiSTR pca - Principal Component Analysis](#inquistr-pca---principal-component-analysis)
  - [inquiSTR association - Statistical Association Testing](#inquistr-association---statistical-association-testing)
- [Legacy R Script Usage](#legacy-r-script-usage)
- [Development](#️-development)
  - [Development Setup](#development-setup)
  - [Code Quality](#code-quality)
  - [Contributing](#contributing)

## Features

- **STR length genotyping**: Determine repeat lengths at specific genomic loci
- **Phased analysis**: Use HP tags from phased BAM/CRAM files to analyze haplotypes separately
- **Multi-sample analysis**: Combine results across samples for cohort studies
- **Association testing**: Perform statistical analysis of STR variations
- **Remote file support**: Seamless access to HTTP/HTTPS/FTP/S3 URLs with automatic index caching
- **Parallelized code**: Efficient use of multi-core systems

## 📦 Installation

### Pre-built Binaries (Recommended)

Pre-built binaries for Linux and macOS can be downloaded from the [releases page](https://github.com/wdecoster/inquiSTR/releases). Make sure the binary is in your $PATH or use the full path to execute.

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

#### System dependencies (Linux)

Building from source requires a C/C++ toolchain and common compression/network libraries used by HTSlib. On Ubuntu/Debian, install:

```bash
sudo apt-get update
sudo apt-get install -y build-essential pkg-config clang libclang-dev \
  zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
  libdeflate-dev libzstd-dev libssl-dev
```

These provide standard headers for clang/bindgen (fixes errors like "stddef.h not found") and libraries used by BAM/CRAM support.

#### Static MUSL build (fully portable Linux binary)

To build a binary that does not depend on your system glibc (useful for older servers or containers):

```bash
# 1) Install musl toolchain support
rustup target add x86_64-unknown-linux-musl

# Option A (recommended): use cross for a reproducible build
cargo install cross   # once
OPENSSL_STATIC=1 LIBZ_SYS_STATIC=1 BZIP2_STATIC=1 ZSTD_STATIC=1 LZMA_API_STATIC=1 CURL_STATIC=1 \
  cross build --release --target x86_64-unknown-linux-musl

# Option B: use cargo directly (requires musl-tools on host)
sudo apt-get install -y musl-tools
OPENSSL_STATIC=1 LIBZ_SYS_STATIC=1 BZIP2_STATIC=1 ZSTD_STATIC=1 LZMA_API_STATIC=1 CURL_STATIC=1 \
  cargo build --release --target x86_64-unknown-linux-musl

# Result
ls -lh target/x86_64-unknown-linux-musl/release/inquiSTR
```

Alternatively, use the helper:

```bash
scripts/build_musl.sh
# or via Makefile
make build-musl
```

Notes:

- Building with MUSL avoids glibc version errors (e.g., GLIBC_2.38 not found).
- The environment variables encourage static linking for transitive C libraries (OpenSSL, zlib, bzip2, lzma, zstd, curl).

### Using Cargo

```bash
cargo install inquiSTR
```

## Usage

The inquiSTR tool has several subcommands, as detailed below. All commands write to stdout.

```text
Usage: inquiSTR <COMMAND>

Commands:
  call         Call lengths
  combine      Combine STR calls or kmer frequencies from multiple samples to a TSV
  outlier      Find outliers from combined STR or kmer data
  query        Lookup genotypes and display
  histogram    Generate histograms for specific repeats
  plot         Show a histogram with multiple groups for a specific repeat
  pca          Perform Principal Component Analysis on combined STR data
  association  Perform statistical association testing for STRs
  unmapped     Count kmer frequencies in unmapped reads
  help         Print this message or the help of the given subcommand(s)

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

#### CRAM inputs: reference FASTA and FAI index

When using CRAM files, provide the matching reference with `--reference` and ensure a FASTA index (`.fai`) exists next to the FASTA. inquiSTR prefers reading contig names and lengths from the reference `.fai` rather than opening the CRAM early; this is faster and avoids edge-case crashes observed in some CRAM/reference setups.

- Required for CRAM: `--reference genome.fa`
- Recommended: a corresponding `genome.fa.fai` in the same directory
- Fallback: if `.fai` is missing, inquiSTR will fall back to contig lengths from the CRAM header (slower and less robust). We strongly recommend providing the `.fai`.

Create the index once with samtools:

```bash
samtools faidx genome.fa
```

Notes:

- The `.fai` must match the exact FASTA used to generate your CRAM.
- Keep the reference FASTA and its `.fai` side by side on local storage. Remote references are not supported for `.fai` loading.
- For extra diagnostics, you can enable debug logs: `RUST_LOG=debug`.

Examples with CRAM + reference:

```bash
# Call STRs in CRAM with provided reference (FAI required for best performance/stability)
inquiSTR call sample.cram --reference genome.fa -R regions.bed

# Single-region call on CRAM
inquiSTR call sample.cram --reference genome.fa -r chr1:1000-1100
```

#### Remote File Access

inquiSTR automatically handles remote BAM/CRAM files with index caching:

- **Auto-caching**: Downloaded indexes are cached at `~/.cache/inquistr/` for instant reuse
- **Local-first**: Checks current directory for index files before downloading
- **Auto-cleanup**: Removes cached files older than 30 days to prevent bloat
- **User control**: Set `INQUISTR_NO_CACHE=1` to disable caching for one-time analyses

See [INDEX_CACHING.md](INDEX_CACHING.md) for detailed configuration options.

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
  -s, --sample <SAMPLE>    sample(s) to consider: can be a single sample name, comma-separated sample names, 
                           or a file path containing sample names (one per line)
  -t, --threads <THREADS>  Number of threads to use for parallel processing [default: 1]
  -h, --help               Print help
```

**Sample Specification:**

The `--sample` option accepts three formats:

1. **Single sample name**: `--sample "Sample123"`
2. **Comma-separated names**: `--sample "Sample1,Sample2,Sample3"`
3. **File path**: `--sample samples.txt` (one sample name per line)

All specified samples are automatically validated against the combined file to ensure they exist.

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

# Find outliers for multiple samples (comma-separated)
inquiSTR outlier combined.tsv --sample "Sample1,Sample2,Sample3"

# Find outliers for samples listed in a file
inquiSTR outlier combined.tsv --sample samples.txt

# Custom minimum expansion size with sample filtering
inquiSTR outlier combined.tsv --minsize 15 --sample "PatientX"

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
  -k, --klength <KLENGTH>                Maximum kmer length to count (counts all sizes from 2 to klength) [default: 6]
      --sample-name <SAMPLE_NAME>        Sample name to use in output header
      --reference <REFERENCE>            Reference fasta for CRAM decoding
  -t, --threads <THREADS>                Number of parallel threads to use [default: 1]
      --target-kmer <TARGET_KMER>        Target kmer to specifically quantify (optional, can be any length). Supports shorthand notation: (CT)4 = CTCTCTCT
      --combine-revcomp                  Combine kmers with their reverse complements (e.g., CTCTCT and AGAGAG counted together)
  -h, --help                             Print help
```

**Key Features:**

- **Comprehensive kmer analysis**: Counts all kmer sizes from 2 to `--klength` (default 6)
- **Canonical representation**: All rotations of the same kmer are represented by the lexicographically smallest form
- **Target kmer search**: Optionally search for a specific kmer pattern with `--target-kmer`
- **Shorthand notation**: Specify repeating kmers easily: `(CT)4` = `CTCTCTCT`, `(CAG)10` = `CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG`
- **Reverse complement combining**: Use `--combine-revcomp` to count forward and reverse complement kmers together (e.g., CTCTCT and AGAGAG)
- **Complete output**: Includes all possible kmers (including zeros) sorted alphabetically by k-size
- **Normalized counts**: Results are normalized by total read count in the entire BAM/CRAM file
- **Parallel processing**: Batch-based parallel processing for efficient analysis of large files
- **Optimized search**: Uses hash-based lookup for O(1) kmer matching (especially fast for target kmer searches)
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

# Use multiple threads for faster processing
inquiSTR unmapped sample.bam --klength 5 --threads 8 --sample-name MySample

# CRAM file with reference
inquiSTR unmapped sample.cram --reference genome.fasta --klength 3 --sample-name CramSample

# Remote file analysis
inquiSTR unmapped https://example.com/sample.bam --threads 4 --sample-name RemoteSample

# Target specific kmer (traditional format)
inquiSTR unmapped sample.bam --target-kmer CTCTCTCT --sample-name MySample

# Target specific kmer with shorthand notation (equivalent to above)
inquiSTR unmapped sample.bam --target-kmer "(CT)4" --sample-name MySample

# Search for CAG repeats (common in neurodegenerative diseases)
inquiSTR unmapped sample.bam --target-kmer "(CAG)10" --threads 4

# Combine reverse complements (CTCTCT and AGAGAG counted together)
inquiSTR unmapped sample.bam --target-kmer "(CT)4" --combine-revcomp

# Full catalog with reverse complement combining
inquiSTR unmapped sample.bam --klength 6 --combine-revcomp --threads 8

# Test with provided example data
inquiSTR unmapped test-data/unmapped.bam --klength 3 --sample-name test_sample
```

CRAM note: as with `inquiSTR call`, providing `--reference genome.fa` with a matching `genome.fa.fai` is strongly recommended. Create the index with `samtools faidx genome.fa`.

**Shorthand Notation:**

For repeating sequences, use parentheses with a repeat count:

- `(CT)4` → `CTCTCTCT`
- `(AG)3` → `AGAGAG`
- `(CAG)5` → `CAGCAGCAGCAGCAG`
- `(ATCG)2` → `ATCGATCG`

Limits: Repeat count must be 1-100, unit must contain only A, C, G, T

**Example Output (Full Catalog):**

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

**Example Output (Target Kmer):**

```tsv
Sample          Target_Kmer     Canonical_Kmer  Kmer_Length     Count   Total_Reads     Frequency
MySample        CTCTCTCT        CTCTCT          8               1523    125000          0.012184
```

When `--combine-revcomp` is used, the canonical form represents the lexicographically smallest among all rotations of both the forward strand and reverse complement.

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

### `inquiSTR association` - Statistical Association Testing

Perform statistical association testing for STRs using the embedded R script. This subcommand provides a convenient wrapper around sophisticated statistical analysis capabilities without requiring separate script management.

```text
Usage: inquiSTR association [OPTIONS] --input <INPUT> --phenocovar <PHENOCOVAR> --phenotype <PHENOTYPE> --out <OUT> --str-mode <STR_MODE> --outcometype <OUTCOMETYPE>

Arguments:
  -i, --input <INPUT>              Combined STR file from inquiSTR combine command
  -p, --phenocovar <PHENOCOVAR>    Phenotype and covariate file with header, first column is individual ID
      --phenotype <PHENOTYPE>      Column name of phenotype in phenocovar file
  -o, --out <OUT>                  Output file name for association results
      --str-mode <STR_MODE>        STR mode: MEAN, MAX, or MIN for H1/H2 combination
      --outcometype <OUTCOMETYPE>  Outcome type: binary or continuous

Options:
      --covnames <COVNAMES>              Covariate names, comma separated (optional)
      --missing-cutoff <MISSING_CUTOFF>  Call rate cutoff for variants [default: 0.8]
      --minimal-length <MINIMAL_LENGTH>  Minimum maximum STR length across samples for inclusion
  -t, --threads <THREADS>                Number of threads for parallel processing [default: 1]
      --chunk-size <CHUNK_SIZE>          Number of variants to process in each chunk [default: 1000]
      --binary-order <BINARY_ORDER>      Binary phenotype order (e.g., Control,Patient) - required for binary outcomes
      --quiet                            Do not print progress messages
  -h, --help                             Print help
```

#### R Environment Setup

The association testing functionality requires R with specific packages. inquiSTR automatically checks your R environment and provides setup instructions if needed.

**Required R packages:**

- `data.table` - Fast data manipulation
- `argparser` - Command-line argument parsing
- `parallel` - Parallel processing (included with base R)

**Installation:**

```bash
# Install R packages
Rscript -e "install.packages(c('data.table', 'argparser'), repos='https://cran.rstudio.com/')"

# Or install interactively in R
R
> install.packages(c('data.table', 'argparser'))
```

#### Phenotype File Format

The phenotype file should be a tab-separated or comma-separated file with:

- **Header row** with column names
- **First column**: Individual IDs matching those in the combined STR file
- **Phenotype column**: The outcome variable
- **Additional columns**: Optional covariates (age, sex, population, etc.)

**Example for continuous outcome:**

```csv
individual_id,height,age,sex
sample1,175.2,25,M
sample2,168.1,30,F
sample3,182.5,28,M
```

**Example for binary outcome:**

```csv
individual_id,disease_status,age,population
patient1,Case,45,EUR
control1,Control,47,EUR
patient2,Case,52,AFR
```

#### STR Mode Options

- **MEAN**: Average of H1 and H2 allele lengths
- **MAX**: Maximum of H1 and H2 allele lengths (good for expansion detection)
- **MIN**: Minimum of H1 and H2 allele lengths

#### Examples

**Continuous phenotype analysis:**

```bash
inquiSTR association \
  --input combined_strs.tsv \
  --phenocovar phenotypes.csv \
  --phenotype height \
  --out height_association.tsv \
  --str-mode MAX \
  --outcometype continuous \
  --covnames age,sex \
  --threads 8
```

**Binary case-control analysis:**

```bash
inquiSTR association \
  --input combined_strs.tsv \
  --phenocovar case_control.csv \
  --phenotype disease_status \
  --out disease_association.tsv \
  --str-mode MEAN \
  --outcometype binary \
  --binary-order Control,Case \
  --covnames age,sex,population \
  --missing-cutoff 0.9 \
  --threads 4
```

**Advanced filtering:**

```bash
inquiSTR association \
  --input combined_strs.tsv \
  --phenocovar phenotypes.csv \
  --phenotype trait \
  --out filtered_association.tsv \
  --str-mode MAX \
  --outcometype continuous \
  --minimal-length 10 \
  --missing-cutoff 0.95 \
  --chunk-size 500 \
  --quiet
```

## Legacy R Script Usage

For advanced users who prefer direct script access, the association testing functionality is also available as a standalone R script `STR_regression.R` in the scripts folder. Detailed examples are provided in [STR_regression_examples.md](STR_regression_examples.md).

## 🛠️ Development

### Development Setup

```bash
# Clone the repository
git clone https://github.com/wdecoster/inquiSTR.git
cd inquiSTR

# Install development tools
make setup

# Install git hooks for automated code quality checks
make install-hooks
```

### Code Quality

This project uses automated code quality checks to maintain consistency and prevent CI failures:

- **Formatting**: `cargo fmt` automatically formats code according to Rust standards
- **Linting**: `cargo clippy` checks for common mistakes and suggests improvements
- **Testing**: `cargo test` runs the test suite

#### Git Hooks (Recommended)

The project includes pre-commit and pre-push hooks that automatically run quality checks:

```bash
# Install hooks (run once after cloning)
make install-hooks

# The hooks will now automatically:
# - Format code on commit (pre-commit)
# - Run clippy and formatting checks before push (pre-push)
```

#### Manual Quality Checks

You can also run quality checks manually:

```bash
# Format code
make fmt

# Run clippy linter
make clippy

# Run all pre-push checks
make pre-push

# Run all CI checks (formatting, linting, tests)
make ci
```

#### Available Make Targets

```bash
make help          # Show all available targets
make fmt           # Format code with cargo fmt
make clippy        # Run clippy linter with warnings as errors
make test          # Run tests
make build         # Build in release mode
make ci            # Run all CI checks (fmt-check, clippy, test)
make pre-push      # Run pre-push checks (fmt, clippy)
make install-hooks # Install git hooks for automated checks
```

### Why Use Git Hooks?

The most common cause of CI failures is code formatting and clippy warnings. The git hooks prevent these issues by:

1. **Pre-commit hook**: Automatically formats your code before each commit
2. **Pre-push hook**: Ensures code passes formatting and clippy checks before pushing

This saves time by catching issues locally instead of discovering them in CI.

### Contributing

1. **Fork and clone** the repository
2. **Install hooks**: `make install-hooks`
3. **Make changes** and commit (hooks will auto-format)
4. **Push changes** (hooks will run quality checks)
5. **Create a pull request**

The automated hooks ensure your contributions meet the project's quality standards before they reach CI.
