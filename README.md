# inquiSTR

[![CI](https://github.com/wdecoster/inquiSTR/actions/workflows/test.yml/badge.svg)](https://github.com/wdecoster/inquiSTR/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A toolkit for lightning-fast Tandem Repeat (TR) length genotyping and downstream analysis from long-read sequencing data. inquiSTR (/ɪnˈkwɪzɪtər/, pronounced like "inquisitor") works with Oxford Nanopore Technologies and PacBio BAM/CRAM files and supports both phased and unphased data. Additional subcommands are provided for combining results across samples for cohort studies with statistical association testing, outlier detection, relatedness assssment and visualization. inquiSTR provides seamless access to remote files via HTTP/HTTPS/FTP/S3 URLs and efficiently leverages multi-core systems through parallelized code.

## Quick Start

```bash
# 1. Genotype a single sample at pathogenic STR loci
inquiSTR call sample.bam --preset pathogenic --output sample.inq

# 2. Combine results for cohort analysis
inquiSTR combine *.inq --output combined.tsv
```

For detailed options, presets and additional commands, see the full documentation below.

## Table of Contents

- [Quick Start](#quick-start)
- [Installation](#-installation)
- [Usage](#usage)
  - [inquiSTR call - STR Genotyping](#inquistr-call---str-genotyping)
  - [inquiSTR combine - Multi-sample Analysis](#inquistr-combine---multi-sample-analysis)
  - [inquiSTR batch - Batch Sample Processing](#inquistr-batch---batch-sample-processing)
  - [inquiSTR convert - VCF to inquiSTR Format](#inquistr-convert---vcf-to-inquistr-format)
  - [inquiSTR filter - Filter STR Data](#inquistr-filter---filter-str-data)
  - [inquiSTR query - Genotype Lookup](#inquistr-query---genotype-lookup)
  - [inquiSTR outlier - Outlier Detection](#inquistr-outlier---outlier-detection)
  - [inquiSTR histogram - Data Visualization](#inquistr-histogram---data-visualization)
  - [inquiSTR plot - Group Comparison Plots](#inquistr-plot---group-comparison-plots)
  - [inquiSTR unmapped - Kmer Frequency Analysis](#inquistr-unmapped---kmer-frequency-analysis)
  - [inquiSTR benchmark - Validate STR Calls](#inquistr-benchmark---validate-str-calls)
  - [inquiSTR pca - Principal Component Analysis](#inquistr-pca---principal-component-analysis)
  - [inquiSTR relate - Compute Sample Relatedness](#inquistr-relate---compute-sample-relatedness)
  - [inquiSTR association - Statistical Association Testing](#inquistr-association---statistical-association-testing)
  - [inquiSTR optimize-call - Parameter Optimization](#inquistr-optimize-call---parameter-optimization)
- [Contributing](#-contributing)

## 📦 Installation

### Pre-built Binaries (Recommended)

Pre-built binaries for Linux and macOS can be downloaded from the [releases page](https://github.com/wdecoster/inquiSTR/releases). Make sure the binary is in your $PATH or use the full path to execute.

```bash
# Linux (musl - static binary)  
curl -L https://github.com/wdecoster/inquiSTR/releases/latest/download/inquiSTR-linux-musl -o inquiSTR
chmod +x inquiSTR

# Linux (glibc)
curl -L https://github.com/wdecoster/inquiSTR/releases/latest/download/inquiSTR-linux -o inquiSTR
chmod +x inquiSTR

# macOS
curl -L https://github.com/wdecoster/inquiSTR/releases/latest/download/inquiSTR-macos -o inquiSTR
chmod +x inquiSTR
```

### From Source

Building from source requires a C/C++ toolchain and common compression/network libraries used by HTSlib. Installation instructions are provided for Ubuntu/Debian, these provide standard headers for clang/bindgen (fixes errors like "stddef.h not found") and libraries used by BAM/CRAM support. The binary will be in target/release/inquiSTR after building.

```bash
sudo apt-get update
sudo apt-get install -y build-essential pkg-config clang libclang-dev \
  zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
  libdeflate-dev libzstd-dev libssl-dev
git clone https://github.com/wdecoster/inquiSTR.git
cd inquiSTR  
cargo build --release
```

To build a binary that does not depend on your system glibc (useful for older servers or containers), you can build with MUSL. The binary will be in target/x86_64-unknown-linux-musl/release/inquiSTR after building.

```bash
rustup target add x86_64-unknown-linux-musl
sudo apt-get install -y musl-tools
OPENSSL_STATIC=1 LIBZ_SYS_STATIC=1 BZIP2_STATIC=1 ZSTD_STATIC=1 LZMA_API_STATIC=1 CURL_STATIC=1 \
  cargo build --release --target x86_64-unknown-linux-musl
```

## Usage

The inquiSTR tool has several subcommands, as detailed below. All commands write to stdout.

```text
Usage: inquiSTR <COMMAND>

Commands:
  call         Call lengths
  combine      Combine STR calls or kmer frequencies from multiple samples to a TSV
  batch        Process multiple samples in batch and combine results
  convert      Convert VCF files to inquiSTR format
  filter       Filter inquiSTR output by various criteria
  outlier      Find outliers from combined STR or kmer data
  query        Lookup genotypes and display
  histogram    Generate histograms for specific repeats
  plot         Show a histogram with multiple groups for a specific repeat
  pca          Perform Principal Component Analysis on combined STR data
  relate       Compute relatedness between samples
  association  Perform statistical association testing for STRs
  optimize-call Optimize batch_size and thread count for your system and dataset
  unmapped     Count kmer frequencies in unmapped reads
  benchmark    Benchmark inquiSTR calls against truth VCF, BED, or another inquiSTR call file
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
      --preset <PRESET>            Use a predefined TR catalog: pathogenic, adotto, trexplorer, or codis
  -m, --minlen <MINLEN>            minimal length of insertion/deletion operation [default: 5]
  -s, --support <SUPPORT>          minimal number of supporting reads [default: 3]
  -t, --threads <THREADS>          Number of parallel threads to use [default: 1]
  -u, --unphased                   If reads have to be considered unphased
      --require-spanning           Only report genotypes supported by spanning reads; loci covered only by soft-clipped reads are reported as missing
      --sample-name <SAMPLE_NAME>  sample name to use in output
      --reference <REFERENCE>      reference fasta for cram decoding
      --max-locus <MAX_LOCUS>      maximum locus size to consider (intervals larger than this will be filtered out)
      --batch-size <BATCH_SIZE>    Batch size in KB for grouping nearby STR targets [default: 50]
      --vcf <VCF>                  Output VCF file path (optional, TSV still written to stdout)
  -h, --help                       Print help
```

When using CRAM files, provide the matching reference with `--reference genome.fa` and preferably have a FASTA index (`genome.fa.fai`) in the same directory. inquiSTR can more efficiently read contig names and lengths from the `genome.fa.fai` rather than opening the CRAM early. A fai index can be created with `samtools faidx genome.fa`

**Performance Note:**

If you are not sure what `--threads` and `--batch-size` values are optimal for your system you can use [`inquiSTR optimize-call`](#inquistr-optimize-call---parameter-optimization) to empirically determine the best configuration through automated benchmarking. inquiSTR typically scales well up to **4-6 threads** for single-sample processing. Beyond this, performance plateaus due to I/O bottlenecks (BAM reading and decompression). Using more than 6 threads per sample will most often not significantly improve speed.

**Examples:**

```bash
# Single region genotyping
inquiSTR call sample.bam -r chr1:1000-1100

# Multiple regions from BED file
inquiSTR call sample.bam -R regions.bed

# Use the predefined STRchive TR catalog (automatically downloads and caches, see below)
inquiSTR call sample.bam --preset pathogenic
```

#### Predefined TR Catalogs

The `--preset` option provides quick access to well-known TR catalogs without manually downloading BED files. All preset catalogs are currently for the **GRCh38/hg38** reference genome. Catalogs are automatically downloaded on first use and cached locally for 7 days in `~/.cache/inquistr/`. If a download fails but a cached version exists (even if expired), inquiSTR will use the cached version with a warning. Adding new preset catalogs is straightforward - please open an issue here on GitHub or see the [developer documentation](src/repeats.rs) for details on extending the `TRPreset` enum.

Available presets:

- **pathogenic**: STRchive pathogenic disease-associated STRs - curated database of STRs linked to human diseases (**75 loci**). See also [Hiatt et al., 2025](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-025-01454-4).
- **adotto**: Adotto TR regions catalog v1.2.1 - comprehensive TR regions from the Adotto project and benchmark (**1,784,804 loci**). See also [English et al., 2025](https://www.nature.com/articles/s41587-024-02225-z).
- **trexplorer**: Broad Institute TR Explorer catalog - genome-wide TR catalog covering 1-1000bp motifs (**4,863,041 loci**). See also [Weisburd et al., 2025](https://www.biorxiv.org/content/10.1101/2024.10.04.615514v2).
- **codis**: CODIS forensic STR markers from the USAT catalog - standard forensic STR markers used for human identification (**20 loci**). See also [Wang et al., 2024](https://link.springer.com/article/10.1186/s12859-022-05021-1).

### `inquiSTR combine` - Multi-sample Analysis

Combine data from multiple samples with `inquiSTR combine`. This command supports either STR call files (from [`inquiSTR call`](#inquistr-call---str-genotyping)) or kmer frequency files (from [`inquiSTR unmapped`](#inquistr-unmapped---kmer-frequency-analysis)), automatically detecting the input format. You can also add new samples to an existing combined file by providing both the combined file and new individual files. This enables efficient cohort expansion without reprocessing all samples. While the input files are validated, the goal is to use inquiSTR combine with files generated for the same loci or motifs.

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

# Add new samples to an existing combined file (incremental cohort building)
inquiSTR combine str_combined.tsv sample4.inq sample5.inq > expanded_cohort.tsv

# Merge multiple combined files together
inquiSTR combine cohort1_combined.tsv cohort2_combined.tsv > merged_cohorts.tsv
```

### `inquiSTR batch` - Batch Sample Processing

Process multiple samples in batch and automatically combine the results, eliminating the need to manually run [`inquiSTR call`](#inquistr-call---str-genotyping) (or [`inquiSTR unmapped`](#inquistr-unmapped---kmer-frequency-analysis)) followed by [`inquiSTR combine`](#inquistr-combine---multi-sample-analysis). STR genotyping mode is the default, using the same options as [`inquiSTR call`](#inquistr-call---str-genotyping). InquiSTR unmapped mode can be activated with the `--unmapped` flag, using the same options as [`inquiSTR unmapped`](#inquistr-unmapped---kmer-frequency-analysis). Temporary individual sample files are stored in the directory specified with `--tmpdir`, specified by the `$TMPDIR` environmental variable, or the current directory and are cleaned up unless `--save-individual` is specified.

```text
Usage: inquiSTR batch [OPTIONS] --output <OUTPUT> <MANIFEST>

Arguments:
  <MANIFEST>  TSV manifest file with bam_path (required) and sample_name (optional) columns

Common Options (both modes):
  -o, --output <OUTPUT>                  Output file for combined results
  -t, --threads <THREADS>                Number of parallel threads to use (per sample) [default: 1]
      --reference <REFERENCE>            Reference fasta for CRAM decoding (applies to all samples)
      --save-individual <DIR>            Save individual sample files to this directory (optional)
      --tmpdir <TMPDIR>                  Temporary directory for intermediate files
      --resume                           Skip already processed samples by checking the output file
      --dry-run                          Validate manifest and preview processing without running

Mode Selection:
      --unmapped                         Process unmapped reads instead of genotyping STRs

Unmapped Kmer Mode Options (only with --unmapped):
  -k, --klength <KLENGTH>                Maximum kmer length to count [default: 6]
      --target-kmer <TARGET_KMER>        Target kmer to quantify (optional, supports shorthand notation)
      --combine-revcomp                  Combine kmers with reverse complements

STR Genotyping Mode Options (default mode, without --unmapped):
  -r, --region <REGION>                  Region string to genotype expansion in
  -R, --region-file <REGION_FILE>        Bed file with region(s) to genotype expansion(s) in
      --preset <PRESET>                  Use a predefined TR catalog (pathogenic, adotto, trexplorer, or codis)
  -m, --minlen <MINLEN>                  Minimal length of insertion/deletion operation [default: 5]
  -s, --support <SUPPORT>                Minimal number of supporting reads [default: 3]
  -u, --unphased                         If reads have to be considered unphased
      --require-spanning                 Only report genotypes supported by spanning reads; loci covered only by soft-clipped reads are reported as missing
      --max-locus <MAX_LOCUS>            Maximum locus size to consider
      --batch-size <BATCH_SIZE>          Batch size in KB for grouping nearby STR targets [default: 50]

  -h, --help                             Print help
```

The manifest file is a tab-separated values (TSV) file with a header row. The bam_path is required, can be a local path or remote URL, and sample_name is optional (if omitted, it will be extracted from the BAM filename).

```tsv
bam_path    sample_name
/path/to/sample1.bam    Patient01
/path/to/sample2.bam    Patient02
/path/to/sample3.cram    Control01
```

**Examples:**

```bash
# Basic batch processing with pathogenic STRs
inquiSTR batch samples.tsv --preset pathogenic --output cohort_results.tsv

# Batch processing with custom regions and save individual files
inquiSTR batch samples.tsv -R regions.bed --output combined.tsv --save-individual results/

# Use custom temporary directory
inquiSTR batch samples.tsv --preset adotto --tmpdir /scratch/tmp --output results.tsv

# Basic unmapped kmer analysis across cohort
inquiSTR batch samples.tsv --unmapped --output kmer_combined.tsv
```

**Workflow Features:**

The `--resume` and `--dry-run` flags provide workflow management capabilities:

- `--dry-run`: Validates the manifest file, checks for missing BAM files, and reports what would be processed without actually running analysis.
- `--resume`: Intelligently resumes interrupted batch jobs by checking the output file for already-processed samples:
and can be used to restart a failed workflow or add new samples to an existing cohort incrementally.
- Use both `--dry-run` and `--resume` together to preview which samples would be skipped and which would be processed, without running anything. This combination is useful for planning before resuming an interrupted job.

### `inquiSTR convert` - VCF to inquiSTR Format

Convert VCF files from other STR genotyping tools to inquiSTR's TSV format. This enables use of inquiSTR's downstream analysis tools ([`inquiSTR outlier`](#inquistr-outlier---outlier-detection), [`inquiSTR association`](#inquistr-association---statistical-association-testing), [`inquiSTR pca`](#inquistr-pca---principal-component-analysis), etc.) with data from other callers. Please let me know if the tool does not work for a specific type of VCF, and I will look into it. Allele lengths are calculated as length(ALT) - length(REF), therefore negative values indicate contractions. Missing genotypes (./.) are converted to "NA".

```text
Usage: inquiSTR convert <VCF>...

Arguments:
  <VCF>...  VCF file(s) to convert (can be compressed). Single sample VCF produces 
            individual_call file, multiple samples or multiple VCFs produce combined_call file

Options:
  -h, --help  Print help
```

**Examples:**

```bash
# Convert a single-sample VCF to inquiSTR format
inquiSTR convert sample1.vcf > sample1.inq

# Convert a multi-sample VCF to combined format
inquiSTR convert cohort.vcf > cohort_combined.tsv

# Combine multiple single-sample VCFs into one file
inquiSTR convert sample1.vcf sample2.vcf sample3.vcf > combined.tsv
```

### `inquiSTR optimize-call` - Parameter Optimization

Empirically determine the optimal `--batch-size` and `--threads` configuration for running [`inquiSTR call`](#inquistr-call---str-genotyping) on your specific system, input file, and repeat catalog. This subcommand provides data-driven recommendations with visualizations. Each test configuration typically completes in 2-5 minutes, and the total optimization time is about ~1-2 hours for comprehensive testing of combinations of parameters. For faster testing, provide a BED file with only one chromosome, rather than random downsampling targets to keep the batching approach comparable to a real dataset.

```text
Usage: inquiSTR optimize-call [OPTIONS] <BAM>

Arguments:
  <BAM>  BAM/CRAM file to optimize parameters for

Options:
  -R, --region-file <REGION_FILE>          Bed file with region(s) to test
      --preset <PRESET>                     Use a predefined TR catalog (pathogenic, adotto, trexplorer, or codis)
      --reference <REFERENCE>               Reference fasta for CRAM decoding
      --min-threads <MIN_THREADS>           Minimum number of threads to test [default: 1]
      --max-threads <MAX_THREADS>           Maximum number of threads to test [default: 16]
      --batch-sizes <BATCH_SIZES>           Batch sizes to test (comma-separated, in KB) [default: 5,10,20,30,50]
      --repeats <REPEATS>                   Number of repetitions per configuration [default: 3]
  -o, --output <OUTPUT>                     Output directory for results and plots [default: optimize_results]
  -h, --help                                Print help
```

**Output Files:**

- `benchmark_results.tsv` - Raw timing data for all test runs
- `aggregated_stats.tsv` - Statistical summaries (mean, std dev, min times)
- `recommendation.txt` - Human-readable recommendation with rationale
- `wall_time_heatmap.html` - Interactive heatmap showing performance across configurations
- `optimization_analysis.html` - 4-panel analysis (speedup, efficiency, batch size trends)

**Examples:**

```bash
# Basic optimization
inquiSTR optimize-call sample.cram \
    -R regions.bed \
    --reference genome.fa

# Comprehensive test with custom output directory
inquiSTR optimize-call sample.bam \
    --preset adotto \
    --min-threads 1 \
    --max-threads 12 \
    --batch-sizes 5,10,15,20,30,50 \
    --repeats 3 \
    --output my_optimization_results
```

**Example Output:**

```text
===========================================
OPTIMIZATION RESULTS
===========================================

🎯 RECOMMENDED CONFIGURATION:
   --threads 4  --batch-size 10
   Expected runtime: 88.3s
   Confidence: High

   Fastest configuration: 88.3s wall time. Uses 4/12 CPUs.
   More threads showed diminishing returns. Small batch size
   provides fine-grained parallelism.

📊 TOP 5 CONFIGURATIONS:
   Rank  Threads  Batch   Time(s)
   ----  -------  -----   -------
      1        4   10KB      88.3
      2        8   10KB     115.6
      3        4   20KB     161.9
      4        8   20KB     163.0
      5        2   10KB     168.2

✓ Optimization complete!

Recommended command:
  inquiSTR call <input> -R <regions> --threads 4 --batch-size 10

Visualization files saved in: optimize_results
```

### `inquiSTR query` - Genotype Lookup

Querying genotypes from a combined file (from [`inquiSTR combine`](#inquistr-combine---multi-sample-analysis)) can be done with `inquiSTR query`, taking a region string or a file with regions to query. When querying a single locus, output will be sorted by repeat length (longest expansions first).

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
inquiSTR query combined.tsv chr1:1000-1200

# Query multiple regions from a BED file (multi-locus table format)
inquiSTR query combined.tsv regions.bed
```

### `inquiSTR filter` - Filter STR Data

Filter inquiSTR output files by various criteria to focus on variants of interest. This command works with output from [`inquiSTR call`](#inquistr-call---str-genotyping) or [`inquiSTR combine`](#inquistr-combine---multi-sample-analysis).

```text
Usage: inquiSTR filter [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Input file from inquiSTR call or inquiSTR combine

Options:
      --minlen <MINLEN>            Minimum allele length (expansion/insertion only, ignores negative values)
      --minchange <MINCHANGE>      Minimum absolute allele length change (uses absolute value)
      --bed <BED>                  BED file to filter by overlap (requires both files to be sorted)
      --call-rate <CALL_RATE>      Minimum call rate (fraction 0.0-1.0) for combined files
      --min-cv <MIN_CV>            Minimum coefficient of variation (only for combined files)
  -s, --samples <SAMPLES>          Sample(s) to keep: can be a single sample name, comma-separated sample names, 
                                   or a file path containing sample names (one per line)
  -d, --drop-samples <DROP_SAMPLES>
                                   Sample(s) to drop: can be a single sample name, comma-separated sample names, 
                                   or a file path containing sample names (one per line)
  -h, --help                       Print help
```

**Examples:**

```bash
# Filter for STR expansions of at least 20 bp
inquiSTR filter combined.tsv --minlen 20 > expansions_20bp.tsv

# Filter for variants with absolute change of at least 10 bp within exons defined in exons.bed
inquiSTR filter combined.tsv --minchange 10 --bed exons.bed > changes_10bp_coding.tsv

# Combine multiple filters, e.g. expansions of at least 15 bp with call rate of at least 90% and coefficient of variation of at least 0.05
inquiSTR filter combined.tsv --minlen 15 --call-rate 0.9 --min-cv 0.05 > filtered.tsv
```

### `inquiSTR outlier` - Outlier Detection

Identify outliers from combined data using `inquiSTR outlier`. This command works with both STR call data (from [`inquiSTR call`](#inquistr-call---str-genotyping)) and kmer frequency data (from [`inquiSTR unmapped`](#inquistr-unmapped---kmer-frequency-analysis)) (automatically detecting the format), using either z-scores or DBSCAN algorithms.

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

All specified samples are automatically validated against the combined file to ensure they exist. The `--sample` option accepts three formats either a specific single sample name, a comma-separated list of sample names, or a file path containing sample names (one per line).

**Examples:**

```bash
# Find outliers in STR call data using default z-score method with a minimum expansion size
inquiSTR outlier str_combined.tsv --minsize 100

# Find outliers in kmer frequency data using a z-score cutoff of 2.5 for a few specific samples only
inquiSTR outlier kmer_combined.tsv --zscore 2.5 --sample "Sample1,Sample2,Sample3"
```

### `inquiSTR histogram` - Data Visualization

Generate histograms for STR length distributions at specific loci from combined data (from [`inquiSTR combine`](#inquistr-combine---multi-sample-analysis)).

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
inquiSTR histogram combined.tsv chr4:3074877-3074933
```

### `inquiSTR plot` - Group Comparison Plots  

Show a histogram with multiple groups for a specific repeat (from [`inquiSTR combine`](#inquistr-combine---multi-sample-analysis)), useful for comparing cohorts.

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
inquiSTR plot combined.tsv metadata.tsv chr1:1000-2000 --condition "group:PAT,CON" --output patient_control.html
```

### `inquiSTR unmapped` - Kmer Frequency Analysis

Count kmer frequencies in unmapped reads from BAM/CRAM files. This analysis extracts all unmapped reads and counts occurrences of kmers of all sizes from 2 to a specified maximum length (default: 6). All kmer rotations (e.g., CAG, AGC, GCA) are canonicalized to a consistent form (AGC, the lexicographically first rotation). You can also search for a specific kmer pattern with `--target-kmer`, using shorthand notation: `(CT)4 = CTCTCTCT`, `(CAG)10 = CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG`. Use `--combine-revcomp` to count forward and reverse complement kmers together (e.g., CTCTCT and AGAGAG). The output is a TSV file with a kmer column and a column with normalized kmer frequency counts for the specified sample (normalized by total read count in the entire BAM/CRAM file).

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

**Examples:**

```bash
# Basic kmer analysis with default settings (k=6, 1 thread), combining rotations but not reverse complements
inquiSTR unmapped sample.bam > kmers.tsv

# Full catalog with combining the reverse complement, using 8 threads for faster processing
inquiSTR unmapped sample.bam --combine-revcomp --threads 8

# Include longer kmer lengths and a sample name, from cram format
inquiSTR unmapped sample.cram --klength 12 --sample-name MySample --reference genome.fasta > kmers_k12.tsv

# Target specific kmer with shorthand notation and combine with the reverse complement 
inquiSTR unmapped sample.bam --target-kmer "(CT)4" --sample-name MySample --combine-revcomp
```

**Shorthand Notation:**

For repeating sequences, use parentheses with a repeat count. Repeat count must be 1-100, unit must contain only A, C, G, T characters. Examples:

- `(CT)4` → `CTCTCTCT`
- `(AG)3` → `AGAGAG`
- `(CAG)5` → `CAGCAGCAGCAGCAG`
- `(ATCG)2` → `ATCGATCG`

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
...     ...
```

Note that only canonical forms are shown. Rotations like CA, GA, GC, TA, TC, TG are represented by their canonical forms (AC, AG, CG, AT, CT, GT respectively) and their counts are combined.

### `inquiSTR benchmark` - Validate STR Calls

Benchmark inquiSTR genotyping results (from [`inquiSTR call`](#inquistr-call---str-genotyping)) against truth data. Supported truth formats:

- **VCF** (`.vcf` or `.vcf.gz`): standard variant call format
- **inquiSTR call file** (`.inq`): another file produced by `inquiSTR call` or `inquiSTR convert`, detected via metadata
- **BED**: 9-column adotto-style BED file where the last two columns are haplotype lengths

This command generates correlation statistics and an optional interactive HTML plot.

```text
Usage: inquiSTR benchmark [OPTIONS] --test <TEST> --truth <TRUTH>

Options:
      --test <TEST>                       inquiSTR call output file (.inq format) to benchmark
      --truth <TRUTH>                     Truth file (VCF, BED, or inquiSTR call file — format auto-detected)
  -m, --mode <MODE>                       Mode for selecting alleles: MAX (default) or MIN [default: MAX]
  -p, --plot <PLOT>                       Output file for correlation plot (HTML)
      --max-plot-length <MAX_PLOT_LENGTH>  Maximum allele length to display on plot [default: 5000]
      --tier1                             Only use Tier1 variants from BED file
  -d, --diff-out <DIFF_OUT>               Output file listing the largest discrepancies
      --max-locus <MAX_LOCUS>             Exclude loci larger than this size (bp) from truth data
      --nonzero                           Exclude zero-zero pairs (loci unchanged in both) from correlation
      --tolerance <TOLERANCE>             Tolerance in bp for counting calls as matching [default: 5]
  -h, --help                              Print help
```

**Output:**

- Correlation statistics (R², exact match rate, within-tolerance rate) printed to stdout
- Optional interactive HTML plot (`--plot`)
- Optional TSV of the largest discrepancies (`--diff-out`)

### `inquiSTR pca` - Principal Component Analysis

Perform Principal Component Analysis (PCA) on combined STR data (from [`inquiSTR combine`](#inquistr-combine---multi-sample-analysis)) to identify population structure and sample relationships. This is particularly useful for large-scale genomic analyses and quality control. The `--aggregation` parameter controls how the two allele lengths (H1/H2) are combined for each sample. You can choose from the max, min or sum of the two alleles.

```text
Usage: inquiSTR pca [OPTIONS] <COMBINED>

Arguments:
  <COMBINED>  Combined file from inquiSTR combine (supports both STR calls and kmer frequencies)

Options:
  -o, --output <OUTPUT>            HTML output file name for interactive PCA plot [default: pca_plot.html]
  -c, --components <COMPONENTS>    Number of principal components to compute (currently only first 2 are plotted) [default: 10]
  -t, --threads <THREADS>          Number of threads to use for parallel processing [default: 1]
  -a, --aggregation <AGGREGATION>  Method for aggregating H1/H2 allele lengths: max (default), min, or sum (only applies to STR files, ignored for kmer files) [default: max]
      --scores <SCORES>            Output file for PC scores (tab-separated: sample, PC1, PC2, ...) to use as covariates
  -h, --help                       Print help
```

**Examples:**

```bash
# Basic PCA analysis
inquiSTR pca combined.tsv

# Compute more principal components (only first 2 are plotted) 
# and save PC scores to a file for use as covariates in association testing
inquiSTR pca combined.tsv --components 20 --output pca.html --scores pc_scores.tsv

# PCA analysis on kmer frequencies from unmapped reads
inquiSTR pca inquiSTR_unmapped.tsv --output kmer_pca.html
```

**PC Scores Output:**

When using the `--scores` option, PC scores (sample coordinates in PC space) are written to a tab-separated file. These PC scores can be used as covariates in association testing to correct for population structure, or for further analysis in external tools.

```text
sample     PC1         PC2         PC3         ...
sample1    0.123456   -0.234567    0.345678   ...
sample2   -0.456789    0.567890   -0.678901   ...
```

### `inquiSTR relate` - Compute Sample Relatedness

Compute relatedness (kinship coefficient) between all pairs of samples using STR genotype data (from [`inquiSTR combine`](#inquistr-combine---multi-sample-analysis)). This command identifies sample relationships such as duplicates, parent-child, siblings, and more distant relatives by calculating the proportion of shared alleles across all STR loci, enabling users to validate pedigrees or identify cryptic relatedness in a cohort.

```text
Usage: inquiSTR relate [OPTIONS] --output <OUTPUT> <COMBINED>

Arguments:
  <COMBINED>  Combined file of STR calls from inquiSTR combine

Options:
  -o, --output <OUTPUT>    Output file for relatedness matrix
  -t, --threads <THREADS>  Number of threads to use for parallel processing [default: 1]
  -h, --help               Print help
```

**Method:**

The relatedness coefficient is computed using Identity-by-State (IBS):

- For each locus, count matching alleles between two samples (0, 1, or 2 matches)
- Relatedness = (IBS2 + 0.5 × IBS1) / n_loci
- Only informative loci (where both samples have valid calls) are used

**Expected Relatedness Values:**

- **1.0**: Identical twins / duplicate samples
- **0.5**: Parent-child / full siblings
- **0.25**: Half-siblings / grandparent-grandchild
- **0.125**: First cousins
- **0.0**: Unrelated individuals

**Output Format:**

TSV file with columns:

- `sample1`, `sample2`: Sample pair names
- `relatedness`: Coefficient (0-1), sorted descending
- `n_loci`: Number of informative loci used
- `ibs0`, `ibs1`, `ibs2`: Counts of loci with 0, 1, or 2 shared alleles

**Example:**

```bash
inquiSTR relate combined.tsv --output relatedness.tsv
```

### `inquiSTR association` - Statistical Association Testing

Perform statistical association testing for STRs from a combined file (from [`inquiSTR combine`](#inquistr-combine---multi-sample-analysis)).

```text
Usage: inquiSTR association [OPTIONS] --input <INPUT> --phenocovar <PHENOCOVAR> --phenotype <PHENOTYPE> --out <OUT> --str-mode <STR_MODE> --outcometype <OUTCOMETYPE>

Arguments:
  -i, --input <INPUT>              Combined STR file from inquiSTR combine
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
      --plot <PREFIX>                    Generate QQ plot and Manhattan plot with custom filename prefix (optional)
  -h, --help                             Print help
```

#### R Environment Setup

The association testing functionality requires R with specific packages (`data.table`, `argparser`, and `qqman` for `--plot`). inquiSTR automatically checks your R environment and provides setup instructions if needed.

**Installation:**

```bash
# Install R packages (including optional qqman for plotting)
Rscript -e "install.packages(c('data.table', 'argparser', 'qqman'), repos='https://cran.rstudio.com/')"

# Or install interactively in R
R
> install.packages(c('data.table', 'argparser', 'qqman'))

# Minimal installation (without plotting capability)
Rscript -e "install.packages(c('data.table', 'argparser'), repos='https://cran.rstudio.com/')"
```

#### Phenotype File Format

The phenotype file should be a tab-separated file (TSV, though comma-separated files are also supported) with a header row. The first column must contain individual IDs that match those in the combined STR file. One column should be designated as the phenotype (outcome variable), and additional columns can be included for optional covariates such as age, sex, population, etc.

**Example for continuous outcome:**

```tsv
individual_id	height	age	sex
sample1	175.2	25	M
sample2	168.1	30	F
sample3	182.5	28	M
```

**Example for binary outcome:**

```tsv
individual_id	disease_status	age	population
patient1	Case	45	EUR
control1	Control	47	EUR
patient2	Case	52	AFR
```

#### STR Mode Options

- **MEAN**: Average of H1 and H2 allele lengths
- **MAX**: Maximum of H1 and H2 allele lengths
- **MIN**: Minimum of H1 and H2 allele lengths

#### Examples

**Binary case-control analysis:**

```bash
inquiSTR association \
  --input combined_strs.tsv \
  --phenocovar case_control.tsv \
  --phenotype disease_status \
  --out disease_association.tsv \
  --str-mode MEAN \
  --outcometype binary \
  --binary-order Control,Case \
  --covnames age,sex,population \
  --missing-cutoff 0.9 \
  --threads 4
```

**Continuous phenotype analysis:**

```bash
inquiSTR association \
  --input combined_strs.tsv \
  --phenocovar phenotypes.tsv \
  --phenotype height \
  --out height_association.tsv \
  --str-mode MAX \
  --outcometype continuous \
  --covnames age,sex \
  --threads 8
```

**Advanced filtering with visualization:**

When `--plot` is used, inquiSTR will generate two additional files, a Manhattan plot showing -log10(p-value) across chromosomes and a QQ plot to assess genomic inflation. If no prefix is provided to `--plot`, it defaults to the output filename stem (e.g., `results.tsv` → `results_manhattan.png`)

```bash
inquiSTR association \
  --input combined_strs.tsv \
  --phenocovar phenotypes.tsv \
  --phenotype trait \
  --out filtered_association.tsv \
  --str-mode MAX \
  --outcometype continuous \
  --minimal-length 10 \
  --missing-cutoff 0.95 \
  --chunk-size 500 \
  --quiet \
  --plot myproject_disease
```

## 🤝 Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for development setup, code quality guidelines, and instructions on how to contribute.
