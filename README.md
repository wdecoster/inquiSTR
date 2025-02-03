# inquiSTR

This repository contains Rust code for inquiSTR, a toolset to genotype and analyze STRs from long read sequencing data, and has been tested with ONT data.

## Installation

Preferably, for most users, download a ready-to-use binary for your system to add directory on your $PATH from the [releases](https://github.com/wdecoster/inquiSTR/releases).  
You may have to change the file permissions to execute it with chmod +x inquiSTR

Alternatively, you can install the tool using cargo:

```bash
git clone https://github.com/wdecoster/inquiSTR.git
cd STRdust
cargo build --release
```

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

Variants from multiple samples can be combined with `inquiSTR combine`.

```text
inquiSTR combine <CALLS>...

Arguments:
  <CALLS>...  files from inquiSTR call

Options:
  -h, --help  Print help
```

Querying genotypes from a combined file can be done with `inquiSTR query`, taking a region or a file with regions to query.

```text
Usage: inquiSTR query <COMBINED> <REGION>

Arguments:
  <COMBINED>  combined file of calls
  <REGION>    region to query or file with regions to query

Options:
  -h, --help  Print help
```

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

## Usage for Association Testing

This repository additionally contains `str_assoc.R`, code to perform association testing of STRs. The code is written in R and can be found in the scripts folder.

Below are some worked usage examples for "MAX" STRmode for binary phenotypes, without covariates.

In the input test file, number of samples were 268, total number of variants were > 654K (about 279K passing quality filters).

### Full (genome-wide) run  

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run full --out full_genome_wide_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```

Though, this may take days, on our testing about 20% of variants were tested in a day. It would be more optimal to take advantage of processes like GNU parallel in combination with --run chromosome mode, for instance:

```bash
seq 1 1 22 | parallel --tmpdir ~/tmp/ --noswap --progress --eta -j 7 'Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run chromosome --chr chr{} --out chr{}.genome_wide_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT'
```

This time it took about 20-25 minutes to run all 22 chromosomes.

### Chromosome-wide run  

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run chromosome --chr chr15 --out chr15_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT'
```

Time:

| type | time       |
|------|------------|
| real | 3m8.689s   |
| user | 13m21.490s |
| sys  | 0m12.138s  |

### Chromosome interval run

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run chr_interval --chr chr15 --chr_begin 34419410 --chr_end 34419465 --out chr15_34419410_34419465_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```

Time:

| type | time       |
|------|------------|
| real | 0m33.452s  |
| user | 12m18.276s |
| sys  | 0m10.538s  |

### Bed interval run

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run bed_interval --bed chr15_roi.bed --out bed_chr15_roi_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```

Time:

| type | time       |
|------|------------|
| real | 0m29.313s  |
| user | 11m19.566s |
| sys  | 0m9.845s   |

### Single variant (Expanded Allele) run

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run single_variant --single_variant chr15_34419414_34419461 --expandedAllele 201 --out singleVariant_chr15_34419414_34419461_expandedAllele201_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```

or

```bash
Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run single_variant --single_variant chr15:34419414-34419461 --expandedAllele 201 --out singleVariant_chr15_34419414_34419461_expandedAllele201_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,PAT
```

Time:

| type | time       |
|------|------------|
| real | 0m37.876s  |
| user | 11m44.990s |
| sys  | 0m12.291s  |
