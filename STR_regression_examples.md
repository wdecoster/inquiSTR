# STR regression examples

Below are some examples for "MAX" STRmode for binary phenotypes, without covariates.

## Full (genome-wide) run  

```bash
Rscript STR_regression.R \
 --input combined.inq.gz \
 --phenocovar inquistr-samples.tsv \
 --phenotype group \
 --run full 
 --out full_genome_wide_testResults.tsv \
 --STRmode MAX \
 --outcometype binary \
 --binaryOrder CON,PAT
```

## Chromosome-wide run  

```bash
Rscript STR_regression.R \
 --input combined.inq.gz \
 --phenocovar inquistr-samples.tsv \
 --phenotype group \
 --run chromosome \
 --chr chr15 \
 --out chr15_testResults.tsv \
 --STRmode MAX \
 --outcometype binary \
 --binaryOrder CON,PAT
```

## Chromosome interval run

```bash
Rscript STR_regression.R \
 --input combined.inq.gz \
 --phenocovar inquistr-samples.tsv \
 --phenotype group \
 --run chr_interval \
 --chr chr15 \
 --chr_begin 34419410 \
 --chr_end 34419465 \
 --out chr15_34419410_34419465_testResults.tsv \
 --STRmode MAX \
 --outcometype binary \
 --binaryOrder CON,PAT
```

## Bed interval run

```bash
Rscript STR_regression.R \
 --input combined.inq.gz \
 --phenocovar inquistr-samples.tsv \
 --phenotype group \
 --run bed_interval \
 --bed chr15_roi.bed \
 --out bed_chr15_roi_testResults.tsv \
 --STRmode MAX \
 --outcometype binary \
 --binaryOrder CON,PAT
```

## Single variant (Expanded Allele) run

```bash
Rscript STR_regression.R \
 --input combined.inq.gz \
 --phenocovar inquistr-samples.tsv \
 --phenotype group \
 --run single_variant \
 --single_variant chr15_34419414_34419461 \
 --expandedAllele 201 \
 --out singleVariant_chr15_34419414_34419461_expandedAllele201_testResults.tsv \
 --STRmode MAX \
 --outcometype binary \
 --binaryOrder CON,PAT
```

or

```bash
Rscript STR_regression.R \
 --input combined.inq.gz \
 --phenocovar inquistr-samples.tsv \
 --phenotype group \
 --run single_variant \
 --single_variant chr15:34419414-34419461 \
 --expandedAllele 201 \
 --out singleVariant_chr15_34419414_34419461_expandedAllele201_testResults.tsv \
 --STRmode MAX \
 --outcometype binary \
 --binaryOrder CON,PAT
```
