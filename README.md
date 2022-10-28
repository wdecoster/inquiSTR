# inquiSTR

Genotyping of STRs with long reads

This repository contains Rust code for inquiSTR, a tool to genotype STRs from long read sequencing data, and has been tested with ONT data.

## Usage for Association Testing

Below are some worked usage examples for "MAX" STRmode for binary phenotypes, without covariates.

In the input test file, number of samples were 268, total number of variants were > 654K (about 279K passing quality filters).

### Full (genome-wide) run  

'Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run full --out full_genome_wide_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,aFTLD-U'

Though, this may take days, on our testing about 20% of variants were tested in a day. It would be more optimal to take advantage of processes like GNU parallel in combination with --run chromosome mode, for instance:

''seq 1 1 22 | parallel --tmpdir ~/tmp/ --noswap --progress --eta -j 7 'Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run chromosome --chr chr{} --out chr{}.genome_wide_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,aFTLD-U''

This time it took about 20-25 minutes to run all 22 chromosomes.

### Chromosome-wide run  

'Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run chromosome --chr chr15 --out chr15_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,aFTLD-U'

Time:

| type | time       |
|------|------------|
| real | 3m8.689s   |
| user | 13m21.490s |
| sys  | 0m12.138s  |

### Chromosome interval run

'Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run chr_interval --chr chr15 --chr_begin 34419410 --chr_end 34419465 --out chr15_34419410_34419465_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,aFTLD-U'

Time:

| type | time       |
|------|------------|
| real | 0m33.452s  |
| user | 12m18.276s |
| sys  | 0m10.538s  |

### Bed interval run

'Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run bed_interval --bed chr15_roi.bed --out bed_chr15_roi_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,aFTLD-U'

Time:

| type | time       |
|------|------------|
| real | 0m29.313s  |
| user | 11m19.566s |
| sys  | 0m9.845s   |

### Single variant (Expanded Allele) run

'Rscript STR_regression.R --input combined.inq.gz --phenocovar inquistr-samples.tsv --phenotype group --run single_variant --single_variant chr15_34419414_34419461 --expandedAllele 201 --out singleVariant_chr15_34419414_34419461_expandedAllele201_testResults.tsv --STRmode MAX --outcometype binary --binaryOrder CON,aFTLD-U'

Time:

| type | time       |
|------|------------|
| real | 0m37.876s  |
| user | 11m44.990s |
| sys  | 0m12.291s  |
