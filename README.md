# inquiSTR
Genotyping of STRs with long reads

This repository contains a snakefile and scripts to genotype specific STRs (in the data folder) from long read sequencing data, and has been optimized for ONT data.

## USAGE
```
snakemake -s /~path~/multi_repeat_typer.smk --cores 24 --use-conda
```