# Output File Formats

inquiSTR produces tab-separated text files. The two primary formats are the **individual call file** (one sample, produced by `inquiSTR call`) and the **combined call file** (multiple samples, produced by `inquiSTR combine` or `inquiSTR batch`).

## Individual call file (`.inq`)

Each file starts with metadata comment lines, followed by a column header and the data rows.

### Metadata header

```
# file_type=individual_call
# version=<inquiSTR version>
# command=call
# sample=<sample name>
# minlen=<minimum STR length threshold>
# support=<minimum read support threshold>
# unphased=<true|false>
# haploid=<comma-separated chromosome names, or None>
# catalog=<preset name or catalog file name>
# catalog_sha256=<SHA-256 of the decompressed catalog>
# genome_sha256=<SHA-256 of the chromosome length map>
```

The last three identify *what* was genotyped, by content rather than by name. `catalog_sha256`
is taken over the decompressed catalog, so a recompressed copy of the same file still matches
and the value can be checked by hand with `zcat catalog.bed.gz | sha256sum`. `genome_sha256` is
taken over the sorted chromosome length map from the BAM header, which distinguishes reference
builds without needing the FASTA. Files produced before inquiSTR recorded provenance carry none
of the three.

### Columns

```
chromosome  begin  end  info  <sample>_H1  <sample>_H2
```

| Column | Description |
|--------|-------------|
| `chromosome` | Chromosome name |
| `begin` | 0-based start coordinate (from the input BED file) |
| `end` | 0-based end coordinate (from the input BED file) |
| `info` | Fourth column of the input BED file, or `.` if absent |
| `<sample>_H1` | STR length change relative to reference for haplotype 1 (bp) |
| `<sample>_H2` | STR length change relative to reference for haplotype 2 (bp) |

The allele length values (`H1`, `H2`) are the median length difference from the reference genome in base pairs, computed from all supporting reads assigned to that haplotype. Positive values indicate expansion, negative values indicate contraction. A value of `NaN` means insufficient reads met the support threshold, no phased reads were available for that haplotype, or the locus is on a haploid chromosome (in which case `H1` carries the single allele value and `H2` is always `NaN`).

## Combined call file

The combined file merges results from multiple individual call files.

### Metadata header

```
# file_type=combined_call
# version=<inquiSTR version>
# command=combine
# catalog=<preset name or catalog file name>
# catalog_sha256=<SHA-256 of the decompressed catalog>
# genome_sha256=<SHA-256 of the chromosome length map>
```

`combine` refuses to merge inputs whose `catalog_sha256` or `genome_sha256` disagree — doing so
would align rows that do not describe the same loci — and carries the agreed values forward.
Inputs lacking the fields are warned about rather than rejected.

### Columns

```
chromosome  begin  end  info  <sample1>_H1  <sample1>_H2  <sample2>_H1  <sample2>_H2  ...
```

The first four columns are the same as in the individual call file. Each subsequent pair of columns (`_H1`, `_H2`) contains the allele lengths for one sample. Column order matches the order in which samples were provided to `inquiSTR combine`.
