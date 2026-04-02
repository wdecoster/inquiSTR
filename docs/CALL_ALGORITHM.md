# inquiSTR call — Algorithm

This document describes how `inquiSTR call` genotypes tandem repeat (STR/TR) loci from a BAM or CRAM file.

## Overview

For each target locus, inquiSTR measures the net insertion/deletion length across the locus by inspecting CIGAR strings of overlapping reads, then computes a per-haplotype median to produce the final genotype call.

## Step 1: Target batching

STR target intervals (from `--region`, `--region-file`, or `--preset`) are grouped into *batches*. Targets on the same chromosome within `--batch-size` kilobases of each other are placed in the same batch. Batching reduces the number of indexed BAM fetches and is the primary lever for I/O performance tuning.

Each batch is padded by 10 bp on both ends to ensure reads that overlap locus boundaries are captured.

## Step 2: Read filtering

For each batch, all reads overlapping the batch region are retrieved from the BAM/CRAM file. Reads with MAPQ ≤ 10 are discarded immediately.

## Step 3: CIGAR-based length measurement

For each read–locus pair, inquiSTR walks the CIGAR string and sums the lengths of INDEL/soft-clip events that:

1. fall within the locus interval (± 10 bp padding, unless `--noextend` is set), and
2. exceed the minimum length threshold (`--minlen`, default 1 bp).

The signed sum is the *call value* for that read at that locus:

- Insertions (`I`) and soft-clips (`S`) add to the value (expansions relative to reference).
- Deletions (`D`) subtract from the value (contractions relative to reference).

Calls are labelled as **spanning** (no soft-clip contribution) or **clipped** (soft-clip contributed). Spanning calls are preferred because they represent reads that fully cross the locus; clipped calls are used as a fallback.

ONT accidental 2D reads are identified via the `SA` supplementary-alignment tag and their soft-clip signal is excluded to avoid double-counting.

## Step 4: Haplotype assignment

### Phased mode (default)

Reads are separated by the `HP` auxiliary tag (`HP=1` → haplotype 1, `HP=2` → haplotype 2). Reads without an `HP` tag are set aside as unphased.

On the first batch, inquiSTR scans up to 10 000 reads to confirm that HP tags are present; if none are found it warns that the data appear unphased and recommends `--unphased`.

If reads with `HP=1` and `HP=2` are both present, each haplotype is genotyped independently. If only one phase has reads (e.g. partially phased data), all reads are pooled and split by rank order as in unphased mode. If no HP tags at all are present at a locus, both haplotype values are reported as `NaN`.

### Unphased mode (`--unphased`)

All reads are sorted by call value and split at the midpoint into a lower half (H1) and an upper half (H2). With `--imbalance <FRACTION>`, the split point is adjusted so that the top `FRACTION` of reads (by call length) are assigned to H2 and the remainder to H1. For example, `--imbalance 0.2` puts the longest 20% of calls into H2 and the shortest 80% into H1. This is useful for target capture data where allele coverage may be unequal.

## Step 5: Summary statistic per haplotype

The genotype value for each haplotype is the **median** call value across its assigned reads, computed as follows:

1. If the number of reads is below `--support` (default 3), return `NaN`.
2. If enough spanning reads are available (≥ `--support`), use only spanning reads.
3. If `--require-spanning` is set and spanning reads are insufficient, return `NaN`.
4. Otherwise, supplement spanning reads with the largest soft-clip calls until `--support` is reached.
5. Return the median of the selected values.

The output values are in base pairs relative to the reference. A value of `0` means no net length change; positive values are expansions; negative values are contractions; `NaN` means the locus could not be called (insufficient support or no phased reads).

## Step 6: Output

Results are sorted by chromosomal position and written to stdout. See [Output File Formats](OUTPUT_FORMATS.md) for the column specification.

With `--vcf`, the TSV intermediate is converted to a VCF file instead (see [docs/FILE_TYPE_METADATA.md](FILE_TYPE_METADATA.md) for metadata details).
