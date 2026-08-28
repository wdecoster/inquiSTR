# File Type Metadata System

## Overview

InquiSTR now includes a metadata system to explicitly identify file types, replacing the previous heuristic-based detection. All output files now include a comment header line that specifies the file type.

## Metadata Format

Files begin with one or more metadata lines (comments starting with `#`):
```
# file_type=<type>
# version=<inquiSTR_version>
# command=<command_name>
# <parameter>=<value>
```

Metadata lines must appear before any data or column headers. The system is designed to be extensible - new metadata types can be added in the future without breaking compatibility.

### Standard Metadata Fields

All output files include these standard metadata fields:

- **file_type**: One of `individual_call`, `combined_call`, `individual_kmer`, `combined_kmer`,
  `target_kmer`, `locus_stats`, `outlier_calls`, or `prioritized_calls`
- **version**: inquiSTR version (from Cargo.toml)
- **command**: The subcommand that generated the file (`call`, `combine`, or `kmer`)

### Command-Specific Metadata

Different commands add additional metadata relevant to their operation:

**call command:**
- `sample`: Sample name
- `minlen`: Minimum STR length threshold
- `support`: Minimum read support threshold
- `unphased`: Whether unphased mode was used (true/false)

**call command (provenance):**
- `catalog`: preset name, or the catalog file's name
- `catalog_sha256`: SHA-256 of the **decompressed** catalog. Decompressed so that a recompressed
  copy of the same catalog still matches, and so the value is reproducible by hand with
  `zcat catalog.bed.gz | sha256sum`
- `genome_sha256`: SHA-256 over the sorted, `chr`-normalised chromosome length map from the BAM
  header. Chromosome lengths differ between GRCh38, T2T and hg19, so this identifies the build
  without needing the reference FASTA

These record *what* was genotyped by content rather than by name, since catalogs are updated in
place at the same URL and can be edited locally. `combine` refuses inputs whose digests conflict
and propagates the agreed values. Both digests are memoised on disk, so a per-sample workflow
computes each one once rather than once per sample.

**varindex command:**
- `length_axis`: always `absolute_bp`, recording the conversion applied to genotypes
- `percentile_method`: how percentiles were computed (`nearest_rank`, never interpolated, so
  every reported value is a length that actually occurs)
- `motif_bins`: the motif-length groups loci were ranked within - a rank is uninterpretable
  without knowing what it was ranked against
- `min_alleles`, `min_bin_size`: the thresholds applied

**outlier command:**
- `method`, `zygosity`, `reference`, `leave_one_out`: how alleles were scored
- `score_units`: always `repeat_units`
- `max_exceedance`, `min_units`: the reporting gate

**kmer command:**
- `klength`: Maximum kmer length analyzed
- `total_reads`: Total number of reads processed
- `target`: Target kmer sequence (only for --target-kmer mode)

## File Types

The following file types are defined:

1. **individual_call** - STR genotyping results for a single sample
   - Format: `chromosome\tbegin\tend\tH1\tH2`
   - Example: `chr1\t1000\t1100\t10\t12`

2. **combined_call** - Combined STR genotyping results for multiple samples
   - Format: `chromosome\tbegin\tend\tsample1_H1\tsample1_H2\tsample2_H1\tsample2_H2\t...`
   - Example: `chr1\t1000\t1100\t10\t12\t11\t13`

3. **individual_kmer** - Kmer frequency counts for a single sample
   - Format: `kmer\tcount`
   - Example: `AAAGGG\t42`

4. **combined_kmer** - Combined kmer frequency counts for multiple samples
   - Format: `kmer\tsample1\tsample2\t...`
   - Example: `AAAGGG\t42\t38\t45`

5. **target_kmer** - Target kmer counts for multiple samples (from --target-kmer mode)
   - Format: `Sample\tcount\tkmers\tsequence`
   - Example: `sample1\t42\t7\tAAGGGAAGGGAAGGG...`

6. **locus_stats** - Per-locus allele-length distribution for a cohort, written by `varindex`
   and also the format in which published population tables are read
   - Format: `chromosome\tbegin\tend\tmotif\tN\t<percentiles>\tmean\tstd\tmad\tmode\tuniqueLengths\tmotif_bin\tlvi_length_internal`
   - Lengths are **absolute base pairs**, converted from inquiSTR's reference-relative genotypes
     as `(end - begin) + value`, so that local and published tables are directly comparable
   - Published tables key rows by a compound `TRID` (`chrom-start-end-motif`) instead of
     coordinate columns; both key forms are accepted, and the variability index column is
     optional since published tables do not carry one

7. **outlier_calls** - Scored allele-level outlier calls, written by `outlier`
   - Format: `chromosome\tbegin\tend\tmotif\tsample\tevent\tallele\tlength\talleles_passing\tscore\tsignificance\tmethod\tevidence`
   - A sample is reported once per locus, on its more extreme allele; `alleles_passing` records
     how many of its alleles cleared the gate, so per-allele counts can be reconstructed
   - `event` is `outlier` or `dropout`; `score` is in repeat units so that it is comparable
     across loci; `evidence` carries method-specific diagnostics such as the reference size and
     how many reference alleles were exceeded

8. **prioritized_calls** - Per-sample ranked candidates, written by `prioritize`
   - Format: `sample\trank\tchromosome\tbegin\tend\tmotif\tevent\tallele\tlength\tscore\tsignificance\tlocus_index\tlocus_alleles\texternal_index\tsample_flag\tevidence`

## Implementation Details

### Writing Metadata

Metadata is written by output functions in:
- `src/call.rs` - Writes `individual_call` metadata
- `src/unmapped.rs` - Writes `individual_kmer` and `target_kmer` metadata
- `src/combine.rs` - Writes `combined_call` and `combined_kmer` metadata

### Reading Metadata

The function `read_file_type_metadata()` in `src/combine.rs` reads and parses metadata headers:

```rust
pub fn read_file_type_metadata(file_path: &Path) -> Option<FileType>
```

This function scans through **all** lines starting with `#` to find the `file_type` metadata, making it extensible for future metadata types.

Returns:

- `Some(FileType::*)` if metadata header is present and valid
- `None` if no metadata header is found (backward compatibility)

The helper function `skip_metadata_lines()` in `src/utils.rs` skips all metadata lines when reading files:

```rust
pub fn skip_metadata_lines(first_line: String, lines: &mut dyn Iterator<...>) -> String
```

This function keeps skipping lines while they start with `#`, ensuring all future metadata types are properly handled.

### Backward Compatibility

All file detection functions maintain backward compatibility by:
1. First attempting to read metadata using `read_file_type_metadata()`
2. If no metadata is found, falling back to heuristic detection (checking column headers, column counts, etc.)
3. Properly handling files that may start with `#` comments

This ensures that:
- Old files without metadata continue to work
- New files with metadata are detected more reliably
- Mixed environments work seamlessly

### File Detection Functions

The following functions use the metadata system:

1. **is_kmer_file()** - Determines if a file contains kmer data
   - Checks metadata first, falls back to column header inspection
   
2. **is_combined_str_file()** - Determines if a file contains combined STR data
   - Checks metadata first, falls back to column count heuristics

3. **get_samples_from_combined_str_file()** - Extracts sample names from combined STR file
   - Skips metadata line before reading column headers
   
4. **get_samples_from_combined_kmer_file()** - Extracts sample names from combined kmer file
   - Skips metadata line before reading column headers

## Example Files

### Individual STR Call
```
# file_type=individual_call
# version=0.19.0
# command=call
# sample=sample1
# minlen=10
# support=3
# unphased=false
chromosome	begin	end	H1	H2
chr1	1000	1100	10	12
chr2	2000	2200	15	15
```

### Combined STR Call
```
# file_type=combined_call
# version=0.19.0
# command=combine
chromosome	begin	end	sample1_H1	sample1_H2	sample2_H1	sample2_H2
chr1	1000	1100	10	12	11	13
chr2	2000	2200	15	15	14	16
```

### Individual Kmer
```
# file_type=individual_kmer
# version=0.19.0
# command=kmer
# klength=6
# total_reads=10000
kmer	sample1
AAAGGG	0.004200
CCCGGG	0.003800
```

### Combined Kmer
```
# file_type=combined_kmer
# version=0.19.0
# command=combine
kmer	sample1	sample2	sample3
AAAGGG	0.004200	0.003800	0.004500
CCCGGG	0.003100	0.002900	0.003300
```

### Target Kmer
```
# file_type=target_kmer
# version=0.19.0
# command=kmer
# target=AAAGGG
Sample	Target_Kmer	Canonical_Kmer	Kmer_Length	Count	Total_Reads	Frequency
sample1	AAAGGG	AAAGGG	6	42	10000	0.004200
sample2	AAAGGG	AAAGGG	6	38	9500	0.004000
```

## Command Validation

Several commands now validate that they receive the correct file type:

- **pca** - Requires `combined_call` or `combined_kmer` (rejects individual files)
- **query** - Requires `combined_call` or `combined_kmer` (rejects individual files)
- **plot** - Requires `combined_call` or `combined_kmer` (rejects individual files)
- **histogram** - Requires `combined_call` or `combined_kmer` (rejects individual files)
- **outlier** - Requires `combined_call` or `combined_kmer` (rejects individual files)

If an invalid file type is provided, these commands will display:

```
ERROR: <Command> requires a combined file (combined_call or combined_kmer).
The provided file appears to be: IndividualCall

Please use 'inquiSTR combine' to merge individual sample files first.
```

This prevents common user errors and provides clear guidance.

## Benefits

1. **Robustness** - Eliminates ambiguity in file type detection
2. **Maintainability** - Easier to add new file formats in the future
3. **Extensibility** - Multiple metadata lines supported (not just file_type)
4. **Debugging** - Files are self-documenting with explicit type information
5. **Backward Compatibility** - Old files without metadata continue to work
6. **Performance** - Metadata check is faster than heuristic inspection
7. **User Safety** - Commands validate input file types and provide helpful error messages

## Future Enhancements

The metadata system is designed to be extensible. Potential additions:

- Version information: `# version=0.19.0`
- Processing parameters: `# klength=6`
- Timestamps: `# created=2024-01-15T10:30:00Z`
- Reference genome: `# reference=hg38`
- Processing history: `# created_by=inquiSTR batch`

All functions use `skip_metadata_lines()` which skips **all** lines starting with `#`, ensuring forward compatibility with any future metadata additions.
