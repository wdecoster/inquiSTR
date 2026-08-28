# Porting the PLVI framework to inquiSTR

Design specification.

**Status:** all four steps are implemented — `varindex`, the outlier framework, `annotate` and
`prioritize` — along with the shared `locus_stats` reader. The pipeline runs end to end and has
been validated on a 434-genome ONT cohort.

## Scope

Danzi et al. (bioRxiv 2025.01.06.631535v3) define the Pure Length Variability Index (PLVI)
and pair it with a per-genome length-outlier scan to prioritise candidate repeat expansion
loci. This document specifies how that framework is adapted to inquiSTR.

The adaptation is not a straight port. inquiSTR genotypes **total allele length** relative to
the reference; PLVI is derived from the **longest pure segment (LPS)**, which requires
sequence-resolved alleles that inquiSTR does not produce. These are different axes. Anything
computed internally from length is named as a length quantity, and LPS-derived values imported
from the publication are never substituted for it.

Target scale: cohorts of ~500 individuals (1,000 alleles), catalogs up to 4.4M loci.

**Scope boundary.** This pipeline finds rare, high-effect expansions where the carrier is
extreme within the cohort. Common risk alleles carried by many people are the job of
`inquiSTR association`, and users with phenotype data should start there — an outlier scan
reports only the most extreme allele at a locus, so multiple carriers mask each other. See
*Start with association, not with outliers* in `docs/PRIORITIZATION.md`.

## Background: what the publication defines

- **PLV** — per-locus standard deviation of LPS length, using only the locus's most common
  LPS motif.
- **PLVI** — loci grouped by motif length, then ranked within each group by PLV.
- **Outlier** — an allele longer than the 99.9th percentile of a reference cohort. The 99.9th
  percentile was chosen to tolerate late onset and incomplete penetrance.
- **Funnel** — a median of 3,654 length outliers per genome, 2,251 when measured on LPS, and
  155 after restricting to loci in the top 5% of PLVI.

The number that matters most for this port is in Figure 3B. Measured as generalized odds ratio
for enrichment of the 42 known RED loci:

| measure | OR |
| --- | --- |
| SD of total allele length | 34.1 |
| SD of LPS length (PLV) | 68.6 |
| PLV ranked within motif-length groups (PLVI) | 243.9 |

A length-only index is therefore a real effect and roughly an order of magnitude weaker than
full PLVI. 85.7% (36/42) of RED loci fall in the top 5% by allele-length SD alone. That is the
honest framing for user-facing documentation: worth computing, not equivalent to PLVI.

## Evidence gathered (2026-08-27)

Findings from inspecting the Zenodo release and the manuscript's own analysis code. These
constrain the design and are recorded so the reasoning survives.

### The Zenodo release contains percentile tables, not distributions

DOI 10.5281/zenodo.19239224 currently resolves to record 19895393. Every locus-level table
shares one schema:

```text
TRID  N  0th 1st 5th 10th 15th 20th 25th 30th 35th 40th 45th 50th 55th 60th 65th 70th
75th 80th 85th 90th 95th 99th 99.9th 100thPercentile  mean  std  mad  mode  uniqueLengths
```

plus `uniqueAlleleSeqs` in the `alleleLengthStats` files and `CPS` (composition polymorphism
score) in the `lpsPerLocus` files.

There are no per-length counts and no histograms. Consequences:

- There is no external raw tail, so a generalized Pareto fit has nothing external to fit
  against. It could only ever fit the user's own thin 1,000-allele tail — the exact weakness it
  was proposed to solve. Dropped; see *Out of scope*.
- No memory-mapped or columnar store is needed. Only five columns are ever read
  (`N`, `99thPercentile`, `99.9thPercentile`, `100thPercentile`, `std`); 4.4M rows of those
  parse into a modest in-memory table.

The cohorts available, all GRCh38:

| file prefix | cohort | technology | caller | catalogs |
| --- | --- | --- | --- | --- |
| `DiscoveryCohort_` | 543 genomes, median 32x | HiFi | TRGT v5.0.0 | TR-Explorer v1.0.1, Adotto |
| `ValidationCohort_` | 2,102 genomes, median 12x | HiFi | TRGT v5.0.0 | both, plus AFR/AMR/EAS/EUR/SAS subsets |
| `ReplicationCohort_` | 500 genomes, median 35x | **ONT** (1KGP) | Medaka Tandem | TR-Explorer v1.0.1 |

`ReplicationCohort_alleleLengthStats.TRExplorer_1.0.1.txt.gz` is both technology-matched and
axis-matched for an ONT length-based tool, and is the recommended default external reference.

### The external tables give an index on inquiSTR's own axis

`alleleLengthStats` reports `std` of **total allele length** per locus. That is the same
quantity as the internal index, computed on a population cohort. So three distinct columns
exist, and each is named for its axis and its source:

| column | axis | source |
| --- | --- | --- |
| `lvi_length_internal` | total allele length | the user's own cohort |
| `lvi_length_external` | total allele length | Zenodo `alleleLengthStats` |
| `plvi_lps_external` | LPS length | Zenodo `lpsPerLocus` (`std` is PLV by definition) |

The first two are directly comparable, which is what makes the internal/external disagreement
check meaningful. The third is an annotation on a different axis and must never be averaged
with or substituted for either of the others.

### TRID is coordinate-derived, so the join is exact

TRIDs take the form `1-160919121-160919141-A`, that is
`chromosome (no chr prefix) - start - end - motif`. The coordinates match the TR-Explorer BED
intervals directly: `end - start` equals the modal allele length in bp on sampled rows
(`1-160919121-160919141-A` has mode 20, interval width 20; `1-12205582-12205592-GCA` has mode
10, interval width 10).

For a run against the `trexplorer` preset the join is therefore an exact key match. Reciprocal
overlap and motif-compatibility checking are only required when joining to a different catalog.

### Unit mismatch between inquiSTR and the reference tables

inquiSTR reports **length change relative to the reference**; the Zenodo tables report
**absolute length**. The bridge is:

```text
absolute_length = (end - begin) + inquiSTR_value
```

verified against the sampled rows above.

Note also that inquiSTR genotypes are not integers. `median_str_length` returns `f64`, and an
even number of supporting reads yields half-integer medians. Keeping allele lengths as `u32`
end to end is therefore not achievable at the inquiSTR end. Exceedance counting works on
floats; the rounding rule applied before comparison against an integer reference threshold must
be documented where it happens.

### The 99.9th percentile threshold is saturated at this cohort size

In `ReplicationCohort_alleleLengthStats` (N = 998 alleles), across the first 10,487 loci:

- **79.7%** of loci have `99.9thPercentile == 100thPercentile`
- **78.0%** have `99thPercentile == 100thPercentile`

At four loci in five, "above the 99.9th percentile" therefore means "longer than anything ever
observed". The same applies to an internal 500-sample reference, which is the same 1,000-allele
resolution.

This is why the exceedance method reports `k` and `N` on every call rather than a nominal
percentile: without them the threshold silently means different things at different loci. It is
also why an empirical exceedance rank is preferred over a p-value that cannot be resolved.

### How the publication actually counted outliers

From `scripts/polars_streaming_outlier_count.py` in the manuscript repository:

- Allele 1 and allele 2 are tested **independently** and the counts **summed**.
- Comparison is strict `>` against the reference percentile column.
- No-calls (`.`) are coerced to `0` and therefore never count as outliers.
- A `(TRID, LPS_Motif)` combination absent from the reference counts as an outlier if length
  > 0. This branch has no length-axis analogue and does not apply here.

The published funnel figures are consequently per-allele sums with dropout excluded by
construction.

### Validated against a real cohort file

Measured on a 434-sample ONT combined file called against the TR-Explorer catalog
(`file_type=combined_call`, 4,863,041 loci, 872 columns, 910 MB gzipped and 9.5 GB
uncompressed — about 1,950 bytes per locus), sampling the first 400,000 loci:

| observation | value | consequence |
| --- | --- | --- |
| `info` parses as a clean motif | 400,000 / 400,000 | motif-length grouping works straight from the combined file for this catalog |
| mean reference locus length | 13.5 bp | most loci are short; a large positive delta is a large *relative* change |
| no-call rate | 2.34% | call rate is high, but dropout is not negligible at ~1 cell in 43 |
| negative (contraction) genotypes | 10.49% | one genotype in ten is a contraction |
| largest observed delta | +93,647 bp | on a mean 13.5 bp locus — certainly artifactual |

Two of these change the design rather than merely confirming it.

**Contractions are 10% of the data.** The existing `streaming_stats` restriction to values
`> 0.0` therefore discards about a tenth of all genotypes from every reference distribution it
builds, on top of conflating no-calls with reference-length alleles. Under absolute lengths the
restriction becomes meaningless, which is the right outcome, but the size of the effect is worth
knowing when old and new results are compared.

**Extreme deltas are real in the data and are not real biology.** A +93 kb call at a 13.5 bp
locus is a soft-clip artifact, not a repeat. Since `prioritize` ranks by score across loci,
artifacts of this magnitude will otherwise occupy the top of every shortlist. An upper
plausibility bound — absolute, or relative to the locus reference length — is required, and
should be reported in the funnel as its own filtered category rather than silently dropped.

## Design decisions

1. **Motif resolution is catalog-aware for presets, with a fallback for custom catalogs.**
   Motif length is needed to group loci for ranking. inquiSTR keeps only BED column 4, as an
   opaque `info` string, and that is the motif for `trexplorer`, a comma-separated motif list
   for `pathogenic` (STRchive), and a bare integer for `adotto` — whose motif sits in a JSON
   blob in column 18 that inquiSTR never reads.

   `combine` preserves `info` from its inputs (`src/combine.rs:354`), so for a `trexplorer` run
   the motif is already present in the combined file and no catalog re-read is needed. Only
   `adotto` and custom catalogs require joining back to the original BED by coordinate. Rules:

   | catalog | rule |
   |---|---|
   | `trexplorer` | column 4 verbatim |
   | `pathogenic` (STRchive) | first of the comma-separated motifs; warn if the listed motifs differ in length |
   | `adotto` | parse the JSON in column 18, take the motif of the highest-scoring span |
   | `codis` | column 4 verbatim |
   | custom BED | column 4 if it parses as a motif; otherwise a single ungrouped ranking with a prominent warning |

   Only motif *length* is needed, not motif identity, which makes the STRchive multi-motif case
   benign in practice.

   Verified on a real 434-sample TR-Explorer combined file: `info` parsed as a clean motif in
   400,000 of 400,000 sampled loci, so for this catalog the motif is available directly from the
   combined file with no catalog re-read.

2. **Motif-length groups are binned to a minimum size, and no locus is excluded.** Ranking
   within a group of one is not a ranking: that locus is rank 1 of 1, sits in the top percentile
   by construction, and passes any "top 5%" filter automatically — a false-positive generator
   placed exactly where the tool claims to be selective, and one that would recur in every
   sample's shortlist rather than appearing at random.

   Full counts from the 434-sample file (all 4,863,041 loci, not a sample):

   | motif length | loci | share |
   | --- | --- | --- |
   | 1 bp | 1,567,337 | 32.2% |
   | 2 bp | 978,972 | 20.1% |
   | 3 bp | 1,432,117 | 29.4% |
   | 4 bp | 590,787 | 12.1% |
   | 5 bp | 177,422 | 3.6% |
   | 6 bp | 56,731 | 1.2% |
   | 7 bp and above | 59,675 | 1.2% |

   There are **208 distinct motif lengths**, the longest motif being 833 bp. The tail is thin but
   not negligible, and it is almost entirely unrankable as-is: 111 groups hold fewer than 10 loci
   each (254 loci in total), 155 groups hold fewer than 100 (1,694 loci), and 186 groups hold
   fewer than 1,000 (16,058 loci). Only 22 motif lengths have 1,000 or more loci.

   **Bin adjacent motif lengths until each bin reaches a floor**, rather than excluding anything.
   Grouping by exact motif length is only a proxy for comparing like with like; a contiguous range
   of motif lengths serves that purpose equally well and keeps every locus rankable. A greedy
   merge over ascending motif length with a floor of 1,000 loci yields 30 bins on this catalog:

   ```text
   1  2  3  4  5  6  7  8  9  10  11  12  13-14  15  16  17  18  19-20  21-22
   23-24  25-26  27-28  29-30  31-32  33-34  35-36  37-39  40-44  45-51  52-833
   ```

   Motif lengths 1 through 12 remain their own groups because the data supports it; only the
   sparse tail is merged. The floor is configurable and the binning is derived from the catalog
   at runtime, so no motif-length constants are hard-coded and other catalogs bin themselves.
   Report the realised bin edges in the output metadata — a rank means nothing without knowing
   what it was ranked against.

   **On excluding VNTRs.** An earlier draft of this document restricted the index to motif
   lengths 1–6 bp, citing the publication's statement that *"we excluded variable number tandem
   repeats (VNTRs) from our analysis. VNTRs are tandem repeats with motif length over 6 base
   pairs."* That reading does not hold up. The same paper reports, of Figure 2E, *"a larger degree
   of LPS length variation in longer tandem repeat motifs (12+bp) compared to shorter motifs
   (generalized odds ratio: 26.35, p<1e-10)"* — an analysis that requires motifs well past 6 bp.
   The published Zenodo tables confirm it: motif lengths of 7, 8, 12, 24, 51 and 62 bp all appear
   within the first ~10,000 loci sampled. The VNTR sentence most plausibly refers to not adopting
   a dedicated VNTR catalog, not to filtering the motif column. Binning therefore replaces
   exclusion, and the design takes no position on where STRs end and VNTRs begin.

   Homopolymers deserve their own note: 1 bp is the single largest group at 32% of the catalog,
   and ONT length estimates are least reliable there. They stay in scope, but the documentation
   should say plainly that a high-ranking homopolymer locus is more likely to be a sequencing
   artifact than a candidate.

3. **A per-locus summary sidecar is a first-class artifact, in the Zenodo schema, in absolute
   length.** A shared schema is pointless if the units differ, so the emitter applies the
   `(end - begin) + value` conversion on the way out and records it in the metadata header. An
   imported Zenodo table and a locally computed one are then the same kind of file, read by one
   code path. This also avoids repeated scans of a combined file that, at 500 samples by 4.4M
   loci, is around 10 GB of TSV (measured: 9.5 GB for 434 samples × 4.86M loci).

4. **Zygosity: both alleles are tested, results collapse per sample.** A sample is flagged at a
   locus if either allele exceeds the threshold. This departs from the publication's per-allele
   sum, so the funnel diagnostic reports **both** the per-sample count and the per-allele sum,
   keeping the published figures usable as a sanity check. Which allele triggered the call is
   carried in the evidence field. `min` and `both` remain selectable for biallelic disorders
   such as RFC1 and FXN, where defaulting to the longer allele would silently lose the signal.

5. **The external reference ships as a downloadable, cached preset.** This reuses the existing
   `TRPreset` machinery in `repeats.rs`: XDG cache directory, lock file, atomic replace. The
   Zenodo record has already moved once (19239224 → 19895393), so the preset resolves through
   the DOI rather than a pinned record ID, and validates the expected header on load. The new
   URLs are added to the existing weekly `preset-urls.yml` workflow.

6. **Provenance metadata: catalog identity is recorded by content, not by name.** A catalog name
   alone is not enough. Catalogs are updated in place at the same URL, and a user can pass a
   locally edited BED, so `# catalog=trexplorer` proves nothing about which loci were actually
   genotyped. `call` gains three metadata fields, written in the block at `src/call.rs:153`:

   | field | value | answers |
   |---|---|---|
   | `catalog` | preset name, or the catalog file's basename | which catalog *version* this is, for joining external tables |
   | `catalog_sha256` | SHA-256 of the **decompressed** catalog bytes | are two inquiSTR files talking about the same loci |
   | `genome_sha256` | SHA-256 over the sorted `(chromosome, length)` pairs from the BAM header | which reference build |

   Hashing decompressed rather than compressed bytes matters: `adotto` and `trexplorer` ship
   gzipped, and gzip output is not reproducible across compressors, so hashing the `.gz` would
   make an identical catalog look different after any recompression. Decompressed content is
   also user-verifiable with standard tools — `zcat catalog.bed.gz | sha256sum` — which is
   essential when the hash is what causes a refusal and the user has to diagnose it.

   The genome fingerprint is free: `get_chrom_lengths_from_bam_header` already loads the
   chromosome length map in `get_targets`. GRCh38, T2T and hg19 differ in chromosome lengths, so
   this identifies the build robustly without needing the FASTA, and is what enforces the
   build-mismatch refusal.

   The two hashes answer different questions and are enforced differently:

   - **Comparing or combining inquiSTR files** needs the same target set. `combine` verifies
     that all inputs agree on `catalog_sha256` and `genome_sha256`, refuses on mismatch, and
     propagates both into the combined header.
   - **Joining an external table** needs the same catalog *version*, since the Zenodo tables
     ship no hash of the BED. That check is on the `catalog` name plus the build, and is the
     refusal described in the annotation step.

   Backward compatibility: files written before this change have none of these fields. A missing
   field warns and proceeds — it can never hard-fail, or every existing `.inq` becomes
   unreadable. Only a *present and mismatched* field refuses.

   **Cost.** Measured on a machine with SHA-NI: hashing the decompressed catalog runs at about
   1.8 GB/s, adding 6–9% to the decompress-and-parse pass the catalog already requires (332 ms
   for TR-Explorer, 910 ms for Adotto as a standalone pass). That cost would otherwise be paid
   once per sample, so the digest is memoised on disk under the cache directory, keyed by path,
   size and modification time, and written atomically. A warm lookup costs **14–31 µs**, four
   orders of magnitude less. On disk rather than in process deliberately: a workflow manager
   runs each `call` as its own process, so an in-memory cache would never hit.

7. **All `outlier` output moves to the new file type.** The current output — `chrom begin end`
   plus a comma-joined sample list, untyped, no metadata header — cannot express a score, and
   ranking is the entire point of the framework. Every method, including the existing z-score
   and DBSCAN, emits the new typed schema. This is a breaking change to anything parsing the
   current output, and is called out in the release notes.

## Existing behaviour to be aware of

The current z-score path in `src/outlier.rs` has properties the new code must not inherit:

- `get_repeat_lengths` (`src/outlier.rs:491`) converts `NaN` and empty fields to `0.0`, making
  a no-call indistinguishable from a genuine reference-length allele.
- `streaming_stats` (`src/outlier.rs:31`) accumulates only values `> 0.0`, so the reference
  mean and SD are computed from expanded alleles alone — contractions and reference-length
  alleles are excluded.
- `z_score_outliers` (`src/outlier.rs:564`) then scores **every** entry, zeros included,
  against that positives-only mean and SD.

Under the absolute-length framing this class of problem disappears, since absolute lengths are
always positive and a no-call is representable as missing. The discrepancy should be expected
when comparing new results against the existing z-score method, and is worth a separate issue.

`--min-cv` in `src/filter.rs` already offers a per-locus coefficient of variation. It is
related to but distinct from the variability index: CV normalises by the mean, the index ranks
raw SD within motif-length groups. Both should not be described as the same thing.

## How this fits the existing pipeline

Today, with `batch` being `call` + `combine` fused over a manifest:

```text
BAM/CRAM ──► call ──► sample.inq ──► combine ──► combined.inq ──► outlier ──► stdout
                  individual_call             combined_call     z-score / DBSCAN
```

After. Two axes meet at the end: `outlier` produces **allele-level** evidence (which sample is
extreme where), `varindex` and `annotate` produce **locus-level** evidence (which loci are worth
caring about), and `prioritize` is the join that turns both into a per-sample ranked shortlist.

```text
BAM/CRAM ─► call ─► sample.inq ─► combine ─► combined.inq   (~10 GB at 500 × 4.9M)
                                                   │
                            ┌──────────────────────┴──────────────────────┐
                            │  locus axis                    allele axis  │
                            ▼                                             ▼
                    varindex [--catalog]                              outlier
                            │                                             │
                            ▼                                             ▼
                      cohort.stats                                 outlier_calls
                   (file_type=locus_stats)                    (sample × locus, scored)
                            │                                             │
                            ├──► annotate --plvi ──► mapping.tsv ──┐      │
                            │                                      │      │
                            └──────────────┬───────────────────────┴──────┘
                                           ▼
                                       prioritize
                                           │
                                           ▼
                              per-sample ranked shortlist
                                    + funnel report

Zenodo table ──preset download + cache──► external.stats
                (TRPreset machinery)      same locus_stats schema
                                              │
                                              ├──► outlier      (external reference)
                                              └──► prioritize   (external index)
```

The one arrow that is easy to misread: `cohort.stats` is **not** an input to `outlier` in the
default internal mode. Internal leave-one-out reads alleles straight from the combined-file
line, so the two branches run independently and can run concurrently. A `locus_stats` table
feeds `outlier` only when it is serving as an *external* reference — a different cohort, or a
Zenodo import.

Per command:

- **`call`** — genotyping and the `.inq` schema are unchanged. None of this implies re-running
  BAMs. The only addition is the three provenance fields from decision 6.
- **`combine`** — already carries `info` through, which is what makes motifs free for
  `trexplorer`. Gains verification and propagation of the provenance fields.
- **`outlier`** — refactored in place, not replaced. The streaming skeleton stays:
  `outlier_str_analysis` (`src/outlier.rs:187`) chunking at 1,000 lines, `process_chunk`
  (`:351`) fanning out over rayon. `ProcessParams` (`:341`) grows a reference source; the
  `match params.method` becomes the enum dispatch; `get_repeat_lengths` (`:491`) stops mapping
  `NaN` to `0.0`; output becomes the typed scored table. It gains no knowledge of the index —
  that stays downstream, per decision C1's rule against silent lookups.
- **`varindex`**, **`annotate`** and **`prioritize`** are new. The first two read the combined
  file; `prioritize` reads only the small derived tables.

The combined file remains the hub — `relate`, `pca`, `assoc`, `filter` and `query` are
unaffected. The practical payoff of the shared `locus_stats` artifact is that the ~10 GB file is
scanned twice at most, once per axis, and everything after that reads a table an order of
magnitude smaller (measured: 744 MB uncompressed, 210 MB gzipped, from a 9.5 GB combined file).
Re-ranking at a different index cutoff costs seconds, not a rescan.

Incidental cleanups the work passes through:

- `get_sample_names_from_file` (`src/outlier.rs:56`) skips 3 columns where
  `STR_FIXED_COLUMNS` is 4, so `info` is counted as a sample. Cosmetic today — it only feeds
  `--sample` validation and the "available samples" listing — but this function is touched.
- Sample-name cleaning exists three times: `outlier.rs:16`, `filter.rs:213`, and
  `locus_search.rs:228` (`extract_clean_sample_names`). The third is the one to keep.

## Build plan

### 1. `varindex` — per-locus summary and variability index *(implemented)*

```text
inquiSTR varindex combined.inq [--catalog <bed|preset>] -o cohort.stats
```

One command, one pass over the combined file, producing one artifact. It emits per-locus `N`,
the percentile grid, `mean`, `std`, `mad`, `mode` and `uniqueLengths` in the Zenodo column
layout and in absolute length, and — in the same pass — resolves motifs and emits the ranked
`lvi_length_internal` column. Splitting the summary from the index would mean two scans of a
~10 GB file to compute quantities that share the same running sums.

**Naming.** The command is `varindex` because the variability index is the reason a user runs
it. The *file type* is the more general `locus_stats`, because the identical schema arrives from
two directions: a `varindex` run, and an imported Zenodo table that has percentiles but no index
column. The index column is therefore optional in the schema, and its presence is flagged in the
metadata header so readers do not have to guess. One reader serves both.

`file_type=locus_stats` is registered in `src/filetype.rs`, with the metadata header recording
the source combined file, the catalog provenance fields from decision 6, the length conversion,
and whether the index column is present.

This is where the leave-one-out machinery lives. Leave-one-out standard deviation is computed
closed-form from running sums, never by recomputation per sample:

```text
n' = n - m                      m = the sample's non-missing alleles at this locus
S' = S - sum(sample alleles)
Q' = Q - sum(sample alleles squared)
mean' = S' / n'
var'  = max(0, Q'/n' - mean'^2)
```

A naive implementation performs N recomputations per locus and will dominate runtime on a 4.4M
locus catalog. Per-locus work is embarrassingly parallel; rayon over loci is the natural
structure.

Motif resolution follows decisions 1 and 2: `info` from the combined file where it carries the motif,
the original catalog by coordinate join otherwise, `--catalog` required only in that second
case. Ranking is within motif-length groups; where no motif is resolvable the ranking is
ungrouped and says so loudly.

Documentation must carry two warnings prominently:

- **Ascertainment.** The published reference is a population sample. A user's 500 may be
  phenotype-ascertained, in which case a causal locus's variability is inflated by the very
  signal being prioritised, and internal ranking partially cancels it. External annotation is
  recommended for phenotype-enriched cohorts.
- **Cohort size is not the concern.** 500 individuals is 1,000 alleles, the same order as the
  reference cohorts: the observed per-locus `N` is 1,086 in the discovery table and 998 in the
  ONT replication table. An SD converges quickly at that scale. What ~1,000 alleles supports
  poorly is the *tail*, which is the outlier test, not the index — and the saturation figures
  above quantify how poorly.

### 2. Outlier framework *(implemented, except the external reference mode)*

An `OutlierMethod` trait, dispatched through an enum so the per-locus hot loop stays
monomorphised, returning:

```rust
(score, significance, evidence)
```

- `score` — monotone in how extreme, comparable **across loci**, independent of motif length.
  Work in repeat units, or normalise bp by motif length. Cross-locus comparability is the
  point: a genome yields thousands of nominal outliers and they need ranking, not
  classification.
- `significance` — p-value, exceedance probability, or empirical rank, per method.
- `evidence` — a typed enum with per-method variants, not a stringly-typed map, so diagnostics
  stay typed through to the output writer. Carries reference N, effective rank, fitted
  parameters, cluster assignment, and which allele triggered the call.

`OutlierResult` has no boolean field, which enforces by construction that thresholding is a
downstream, user-configurable step.

The resulting method set on `--method`:

| method | status |
| --- | --- |
| `zscore` | existing, moved behind the interface unchanged in behaviour |
| `dbscan` | existing, moved behind the interface — but see the note below |
| `percentile` | **new, and the new default** |
| `robustz` | optional fourth, cheap to add alongside `percentile` |

Changing the default from `zscore` to `percentile` is a behaviour change for existing
invocations that do not pass `--method`, and belongs in the release notes next to the output
schema change.

All four emit the typed output described in decision 7 — no method returns a bare sample list.

`dbscan`'s neighbourhood radius, `max(2 x mode, 10)`, is unchanged in form but is now computed
on **absolute** lengths rather than reference-relative deltas. The mode was previously almost
always zero, making the radius a constant 10; it now scales with locus length. The formula was
carried over, its behaviour was not — worth revisiting rather than assuming it still holds.

- **Exceedance / percentile method**, the recommended default. Alleles are pooled (2N per
  locus, minus no-calls), never samples. The threshold is defined by **exceedance count** —
  "longer than the k-th largest reference allele" — not by an interpolated quantile, since
  lengths are heavily tied and interpolation produces thresholds corresponding to no observable
  value. `k` and `N` are emitted alongside every call so a degenerate max-based rule is visible.
  One-sided (expansion) by default, two-sided opt-in for contractions.
- **Robust z (median/MAD)** as a cheap improvement on the current z-score. Both carry an
  explicit documented caveat: a z-score divides by the locus standard deviation, which is the
  same quantity the variability index rewards, so sensitivity degrades at exactly the
  high-variability loci the index prioritises — the two filters pull against each other.
  Rank- and tail-based methods are scale-free within locus and do not have this property.
  Repeat-length distributions are also heavy-tailed and frequently multimodal, so the SD in the
  denominator is partly inflated by the tail being tested. The same tension appears in DBSCAN
  via `eps`: an absolute `eps` fragments high-variability loci into multiple clusters, while an
  `eps` scaled to locus spread reimports the z-score behaviour. Worth checking empirically
  whether the noise-point rate correlates with the variability index across loci; if it does,
  the two signals are not independent.

Reference distribution modes:

- **Internal** (default) — the cohort in the current run, leave-one-out, so a sample never
  contributes to the null it is tested against. This mode needs **no precomputed table**: the
  combined-file line already holds every allele at the locus, so exclusion and ranking happen in
  place, exactly where `process_chunk` already operates. That matters most for the exceedance
  method, where an expansion carrier holding the longest allele would otherwise define the very
  threshold they are tested against and score as unremarkable. Must also support a relatedness
  exclusion list, since a sib pair sharing an expanded allele mutually suppresses each other's
  significance.
- **External** — *not implemented.* A `locus_stats` table (another cohort's `varindex` output,
  or an imported published table) would supply the tail from its percentile columns, with no
  leave-one-out applied. `ScanConfig` has no field for it and there is no CLI flag; the scan
  writes `# reference=internal_leave_one_out` unconditionally. This is the largest remaining
  gap in the design, and the one that matters most for phenotype-ascertained cohorts.

Note that a `locus_stats` table is sufficient for closed-form moment-based statistics even
though it stores no sums — `S = N·mean` and `Q = N·(std² + mean²)` recover them exactly. It is
*not* sufficient for exact tail ranking, which is why internal mode reads alleles directly
rather than round-tripping through `varindex`.

**How much does `outlier` emit?** Decision A1 says thresholding belongs downstream, but scoring
every locus in every sample is 4.4M × 500 ≈ 2.2 billion rows, which is not a file anyone wants.
The resolution: `outlier` applies a deliberately *permissive* gate — by default, alleles
exceeding the reference maximum — and emits a scored row for everything that passes, with no
decision thresholds applied. At the published rate of ~3,654 outliers per genome that is roughly
1.8M rows for a 500-sample cohort, a file of a few hundred MB. Decision thresholds, index
filtering and ranking all live in `prioritize`, which can then be re-run at many settings without
touching the combined file again.

### 3. Annotation join *(implemented)*

```text
inquistr annotate --plvi <table> --catalog <bed> --out <mapping>
```

A separate, inspectable step producing an auditable mapping file rather than a silent lookup
inside the outlier command. Requirements:

- Report the match rate up front. A user whose catalog largely fails to map needs to know
  before reading any rankings.
- Exact TRID key match for TR-Explorer v1.0.1. For other catalogs, reciprocal overlap at a
  configurable threshold **plus** a motif-compatibility check: the annotation's motif should be
  rotation- or reverse-complement-equivalent to the catalog's.
- Emit a per-locus join quality: exact / high-confidence / ambiguous / unmapped, where
  ambiguous covers many-to-one and one-to-many overlaps, i.e. the compound-locus case.
- **Never resolve an ambiguous many-to-one match by taking the maximum PLVI.** That
  systematically inflates exactly the messy loci.
- Reference build and catalog version mismatches are **refusals, not warnings**. The published
  tables are GRCh38 / TR-Explorer v1.0.1; joining to a T2T or hg19 catalog fails loudly.
- Report both internal and external indices where available, and flag large rank
  discrepancies rather than averaging or silently preferring one. A locus ranking far higher
  internally is either an ancestry-specific difference or a cohort-composition artifact, and
  both merit attention. *(Partly implemented: `prioritize --annotation` reports the external
  index in its own column, but ranking and filtering use the internal index only, and no
  discrepancy flag is computed.)*

PLV is boundary-sensitive: LPS depends on where a locus is cut and whether adjacent repeats are
merged, which is what TR-Explorer's variation-cluster system arbitrates. The published Adotto
replication shows the metric is robust when *recomputed* on another catalog, not that values
transfer by coordinate overlap. The annotation is a convenience, not a substitute for
recomputation.

### 4. `prioritize` — the join, and the funnel *(implemented)*

```text
inquiSTR prioritize outlier_calls --index cohort.stats [--annotation mapping.tsv]
                                  [--top-index-pct 5] -o shortlist.tsv
```

This is the command that answers the question a user actually has: *for this genome, which
repeats should I look at first?* It is the operation behind Figure 4C → 4D, where filtering an
unranked outlier list by locus index moves the pathogenic HTT allele from rank 143 to rank 2.

Neither input can do this alone. `outlier` knows which alleles are extreme but nothing about
which loci matter; `varindex` knows which loci are unstable but nothing about any individual.
`prioritize` joins the allele axis to the locus axis and ranks the result per sample:

1. Join `outlier_calls` to the locus-level index by coordinate.
2. Apply the locus filter — by default the top 5% by index, matching the published threshold.
3. Apply QC (below), dropping or flagging as appropriate.
4. Rank within each sample by the outlier score, which is why the score must be comparable
   across loci rather than merely within one.
5. Emit a per-sample ranked shortlist, plus the funnel report.

It reads only the small derived tables, never the combined file, so re-running it at a different
cutoff or with an annotation obtained later is cheap. That is the reason the pipeline splits here
rather than folding the index into `outlier`.

Where both an internal and an external index are available, report both and flag large rank
discrepancies rather than averaging or silently preferring one, per decision B4.

**Funnel report.** Reproduce the published filtering funnel as a standard part of the output:
outliers per sample before filtering, after the locus-index filter, and in the final ranked
list, reported both per-sample and as a per-allele sum. Published reference figures for
orientation: ~3,654 length outliers per genome, 2,251 on LPS, 155 after restricting to the top
5% of PLVI. This is probably the single most useful diagnostic in the tool — it tells a user
immediately whether their thresholds are sane.

QC will dominate the false-positive rate more than the choice of statistic:

- Per-locus reference call rate, with refusal to test below a threshold. A locus where controls
  frequently fail to genotype has an artificially truncated tail.
- Per-sample outlier count as a QC metric. A sample producing many times the cohort median is a
  coverage or alignment failure, not a patient.
- **Genotype failure at a high-index locus is its own flagged event**, not missing data. In
  length genotyping, dropout of a very long allele is a plausible expansion signature and is
  invisible to any method that only tests successfully called lengths.

## Out of scope

- **Generalized Pareto tail fitting.** No external raw tail exists to fit (see *Evidence*), and
  Rust has no maximum-likelihood GPD estimator off the shelf — `statrs` provides distributions
  but not the fit. Implementing it means Grimshaw's method or a bounded solve on the profile
  likelihood plus convergence guards for the near-zero-shape case: roughly a day of work and a
  numerical-correctness surface, not a dependency line. The exceedance method is sufficient on
  its own. Revisit only if tail saturation is demonstrably costing discriminative power on real
  data.
- **Memory-mapped or columnar reference store.** Unnecessary once the reference is five columns
  of percentile data rather than full histograms.

## Validation

- **Known outlier alleles already identified in the local cohort are the primary test.** These
  are true positives on the exact technology, catalog and pipeline the tool will be used with,
  which makes them stronger evidence than any simulation. Each method should recover them, and
  the rank they receive in the `prioritize` shortlist is the metric to track across changes —
  recovering an allele but burying it at rank 3,000 is a failure. Record the loci and carrier
  samples here once agreed, and pin them as a regression fixture. The published benchmark for
  this is a median rank of 2 across six known-RED genomes, with all six in the top 10.
- Spike-in of long alleles is the easy test and exercises the least interesting failure mode.
- **Simulate allele dropout at expanded loci.** This is the failure that will actually bite,
  and it connects directly to the requirement that genotype failure be a flagged event.
- Report sensitivity as a function of coverage. ONT users run a wide range of depths, and the
  published cohorts span ~8x to ~32x with materially different behaviour in the upper length
  percentiles.
- Regression-test the funnel counts so threshold or method changes surface as visible shifts.

## References

Danzi MC, Xu IRL, et al. *Population-scale variability at short tandem repeat loci reveals
pathogenicity signature.* bioRxiv 2025.01.06.631535v3, posted 5 June 2026.

- Analysis code: <https://github.com/ZuchnerLab/AoU_TRs_PLVI_Manuscript>
- Locus-level summary tables: <https://doi.org/10.5281/zenodo.19239224>
- TR-Explorer catalog: Weisburd B, Dolzhenko E, et al., bioRxiv 2024.10.04.615514
