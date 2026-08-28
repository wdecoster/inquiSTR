# Prioritizing repeat expansion candidates

> **Status: implemented.** `varindex`, the reworked `outlier` scan, `annotate` and `prioritize`
> all run. See `docs/PLVI_PORT.md` for the implementation specification and the design
> rationale.

A genome contains millions of tandem repeat loci, and thousands of them will carry an allele
longer than anything seen in a reference cohort. Listing those outliers is not useful on its
own — the pathogenic one is somewhere in the list, unranked, among thousands of benign ones.

Prioritization narrows that list by combining two independent kinds of evidence:

- **Which alleles are unusual** in this individual, relative to a population.
- **Which loci are unstable** in the population at all, regardless of any individual.

Neither is sufficient. The second is the useful and less obvious half: repeat loci that cause
disease are, as a class, among the most length-variable loci in the general population for their
motif size. That means a population cohort with no patients in it can tell you which loci are
worth looking at — and filtering an outlier list by that property moves genuine pathogenic
expansions from the middle of a list of thousands into the top handful.

## Start with association, not with outliers

If you have case/control labels or a quantitative phenotype, **run `inquiSTR association`
first**, and come here only if it returns nothing.

The two pipelines answer different questions and are good at different things:

| | `association` | `outlier` + `prioritize` |
| --- | --- | --- |
| finds | alleles whose length tracks a phenotype across the cohort | alleles that are extreme in one individual |
| suits | common risk alleles, moderate effects, many carriers | rare high-effect expansions, often a single carrier |
| needs | phenotype data | nothing beyond the genotypes |
| evidence | a p-value against the phenotype | a rank against the cohort, weighted by locus instability |

An outlier scan reports only the most extreme allele at a locus, so a variant carried by tens
of people is largely invisible to it — each carrier hides the others. That is exactly the shape
of a common risk allele, and exactly what association testing is built for.

The *GOLGA8A* repeat is the worked example. A CT-rich expansion at chr15:34,419,425–34,419,451
is a major risk factor for atypical FTLD-U, carried by around 60% of cases; the associated
haplotype is present in 6% of controls. It was found by testing repeat length against
phenotype, and a second significant locus nearby — an intergenic repeat between *GOLGA8A* and
*GOLGA8B*, on average 8 bp longer on the risk haplotype — has no expanded alleles at all.
Neither would surface from an outlier scan, and the second would not clear a locus-instability
filter either. Both are unambiguous association findings.

So: phenotype in hand, start with association. Turn to prioritisation when the case/control
comparison comes back empty, when the cohort is too small for association to have power, or
when you are looking at a single genome with no comparison group at all.

One caveat that applies to both. inquiSTR measures **length only**. Where pathogenicity is
defined partly by sequence composition — *GOLGA8A* again, where the published criterion is
longer than 450 bp **and** more than 80% CT content — a length-based tool can identify
candidates but cannot confirm them. Resolving the motif needs a sequence-resolved genotyper.

## How to run it

The pipeline splits along those two kinds of evidence and rejoins at the end.

```text
BAM/CRAM ─► call ─► sample.inq ─► combine ─► combined.inq
                                                   │
                            ┌──────────────────────┴──────────────────────┐
                            │  locus axis                    allele axis  │
                            ▼                                             ▼
                    varindex [--catalog]                              outlier
                            │                                             │
                            ▼                                             ▼
                      cohort.stats                                 outlier_calls
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
```

Starting from a combined file produced by `inquiSTR combine` or `inquiSTR batch`:

```bash
# 1. Describe the loci: per-locus length distribution and variability index
inquiSTR varindex combined.inq -o cohort.stats

# 2. Find unusual alleles: each sample tested against the rest of the cohort
inquiSTR outlier combined.inq --method percentile -o outlier_calls.tsv

# 3. Rank per sample, keeping only the most variable loci
inquiSTR prioritize outlier_calls.tsv --index cohort.stats --top-index-pct 5 -o shortlist.tsv
```

Steps 1 and 2 are independent and can run at the same time. Step 3 reads only the two small
tables, so re-running it at a different cutoff takes seconds.

Optionally, annotate your catalog with the published PLVI values before ranking. The external
index is reported alongside each candidate for comparison; **ranking and filtering still use the
locally computed index**, since the two measure different things and averaging them would hide
exactly the disagreements worth looking at:

```bash
inquiSTR annotate --plvi <reference-table> --catalog <your-catalog.bed> --out mapping.tsv
```

### Choosing a reference

`outlier` compares each allele against the rest of your own cohort, with the tested sample's
alleles removed. That exclusion is what lets a carrier be seen at all: without it, an individual
holding the longest allele at a locus is measured against their own allele and never stands out.

> **Not yet available:** scoring against an external population table instead of your own cohort
> is specified in `docs/PLVI_PORT.md` but is **not implemented**. It matters most for small
> cohorts and for cohorts ascertained on a phenotype, where the cohort's own tail is a poor null.
> Until it exists, `varindex` output from a larger cohort can serve as an external *index* for
> `prioritize`, but the outlier threshold is always internal.

### Resource expectations

A combined file for 500 samples against a 4.4M-locus catalog is on the order of 10 GB of text
(measured: 9.5 GB for 434 samples against 4.86M loci, 910 MB gzipped).
`varindex` and `outlier` each stream it once and are parallel across loci; everything downstream
works on tables that are three to four orders of magnitude smaller. Plan for one pass per axis,
and expect re-ranking to be effectively free.

## How it works

### The variability index

For each locus, `varindex` computes the standard deviation of allele length across the cohort,
then ranks loci **within groups of comparable motif length**. The grouping matters: a
dinucleotide repeat and a pentanucleotide repeat have systematically different length
variability, and comparing them directly just sorts loci by motif size. Ranking within
motif-length groups removes that and leaves the property of interest — how unstable this locus
is *for its kind*.

Motif lengths are not evenly populated. In a typical genome-wide catalog, motifs of 1 to 6 bp
account for around 99% of loci, while the remainder is spread across a long tail of rare motif
lengths — often only a handful of loci each. A group of one cannot be ranked: that locus is
automatically in its own top percentile and would pass any cutoff regardless of its behaviour.
Adjacent motif lengths are therefore merged until each group is large enough to rank
meaningfully. Common motif lengths keep their own group; only the sparse tail is combined. The
group boundaries actually used are recorded in the output, since a rank is uninterpretable
without knowing what it was ranked against.

The output is a per-locus table: allele-length percentiles, mean, standard deviation, median
absolute deviation, mode, and the index itself. That table serves two purposes. It is the locus
evidence used for ranking, and it is a reference distribution that can be handed to `outlier`
for a different cohort.

### Finding outlier alleles

`outlier` tests every allele of every sample against the reference distribution at that locus.
Four methods are available:

| method | basis |
| --- | --- |
| `percentile` (default) | how many reference alleles are longer than this one |
| `robustz` | distance from the locus median, scaled by median absolute deviation |
| `zscore` | distance from the locus mean, scaled by standard deviation |
| `dbscan` | density-based clustering, outliers are points in no cluster |

`percentile` is the default because it makes no assumption about the shape of the distribution.
Repeat length distributions are frequently multimodal and heavy-tailed, which is exactly where
mean-and-standard-deviation methods behave worst.

There is a deeper reason to prefer it here. A z-score divides by the locus standard deviation —
the same quantity the variability index rewards. Sensitivity therefore drops at precisely the
high-variability loci that prioritization is designed to highlight, so the two steps work against
each other. Rank- and tail-based methods do not have this problem. The same tension appears in
`dbscan` through its neighbourhood radius.

Every method returns a **score**, a **significance**, and **evidence** — never a bare yes/no.
Scores are comparable across loci, which is what makes ranking a whole genome meaningful.

### Leave-one-out, and why it matters

When the reference is your own cohort, a sample is never included in the distribution it is
tested against. This is not a technicality. If a carrier holds the longest allele at a locus,
including them makes their allele the threshold, and they score as unremarkable — the one
individual you were looking for is the one the test cannot see.

For the same reason, related individuals should be excluded from each other's reference. A pair
of siblings sharing an expanded allele will each suppress the other's significance.

### Ranking

`prioritize` joins the outlier calls to the locus table, keeps loci above an index cutoff
(default: the top 5%), applies quality control, and ranks what remains within each sample.

It also reports the **funnel** — how many candidates survive each stage. This is the fastest way
to tell whether your thresholds are sensible: if a stage removes almost nothing, or almost
everything, the numbers say so immediately.

### Quality control

Quality control affects the false positive rate more than the choice of statistic does.

- Loci where the reference cohort frequently fails to genotype are not tested, because a locus
  with many failed calls has an artificially truncated tail.
- A sample producing far more outliers than the cohort median is flagged. That pattern is
  usually a coverage or alignment problem, not biology.
- **A failed genotype at a highly variable locus is reported as its own event**, not silently
  dropped. In length-based genotyping, failure to call a very long allele is itself a plausible
  expansion signature, and it is invisible to any method that only examines successful calls.

## How this differs from the PLVI paper

This implementation follows Danzi et al. (bioRxiv 2025.01.06.631535) but is not a reimplementation
of it. The differences below are deliberate and affect interpretation.

### Length, not longest pure segment

This is the fundamental difference. The published index is built on the **longest pure segment**
(LPS) — the longest run of uninterrupted repeat motif within an allele. inquiSTR measures total
allele length from alignment, and does not resolve allele sequence, so it cannot compute LPS.

The published work quantifies exactly what this costs. Measured as enrichment for known
disease-associated loci:

| measure | odds ratio |
| --- | --- |
| variability of total allele length | 34.1 |
| variability of LPS length (PLV) | 68.6 |
| PLV ranked within motif-length groups (PLVI) | 243.9 |

A length-based index is a real and useful signal — 36 of 42 known disease loci fall in the top
5% by length variability alone — but it is roughly an order of magnitude weaker than PLVI.
Expect a longer candidate list for the same sensitivity.

Because of this, the length-based index computed here is never called PLVI. Output columns state
both what was measured and where it came from:

| column | what it is |
| --- | --- |
| `lvi_length_internal` | length variability, computed from your cohort |
| `lvi_length_external` | length variability, from a published population cohort |
| `plvi_lps_external` | the published LPS-based PLVI, imported as an annotation |

An imported LPS-derived value is never substituted for a locally computed length-based one, and
the two are never averaged. They measure different things.

### Interruptions are invisible

Repeat interruptions strongly influence expansion instability and clinical presentation, and are
a large part of why LPS outperforms raw length. A length-based analysis cannot see them. An
allele that is long but heavily interrupted and one that is long and pure look identical here.

### A reference cohort of your own

The published analysis used a dedicated 543-genome reference cohort. Here the default reference
is your own cohort, leave-one-out. This is convenient but carries an assumption worth stating:
if a substantial fraction of your cohort shares a disease, the causal locus's variability is
inflated by the very signal you are trying to find, and ranking within your own cohort partly
cancels it.

**For phenotype-ascertained cohorts an external reference would be the right answer, but
scoring against one is not implemented yet** (see *Choosing a reference* above); until it is,
treat internally computed rankings on such a cohort with corresponding caution. Cohort size is not the concern
— 500 individuals is 1,000 alleles, comparable to the published reference cohorts. Cohort
composition is.

### Thresholds are reported, not just applied

The published outlier definition is "longer than the 99.9th percentile of the reference". At
realistic cohort sizes that threshold is frequently saturated: in the published ONT reference of
998 alleles, around 80% of loci have a 99.9th percentile equal to the longest allele ever
observed there. At those loci the rule is really "longer than anything we have seen", and no
threshold can be finer.

That is not a flaw so much as a limit of the data, but it should be visible rather than implied.
Every call therefore reports how many reference alleles it exceeded and how many there were, so
a nominal percentile that has degenerated into a maximum is apparent.

### Counting per sample, not per allele

The published funnel counts each allele separately and sums them, so a locus where both alleles
are outliers contributes twice. Here a sample is counted once per locus. Funnel output reports
both figures so the published numbers remain comparable.

By default both alleles are tested. For biallelic conditions such as RFC1 or FXN, where the
*shorter* allele being long is the signal, the allele selection is configurable — testing only
the longer allele silently loses recessive disease.

### Failed genotypes

In the published analysis, uncalled genotypes are treated as zero and therefore never appear as
outliers. Here they are tracked and reported. This produces different counts than the published
figures, deliberately.

## References

Danzi MC, Xu IRL, et al. *Population-scale variability at short tandem repeat loci reveals
pathogenicity signature.* bioRxiv 2025.01.06.631535.

For the *GOLGA8A* example and the association-first argument: *A repeat expansion in GOLGA8A is
a major risk factor for atypical frontotemporal lobar degeneration with ubiquitin-positive
inclusions*, Nature Genetics 58, 726–736 (2026), doi:10.1038/s41588-026-02537-7 — which used
inquiSTR for the genome-wide repeat-length association and a sequence-resolved genotyper to
characterise the motif composition at the associated locus.
