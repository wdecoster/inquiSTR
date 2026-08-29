# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

inquiSTR is a single-binary Rust CLI (edition 2024, MSRV 1.88) for genotyping and analyzing short tandem repeats (STRs / tandem repeats) from long-read sequencing data (BAM/CRAM, ONT-focused). Everything lives in `src/` as flat modules declared in `src/main.rs`; there is no `lib.rs`, so most tests are inline `#[cfg(test)]` modules. `tests/` holds the exception: end-to-end tests that drive the compiled binary via `env!("CARGO_BIN_EXE_inquiSTR")`, for behaviour an inline test cannot reach (a subcommand writing to stdout, for instance). Prefer an inline test unless the behaviour only exists at the process boundary.

## Commands

```bash
make ci                 # fmt-check + clippy + test — run this before considering a change done
make fmt                # cargo fmt
make clippy             # cargo clippy --all-targets --all-features -- -D warnings
make test               # cargo test
make build              # cargo build --release
make musl               # static x86_64-unknown-linux-musl binary (uses cross if available)
make install-hooks      # install .githooks/{pre-commit,pre-push} into .git/hooks

cargo test <name>                              # run a single test, e.g. cargo test test_unphased
cargo test -- --ignored --nocapture            # run network-dependent tests (preset URL downloads)
cargo test test_preset_urls -- --ignored       # verify the preset catalog URLs still resolve
```

Tests use relative paths into `test-data/`, so run them from the repository root. CI runs both a native and a musl (`cross`) build, and treats clippy warnings as errors — formatting/clippy is the most common CI failure. Clippy thresholds are relaxed in `.config/clippy.toml` (100-line functions, 10 arguments, cognitive complexity 30).

Note: `.rustfmt.toml` declares `edition = "2021"` while `Cargo.toml` is `edition = "2024"`; leave that alone unless fixing it deliberately — changing it reformats the tree.

## Architecture

### Subcommand dispatch
`src/main.rs` defines the entire clap `Commands` enum and dispatches to one module per subcommand. Shared clap arg groups (`BatchCommonArgs`, mode selection) also live there. Adding a subcommand means: new module, `pub mod` in main.rs, new `Commands` variant, and a match arm.

### The genotyping pipeline (`call` / `batch`)
The core path is `repeats.rs` → `batching.rs` → `genotype_batch.rs` → `call.rs`, with `batch_process.rs` wrapping it for many samples:

1. **`repeats.rs`** resolves targets from `--region`, `--region-file`, or `--preset`. `TRPreset` (pathogenic/adotto/trexplorer/codis) maps to a remote catalog URL that is downloaded and cached under `$XDG_CACHE_HOME` (or `~/.cache`) with a lock file to make concurrent runs safe. `ChromosomeMapper` reconciles `chr`-prefixed vs unprefixed contig naming between the catalog and the BAM header. Adding a preset = new enum variant + entry in `TRPreset::metadata()`.
2. **`batching.rs`** groups targets on the same chromosome within `--batch-size` KB into a `Batch`, so one indexed BAM fetch serves many loci. This is the main I/O performance lever (`optimize.rs` / `optimize-call` tunes it empirically).
3. **`genotype_batch.rs`** does the actual calling: walks CIGARs of reads overlapping the batch, records signed indel/soft-clip length per read as a `Call` (`CallType` spanning vs clipped, `Phase` from the `HP` tag), then `median_str_length` produces a per-haplotype median. `docs/CALL_ALGORITHM.md` is the authoritative description of this logic — keep it in sync with changes here.
4. **`bam_pool.rs`** hands each rayon worker its own BAM reader instead of locking a shared one, while respecting file-descriptor limits.
5. **`batch_process.rs`** runs the manifest (`bam_path` + optional `sample_name` TSV) across samples. Default threading is *one thread per sample* (CRAM-friendly); `--parallel-samples` overrides the split. Per-sample panics are caught (`catch_sample_panic`) so `--keep-going` and `--resume` work; intermediate per-sample files go to a tmpdir and are merged via `combine`.

### File formats and type dispatch
All output files carry a `#`-prefixed metadata header (`file_type`, `version`, `command`, plus command-specific keys). `filetype.rs` reads that header to decide what a file is (`individual_call`, `combined_call`, `individual_kmer`, `combined_kmer`, `target_kmer`) instead of guessing from content — downstream commands validate with it and fail early on a mismatched input. `docs/FILE_TYPE_METADATA.md` and `docs/OUTPUT_FORMATS.md` document the formats; both must be updated when a column or metadata key changes.

Genotype values are the **median length difference from reference in bp** (positive = expansion), with `NaN` for uncallable loci — not repeat-unit counts.

### Downstream analysis
`combine.rs` (merge/extend cohort files, gzip-transparent readers), then `filter.rs`, `outlier.rs`, `query.rs`, `histogram.rs`/`plot.rs`, `pca.rs`, `assoc.rs` (native Rust linear/logistic regression via nalgebra), `relate.rs`, `benchmark.rs`, `convert.rs` (VCF → inquiSTR). `locus_search.rs` is the shared "find this region in a combined file" helper used by query/plot/histogram; `sample_info.rs` parses phenotype/group files (distinct from the `#` header metadata).

`unmapped.rs` is a separate track: kmer frequency counting over unmapped reads, producing kmer-type files that flow through the same combine/pca/assoc/outlier machinery.

### Dependency notes
- `rust-htslib` and `kuva` (plotting) are pinned to **git forks** under `wdecoster/`. `kuva` is on a feature branch pending upstream merge — revert to crates.io once merged.
- `histo.rs` is a vendored, dependency-free replacement for the abandoned `histo_fp` crate; `histogram.rs` is the subcommand that uses it. Don't confuse the two.
- Release profile uses `panic = "abort"`, so panic-catching in batch processing only works in dev/test builds — user-facing errors should go through `errors.rs` (`InquiSTRError`), not panics.

## Conventions

From `.github/copilot-instructions.md`: after code changes, run `cargo fmt`, `cargo clippy`, and `cargo test` and fix issues before finishing. New features must also update doc comments, README, and clap help strings. Keep documentation user-facing — what it does and when to use it, not internals — with at most one concise example per feature.
