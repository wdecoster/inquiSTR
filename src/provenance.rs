//! Content fingerprints for the catalog and the reference genome.
//!
//! A catalog *name* proves nothing: catalogs are updated in place at the same URL, and a user
//! can pass a locally edited BED. Recording a hash of the content lets `combine` refuse to
//! merge files that were called against different loci, and lets an external join refuse a
//! mismatched genome build, instead of silently producing coordinates that do not correspond.
//!
//! Two fingerprints, answering different questions:
//!
//! - `catalog_sha256` — are two inquiSTR files talking about the same loci?
//! - `genome_sha256` — which reference build? Chromosome lengths differ between GRCh38, T2T
//!   and hg19, so the length map identifies the build without needing the FASTA.
//!
//! **Cost.** Hashing is memoised. The digest is computed once per catalog and cached beside it,
//! keyed by size and modification time, because `call` runs once per sample and a batch would
//! otherwise re-hash an identical file for every one of them. On a machine with SHA-NI the
//! digest itself adds roughly 6-9% to the decompress-and-read the catalog already requires;
//! after the first run it costs nothing.

use sha2::{Digest, Sha256};

/// Lowercase hex, so a digest printed here matches what `sha256sum` prints
fn hex(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        use std::fmt::Write;
        let _ = write!(out, "{:02x}", b);
    }
    out
}
use std::collections::HashMap;
use std::io::Read;
use std::path::{Path, PathBuf};

/// Hash of the **decompressed** catalog content.
///
/// Decompressed rather than raw: catalogs ship gzipped, gzip output is not reproducible across
/// compressors, and hashing the container would make an identical catalog look changed after
/// any recompression. Decompressed content is also verifiable by hand with
/// `zcat catalog.bed.gz | sha256sum`, which matters when a hash is what causes a refusal.
pub fn catalog_sha256(path: &Path) -> Option<String> {
    // Provenance is an optional field, so a missing catalog must degrade to "unknown" rather
    // than abort: `utils::reader` exits the process when it cannot open a file.
    if !path.is_file() {
        return None;
    }
    if let Some(cached) = read_cached(path) {
        return Some(cached);
    }
    let before = stamp(path)?;
    let digest = hash_decompressed(path)?;
    write_cached(path, &digest, before);
    Some(digest)
}

fn hash_decompressed(path: &Path) -> Option<String> {
    let mut reader = crate::utils::reader(&path.to_string_lossy());
    let mut hasher = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        match reader.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => hasher.update(&buf[..n]),
            Err(e) => {
                eprintln!("Warning: could not hash {}: {}", path.display(), e);
                return None;
            }
        }
    }
    Some(hex(&hasher.finalize()))
}

/// Fingerprint of the reference build, from the chromosome length map.
///
/// The map is already loaded to validate target coordinates, so this costs nothing measurable.
/// Names are normalised and sorted so that `chr1` and `1`, or a differently ordered header,
/// still identify the same build.
pub fn genome_sha256(chrom_lengths: &HashMap<String, u64>) -> String {
    let mut entries: Vec<(String, u64)> = chrom_lengths
        .iter()
        .map(|(k, v)| (crate::locus_stats::normalise_chrom(k), *v))
        .collect();
    entries.sort();
    let mut hasher = Sha256::new();
    for (name, len) in entries {
        hasher.update(name.as_bytes());
        hasher.update(b"\t");
        hasher.update(len.to_string().as_bytes());
        hasher.update(b"\n");
    }
    hex(&hasher.finalize())
}

/// Where the memoised digest for `path` lives
fn cache_path(path: &Path) -> Option<PathBuf> {
    // Reuses the catalog cache location rather than re-deriving it: a second copy of that logic
    // drifted from the original (no macOS branch, no absolute-path filter) and scattered
    // fingerprint directories into whatever the working directory happened to be.
    let dir = crate::repeats::preset_cache_dir().join("fingerprints");
    std::fs::create_dir_all(&dir).ok()?;
    // Key on the full path so two catalogs with the same basename cannot collide
    let mut hasher = Sha256::new();
    hasher.update(
        path.canonicalize()
            .unwrap_or_else(|_| path.to_path_buf())
            .to_string_lossy()
            .as_bytes(),
    );
    Some(dir.join(format!("{}.sha256", hex(&hasher.finalize()))))
}

/// File identity: size, nanosecond mtime, inode and device
#[derive(PartialEq, Eq, Clone, Copy)]
struct Stamp {
    size: u64,
    mtime: u128,
    ino: u64,
    dev: u64,
}

impl Stamp {
    fn encode(&self) -> String {
        format!("{}\t{}\t{}\t{}", self.size, self.mtime, self.ino, self.dev)
    }
    fn decode(fields: &mut std::str::Split<'_, char>) -> Option<Self> {
        Some(Stamp {
            size: fields.next()?.parse().ok()?,
            mtime: fields.next()?.parse().ok()?,
            ino: fields.next()?.parse().ok()?,
            dev: fields.next()?.parse().ok()?,
        })
    }
}

/// Identity stamp used to notice that a catalog changed under a cached digest.
///
/// Nanosecond mtime, not whole seconds: a pipeline that regenerates a BED at a fixed path can
/// rewrite it to the same size within one second, and a second-resolution stamp would then serve
/// the previous file's digest. Inode and device are included too, because `rsync -a`, `cp -p`
/// and `tar -x` restore the *source* mtime, so a same-size replacement can arrive carrying an
/// arbitrarily old timestamp.
fn stamp(path: &Path) -> Option<Stamp> {
    use std::os::unix::fs::MetadataExt;
    let meta = std::fs::metadata(path).ok()?;
    let mtime = meta
        .modified()
        .ok()?
        .duration_since(std::time::UNIX_EPOCH)
        .ok()?
        .as_nanos();
    Some(Stamp { size: meta.len(), mtime, ino: meta.ino(), dev: meta.dev() })
}

fn read_cached(path: &Path) -> Option<String> {
    let current = stamp(path)?;
    let contents = std::fs::read_to_string(cache_path(path)?).ok()?;
    let mut fields = contents.trim().split('\t');
    let digest = fields.next()?;
    let cached = Stamp::decode(&mut fields)?;
    // Validate the shape as well as the stamp: a corrupt memo must cause a recompute, never a
    // "digest" that was never a digest.
    let looks_like_a_digest = digest.len() == 64 && digest.bytes().all(|b| b.is_ascii_hexdigit());
    if cached == current && looks_like_a_digest {
        Some(digest.to_string())
    } else {
        None
    }
}

fn write_cached(path: &Path, digest: &str, before: Stamp) {
    // The stamp is taken before hashing and re-checked here. If the catalog changed while it was
    // being read, the digest describes content no longer at this path, and pairing it with the
    // new stamp would serve a wrong answer indefinitely.
    if stamp(path) != Some(before) {
        return;
    }
    let Some(cache) = cache_path(path) else {
        return;
    };
    // Written atomically. A workflow manager will start many `call` processes at once, all of
    // them missing the cache and all of them writing it; a torn file would then be read back as
    // a miss for every later run. Write to a unique temporary name and rename over the target,
    // which is atomic within a filesystem.
    let tmp = cache.with_extension(format!("tmp.{}", std::process::id()));
    if std::fs::write(&tmp, format!("{}\t{}\n", digest, before.encode())).is_ok()
        && std::fs::rename(&tmp, &cache).is_err()
    {
        let _ = std::fs::remove_file(&tmp);
    }
}

/// The provenance fields recorded in an inquiSTR output header.
///
/// Captured once and shared, so that every path that writes a call file records the same thing.
#[derive(Debug, Clone, Default)]
pub struct Provenance {
    pub catalog: Option<String>,
    pub catalog_sha256: Option<String>,
    pub genome_sha256: Option<String>,
}

impl Provenance {
    /// Capture provenance for a run.
    ///
    /// Must be called *after* the targets are resolved: a preset catalog is downloaded or
    /// refreshed on demand, and hashing beforehand would digest a stale copy or a file that
    /// does not exist yet.
    pub fn capture(
        targets: &crate::repeats::TargetConfig,
        bam: &str,
        reference: &Option<String>,
    ) -> Self {
        Provenance {
            catalog: targets.catalog_name(),
            catalog_sha256: targets.catalog_path().as_deref().and_then(catalog_sha256),
            genome_sha256: crate::bam_utils::get_chrom_lengths_from_bam_header(
                bam.to_string(),
                reference,
            )
            .ok()
            .map(|m| genome_sha256(&m)),
        }
    }

    /// Write the fields that are known, as `# key=value` metadata lines
    pub fn write<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()> {
        for (key, value) in [
            ("catalog", &self.catalog),
            ("catalog_sha256", &self.catalog_sha256),
            ("genome_sha256", &self.genome_sha256),
        ] {
            if let Some(v) = value {
                writeln!(w, "# {}={}", key, v)?;
            }
        }
        Ok(())
    }
}

/// Compare a provenance field between two files.
///
/// A missing field is not a mismatch: files written before provenance was recorded have none,
/// and refusing them would make every existing `.inq` unreadable. Only a present and differing
/// value is an error.
pub fn fields_conflict(a: Option<&str>, b: Option<&str>) -> bool {
    matches!((a, b), (Some(x), Some(y)) if x != y)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_genome_fingerprint_ignores_naming_and_order() {
        let mut a = HashMap::new();
        a.insert("chr1".to_string(), 248_956_422u64);
        a.insert("chr2".to_string(), 242_193_529u64);
        let mut b = HashMap::new();
        b.insert("2".to_string(), 242_193_529u64);
        b.insert("1".to_string(), 248_956_422u64);
        assert_eq!(genome_sha256(&a), genome_sha256(&b), "chr prefix and order must not matter");

        // A different build has different chromosome lengths
        let mut c = HashMap::new();
        c.insert("chr1".to_string(), 249_250_621u64); // hg19 chr1
        c.insert("chr2".to_string(), 242_193_529u64);
        assert_ne!(genome_sha256(&a), genome_sha256(&c), "a different build must differ");
    }

    #[test]
    fn test_catalog_hash_is_of_content_not_container() {
        let mut plain = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        plain.write_all(b"chr1\t100\t130\tCAG\n").unwrap();
        plain.flush().unwrap();
        let expected = hash_decompressed(plain.path()).unwrap();

        // The same content, gzipped, must hash identically
        let gz = tempfile::Builder::new()
            .suffix(".bed.gz")
            .tempfile()
            .unwrap();
        let mut enc = flate2::write::GzEncoder::new(
            std::fs::File::create(gz.path()).unwrap(),
            flate2::Compression::default(),
        );
        enc.write_all(b"chr1\t100\t130\tCAG\n").unwrap();
        enc.finish().unwrap();
        assert_eq!(
            hash_decompressed(gz.path()).unwrap(),
            expected,
            "compression must not change the digest"
        );

        // And it must match what a user gets from the shell
        assert_eq!(expected.len(), 64);
    }

    #[test]
    fn test_memoised_digest_is_reused_and_invalidated() {
        // Redirect the memo into a temporary cache so the test never writes to the user's own
        let home = tempfile::tempdir().unwrap();
        unsafe { std::env::set_var("XDG_CACHE_HOME", home.path()) };

        let mut f = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        f.write_all(b"chr1\t100\t130\tCAG\n").unwrap();
        f.flush().unwrap();
        let first = catalog_sha256(f.path()).unwrap();
        assert_eq!(catalog_sha256(f.path()).unwrap(), first, "second call must hit the cache");

        // A same-size rewrite within the same second must still invalidate the memo. With a
        // whole-second mtime this returned the previous file's digest.
        std::fs::write(f.path(), b"chr1\t100\t130\tCAC\n").unwrap();
        let same_second = catalog_sha256(f.path()).unwrap();
        assert_ne!(same_second, first, "same-size rewrite must not serve a stale digest");
        assert_eq!(same_second, hash_decompressed(f.path()).unwrap());

        // A restored file carrying an old mtime must invalidate too, via inode identity
        std::fs::write(f.path(), b"chr1\t100\t130\tCAG\n").unwrap();
        assert_eq!(catalog_sha256(f.path()).unwrap(), first);
    }

    #[test]
    fn test_corrupt_memo_recomputes_rather_than_returning_junk() {
        let home = tempfile::tempdir().unwrap();
        unsafe { std::env::set_var("XDG_CACHE_HOME", home.path()) };
        let mut f = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        f.write_all(b"chr1\t100\t130\tCAG\n").unwrap();
        f.flush().unwrap();
        let good = catalog_sha256(f.path()).unwrap();

        // A 64-character non-hex first field must not be handed back as a digest
        let memo = cache_path(f.path()).unwrap();
        let st = stamp(f.path()).unwrap();
        std::fs::write(&memo, format!("{}\t{}\n", "Z".repeat(64), st.encode())).unwrap();
        assert_eq!(catalog_sha256(f.path()).unwrap(), good, "junk memo must be recomputed");

        for junk in ["", "not a memo", "abc\t1\t2"] {
            std::fs::write(&memo, junk).unwrap();
            assert_eq!(catalog_sha256(f.path()).unwrap(), good);
        }
    }

    #[test]
    fn test_missing_fields_are_not_a_conflict() {
        assert!(!fields_conflict(None, Some("abc")), "an older file has no field to conflict");
        assert!(!fields_conflict(Some("abc"), None));
        assert!(!fields_conflict(Some("abc"), Some("abc")));
        assert!(fields_conflict(Some("abc"), Some("def")), "present and different is an error");
    }
}

#[cfg(test)]
mod bench {
    use super::*;
    use std::path::Path;

    #[test]
    #[ignore] // timing probe, run with: cargo test --release bench_ -- --ignored --nocapture
    fn bench_catalog_fingerprint() {
        for p in [
            "repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz",
            concat!(env!("HOME"), "/.cache/inquistr/adotto_TRregions_v1.2.1.bed.gz"),
        ] {
            let path = Path::new(p);
            if !path.exists() {
                continue;
            }
            let t = std::time::Instant::now();
            let cold = hash_decompressed(path);
            let cold_ms = t.elapsed().as_millis();
            let t = std::time::Instant::now();
            let warm = catalog_sha256(path);
            let warm_ms = t.elapsed().as_millis();
            let t = std::time::Instant::now();
            let warm2 = catalog_sha256(path);
            let warm2_us = t.elapsed().as_micros();
            assert_eq!(cold, warm);
            assert_eq!(cold, warm2);
            println!(
                "{:<48} cold {:>5} ms | first {:>5} ms | memoised {:>6} us",
                path.file_name().unwrap().to_string_lossy(),
                cold_ms,
                warm_ms,
                warm2_us
            );
        }
    }
}
