//! Reader for per-locus statistics tables.
//!
//! One format, two origins: a table written by `varindex` from a local cohort, and a published
//! population table. They carry the same statistics but key their rows differently — inquiSTR
//! uses explicit coordinate columns, the published tables use a compound `TRID` of the form
//! `chrom-start-end-motif`. Both are accepted here so that everything downstream has a single
//! way to read a reference distribution, whatever produced it.
//!
//! Columns are located by name rather than position, so extra or reordered columns are
//! tolerated and a missing required column is an error rather than a silent misread.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

/// Locus identity, normalised so that `chr1` and `1` compare equal
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct LocusKey {
    pub chrom: String,
    pub begin: u32,
    pub end: u32,
}

impl LocusKey {
    pub fn new(chrom: &str, begin: u32, end: u32) -> Self {
        LocusKey { chrom: normalise_chrom(chrom), begin, end }
    }
}

/// Drop a `chr` prefix so that catalogs and published tables can be joined.
pub fn normalise_chrom(chrom: &str) -> String {
    chrom.strip_prefix("chr").unwrap_or(chrom).to_string()
}

/// One row of a locus statistics table
#[derive(Clone, Debug)]
pub struct LocusStatsRow {
    pub key: LocusKey,
    pub motif: String,
    /// Called alleles behind the statistics
    pub n: usize,
    pub std: f64,
    /// The variability index, when the table carries one
    pub index: Option<f64>,
    /// Percentile values, in the order given by [`StatsHeader::percentiles`]
    pub percentiles: Vec<f64>,
}

impl LocusStatsRow {
    /// Value of the named percentile, e.g. 99.9
    pub fn percentile(&self, header: &StatsHeader, p: f64) -> Option<f64> {
        header
            .percentiles
            .iter()
            .position(|(q, _)| (*q - p).abs() < 1e-9)
            .and_then(|i| self.percentiles.get(i).copied())
    }
}

/// Metadata and column layout of a locus statistics table
#[derive(Clone, Debug, Default)]
pub struct StatsHeader {
    pub samples: Option<usize>,
    pub index_column: Option<String>,
    pub motif_bins: Option<String>,
    pub length_axis: Option<String>,
    pub source: Option<String>,
    /// (percentile, column position)
    pub percentiles: Vec<(f64, usize)>,
    col_chrom: Option<usize>,
    col_begin: Option<usize>,
    col_end: Option<usize>,
    col_trid: Option<usize>,
    col_motif: Option<usize>,
    col_n: Option<usize>,
    col_std: Option<usize>,
    col_index: Option<usize>,
}

/// Parse `99.9thPercentile` into 99.9
fn percentile_of(name: &str) -> Option<f64> {
    let stem = name.strip_suffix("Percentile")?;
    let stem = stem
        .strip_suffix("th")
        .or_else(|| stem.strip_suffix("st"))
        .or_else(|| stem.strip_suffix("nd"))
        .or_else(|| stem.strip_suffix("rd"))
        .unwrap_or(stem);
    stem.parse::<f64>().ok()
}

/// Split a published `TRID` such as `1-160919121-160919141-A`.
///
/// The motif may itself be long, and coordinates are plain integers, so the identifier is split
/// from the left on the first two separators and the remainder is taken as the motif.
pub fn parse_trid(trid: &str) -> Option<(String, u32, u32, String)> {
    let mut it = trid.splitn(4, '-');
    let chrom = it.next()?;
    let begin: u32 = it.next()?.parse().ok()?;
    let end: u32 = it.next()?.parse().ok()?;
    let motif = it.next().unwrap_or(".").to_string();
    Some((normalise_chrom(chrom), begin, end, motif))
}

/// Read the metadata and column header of a locus statistics table.
/// Returns the header plus the reader positioned at the first data line.
pub fn open(path: &Path) -> (StatsHeader, Box<dyn BufRead>) {
    let mut reader: Box<dyn BufRead> = Box::new(crate::utils::reader(&path.to_string_lossy()));
    let mut header = StatsHeader::default();
    let column_line: String;

    loop {
        let mut line = String::new();
        match reader.read_line(&mut line) {
            Ok(0) => {
                eprintln!("ERROR: {} contains no data", path.display());
                std::process::exit(1);
            }
            Ok(_) => {}
            Err(e) => {
                eprintln!("ERROR: Failed reading {}: {}", path.display(), e);
                std::process::exit(1);
            }
        }
        let trimmed = line.trim_end().to_string();
        if let Some(rest) = trimmed.strip_prefix("# ") {
            if let Some((k, v)) = rest.split_once('=') {
                match k {
                    "samples" => header.samples = v.parse().ok(),
                    "index_column" => header.index_column = Some(v.to_string()),
                    "motif_bins" => header.motif_bins = Some(v.to_string()),
                    "length_axis" => header.length_axis = Some(v.to_string()),
                    "source" => header.source = Some(v.to_string()),
                    _ => {}
                }
            }
            continue;
        }
        column_line = trimmed;
        break;
    }

    for (i, name) in column_line.split('\t').enumerate() {
        match name {
            "chromosome" | "chrom" => header.col_chrom = Some(i),
            "begin" | "start" => header.col_begin = Some(i),
            "end" => header.col_end = Some(i),
            "TRID" => header.col_trid = Some(i),
            "motif" | "LPS_Motif" => header.col_motif = Some(i),
            "N" => header.col_n = Some(i),
            "std" => header.col_std = Some(i),
            _ => {
                if let Some(p) = percentile_of(name) {
                    header.percentiles.push((p, i));
                } else if name.starts_with("lvi_") || name.starts_with("plvi_") {
                    // A table may carry several index columns (an internal and an external one
                    // side by side). Prefer the one the metadata header names; otherwise take
                    // the first. Taking the last would silently filter on a different column
                    // from the one reported downstream.
                    let named = header.index_column.as_deref() == Some(name);
                    if named || header.col_index.is_none() {
                        header.col_index = Some(i);
                        if header.index_column.is_none() {
                            header.index_column = Some(name.to_string());
                        }
                    }
                }
            }
        }
    }

    // If the metadata names an index column the table does not actually contain, report the
    // column that will really be read rather than the name that was asked for.
    if let Some(i) = header.col_index {
        let actual = column_line.split('\t').nth(i).map(str::to_string);
        if header.index_column != actual {
            header.index_column = actual;
        }
    } else {
        header.index_column = None;
    }

    let keyed_by_coordinates =
        header.col_chrom.is_some() && header.col_begin.is_some() && header.col_end.is_some();
    if !keyed_by_coordinates && header.col_trid.is_none() {
        eprintln!(
            "ERROR: {} has neither chromosome/begin/end columns nor a TRID column; \
             it does not look like a locus statistics table.",
            path.display()
        );
        std::process::exit(1);
    }

    (header, reader)
}

/// Parse one data line against a header
pub fn parse_row(header: &StatsHeader, line: &str) -> Option<LocusStatsRow> {
    let f: Vec<&str> = line.split('\t').collect();
    let get = |i: Option<usize>| -> Option<&str> { i.and_then(|i| f.get(i).copied()) };

    let (chrom, begin, end, motif) = if let Some(trid) = get(header.col_trid) {
        let (c, b, e, m) = parse_trid(trid)?;
        // An explicit motif column wins over the one embedded in the identifier
        let m = get(header.col_motif)
            .filter(|s| *s != "." && !s.is_empty())
            .unwrap_or(&m)
            .to_string();
        (c, b, e, m)
    } else {
        (
            normalise_chrom(get(header.col_chrom)?),
            get(header.col_begin)?.parse().ok()?,
            get(header.col_end)?.parse().ok()?,
            get(header.col_motif).unwrap_or(".").to_string(),
        )
    };

    let parse_f =
        |s: Option<&str>| -> f64 { s.and_then(|v| v.parse::<f64>().ok()).unwrap_or(f64::NAN) };

    let index = header
        .col_index
        .and_then(|i| f.get(i))
        .and_then(|v| v.parse::<f64>().ok());
    let percentiles: Vec<f64> = header
        .percentiles
        .iter()
        .map(|(_, i)| parse_f(f.get(*i).copied()))
        .collect();

    Some(LocusStatsRow {
        key: LocusKey { chrom, begin, end },
        motif,
        n: get(header.col_n).and_then(|v| v.parse().ok()).unwrap_or(0),
        std: parse_f(get(header.col_std)),
        index: index.filter(|v| v.is_finite()),
        percentiles,
    })
}

/// Stream a table, keeping only rows whose locus appears in `wanted`.
///
/// Filtering during the scan keeps memory proportional to the loci actually needed rather than
/// to the size of the catalog, which matters when the table has millions of rows and only a
/// few hundred thousand are of interest.
pub fn load_subset(
    path: &Path,
    wanted: &std::collections::HashSet<LocusKey>,
) -> (StatsHeader, HashMap<LocusKey, LocusStatsRow>) {
    let (header, reader) = open(path);
    let mut map = HashMap::with_capacity(wanted.len());
    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("ERROR: Failed reading {}: {}", path.display(), e);
                std::process::exit(1);
            }
        };
        if let Some(row) = parse_row(&header, &line)
            && wanted.contains(&row.key)
        {
            map.insert(row.key.clone(), row);
        }
    }
    (header, map)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(content: &str, suffix: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn test_percentile_column_names() {
        assert_eq!(percentile_of("0thPercentile"), Some(0.0));
        assert_eq!(percentile_of("1stPercentile"), Some(1.0));
        assert_eq!(percentile_of("99.9thPercentile"), Some(99.9));
        assert_eq!(percentile_of("100thPercentile"), Some(100.0));
        assert_eq!(percentile_of("mean"), None);
    }

    #[test]
    fn test_parse_trid_handles_long_motifs() {
        let (c, b, e, m) = parse_trid("1-160919121-160919141-A").unwrap();
        assert_eq!((c.as_str(), b, e, m.as_str()), ("1", 160919121, 160919141, "A"));
        // A motif containing no separators of its own, but long
        let (_, _, _, m) = parse_trid("chr7-1000-2000-CCTCATGGTGGTGGCTGGGGGCAG").unwrap();
        assert_eq!(m, "CCTCATGGTGGTGGCTGGGGGCAG");
    }

    #[test]
    fn test_chromosome_naming_is_reconciled() {
        assert_eq!(LocusKey::new("chr1", 10, 20), LocusKey::new("1", 10, 20));
    }

    #[test]
    fn test_reads_inquistr_native_layout() {
        let content = "# file_type=locus_stats\n# samples=4\n# index_column=lvi_length_internal\n\
            chromosome\tbegin\tend\tmotif\tN\t50thPercentile\t99.9thPercentile\tstd\tlvi_length_internal\n\
            chr1\t100\t130\tCAG\t8\t30\t90\t4.5\t0.987\n";
        let f = write_tmp(content, ".tsv");
        let (header, reader) = open(f.path());
        assert_eq!(header.samples, Some(4));
        assert_eq!(header.percentiles.len(), 2);
        let line = reader.lines().next().unwrap().unwrap();
        let row = parse_row(&header, &line).unwrap();
        assert_eq!(row.key, LocusKey::new("1", 100, 130));
        assert_eq!(row.motif, "CAG");
        assert_eq!(row.n, 8);
        assert_eq!(row.index, Some(0.987));
        assert_eq!(row.percentile(&header, 99.9), Some(90.0));
    }

    #[test]
    fn test_reads_published_trid_layout_without_an_index() {
        // The published tables have no index column; that must not be an error, and must not
        // be mistaken for an index of zero.
        let content = "TRID\tN\t50thPercentile\t99.9thPercentile\tstd\n\
            1-100-130-CAG\t1086\t30\t90\t4.5\n";
        let f = write_tmp(content, ".tsv");
        let (header, reader) = open(f.path());
        assert!(header.col_index.is_none());
        let line = reader.lines().next().unwrap().unwrap();
        let row = parse_row(&header, &line).unwrap();
        assert_eq!(row.key, LocusKey::new("chr1", 100, 130));
        assert_eq!(row.motif, "CAG");
        assert_eq!(row.n, 1086);
        assert_eq!(row.index, None, "absent index must stay absent");
    }

    #[test]
    fn test_named_index_column_wins_over_others_present() {
        let content = "# index_column=lvi_length_internal\n\
            chromosome\tbegin\tend\tmotif\tN\tstd\tplvi_lps_external\tlvi_length_internal\n\
            chr1\t100\t130\tCAG\t8\t4.5\t0.01\t0.99\n";
        let f = write_tmp(content, ".tsv");
        let (header, reader) = open(f.path());
        let line = reader.lines().next().unwrap().unwrap();
        let row = parse_row(&header, &line).unwrap();
        assert_eq!(header.index_column.as_deref(), Some("lvi_length_internal"));
        assert_eq!(row.index, Some(0.99), "must read the column the header names");
    }

    #[test]
    fn test_columns_located_by_name_not_position() {
        // Same data, columns reordered and an unknown column added
        let content = "# samples=4\n\
            extra\tstd\tend\tN\tbegin\tmotif\tchromosome\tlvi_length_internal\n\
            junk\t4.5\t130\t8\t100\tCAG\tchr1\t0.5\n";
        let f = write_tmp(content, ".tsv");
        let (header, reader) = open(f.path());
        let line = reader.lines().next().unwrap().unwrap();
        let row = parse_row(&header, &line).unwrap();
        assert_eq!(row.key, LocusKey::new("1", 100, 130));
        assert_eq!(row.std, 4.5);
        assert_eq!(row.index, Some(0.5));
    }
}
