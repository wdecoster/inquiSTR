//! Text histogram with a fixed number of buckets.
//!
//! This is a vendored, dependency-free replacement for `histo_fp` 0.2.1, whose
//! last release was in 2020 and which pulled in `rand` 0.7.3 (RUSTSEC-2026-0097)
//! along with eight other crates nothing else in inquiSTR used. It reproduces
//! that crate's bucketing, statistics and `Display` byte-for-byte;
//! `tests::golden_output` pins the format.
//!
//! Two upstream quirks are deliberately preserved because they are visible in
//! the output:
//!
//! - the bucket range is widened by `range as u64 % num_buckets` before being
//!   divided, a leftover from the integer-based crate this was forked from;
//! - both range columns are padded to the width of the widest *end* value, so a
//!   wider start value is left unpadded rather than widening the column.
//!
//! `histo_fp::Histogram::with_buckets` also took a `precision` argument. It was
//! stored and never read, so it is not reproduced here.

use std::cmp;
use std::collections::BTreeMap;
use std::fmt;

/// An `f64` usable as a `BTreeMap` key.
///
/// `histo_fp` used `noisy_float::R64`, which orders through
/// `partial_cmp().unwrap()` and so panics on NaN. Same contract here: callers
/// filter NaN before adding.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Key(f64);

impl Eq for Key {}

impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Key {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.0
            .partial_cmp(&other.0)
            .expect("NaN samples must be filtered before adding to a Histogram")
    }
}

/// A collection of samples, sorted into a fixed number of buckets.
#[derive(Debug, Clone)]
pub struct Histogram {
    num_buckets: u64,
    samples: BTreeMap<Key, u64>,
    min: Option<f64>,
    max: Option<f64>,
    // Welford accumulators, matching streaming-stats' OnlineStats term for term
    // so that the printed mean/stddev/variance keep their exact f64 rounding.
    size: u64,
    mean: f64,
    variance: f64,
}

impl Histogram {
    /// Construct a histogram with the given number of buckets.
    ///
    /// # Panics
    ///
    /// Panics if `num_buckets` is zero.
    pub fn with_buckets(num_buckets: u64) -> Histogram {
        assert!(num_buckets > 0, "a histogram needs at least one bucket");
        Histogram {
            num_buckets,
            samples: BTreeMap::new(),
            min: None,
            max: None,
            size: 0,
            mean: 0.0,
            variance: 0.0,
        }
    }

    /// Add a sample. NaN must be filtered out by the caller.
    pub fn add(&mut self, sample: f64) {
        *self.samples.entry(Key(sample)).or_insert(0) += 1;

        self.min = Some(
            self.min
                .map_or(sample, |m| if sample < m { sample } else { m }),
        );
        self.max = Some(
            self.max
                .map_or(sample, |m| if sample > m { sample } else { m }),
        );

        let oldmean = self.mean;
        let prevq = self.variance * self.size as f64;
        self.size += 1;
        self.mean += (sample - oldmean) / self.size as f64;
        self.variance = (prevq + (sample - oldmean) * (sample - self.mean)) / self.size as f64;
    }

    fn buckets(&self) -> Buckets<'_> {
        Buckets { histogram: self, index: 0 }
    }
}

/// One bucket: a half-open range and the number of samples in it.
#[derive(Debug, Clone, Copy)]
struct Bucket {
    start: f64,
    end: f64,
    count: u64,
}

struct Buckets<'a> {
    histogram: &'a Histogram,
    index: u64,
}

impl Iterator for Buckets<'_> {
    type Item = Bucket;

    fn next(&mut self) -> Option<Bucket> {
        if self.index >= self.histogram.num_buckets {
            return None;
        }
        let (min, max) = match (self.histogram.min, self.histogram.max) {
            (Some(min), Some(max)) => (min, max),
            _ => return None,
        };

        let range = max - min;
        // Upstream quirk, kept for identical bucket edges: widen the range by
        // its integer remainder over the bucket count before dividing.
        let range = range + (range as u64 % self.histogram.num_buckets) as f64;
        let bucket_size = range / self.histogram.num_buckets as f64;

        let start = min + self.index as f64 * bucket_size;
        let end = min + (self.index + 1) as f64 * bucket_size;

        self.index += 1;

        // The final bucket is unbounded above so that the maximum sample lands
        // in it rather than falling off the end of the half-open range.
        let count = if self.index == self.histogram.num_buckets {
            self.histogram
                .samples
                .range(Key(start)..)
                .map(|(_, c)| c)
                .sum()
        } else {
            self.histogram
                .samples
                .range(Key(start)..Key(end))
                .map(|(_, c)| c)
                .sum()
        };

        Some(Bucket { start, end, count })
    }
}

impl fmt::Display for Histogram {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let num_samples: u64 = self.samples.values().sum();
        writeln!(f, "# Number of samples = {num_samples}")?;
        if num_samples == 0 {
            return Ok(());
        }

        let min = self.min.expect("min is set whenever samples exist");
        let max = self.max.expect("max is set whenever samples exist");
        writeln!(f, "# Min = {min}")?;
        writeln!(f, "# Max = {max}")?;
        writeln!(f, "#")?;
        writeln!(f, "# Mean = {}", self.mean)?;
        writeln!(f, "# Standard deviation = {}", self.variance.sqrt())?;
        writeln!(f, "# Variance = {}", self.variance)?;
        writeln!(f, "#")?;

        const WIDTH: u64 = 50;
        let max_bucket_count = self.buckets().map(|b| b.count).fold(0, cmp::max);
        let count_per_char = cmp::max(max_bucket_count / WIDTH, 1);

        writeln!(f, "# Each ∎ is a count of {count_per_char}")?;
        writeln!(f, "#")?;

        let widest_count = self
            .buckets()
            .fold(0, |n, b| cmp::max(n, b.count.to_string().len()));
        let widest_range = self
            .buckets()
            .fold(0, |n, b| cmp::max(n, b.end.to_string().len()));

        for bucket in self.buckets() {
            write!(
                f,
                "{start:>range_width$} .. {end:>range_width$} [ {count:>count_width$} ]: ",
                start = bucket.start,
                end = bucket.end,
                count = bucket.count,
                range_width = widest_range,
                count_width = widest_count,
            )?;
            for _ in 0..bucket.count / count_per_char {
                write!(f, "∎")?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn render(num_buckets: u64, samples: &[f64]) -> String {
        let mut h = Histogram::with_buckets(num_buckets);
        for &s in samples {
            h.add(s);
        }
        h.to_string()
    }

    /// Pins the exact format: bucket edges, column padding and bar scaling.
    ///
    /// This is the output histo_fp 0.2.1 produced for the same input, captured
    /// while a differential test ran both implementations side by side over
    /// ~900 cases. Note it does NOT match the example in histo_fp's own crate
    /// docs, which was never updated after the crate moved from integer to
    /// floating-point bucket edges.
    #[test]
    fn golden_output() {
        let samples: Vec<f64> = (0..100).flat_map(|i| [i as f64, (i * i) as f64]).collect();
        let expected = "\
# Number of samples = 200
# Min = 0
# Max = 9801
#
# Mean = 1666.5000000000005
# Standard deviation = 2641.2281518263426
# Variance = 6976086.1499999985
#
# Each ∎ is a count of 2
#
                 0 ..              980.2 [ 132 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
             980.2 ..             1960.4 [  13 ]: ∎∎∎∎∎∎
            1960.4 .. 2940.6000000000004 [  10 ]: ∎∎∎∎∎
2940.6000000000004 ..             3920.8 [   8 ]: ∎∎∎∎
            3920.8 ..               4901 [   8 ]: ∎∎∎∎
              4901 ..  5881.200000000001 [   6 ]: ∎∎∎
 5881.200000000001 ..  6861.400000000001 [   6 ]: ∎∎∎
 6861.400000000001 ..             7841.6 [   6 ]: ∎∎∎
            7841.6 ..  8821.800000000001 [   5 ]: ∎∎
 8821.800000000001 ..               9802 [   6 ]: ∎∎∎
";
        assert_eq!(render(10, &samples), expected);
    }

    #[test]
    fn empty_histogram_prints_only_the_count() {
        assert_eq!(render(100, &[]), "# Number of samples = 0\n");
    }

    #[test]
    fn single_sample_lands_in_the_last_bucket() {
        let out = render(4, &[7.5]);
        assert!(out.contains("# Number of samples = 1"), "{out}");
        assert!(out.contains("# Min = 7.5"), "{out}");
        assert!(out.contains("# Max = 7.5"), "{out}");
        // zero range: every bucket edge is the sample itself, and the
        // unbounded final bucket collects it
        assert!(out.contains("7.5 .. 7.5 [ 1 ]: ∎"), "{out}");
    }

    #[test]
    fn identical_samples_have_zero_variance() {
        let out = render(10, &[3.0; 25]);
        assert!(out.contains("# Mean = 3"), "{out}");
        assert!(out.contains("# Standard deviation = 0"), "{out}");
        assert!(out.contains("# Variance = 0"), "{out}");
    }

    #[test]
    fn every_sample_is_counted_exactly_once() {
        for &n in &[1u64, 2, 7, 100] {
            let samples: Vec<f64> = (0..250).map(|i| f64::from(i) * 0.37 - 20.0).collect();
            let mut h = Histogram::with_buckets(n);
            for &s in &samples {
                h.add(s);
            }
            let total: u64 = h.buckets().map(|b| b.count).sum();
            assert_eq!(total, samples.len() as u64, "with {n} buckets");
        }
    }

    /// histo_fp 0.2.1 panicked with "attempt to subtract with overflow" on this
    /// shape (kmer frequencies in [0, 1] at inquiSTR's 100 buckets) because it
    /// padded start values to the width of the widest *end* value and underflowed
    /// when a start was wider. Here the width saturates instead.
    #[test]
    fn wide_start_values_do_not_underflow_the_padding() {
        let mut state: u64 = 42;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let samples: Vec<f64> = (0..500).map(|_| next()).collect();
        let out = render(100, &samples);
        assert_eq!(out.lines().count(), 110, "8 header lines + 100 buckets + trailing");
        assert!(out.contains("# Number of samples = 500"), "{out}");
    }

    #[test]
    fn negative_values_render() {
        let samples: Vec<f64> = (-50..50).map(f64::from).collect();
        let out = render(5, &samples);
        assert!(out.contains("# Min = -50"), "{out}");
        assert!(out.contains("# Max = 49"), "{out}");
    }
}
