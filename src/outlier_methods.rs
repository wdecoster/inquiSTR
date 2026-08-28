//! Scoring methods for the outlier scan.
//!
//! Every method answers the same question — how extreme is this allele, given a reference
//! distribution — and returns a score, a significance, and typed evidence. None of them
//! returns a boolean: deciding what counts as a call is a downstream, user-configurable step,
//! and is enforced here by [`OutlierResult`] having no such field.
//!
//! Scores are expressed in **repeat units** (base pairs divided by motif length) so that they
//! are comparable across loci. Ranking a whole genome depends on that: a score that is only
//! meaningful within one locus produces a meaningless shortlist.

use clap::ValueEnum;

/// Available scoring methods
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Method {
    /// How many reference alleles are longer. Makes no assumption about distribution shape.
    Percentile,
    /// Distance from the locus median, scaled by median absolute deviation
    Robustz,
    /// Distance from the locus mean, scaled by standard deviation
    Zscore,
    /// Density-based clustering; outliers are points belonging to no cluster
    Dbscan,
}

impl Method {
    pub fn as_str(&self) -> &'static str {
        match self {
            Method::Percentile => "percentile",
            Method::Robustz => "robustz",
            Method::Zscore => "zscore",
            Method::Dbscan => "dbscan",
        }
    }

    /// Whether leave-one-out exclusion changes what this method sees.
    ///
    /// Mean and standard deviation are moved by the very allele under test, so excluding it
    /// matters. Median and MAD are not — resisting a single extreme value is what makes them
    /// robust — and density clustering classifies all points jointly.
    pub fn uses_leave_one_out(&self) -> bool {
        matches!(self, Method::Percentile | Method::Zscore)
    }
}

/// Method-specific diagnostics, kept typed all the way to the output writer
#[derive(Clone, Debug, PartialEq)]
pub enum Evidence {
    /// Empirical tail position
    Exceedance {
        /// Reference alleles considered, after any exclusions
        ref_n: usize,
        /// Reference alleles strictly longer than the query
        greater: usize,
        /// Reference alleles at least as long as the query
        at_least: usize,
        /// Longest reference allele
        ref_max: f64,
        /// Finest significance this reference size can resolve, 1/(ref_n + 1)
        alpha_floor: f64,
        /// True when the reference maximum equals the nominal high percentile, meaning a
        /// percentile threshold has degenerated into "longer than anything observed"
        saturated: bool,
    },
    /// Centre-and-scale statistics
    Moment {
        ref_n: usize,
        center: f64,
        scale: f64,
        robust: bool,
    },
    /// Density clustering assignment
    Cluster {
        ref_n: usize,
        eps: f64,
        min_points: usize,
        noise: bool,
    },
}

impl Evidence {
    /// Render as a compact `key=value;...` field. The typing above is what matters; this is
    /// only the serialisation at the boundary.
    pub fn render(&self) -> String {
        match self {
            Evidence::Exceedance { ref_n, greater, at_least, ref_max, alpha_floor, saturated } => {
                format!(
                    "ref_n={};greater={};at_least={};ref_max={};alpha_floor={:.3e}{}",
                    ref_n,
                    greater,
                    at_least,
                    ref_max,
                    alpha_floor,
                    if *saturated { ";saturated=true" } else { "" }
                )
            }
            Evidence::Moment { ref_n, center, scale, robust } => {
                format!("ref_n={};center={:.4};scale={:.4};robust={}", ref_n, center, scale, robust)
            }
            Evidence::Cluster { ref_n, eps, min_points, noise } => {
                format!("ref_n={};eps={};min_points={};noise={}", ref_n, eps, min_points, noise)
            }
        }
    }

    pub fn ref_n(&self) -> usize {
        match self {
            Evidence::Exceedance { ref_n, .. }
            | Evidence::Moment { ref_n, .. }
            | Evidence::Cluster { ref_n, .. } => *ref_n,
        }
    }
}

/// The outcome of scoring one allele against one reference distribution.
///
/// Deliberately has no boolean field. Whether a result is "an outlier" is a thresholding
/// decision that belongs to the caller and, ultimately, to the user.
#[derive(Clone, Debug)]
pub struct OutlierResult {
    /// Monotone in how extreme the allele is, in repeat units, comparable across loci
    pub score: f64,
    /// P-value, exceedance probability, or empirical rank, depending on method
    pub significance: f64,
    pub evidence: Evidence,
}

/// Reference alleles at one locus, plus the running sums needed for closed-form exclusion.
pub struct LocusReference<'a> {
    /// All called alleles at the locus, ascending
    pub sorted: &'a [f64],
    pub sum: f64,
    pub sumsq: f64,
    /// Motif length in bp; scores are divided by this to give repeat units
    pub motif_len: f64,
}

impl LocusReference<'_> {
    /// Largest reference allele once `excluded` values have been removed, one occurrence each.
    fn max_excluding(&self, excluded: &[f64]) -> f64 {
        let mut skip: Vec<f64> = excluded.to_vec();
        for &v in self.sorted.iter().rev() {
            if let Some(pos) = skip.iter().position(|&s| s == v) {
                skip.swap_remove(pos);
                continue;
            }
            return v;
        }
        f64::NAN
    }
}

/// Count of reference alleles strictly greater than, and at least as large as, `query`,
/// after removing one occurrence of each value in `excluded`.
fn tail_counts(sorted: &[f64], query: f64, excluded: &[f64]) -> (usize, usize) {
    let n = sorted.len();
    let greater = n - sorted.partition_point(|&x| x <= query);
    let at_least = n - sorted.partition_point(|&x| x < query);
    let own_greater = excluded.iter().filter(|&&v| v > query).count();
    let own_at_least = excluded.iter().filter(|&&v| v >= query).count();
    (greater.saturating_sub(own_greater), at_least.saturating_sub(own_at_least))
}

/// Score an allele by its position in the empirical tail.
///
/// `excluded` holds the alleles belonging to the sample under test, which are removed from the
/// reference so that a carrier cannot define the threshold they are measured against. Without
/// this, an individual holding the longest allele at a locus is never remarkable — which is
/// precisely the individual the scan exists to find.
pub fn score_exceedance(reference: &LocusReference, query: f64, excluded: &[f64]) -> OutlierResult {
    let ref_n = reference.sorted.len().saturating_sub(excluded.len());
    let (greater, at_least) = tail_counts(reference.sorted, query, excluded);
    let ref_max = reference.max_excluding(excluded);
    let alpha_floor = 1.0 / (ref_n as f64 + 1.0);
    let significance = (at_least as f64 + 1.0) / (ref_n as f64 + 1.0);
    let score = if ref_max.is_finite() {
        (query - ref_max) / reference.motif_len
    } else {
        f64::NAN
    };
    OutlierResult {
        score,
        significance,
        evidence: Evidence::Exceedance {
            ref_n,
            greater,
            at_least,
            ref_max,
            alpha_floor,
            saturated: false,
        },
    }
}

/// Standard normal upper-tail probability for a given z, used to turn a z-score cutoff into
/// the significance gate applied by the scan.
pub fn normal_sf_public(z: f64) -> f64 {
    normal_sf(z)
}

/// Standard normal upper-tail probability, via the complementary error function.
fn normal_sf(z: f64) -> f64 {
    // Infinite z is a real result, not bad input: it is what a zero scale estimate produces,
    // and the tail probability there is exactly zero (or one). Only NaN is undefined.
    if z.is_nan() {
        return f64::NAN;
    }
    if z == f64::INFINITY {
        return 0.0;
    }
    if z == f64::NEG_INFINITY {
        return 1.0;
    }
    0.5 * erfc(z / std::f64::consts::SQRT_2)
}

/// Complementary error function (Numerical Recipes rational approximation, ~1e-7 accuracy)
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 2.0 / (2.0 + z);
    let ty = 4.0 * t - 2.0;
    const COF: [f64; 28] = [
        -1.3026537197817094,
        6.419_697_923_564_902e-1,
        1.9476473204185836e-2,
        -9.561_514_786_808_63e-3,
        -9.46595344482036e-4,
        3.66839497852761e-4,
        4.2523324806907e-5,
        -2.0278578112534e-5,
        -1.624290004647e-6,
        1.303655835580e-6,
        1.5626441722e-8,
        -8.5238095915e-8,
        6.529054439e-9,
        5.059343495e-9,
        -9.91364156e-10,
        -2.27365122e-10,
        9.6467911e-11,
        2.394038e-12,
        -6.886027e-12,
        8.94487e-13,
        3.13092e-13,
        -1.12708e-13,
        3.81e-16,
        7.106e-15,
        -1.523e-15,
        -9.4e-17,
        1.21e-16,
        -2.8e-17,
    ];
    let (mut d, mut dd) = (0.0f64, 0.0f64);
    for &c in COF.iter().rev().take(COF.len() - 1) {
        let tmp = d;
        d = ty * d - dd + c;
        dd = tmp;
    }
    let ans = t * (-z * z + 0.5 * (COF[0] + ty * d) - dd).exp();
    if x >= 0.0 { ans } else { 2.0 - ans }
}

/// Score an allele by distance from a centre, in units of a scale estimate.
///
/// With `robust`, centre and scale are the median and median absolute deviation; otherwise the
/// mean and standard deviation, with the tested sample's own alleles removed in closed form
/// from the running sums.
///
/// Both carry the same caveat, which belongs in user-facing documentation: dividing by a
/// measure of locus spread is dividing by the quantity the variability index rewards, so
/// sensitivity is worst at exactly the high-variability loci prioritisation is meant to
/// surface.
pub fn score_moment(
    reference: &LocusReference,
    query: f64,
    excluded: &[f64],
    robust: bool,
    median: f64,
    mad: f64,
) -> OutlierResult {
    let (center, scale, ref_n) = if robust {
        // Median and MAD are unmoved by a single extreme value, so no exclusion is applied.
        (median, mad, reference.sorted.len())
    } else {
        let n = reference.sorted.len().saturating_sub(excluded.len()) as f64;
        let s = reference.sum - excluded.iter().sum::<f64>();
        let q = reference.sumsq - excluded.iter().map(|v| v * v).sum::<f64>();
        if n <= 1.0 {
            (f64::NAN, f64::NAN, 0)
        } else {
            let mean = s / n;
            let var = (q / n - mean * mean).max(0.0);
            (mean, var.sqrt(), n as usize)
        }
    };

    // 1.4826 rescales MAD to be comparable with a standard deviation under normality
    let scaled = if robust { scale * 1.4826 } else { scale };
    // A zero scale estimate is common rather than exceptional - MAD is exactly zero at most
    // loci, because more than half the cohort usually carries the same allele - so returning
    // NaN there discards the cleanest signals in the data. Treating it as infinite
    // significance is the opposite error: every one-base deviation would then qualify.
    // Instead the scale falls back to one repeat unit, the finest difference that means
    // anything at a tandem repeat. The result stays finite and proportional to magnitude.
    let effective_scale = if scaled > 0.0 {
        scaled
    } else {
        reference.motif_len.max(1.0)
    };
    let z = (query - center) / effective_scale;
    OutlierResult {
        score: if center.is_finite() {
            (query - center) / reference.motif_len
        } else {
            f64::NAN
        },
        significance: normal_sf(z),
        evidence: Evidence::Moment { ref_n, center, scale, robust },
    }
}

/// Build a result for a point classified by density clustering.
pub fn dbscan_result(
    reference: &LocusReference,
    query: f64,
    eps: f64,
    min_points: usize,
    noise: bool,
    nearest_core: f64,
) -> OutlierResult {
    OutlierResult {
        score: if nearest_core.is_finite() {
            (query - nearest_core) / reference.motif_len
        } else {
            f64::NAN
        },
        significance: f64::NAN,
        evidence: Evidence::Cluster { ref_n: reference.sorted.len(), eps, min_points, noise },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn refr(v: &[f64]) -> (Vec<f64>, f64, f64) {
        let mut s = v.to_vec();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let sum = s.iter().sum();
        let sumsq = s.iter().map(|x| x * x).sum();
        (s, sum, sumsq)
    }

    #[test]
    fn test_carrier_does_not_define_its_own_threshold() {
        // One sample holds the longest allele by a wide margin. Without exclusion it is tied
        // with the reference maximum and looks unremarkable; with exclusion it stands out.
        let (sorted, sum, sumsq) = refr(&[10.0, 10.0, 11.0, 12.0, 12.0, 90.0]);
        let r = LocusReference { sorted: &sorted, sum, sumsq, motif_len: 3.0 };

        let without = score_exceedance(&r, 90.0, &[]);
        assert_eq!(without.score, 0.0, "self-inclusion makes the carrier look ordinary");

        let with = score_exceedance(&r, 90.0, &[90.0]);
        assert!(with.score > 0.0, "exclusion must reveal the carrier");
        assert_eq!(with.score, (90.0 - 12.0) / 3.0);
        if let Evidence::Exceedance { greater, at_least, ref_max, ref_n, .. } = with.evidence {
            assert_eq!(greater, 0);
            assert_eq!(at_least, 0);
            assert_eq!(ref_max, 12.0);
            assert_eq!(ref_n, 5);
        } else {
            panic!("wrong evidence variant");
        }
    }

    #[test]
    fn test_both_alleles_of_the_tested_sample_are_excluded() {
        let (sorted, sum, sumsq) = refr(&[5.0, 5.0, 6.0, 40.0, 41.0]);
        let r = LocusReference { sorted: &sorted, sum, sumsq, motif_len: 1.0 };
        // A sibling pair effect in miniature: the sample carries both long alleles
        let res = score_exceedance(&r, 41.0, &[40.0, 41.0]);
        if let Evidence::Exceedance { ref_max, ref_n, .. } = res.evidence {
            assert_eq!(ref_max, 6.0, "the sample's other long allele must not remain");
            assert_eq!(ref_n, 3);
        } else {
            panic!("wrong evidence variant");
        }
    }

    #[test]
    fn test_monomorphic_locus_yields_no_margin() {
        let (sorted, sum, sumsq) = refr(&[20.0; 8]);
        let r = LocusReference { sorted: &sorted, sum, sumsq, motif_len: 4.0 };
        let res = score_exceedance(&r, 20.0, &[20.0, 20.0]);
        assert_eq!(res.score, 0.0, "everyone tied at the top must score zero, not top rank");
    }

    #[test]
    fn test_alpha_floor_reports_resolution_limit() {
        let (sorted, sum, sumsq) = refr(&(0..100).map(|i| i as f64).collect::<Vec<_>>());
        let r = LocusReference { sorted: &sorted, sum, sumsq, motif_len: 1.0 };
        let res = score_exceedance(&r, 500.0, &[]);
        if let Evidence::Exceedance { alpha_floor, ref_n, .. } = res.evidence {
            assert_eq!(ref_n, 100);
            assert!((alpha_floor - 1.0 / 101.0).abs() < 1e-12);
            // Nothing can be more significant than the floor, however extreme it is
            assert!((res.significance - alpha_floor).abs() < 1e-12);
        } else {
            panic!("wrong evidence variant");
        }
    }

    #[test]
    fn test_score_is_in_repeat_units() {
        // The same 30 bp excess is 10 units on a trinucleotide and 6 on a pentanucleotide,
        // which is what makes scores comparable between loci.
        let (sorted, sum, sumsq) = refr(&[10.0, 10.0, 10.0]);
        let tri = LocusReference { sorted: &sorted, sum, sumsq, motif_len: 3.0 };
        let pen = LocusReference { sorted: &sorted, sum, sumsq, motif_len: 5.0 };
        assert_eq!(score_exceedance(&tri, 40.0, &[]).score, 10.0);
        assert_eq!(score_exceedance(&pen, 40.0, &[]).score, 6.0);
    }

    #[test]
    fn test_zscore_excludes_query_from_its_own_mean() {
        let (sorted, sum, sumsq) = refr(&[10.0, 10.0, 10.0, 10.0, 100.0]);
        let r = LocusReference { sorted: &sorted, sum, sumsq, motif_len: 1.0 };
        let res = score_moment(&r, 100.0, &[100.0], false, 10.0, 0.0);
        if let Evidence::Moment { center, scale, ref_n, .. } = res.evidence {
            assert_eq!(ref_n, 4);
            assert!((center - 10.0).abs() < 1e-9, "the outlier must not inflate its own mean");
            assert!(scale.abs() < 1e-9);
        } else {
            panic!("wrong evidence variant");
        }
    }

    #[test]
    fn test_normal_sf_matches_known_values() {
        assert!((normal_sf(0.0) - 0.5).abs() < 1e-6);
        assert!((normal_sf(1.0) - 0.158655).abs() < 1e-5);
        assert!((normal_sf(1.96) - 0.025).abs() < 1e-4);
        assert!((normal_sf(3.0) - 0.001349).abs() < 1e-6);
    }

    #[test]
    fn test_result_carries_no_verdict() {
        // A compile-time guarantee in spirit: thresholding lives downstream, so the result
        // exposes only continuous quantities and diagnostics.
        let (sorted, sum, sumsq) = refr(&[1.0, 2.0, 3.0]);
        let r = LocusReference { sorted: &sorted, sum, sumsq, motif_len: 1.0 };
        let res = score_exceedance(&r, 99.0, &[]);
        assert!(res.score.is_finite() && res.significance.is_finite());
    }
}

#[cfg(test)]
mod scale_zero_tests {
    use super::*;

    #[test]
    fn test_zero_scale_is_infinitely_significant() {
        assert_eq!(normal_sf(f64::INFINITY), 0.0, "erfc must saturate to zero at infinity");
        let sorted: Vec<f64> = std::iter::repeat_n(30.0, 19).chain([330.0]).collect();
        let sum: f64 = sorted.iter().sum();
        let sumsq: f64 = sorted.iter().map(|v| v * v).sum();
        let r = LocusReference { sorted: &sorted, sum, sumsq, motif_len: 3.0 };
        let res = score_moment(&r, 330.0, &[330.0, 30.0], false, 30.0, 0.0);
        if let Evidence::Moment { center, scale, ref_n, .. } = res.evidence {
            assert_eq!(ref_n, 18);
            assert_eq!(center, 30.0);
            assert_eq!(scale, 0.0, "every remaining allele is identical");
        } else {
            panic!("wrong variant");
        }
        assert!(res.significance.is_finite(), "a zero scale must not yield NaN significance");
        assert!(res.significance < 1e-6, "a 300 bp excess over an invariant locus is extreme");
        assert_eq!(res.score, 100.0, "score must still be defined in repeat units");

        // The fallback must not make trivial deviations look significant
        let near = score_moment(&r, 33.0, &[33.0, 30.0], false, 30.0, 0.0);
        assert!(near.significance > 0.05, "a one-unit difference is not an outlier");
    }
}
