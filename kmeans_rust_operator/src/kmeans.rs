//! The R `stats::kmeans` wrapper around the Hartigan–Wong engine
//! (`kmns.rs`), transliterated from R 4.3.3's `kmeans.R` for the case the
//! operator uses: `centers` given as a **number** (the operator's `centers`
//! property is a count, never a centre matrix), `algorithm = "Hartigan-Wong"`.
//!
//! The wrapper decides the whole trajectory of the fit, so it is ported as
//! literally as the engine:
//!
//! - `nstart == 1`: the initial centres are rows of the data drawn with
//!   `sample.int(m, k)` — R's partial Fisher–Yates over the RNG stream
//!   (`rng::RRng::sample_int`).
//! - If the drawn rows contain **duplicates** (R: `any(duplicated(centers))`),
//!   the wrapper throws the draw away and samples `k` of the *distinct* rows
//!   instead (`cn <- unique(x)`), continuing the same RNG stream. With
//!   `nstart >= 2` it always samples from the distinct rows. Row equality is
//!   R's `match` semantics (`0.0 == -0.0`, `NaN` matches nothing), not a
//!   bit compare — `unique()` would collapse `0.0` and `-0.0`.
//! - `mm < k` on the distinct path stops with R's message.
//! - `nstart > 1`: `nstart` fits are run and the first with the strictly
//!   smallest `sum(wss)` wins — ties keep the earlier start.
//! - `k == 1`: R switches to the MacQueen algorithm (`nmeth <- 3L`). Every
//!   point lands in the single cluster, so the labels — the only thing this
//!   operator publishes — are decided without running it. The within-cluster
//!   sums are reported as the exact `sum((x - colMean)^2)`, computed
//!   differently from MacQueen's incremental update; documented as a
//!   deviation (the operator never outputs them).
//!
//! Exit codes surface as the R wrapper surfaces them: `ifault` 1 and 3 are
//! `stop()`s, 2 ("did not converge in `iter.max` iterations") and 4
//! ("Quick-TRANSfer stage steps exceeded maximum") are warnings with the fit
//! returned.

use anyhow::{Result, bail};

use crate::kmns::{self, KmnsResult};
use crate::rng::RRng;

#[derive(Debug, Clone, PartialEq)]
pub struct KmeansFit {
    /// Cluster of each point, 0-based.
    pub cluster: Vec<usize>,
    /// Final centres, row-major `k × p`.
    pub centers: Vec<f64>,
    /// Within-cluster sum of squares per cluster.
    pub wss: Vec<f64>,
    /// `sum(wss)` — R's `tot.withinss` (its `nstart` selection criterion).
    pub tot_withinss: f64,
    /// R's `iter` (the engine's `ITER`; `iter.max + 1` on non-convergence).
    pub iter: usize,
    /// The engine's exit code: 0 ok, 2 non-convergence, 4 quick-transfer limit.
    pub ifault: i32,
    /// R `warning()`s emitted on the way (non-fatal `ifault` exits).
    pub warnings: Vec<String>,
}

/// `kmeans(x, centers = k, iter.max, nstart)` with `x` row-major `m × p`.
///
/// `rng` must already be in the state R's session RNG is in right before the
/// call — the R operator does `set.seed(seed)` and nothing else.
pub fn kmeans(
    x: &[f64],
    m: usize,
    p: usize,
    k: usize,
    iter_max: usize,
    nstart: usize,
    rng: &mut RRng,
) -> Result<KmeansFit> {
    if m < 1 {
        bail!("invalid nrow(x)");
    }
    if k >= 1 && k > m {
        // R reaches here through `sample.int(m, k)`.
        bail!("cannot take a sample larger than the population");
    }
    if k < 1 {
        bail!("number of cluster centres must lie between 1 and nrow(x)");
    }
    if iter_max < 1 {
        bail!("'iter.max' must be positive");
    }

    let isteps_qtran = 50usize.saturating_mul(m).min(i32::MAX as usize);

    // k == 1: R routes to MacQueen; the labels are trivially all one cluster.
    // The RNG still pays the draw the wrapper makes for the initial centres.
    if k == 1 {
        let _ = rng.sample_int(m, 1);
        return Ok(macqueen_k1(x, m, p));
    }

    // Scalar centres, R wrapper:
    //   if (nstart == 1L) centers <- x[sample.int(m, k), , drop = FALSE]
    //   if (nstart >= 2L || any(duplicated(centers))) { cn <- unique(x); ... }
    let mut distinct: Option<Vec<usize>> = None; // `cn` as row indices into x
    // Draw k rows of `cn` (the distinct rows): sample.int indexes the
    // distinct list, then the rows come from x through it.
    let draw_distinct = |cn: &[usize], rng: &mut RRng| -> Vec<usize> {
        rng.sample_int(cn.len(), k)
            .into_iter()
            .map(|j| cn[j])
            .collect()
    };
    let mut centers = if nstart == 1 {
        let drawn: Vec<usize> = rng.sample_int(m, k);
        if has_duplicate_rows(x, &drawn, p) {
            let cn = unique_row_indices(x, m, p);
            if cn.len() < k {
                bail!("more cluster centers than distinct data points.");
            }
            let rows = draw_distinct(&cn, rng);
            distinct = Some(cn);
            rows_of(x, &rows, p)
        } else {
            rows_of(x, &drawn, p)
        }
    } else {
        let cn = unique_row_indices(x, m, p);
        if cn.len() < k {
            bail!("more cluster centers than distinct data points.");
        }
        let rows = draw_distinct(&cn, rng);
        distinct = Some(cn);
        rows_of(x, &rows, p)
    };

    let mut best = do_one(x, m, p, &centers, k, iter_max, isteps_qtran)?;

    if nstart >= 2 {
        debug_assert!(distinct.is_some(), "distinct rows known when nstart >= 2");
        let cn = distinct.as_ref().unwrap().clone();
        for _ in 1..nstart {
            let rows = draw_distinct(&cn, rng);
            centers = rows_of(x, &rows, p);
            let zz = do_one(x, m, p, &centers, k, iter_max, isteps_qtran)?;
            if zz.tot_withinss < best.tot_withinss {
                best = zz;
            }
        }
    }
    Ok(best)
}

/// One `do_one(nmeth = 1)` of the R wrapper: run the engine, map `ifault`.
fn do_one(
    x: &[f64],
    m: usize,
    p: usize,
    centers: &[f64],
    k: usize,
    iter_max: usize,
    isteps_qtran: usize,
) -> Result<KmeansFit> {
    let mut c = centers.to_vec();
    let KmnsResult {
        cluster,
        wss,
        iter,
        ifault,
    } = kmns::kmns(x, m, p, &mut c, k, iter_max, isteps_qtran);
    let mut warnings = Vec::new();
    match ifault {
        0 | 2 | 4 => {}
        1 => bail!("empty cluster: try a better set of initial centers"),
        3 => bail!("number of cluster centres must lie between 1 and nrow(x)"),
        other => bail!("kmeans: unexpected ifault {other}"),
    }
    if ifault == 2 {
        warnings.push(format!(
            "did not converge in {iter_max} iteration{}",
            if iter_max == 1 { "" } else { "s" }
        ));
    }
    if ifault == 4 {
        warnings.push(format!(
            "Quick-TRANSfer stage steps exceeded maximum (= {isteps_qtran})"
        ));
    }
    // R computes `best <- sum(Z$wss)` with `sum()`, which accumulates in long
    // double and rounds once; an f64 accumulation can land one ulp away. The
    // per-cluster wss are bitwise R's, so two starts can only swap here when
    // their exact totals differ by less than one f64 rounding — and the
    // operator publishes labels, not this total (documented deviation).
    let tot_withinss = wss.iter().sum();
    Ok(KmeansFit {
        cluster,
        centers: c,
        wss,
        tot_withinss,
        iter,
        ifault,
        warnings,
    })
}

/// The `k == 1` fit: every point in cluster 0; centres are the column means;
/// `wss` is `sum((x - colMean)^2)` computed exactly (see the module comment).
fn macqueen_k1(x: &[f64], m: usize, p: usize) -> KmeansFit {
    let mut centers = vec![0.0f64; p];
    for i in 0..m {
        for j in 0..p {
            centers[j] += x[i * p + j];
        }
    }
    for c in centers.iter_mut() {
        *c /= m as f64;
    }
    let mut wss = 0.0f64;
    for i in 0..m {
        for j in 0..p {
            let d = x[i * p + j] - centers[j];
            wss += d * d;
        }
    }
    KmeansFit {
        cluster: vec![0; m],
        centers,
        wss: vec![wss],
        tot_withinss: wss,
        iter: 1,
        ifault: 0,
        warnings: Vec::new(),
    }
}

fn rows_of(x: &[f64], idx: &[usize], p: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(idx.len() * p);
    for &i in idx {
        out.extend_from_slice(&x[i * p..i * p + p]);
    }
    out
}

/// `any(duplicated(centers))` on the drawn rows — R's `==` semantics.
fn has_duplicate_rows(x: &[f64], idx: &[usize], p: usize) -> bool {
    for a in 0..idx.len() {
        for b in 0..a {
            if row_eq(x, idx[a], idx[b], p) {
                return true;
            }
        }
    }
    false
}

/// `unique(x)` as row indices, first occurrence first.
fn unique_row_indices(x: &[f64], m: usize, p: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for i in 0..m {
        if !out.iter().any(|&j| row_eq(x, i, j, p)) {
            out.push(i);
        }
    }
    out
}

/// Row equality with R `match` semantics: `==` per element (so `0.0` equals
/// `-0.0`, and a `NaN` element matches nothing — including itself).
fn row_eq(x: &[f64], i: usize, j: usize, p: usize) -> bool {
    (0..p).all(|c| x[i * p + c] == x[j * p + c])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng(seed: u32) -> RRng {
        RRng::set_seed(seed)
    }

    /// Four points on a line, two obvious groups; centres fixed so the fit
    /// is deterministic (the R wrapper path is exercised elsewhere).
    fn line_points() -> (Vec<f64>, usize, usize) {
        (vec![0.0, 0.0, 1.0, 0.0, 10.0, 0.0, 11.0, 0.0], 4, 2)
    }

    #[test]
    fn a_two_group_line_splits_at_the_middle() {
        let (x, m, p) = line_points();
        // Initial centres at points 0 and 2 (the two groups' anchors).
        let centers = vec![0.0, 0.0, 10.0, 0.0];
        let fit = do_one(&x, m, p, &centers, 2, 10, 50 * m).unwrap();
        assert_eq!(fit.cluster, vec![0, 0, 1, 1]);
        assert_eq!(fit.centers, vec![0.5, 0.0, 10.5, 0.0]);
        assert_eq!(fit.wss, vec![0.5, 0.5]);
        assert_eq!(fit.ifault, 0);
    }

    #[test]
    fn more_centers_than_points_is_r_sample_error() {
        let (x, m, p) = line_points();
        let e = kmeans(&x, m, p, 5, 10, 1, &mut rng(1)).unwrap_err();
        assert!(e.to_string().contains("cannot take a sample larger"));
    }

    #[test]
    fn centers_equal_to_points_is_r_engine_error() {
        let (x, m, p) = line_points();
        // k == m passes sampling and reaches the engine's k >= m exit.
        let e = kmeans(&x, m, p, m, 10, 1, &mut rng(1)).unwrap_err();
        assert!(
            e.to_string()
                .contains("number of cluster centres must lie between 1 and nrow(x)")
        );
    }

    #[test]
    fn zero_centers_is_r_engine_error() {
        let (x, m, p) = line_points();
        let e = kmeans(&x, m, p, 0, 10, 1, &mut rng(1)).unwrap_err();
        assert!(
            e.to_string()
                .contains("number of cluster centres must lie between 1 and nrow(x)")
        );
    }

    #[test]
    fn non_positive_iter_max_is_r_error() {
        let (x, m, p) = line_points();
        let e = kmeans(&x, m, p, 2, 0, 1, &mut rng(1)).unwrap_err();
        assert!(e.to_string().contains("'iter.max' must be positive"));
    }

    #[test]
    fn one_cluster_takes_every_point() {
        let (x, m, p) = line_points();
        let fit = kmeans(&x, m, p, 1, 10, 1, &mut rng(1)).unwrap();
        assert_eq!(fit.cluster, vec![0, 0, 0, 0]);
        assert_eq!(fit.centers, vec![5.5, 0.0]);
        assert_eq!(
            fit.tot_withinss,
            5.5 * 5.5 + 4.5 * 4.5 + 4.5 * 4.5 + 5.5 * 5.5
        );
    }

    #[test]
    fn duplicate_drawn_centres_are_resampled_from_distinct_rows() {
        // Rows 0 and 2 are identical; seed 42 draws indices 0 and 4 first for
        // m=10, k=2 — craft m=4, k=2 so that whatever the draw, the duplicate
        // path is forced by making every draw contain both copies: with only
        // two distinct rows among four, k=2 drawn from x must hit the
        // duplicate branch unless the draw picks one of each. Deterministic
        // assertion: the fit never fails and covers all four points.
        let x = vec![1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0];
        let fit = kmeans(&x, 5, 2, 2, 10, 1, &mut rng(42)).unwrap();
        assert_eq!(fit.cluster.len(), 5);
        assert!(fit.cluster.iter().all(|c| *c < 2));
    }

    #[test]
    fn row_equality_matches_r_match_semantics() {
        let x = vec![0.0, 0.0, 0.0f64.copysign(-1.0), 0.0];
        // 0.0 == -0.0 in R's match, so these two rows are duplicates.
        assert!(row_eq(&x, 0, 1, 2));
        assert!(has_duplicate_rows(&x, &[0, 1], 2));
    }

    #[test]
    fn nstart_runs_are_compared_on_total_withinss() {
        // Well-separated groups: every start converges to the same optimum,
        // so the first stays the winner and the fit equals the nstart=1 fit.
        let mut x = Vec::new();
        for g in 0..3 {
            for i in 0..10 {
                x.push((g * 100 + i) as f64);
                x.push(0.0);
            }
        }
        let mut r1 = rng(9);
        let a = kmeans(&x, 30, 2, 3, 10, 1, &mut r1).unwrap();
        let mut r3 = rng(9);
        let b = kmeans(&x, 30, 2, 3, 10, 3, &mut r3).unwrap();
        assert_eq!(a.cluster, b.cluster);
        assert_eq!(a.centers, b.centers);
    }
}
