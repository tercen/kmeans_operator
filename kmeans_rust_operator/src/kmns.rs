//! `stats::kmeans`'s Hartigan–Wong engine, transliterated from R's Fortran.
//!
//! R runs the default algorithm through `KMNS` / `OPTRA` / `QTRAN`
//! (`src/library/stats/src/kmns.f`, a minimal double-precision modification
//! of Applied Statistics AS 136). A clustering has no tolerance — a tie broken
//! the other way is a different cluster for every point in the node — so this
//! is a statement-for-statement transliteration, including the parts that look
//! like mistakes:
//!
//! - the early-exit distance loops (`IF (DB .GE. DT(2)) GO TO 50`) **skip**
//!   the update check on exit, they do not fall through to it;
//! - every reduction adds the squared terms in the same order as the Fortran
//!   `DO J` loops, so the last bits of every centre are R's;
//! - `BIG` is `1e30`, not `f64::MAX`;
//! - `QTRAN`'s step counter and `NCP` bookkeeping keep R's 1-based convention.
//!
//! `ifault` carries R's exit codes: 0 ok, 1 empty cluster at the initial
//! assignment, 2 `iter.max` exceeded, 4 quick-transfer step limit exceeded.
//! (3 — k out of range — is decided by the caller, which mirrors the R
//! wrapper's own checks.)

pub const BIG: f64 = 1.0e30;

#[derive(Debug, Clone, PartialEq)]
pub struct KmnsResult {
    /// Cluster of each point, 0-based (Fortran `IC1 - 1`).
    pub cluster: Vec<usize>,
    /// Within-cluster sum of squares per cluster.
    pub wss: Vec<f64>,
    /// Fortran `ITER` after the run (`iter.max + 1` on the non-convergence exit).
    pub iter: usize,
    pub ifault: i32,
}

/// `KMNS`: divide `m` points (row-major `a`, `m × p`) into `k` clusters.
///
/// `centers` (row-major `k × p`) is the initial-centre matrix the R wrapper
/// samples; `iter_max` is `iter.max`; `isteps_qtran` is the R wrapper's
/// `min(.Machine$integer.max, 50 * m)`.
pub fn kmns(
    a: &[f64],
    m: usize,
    p: usize,
    centers: &mut [f64],
    k: usize,
    iter_max: usize,
    isteps_qtran: usize,
) -> KmnsResult {
    let mut ic1 = vec![0usize; m];
    let mut ic2 = vec![0usize; m];
    let mut nc = vec![0usize; k];
    let mut an1 = vec![0.0f64; k];
    let mut an2 = vec![0.0f64; k];
    let mut ncp = vec![0i32; k];
    let mut d = vec![0.0f64; m];
    // ITRAN is ITRAN(K+1) in the Fortran: entry 0 is iMaxQtr, 1..k the flags.
    let mut itr = vec![0i32; k + 1];
    let mut live = vec![0i32; k];
    let mut wss = vec![0.0f64; k];
    let mut dt = [0.0f64; 2];

    let mut ifault;

    // IFAULT = 3; IF (K .LE. 1 .OR. K .GE. M) RETURN — the caller never sends
    // k <= 1 (the R wrapper routes k == 1 to MacQueen) but k >= m is reachable.
    if k <= 1 || k >= m {
        return KmnsResult {
            cluster: ic1,
            wss,
            iter: 0,
            ifault: 3,
        };
    }
    ifault = 0;

    // For each point I, find its two closest centres, IC1(I) and IC2(I).
    for i in 0..m {
        ic1[i] = 0;
        ic2[i] = 1;
        for (il, dtil) in dt.iter_mut().enumerate() {
            *dtil = 0.0;
            for j in 0..p {
                let da = a[i * p + j] - centers[il * p + j];
                *dtil += da * da;
            }
        }
        if dt[0] > dt[1] {
            ic1[i] = 1;
            ic2[i] = 0;
            dt.swap(0, 1);
        }
        for l in 2..k {
            let mut db = 0.0f64;
            let mut early = false;
            for j in 0..p {
                let dc = a[i * p + j] - centers[l * p + j];
                db += dc * dc;
                if db >= dt[1] {
                    early = true;
                    break;
                }
            }
            if early {
                continue;
            }
            if db >= dt[0] {
                dt[1] = db;
                ic2[i] = l;
            } else {
                dt[1] = dt[0];
                ic2[i] = ic1[i];
                dt[0] = db;
                ic1[i] = l;
            }
        }
    }

    // Update cluster centres to be the average of points contained in them.
    for l in 0..k {
        nc[l] = 0;
        for j in 0..p {
            centers[l * p + j] = 0.0;
        }
    }
    for i in 0..m {
        let l = ic1[i];
        nc[l] += 1;
        for j in 0..p {
            centers[l * p + j] += a[i * p + j];
        }
    }
    for l in 0..k {
        if nc[l] == 0 {
            return KmnsResult {
                cluster: ic1,
                wss,
                iter: 0,
                ifault: 1,
            };
        }
        let aa = nc[l] as f64;
        for j in 0..p {
            centers[l * p + j] /= aa;
        }
        an2[l] = aa / (aa + 1.0);
        an1[l] = if aa > 1.0 { aa / (aa - 1.0) } else { BIG };
        itr[l + 1] = 1;
        ncp[l] = -1;
    }

    let mut indx: usize = 0;
    let mut ij: usize = 0;
    let mut stopped = false;
    for it in 1..=iter_max {
        ij = it;

        // Optimal-transfer stage: one pass over the data.
        optra(
            a, m, p, centers, k, &mut ic1, &mut ic2, &mut nc, &mut an1, &mut an2, &mut ncp, &mut d,
            &mut itr, &mut live, &mut indx,
        );

        // Stop if no transfer took place in the last M optimal transfer steps.
        if indx == m {
            stopped = true;
            break;
        }

        // Quick-transfer stage.
        let mut imaxqtr = isteps_qtran as i64;
        qtran(
            a,
            m,
            p,
            centers,
            &mut ic1,
            &mut ic2,
            &mut nc,
            &mut an1,
            &mut an2,
            &mut ncp,
            &mut d,
            &mut itr,
            &mut indx,
            &mut imaxqtr,
        );
        if imaxqtr < 0 {
            ifault = 4;
            stopped = true;
            break;
        }

        // If there are only two clusters there is no need to re-enter OPTRA.
        if k == 2 {
            stopped = true;
            break;
        }

        for v in ncp.iter_mut() {
            *v = 0;
        }
    }
    if !stopped {
        // Normal completion of `DO IJ = 1, ITER` leaves IJ = ITER + 1 in the
        // Fortran, which is exactly how R's "did not converge" run reports
        // one more iteration than allowed.
        ifault = 2;
        ij = iter_max + 1;
    }

    // Compute within-cluster sum of squares for each cluster (and the final
    // centres as their means), in the Fortran's loop order.
    for l in 0..k {
        wss[l] = 0.0;
        for j in 0..p {
            centers[l * p + j] = 0.0;
        }
    }
    for i in 0..m {
        let ii = ic1[i];
        for j in 0..p {
            centers[ii * p + j] += a[i * p + j];
        }
    }
    for j in 0..p {
        for l in 0..k {
            centers[l * p + j] /= nc[l] as f64;
        }
        for i in 0..m {
            let ii = ic1[i];
            let da = a[i * p + j] - centers[ii * p + j];
            wss[ii] += da * da;
        }
    }

    KmnsResult {
        cluster: ic1,
        wss,
        iter: ij,
        ifault,
    }
}

/// `OPTRA` — Algorithm AS 136.1: re-allocate each point to the cluster that
/// induces the maximum reduction in within-cluster sum of squares.
#[allow(clippy::too_many_arguments)]
fn optra(
    a: &[f64],
    m: usize,
    p: usize,
    c: &mut [f64],
    k: usize,
    ic1: &mut [usize],
    ic2: &mut [usize],
    nc: &mut [usize],
    an1: &mut [f64],
    an2: &mut [f64],
    ncp: &mut [i32],
    d: &mut [f64],
    itr: &mut [i32],
    live: &mut [i32],
    indx: &mut usize,
) {
    for l in 0..k {
        if itr[l + 1] == 1 {
            live[l] = m as i32 + 1;
        }
    }

    for i in 0..m {
        let i1 = i as i32 + 1; // Fortran's I, 1-based
        *indx += 1;
        let l1 = ic1[i];
        let mut l2 = ic2[i];
        let ll = l2;

        // If point I is the only member of cluster L1, no transfer.
        if nc[l1] == 1 {
            if *indx == m {
                return;
            }
            continue;
        }

        // If L1 has not yet been updated in this stage, no need to recompute D(I).
        if ncp[l1] != 0 {
            let mut de = 0.0f64;
            for j in 0..p {
                let df = a[i * p + j] - c[l1 * p + j];
                de += df * df;
            }
            d[i] = de * an1[l1];
        }

        // Find the cluster with minimum R2.
        let mut da = 0.0f64;
        for j in 0..p {
            let db = a[i * p + j] - c[l2 * p + j];
            da += db * db;
        }
        let mut r2 = da * an2[l2];
        for l in 0..k {
            // Only clusters in the live set can take the point, unless L1 is
            // still live itself; L1 and the second-best L2 are skipped.
            if (i1 >= live[l1] && i1 >= live[l]) || l == l1 || l == ll {
                continue;
            }
            let rr = r2 / an2[l];
            let mut dc = 0.0f64;
            let mut early = false;
            for j in 0..p {
                let dd = a[i * p + j] - c[l * p + j];
                dc += dd * dd;
                if dc >= rr {
                    early = true;
                    break;
                }
            }
            if early {
                continue;
            }
            r2 = dc * an2[l];
            l2 = l;
        }
        if r2 >= d[i] {
            // No transfer: L2 is the new IC2(I).
            ic2[i] = l2;
        } else {
            // Update centres, LIVE, NCP, AN1 & AN2 for clusters L1 and L2.
            *indx = 0;
            live[l1] = m as i32 + i1;
            live[l2] = m as i32 + i1;
            ncp[l1] = i1;
            ncp[l2] = i1;
            let al1 = nc[l1] as f64;
            let alw = al1 - 1.0;
            let al2 = nc[l2] as f64;
            let alt = al2 + 1.0;
            for j in 0..p {
                c[l1 * p + j] = (c[l1 * p + j] * al1 - a[i * p + j]) / alw;
                c[l2 * p + j] = (c[l2 * p + j] * al2 + a[i * p + j]) / alt;
            }
            nc[l1] -= 1;
            nc[l2] += 1;
            an2[l1] = alw / al1;
            an1[l1] = if alw > 1.0 { alw / (alw - 1.0) } else { BIG };
            an1[l2] = alt / al2;
            an2[l2] = alt / (alt + 1.0);
            ic1[i] = l2;
            ic2[i] = l1;
        }

        // If no re-allocation took place in the last M steps, return.
        if *indx == m {
            return;
        }
    }

    for l in 0..k {
        itr[l + 1] = 0; // before entering QTRAN
        live[l] -= m as i32; // LIVE has to be decreased by M before re-entering OPTRA
    }
}

/// `QTRAN` — Algorithm AS 136.2: keep testing each point in turn against its
/// second-best cluster, updating the centres after every move, until a full
/// pass over the data moves nothing.
#[allow(clippy::too_many_arguments)]
fn qtran(
    a: &[f64],
    m: usize,
    p: usize,
    c: &mut [f64],
    ic1: &mut [usize],
    ic2: &mut [usize],
    nc: &mut [usize],
    an1: &mut [f64],
    an2: &mut [f64],
    ncp: &mut [i32],
    d: &mut [f64],
    itr: &mut [i32],
    indx: &mut usize,
    imaxqtr: &mut i64,
) {
    let mut icoun: usize = 0;
    let mut istep: i64 = 0;
    loop {
        for i in 0..m {
            icoun += 1;
            istep += 1;
            if istep >= *imaxqtr {
                *imaxqtr = -1;
                return;
            }
            let l1 = ic1[i];
            let l2 = ic2[i];

            // If point I is the only member of cluster L1, no transfer.
            if nc[l1] == 1 {
                if icoun == m {
                    return;
                }
                continue;
            }

            // If ISTEP > NCP(L1), no need to recompute the distance to L1.
            // (L1 updated exactly M steps ago still needs it — hence `<=`.)
            if istep <= ncp[l1] as i64 {
                let mut da = 0.0f64;
                for j in 0..p {
                    let db = a[i * p + j] - c[l1 * p + j];
                    da += db * db;
                }
                d[i] = da * an1[l1];
            }

            // If ISTEP >= both NCP(L1) & NCP(L2) there will be no transfer here.
            if istep < ncp[l1] as i64 || istep < ncp[l2] as i64 {
                let r2 = d[i] / an2[l2];
                let mut dd = 0.0f64;
                let mut early = false;
                for j in 0..p {
                    let de = a[i * p + j] - c[l2 * p + j];
                    dd += de * de;
                    if dd >= r2 {
                        early = true;
                        break;
                    }
                }
                if !early {
                    // Move I from L1 to L2; centres update after every step.
                    icoun = 0;
                    *indx = 0;
                    itr[l1 + 1] = 1;
                    itr[l2 + 1] = 1;
                    let nm = istep + m as i64;
                    ncp[l1] = nm as i32;
                    ncp[l2] = nm as i32;
                    let al1 = nc[l1] as f64;
                    let alw = al1 - 1.0;
                    let al2 = nc[l2] as f64;
                    let alt = al2 + 1.0;
                    for j in 0..p {
                        c[l1 * p + j] = (c[l1 * p + j] * al1 - a[i * p + j]) / alw;
                        c[l2 * p + j] = (c[l2 * p + j] * al2 + a[i * p + j]) / alt;
                    }
                    nc[l1] -= 1;
                    nc[l2] += 1;
                    an2[l1] = alw / al1;
                    an1[l1] = if alw > 1.0 { alw / (alw - 1.0) } else { BIG };
                    an1[l2] = alt / al2;
                    an2[l2] = alt / (alt + 1.0);
                    ic1[i] = l2;
                    ic2[i] = l1;
                }
            }

            // If no re-allocation took place in the last M steps, return.
            if icoun == m {
                return;
            }
        }
    }
}
