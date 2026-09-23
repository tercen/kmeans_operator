//! Reading the crosstab and gathering it into the matrix `kmeans` clusters.
//!
//! The projection (R `main.R`): rows are variables, columns are observations,
//! y is the measurement. R builds its data with
//! `reshape2::acast(.ci ~ .ri, value.var = ".y", fill = NaN,
//! fun.aggregate = mean)` — one point per **column** factor value, one
//! dimension per **row** factor value, duplicates averaged, missing cells
//! `NaN` — and clusters the resulting rows.
//!
//! **This gather is the memory pay-off of the port.** R materialises the whole
//! long-format cell table (`.ci`, `.ri`, `.y` — 24 B/cell before R's per-element
//! overhead) and then the dense matrix on top of it. Here the cell table is
//! read in bounded chunks (`chunked::for_each_chunk`, counted against
//! `Schema.nRows`) and folded straight into two dense buffers — value sums
//! (8 B/cell) and occurrence counts (4 B/cell) sized from the row/column
//! factor schemas, which the task already carries. Nothing of the cell table
//! is ever held whole; peak working memory is 12 B/cell plus one chunk, and
//! after the mean the counts are freed and `kmeans` runs on the sums buffer
//! itself (8 B/cell). See `memory_model.json` and `CLAUDE.md`.
//!
//! Without the factor schemas (no row/column hash on the task) the gather
//! falls back to a one-pass map of observed cells, as `clusterx_rust_operator`
//! does; the dense path is the normal one.

use anyhow::{Result, bail};
use std::collections::HashMap;

use tercen_rs::context::ContextBase;

use crate::chunked::{cell_count, for_each_chunk};

/// The gathered crosstab, in the layout the R operator hands to `kmeans`:
/// one point per `.ci` (rows of the matrix), one dimension per `.ri`.
pub struct CrosstabMatrix {
    /// Row-major `n_points × n_vars`; entry `(i, k)` is the mean y of the
    /// cells whose column is `ci_values[i]` and whose row is `ri_values[k]`.
    pub data: Vec<f64>,
    /// The `.ci` values in matrix row order (ascending).
    pub ci_values: Vec<i32>,
    /// The `.ri` values in matrix column order (ascending).
    pub ri_values: Vec<i32>,
}

impl CrosstabMatrix {
    pub fn n_points(&self) -> usize {
        self.ci_values.len()
    }
    pub fn n_vars(&self) -> usize {
        self.ri_values.len()
    }
}

/// Stream the whole cell table and aggregate it exactly like
/// `acast(.ci ~ .ri, fun.aggregate = mean, fill = NaN)`.
pub async fn gather_matrix(ctx: &ContextBase) -> Result<CrosstabMatrix> {
    let n_cells = cell_count(ctx).await?;
    match factor_cardinalities(ctx).await? {
        Some((n_vars, n_points)) => gather_dense(ctx, n_cells, n_vars, n_points).await,
        None => gather_observed(ctx, n_cells).await,
    }
}

/// Distinct factor counts from the row and column factor schemas:
/// `(n_rows, n_cols)` of the crosstab — the dimensions of the matrix.
///
/// `.ri` indexes the row factor table and `.ci` the column factor table, both
/// zero-based and dense, so the schema row counts size the gather exactly.
/// `stream_tson` cannot be asked for this (the cell table does not carry it),
/// which is why the R port leans on `reshape2` instead.
async fn factor_cardinalities(ctx: &ContextBase) -> Result<Option<(usize, usize)>> {
    if ctx.row_hash().is_empty() || ctx.column_hash().is_empty() {
        return Ok(None);
    }
    let row_n = schema_n_rows(ctx, ctx.row_hash()).await?;
    let col_n = schema_n_rows(ctx, ctx.column_hash()).await?;
    match (row_n, col_n) {
        (Some(r), Some(c)) if r > 0 && c > 0 => Ok(Some((r, c))),
        _ => Ok(None),
    }
}

async fn schema_n_rows(ctx: &ContextBase, hash: &str) -> Result<Option<usize>> {
    let schema = match ctx.streamer().get_schema(hash).await {
        Ok(s) => s,
        Err(e) => {
            tracing::warn!(hash, error = %e, "could not fetch factor schema; using the observed-cell gather");
            return Ok(None);
        }
    };
    use tercen_rs::client::proto::e_schema;
    let n = match schema.object {
        Some(e_schema::Object::Schema(s)) => s.n_rows,
        Some(e_schema::Object::Tableschema(s)) => s.n_rows,
        Some(e_schema::Object::Computedtableschema(s)) => s.n_rows,
        Some(e_schema::Object::Cubequerytableschema(s)) => s.n_rows,
        None => bail!("factor schema {hash} has no object"),
    };
    Ok(usize::try_from(n).ok())
}

/// The dense gather: one pass over the cell table into sums and counts.
async fn gather_dense(
    ctx: &ContextBase,
    n_cells: usize,
    n_vars: usize,
    n_points: usize,
) -> Result<CrosstabMatrix> {
    // The cell table should be exactly the projection's cells; if it is not
    // (duplicated (Variable, Observation) pairs), the sums below still mean
    // what acast's `fun.aggregate = mean` means, and the empty-cell check at
    // the end decides completeness.
    if n_cells != n_points.saturating_mul(n_vars) {
        tracing::warn!(
            n_cells,
            expected = n_points * n_vars,
            "cell count differs from the factor-grid size; continuing (duplicates are averaged)"
        );
    }
    let mut sums = vec![0.0f64; n_points * n_vars];
    let mut counts = vec![0u32; n_points * n_vars];
    for_each_chunk(ctx, &[".ri", ".ci", ".y"], n_cells, CHUNK, |c| {
        fold_chunk(&mut sums, &mut counts, n_vars, n_points, &c)
    })
    .await?;

    // Mean in place: the sums buffer becomes the matrix `kmeans` clusters.
    let mut missing = 0usize;
    for cell in 0..n_points * n_vars {
        if counts[cell] == 0 {
            missing += 1;
        } else {
            sums[cell] /= counts[cell] as f64;
        }
    }
    drop(counts);
    if missing > 0 {
        bail!(
            "the crosstab is missing {missing} of {} measurements (an observation × variable \
             pair has no value); every cell of the projection must have a value — the R \
             reference fills these with NaN and the clustering is meaningless",
            n_points * n_vars
        );
    }
    Ok(CrosstabMatrix {
        data: sums,
        ci_values: (0..n_points as i32).collect(),
        ri_values: (0..n_vars as i32).collect(),
    })
}

/// Fallback gather for a task without row/column factor schemas: one pass
/// over observed cells into a map, then the dense matrix (the
/// `clusterx_rust_operator` gather). Peak is the map (~50 B/cell), not 12.
async fn gather_observed(ctx: &ContextBase, n_cells: usize) -> Result<CrosstabMatrix> {
    let mut sums: HashMap<(i32, i32), (f64, u32)> = HashMap::with_capacity(n_cells);
    for_each_chunk(ctx, &[".ri", ".ci", ".y"], n_cells, CHUNK, |c| {
        for e in 0..c.len() {
            let y = c.y[e];
            if !y.is_finite() {
                bail!(
                    "the crosstab has a non-finite measurement (cell (.ci={}, .ri={}) is {y})",
                    c.ci[e],
                    c.ri[e]
                );
            }
            let entry = sums.entry((c.ci[e], c.ri[e])).or_insert((0.0, 0));
            entry.0 += y;
            entry.1 += 1;
        }
        Ok(())
    })
    .await?;

    let mut ci_values: Vec<i32> = sums.keys().map(|k| k.0).collect();
    ci_values.sort_unstable();
    ci_values.dedup();
    let mut ri_values: Vec<i32> = sums.keys().map(|k| k.1).collect();
    ri_values.sort_unstable();
    ri_values.dedup();
    let n_points = ci_values.len();
    let n_vars = ri_values.len();
    if n_points < 1 {
        bail!("the projection has no observation column factor values");
    }
    if n_vars < 1 {
        bail!("the projection has no variable row factor values");
    }

    let mut ci_rank: HashMap<i32, usize> = HashMap::with_capacity(n_points);
    for (i, &v) in ci_values.iter().enumerate() {
        ci_rank.insert(v, i);
    }
    let mut ri_rank: HashMap<i32, usize> = HashMap::with_capacity(n_vars);
    for (k, &v) in ri_values.iter().enumerate() {
        ri_rank.insert(v, k);
    }

    let mut data = vec![f64::NAN; n_points * n_vars];
    for ((ci, ri), (sum, count)) in &sums {
        let (i, k) = (ci_rank[ci], ri_rank[ri]);
        data[i * n_vars + k] = sum / *count as f64;
    }
    Ok(CrosstabMatrix {
        data,
        ci_values,
        ri_values,
    })
}

/// Cells fetched per gRPC round trip (same measurement as asinh's).
pub const CHUNK: usize = 200_000;

/// Fold one decoded chunk into the sums/counts buffers (`acast`'s
/// `fun.aggregate = mean`, made explicit). Pure, so the aggregation
/// semantics are unit-tested without a server.
fn fold_chunk(
    sums: &mut [f64],
    counts: &mut [u32],
    n_vars: usize,
    n_points: usize,
    c: &crate::chunked::Chunk,
) -> Result<()> {
    for e in 0..c.len() {
        let (ri, ci) = (c.ri[e], c.ci[e]);
        if ri < 0 || ci < 0 {
            bail!("the crosstab has a negative factor index (.ri={ri}, .ci={ci})");
        }
        let (ri, ci) = (ri as usize, ci as usize);
        if ri >= n_vars || ci >= n_points {
            bail!(
                "the crosstab carries cell (.ri={ri}, .ci={ci}) but the factor schemas say \
                 {n_vars} rows × {n_points} columns; the projection and its schemas disagree"
            );
        }
        let y = c.y[e];
        if !y.is_finite() {
            bail!(
                "the crosstab has a non-finite measurement (cell (.ci={ci}, .ri={ri}) is \
                 {y}); every value must be finite"
            );
        }
        let cell = ci * n_vars + ri;
        sums[cell] += y;
        counts[cell] += 1;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chunked::Chunk;

    #[test]
    fn duplicate_cells_are_averaged_like_acast() {
        // acast(.ci ~ .ri, fun.aggregate = mean): two values for the same
        // (ci, ri) collapse to their mean.
        let c = Chunk {
            ri: vec![0, 1, 0],
            ci: vec![0, 0, 0],
            y: vec![1.0, 5.0, 3.0],
        };
        let mut sums = vec![0.0; 2];
        let mut counts = vec![0u32; 2];
        fold_chunk(&mut sums, &mut counts, 2, 1, &c).unwrap();
        let means: Vec<f64> = sums
            .iter()
            .zip(&counts)
            .map(|(s, n)| s / *n as f64)
            .collect();
        assert_eq!(means, vec![2.0, 5.0]);
    }

    #[test]
    fn an_index_outside_the_factor_schemas_is_an_error() {
        let c = Chunk {
            ri: vec![5],
            ci: vec![0],
            y: vec![1.0],
        };
        let mut sums = vec![0.0; 5];
        let mut counts = vec![0u32; 5];
        assert!(fold_chunk(&mut sums, &mut counts, 5, 1, &c).is_err());
    }

    #[test]
    fn a_non_finite_measurement_is_an_error() {
        let c = Chunk {
            ri: vec![0],
            ci: vec![0],
            y: vec![f64::INFINITY],
        };
        let mut sums = vec![0.0; 1];
        let mut counts = vec![0u32; 1];
        let err = fold_chunk(&mut sums, &mut counts, 1, 1, &c).unwrap_err();
        assert!(err.to_string().contains("non-finite"));
    }
}
