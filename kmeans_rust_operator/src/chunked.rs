//! Reading the crosstab cell table in chunks (the shared plumbing, copied
//! from `asinh_rust_operator/src/input.rs`).
//!
//! The cell table is never held whole by the reader: `stream_tson` hands back
//! TSON bytes and the caller folds them. A chunk is decoded with `rustson`
//! rather than polars, so a chunk of 1 M cells costs the three column vectors
//! and nothing else.
use anyhow::{Result, anyhow, bail};
use std::collections::HashMap;

use tercen_rs::context::ContextBase;

/// One decoded chunk of the crosstab: the columns asked for, in row order.
#[derive(Debug, Default)]
pub struct Chunk {
    pub ri: Vec<i32>,
    pub ci: Vec<i32>,
    pub y: Vec<f64>,
}

impl Chunk {
    /// Rows in this chunk, whichever columns were requested.
    pub fn len(&self) -> usize {
        self.y.len().max(self.ri.len()).max(self.ci.len())
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Read the cell table in chunks, stopping at the row count the schema declared.
///
/// `TableStreamer::stream_table_chunked` cannot be used for this: it decides it has reached the
/// end when a chunk comes back as **zero bytes**, and a TSON table with no rows is not zero
/// bytes. On a 1,000-cell crosstab it therefore kept asking for the next million rows until the
/// server rejected an offset above `i32::MAX` ("Illegal to set field start … out of range for
/// signed 32-bit int"). Counting rows against `Schema.nRows` is both correct and bounded.
pub async fn for_each_chunk<F>(
    ctx: &ContextBase,
    columns: &[&str],
    n_cells: usize,
    chunk: usize,
    mut f: F,
) -> Result<()>
where
    F: FnMut(Chunk) -> Result<()>,
{
    let want = |n: &str| columns.contains(&n);
    let cols: Vec<String> = columns.iter().map(|c| c.to_string()).collect();
    let mut offset = 0usize;
    while offset < n_cells {
        let want_rows = chunk.min(n_cells - offset);
        let bytes = ctx
            .streamer()
            .stream_tson(
                ctx.qt_hash(),
                Some(cols.clone()),
                offset as i64,
                want_rows as i64,
            )
            .await
            .map_err(|e| anyhow!("read cells {offset}..{}: {e}", offset + want_rows))?;
        let c = decode_chunk(&bytes, want(".ri"), want(".ci"), want(".y"))?;
        if c.is_empty() {
            break;
        }
        offset += c.len();
        f(c)?;
    }
    if offset != n_cells {
        bail!(
            "the crosstab returned {offset} cells where its schema says {n_cells}; the table \
             changed under the operator"
        );
    }
    Ok(())
}

/// Number of cells in the projected crosstab, from the table schema (`Schema.nRows`).
pub async fn cell_count(ctx: &ContextBase) -> Result<usize> {
    let schema = ctx
        .streamer()
        .get_schema(ctx.qt_hash())
        .await
        .map_err(|e| anyhow!("get schema {}: {e}", ctx.qt_hash()))?;
    use tercen_rs::client::proto::e_schema;
    let n = match schema.object {
        Some(e_schema::Object::Schema(s)) => s.n_rows,
        Some(e_schema::Object::Tableschema(s)) => s.n_rows,
        Some(e_schema::Object::Computedtableschema(s)) => s.n_rows,
        Some(e_schema::Object::Cubequerytableschema(s)) => s.n_rows,
        None => bail!("schema {} has no object", ctx.qt_hash()),
    };
    usize::try_from(n).map_err(|_| anyhow!("schema reports {n} rows"))
}

/// Decode one TSON chunk of the cell table into the requested columns.
///
/// A response is **not one document**. `stream_tson` concatenates every gRPC message it
/// receives, and the server pages its answer — about 15,000 rows per page here — so the buffer
/// holds one complete TSON document per page, back to back. Decoding only the first is what
/// made a large request pathological: asking for a million rows transferred a million rows and
/// used fifteen thousand of them, then asked again from a slightly later offset. Reading every
/// document in the buffer turns the same request into 775,000 cells a second instead of 30,000.
pub fn decode_chunk(bytes: &[u8], want_ri: bool, want_ci: bool, want_y: bool) -> Result<Chunk> {
    let mut out = Chunk::default();
    let mut cur = std::io::Cursor::new(bytes);
    while (cur.position() as usize) < bytes.len() {
        let before = cur.position();
        let v = rustson::decode(cur.clone()).map_err(|e| anyhow!("decode chunk: {e:?}"))?;
        // `rustson::decode` takes the cursor by value, so advance our own by the document's
        // length: re-encoding is not an option, and the decoder does not report what it read.
        let consumed = rustson::encode(&v)
            .map_err(|e| anyhow!("measure chunk: {e:?}"))?
            .len();
        cur.set_position(before + consumed as u64);
        decode_one(&v, want_ri, want_ci, want_y, &mut out)?;
        if consumed == 0 {
            break;
        }
    }
    Ok(out)
}

fn decode_one(
    v: &rustson::Value,
    want_ri: bool,
    want_ci: bool,
    want_y: bool,
    out: &mut Chunk,
) -> Result<()> {
    let rustson::Value::MAP(m) = v else {
        bail!("chunk is not a map")
    };
    for (name, want) in [(".ri", want_ri), (".ci", want_ci), (".y", want_y)] {
        if !want {
            continue;
        }
        let values = find_column(m, name).ok_or_else(|| {
            anyhow!(
                "chunk has no column '{name}' (it has {:?})",
                column_names(m)
            )
        })?;
        match name {
            ".ri" => out.ri.extend_from_slice(&as_i32(values, name)?),
            ".ci" => out.ci.extend_from_slice(&as_i32(values, name)?),
            _ => out.y.extend_from_slice(&as_f64(values, name)?),
        }
    }
    Ok(())
}

fn as_i32(v: &rustson::Value, name: &str) -> Result<Vec<i32>> {
    Ok(match v {
        rustson::Value::LSTI32(v) => v.clone(),
        rustson::Value::LSTF64(v) => v.iter().map(|x| *x as i32).collect(),
        rustson::Value::LSTU8(v) => v.iter().map(|x| *x as i32).collect(),
        other => bail!("column '{name}' is not an integer list ({other:?})"),
    })
}

fn as_f64(v: &rustson::Value, name: &str) -> Result<Vec<f64>> {
    Ok(match v {
        rustson::Value::LSTF64(v) => v.clone(),
        rustson::Value::LSTI32(v) => v.iter().map(|x| *x as f64).collect(),
        rustson::Value::LSTU8(v) => v.iter().map(|x| *x as f64).collect(),
        other => bail!("column '{name}' is not a numeric list ({other:?})"),
    })
}

/// Find a named column's values in a decoded TSON table.
///
/// Tercen uses **two** layouts and they are easy to confuse. A table streamed out of the server
/// is `{"cols": [{"name", "type", "data"}]}`; the `OperatorResult` an operator writes back is
/// `{"tables": [{"columns": [{"name", "values"}]}]}`. Reading a streamed chunk with the second
/// shape fails with "chunk has no columns", which says nothing about the real cause, so both are
/// accepted here.
fn find_column<'a>(
    m: &'a HashMap<String, rustson::Value>,
    name: &str,
) -> Option<&'a rustson::Value> {
    for (list_key, value_key) in [("cols", "data"), ("columns", "values")] {
        let Some(rustson::Value::LST(cols)) = m.get(list_key) else {
            continue;
        };
        for c in cols {
            let rustson::Value::MAP(c) = c else { continue };
            let Some(rustson::Value::STR(n)) = c.get("name") else {
                continue;
            };
            if n == name {
                return c.get(value_key);
            }
        }
    }
    None
}

/// Every column name a decoded TSON table carries, for error messages.
fn column_names(m: &HashMap<String, rustson::Value>) -> Vec<String> {
    let mut out = Vec::new();
    for list_key in ["cols", "columns"] {
        if let Some(rustson::Value::LST(cols)) = m.get(list_key) {
            for c in cols {
                if let rustson::Value::MAP(c) = c
                    && let Some(rustson::Value::STR(n)) = c.get("name")
                {
                    out.push(n.clone());
                }
            }
        }
    }
    out
}
