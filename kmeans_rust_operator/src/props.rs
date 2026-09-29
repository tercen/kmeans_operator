//! Operator properties, with the R operator's defaults
//! (`kmeans_operator/main.R`): `centers` (10), `iter.max` (10), `nstart` (1)
//! and `seed` ("NULL" = R's `set.seed(NULL)`, a time-seeded, unreproducible
//! session — the port uses OS entropy there).
use anyhow::{Result, bail};
use tercen_rs::PropertyReader;
use tercen_rs::context::ContextBase;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Settings {
    pub centers: usize,
    pub iter_max: usize,
    pub nstart: usize,
    /// `Some(s)` = R `set.seed(s)`; `None` = R `set.seed(NULL)` (entropy here).
    pub seed: Option<i32>,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            centers: 10,
            iter_max: 10,
            nstart: 1,
            seed: None,
        }
    }
}

/// Read the properties off the task's `CubeQueryTask` snapshot.
pub fn settings_from_ctx(ctx: &ContextBase) -> Result<Settings> {
    let pr = PropertyReader::from_operator_settings(ctx.operator_settings());
    settings_from(|name, default| pr.get_string(name, default))
}

/// Parse the properties the R operator reads in `main.R`, through a getter
/// returning the raw string (so tests can feed property tables directly).
pub fn settings_from(get: impl Fn(&str, &str) -> String) -> Result<Settings> {
    // Tercen serialises a numeric property as e.g. "2.0" even where the
    // operator wants an integer, so parse every number as f64 and cast
    // (create-rust-operator §2). R reads them with `as.integer`, which
    // truncates toward zero.
    let int = |name: &str, default: f64| -> Result<f64> {
        let raw = get(name, &default.to_string());
        let v: f64 = raw
            .trim()
            .parse()
            .map_err(|_| anyhow::anyhow!("property '{name}' is not a number: '{raw}'"))?;
        if !v.is_finite() {
            bail!("property '{name}' must be a finite number, got {v}");
        }
        if !(-2_147_483_648.0..=2_147_483_647.0).contains(&v) {
            bail!("property '{name}' does not fit in an R integer: {v}");
        }
        Ok(v.trunc())
    };
    let centers = int("centers", 10.0)? as i32;
    let iter_max = int("iter.max", 10.0)? as i32;
    let nstart = int("nstart", 1.0)? as i32;
    let seed = match get("seed", "NULL").trim() {
        "NULL" => None,
        raw => {
            let v: f64 = raw.parse().map_err(|_| {
                anyhow::anyhow!("property 'seed' is not a number or 'NULL': '{raw}'")
            })?;
            if !v.is_finite() || !(-2_147_483_648.0..=2_147_483_647.0).contains(&v) {
                bail!("property 'seed' must be a finite R integer or 'NULL', got '{raw}'");
            }
            Some(v.trunc() as i32)
        }
    };
    if nstart < 1 {
        // R errors here too, but with an opaque "argument is of length zero"
        // from its own wrapper (a scalar `centers` never becomes a matrix).
        bail!("property 'nstart' must be at least 1, got {nstart}");
    }
    Ok(Settings {
        centers: usize::try_from(centers)?,
        iter_max: usize::try_from(iter_max)?,
        nstart: usize::try_from(nstart)?,
        seed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn table<'a>(pairs: &'a [(&'a str, &'a str)]) -> impl Fn(&str, &str) -> String + 'a {
        let m: HashMap<&str, &str> = pairs.iter().copied().collect();
        move |name, default| {
            m.get(name)
                .map(|s| s.to_string())
                .unwrap_or_else(|| default.to_string())
        }
    }

    #[test]
    fn defaults_match_the_r_operator() {
        let s = settings_from(|_, default| default.to_string()).unwrap();
        assert_eq!(s, Settings::default());
    }

    #[test]
    fn a_double_property_written_as_2_point_0_parses_as_2() {
        // The UI stores a DoubleProperty as "2.0"; R's as.integer truncates.
        let s = settings_from(table(&[("centers", "2.0"), ("iter.max", "5.0")])).unwrap();
        assert_eq!(s.centers, 2);
        assert_eq!(s.iter_max, 5);
        assert_eq!(s.nstart, 1);
        assert_eq!(s.seed, None);
    }

    #[test]
    fn a_numeric_seed_is_the_parity_path() {
        let s = settings_from(table(&[("seed", "42")])).unwrap();
        assert_eq!(s.seed, Some(42));
    }

    #[test]
    fn nstart_below_one_is_refused() {
        assert!(settings_from(table(&[("nstart", "0")])).is_err());
    }

    #[test]
    fn a_non_numeric_property_is_an_error_not_a_silent_default() {
        assert!(settings_from(table(&[("centers", "many")])).is_err());
    }
}
