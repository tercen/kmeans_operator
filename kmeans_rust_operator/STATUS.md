# STATUS

## Memory model — measured 2026-09-23

`tests/memory.rs` (release, jemalloc, `MALLOC_CONF=background_thread:true,
dirty_decay_ms:0,muzzy_decay_ms:0,thp:never`), 10 000 000 cells
(2 000 000 observations × 5 variables):

| quantity | value |
|---|---|
| gather peak (sums + counts) | 119 880 kB |
| after in-place mean, counts freed | 119 296 kB |
| final peak incl. Hartigan–Wong | 126 716 kB (**123 MB** ≈ 12.4 B/cell) |
| Hartigan–Wong fit time (k = 5) | 0.4 s |

Booked in `memory_model.json`: `2.0e-05 MB/cell × cells + 150 MB offset`
(= 425 MB for 10 M cells) — about 3× the measured peak, deliberately, since
the measurement excludes the runtime and the upload path. Refit from
`stats_d_actual_ram_peak` once there are real runs on the platform.

## Platform verification not yet done

- `tests/test.json` golden has not run through a Studio install-check.
- No production-path dry run (§8a of create-rust-operator) against a live
  instance yet: the dev binary exists (`src/bin/dev.rs`) but no Studio was
  available while building.
