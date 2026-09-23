//! Peak-memory measurement for the memory model (create-rust-operator §7).
//!
//! Not part of the normal suite: `cargo test --release --ignored -- --nocapture memory_bench`.
//! It reproduces the operator's working set at scale without a server: the
//! dense sums+counts gather buffers (12 B/cell), the in-place mean, and the
//! Hartigan–Wong run on the sums buffer. The production run adds one read
//! chunk (200 000 cells ≈ 5 MB) and the tokio/tonic runtime to this.

use kmeans_rust_operator::kmns::kmns;

fn peak_rss_kb() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .unwrap()
        .lines()
        .find(|l| l.starts_with("VmHWM:"))
        .unwrap()
        .split_whitespace()
        .nth(1)
        .unwrap()
        .parse()
        .unwrap()
}

#[test]
#[ignore]
fn memory_bench_10m_cells() {
    // 2 000 000 observations × 5 variables = 10 000 000 cells.
    let (m, p) = (2_000_000usize, 5usize);
    let cells = m * p;
    let mut sums = vec![0.0f64; cells];
    let mut counts = vec![0u32; cells];
    for i in 0..m {
        for j in 0..p {
            let cell = i * p + j;
            sums[cell] = ((i % 97) as f64) * 0.1 + (j as f64);
            counts[cell] = 1;
        }
    }
    let before_gather = peak_rss_kb();
    for cell in 0..cells {
        sums[cell] /= counts[cell] as f64;
    }
    drop(counts);
    let before_kmeans = peak_rss_kb();

    let t = std::time::Instant::now();
    // Initial centres: the first k points, as the wrapper would hand over.
    let centers: Vec<f64> = (0..5)
        .flat_map(|i| sums[i * p..(i + 1) * p].to_vec())
        .collect();
    let result = kmns(&sums, m, p, &mut centers.clone(), 5, 10, 50 * m);
    let dt = t.elapsed().as_secs_f64();
    let peak = peak_rss_kb();
    println!(
        "10M cells: gather peak {before_gather} kB, after mean {before_kmeans} kB, \
         final peak {peak} kB ({} MB), kmeans {dt:.1}s, {} clusters",
        peak / 1024,
        result
            .cluster
            .iter()
            .collect::<std::collections::HashSet<_>>()
            .len()
    );
}
