//! R's random number generator, exactly.
//!
//! FlowSOM's map depends on the RNG twice — the initial codes are a `sample()` of the data, and
//! every training step draws an observation with `unif_rand()` — so a port that does not
//! reproduce R's stream cannot be compared with R's answer at all. Everything here is a
//! transliteration of `src/main/RNG.c` and `src/main/random.c` in R 4.0.4.
//!
//! Two details matter more than they look:
//!
//! - `set.seed` does not seed Mersenne-Twister with the seed. It scrambles it through 50 rounds
//!   of a linear congruential generator first, then fills all 625 state words from the same LCG.
//! - Since R 3.6 the default `sample.kind` is `"Rejection"`, not the old rounding. `sample()`
//!   draws whole 16-bit blocks, masks them down to the bits it needs and rejects anything out of
//!   range, so it consumes a variable number of uniforms.

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;
const TEMPERING_MASK_B: u32 = 0x9d2c_5680;
const TEMPERING_MASK_C: u32 = 0xefc6_0000;

/// `i2_32m1` in R: the smallest step away from an endpoint that `fixup` uses.
const I2_32M1: f64 = 2.328_306_437_080_797e-10;

/// R's Mersenne-Twister, seeded the way `set.seed` seeds it.
pub struct RRng {
    mt: [u32; N],
    mti: usize,
}

impl RRng {
    /// `set.seed(seed)` with the default `kind` and `sample.kind`.
    pub fn set_seed(seed: u32) -> Self {
        // RNG_Init: initial scrambling, then one LCG step per state word. `i_seed[0]` is the
        // index and is overwritten by FixupSeeds, so it is drawn and discarded here.
        let mut s = seed;
        let mut lcg = || {
            s = s.wrapping_mul(69069).wrapping_add(1);
            s
        };
        for _ in 0..50 {
            lcg();
        }
        lcg(); // i_seed[0], replaced by FixupSeeds below
        let mut mt = [0u32; N];
        for w in mt.iter_mut() {
            *w = lcg();
        }
        // FixupSeeds(MERSENNE_TWISTER, initial = 1): I624 = 624, so the first draw regenerates.
        Self { mt, mti: N }
    }

    /// `MT_genrand`: the raw generator, before `fixup`.
    fn genrand(&mut self) -> f64 {
        if self.mti >= N {
            let mag01 = [0u32, MATRIX_A];
            for kk in 0..N - M {
                let y = (self.mt[kk] & UPPER_MASK) | (self.mt[kk + 1] & LOWER_MASK);
                self.mt[kk] = self.mt[kk + M] ^ (y >> 1) ^ mag01[(y & 1) as usize];
            }
            for kk in N - M..N - 1 {
                let y = (self.mt[kk] & UPPER_MASK) | (self.mt[kk + 1] & LOWER_MASK);
                self.mt[kk] = self.mt[kk + M - N] ^ (y >> 1) ^ mag01[(y & 1) as usize];
            }
            let y = (self.mt[N - 1] & UPPER_MASK) | (self.mt[0] & LOWER_MASK);
            self.mt[N - 1] = self.mt[M - 1] ^ (y >> 1) ^ mag01[(y & 1) as usize];
            self.mti = 0;
        }
        let mut y = self.mt[self.mti];
        self.mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & TEMPERING_MASK_B;
        y ^= (y << 15) & TEMPERING_MASK_C;
        y ^= y >> 18;
        (y as f64) * 2.328_306_436_538_696_3e-10
    }

    /// `unif_rand()`: the generator with R's `fixup`, so the result is strictly inside (0, 1).
    pub fn unif_rand(&mut self) -> f64 {
        let x = self.genrand();
        if x <= 0.0 {
            0.5 * I2_32M1
        } else if 1.0 - x <= 0.0 {
            1.0 - 0.5 * I2_32M1
        } else {
            x
        }
    }

    /// `rbits(bits)`: whole 16-bit blocks, then masked down.
    fn rbits(&mut self, bits: i32) -> f64 {
        let mut v: i64 = 0;
        let mut n = 0;
        while n <= bits {
            let v1 = (self.unif_rand() * 65536.0).floor() as i64;
            v = 65536 * v + v1;
            n += 16;
        }
        (v & ((1i64 << bits) - 1)) as f64
    }

    /// `R_unif_index(dn)` with `sample.kind = "Rejection"`: a uniform integer in `[0, dn)`.
    pub fn unif_index(&mut self, dn: f64) -> f64 {
        if dn <= 0.0 {
            return 0.0;
        }
        let bits = dn.log2().ceil() as i32;
        loop {
            let dv = self.rbits(bits);
            if dn > dv {
                return dv;
            }
        }
    }

    /// `sample(1:n, k, replace = FALSE)`, returning **0-based** indices in draw order.
    ///
    /// R's partial Fisher-Yates: draw a position, take what is there, and move the last element
    /// into the hole. The order of the draws is part of the answer, not an implementation
    /// detail — the codes a map starts from are these rows in this order.
    pub fn sample_int(&mut self, n: usize, k: usize) -> Vec<usize> {
        let mut ix: Vec<usize> = (0..n).collect();
        let mut left = n;
        let mut out = Vec::with_capacity(k);
        for _ in 0..k {
            let j = self.unif_index(left as f64) as usize;
            out.push(ix[j]);
            left -= 1;
            ix[j] = ix[left];
        }
        out
    }
}
