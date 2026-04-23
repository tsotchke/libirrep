/* SPDX-License-Identifier: MIT */
/* NEON (aarch64) batch kernel for `irrep_wigner_d_matrix_batch`.
 *
 * Vectorisation strategy: lane 0 of every `float64x2_t` carries β_b, lane 1
 * carries β_{b+1}. The entire matrix build — trig, powers, Jacobi recurrence,
 * normalisation — runs in lockstep across the two β values. The Jacobi
 * recurrence's three-term update is the hot inner loop, and its coefficients
 * (which depend on α, β_jac, k but NOT on x) are scalar-shared across both
 * lanes; only the `x * c1` term differs per lane.
 *
 * Pair count `n_betas / 2` iterations handle the main loop; the tail
 * (`n_betas & 1`) falls through to the scalar path. */

#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)

#include <arm_neon.h>
#include <math.h>
#include <stddef.h>

#include <irrep/wigner_d.h>

#pragma STDC FP_CONTRACT OFF

/* Shared scalar helpers — same math as `irrep_wigner_d_small_2j` minus the
 * canonicalisation, which we perform once before entering the SIMD body. */

extern double irrep_wigner_d_small(int j, int mp, int m, double beta);

/* Per-iteration `x^n` accumulator. `n` is an integer exponent with `n <= 2j`.
 * Square-and-multiply keeps the operation count logarithmic. */
static inline float64x2_t ipow_v(float64x2_t x, int n) {
    float64x2_t r = vdupq_n_f64(1.0);
    float64x2_t base = x;
    while (n > 0) {
        if (n & 1)
            r = vmulq_f64(r, base);
        n >>= 1;
        if (n)
            base = vmulq_f64(base, base);
    }
    return r;
}

/* Evaluate P_n^{(α,β)}(x_v) via DLMF §18.9.1 recurrence, two x-lanes at once.
 * α, β, n are identical across lanes; only x differs. Returns a vector of
 * two P_n values. */
static float64x2_t jacobi_P_v_(int n, int alpha, int beta, float64x2_t x_v) {
    if (n < 0)
        return vdupq_n_f64(0.0);
    if (n == 0)
        return vdupq_n_f64(1.0);
    double      a = (double)alpha;
    double      b = (double)beta;
    float64x2_t p0 = vdupq_n_f64(1.0);
    double      p1_const = 0.5 * (a - b);
    double      p1_slope = 0.5 * (a + b + 2.0);
    float64x2_t p1 = vaddq_f64(vmulq_n_f64(x_v, p1_slope), vdupq_n_f64(p1_const));
    if (n == 1)
        return p1;
    for (int k = 1; k < n; ++k) {
        double dk = (double)k;
        double two_k_ab = 2.0 * dk + a + b;
        double denom = 2.0 * (dk + 1.0) * (dk + a + b + 1.0) * two_k_ab;
        double slope = (two_k_ab + 1.0) * two_k_ab * (two_k_ab + 2.0);
        double offset = (two_k_ab + 1.0) * (a * a - b * b);
        double c2 = 2.0 * (dk + a) * (dk + b) * (two_k_ab + 2.0);
        /* c1_v = slope · x + offset */
        float64x2_t c1_v = vaddq_f64(vmulq_n_f64(x_v, slope), vdupq_n_f64(offset));
        /* p2 = (c1_v · p1 − c2 · p0) / denom */
        float64x2_t num = vsubq_f64(vmulq_f64(c1_v, p1), vmulq_n_f64(p0, c2));
        float64x2_t p2 = vmulq_n_f64(num, 1.0 / denom);
        p0 = p1;
        p1 = p2;
    }
    return p1;
}

/* Canonicalise (mp, m) to m ≥ |mp|; returns the accumulated (-1)^{m'−m} sign
 * and writes canonicalised (mp, m). Mirror of `canonicalise_` in wigner_d.c
 * but for integer labels only (caller has already ruled out half-integer). */
static double canonicalise_int_(int *mp_io, int *m_io) {
    double sign = 1.0;
    int    mp = *mp_io;
    int    m = *m_io;
    if (m < mp) {
        int delta = mp - m;
        if (delta & 1)
            sign = -sign;
        int tmp = mp;
        mp = m;
        m = tmp;
    }
    if (m < -mp) {
        int tmp = -mp;
        mp = -m;
        m = tmp;
    }
    *mp_io = mp;
    *m_io = m;
    return sign;
}

/* Fill a single (mp, m) entry across both β-lanes. */
static float64x2_t d_small_pair_(int j, int mp_in, int m_in, float64x2_t cb2_v, float64x2_t sb2_v,
                                 float64x2_t x_v) {
    int mp = mp_in, m = m_in;
    if (m > j || m < -j || mp > j || mp < -j)
        return vdupq_n_f64(0.0);
    double sign = canonicalise_int_(&mp, &m);

    int mu = m - mp;       /* ≥ 0 */
    int nu = m + mp;       /* ≥ 0 */
    int n = j - m;         /* ≥ 0 */

    int    jpm = j + m;
    int    jmm = j - m;
    int    jpmp = j + mp;
    int    jmmp = j - mp;
    double lg_norm =
        0.5 * (lgamma((double)jpm + 1.0) + lgamma((double)jmm + 1.0) -
               lgamma((double)jpmp + 1.0) - lgamma((double)jmmp + 1.0));
    double norm = exp(lg_norm);

    float64x2_t cb2_pow = ipow_v(cb2_v, nu);
    float64x2_t sb2_pow = ipow_v(sb2_v, mu);
    float64x2_t Pn = jacobi_P_v_(n, mu, nu, x_v);

    float64x2_t prod = vmulq_f64(cb2_pow, sb2_pow);
    prod = vmulq_f64(prod, Pn);
    prod = vmulq_n_f64(prod, sign * norm);
    return prod;
}

/* Build the full (2j+1)² matrix for a pair of β values. Output layout:
 *   out_A = out for β_A (first β), out_B = out for β_B (second β). */
static void fill_matrix_pair_(int j, double beta_A, double beta_B, double *out_A, double *out_B) {
    int d = 2 * j + 1;
    /* Zero the anti-diagonals as a guard (keeps parity with the scalar
     * `_matrix` which pre-zeroes them). */
    for (int imp = 0; imp < d; ++imp) {
        out_A[imp * d + (d - 1 - imp)] = 0.0;
        out_B[imp * d + (d - 1 - imp)] = 0.0;
    }

    double betas[2] = {beta_A, beta_B};
    (void)betas;
    float64x2_t beta_v = {beta_A, beta_B};
    float64x2_t half_v = vdupq_n_f64(0.5);
    float64x2_t beta_half = vmulq_f64(beta_v, half_v);
    double      cb2_0 = cos(vgetq_lane_f64(beta_half, 0));
    double      cb2_1 = cos(vgetq_lane_f64(beta_half, 1));
    double      sb2_0 = sin(vgetq_lane_f64(beta_half, 0));
    double      sb2_1 = sin(vgetq_lane_f64(beta_half, 1));
    double      x_0 = cos(beta_A);
    double      x_1 = cos(beta_B);
    float64x2_t cb2_v = {cb2_0, cb2_1};
    float64x2_t sb2_v = {sb2_0, sb2_1};
    float64x2_t x_v = {x_0, x_1};

    for (int imp = 0; imp < d; ++imp) {
        int mp = imp - j;
        int im_lo = (mp < 0) ? -mp + j : mp + j;
        for (int im = im_lo; im < d; ++im) {
            int         m = im - j;
            float64x2_t v = d_small_pair_(j, mp, m, cb2_v, sb2_v, x_v);
            double      vA = vgetq_lane_f64(v, 0);
            double      vB = vgetq_lane_f64(v, 1);
            out_A[imp * d + im] = vA;
            out_B[imp * d + im] = vB;
            int imp_a = -m + j, im_a = -mp + j;
            out_A[imp_a * d + im_a] = vA;
            out_B[imp_a * d + im_a] = vB;
            if (mp != m) {
                double signA = ((m - mp) & 1) ? -vA : vA;
                double signB = ((m - mp) & 1) ? -vB : vB;
                out_A[im * d + imp] = signA;
                out_B[im * d + imp] = signB;
                int imp_b = -mp + j, im_b = -m + j;
                if (imp_b != im || im_b != imp) {
                    out_A[imp_b * d + im_b] = signA;
                    out_B[imp_b * d + im_b] = signB;
                }
            }
        }
    }
}

extern void irrep_wigner_d_matrix(int j, double *out, double beta);

void irrep_wigner_d_matrix_batch_neon(int j, size_t n_betas, const double *betas, double *out) {
    int    d = 2 * j + 1;
    size_t stride = (size_t)d * (size_t)d;
    size_t b = 0;
    for (; b + 2 <= n_betas; b += 2) {
        fill_matrix_pair_(j, betas[b], betas[b + 1], out + b * stride, out + (b + 1) * stride);
    }
    if (b < n_betas) {
        irrep_wigner_d_matrix(j, out + b * stride, betas[b]);
    }
}

#endif /* aarch64 */
