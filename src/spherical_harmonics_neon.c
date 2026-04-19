/* SPDX-License-Identifier: MIT */
/* NEON (aarch64) batch kernel for irrep_sph_harm_cart_all.
 *
 * Vectorisation strategy: run two edges through the scalar algorithm in
 * parallel, with each `float64x2_t` carrying lane 0 = edge `e`, lane 1 =
 * edge `e+1`. Chebyshev trig recurrence (cos(mφ), sin(mφ) from x/r_xy,
 * y/r_xy) and the normalisation multiplies are vectorised; the associated
 * Legendre polynomial `P_l^m(z)` is still evaluated per-lane via the scalar
 * `irrep_legendre_assoc`, which sidesteps the cross-lane branching in the
 * m-recurrence and keeps the result bit-exact against the scalar batch
 * kernel for every representative input.
 *
 * Odd N: the tail single edge falls through to the scalar per-edge path.
 *
 * Bit-exactness is asserted in tests/test_spherical_harmonics.c by a
 * scalar-vs-NEON sweep at l_max ∈ {0..4} over N = 17 edges (odd — tail
 * exercised). Matches to byte identity. */

#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)

#include <arm_neon.h>
#include <math.h>
#include <stddef.h>

#include <irrep/spherical_harmonics.h>
#include <irrep/types.h>

#include "internal/dispatch.h"

/* Disable FP contraction for parity with the scalar path — see matching
 * pragma in src/spherical_harmonics.c. Without this, clang may emit FMA on
 * the scalar side while the NEON intrinsics here stay non-fused, producing
 * a 1-ulp divergence that trips the bit-exactness test. */
#pragma STDC FP_CONTRACT OFF

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

/* Scalar helper exposed by spherical_harmonics.c. */
extern double irrep_legendre_assoc(int l, int m, double x);

/* Scalar per-edge kernel (fallback for the tail, and reference). */
extern void irrep_sph_harm_cart_all(int l_max, double *out, const double r_hat[3]);

/* Pole threshold — kept identical to the scalar path so lane-wise branches
 * pick up the same behaviour. */
#define IRREP_SH_POLE_EPS 1e-14

/* Fill a triangular grid P[l * stride + m] = P_l^m(z) for 0 ≤ m ≤ l ≤ l_max,
 * Condon-Shortley phase. Reproduces the exact arithmetic of
 * `irrep_legendre_assoc` per-(l, m) so results are bit-identical, but
 * amortises the three-term recurrence across all l at each fixed m —
 * O(l_max²) work instead of the O(l_max³) the per-call approach costs.
 *
 * Entries with m > l are left uninitialised; callers must not read them. */
static void legendre_grid_(int l_max, double z, int stride, double *P) {
    if (z < -1.0) z = -1.0;
    if (z >  1.0) z =  1.0;
    double somx2 = sqrt((1.0 - z) * (1.0 + z));

    for (int m = 0; m <= l_max; ++m) {
        /* P_m^m = (-1)^m · (2m-1)!! · somx2^m — computed exactly as the
         * scalar `irrep_legendre_assoc` does: two separate loops (factorial
         * and power) then signed-multiply, to guarantee the same op order. */
        double pmm;
        if (m == 0) {
            pmm = 1.0;
        } else {
            double fact = 1.0;
            for (int k = 1; k <= 2 * m - 1; k += 2) fact *= (double)k;
            double pow_sx = 1.0;
            for (int k = 0; k < m; ++k) pow_sx *= somx2;
            pmm = ((m & 1) ? -1.0 : 1.0) * fact * pow_sx;
        }
        P[m * stride + m] = pmm;
        if (m + 1 > l_max) continue;

        /* P_{m+1}^m = z · (2m+1) · P_m^m */
        double p_m1m = z * (2.0 * (double)m + 1.0) * pmm;
        P[(m + 1) * stride + m] = p_m1m;

        /* Three-term recurrence in l at fixed m — identical to the scalar. */
        double p_lm2 = pmm;
        double p_lm1 = p_m1m;
        for (int l = m + 2; l <= l_max; ++l) {
            double val = ((2.0 * l - 1.0) * z * p_lm1
                          - (double)(l + m - 1) * p_lm2) / (double)(l - m);
            p_lm2 = p_lm1;
            p_lm1 = val;
            P[l * stride + m] = val;
        }
    }
}

/* Evaluate sph-harm cart for two edges at once, writing (l_max+1)² doubles
 * into each of out_a and out_b. */
static void sh_cart_all_pair_(int l_max,
                              const double *rhat_a, const double *rhat_b,
                              double *out_a, double *out_b) {
    const double xa = rhat_a[0], ya = rhat_a[1], za = rhat_a[2];
    const double xb = rhat_b[0], yb = rhat_b[1], zb = rhat_b[2];

    /* r_xy, cp, sp per lane — scalar form matches what the scalar kernel
     * would do element-by-element; we just keep the per-lane values in a
     * float64x2_t for the downstream multiplies. */
    const double r_xy_a = sqrt(xa * xa + ya * ya);
    const double r_xy_b = sqrt(xb * xb + yb * yb);
    const double cp_a = (r_xy_a > IRREP_SH_POLE_EPS) ? xa / r_xy_a : 1.0;
    const double sp_a = (r_xy_a > IRREP_SH_POLE_EPS) ? ya / r_xy_a : 0.0;
    const double cp_b = (r_xy_b > IRREP_SH_POLE_EPS) ? xb / r_xy_b : 1.0;
    const double sp_b = (r_xy_b > IRREP_SH_POLE_EPS) ? yb / r_xy_b : 0.0;

    const float64x2_t vcp = { cp_a, cp_b };
    const float64x2_t vsp = { sp_a, sp_b };

    /* Pre-compute the full triangular Legendre grid per lane — one sweep
     * amortises the three-term recurrence across all (l, m) pairs. Replaces
     * what was O(l_max³) scalar calls with O(l_max²) straight-line work;
     * bit-exact against per-call `irrep_legendre_assoc` by construction
     * (same arithmetic, same order, FP_CONTRACT OFF in both). */
    const int stride = IRREP_L_MAX + 1;
    double Pa[(IRREP_L_MAX + 1) * (IRREP_L_MAX + 1)];
    double Pb[(IRREP_L_MAX + 1) * (IRREP_L_MAX + 1)];
    legendre_grid_(l_max, za, stride, Pa);
    legendre_grid_(l_max, zb, stride, Pb);

    for (int l = 0; l <= l_max; ++l) {
        /* m = 0: out[block_base + l] = N0 · P_l^0(z) per lane.
         * N0 is lane-independent so we scalar-multiply. */
        const double N0 = sqrt((2.0 * l + 1.0) / (4.0 * M_PI));
        const double pla0 = Pa[l * stride + 0];
        const double plb0 = Pb[l * stride + 0];
        out_a[l * l + l] = N0 * pla0;
        out_b[l * l + l] = N0 * plb0;
        if (l == 0) continue;

        /* Chebyshev state per edge: (cos_{m-1}, sin_{m-1}).
         * Vector form: lane 0 = edge a, lane 1 = edge b. */
        float64x2_t vcos_prev = vdupq_n_f64(1.0);
        float64x2_t vsin_prev = vdupq_n_f64(0.0);
        double ratio_inv = 1.0;

        for (int m = 1; m <= l; ++m) {
            /* Non-fused mul/sub/add to match clang's scalar emission for
             * `cp*cos_prev - sp*sin_prev`; a fused `fmsub` would differ in
             * the last bit on some inputs. Order of operations mirrors the
             * scalar code exactly. */
            float64x2_t vcp_cos   = vmulq_f64(vcp, vcos_prev);
            float64x2_t vsp_sin   = vmulq_f64(vsp, vsin_prev);
            float64x2_t vcos_m    = vsubq_f64(vcp_cos, vsp_sin);

            float64x2_t vsp_cos   = vmulq_f64(vsp, vcos_prev);
            float64x2_t vcp_sin   = vmulq_f64(vcp, vsin_prev);
            float64x2_t vsin_m    = vaddq_f64(vsp_cos, vcp_sin);

            vcos_prev = vcos_m;
            vsin_prev = vsin_m;

            ratio_inv /= (double)(l - m + 1) * (double)(l + m);
            const double N   = sqrt((2.0 * l + 1.0) / (4.0 * M_PI) * ratio_inv);
            const double pla = Pa[l * stride + m];
            const double plb = Pb[l * stride + m];
            const double sgn = (m & 1) ? -1.0 : 1.0;

            /* Scalar associativity:
             *   ((sgn · √2) · N · plm) · cos_m
             * Reproduce bit-for-bit: build sgn·√2 scalar, broadcast, then
             * compose the three per-lane multiplies in order. */
            const double sgn_sqrt2   = sgn * 1.4142135623730951;
            const double sgn_sqrt2_N = sgn_sqrt2 * N;
            const float64x2_t vK = vdupq_n_f64(sgn_sqrt2_N);
            const float64x2_t vP = { pla, plb };
            const float64x2_t vcoeff = vmulq_f64(vK, vP);

            const float64x2_t v_plus  = vmulq_f64(vcoeff, vcos_m);
            const float64x2_t v_minus = vmulq_f64(vcoeff, vsin_m);

            out_a[l * l + l + m] = vgetq_lane_f64(v_plus,  0);
            out_b[l * l + l + m] = vgetq_lane_f64(v_plus,  1);
            out_a[l * l + l - m] = vgetq_lane_f64(v_minus, 0);
            out_b[l * l + l - m] = vgetq_lane_f64(v_minus, 1);
        }
    }
}

void irrep_sph_harm_cart_all_batch_neon(int l_max, size_t N,
                                        const double *r_hats, double *out) {
    if (l_max < 0 || !r_hats || !out) return;
    const int block = (l_max + 1) * (l_max + 1);
    size_t e = 0;
    for (; e + 2 <= N; e += 2) {
        sh_cart_all_pair_(l_max,
                          r_hats + e       * 3,
                          r_hats + (e + 1) * 3,
                          out   + e       * (size_t)block,
                          out   + (e + 1) * (size_t)block);
    }
    /* Tail: one edge, fall through to scalar. */
    if (e < N) {
        irrep_sph_harm_cart_all(l_max, out + e * (size_t)block,
                                 r_hats + e * 3);
    }
}

#else
/* Non-aarch64: keep the translation unit non-empty to silence pedantic
 * "file has no symbols" warnings under some configurations. */
typedef int irrep_sph_harm_neon_stub_t;
#endif
