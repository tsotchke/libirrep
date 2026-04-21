/* SPDX-License-Identifier: MIT */
/* AVX2 (x86_64) batch kernel for irrep_sph_harm_cart_all.
 *
 * Mirrors src/spherical_harmonics_neon.c, but four edges per vector
 * instead of two: each __m256d carries lanes {e, e+1, e+2, e+3} and
 * the Chebyshev trig recurrence for cos(mφ) / sin(mφ) runs across them
 * in lockstep. The associated-Legendre grid is evaluated per-lane by
 * the scalar kernel (sidesteps cross-lane branching in the m-recurrence
 * and keeps the output bit-exact against the scalar reference).
 *
 * Tail: N mod 4 ∈ {1, 2, 3} falls through to the scalar per-edge path.
 *
 * Bit-exactness is asserted in tests/test_spherical_harmonics.c by a
 * scalar-vs-AVX2 sweep at l_max ∈ {0..4} over N = 19 edges (tail 3
 * exercised). */

#if (defined(__x86_64__) || defined(_M_X64)) && defined(__AVX2__) && defined(__FMA__)

#include <immintrin.h>
#include <math.h>
#include <stddef.h>

#include <irrep/spherical_harmonics.h>
#include <irrep/types.h>

#include "internal/dispatch.h"

/* Disable FP contraction so the AVX2 path and the scalar path use the same
 * unfused mul/add/sub sequence — otherwise clang's scalar side sometimes
 * emits vfmadd while these intrinsics stay unfused, giving a 1-ulp
 * divergence that trips the bit-exactness test. */
#pragma STDC FP_CONTRACT OFF

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern double irrep_legendre_assoc(int l, int m, double x);
extern void   irrep_sph_harm_cart_all(int l_max, double *out, const double r_hat[3]);

#define IRREP_SH_POLE_EPS 1e-14

/* Associated-Legendre triangle grid P[l * stride + m] = P_l^m(z) for
 * 0 ≤ m ≤ l ≤ l_max. Identical arithmetic to the NEON variant — kept as
 * a free-standing function because each lane needs its own grid. */
static void legendre_grid_(int l_max, double z, int stride, double *P) {
    if (z < -1.0)
        z = -1.0;
    if (z > 1.0)
        z = 1.0;
    double somx2 = sqrt((1.0 - z) * (1.0 + z));

    for (int m = 0; m <= l_max; ++m) {
        double pmm;
        if (m == 0) {
            pmm = 1.0;
        } else {
            double fact = 1.0;
            for (int k = 1; k <= 2 * m - 1; k += 2)
                fact *= (double)k;
            double pow_sx = 1.0;
            for (int k = 0; k < m; ++k)
                pow_sx *= somx2;
            pmm = ((m & 1) ? -1.0 : 1.0) * fact * pow_sx;
        }
        P[m * stride + m] = pmm;
        if (m + 1 > l_max)
            continue;

        double p_m1m = z * (2.0 * (double)m + 1.0) * pmm;
        P[(m + 1) * stride + m] = p_m1m;

        double p_lm2 = pmm;
        double p_lm1 = p_m1m;
        for (int l = m + 2; l <= l_max; ++l) {
            double val =
                ((2.0 * l - 1.0) * z * p_lm1 - (double)(l + m - 1) * p_lm2) / (double)(l - m);
            p_lm2 = p_lm1;
            p_lm1 = val;
            P[l * stride + m] = val;
        }
    }
}

/* Evaluate sph-harm cart for four edges at once, writing (l_max+1)²
 * doubles into each of out_0..out_3. */
static void sh_cart_all_quad_(int l_max, const double *rhat_0, const double *rhat_1,
                              const double *rhat_2, const double *rhat_3, double *out_0,
                              double *out_1, double *out_2, double *out_3) {
    const double x0 = rhat_0[0], y0 = rhat_0[1], z0 = rhat_0[2];
    const double x1 = rhat_1[0], y1 = rhat_1[1], z1 = rhat_1[2];
    const double x2 = rhat_2[0], y2 = rhat_2[1], z2 = rhat_2[2];
    const double x3 = rhat_3[0], y3 = rhat_3[1], z3 = rhat_3[2];

    const double r_xy_0 = sqrt(x0 * x0 + y0 * y0);
    const double r_xy_1 = sqrt(x1 * x1 + y1 * y1);
    const double r_xy_2 = sqrt(x2 * x2 + y2 * y2);
    const double r_xy_3 = sqrt(x3 * x3 + y3 * y3);

    const double cp_0 = (r_xy_0 > IRREP_SH_POLE_EPS) ? x0 / r_xy_0 : 1.0;
    const double sp_0 = (r_xy_0 > IRREP_SH_POLE_EPS) ? y0 / r_xy_0 : 0.0;
    const double cp_1 = (r_xy_1 > IRREP_SH_POLE_EPS) ? x1 / r_xy_1 : 1.0;
    const double sp_1 = (r_xy_1 > IRREP_SH_POLE_EPS) ? y1 / r_xy_1 : 0.0;
    const double cp_2 = (r_xy_2 > IRREP_SH_POLE_EPS) ? x2 / r_xy_2 : 1.0;
    const double sp_2 = (r_xy_2 > IRREP_SH_POLE_EPS) ? y2 / r_xy_2 : 0.0;
    const double cp_3 = (r_xy_3 > IRREP_SH_POLE_EPS) ? x3 / r_xy_3 : 1.0;
    const double sp_3 = (r_xy_3 > IRREP_SH_POLE_EPS) ? y3 / r_xy_3 : 0.0;

    /* AVX2 layout: lane 0..3 = edges 0..3. */
    const __m256d vcp = _mm256_setr_pd(cp_0, cp_1, cp_2, cp_3);
    const __m256d vsp = _mm256_setr_pd(sp_0, sp_1, sp_2, sp_3);

    const int     stride = IRREP_L_MAX + 1;
    double        P0[(IRREP_L_MAX + 1) * (IRREP_L_MAX + 1)];
    double        P1[(IRREP_L_MAX + 1) * (IRREP_L_MAX + 1)];
    double        P2[(IRREP_L_MAX + 1) * (IRREP_L_MAX + 1)];
    double        P3[(IRREP_L_MAX + 1) * (IRREP_L_MAX + 1)];
    legendre_grid_(l_max, z0, stride, P0);
    legendre_grid_(l_max, z1, stride, P1);
    legendre_grid_(l_max, z2, stride, P2);
    legendre_grid_(l_max, z3, stride, P3);

    for (int l = 0; l <= l_max; ++l) {
        const double N0 = sqrt((2.0 * l + 1.0) / (4.0 * M_PI));
        out_0[l * l + l] = N0 * P0[l * stride + 0];
        out_1[l * l + l] = N0 * P1[l * stride + 0];
        out_2[l * l + l] = N0 * P2[l * stride + 0];
        out_3[l * l + l] = N0 * P3[l * stride + 0];
        if (l == 0)
            continue;

        __m256d vcos_prev = _mm256_set1_pd(1.0);
        __m256d vsin_prev = _mm256_set1_pd(0.0);
        double  ratio_inv = 1.0;

        for (int m = 1; m <= l; ++m) {
            /* cos(mφ) = cp · cos((m-1)φ) − sp · sin((m-1)φ)
             * sin(mφ) = sp · cos((m-1)φ) + cp · sin((m-1)φ)
             * Non-fused sequence for bit-exactness against the scalar path. */
            __m256d vcp_cos = _mm256_mul_pd(vcp, vcos_prev);
            __m256d vsp_sin = _mm256_mul_pd(vsp, vsin_prev);
            __m256d vcos_m = _mm256_sub_pd(vcp_cos, vsp_sin);

            __m256d vsp_cos = _mm256_mul_pd(vsp, vcos_prev);
            __m256d vcp_sin = _mm256_mul_pd(vcp, vsin_prev);
            __m256d vsin_m = _mm256_add_pd(vsp_cos, vcp_sin);

            vcos_prev = vcos_m;
            vsin_prev = vsin_m;

            ratio_inv /= (double)(l - m + 1) * (double)(l + m);
            const double  N = sqrt((2.0 * l + 1.0) / (4.0 * M_PI) * ratio_inv);
            const double  sgn = (m & 1) ? -1.0 : 1.0;
            const double  sgn_sqrt2_N = sgn * 1.4142135623730951 * N;

            const __m256d vK = _mm256_set1_pd(sgn_sqrt2_N);
            const __m256d vP = _mm256_setr_pd(P0[l * stride + m], P1[l * stride + m],
                                              P2[l * stride + m], P3[l * stride + m]);
            const __m256d vcoeff = _mm256_mul_pd(vK, vP);
            const __m256d v_plus = _mm256_mul_pd(vcoeff, vcos_m);
            const __m256d v_minus = _mm256_mul_pd(vcoeff, vsin_m);

            double        plus_buf[4], minus_buf[4];
            _mm256_storeu_pd(plus_buf, v_plus);
            _mm256_storeu_pd(minus_buf, v_minus);

            out_0[l * l + l + m] = plus_buf[0];
            out_1[l * l + l + m] = plus_buf[1];
            out_2[l * l + l + m] = plus_buf[2];
            out_3[l * l + l + m] = plus_buf[3];

            out_0[l * l + l - m] = minus_buf[0];
            out_1[l * l + l - m] = minus_buf[1];
            out_2[l * l + l - m] = minus_buf[2];
            out_3[l * l + l - m] = minus_buf[3];
        }
    }
}

void irrep_sph_harm_cart_all_batch_avx2(int l_max, size_t N, const double *r_hats, double *out) {
    if (l_max < 0 || !r_hats || !out)
        return;
    const int block = (l_max + 1) * (l_max + 1);
    size_t    e = 0;
    for (; e + 4 <= N; e += 4) {
        sh_cart_all_quad_(l_max, r_hats + e * 3, r_hats + (e + 1) * 3, r_hats + (e + 2) * 3,
                          r_hats + (e + 3) * 3, out + e * (size_t)block,
                          out + (e + 1) * (size_t)block, out + (e + 2) * (size_t)block,
                          out + (e + 3) * (size_t)block);
    }
    /* Tail (1..3 edges): scalar per-edge path. */
    for (; e < N; ++e) {
        irrep_sph_harm_cart_all(l_max, out + e * (size_t)block, r_hats + e * 3);
    }
}

#else
/* Non-AVX2 build: keep the translation unit non-empty. */
typedef int irrep_sph_harm_avx2_stub_t;
#endif
