/* SPDX-License-Identifier: MIT */
/* M3: spherical harmonics — complex, real, cartesian, gradient, basis change.
 *
 * Conventions: Condon-Shortley phase, orthonormal on S² under ∫ dΩ. Real SH
 * follow the e3nn / Wikipedia sign convention (Y_{l,+1} ∝ +x at the equator).
 * Associated Legendre via the stable three-term recurrence (Numerical Recipes
 * 3e §6.7; Limpanuparb & Milthorpe 2014).
 *
 * M3 gradient is a centered finite-difference on the normalized input. The
 * analytic solid-harmonic form lands in M10 alongside the SIMD work.
 */

#include <complex.h>
#include <math.h>
#include <stddef.h>

#include <irrep/spherical_harmonics.h>
#include <irrep/types.h>

#include "internal/dispatch.h"

/* The NEON batch kernel reproduces the scalar algorithm lane-by-lane and
 * compares bit-exactly. Disabling fused-multiply-add contraction here keeps
 * clang from pattern-matching `a*b - c*d` to `fmsub`, which would give a
 * last-bit difference from the NEON path's separate mul/sub ops. */
#pragma STDC FP_CONTRACT OFF

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define UNUSED(x) ((void)(x))

/* -------------------------------------------------------------------------- *
 * Associated Legendre polynomial (Condon-Shortley phase)                     *
 * -------------------------------------------------------------------------- */

double irrep_legendre_assoc(int l, int m, double x) {
    if (l < 0)
        return 0.0;
    int am = m < 0 ? -m : m;
    if (am > l)
        return 0.0;
    if (x < -1.0)
        x = -1.0;
    if (x > 1.0)
        x = 1.0;

    /* P_m^m(x) = (-1)^m (2m-1)!! (1 - x²)^{m/2}  (Condon-Shortley). */
    double somx2 = sqrt((1.0 - x) * (1.0 + x)); /* sin θ, ≥ 0 */
    double pmm = 1.0;
    if (am > 0) {
        double fact = 1.0;
        for (int k = 1; k <= 2 * am - 1; k += 2)
            fact *= (double)k;
        double pow_sx = 1.0;
        for (int k = 0; k < am; ++k)
            pow_sx *= somx2;
        pmm = ((am & 1) ? -1.0 : 1.0) * fact * pow_sx;
    }

    double plm;
    if (l == am) {
        plm = pmm;
    } else {
        /* P_{m+1}^m = x (2m + 1) P_m^m. */
        double p_m1m = x * (2.0 * am + 1.0) * pmm;
        if (l == am + 1) {
            plm = p_m1m;
        } else {
            /* Three-term: P_l^m = ((2l − 1) x P_{l−1}^m − (l + m − 1) P_{l−2}^m) / (l − m). */
            double p_lm2 = pmm;
            double p_lm1 = p_m1m;
            for (int ll = am + 2; ll <= l; ++ll) {
                double val = ((2.0 * ll - 1.0) * x * p_lm1 - (double)(ll + am - 1) * p_lm2) /
                             (double)(ll - am);
                p_lm2 = p_lm1;
                p_lm1 = val;
            }
            plm = p_lm1;
        }
    }

    if (m < 0) {
        /* P_l^{−|m|} = (−1)^{|m|} (l − |m|)! / (l + |m|)! · P_l^{|m|}. */
        double ratio = 1.0;
        for (int k = l - am + 1; k <= l + am; ++k)
            ratio *= (double)k;
        double sign = (am & 1) ? -1.0 : 1.0;
        plm = sign * plm / ratio;
    }
    return plm;
}

/* -------------------------------------------------------------------------- *
 * Complex and real-polar spherical harmonics                                 *
 * -------------------------------------------------------------------------- */

static double norm_lm_(int l, int am) {
    double ratio = 1.0;
    for (int k = l - am + 1; k <= l + am; ++k)
        ratio *= (double)k;
    return sqrt((2.0 * l + 1.0) / (4.0 * M_PI) / ratio);
}

double _Complex irrep_sph_harm(int l, int m, double theta, double phi) {
    int am = m < 0 ? -m : m;
    if (l < 0 || am > l)
        return 0.0;
    double p = irrep_legendre_assoc(l, am, cos(theta));
    double N = norm_lm_(l, am);
    double c = cos((double)am * phi), s = sin((double)am * phi);
    double _Complex base = N * p * (c + I * s);
    if (m >= 0)
        return base;
    double sign = (am & 1) ? -1.0 : 1.0;
    return sign * conj(base);
}

double irrep_sph_harm_real(int l, int m, double theta, double phi) {
    int am = m < 0 ? -m : m;
    if (l < 0 || am > l)
        return 0.0;
    double p = irrep_legendre_assoc(l, am, cos(theta));
    double N = norm_lm_(l, am);
    if (m == 0)
        return N * p;
    double       sign = (am & 1) ? -1.0 : 1.0;
    const double sqrt2 = 1.4142135623730951;
    if (m > 0)
        return sign * sqrt2 * N * p * cos((double)am * phi);
    /*      */ return sign * sqrt2 * N * p * sin((double)am * phi);
}

/* -------------------------------------------------------------------------- *
 * Cartesian real spherical harmonics (hot path)                              *
 * -------------------------------------------------------------------------- */

void irrep_sph_harm_cart(int l, double *out, const double r_hat[3]) {
    if (l < 0 || !out || !r_hat)
        return;
    double x = r_hat[0], y = r_hat[1], z = r_hat[2];

    /* m = 0: Y_{l,0} = N(l,0) P_l(z). */
    double N0 = sqrt((2.0 * l + 1.0) / (4.0 * M_PI));
    out[l] = N0 * irrep_legendre_assoc(l, 0, z);
    if (l == 0)
        return;

    /* Trig recurrence: cos(mφ), sin(mφ) from (x, y).
     * At the pole (r_xy ≈ 0) φ is undefined; any Y_{l, m ≠ 0} contains a
     * (sin θ)^|m| prefactor that vanishes at the pole anyway, so the choice
     * (cp, sp) = (1, 0) gives the correct limit. The threshold below
     * IRREP_SH_POLE_EPS isolates the divide and is set well inside the
     * region where double-precision loss on x/r_xy would matter. */
    const double IRREP_SH_POLE_EPS = 1e-14;
    double       r_xy = sqrt(x * x + y * y);
    double       cp, sp;
    if (r_xy > IRREP_SH_POLE_EPS) {
        cp = x / r_xy;
        sp = y / r_xy;
    } else {
        cp = 1.0;
        sp = 0.0;
    }

    double       cos_prev = 1.0, sin_prev = 0.0;
    const double sqrt2 = 1.4142135623730951;
    double       ratio_inv = 1.0; /* (l − m)! / (l + m)! */

    for (int m = 1; m <= l; ++m) {
        double cos_m = cp * cos_prev - sp * sin_prev;
        double sin_m = sp * cos_prev + cp * sin_prev;
        cos_prev = cos_m;
        sin_prev = sin_m;

        ratio_inv /= (double)(l - m + 1) * (double)(l + m);

        double plm = irrep_legendre_assoc(l, m, z);
        double N = sqrt((2.0 * l + 1.0) / (4.0 * M_PI) * ratio_inv);
        double sgn = (m & 1) ? -1.0 : 1.0;

        out[l + m] = sgn * sqrt2 * N * plm * cos_m; /* +m stored at index l + m */
        out[l - m] = sgn * sqrt2 * N * plm * sin_m; /* −m stored at index l − m */
    }
}

void irrep_sph_harm_cart_all(int l_max, double *out, const double r_hat[3]) {
    if (l_max < 0 || !out || !r_hat)
        return;
    int offset = 0;
    for (int l = 0; l <= l_max; ++l) {
        irrep_sph_harm_cart(l, out + offset, r_hat);
        offset += 2 * l + 1;
    }
}

/* -------------------------------------------------------------------------- *
 * Gradient of cartesian real surface SH.                                     *
 * Layout: out[axis * (2l+1) + (m + l)] for axis ∈ {0:x, 1:y, 2:z}.           *
 *                                                                            *
 * Analytic via solid harmonics: on the unit sphere,                          *
 *                                                                            *
 *   ∇_i Y_{l,m}(r̂) = ∇_i R_{l,m}(r)|_{|r|=1} − l · r̂_i · Y_{l,m}(r̂).        *
 *                                                                            *
 * Both R and ∇R are computed analytically by `irrep/solid_harmonics.h`.      *
 * -------------------------------------------------------------------------- */

#include <irrep/solid_harmonics.h>

void irrep_sph_harm_cart_grad(int l, double *out, const double r_hat[3]) {
    if (l < 0 || !out)
        return;
    int    d = 2 * l + 1;

    double Y[2 * IRREP_L_MAX + 1];
    irrep_sph_harm_cart(l, Y, r_hat);

    double dR[3 * (2 * IRREP_L_MAX + 1)];
    irrep_solid_harm_cart_grad(l, dR, r_hat);

    for (int axis = 0; axis < 3; ++axis) {
        for (int i = 0; i < d; ++i) {
            out[axis * d + i] = dR[axis * d + i] - (double)l * r_hat[axis] * Y[i];
        }
    }
}

/* -------------------------------------------------------------------------- *
 * Complex → real basis change matrix                                         *
 *                                                                            *
 * Y_real = U · Y_complex, with rows and columns both indexed m + l (so       *
 * m = −l sits at index 0, m = +l at index 2l). The convention matches e3nn  *
 * and Wikipedia "Real-valued spherical harmonics" (Condon-Shortley Y_l^m):   *
 *                                                                            *
 *   Y_{l, m}^{real} = Y_l^0                           if m = 0               *
 *                   = (Y_l^{-m} + (-1)^m Y_l^m) / √2   if m > 0              *
 *                   = i (Y_l^{-|m|} − (-1)^{|m|} Y_l^{|m|}) · (−1) / √2      *
 *                     for m < 0; rearranged below.                           *
 * -------------------------------------------------------------------------- */

void irrep_sph_harm_complex_to_real(int l, double _Complex *out) {
    if (l < 0)
        return;
    const int d = 2 * l + 1;
    for (int i = 0; i < d * d; ++i)
        out[i] = 0.0;
    out[l * d + l] = 1.0; /* m_r = 0 */

    const double s2i = 1.0 / 1.4142135623730951;
    for (int m = 1; m <= l; ++m) {
        double sign = (m & 1) ? -1.0 : 1.0;

        /* Y_{l,+m}^real = (-1)^m / √2 · Y_l^m + 1 / √2 · Y_l^{-m}. */
        out[(l + m) * d + (l + m)] = sign * s2i;
        out[(l + m) * d + (l - m)] = s2i;

        /* Y_{l,-m}^real = i · (-1)^{m+1} / √2 · Y_l^m + i / √2 · Y_l^{-m}. */
        out[(l - m) * d + (l + m)] = -I * sign * s2i;
        out[(l - m) * d + (l - m)] = I * s2i;
    }
}

/* -------------------------------------------------------------------------- *
 * Single-precision wrappers (real SIMD kernels in M10/M11)                   *
 * -------------------------------------------------------------------------- */

void irrep_sph_harm_cart_f32(int l, float *out, const float r_hat[3]) {
    if (l < 0)
        return;
    double rd[3] = {(double)r_hat[0], (double)r_hat[1], (double)r_hat[2]};
    double tmp[2 * IRREP_L_MAX + 3];
    irrep_sph_harm_cart(l, tmp, rd);
    for (int i = 0; i < 2 * l + 1; ++i)
        out[i] = (float)tmp[i];
}

void irrep_sph_harm_cart_all_f32(int l_max, float *out, const float r_hat[3]) {
    if (l_max < 0)
        return;
    double rd[3] = {(double)r_hat[0], (double)r_hat[1], (double)r_hat[2]};
    double tmp[(IRREP_L_MAX + 1) * (IRREP_L_MAX + 1)];
    irrep_sph_harm_cart_all(l_max, tmp, rd);
    int total = (l_max + 1) * (l_max + 1);
    for (int i = 0; i < total; ++i)
        out[i] = (float)tmp[i];
}

void irrep_sph_harm_cart_all_batch_scalar(int l_max, size_t N, const double *r_hats, double *out) {
    if (l_max < 0)
        return;
    int block = (l_max + 1) * (l_max + 1);
    for (size_t i = 0; i < N; ++i) {
        irrep_sph_harm_cart_all(l_max, out + i * (size_t)block, r_hats + i * 3);
    }
}

void irrep_sph_harm_cart_all_batch(int l_max, size_t N, const double *r_hats, double *out) {
    irrep_dispatch_get()->sph_harm_cart_all_batch(l_max, N, r_hats, out);
}

/* Batched gradient: for each edge, for each l = 0..l_max, the gradient block
 * is 3 · (2l+1) reals. The per-edge output is
 *
 *     out[edge, axis, l, m + l]
 *       → out[(edge · 3 + axis) · (l_max+1)² + (l² + (m + l))]
 *
 * where each `l`-group of `2l + 1` entries is contiguous within an axis plane.
 * This matches the layout produced by stacking irrep_sph_harm_cart_grad
 * across l (per axis). */
void irrep_sph_harm_cart_all_grad_batch(int l_max, size_t N, const double *r_hats, double *out) {
    if (l_max < 0 || !out || !r_hats)
        return;
    int    block = (l_max + 1) * (l_max + 1);
    double tmp_grad[3 * (2 * IRREP_L_MAX + 1)];
    for (size_t e = 0; e < N; ++e) {
        const double *rhat = r_hats + e * 3;
        double       *base = out + e * 3 * (size_t)block;
        int           off = 0;
        for (int l = 0; l <= l_max; ++l) {
            int d = 2 * l + 1;
            irrep_sph_harm_cart_grad(l, tmp_grad, rhat);
            for (int axis = 0; axis < 3; ++axis) {
                for (int i = 0; i < d; ++i) {
                    base[axis * (size_t)block + off + i] = tmp_grad[axis * d + i];
                }
            }
            off += d;
        }
    }
}
