/* SPDX-License-Identifier: MIT */
/* M8: radial basis functions and cutoff envelopes for equivariant GNNs. */

#include <math.h>
#include <stddef.h>

#include <irrep/radial.h>

#include "internal/dispatch.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define UNUSED(x) ((void)(x))

/* DimeNet Bessel RBF (Klicpera et al. 2020):
 *
 *   φ_n(r) = √(2/r_cut) · sin(n π r / r_cut) / r        for 0 < r < r_cut
 *          = √(2/r_cut) · (n π / r_cut)                 at r = 0 (sinc limit)
 *          = 0                                           for r ≥ r_cut.
 *
 * Orthonormal on [0, r_cut] with the r² volume element. */
double irrep_rbf_bessel(int n, double r, double r_cut) {
    if (n <= 0 || r_cut <= 0.0 || r < 0.0 || r >= r_cut)
        return 0.0;
    double norm = sqrt(2.0 / r_cut);
    double arg = (double)n * M_PI * r / r_cut;
    if (r < 1e-12)
        return norm * (double)n * M_PI / r_cut;
    return norm * sin(arg) / r;
}

void irrep_rbf_bessel_all(int n_max, double *out, double r, double r_cut) {
    if (!out || n_max <= 0)
        return;
    for (int n = 1; n <= n_max; ++n) {
        out[n - 1] = irrep_rbf_bessel(n, r, r_cut);
    }
}

/* Derivative d/dr φ_n(r).
 *
 *   φ_n(r) = C · sin(a r) / r        with C = √(2/r_cut), a = n π / r_cut,
 *
 * so
 *
 *   φ_n'(r) = C · [ a cos(a r) / r − sin(a r) / r² ]
 *
 * at small r the two terms are the leading entries of a Taylor series that
 * cancels to ~(a r)² / 6 — by r ≈ 1e-7 · a ~ 1, the naive form has lost
 * all double-precision significance to cancellation. Switch to a
 * three-term Taylor expansion
 *
 *   φ_n'(r) = C · (−a³ r / 3) · (1 − (ar)² / 10 + (ar)⁴ / 280 − …)
 *
 * whenever |a r| < 1e-3. The truncated fourth term there is ≲ 1e-14, below
 * double-precision round-off on the leading coefficient. At |a r| ≥ 1e-3
 * the full formula loses at most ~6 significant digits to cancellation,
 * still well above our 1e-10 tolerance. */
double irrep_rbf_bessel_d(int n, double r, double r_cut) {
    if (n <= 0 || r_cut <= 0.0 || r < 0.0 || r >= r_cut)
        return 0.0;
    double C = sqrt(2.0 / r_cut);
    double a = (double)n * M_PI / r_cut;
    double ar = a * r;
    if (ar < 1e-3) {
        double ar2 = ar * ar;
        return -C * a * a * a * r / 3.0 * (1.0 - ar2 / 10.0 + ar2 * ar2 / 280.0);
    }
    return C * (a * cos(ar) / r - sin(ar) / (r * r));
}

void irrep_rbf_bessel_d_all(int n_max, double *out, double r, double r_cut) {
    if (!out || n_max <= 0)
        return;
    for (int n = 1; n <= n_max; ++n) {
        out[n - 1] = irrep_rbf_bessel_d(n, r, r_cut);
    }
}

void irrep_rbf_bessel_d_batch(int n, size_t N, const double *r, double r_cut, double *out) {
    for (size_t i = 0; i < N; ++i)
        out[i] = irrep_rbf_bessel_d(n, r[i], r_cut);
}

/* Gaussian RBF: φ(r) = exp(−(r − μ)² / (2σ²)). */
double irrep_rbf_gaussian(double r, double mu, double sigma) {
    if (sigma <= 0.0)
        return 0.0;
    double d = r - mu;
    return exp(-0.5 * d * d / (sigma * sigma));
}

void irrep_rbf_gaussian_grid(int n, double *out, double r, double r_min, double r_max,
                             double sigma) {
    if (!out || n <= 0 || sigma <= 0.0)
        return;
    if (n == 1) {
        out[0] = irrep_rbf_gaussian(r, r_min, sigma);
        return;
    }
    double step = (r_max - r_min) / (double)(n - 1);
    for (int i = 0; i < n; ++i) {
        double mu = r_min + (double)i * step;
        out[i] = irrep_rbf_gaussian(r, mu, sigma);
    }
}

/* Cosine cutoff: c(r) = (1/2)(1 + cos(π r / r_cut)) for r < r_cut, else 0. */
double irrep_cutoff_cosine(double r, double r_cut) {
    if (r_cut <= 0.0 || r < 0.0 || r >= r_cut)
        return 0.0;
    return 0.5 * (1.0 + cos(M_PI * r / r_cut));
}

double irrep_cutoff_cosine_d(double r, double r_cut) {
    if (r_cut <= 0.0 || r < 0.0 || r >= r_cut)
        return 0.0;
    return -0.5 * M_PI / r_cut * sin(M_PI * r / r_cut);
}

/* NequIP polynomial cutoff (Batzner et al. 2022):
 *
 *   f_p(u) = 1 − (p+1)(p+2)/2 · u^p
 *              + p (p+2)    · u^{p+1}
 *              − p (p+1)/2  · u^{p+2},        u = r / r_cut,
 *
 * smooth at u = 1 through order ~ p. Derivative simplifies to
 *
 *   f_p'(r) = −p(p+1)(p+2) / (2 r_cut) · u^{p−1} · (1 − u)². */
double irrep_cutoff_polynomial(double r, double r_cut, int p) {
    if (r_cut <= 0.0 || r < 0.0 || r >= r_cut || p < 1)
        return 0.0;
    double u = r / r_cut;
    /* u^p via repeated multiply — integer p, so this is both faster and
     * bit-reproducible across the scalar and SIMD batched kernels. */
    double a = u;
    for (int k = 1; k < p; ++k)
        a *= u;
    double b = a * u; /* u^{p+1} */
    double c = b * u; /* u^{p+2} */
    double c1 = 0.5 * (double)(p + 1) * (double)(p + 2);
    double c2 = (double)p * (double)(p + 2);
    double c3 = 0.5 * (double)p * (double)(p + 1);
    return 1.0 - c1 * a + c2 * b - c3 * c;
}

double irrep_cutoff_polynomial_d(double r, double r_cut, int p) {
    if (r_cut <= 0.0 || r < 0.0 || r >= r_cut || p < 1)
        return 0.0;
    double u = r / r_cut;
    double om = 1.0 - u;
    double pre = -(double)p * (double)(p + 1) * (double)(p + 2) / (2.0 * r_cut);
    double upm1 = 1.0;
    for (int k = 0; k < p - 1; ++k)
        upm1 *= u;
    return pre * upm1 * om * om;
}

/* -------------------------------------------------------------------------- *
 * Batched variants. Public entry points trampoline through the runtime SIMD  *
 * dispatch table; scalar implementations below are the fallback.             *
 * -------------------------------------------------------------------------- */

void irrep_rbf_bessel_batch_scalar(int n, size_t N, const double *r, double r_cut, double *out) {
    for (size_t i = 0; i < N; ++i)
        out[i] = irrep_rbf_bessel(n, r[i], r_cut);
}

void irrep_cutoff_cosine_batch_scalar(size_t N, const double *r, double r_cut, double *out) {
    for (size_t i = 0; i < N; ++i)
        out[i] = irrep_cutoff_cosine(r[i], r_cut);
}

void irrep_cutoff_polynomial_batch_scalar(size_t N, const double *r, double r_cut, int p,
                                          double *out) {
    for (size_t i = 0; i < N; ++i)
        out[i] = irrep_cutoff_polynomial(r[i], r_cut, p);
}

void irrep_cutoff_cosine_d_batch_scalar(size_t N, const double *r, double r_cut, double *out) {
    for (size_t i = 0; i < N; ++i)
        out[i] = irrep_cutoff_cosine_d(r[i], r_cut);
}

void irrep_cutoff_polynomial_d_batch_scalar(size_t N, const double *r, double r_cut, int p,
                                            double *out) {
    for (size_t i = 0; i < N; ++i)
        out[i] = irrep_cutoff_polynomial_d(r[i], r_cut, p);
}

void irrep_rbf_bessel_batch(int n, size_t N, const double *r, double r_cut, double *out) {
    irrep_dispatch_get()->rbf_bessel_batch(n, N, r, r_cut, out);
}

void irrep_cutoff_cosine_batch(size_t N, const double *r, double r_cut, double *out) {
    irrep_dispatch_get()->cutoff_cosine_batch(N, r, r_cut, out);
}

void irrep_cutoff_polynomial_batch(size_t N, const double *r, double r_cut, int p, double *out) {
    irrep_dispatch_get()->cutoff_polynomial_batch(N, r, r_cut, p, out);
}

void irrep_cutoff_cosine_d_batch(size_t N, const double *r, double r_cut, double *out) {
    irrep_dispatch_get()->cutoff_cosine_d_batch(N, r, r_cut, out);
}

void irrep_cutoff_polynomial_d_batch(size_t N, const double *r, double r_cut, int p, double *out) {
    irrep_dispatch_get()->cutoff_polynomial_d_batch(N, r, r_cut, p, out);
}
