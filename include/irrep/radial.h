/* SPDX-License-Identifier: MIT */
/** @file radial.h
 *  @brief Radial basis functions and smooth cutoff envelopes for equivariant
 *         GNNs. Paired with @ref spherical_harmonics.h in every NequIP /
 *         MACE / Allegro message-passing layer.
 */
#ifndef IRREP_RADIAL_H
#define IRREP_RADIAL_H

#include <stddef.h>

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Bessel RBF of degree @p n on `[0, r_cut]`:
 *         `φ_n(r) = √(2 / r_cut) · sin(n π r / r_cut) / r`.
 *  Returns 0 for `r ≥ r_cut` and the regular sinc limit at `r → 0`.
 *  Orthonormal with the `r²` volume element. */
IRREP_API double irrep_rbf_bessel    (int n, double r, double r_cut);

/** @brief Fills `out[n-1] = irrep_rbf_bessel(n, r, r_cut)` for `n = 1..n_max`. */
IRREP_API void   irrep_rbf_bessel_all(int n_max, double *out, double r, double r_cut);

/** @brief Derivative `d/dr φ_n(r)` of the Bessel RBF. Zero at `r = 0`
 *         (the RBF is smooth through the origin with `φ_n'(0) = 0`) and
 *         zero for `r ≥ r_cut`. Small-r guarded by a Taylor expansion so
 *         the naive `a cos(ar)/r − sin(ar)/r²` doesn't catastrophically
 *         cancel. */
IRREP_API double irrep_rbf_bessel_d   (int n, double r, double r_cut);

/** @brief Fills `out[n-1] = irrep_rbf_bessel_d(n, r, r_cut)` for `n = 1..n_max`. */
IRREP_API void   irrep_rbf_bessel_d_all(int n_max, double *out, double r, double r_cut);

/** @brief Batched derivative — @p N radii in one call. */
IRREP_API void   irrep_rbf_bessel_d_batch(int n, size_t N, const double *r,
                                          double r_cut, double *out);

/** @brief Gaussian RBF: `φ(r) = exp(−(r − μ)² / (2 σ²))`. */
IRREP_API double irrep_rbf_gaussian     (double r, double mu, double sigma);

/** @brief Fills @p out with @p n Gaussians whose centres are evenly spaced
 *         on `[r_min, r_max]` and whose common width is @p sigma. */
IRREP_API void   irrep_rbf_gaussian_grid(int n, double *out, double r,
                                         double r_min, double r_max, double sigma);

/** @brief Cosine cutoff: `(1 + cos(π r / r_cut)) / 2` on `[0, r_cut]`, else 0. */
IRREP_API double irrep_cutoff_cosine    (double r, double r_cut);

/** @brief NequIP polynomial cutoff (smooth to order `p` at `r = r_cut`). */
IRREP_API double irrep_cutoff_polynomial(double r, double r_cut, int p);

/** @brief Derivative of the cosine cutoff (for force evaluation). */
IRREP_API double irrep_cutoff_cosine_d    (double r, double r_cut);

/** @brief Derivative of the polynomial cutoff. */
IRREP_API double irrep_cutoff_polynomial_d(double r, double r_cut, int p);

/** @name Batched variants — runtime SIMD dispatch.
 *  Each evaluates an @c N -vector of radii in one call. NEON / AVX2 kernels
 *  (where available) are bit-exact against the single-radius functions for
 *  representative inputs; associativity is preserved inside safe margins.
 *  @{ */
IRREP_API void irrep_rbf_bessel_batch    (int n, size_t N, const double *r,
                                          double r_cut, double *out);
IRREP_API void irrep_cutoff_cosine_batch (size_t N, const double *r,
                                          double r_cut, double *out);
IRREP_API void irrep_cutoff_polynomial_batch(size_t N, const double *r,
                                             double r_cut, int p, double *out);
IRREP_API void irrep_cutoff_cosine_d_batch   (size_t N, const double *r,
                                              double r_cut, double *out);
IRREP_API void irrep_cutoff_polynomial_d_batch(size_t N, const double *r,
                                               double r_cut, int p, double *out);
/** @} */

#ifdef __cplusplus
}
#endif

#endif /* IRREP_RADIAL_H */
