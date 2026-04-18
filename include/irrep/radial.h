/* SPDX-License-Identifier: MIT */
#ifndef IRREP_RADIAL_H
#define IRREP_RADIAL_H

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Bessel RBF of degree n on [0, r_cut]:
 *   phi_n(r) = sqrt(2 / r_cut) * sin(n * pi * r / r_cut) / r
 * Returns 0 for r > r_cut. */
IRREP_API double irrep_rbf_bessel    (int n, double r, double r_cut);
IRREP_API void   irrep_rbf_bessel_all(int n_max, double *out, double r, double r_cut);

/* Gaussian RBF: phi(r) = exp(-((r - mu)^2) / (2 * sigma^2)). */
IRREP_API double irrep_rbf_gaussian     (double r, double mu, double sigma);
IRREP_API void   irrep_rbf_gaussian_grid(int n, double *out, double r,
                                         double r_min, double r_max, double sigma);

/* Smooth cutoffs on [0, r_cut], zero for r > r_cut. */
IRREP_API double irrep_cutoff_cosine    (double r, double r_cut);
IRREP_API double irrep_cutoff_polynomial(double r, double r_cut, int p);

/* Cutoff derivatives (for force evaluation). */
IRREP_API double irrep_cutoff_cosine_d    (double r, double r_cut);
IRREP_API double irrep_cutoff_polynomial_d(double r, double r_cut, int p);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_RADIAL_H */
