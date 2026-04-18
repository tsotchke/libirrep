/* SPDX-License-Identifier: MIT */
#ifndef IRREP_SPHERICAL_HARMONICS_H
#define IRREP_SPHERICAL_HARMONICS_H

#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Polar-coordinate forms (Condon-Shortley phase, orthonormal on S^2). */
IRREP_API double _Complex irrep_sph_harm     (int l, int m, double theta, double phi);
IRREP_API double          irrep_sph_harm_real(int l, int m, double theta, double phi);

/* Cartesian: writes (2l+1) real values for a single l; `r_hat` must be unit. */
IRREP_API void irrep_sph_harm_cart(int l, double *out, const double r_hat[3]);

/* Batched cartesian: writes sum_{l=0..l_max} (2l+1) = (l_max+1)^2 real values,
 * laid out as [l=0; l=1 (-1,0,+1); l=2 (-2..+2); ...]. */
IRREP_API void irrep_sph_harm_cart_all(int l_max, double *out, const double r_hat[3]);

/* Gradient of cartesian real SH: writes 3 * (2l+1) reals (dx, dy, dz planes). */
IRREP_API void irrep_sph_harm_cart_grad(int l, double *out, const double r_hat[3]);

/* Associated Legendre polynomial P_l^m(x), x = cos(theta), Condon-Shortley phase. */
IRREP_API double irrep_legendre_assoc(int l, int m, double x);

/* Change-of-basis matrix U s.t. Y_real = U * Y_complex, (2l+1) x (2l+1), complex. */
IRREP_API void irrep_sph_harm_complex_to_real(int l, double _Complex *out);

/* Single-precision variants (SIMD hot paths). */
IRREP_API void irrep_sph_harm_cart_f32    (int l,     float *out, const float r_hat[3]);
IRREP_API void irrep_sph_harm_cart_all_f32(int l_max, float *out, const float r_hat[3]);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SPHERICAL_HARMONICS_H */
