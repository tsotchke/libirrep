/* SPDX-License-Identifier: MIT */
/** @file solid_harmonics.h
 *  @brief Real **regular solid harmonics** — homogeneous polynomials
 *         `R_{l,m}(r) = |r|^l · Y^real_{l,m}(r / |r|)`.
 *
 *  Inputs are un-normalised 3-vectors; outputs use the same
 *  `m = −l .. +l` layout as @ref irrep_sph_harm_cart. Normalisation matches
 *  the surface spherical harmonics, so `R_{l,m}(r̂) = Y^real_{l,m}(r̂)` for
 *  unit `|r̂| = 1`.
 *
 *  For `l ≤ 4` the library uses explicit polynomial expressions; higher `l`
 *  runs the Limpanuparb–Milthorpe recurrence with gradients co-evolved
 *  through each step — both paths are analytic and accurate to machine
 *  precision in double.
 */
#ifndef IRREP_SOLID_HARMONICS_H
#define IRREP_SOLID_HARMONICS_H

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Maximum order of solid harmonics supported (same as IRREP_L_MAX). */
#define IRREP_SOLID_L_MAX 16

/** @brief Write `2l + 1` solid harmonic values for a single @p l. */
IRREP_API void irrep_solid_harm_cart     (int l, double *out, const double r[3]);

/** @brief Write `3 · (2l + 1)` gradient components,
 *         `out[axis · (2l+1) + (m+l)] = ∂R_{l,m}/∂r_axis`. */
IRREP_API void irrep_solid_harm_cart_grad(int l, double *out, const double r[3]);

/** @brief Write `R_{l, m}(r)` for `l = 0 .. l_max`, packed as
 *         `[l=0 (1); l=1 (3); l=2 (5); …]`, total `(l_max + 1)²` values. */
IRREP_API void irrep_solid_harm_cart_all (int l_max, double *out, const double r[3]);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SOLID_HARMONICS_H */
