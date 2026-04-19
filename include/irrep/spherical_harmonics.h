/* SPDX-License-Identifier: MIT */
/** @file spherical_harmonics.h
 *  @brief Spherical harmonics on S² in polar, real-polar, and cartesian form,
 *         plus gradients and the complex↔real basis change.
 *
 *  Conventions (single source of truth):
 *  - **Condon-Shortley** phase `(−1)^m` applied once in the associated
 *    Legendre polynomial, not repeated downstream.
 *  - **Orthonormal** on the 2-sphere: `∫ Yₗᵐ · Yₗ′ᵐ′* dΩ = δ_{ll′} δ_{mm′}`.
 *  - **e3nn real-SH sign**: `Y_{l,+1} ∝ +x` at the equator (matches Wikipedia
 *    and e3nn, differs from some quantum-chemistry packages).
 *  - Storage layout for a single `l`: `out[m + l]`, so `m = −l` is index 0
 *    and `m = +l` is index `2l`.
 *  - Cartesian entry points take a unit vector `r̂`; no (θ, φ) conversion.
 */
#ifndef IRREP_SPHERICAL_HARMONICS_H
#define IRREP_SPHERICAL_HARMONICS_H

#include <complex.h>
#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Complex spherical harmonic `Yₗᵐ(θ, φ)` (Condon-Shortley, orthonormal). */
IRREP_API double _Complex irrep_sph_harm     (int l, int m, double theta, double phi);

/** @brief Real spherical harmonic `Y_{l,m}^{real}(θ, φ)` (e3nn sign convention). */
IRREP_API double          irrep_sph_harm_real(int l, int m, double theta, double phi);

/** @brief Cartesian real SH: writes `2l + 1` values into @p out at indices
 *         `m + l`. @p r_hat must be unit-norm. */
IRREP_API void irrep_sph_harm_cart(int l, double *out, const double r_hat[3]);

/** @brief All cartesian real SH up to @p l_max, packed as `[l=0; l=1(−1,0,+1);
 *         l=2(−2..+2); …]`. Total output size: `(l_max + 1)²`. */
IRREP_API void irrep_sph_harm_cart_all(int l_max, double *out, const double r_hat[3]);

/** @brief Gradient of cartesian real SH: writes `3 · (2l + 1)` reals
 *         laid out as `out[axis · (2l + 1) + (m + l)]`, axis ∈ {0:x, 1:y, 2:z}. */
IRREP_API void irrep_sph_harm_cart_grad(int l, double *out, const double r_hat[3]);

/** @brief Associated Legendre polynomial `Pₗᵐ(x)` at `x = cos θ`, Condon-Shortley
 *         phase, valid for `|m| ≤ l`. */
IRREP_API double irrep_legendre_assoc(int l, int m, double x);

/** @brief Change-of-basis matrix `U` such that `Y_real = U · Y_complex`.
 *         @p out is row-major `(2l+1) × (2l+1)` complex; rows and columns are
 *         indexed `m + l`. */
IRREP_API void irrep_sph_harm_complex_to_real(int l, double _Complex *out);

/** @brief Single-precision wrapper — delegates to the double kernel and casts. */
IRREP_API void irrep_sph_harm_cart_f32    (int l,     float *out, const float r_hat[3]);
/** @brief Single-precision wrapper for #irrep_sph_harm_cart_all. */
IRREP_API void irrep_sph_harm_cart_all_f32(int l_max, float *out, const float r_hat[3]);

/** @brief Batched cartesian SH: fills `N · (l_max+1)²` values, one contiguous
 *         block per unit-vector input @p r_hats. Routes through the runtime
 *         SIMD dispatch table. */
IRREP_API void irrep_sph_harm_cart_all_batch(int l_max, size_t N,
                                             const double *r_hats,
                                             double *out);

/** @brief Batched gradient of cartesian real SH up to @p l_max.
 *
 *  For each of @p N unit-vector inputs, writes `3 · (l_max+1)²` doubles
 *  packed as `out[e, axis, l, m] → out[(e · 3 + axis) · (l_max+1)² + (l² + m + l)]`
 *  — axis ∈ {0:x, 1:y, 2:z}. This is the force-evaluation hot path for
 *  NequIP / MACE / Allegro. Routes through the runtime SIMD dispatch. */
IRREP_API void irrep_sph_harm_cart_all_grad_batch(int l_max, size_t N,
                                                  const double *r_hats,
                                                  double *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SPHERICAL_HARMONICS_H */
