/* SPDX-License-Identifier: MIT */
/** @file quadrature.h
 *  @brief Spherical (Lebedev, tensor-product) and 1-D (Gauss-Legendre)
 *         quadrature rules, used internally for the SH orthonormality tests
 *         and exposed for general solid-angle integration.
 */
#ifndef IRREP_QUADRATURE_H
#define IRREP_QUADRATURE_H

#include <stdbool.h>

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Lebedev rule size (points) for @p order ∈
 *         {3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41}.
 *  Returns 0 for unknown orders. */
IRREP_API int irrep_lebedev_size(int order);

/** @brief Fill a `size × 4` array `[x, y, z, weight]` for the Lebedev rule of
 *         @p order. Weights sum to 1. Returns @c false if @p order is unknown. */
IRREP_API bool irrep_lebedev_fill(int order, double *xyz_weights);

/** @brief Gauss-Legendre nodes and weights on `[−1, 1]`. Exact for polynomials
 *         of degree `≤ 2n − 1`. */
IRREP_API bool irrep_gauss_legendre(int n, double *nodes, double *weights);

/** @brief Tensor-product S² rule: Gauss-Legendre in `cos θ` × uniform in `φ`.
 *         Exact for real SH up to @p exactness_deg. Covers any order; use
 *         this when the Lebedev tabulated order is missing.
 *  @return number of sample points `N = (exactness_deg/2 + 1) · (exactness_deg + 1)`. */
IRREP_API int irrep_quadrature_sphere_size(int exactness_deg);

/** @brief Fill a `N × 4` `[x, y, z, weight]` array for the tensor-product rule.
 *         Weights sum to 1. */
IRREP_API bool irrep_quadrature_sphere_fill(int exactness_deg, double *xyz_weights);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_QUADRATURE_H */
