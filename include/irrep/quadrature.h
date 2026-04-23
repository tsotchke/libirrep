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
 *  Returns 0 for unknown orders. Orders 3 / 5 / 7 are hard-coded; all
 *  higher orders become known at runtime after a call to
 *  @ref irrep_lebedev_register_rule (see `scripts/fetch_lebedev_tables.sh`
 *  for the Lebedev-Laikov 1999 tables in public-domain form, and
 *  `examples/register_lebedev.c` for the loader). */
IRREP_API int irrep_lebedev_size(int order);

/** @brief Fill a `size × 4` array `[x, y, z, weight]` for the Lebedev rule of
 *         @p order. Weights sum to 1. Returns @c false if @p order is unknown. */
IRREP_API bool irrep_lebedev_fill(int order, double *xyz_weights);

/** @brief Register a Lebedev rule at runtime. Copies the caller's
 *         @p xyz_weights buffer into library-owned storage; subsequent
 *         calls to @ref irrep_lebedev_size and @ref irrep_lebedev_fill
 *         at this @p order return the registered rule.
 *
 *  Replaces any previously-registered rule at the same order.
 *  Hard-coded orders (3 / 5 / 7) cannot be overridden and return
 *  `IRREP_ERR_PRECONDITION`. The library keeps the registration until
 *  @ref irrep_lebedev_clear_registry or process exit.
 *
 *  This is the path for the Lebedev-Laikov 1999 tables at orders
 *  9–41: the shipped `scripts/fetch_lebedev_tables.sh` downloads the
 *  public-domain data, `examples/register_lebedev.c` parses and
 *  registers; no data file is bundled in the library tree.
 *
 *  @param order          Lebedev rule order (must be odd, ≥ 3).
 *  @param n_points       number of quadrature points in @p xyz_weights.
 *  @param xyz_weights    `n_points × 4` row-major `(x, y, z, w)` array;
 *                        each `(x, y, z)` must be unit-sphere-normalised
 *                        and weights must sum to 1 to within 1e-10.
 *  @return `IRREP_OK`, `IRREP_ERR_INVALID_ARG` (bad args, non-unit
 *          points, weights don't sum to 1), `IRREP_ERR_PRECONDITION`
 *          (order 3/5/7 — hard-coded, cannot override), or
 *          `IRREP_ERR_OUT_OF_MEMORY`. */
#include <irrep/types.h>
IRREP_API irrep_status_t irrep_lebedev_register_rule(int order, int n_points,
                                                     const double *xyz_weights);

/** @brief Discard all runtime-registered Lebedev rules. Hard-coded orders
 *         remain available. Useful for tests or repeated registration
 *         cycles; unnecessary in normal use. */
IRREP_API void irrep_lebedev_clear_registry(void);

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
