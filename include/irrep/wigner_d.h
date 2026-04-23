/* SPDX-License-Identifier: MIT */
/** @file wigner_d.h
 *  @brief Wigner small-d and full-D rotation matrices on SO(3) / SU(2)
 *         irreps, plus the block-diagonal lifting to an @ref irrep_multiset_t.
 *
 *  Full D factors as
 *
 *  @verbatim
 *      D^j_{m'm}(α, β, γ) = e^{−i m' α} · d^j_{m'm}(β) · e^{−i m γ}
 *  @endverbatim
 *
 *  Small-d is evaluated via the Edmonds (4.1.23) Jacobi-polynomial form
 *  with symmetry reduction to the canonical `m ≥ |m'|` region
 *  (Varshalovich §4.4.1); the Jacobi polynomial itself uses the NIST
 *  DLMF §18.9.1 forward three-term recurrence. Measured unitarity is
 *  `≤ 10⁻¹²` through at least `j = 80` and bounded only by the
 *  IEEE-754 `lgamma` overflow limit past that. See `src/wigner_d.c`
 *  for the implementation and `docs/METHODS.md` §3.2 for references.
 */
#ifndef IRREP_WIGNER_D_H
#define IRREP_WIGNER_D_H

#include <complex.h>
#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Small Wigner-d scalar `d^j_{mp, m}(β)` — real, depends on β only. */
IRREP_API double irrep_wigner_d_small(int j, int mp, int m, double beta);
/** @brief Doubled-integer variant (half-spin). */
IRREP_API double irrep_wigner_d_small_2j(int two_j, int two_mp, int two_m, double beta);

/** @brief Full Wigner-D scalar `D^j_{mp, m}(α, β, γ)`. */
IRREP_API double _Complex irrep_wigner_D(int j, int mp, int m, double alpha, double beta,
                                         double gamma);
/** @brief Doubled-integer variant (half-spin). */
IRREP_API double _Complex irrep_wigner_D_2j(int two_j, int two_mp, int two_m, double alpha,
                                            double beta, double gamma);

/** @brief Full `(2j+1) × (2j+1)` complex Wigner-D matrix, row-major,
 *         `m'` indexes rows, `m` indexes columns. */
IRREP_API void irrep_wigner_D_matrix(int j, double _Complex *out, double alpha, double beta,
                                     double gamma);

/** @brief Full real small-d matrix (same layout). */
IRREP_API void irrep_wigner_d_matrix(int j, double *out, double beta);

/** @brief Batched small-d matrix build for many β values.
 *
 *  Output layout: `out[b * (2j+1)² + mp * (2j+1) + m] = d^j_{mp, m}(β_b)`.
 *  Callers invoking @ref irrep_wigner_d_matrix inside a hot loop over
 *  edge angles (e.g. one rotation per graph edge in a NequIP forward
 *  pass) should use this batched variant: it runs pairs of β values
 *  through the Jacobi recurrence in lockstep under a runtime-dispatched
 *  NEON kernel on arm64, giving ~2× throughput on the dominant inner
 *  loop. The scalar fallback is a trivial per-β loop.
 *
 *  @param j        spin quantum (≥ 0).
 *  @param n_betas  number of β values.
 *  @param betas    length-`n_betas` array of β values.
 *  @param out      output buffer of size `n_betas · (2j+1)²`. */
IRREP_API void irrep_wigner_d_matrix_batch(int j, size_t n_betas, const double *betas, double *out);

/** @brief Block-diagonal Wigner-D lifted onto an @ref irrep_multiset_t.
 *         Writes `total_dim × total_dim` complex entries; each irrep term
 *         with multiplicity μ contributes μ identical `(2l+1) × (2l+1)` blocks. */
IRREP_API void irrep_wigner_D_multiset(const irrep_multiset_t *m, double _Complex *out,
                                       double alpha, double beta, double gamma);

/** @brief Partial derivative `∂ d^j_{mp, m} / ∂β` — useful for forces and
 *         for finite-difference sanity checks on rotational equivariance. */
IRREP_API double irrep_wigner_d_small_dbeta(int j, int mp, int m, double beta);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_WIGNER_D_H */
