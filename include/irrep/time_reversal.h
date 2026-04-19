/* SPDX-License-Identifier: MIT */
/** @file time_reversal.h
 *  @brief Time-reversal operator `T` on SO(3) / SU(2) irreps and on a whole
 *         multiset feature vector.
 *
 *  `T` is anti-linear: `T = K` (complex conjugation) on integer-l irreps,
 *  `T = i σ_y · K` on half-integer spin-½ factors. The composite operator
 *  squares to `+1` on even-spin combinations and to `−1` on odd-spin ones
 *  (Kramers degeneracy).
 */
#ifndef IRREP_TIME_REVERSAL_H
#define IRREP_TIME_REVERSAL_H

#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Matrix of `T` on an integer-l irrep in the `m` basis. Writes
 *         `(2l+1) × (2l+1)` complex entries. */
IRREP_API void irrep_time_reversal_integer(int l, double _Complex *out);

/** @brief Matrix of `T` on a half-integer irrep of doubled-spin @p two_j.
 *         Writes `(2j+1) × (2j+1)` complex entries. */
IRREP_API void irrep_time_reversal_half_integer(int two_j, double _Complex *out);

/** @brief Block-diagonal `T` on an @ref irrep_multiset_t; writes
 *         `total_dim × total_dim` complex entries. */
IRREP_API void irrep_time_reversal_multiset(const irrep_multiset_t *m, double _Complex *out);

/** @brief Returns `+1` on any non-empty multiset, `0` on empty input.
 *
 *  Rationale: @ref irrep_multiset_t carries only integer-l labels (the
 *  e3nn-style "Nxle | Nxlo" grammar does not encode half-integer spin), so
 *  a multiset-level query cannot detect Kramers degeneracy. This function
 *  therefore returns `+1` unconditionally on integer-only content; for
 *  half-integer content use `irrep_time_reversal_square_sign_2j`
 *  (header @c irrep/multiset_2j.h ), which takes a doubled-integer
 *  multiset and returns `−1` when any block has odd `two_j`. Callers
 *  that need the per-block matrix itself go through
 *  @ref irrep_time_reversal_half_integer. */
IRREP_API int irrep_time_reversal_square_sign(const irrep_multiset_t *m);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_TIME_REVERSAL_H */
