/* SPDX-License-Identifier: MIT */
#ifndef IRREP_TIME_REVERSAL_H
#define IRREP_TIME_REVERSAL_H

#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Time-reversal operator on an integer-l irrep: written as a (2l+1) x (2l+1)
 * complex matrix (complex conjugation in the m basis is a diagonal phase; the
 * explicit matrix is the identity times (-1)^l in Condon-Shortley). */
IRREP_API void irrep_time_reversal_integer(int l, double _Complex *out);

/* Half-integer: T = i * sigma_y on the spin-1/2 factor, composed per copy.
 * `two_j` is the doubled spin label; output is (2j+1) x (2j+1) complex. */
IRREP_API void irrep_time_reversal_half_integer(int two_j, double _Complex *out);

/* Block-diagonal T on a multiset; writes total_dim x total_dim complex. */
IRREP_API void irrep_time_reversal_multiset(const irrep_multiset_t *m, double _Complex *out);

/* T^2 = +1 (all integer-l components) or -1 (any half-integer component).
 * Returns +1 or -1 (or 0 if the multiset is empty). */
IRREP_API int irrep_time_reversal_square_sign(const irrep_multiset_t *m);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_TIME_REVERSAL_H */
