/* SPDX-License-Identifier: MIT */
/** @file multiset.h
 *  @brief Construction, parsing, simplification, and direct sum of
 *         @ref irrep_multiset_t — e3nn-style direct-sum descriptions of
 *         feature spaces.
 *
 *  Grammar for irrep_multiset_parse():
 *
 *  @verbatim
 *    multiset := ε | term ( '+' term )*
 *    term     := mult 'x' l parity
 *    mult     := positive integer
 *    l        := non-negative integer (≤ IRREP_L_MAX)
 *    parity   := 'e' | 'o'
 *  @endverbatim
 *
 *  Whitespace is ignored. Empty input is accepted and yields an empty
 *  multiset (matching e3nn's convention for the trivial space).
 */
#ifndef IRREP_MULTISET_H
#define IRREP_MULTISET_H

#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Allocate an empty multiset with room for @p capacity terms. */
IRREP_API irrep_multiset_t *irrep_multiset_new(int capacity);

/** @brief Release a multiset built by any constructor in this header. */
IRREP_API void irrep_multiset_free(irrep_multiset_t *m);

/** @brief Append a `(label, multiplicity)` term. Grows the backing storage
 *         when necessary. Returns #IRREP_OK, or an error code on bad input. */
IRREP_API irrep_status_t irrep_multiset_append(irrep_multiset_t *m, irrep_label_t label,
                                               int multiplicity);

/** @brief Parse an e3nn-style spec (e.g. `"1x0e + 2x1o + 1x2e"`).
 *  @return freshly-allocated multiset, or @c NULL on malformed input; the
 *          reason is available via @ref irrep_last_error(). */
IRREP_API irrep_multiset_t *irrep_multiset_parse(const char *spec);

/** @brief Format into @p buf. Writes at most `buflen - 1` chars plus a NUL.
 *  @return the required length (excluding terminator); may exceed @p buflen. */
IRREP_API int irrep_multiset_format(const irrep_multiset_t *m, char *buf, size_t buflen);

/** @brief Merge like terms, sort canonically by `(l ascending, even-before-odd)`.
 *         Idempotent. */
IRREP_API void irrep_multiset_simplify(irrep_multiset_t *m);

/** @brief Direct sum `m1 ⊕ m2`. Returns a newly-allocated multiset. */
IRREP_API irrep_multiset_t *irrep_multiset_direct_sum(const irrep_multiset_t *m1,
                                                      const irrep_multiset_t *m2);

/** @brief Total dimension (synonym for `m->total_dim`). */
IRREP_API int irrep_multiset_dim(const irrep_multiset_t *m);

/** @brief Start offset of term @p term_idx in a feature vector laid out by
 *         this multiset. Returns -1 on out-of-range input. */
IRREP_API int irrep_multiset_block_offset(const irrep_multiset_t *m, int term_idx);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_MULTISET_H */
