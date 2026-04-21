/* SPDX-License-Identifier: MIT */
/** @file multiset_2j.h
 *  @brief Doubled-integer multiset for mixed integer + half-integer spin
 *         content. Complements @ref irrep_multiset_t, which is restricted
 *         to integer `l`.
 *
 *  Grammar (extends the e3nn-style multiset parser):
 *
 *  @verbatim
 *    multiset := ε | term ( '+' term )*
 *    term     := mult 'x' label parity
 *    label    := integer | half-integer
 *    integer  := '0' | '1' | '2' | ...                (e.g. "1x0e", "2x2e")
 *    half-integer := N '/2'  (N positive odd integer)  (e.g. "1x1/2o", "2x3/2e")
 *    parity   := 'e' | 'o'
 *  @endverbatim
 *
 *  Every entry stores `two_j = 2·j`, so spin-½ has `two_j = 1`, spin-1 has
 *  `two_j = 2`, etc. Block dimension is `2j + 1 = two_j + 1`.
 *
 *  Added in 1.3 (advance preview shipped in 1.2.1): lets consumers that
 *  carry half-integer spin features (magnetic moments in spinor
 *  representations, spin-orbit-coupled orbital bases, etc.) obtain correct
 *  Kramers degeneracy information from
 *  @ref irrep_time_reversal_square_sign_2j. */
#ifndef IRREP_MULTISET_2J_H
#define IRREP_MULTISET_2J_H

#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Doubled-spin irrep label: `two_j ∈ {0, 1, 2, …}`. */
typedef struct {
    int two_j;  /**< doubled-integer spin label */
    int parity; /**< +1 (even) or -1 (odd) */
} irrep_label_2j_t;

/** @brief Direct sum of doubled-integer irreps. `total_dim = Σ mult_i · (two_j_i + 1)`. */
typedef struct {
    irrep_label_2j_t *labels;
    int              *multiplicities;
    int               num_terms;
    int               capacity;
    int               total_dim;
} irrep_multiset_2j_t;

/** @brief Allocate an empty multiset with room for @p capacity terms. */
IRREP_API irrep_multiset_2j_t *irrep_multiset_2j_new(int capacity);

/** @brief Release a multiset. */
IRREP_API void irrep_multiset_2j_free(irrep_multiset_2j_t *m);

/** @brief Append `(label, mult)`. Grows backing storage if necessary. */
IRREP_API irrep_status_t irrep_multiset_2j_append(irrep_multiset_2j_t *m, irrep_label_2j_t label,
                                                  int multiplicity);

/** @brief Parse a spec string supporting both integer (e.g. `"1x0e"`)
 *         and half-integer (e.g. `"2x1/2o"`) terms. */
IRREP_API irrep_multiset_2j_t *irrep_multiset_2j_parse(const char *spec);

/** @brief Format into @p buf. Writes at most `buflen - 1` chars + NUL.
 *  @return required length (excluding terminator). */
IRREP_API int irrep_multiset_2j_format(const irrep_multiset_2j_t *m, char *buf, size_t buflen);

/** @brief Total dimension (synonym for `m->total_dim`). */
IRREP_API int irrep_multiset_2j_dim(const irrep_multiset_2j_t *m);

/** @brief Check whether any block has odd `two_j` (half-integer spin).
 *         Returns 1 if yes, 0 if purely integer-spin. */
IRREP_API int irrep_multiset_2j_has_half_integer(const irrep_multiset_2j_t *m);

/** @brief Kramers-sign check: returns `−1` if `two_j_total` (the sum of
 *         `(two_j + 1) · mult` reduced mod 4) falls in the half-integer
 *         regime, `+1` if the multiset is purely integer-spin, `0` on
 *         empty input.
 *
 *  Concretely: `T² = (−1)^{2j}` on a spin-j block (Sakurai §4.4). For a
 *  direct sum, `T²` has a definite sign iff every block has the same
 *  parity of `two_j`; if any half-integer block is present and any
 *  integer block is also present, the return is still `−1` by Kramers'
 *  argument — the sign is carried by the half-integer sector, the
 *  integer one just inherits `+1` and the overall `T²` is the direct
 *  sum. Callers wanting the per-block sign should iterate the multiset
 *  directly. */
IRREP_API int irrep_time_reversal_square_sign_2j(const irrep_multiset_2j_t *m);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_MULTISET_2J_H */
