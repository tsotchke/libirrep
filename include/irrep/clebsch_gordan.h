/* SPDX-License-Identifier: MIT */
/** @file clebsch_gordan.h
 *  @brief Clebsch-Gordan coefficients and Wigner-3j symbols for both integer
 *         and half-integer angular momentum, plus a cached table builder.
 *
 *  Coefficients satisfy the Condon-Shortley phase convention and the standard
 *  orthonormality relations. The kernel is Miller two-directional iteration
 *  over the Schulten–Gordon three-term recurrence for Wigner-3j (Luscombe–
 *  Luban, *Phys. Rev. E* **57**, 7274, 1998); Clebsch–Gordan is derived from
 *  3j via `⟨j₁ m₁ j₂ m₂ | J M⟩ = (−1)^{j₁−j₂+M} √(2J+1) (j₁ j₂ J; m₁ m₂ −M)`.
 *  Machine precision through at least `j = 200` in regression tests.
 *
 *  Every coefficient-returning function returns exactly `0.0` when the
 *  selection rules
 *  `|j1 − j2| ≤ J ≤ j1 + j2`, `m1 + m2 = M`, `|m_i| ≤ j_i`
 *  are violated — so a zero result is indistinguishable from an invalid one.
 *  (Builder functions signal failure via @c NULL + @ref irrep_last_error.)
 */
#ifndef IRREP_CLEBSCH_GORDAN_H
#define IRREP_CLEBSCH_GORDAN_H

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Integer-j Clebsch-Gordan coefficient `⟨j1 m1; j2 m2 | J M⟩`. */
IRREP_API double irrep_cg(int j1, int m1, int j2, int m2, int J, int M);

/** @brief Doubled-integer (half-spin) variant. Every argument is `2 ·` its
 *         physical value; e.g. `j = 1/2` → `two_j = 1`. */
IRREP_API double irrep_cg_2j(int two_j1, int two_m1,
                             int two_j2, int two_m2,
                             int two_J,  int two_M);

/** @brief Wigner 3j symbol — same information as CG with a cleaner column
 *         symmetry, related by
 *         `⟨j1 m1; j2 m2 | j3 m3⟩ = (−1)^{j1−j2+m3} √(2j3+1) · (j1 j2 j3 / m1 m2 −m3)`. */
IRREP_API double irrep_wigner_3j   (int j1, int m1, int j2, int m2, int j3, int m3);
/** @brief Doubled-integer Wigner 3j symbol. */
IRREP_API double irrep_wigner_3j_2j(int two_j1, int two_m1,
                                    int two_j2, int two_m2,
                                    int two_j3, int two_m3);

/** @brief Opaque pre-built CG lookup table. */
typedef struct cg_table cg_table_t;

/** @brief Build a CG lookup table covering integer `j1 ≤ j1_max`, `j2 ≤ j2_max`. */
IRREP_API cg_table_t *irrep_cg_table_build   (int j1_max, int j2_max);
/** @brief Doubled-integer variant — covers half-spins too. */
IRREP_API cg_table_t *irrep_cg_table_build_2j(int two_j1_max, int two_j2_max);
/** @brief Release a table built by either #irrep_cg_table_build variant. */
IRREP_API void        irrep_cg_table_free    (cg_table_t *table);

/** @brief O(1) lookup; returns `0.0` for out-of-range queries (same contract
 *         as the direct computation). */
IRREP_API double      irrep_cg_lookup   (const cg_table_t *table,
                                         int j1, int m1, int j2, int m2, int J, int M);
/** @brief Doubled-integer lookup. */
IRREP_API double      irrep_cg_lookup_2j(const cg_table_t *table,
                                         int two_j1, int two_m1,
                                         int two_j2, int two_m2,
                                         int two_J,  int two_M);

/** @brief Fill a flat `(2j1+1) · (2j2+1) · (2J+1)` block.
 *  Index layout: `[(m1 + j1) · (2j2+1) + (m2 + j2)] · (2J+1) + (M + J)`. */
IRREP_API void irrep_cg_block(int j1, int j2, int J, double *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_CLEBSCH_GORDAN_H */
