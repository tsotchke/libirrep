/* SPDX-License-Identifier: MIT */
#ifndef IRREP_CLEBSCH_GORDAN_H
#define IRREP_CLEBSCH_GORDAN_H

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Integer j: Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M>.
 * Returns 0.0 if the selection rules are violated (|j1 - j2| <= J <= j1 + j2,
 * m1 + m2 = M, |m_i| <= j_i). */
IRREP_API double irrep_cg(int j1, int m1, int j2, int m2, int J, int M);

/* Doubled-integer (half-spin) variant: every argument is 2 * the physical value. */
IRREP_API double irrep_cg_2j(int two_j1, int two_m1,
                             int two_j2, int two_m2,
                             int two_J,  int two_M);

/* Wigner 3j symbol (related to CG by a phase and a 1/sqrt(2J+1)). */
IRREP_API double irrep_wigner_3j   (int j1, int m1, int j2, int m2, int j3, int m3);
IRREP_API double irrep_wigner_3j_2j(int two_j1, int two_m1,
                                    int two_j2, int two_m2,
                                    int two_j3, int two_m3);

/* Cached CG table. Opaque to callers. */
typedef struct cg_table cg_table_t;

IRREP_API cg_table_t *irrep_cg_table_build   (int j1_max, int j2_max);
IRREP_API cg_table_t *irrep_cg_table_build_2j(int two_j1_max, int two_j2_max);
IRREP_API void        irrep_cg_table_free    (cg_table_t *table);

IRREP_API double      irrep_cg_lookup   (const cg_table_t *table,
                                         int j1, int m1, int j2, int m2, int J, int M);
IRREP_API double      irrep_cg_lookup_2j(const cg_table_t *table,
                                         int two_j1, int two_m1,
                                         int two_j2, int two_m2,
                                         int two_J,  int two_M);

/* Fill a flat (2j1+1)(2j2+1)(2J+1) block with CG coefficients; indexing is
 * [(m1 + j1) * (2j2+1) + (m2 + j2)] * (2J+1) + (M + J). */
IRREP_API void irrep_cg_block(int j1, int j2, int J, double *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_CLEBSCH_GORDAN_H */
