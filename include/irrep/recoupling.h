/* SPDX-License-Identifier: MIT */
/** @file recoupling.h
 *  @brief Wigner 6j, 9j, and Racah-W recoupling coefficients.
 *
 *  The 6j symbol measures the two equivalent ways of coupling three angular
 *  momenta; the 9j the four ways of coupling four. All return `0` on
 *  selection-rule violation (triangle conditions failing). Doubled-integer
 *  (`_2j`) variants allow half-integer inputs.
 */
#ifndef IRREP_RECOUPLING_H
#define IRREP_RECOUPLING_H

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Wigner 6j symbol `{j1 j2 j3; j4 j5 j6}`. */
IRREP_API double irrep_wigner_6j   (int j1, int j2, int j3,
                                    int j4, int j5, int j6);
/** @brief Doubled-integer variant (half-spin). */
IRREP_API double irrep_wigner_6j_2j(int two_j1, int two_j2, int two_j3,
                                    int two_j4, int two_j5, int two_j6);

/** @brief Wigner 9j symbol `{j1 j2 j3; j4 j5 j6; j7 j8 j9}`. */
IRREP_API double irrep_wigner_9j   (int j1, int j2, int j3,
                                    int j4, int j5, int j6,
                                    int j7, int j8, int j9);
/** @brief Doubled-integer variant. */
IRREP_API double irrep_wigner_9j_2j(int two_j1, int two_j2, int two_j3,
                                    int two_j4, int two_j5, int two_j6,
                                    int two_j7, int two_j8, int two_j9);

/** @brief Racah W coefficient — related to the 6j symbol by a phase. */
IRREP_API double irrep_racah_w(int j1, int j2, int J, int j3, int j12, int j23);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_RECOUPLING_H */
