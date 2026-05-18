/* SPDX-License-Identifier: MIT */
/** @file color_code_3d.c
 *  @brief Implementation of the `[[8, 3, 2]]` cubic 3D color code. */
#include <irrep/color_code_3d.h>
#include <irrep/types.h>

#include <stddef.h>

/* Face / volume supports.  Cube qubits indexed 0..7 with binary
 * coordinates (x, y, z) and linearisation `4z + 2y + x`. */
static const int kTopFace[4]   = { 4, 5, 6, 7 };  /* z = 1 */
static const int kLeftFace[4]  = { 0, 2, 4, 6 };  /* x = 0 */
static const int kFrontFace[4] = { 0, 1, 4, 5 };  /* y = 0 */

irrep_status_t
irrep_color_3d_cube_8_3_2(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;

    irrep_status_t s = irrep_css_code_new(out, IRREP_COLOR_3D_CUBE_N,
                                          IRREP_COLOR_3D_CUBE_N_X_STABS,
                                          IRREP_COLOR_3D_CUBE_N_Z_STABS);
    if (s != IRREP_OK) return s;

    /* H_X rows: 3 face stabs + 1 volume stab. */
    for (int i = 0; i < 4; ++i) irrep_parity_matrix_set(&out->H_X, 0, kTopFace[i]);
    for (int i = 0; i < 4; ++i) irrep_parity_matrix_set(&out->H_X, 1, kLeftFace[i]);
    for (int i = 0; i < 4; ++i) irrep_parity_matrix_set(&out->H_X, 2, kFrontFace[i]);
    for (int q = 0; q < IRREP_COLOR_3D_CUBE_N; ++q) {
        irrep_parity_matrix_set(&out->H_X, 3, q);  /* X-volume */
    }

    /* H_Z rows: 1 volume stab. */
    for (int q = 0; q < IRREP_COLOR_3D_CUBE_N; ++q) {
        irrep_parity_matrix_set(&out->H_Z, 0, q);  /* Z-volume */
    }

    return IRREP_OK;
}

/* ====================================================================
 * [[15, 1, 3]] Reed-Muller / tetrahedral 3D color code
 *
 * Qubit indexing: qubit q ∈ {0..14} corresponds to PG(3,2) point
 * (q + 1), i.e. the binary 4-tuple `q + 1`. The bit `b ∈ {0..3}` of
 * `(q + 1)` selects an X-stabilizer's "single-bit" support; pairs of
 * bits select a Z-stabilizer's weight-4 "pair" support.
 *
 * The construction uses punctured Reed-Muller codes:
 *   H_X = punctured RM(1, 4) generators (one per bit, weight 8)
 *   H_Z = even-weight subcode of dual punctured RM(2, 4)
 *        = [the 4 weight-8 single-bit rows] + [6 weight-4 pair rows]
 *
 * Commutation: a single-bit and a pair row share exactly 2 or 4 qubits
 * (even), so H_X · H_Z^T = 0 over GF(2). ✓
 * ==================================================================== */

irrep_status_t
irrep_color_3d_rm_15_1_3(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;

    irrep_status_t s = irrep_css_code_new(out, IRREP_COLOR_3D_RM15_N,
                                          IRREP_COLOR_3D_RM15_N_X_STABS,
                                          IRREP_COLOR_3D_RM15_N_Z_STABS);
    if (s != IRREP_OK) return s;

    /* H_X: 4 weight-8 rows. Row b's support = {q : bit b of (q + 1) = 1}. */
    for (int b = 0; b < 4; ++b) {
        for (int q = 0; q < IRREP_COLOR_3D_RM15_N; ++q) {
            int point = q + 1;
            if ((point >> b) & 1) {
                irrep_parity_matrix_set(&out->H_X, b, q);
            }
        }
    }

    /* H_Z rows 0..3: same 4 weight-8 single-bit supports as H_X (since
     * C_X ⊂ C_Z in the punctured-RM CSS construction). */
    for (int b = 0; b < 4; ++b) {
        for (int q = 0; q < IRREP_COLOR_3D_RM15_N; ++q) {
            int point = q + 1;
            if ((point >> b) & 1) {
                irrep_parity_matrix_set(&out->H_Z, b, q);
            }
        }
    }

    /* H_Z rows 4..9: 6 weight-4 pair generators, lexicographic on the
     * pair (b_1, b_2). Row's support = {q : bit b_1 AND bit b_2 of
     * (q + 1) are both set}. */
    int row = 4;
    for (int b1 = 0; b1 < 4; ++b1) {
        for (int b2 = b1 + 1; b2 < 4; ++b2) {
            for (int q = 0; q < IRREP_COLOR_3D_RM15_N; ++q) {
                int point = q + 1;
                if (((point >> b1) & 1) && ((point >> b2) & 1)) {
                    irrep_parity_matrix_set(&out->H_Z, row, q);
                }
            }
            ++row;
        }
    }
    return IRREP_OK;
}
