/* SPDX-License-Identifier: MIT */
/** @file color_code.c
 *  @brief Implementation of triangular color codes.
 *
 *  Currently only the [[7, 1, 3]] Steane variant is exposed; this is
 *  the smallest color code, with three weight-4 X-stabs and three
 *  weight-4 Z-stabs sharing the same support (one per face/colour). */
#include <irrep/color_code.h>
#include <irrep/types.h>

#include <stddef.h>

irrep_status_t
irrep_color_steane(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_css_code_new(out, 7, 3, 3);
    if (s != IRREP_OK) return s;

    /* Hamming-code parity matrix (cols 1..7 in 1-indexed convention,
     * 0..6 in 0-indexed):
     *   row 0 (face R): cols 1, 3, 5, 7  →  0-idx 0, 2, 4, 6
     *   row 1 (face G): cols 2, 3, 6, 7  →  0-idx 1, 2, 5, 6
     *   row 2 (face B): cols 4, 5, 6, 7  →  0-idx 3, 4, 5, 6
     */
    static const int support[3][4] = {
        { 0, 2, 4, 6 },
        { 1, 2, 5, 6 },
        { 3, 4, 5, 6 },
    };
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 4; ++k) {
            irrep_parity_matrix_set(&out->H_X, i, support[i][k]);
            irrep_parity_matrix_set(&out->H_Z, i, support[i][k]);
        }
    }
    return IRREP_OK;
}
