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
