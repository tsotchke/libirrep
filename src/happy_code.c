/* SPDX-License-Identifier: MIT */
/** @file happy_code.c
 *  @brief HaPPY-code [[5, 1, 3]] perfect-tensor primitive. */
#include <irrep/happy_code.h>
#include <irrep/types.h>

#include <stddef.h>

irrep_status_t
irrep_happy_perfect_tensor_5_1_3(irrep_stabilizer_group_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_stabilizer_group_new(out, 5, 4);
    if (s != IRREP_OK) return s;

    /* Stabilizer table for the [[5, 1, 3]] code (Laflamme et al. 1996):
     *   g_1 = X Z Z X I
     *   g_2 = I X Z Z X
     *   g_3 = X I X Z Z
     *   g_4 = Z X I X Z
     * Cyclic shift of the same pattern (X Z Z X I) under rotation. */
    static const irrep_pauli_letter_t letters[4][5] = {
        { IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z,
          IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_I },
        { IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_Z,
          IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_X },
        { IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_X,
          IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z },
        { IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_I,
          IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_Z },
    };
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (letters[i][j] != IRREP_PAULI_LETTER_I) {
                irrep_pauli_set(&out->gens[i], j, letters[i][j]);
            }
        }
    }
    return IRREP_OK;
}
