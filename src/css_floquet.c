/* SPDX-License-Identifier: MIT */
/** @file css_floquet.c
 *  @brief CSS Floquet code construction on the square lattice. */
#include <irrep/css_floquet.h>
#include <irrep/types.h>

#include <stddef.h>

static inline int wrap(int v, int L) {
    int r = v % L;
    if (r < 0) r += L;
    return r;
}

/* Qubit at (x, y, sub): sub ∈ {0=A, 1=B}. Linear index. */
static inline int q_idx(int x, int y, int sub, int Lx)
{
    return 2 * (y * Lx + x) + sub;
}

irrep_status_t
irrep_css_floquet_square_build(int Lx, int Ly, irrep_floquet_code_t *out)
{
    if (out == NULL || Lx <= 0 || Ly <= 0) return IRREP_ERR_INVALID_ARG;
    int n_qubits = 2 * Lx * Ly;
    int n_meas_per_round = Lx * Ly;

    irrep_status_t s = irrep_floquet_code_new(out, n_qubits, 4);
    if (s != IRREP_OK) return s;
    for (int r = 0; r < 4; ++r) {
        s = irrep_floquet_round_alloc(out, r, n_meas_per_round);
        if (s != IRREP_OK) {
            irrep_floquet_code_free(out);
            return s;
        }
    }

    /* Round 0: XX on horizontal A-B edges within each unit cell.
     * Round 1: ZZ on the same horizontal edges. */
    int idx = 0;
    for (int y = 0; y < Ly; ++y) {
        for (int x = 0; x < Lx; ++x) {
            int qa = q_idx(x, y, 0, Lx);
            int qb = q_idx(x, y, 1, Lx);
            irrep_pauli_set(&out->rounds[0].meas[idx], qa, IRREP_PAULI_LETTER_X);
            irrep_pauli_set(&out->rounds[0].meas[idx], qb, IRREP_PAULI_LETTER_X);
            irrep_pauli_set(&out->rounds[1].meas[idx], qa, IRREP_PAULI_LETTER_Z);
            irrep_pauli_set(&out->rounds[1].meas[idx], qb, IRREP_PAULI_LETTER_Z);
            ++idx;
        }
    }

    /* Round 2: XX on vertical edges connecting B(x, y) ↔ A(x, y+1).
     * Round 3: ZZ on the same vertical edges. */
    idx = 0;
    for (int y = 0; y < Ly; ++y) {
        for (int x = 0; x < Lx; ++x) {
            int qb = q_idx(x, y, 1, Lx);
            int qa_up = q_idx(x, wrap(y + 1, Ly), 0, Lx);
            irrep_pauli_set(&out->rounds[2].meas[idx], qb,    IRREP_PAULI_LETTER_X);
            irrep_pauli_set(&out->rounds[2].meas[idx], qa_up, IRREP_PAULI_LETTER_X);
            irrep_pauli_set(&out->rounds[3].meas[idx], qb,    IRREP_PAULI_LETTER_Z);
            irrep_pauli_set(&out->rounds[3].meas[idx], qa_up, IRREP_PAULI_LETTER_Z);
            ++idx;
        }
    }
    return IRREP_OK;
}
