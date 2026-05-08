/* SPDX-License-Identifier: MIT */
/** @file honeycomb_floquet.c
 *  @brief Hastings-Haah honeycomb Floquet code construction. */
#include <irrep/honeycomb_floquet.h>
#include <irrep/types.h>

#include <stddef.h>

/* Linear index of vertex (x, y, sublattice) on the Lx × Ly torus. */
static inline int hc_qubit(int x, int y, int s, int Lx, int Ly)
{
    (void)Ly;
    return 2 * (y * Lx + x) + s;
}

static inline int wrap(int v, int L) {
    int r = v % L;
    if (r < 0) r += L;
    return r;
}

irrep_status_t
irrep_honeycomb_floquet_build(int Lx, int Ly, irrep_floquet_code_t *out)
{
    if (out == NULL || Lx <= 0 || Ly <= 0) return IRREP_ERR_INVALID_ARG;
    int n_qubits = 2 * Lx * Ly;
    int edges_per_color = Lx * Ly;
    irrep_status_t s = irrep_floquet_code_new(out, n_qubits, 3);
    if (s != IRREP_OK) return s;
    for (int t = 0; t < 3; ++t) {
        s = irrep_floquet_round_alloc(out, t, edges_per_color);
        if (s != IRREP_OK) {
            irrep_floquet_code_free(out);
            return s;
        }
    }

    /* One measurement per edge per color. The triple of edges incident
     * to each A-vertex is enumerated below; we walk the unit cells and
     * append one edge of each color per cell, which exhausts the
     * 3·Lx·Ly edges of the lattice exactly. */
    int idx_R = 0, idx_G = 0, idx_B = 0;
    for (int y = 0; y < Ly; ++y) {
        for (int x = 0; x < Lx; ++x) {
            int qa = hc_qubit(x, y, 0, Lx, Ly);

            /* Red edge: A(x,y) ↔ B(x,y) */
            int qbR = hc_qubit(x, y, 1, Lx, Ly);
            irrep_pauli_set(&out->rounds[0].meas[idx_R], qa,  IRREP_PAULI_LETTER_X);
            irrep_pauli_set(&out->rounds[0].meas[idx_R], qbR, IRREP_PAULI_LETTER_X);
            ++idx_R;

            /* Green edge: A(x, y) ↔ B(x-1, y) */
            int qbG = hc_qubit(wrap(x - 1, Lx), y, 1, Lx, Ly);
            irrep_pauli_set(&out->rounds[1].meas[idx_G], qa,  IRREP_PAULI_LETTER_Y);
            irrep_pauli_set(&out->rounds[1].meas[idx_G], qbG, IRREP_PAULI_LETTER_Y);
            ++idx_G;

            /* Blue edge: A(x, y) ↔ B(x, y-1) */
            int qbB = hc_qubit(x, wrap(y - 1, Ly), 1, Lx, Ly);
            irrep_pauli_set(&out->rounds[2].meas[idx_B], qa,  IRREP_PAULI_LETTER_Z);
            irrep_pauli_set(&out->rounds[2].meas[idx_B], qbB, IRREP_PAULI_LETTER_Z);
            ++idx_B;
        }
    }
    return IRREP_OK;
}
