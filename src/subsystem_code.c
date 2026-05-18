/* SPDX-License-Identifier: MIT */
/** @file subsystem_code.c
 *  @brief Implementation of subsystem-code primitives. */
#include <irrep/subsystem_code.h>
#include <irrep/types.h>

#include <stdlib.h>

irrep_status_t
irrep_subsystem_code_new(irrep_subsystem_code_t *out,
                         int n_qubits, int n_gauge)
{
    if (out == NULL || n_qubits <= 0 || n_gauge <= 0) return IRREP_ERR_INVALID_ARG;
    out->n = n_qubits;
    out->n_gauge = n_gauge;
    out->gauge = (irrep_pauli_t *)calloc(n_gauge, sizeof(irrep_pauli_t));
    if (out->gauge == NULL) return IRREP_ERR_OUT_OF_MEMORY;
    for (int i = 0; i < n_gauge; ++i) {
        irrep_status_t s = irrep_pauli_new(&out->gauge[i], n_qubits);
        if (s != IRREP_OK) {
            for (int j = 0; j < i; ++j) irrep_pauli_free(&out->gauge[j]);
            free(out->gauge);
            out->gauge = NULL;
            return s;
        }
    }
    return IRREP_OK;
}

void
irrep_subsystem_code_free(irrep_subsystem_code_t *c)
{
    if (c == NULL) return;
    if (c->gauge) {
        for (int i = 0; i < c->n_gauge; ++i) irrep_pauli_free(&c->gauge[i]);
        free(c->gauge);
        c->gauge = NULL;
    }
    c->n = 0;
    c->n_gauge = 0;
}

int
irrep_subsystem_in_centraliser(const irrep_subsystem_code_t *c,
                               const irrep_pauli_t *p)
{
    if (c == NULL || p == NULL) return 0;
    for (int i = 0; i < c->n_gauge; ++i) {
        if (!irrep_pauli_commute(&c->gauge[i], p)) return 0;
    }
    return 1;
}

/* Bacon-Shor [[9, 1, 3]] subsystem code.
 *
 * Layout: 3×3 grid, qubit q(r, c) = 3r + c.
 *
 * Gauge group (12 weight-2 generators):
 *   - 6 X-gauges (horizontal nearest-neighbour XX in each row):
 *       X(q(r, 0))·X(q(r, 1)) and X(q(r, 1))·X(q(r, 2)) for r = 0, 1, 2.
 *   - 6 Z-gauges (vertical nearest-neighbour ZZ in each column):
 *       Z(q(0, c))·Z(q(1, c)) and Z(q(1, c))·Z(q(2, c)) for c = 0, 1, 2.
 *
 * Stabilizer subgroup (4 weight-6 generators, in the centre of the gauge
 * group; products of 3 gauge generators each):
 *   - S_X1 = X over columns 0+1 = X_0 X_1 X_3 X_4 X_6 X_7
 *   - S_X2 = X over columns 1+2 = X_1 X_2 X_4 X_5 X_7 X_8
 *   - S_Z1 = Z over rows 0+1    = Z_0 Z_1 Z_2 Z_3 Z_4 Z_5
 *   - S_Z2 = Z over rows 1+2    = Z_3 Z_4 Z_5 Z_6 Z_7 Z_8
 *
 * Logical operators:
 *   - L_X = X over a single column, e.g. X_0 X_3 X_6.
 *   - L_Z = Z over a single row, e.g. Z_0 Z_1 Z_2.
 */
irrep_status_t
irrep_subsystem_bacon_shor_9_1_3(irrep_subsystem_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_subsystem_code_new(out, 9, 12);
    if (s != IRREP_OK) return s;

    int idx = 0;
    /* X-gauges: row r, qubits 3r+0 -- 3r+1 and 3r+1 -- 3r+2. */
    for (int r = 0; r < 3; ++r) {
        irrep_pauli_set(&out->gauge[idx], 3*r + 0, IRREP_PAULI_LETTER_X);
        irrep_pauli_set(&out->gauge[idx], 3*r + 1, IRREP_PAULI_LETTER_X);
        ++idx;
        irrep_pauli_set(&out->gauge[idx], 3*r + 1, IRREP_PAULI_LETTER_X);
        irrep_pauli_set(&out->gauge[idx], 3*r + 2, IRREP_PAULI_LETTER_X);
        ++idx;
    }
    /* Z-gauges: col c, qubits c -- 3+c and 3+c -- 6+c. */
    for (int c = 0; c < 3; ++c) {
        irrep_pauli_set(&out->gauge[idx], 0 + c, IRREP_PAULI_LETTER_Z);
        irrep_pauli_set(&out->gauge[idx], 3 + c, IRREP_PAULI_LETTER_Z);
        ++idx;
        irrep_pauli_set(&out->gauge[idx], 3 + c, IRREP_PAULI_LETTER_Z);
        irrep_pauli_set(&out->gauge[idx], 6 + c, IRREP_PAULI_LETTER_Z);
        ++idx;
    }
    return IRREP_OK;
}

/* Smallest 3D Bombín gauge color code on a tetrahedral 3-colex.
 *
 * The 4 qubits sit at the 4 vertices of K_4 (the complete graph on
 * 4 vertices, = the boundary of a 3-simplex). The 6 edges are the
 * 6 unordered pairs (i, j) with i < j:
 *   edge 0: {0,1}   edge 1: {0,2}   edge 2: {0,3}
 *   edge 3: {1,2}   edge 4: {1,3}   edge 5: {2,3}
 *
 * Gauge generators: weight-2 X (XX) and weight-2 Z (ZZ) on each edge,
 * for 12 generators total. This is the symmetric "doubled-edge"
 * variant of the Bombín 2010 gauge color code on the smallest
 * tetrahedral 3-colex.
 */
static const int kK4Edges[6][2] = {
    {0, 1}, {0, 2}, {0, 3},
    {1, 2}, {1, 3}, {2, 3},
};

irrep_status_t
irrep_subsystem_bombin_3d_tetrahedron(irrep_subsystem_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_subsystem_code_new(out, /*n=*/4, /*n_gauge=*/12);
    if (s != IRREP_OK) return s;

    int idx = 0;
    /* 6 XX-edges. */
    for (int e = 0; e < 6; ++e) {
        irrep_pauli_set(&out->gauge[idx], kK4Edges[e][0], IRREP_PAULI_LETTER_X);
        irrep_pauli_set(&out->gauge[idx], kK4Edges[e][1], IRREP_PAULI_LETTER_X);
        ++idx;
    }
    /* 6 ZZ-edges. */
    for (int e = 0; e < 6; ++e) {
        irrep_pauli_set(&out->gauge[idx], kK4Edges[e][0], IRREP_PAULI_LETTER_Z);
        irrep_pauli_set(&out->gauge[idx], kK4Edges[e][1], IRREP_PAULI_LETTER_Z);
        ++idx;
    }
    return IRREP_OK;
}
