/* SPDX-License-Identifier: MIT */
/* Tests for the 3D toric code on the cubic lattice (Castelnovo-Chamon 2008).
 *
 * Verifies:
 *  - Counts: n = 3·Lx·Ly·Lz, n_X = Lx·Ly·Lz, n_Z = 3·Lx·Ly·Lz.
 *  - Each vertex stab has weight 6 (six incident edges).
 *  - Each face stab has weight 4.
 *  - CSS-orthogonality across all (vertex, face) pairs (each pair shares
 *    0 or 2 edges).
 *  - Stabilizer-group materialisation passes pairwise commutativity.
 *  - For L=2 (smallest non-trivial PBC torus), product of all vertex
 *    stabilizers acts as the identity (each edge appears in exactly two
 *    incident-vertex stars).
 */
#include "harness.h"
#include <irrep/css_code.h>
#include <irrep/stabilizer_group.h>
#include <irrep/toric3d.h>
#include <stdio.h>

static int row_weight(const irrep_parity_matrix_t *m, int row) {
    int w = 0;
    for (int c = 0; c < m->n_cols; ++c) if (irrep_parity_matrix_get(m, row, c)) ++w;
    return w;
}

static int test_counts(int Lx, int Ly, int Lz) {
    irrep_toric3d_params_t p;
    if (irrep_toric3d_init(&p, Lx, Ly, Lz) != IRREP_OK) return 1;
    int N = Lx * Ly * Lz;
    return (p.n_qubits == 3*N && p.n_vertex_stabs == N && p.n_face_stabs == 3*N) ? 0 : 1;
}

static int test_weights(int Lx, int Ly, int Lz) {
    irrep_toric3d_params_t p;
    irrep_toric3d_init(&p, Lx, Ly, Lz);
    irrep_css_code_t c;
    if (irrep_toric3d_build(&p, &c) != IRREP_OK) return 1;
    int rc = 0;
    for (int i = 0; i < c.H_X.n_rows; ++i) {
        if (row_weight(&c.H_X, i) != 6) rc = 1;
    }
    for (int i = 0; i < c.H_Z.n_rows; ++i) {
        if (row_weight(&c.H_Z, i) != 4) rc = 1;
    }
    irrep_css_code_free(&c);
    return rc;
}

static int test_css_orth(int Lx, int Ly, int Lz) {
    irrep_toric3d_params_t p;
    irrep_toric3d_init(&p, Lx, Ly, Lz);
    irrep_css_code_t c;
    if (irrep_toric3d_build(&p, &c) != IRREP_OK) return 1;
    int rc = irrep_css_code_verify(&c) == IRREP_OK ? 0 : 1;
    irrep_css_code_free(&c);
    return rc;
}

static int test_stabilizer_group(int Lx, int Ly, int Lz) {
    irrep_toric3d_params_t p;
    irrep_toric3d_init(&p, Lx, Ly, Lz);
    irrep_css_code_t c;
    if (irrep_toric3d_build(&p, &c) != IRREP_OK) return 1;
    irrep_stabilizer_group_t g;
    int rc = 1;
    if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
        if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK) rc = 0;
        irrep_stabilizer_group_free(&g);
    }
    irrep_css_code_free(&c);
    return rc;
}

/* Vertex-stabilizer redundancy: product of all vertex stars on a closed
 * torus equals identity (each edge appears in exactly 2 incident-vertex
 * stars). Equivalently, each column of H_X has even parity sum. */
static int test_vertex_redundancy(int Lx, int Ly, int Lz) {
    irrep_toric3d_params_t p;
    irrep_toric3d_init(&p, Lx, Ly, Lz);
    irrep_css_code_t c;
    if (irrep_toric3d_build(&p, &c) != IRREP_OK) return 1;
    int rc = 0;
    for (int col = 0; col < c.n; ++col) {
        int parity = 0;
        for (int row = 0; row < c.H_X.n_rows; ++row) {
            parity ^= irrep_parity_matrix_get(&c.H_X, row, col);
        }
        if (parity != 0) rc = 1;
    }
    irrep_css_code_free(&c);
    return rc;
}

int main(void) {
    int rc = 0;
    int sizes[][3] = { {2, 2, 2}, {3, 2, 2}, {3, 3, 3}, {4, 3, 2} };
    for (size_t k = 0; k < sizeof(sizes)/sizeof(sizes[0]); ++k) {
        int Lx = sizes[k][0], Ly = sizes[k][1], Lz = sizes[k][2];
        if (test_counts(Lx, Ly, Lz))
            { fprintf(stderr, "FAIL test_counts(%d,%d,%d)\n", Lx, Ly, Lz); rc = 1; }
        if (test_weights(Lx, Ly, Lz))
            { fprintf(stderr, "FAIL test_weights(%d,%d,%d)\n", Lx, Ly, Lz); rc = 1; }
        if (test_css_orth(Lx, Ly, Lz))
            { fprintf(stderr, "FAIL test_css_orth(%d,%d,%d)\n", Lx, Ly, Lz); rc = 1; }
        if (test_stabilizer_group(Lx, Ly, Lz))
            { fprintf(stderr, "FAIL test_stabilizer_group(%d,%d,%d)\n", Lx, Ly, Lz); rc = 1; }
        if (test_vertex_redundancy(Lx, Ly, Lz))
            { fprintf(stderr, "FAIL test_vertex_redundancy(%d,%d,%d)\n", Lx, Ly, Lz); rc = 1; }
    }
    return rc;
}
