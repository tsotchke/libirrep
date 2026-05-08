/* SPDX-License-Identifier: MIT */
/* Tests for the X-cube fracton code (Vijay-Haah-Fu 2016).
 *
 * Verifies:
 *  - Counts: n = 3·N, n_X = 3·N (3 planar stabs/vertex), n_Z = N (cube
 *    stabs), where N = Lx·Ly·Lz.
 *  - Each planar vertex stab has weight 4.
 *  - Each cube stab has weight 12.
 *  - CSS-orthogonality (every (planar-vertex-stab, cube-stab) pair
 *    shares 0 or 2 edges, by the corner-sharing argument).
 *  - Stabilizer-group materialisation passes pairwise commutativity.
 *  - Per-vertex redundancy A_v^{xy} · A_v^{xz} · A_v^{yz} = identity:
 *    the 3 planar stabs at a single vertex sum to zero in F₂ (each
 *    incident edge appears in exactly 2 of the 3 stabs).
 */
#include "harness.h"
#include <irrep/css_code.h>
#include <irrep/stabilizer_group.h>
#include <irrep/xcube_code.h>
#include <stdio.h>

static int row_weight(const irrep_parity_matrix_t *m, int row) {
    int w = 0;
    for (int c = 0; c < m->n_cols; ++c) if (irrep_parity_matrix_get(m, row, c)) ++w;
    return w;
}

static int test_counts(int Lx, int Ly, int Lz) {
    irrep_xcube_params_t p;
    if (irrep_xcube_init(&p, Lx, Ly, Lz) != IRREP_OK) return 1;
    int N = Lx * Ly * Lz;
    return (p.n_qubits == 3*N && p.n_vertex_stabs == 3*N && p.n_cube_stabs == N) ? 0 : 1;
}

static int test_weights(int Lx, int Ly, int Lz) {
    irrep_xcube_params_t p;
    irrep_xcube_init(&p, Lx, Ly, Lz);
    irrep_css_code_t c;
    if (irrep_xcube_build(&p, &c) != IRREP_OK) return 1;
    int rc = 0;
    for (int i = 0; i < c.H_X.n_rows; ++i) if (row_weight(&c.H_X, i) != 4)  rc = 1;
    for (int i = 0; i < c.H_Z.n_rows; ++i) if (row_weight(&c.H_Z, i) != 12) rc = 1;
    irrep_css_code_free(&c);
    return rc;
}

static int test_css_orth(int Lx, int Ly, int Lz) {
    irrep_xcube_params_t p;
    irrep_xcube_init(&p, Lx, Ly, Lz);
    irrep_css_code_t c;
    if (irrep_xcube_build(&p, &c) != IRREP_OK) return 1;
    int rc = irrep_css_code_verify(&c) == IRREP_OK ? 0 : 1;
    irrep_css_code_free(&c);
    return rc;
}

static int test_stabilizer_group(int Lx, int Ly, int Lz) {
    irrep_xcube_params_t p;
    irrep_xcube_init(&p, Lx, Ly, Lz);
    irrep_css_code_t c;
    if (irrep_xcube_build(&p, &c) != IRREP_OK) return 1;
    irrep_stabilizer_group_t g;
    int rc = 1;
    if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
        if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK) rc = 0;
        irrep_stabilizer_group_free(&g);
    }
    irrep_css_code_free(&c);
    return rc;
}

/* Verify per-vertex 3-stab redundancy. The xy-stab at v contains the 4 edges
 * of v in the xy-plane (2 x-edges, 2 y-edges); xz contains 2 x + 2 z;
 * yz contains 2 y + 2 z. Their F₂ sum = 0 since each edge type appears twice.
 * For each vertex, sum of the 3 stab rows must equal 0 over F₂. */
static int test_vertex_three_stab_redundancy(int Lx, int Ly, int Lz) {
    irrep_xcube_params_t p;
    irrep_xcube_init(&p, Lx, Ly, Lz);
    irrep_css_code_t c;
    if (irrep_xcube_build(&p, &c) != IRREP_OK) return 1;
    int N = Lx * Ly * Lz;
    int rc = 0;
    for (int v = 0; v < N; ++v) {
        int row_xy = v;
        int row_xz = N + v;
        int row_yz = 2 * N + v;
        for (int col = 0; col < c.n; ++col) {
            int parity =
                irrep_parity_matrix_get(&c.H_X, row_xy, col) ^
                irrep_parity_matrix_get(&c.H_X, row_xz, col) ^
                irrep_parity_matrix_get(&c.H_X, row_yz, col);
            if (parity != 0) rc = 1;
        }
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
        if (test_vertex_three_stab_redundancy(Lx, Ly, Lz))
            { fprintf(stderr, "FAIL test_vertex_three_stab_redundancy(%d,%d,%d)\n", Lx, Ly, Lz); rc = 1; }
    }
    return rc;
}
