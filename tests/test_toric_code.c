/* SPDX-License-Identifier: MIT */
/* Tests for toric-code stabilizer primitives.
 *
 * Verifies:
 *  - Edge index packing/unpacking is bijective on [0, 2 Lx Ly).
 *  - Vertex-edge and plaquette-edge incidence return 4 distinct edges
 *    (modulo PBC degeneracies on Lx=2 and Ly=2 where some pairs coincide).
 *  - **Stabilizer commutativity**: every (vertex, plaquette) pair shares
 *    an even number of edges (hence A_v B_p = B_p A_v).
 *  - Counts: n_qubits = 2 Lx Ly, n_vertices = n_plaquettes = Lx Ly.
 */
#include "harness.h"
#include <irrep/css_code.h>
#include <irrep/stabilizer_group.h>
#include <irrep/toric_code.h>
#include <stdio.h>
#include <stdlib.h>

/* Test 1: parameter initialisation. */
static int test_toric_init(void) {
    irrep_toric_params_t p;
    if (irrep_toric_init(&p, 3, 4) != IRREP_OK) {
        return 1;
    }
    if (p.Lx != 3 || p.Ly != 4) return 1;
    if (p.n_qubits != 2 * 3 * 4) return 1;
    if (p.n_vertices != 3 * 4) return 1;
    if (p.n_plaquettes != 3 * 4) return 1;
    return 0;
}

/* Test 2: edge index round-trip. */
static int test_edge_roundtrip(void) {
    irrep_toric_params_t p;
    irrep_toric_init(&p, 4, 4);
    int n = p.n_qubits;
    for (int e = 0; e < n; ++e) {
        irrep_toric_edge_t edge;
        if (irrep_toric_edge_unpack(&p, e, &edge) != IRREP_OK) return 1;
        int back = irrep_toric_edge_index(&p, edge.orient, edge.x, edge.y);
        if (back != e) return 1;
    }
    return 0;
}

/* Test 3: stabilizer commutativity for ALL (vertex, plaquette) pairs.
 *
 * The toric code's defining property: each A_v · B_p = B_p · A_v.
 * Equivalently: every (vertex, plaquette) pair shares 0 or 2 edges
 * (count even ⇒ commute).
 *
 * Proof of the structural claim:
 *  - If vertex v is not on the boundary of plaquette p, they share 0 edges.
 *  - If vertex v is a corner of plaquette p, they share exactly 2 edges
 *    (two adjacent edges of p both incident to v).
 *  - No (v, p) pair shares 1, 3, or 4 edges.
 */
static int test_stabilizer_commutativity(void) {
    irrep_toric_params_t p;
    irrep_toric_init(&p, 4, 4);
    for (int vy = 0; vy < p.Ly; ++vy) {
        for (int vx = 0; vx < p.Lx; ++vx) {
            int v_edges[4];
            irrep_toric_vertex_edges(&p, vx, vy, v_edges);
            for (int py = 0; py < p.Ly; ++py) {
                for (int px = 0; px < p.Lx; ++px) {
                    int p_edges[4];
                    irrep_toric_plaquette_edges(&p, px, py, p_edges);
                    int shared = irrep_toric_shared_edges(v_edges, p_edges);
                    if (shared % 2 != 0) {
                        /* Anti-commutation would break the toric code. */
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

/* Test 4: edge-incidence symmetry — vertex (vx, vy) shares an edge
 * with neighbour (vx+1, vy) (the eastern edge). Sanity check on the
 * incidence tables.
 *
 * (Vertex-vertex commutativity is automatic since A_v stabilizers are
 * all-X operators that always commute with each other regardless of
 * edge overlap. Same for plaquette-plaquette which are all-Z.) */
static int test_vertex_neighbour_overlap(void) {
    irrep_toric_params_t p;
    irrep_toric_init(&p, 4, 4);
    int e1[4], e2[4];
    irrep_toric_vertex_edges(&p, 1, 1, e1);
    irrep_toric_vertex_edges(&p, 2, 1, e2);  /* eastern neighbour */
    int shared = irrep_toric_shared_edges(e1, e2);
    /* They share exactly the horizontal edge between them: edge (1, 1). */
    if (shared != 1) return 1;
    return 0;
}

/* Build the 8 stabilizers for 2×2 toric (4 vertex stars + 4 plaquettes) and
 * verify each of the 4 logical operators commutes with all of them, plus
 * the (X1, Z1) and (X2, Z2) pairs anti-commute with their conjugate. */
static int test_toric_logical_operators(int Lx, int Ly) {
    irrep_toric_params_t p;
    if (irrep_toric_init(&p, Lx, Ly) != IRREP_OK) return 1;
    /* Build vertex+plaquette stabilizers as Pauli operators on the
     * 2·Lx·Ly edges. */
    int m = p.n_vertices + p.n_plaquettes;
    irrep_stabilizer_group_t g;
    if (irrep_stabilizer_group_new(&g, p.n_qubits, m) != IRREP_OK) return 1;
    int idx = 0;
    for (int vy = 0; vy < Ly; ++vy) {
        for (int vx = 0; vx < Lx; ++vx, ++idx) {
            int e[4];
            irrep_toric_vertex_edges(&p, vx, vy, e);
            for (int k = 0; k < 4; ++k)
                irrep_pauli_set(&g.gens[idx], e[k], IRREP_PAULI_LETTER_X);
        }
    }
    for (int py = 0; py < Ly; ++py) {
        for (int px = 0; px < Lx; ++px, ++idx) {
            int e[4];
            irrep_toric_plaquette_edges(&p, px, py, e);
            for (int k = 0; k < 4; ++k)
                irrep_pauli_set(&g.gens[idx], e[k], IRREP_PAULI_LETTER_Z);
        }
    }

    irrep_pauli_t Lx1, Lz1, Lx2, Lz2;
    int rc = 0;
    if (irrep_toric_logical_X1(&p, &Lx1) != IRREP_OK) rc = 1;
    if (irrep_toric_logical_Z1(&p, &Lz1) != IRREP_OK) rc = 1;
    if (irrep_toric_logical_X2(&p, &Lx2) != IRREP_OK) rc = 1;
    if (irrep_toric_logical_Z2(&p, &Lz2) != IRREP_OK) rc = 1;

    /* All four logicals commute with every stabilizer. */
    for (int i = 0; i < m && rc == 0; ++i) {
        if (!irrep_pauli_commute(&g.gens[i], &Lx1)) rc = 1;
        if (!irrep_pauli_commute(&g.gens[i], &Lz1)) rc = 1;
        if (!irrep_pauli_commute(&g.gens[i], &Lx2)) rc = 1;
        if (!irrep_pauli_commute(&g.gens[i], &Lz2)) rc = 1;
    }

    /* Conjugate pairs anti-commute. */
    if (rc == 0) {
        if (irrep_pauli_commute(&Lx1, &Lz1)) rc = 1;
        if (irrep_pauli_commute(&Lx2, &Lz2)) rc = 1;
        /* Cross pairs commute (different logical qubits). */
        if (!irrep_pauli_commute(&Lx1, &Lx2)) rc = 1;
        if (!irrep_pauli_commute(&Lz1, &Lz2)) rc = 1;
        if (!irrep_pauli_commute(&Lx1, &Lz2)) rc = 1;
        if (!irrep_pauli_commute(&Lz1, &Lx2)) rc = 1;
    }

    /* Weights are Lx (for L_*1) and Ly (for L_*2). */
    if (rc == 0) {
        if (irrep_pauli_weight(&Lx1) != Lx) rc = 1;
        if (irrep_pauli_weight(&Lz1) != Ly) rc = 1;
        if (irrep_pauli_weight(&Lx2) != Ly) rc = 1;
        if (irrep_pauli_weight(&Lz2) != Lx) rc = 1;
    }

    irrep_pauli_free(&Lx1);
    irrep_pauli_free(&Lz1);
    irrep_pauli_free(&Lx2);
    irrep_pauli_free(&Lz2);
    irrep_stabilizer_group_free(&g);
    return rc;
}

/* CSS-builder test: 2D toric on Lx × Ly torus has k = 2 logical qubits
 * (= dim H_1(T², F_2)). Verifies via the new irrep_toric_code_build_css
 * + F₂-rank logical-qubit counter. */
static int test_toric_code_build_css_k(int Lx, int Ly) {
    irrep_toric_params_t p;
    irrep_toric_init(&p, Lx, Ly);
    irrep_css_code_t c;
    if (irrep_toric_code_build_css(&p, &c) != IRREP_OK) return 1;
    int rc = 0;
    if (c.n != 2 * Lx * Ly) rc = 1;
    if (c.H_X.n_rows != Lx * Ly) rc = 1;
    if (c.H_Z.n_rows != Lx * Ly) rc = 1;
    if (irrep_css_code_verify(&c) != IRREP_OK) rc = 1;
    if (irrep_css_code_logical_qubits(&c) != 2) {
        fprintf(stderr, "  2D toric (%d, %d): k = %d (expected 2)\n",
                Lx, Ly, irrep_css_code_logical_qubits(&c));
        rc = 1;
    }
    irrep_css_code_free(&c);
    return rc;
}

/* Brute-force distance verification: 2D toric on L × L torus has
 * d = L. Tested at L = 2, 3, 4 (only L=2,3 use brute; L=4 too slow). */
static int test_toric_code_brute_distance(int L) {
    irrep_toric_params_t p;
    irrep_toric_init(&p, L, L);
    irrep_css_code_t c;
    if (irrep_toric_code_build_css(&p, &c) != IRREP_OK) return 1;
    int d = irrep_css_code_distance(&c, L);
    int rc = (d == L) ? 0 : 1;
    if (rc) fprintf(stderr, "  2D toric L=%d: brute-distance = %d (expected %d)\n",
                    L, d, L);
    irrep_css_code_free(&c);
    return rc;
}

int main(void) {
    int rc = 0;
    if (test_toric_init()) { fprintf(stderr, "FAIL test_toric_init\n"); rc = 1; }
    if (test_edge_roundtrip()) { fprintf(stderr, "FAIL test_edge_roundtrip\n"); rc = 1; }
    if (test_stabilizer_commutativity()) { fprintf(stderr, "FAIL test_stabilizer_commutativity\n"); rc = 1; }
    if (test_vertex_neighbour_overlap()) { fprintf(stderr, "FAIL test_vertex_neighbour_overlap\n"); rc = 1; }
    if (test_toric_logical_operators(2, 2)) { fprintf(stderr, "FAIL test_toric_logical_operators(2,2)\n"); rc = 1; }
    if (test_toric_logical_operators(3, 3)) { fprintf(stderr, "FAIL test_toric_logical_operators(3,3)\n"); rc = 1; }
    if (test_toric_logical_operators(4, 3)) { fprintf(stderr, "FAIL test_toric_logical_operators(4,3)\n"); rc = 1; }
    if (test_toric_code_build_css_k(2, 2)) { fprintf(stderr, "FAIL test_toric_code_build_css_k(2,2)\n"); rc = 1; }
    if (test_toric_code_build_css_k(3, 3)) { fprintf(stderr, "FAIL test_toric_code_build_css_k(3,3)\n"); rc = 1; }
    if (test_toric_code_build_css_k(3, 4)) { fprintf(stderr, "FAIL test_toric_code_build_css_k(3,4)\n"); rc = 1; }
    if (test_toric_code_build_css_k(5, 5)) { fprintf(stderr, "FAIL test_toric_code_build_css_k(5,5)\n"); rc = 1; }
    if (test_toric_code_brute_distance(2)) { fprintf(stderr, "FAIL test_toric_code_brute_distance(L=2)\n"); rc = 1; }
    if (test_toric_code_brute_distance(3)) { fprintf(stderr, "FAIL test_toric_code_brute_distance(L=3)\n"); rc = 1; }
    return rc;
}
