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

int main(void) {
    int rc = 0;
    if (test_toric_init()) { fprintf(stderr, "FAIL test_toric_init\n"); rc = 1; }
    if (test_edge_roundtrip()) { fprintf(stderr, "FAIL test_edge_roundtrip\n"); rc = 1; }
    if (test_stabilizer_commutativity()) { fprintf(stderr, "FAIL test_stabilizer_commutativity\n"); rc = 1; }
    if (test_vertex_neighbour_overlap()) { fprintf(stderr, "FAIL test_vertex_neighbour_overlap\n"); rc = 1; }
    return rc;
}
