/* SPDX-License-Identifier: MIT */
/* Tests for the rotated surface code [[d², 1, d]].
 *
 * Verifies (for d ∈ {2, 3, 5, 7}):
 *  - Stabilizer counts match the analytical formulas: n_q = d²,
 *    n_X + n_Z = d² - 1 (one logical qubit).
 *  - CSS orthogonality `H_X · H_Z^T = 0` holds.
 *  - Materialised stabilizer group has full pairwise commutativity.
 *  - For d=3 specifically: hand-checks the [[9, 1, 3]] stabilizer set
 *    (4 X + 4 Z, weight 4 bulk + weight 2 boundary).
 *  - For d=5: each row has weight 2 or 4; total set has expected cardinality.
 */
#include "harness.h"
#include <irrep/css_code.h>
#include <irrep/qec_distance.h>
#include <irrep/stabilizer_group.h>
#include <irrep/surface_code.h>
#include <stdio.h>

static int row_weight(const irrep_parity_matrix_t *m, int row) {
    int w = 0;
    for (int c = 0; c < m->n_cols; ++c) if (irrep_parity_matrix_get(m, row, c)) ++w;
    return w;
}

static int test_init_counts(void) {
    struct { int d, n_q, n_x, n_z; } cases[] = {
        { 2,  4,  1,  2 },
        { 3,  9,  4,  4 },
        { 5, 25, 12, 12 },
        { 7, 49, 24, 24 },
        { 9, 81, 40, 40 },
    };
    for (size_t k = 0; k < sizeof(cases)/sizeof(cases[0]); ++k) {
        irrep_surface_params_t p;
        if (irrep_surface_init(&p, cases[k].d) != IRREP_OK) return 1;
        if (p.n_qubits  != cases[k].n_q) return 1;
        if (p.n_X_stabs != cases[k].n_x) return 1;
        if (p.n_Z_stabs != cases[k].n_z) return 1;
        /* Encodes 1 logical qubit: n - (n_X + n_Z) = 1. */
        if (p.n_qubits - p.n_X_stabs - p.n_Z_stabs != 1) return 1;
    }
    return 0;
}

/* For each distance, build the code and verify CSS-orthogonality plus
 * full stabilizer-group commutativity. */
static int test_build_and_verify(int d) {
    irrep_surface_params_t p;
    if (irrep_surface_init(&p, d) != IRREP_OK) return 1;
    irrep_css_code_t c;
    if (irrep_surface_build(&p, &c) != IRREP_OK) return 1;
    int rc = 1;
    if (irrep_css_code_verify(&c) == IRREP_OK) {
        irrep_stabilizer_group_t g;
        if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
            if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK) {
                /* Surface code on a planar patch with rough/smooth
                 * boundaries encodes k = 1 logical qubit at any d ≥ 2. */
                if (irrep_css_code_logical_qubits(&c) == 1) rc = 0;
            }
            irrep_stabilizer_group_free(&g);
        }
    }
    irrep_css_code_free(&c);
    return rc;
}

/* Brute-force distance verification: at small d, the planar surface code
 * actually achieves distance d. For d ∈ {2, 3} the enumeration is fast
 * (n = 4 or 9 qubits); for d=5 it takes a few seconds. */
static int test_brute_distance(int d) {
    irrep_surface_params_t p;
    if (irrep_surface_init(&p, d) != IRREP_OK) return 1;
    irrep_css_code_t c;
    if (irrep_surface_build(&p, &c) != IRREP_OK) return 1;
    int rc = 1;
    irrep_stabilizer_group_t g;
    if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
        int dist = irrep_qec_distance_brute(&g, /*max_weight*/ d);
        if (dist == d) rc = 0;
        else fprintf(stderr, "  surface d=%d: brute-distance = %d\n", d, dist);
        irrep_stabilizer_group_free(&g);
    }
    irrep_css_code_free(&c);
    return rc;
}

static int test_d3_explicit_weights(void) {
    /* For d=3: 4 X-stabs, 4 Z-stabs.
     * Bulk: 2 X + 2 Z weight-4 each.
     * Boundary: 2 X (top/bottom) + 2 Z (left/right) weight-2 each. */
    irrep_surface_params_t p;
    irrep_surface_init(&p, 3);
    irrep_css_code_t c;
    if (irrep_surface_build(&p, &c) != IRREP_OK) return 1;
    int x_w4 = 0, x_w2 = 0, z_w4 = 0, z_w2 = 0;
    for (int i = 0; i < c.H_X.n_rows; ++i) {
        int w = row_weight(&c.H_X, i);
        if (w == 4) ++x_w4;
        else if (w == 2) ++x_w2;
        else { irrep_css_code_free(&c); return 1; }
    }
    for (int i = 0; i < c.H_Z.n_rows; ++i) {
        int w = row_weight(&c.H_Z, i);
        if (w == 4) ++z_w4;
        else if (w == 2) ++z_w2;
        else { irrep_css_code_free(&c); return 1; }
    }
    irrep_css_code_free(&c);
    if (x_w4 != 2 || x_w2 != 2) return 1;
    if (z_w4 != 2 || z_w2 != 2) return 1;
    return 0;
}

/* For d=5: bulk = 4²=16 plaquettes (8 X bulk + 8 Z bulk, all weight 4).
 * Boundary: 4 X + 4 Z (weight 2 each). Total 12 X + 12 Z. */
static int test_d5_weights(void) {
    irrep_surface_params_t p;
    irrep_surface_init(&p, 5);
    irrep_css_code_t c;
    if (irrep_surface_build(&p, &c) != IRREP_OK) return 1;
    int x_w4 = 0, x_w2 = 0, z_w4 = 0, z_w2 = 0;
    for (int i = 0; i < c.H_X.n_rows; ++i) {
        int w = row_weight(&c.H_X, i);
        if (w == 4) ++x_w4; else if (w == 2) ++x_w2;
    }
    for (int i = 0; i < c.H_Z.n_rows; ++i) {
        int w = row_weight(&c.H_Z, i);
        if (w == 4) ++z_w4; else if (w == 2) ++z_w2;
    }
    irrep_css_code_free(&c);
    if (x_w4 != 8 || x_w2 != 4) return 1;
    if (z_w4 != 8 || z_w2 != 4) return 1;
    return 0;
}

/* Logical X̄ and Z̄ must commute with every stabilizer and anti-commute
 * with each other.  */
static int test_logical_operators(int d) {
    irrep_surface_params_t p;
    if (irrep_surface_init(&p, d) != IRREP_OK) return 1;
    irrep_css_code_t c;
    if (irrep_surface_build(&p, &c) != IRREP_OK) return 1;
    irrep_stabilizer_group_t g;
    int rc = 1;
    if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
        irrep_pauli_t Lx, Lz;
        if (irrep_surface_logical_X(&p, &Lx) == IRREP_OK &&
            irrep_surface_logical_Z(&p, &Lz) == IRREP_OK) {
            rc = 0;
            /* Weight checks. */
            if (irrep_pauli_weight(&Lx) != d) rc = 1;
            if (irrep_pauli_weight(&Lz) != d) rc = 1;
            /* Commute with every stabilizer. */
            for (int i = 0; i < g.n_generators; ++i) {
                if (!irrep_pauli_commute(&g.gens[i], &Lx)) rc = 1;
                if (!irrep_pauli_commute(&g.gens[i], &Lz)) rc = 1;
            }
            /* Anti-commute with each other. */
            if (irrep_pauli_commute(&Lx, &Lz)) rc = 1;
        }
        irrep_pauli_free(&Lx);
        irrep_pauli_free(&Lz);
        irrep_stabilizer_group_free(&g);
    }
    irrep_css_code_free(&c);
    return rc;
}

int main(void) {
    int rc = 0;
    if (test_init_counts())          { fprintf(stderr, "FAIL test_init_counts\n"); rc = 1; }
    if (test_build_and_verify(2))    { fprintf(stderr, "FAIL test_build_and_verify(d=2)\n"); rc = 1; }
    if (test_build_and_verify(3))    { fprintf(stderr, "FAIL test_build_and_verify(d=3)\n"); rc = 1; }
    if (test_build_and_verify(5))    { fprintf(stderr, "FAIL test_build_and_verify(d=5)\n"); rc = 1; }
    if (test_build_and_verify(7))    { fprintf(stderr, "FAIL test_build_and_verify(d=7)\n"); rc = 1; }
    if (test_build_and_verify(9))    { fprintf(stderr, "FAIL test_build_and_verify(d=9)\n"); rc = 1; }
    if (test_brute_distance(2))      { fprintf(stderr, "FAIL test_brute_distance(d=2)\n"); rc = 1; }
    if (test_brute_distance(3))      { fprintf(stderr, "FAIL test_brute_distance(d=3)\n"); rc = 1; }
    if (test_d3_explicit_weights())  { fprintf(stderr, "FAIL test_d3_explicit_weights\n"); rc = 1; }
    if (test_d5_weights())           { fprintf(stderr, "FAIL test_d5_weights\n"); rc = 1; }
    if (test_logical_operators(3))   { fprintf(stderr, "FAIL test_logical_operators(d=3)\n"); rc = 1; }
    if (test_logical_operators(5))   { fprintf(stderr, "FAIL test_logical_operators(d=5)\n"); rc = 1; }
    if (test_logical_operators(7))   { fprintf(stderr, "FAIL test_logical_operators(d=7)\n"); rc = 1; }
    return rc;
}
