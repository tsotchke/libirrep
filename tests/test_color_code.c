/* SPDX-License-Identifier: MIT */
/* Tests for the triangular color code (Bombín–Martín-Delgado 2006).
 *
 * Verifies:
 *  - Steane [[7, 1, 3]] has 7 qubits, 3 X-stabs + 3 Z-stabs, weight 4 each.
 *  - CSS-orthogonality holds (every X-stab is F2-orthogonal to every Z-stab).
 *  - Round-trip through stabilizer-group abstraction passes commutativity.
 *  - The X- and Z-supports coincide face-by-face (the defining property
 *    of color codes — every face hosts both an X-stab and a Z-stab on
 *    the same qubits).
 */
#include "harness.h"
#include <irrep/color_code.h>
#include <irrep/css_code.h>
#include <irrep/stabilizer_group.h>
#include <stdio.h>

static int test_steane_shape(void) {
    irrep_css_code_t c;
    if (irrep_color_steane(&c) != IRREP_OK) return 1;
    int rc = (c.n == 7 && c.H_X.n_rows == 3 && c.H_Z.n_rows == 3) ? 0 : 1;
    /* Each row has weight 4 (Hamming-code parity rows). */
    for (int i = 0; i < 3 && rc == 0; ++i) {
        int wx = 0, wz = 0;
        for (int j = 0; j < 7; ++j) {
            if (irrep_parity_matrix_get(&c.H_X, i, j)) ++wx;
            if (irrep_parity_matrix_get(&c.H_Z, i, j)) ++wz;
        }
        if (wx != 4 || wz != 4) rc = 1;
    }
    irrep_css_code_free(&c);
    return rc;
}

static int test_steane_x_z_support_coincide(void) {
    /* Defining property of color codes: each face's X-support equals its
     * Z-support (both stabilizers live on the same boundary). */
    irrep_css_code_t c;
    if (irrep_color_steane(&c) != IRREP_OK) return 1;
    int rc = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 7; ++j) {
            int x = irrep_parity_matrix_get(&c.H_X, i, j);
            int z = irrep_parity_matrix_get(&c.H_Z, i, j);
            if (x != z) { rc = 1; break; }
        }
    }
    irrep_css_code_free(&c);
    return rc;
}

static int test_steane_css_verify(void) {
    irrep_css_code_t c;
    if (irrep_color_steane(&c) != IRREP_OK) return 1;
    irrep_status_t s = irrep_css_code_verify(&c);
    irrep_css_code_free(&c);
    return s == IRREP_OK ? 0 : 1;
}

static int test_steane_stabilizer_group(void) {
    irrep_css_code_t c;
    if (irrep_color_steane(&c) != IRREP_OK) return 1;
    irrep_stabilizer_group_t g;
    int rc = 1;
    if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
        if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK) rc = 0;
        irrep_stabilizer_group_free(&g);
    }
    irrep_css_code_free(&c);
    return rc;
}

/* [[15, 7, 3]] Hamming-CSS code from [15, 11, 3] punctured RM(2, 4):
 * structural counts + CSS orthogonality + k = 7 + brute-force d = 3. */
static int test_hamming_15_7_3(void) {
    irrep_css_code_t cs;
    if (irrep_color_hamming_15_7_3(&cs) != IRREP_OK) return 1;
    int rc = 0;
    if (cs.n != 15)              rc = 1;
    if (cs.H_X.n_rows != 4)      rc = 1;
    if (cs.H_Z.n_rows != 4)      rc = 1;
    if (irrep_css_code_verify(&cs) != IRREP_OK) rc = 1;
    if (irrep_css_code_logical_qubits(&cs) != 7) {
        fprintf(stderr, "  [[15,7,3]] k = %d (expected 7)\n",
                irrep_css_code_logical_qubits(&cs));
        rc = 1;
    }
    if (irrep_css_code_distance(&cs, 3) != 3) {
        fprintf(stderr, "  [[15,7,3]] d = %d (expected 3)\n",
                irrep_css_code_distance(&cs, 3));
        rc = 1;
    }
    irrep_css_code_free(&cs);
    return rc;
}

int main(void) {
    int rc = 0;
    if (test_steane_shape())                 { fprintf(stderr, "FAIL test_steane_shape\n"); rc = 1; }
    if (test_steane_x_z_support_coincide())  { fprintf(stderr, "FAIL test_steane_x_z_support_coincide\n"); rc = 1; }
    if (test_steane_css_verify())            { fprintf(stderr, "FAIL test_steane_css_verify\n"); rc = 1; }
    if (test_steane_stabilizer_group())      { fprintf(stderr, "FAIL test_steane_stabilizer_group\n"); rc = 1; }
    if (test_hamming_15_7_3())               { fprintf(stderr, "FAIL test_hamming_15_7_3\n"); rc = 1; }
    /* Steane logical operators: weight 3, commute with all stabs,
     * not in stab span, mutually anti-commute. */
    {
        irrep_css_code_t cs;
        if (irrep_color_steane(&cs) != IRREP_OK) { rc = 1; }
        else {
            irrep_pauli_t Lx, Lz;
            if (irrep_color_steane_logical_X(&Lx) != IRREP_OK) rc = 1;
            if (irrep_color_steane_logical_Z(&Lz) != IRREP_OK) rc = 1;
            if (irrep_pauli_weight(&Lx) != 3) rc = 1;
            if (irrep_pauli_weight(&Lz) != 3) rc = 1;
            irrep_stabilizer_group_t g;
            if (irrep_css_code_to_stabilizer_group(&cs, &g) != IRREP_OK) rc = 1;
            else {
                for (int i = 0; i < g.n_generators; ++i) {
                    if (!irrep_pauli_commute(&g.gens[i], &Lx)) rc = 1;
                    if (!irrep_pauli_commute(&g.gens[i], &Lz)) rc = 1;
                }
                if (irrep_pauli_in_stabilizer_span(&g, &Lx) != 0) rc = 1;
                if (irrep_pauli_in_stabilizer_span(&g, &Lz) != 0) rc = 1;
                if (irrep_pauli_commute(&Lx, &Lz)) rc = 1; /* anti-commute */
                irrep_stabilizer_group_free(&g);
            }
            irrep_pauli_free(&Lx);
            irrep_pauli_free(&Lz);
            irrep_css_code_free(&cs);
        }
        if (rc) fprintf(stderr, "FAIL Steane logical-operator proof\n");
    }
    return rc;
}
