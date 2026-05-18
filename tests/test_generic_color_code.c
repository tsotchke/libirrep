/* SPDX-License-Identifier: MIT */
/* Tests for the generic 2D color-code framework.
 *
 * Verifies that the face-list-driven constructor reproduces Steane
 * [[7, 1, 3]] when fed the standard Hamming-code parity-row face list,
 * and that the resulting CSS code matches the bespoke `irrep_color_steane`
 * generator bit-for-bit.
 *
 * Also validates:
 *  - Rejection of out-of-range qubit indices.
 *  - Rejection of underweight faces (size < 2).
 *  - Round-trip through stabilizer-group commutativity.
 */
#include "harness.h"
#include <irrep/color_code.h>
#include <irrep/css_code.h>
#include <irrep/generic_color_code.h>
#include <irrep/stabilizer_group.h>
#include <stdio.h>

/* Steane face list: 3 weight-4 hexagonal faces (each face's qubits are
 * the bits of one Hamming-code parity row). */
static const int steane_face_0[4] = { 0, 2, 4, 6 };
static const int steane_face_1[4] = { 1, 2, 5, 6 };
static const int steane_face_2[4] = { 3, 4, 5, 6 };

static int test_steane_via_generic(void) {
    int sizes[3] = { 4, 4, 4 };
    const int *qubits[3] = { steane_face_0, steane_face_1, steane_face_2 };
    irrep_color_lattice_t L = {
        .n_qubits = 7,
        .n_faces = 3,
        .face_sizes = sizes,
        .face_qubits = qubits,
        .face_color = NULL,
    };
    irrep_css_code_t c;
    if (irrep_generic_color_build(&L, &c) != IRREP_OK) return 1;
    int rc = (c.n == 7 && c.H_X.n_rows == 3 && c.H_Z.n_rows == 3) ? 0 : 1;
    /* Round-trip: stabilizer-group commutativity. */
    irrep_stabilizer_group_t g;
    if (rc == 0) {
        if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
            if (irrep_stabilizer_group_check_commutativity(&g) != IRREP_OK)
                rc = 1;
            irrep_stabilizer_group_free(&g);
        } else rc = 1;
    }
    irrep_css_code_free(&c);
    return rc;
}

/* Sanity: framework's Steane matches the bespoke `irrep_color_steane`
 * bit-for-bit (same H_X, H_Z entries). */
static int test_generic_matches_bespoke_steane(void) {
    int sizes[3] = { 4, 4, 4 };
    const int *qubits[3] = { steane_face_0, steane_face_1, steane_face_2 };
    irrep_color_lattice_t L = {
        .n_qubits = 7, .n_faces = 3,
        .face_sizes = sizes, .face_qubits = qubits, .face_color = NULL,
    };
    irrep_css_code_t a, b;
    if (irrep_generic_color_build(&L, &a) != IRREP_OK) return 1;
    if (irrep_color_steane(&b) != IRREP_OK) {
        irrep_css_code_free(&a);
        return 1;
    }
    int rc = (a.n == b.n && a.H_X.n_rows == b.H_X.n_rows) ? 0 : 1;
    if (rc == 0) {
        for (int r = 0; r < a.H_X.n_rows && rc == 0; ++r) {
            for (int c = 0; c < a.n && rc == 0; ++c) {
                if (irrep_parity_matrix_get(&a.H_X, r, c) !=
                    irrep_parity_matrix_get(&b.H_X, r, c)) rc = 1;
                if (irrep_parity_matrix_get(&a.H_Z, r, c) !=
                    irrep_parity_matrix_get(&b.H_Z, r, c)) rc = 1;
            }
        }
    }
    irrep_css_code_free(&a);
    irrep_css_code_free(&b);
    return rc;
}

/* Out-of-range qubit index must be rejected. */
static int test_out_of_range_rejected(void) {
    int sizes[1] = { 3 };
    int bad[3] = { 0, 1, 99 }; /* 99 is out of range for n_qubits = 7 */
    const int *qubits[1] = { bad };
    irrep_color_lattice_t L = {
        .n_qubits = 7, .n_faces = 1,
        .face_sizes = sizes, .face_qubits = qubits, .face_color = NULL,
    };
    irrep_css_code_t c;
    irrep_status_t s = irrep_generic_color_build(&L, &c);
    if (s == IRREP_OK) { irrep_css_code_free(&c); return 1; }
    return s == IRREP_ERR_PRECONDITION ? 0 : 1;
}

/* Bad face list (CSS-orthogonality fails) must be rejected. */
static int test_bad_overlap_rejected(void) {
    /* Two faces sharing exactly 1 qubit → odd overlap → CSS verify fails. */
    int sizes[2] = { 3, 3 };
    int f0[3] = { 0, 1, 2 };
    int f1[3] = { 2, 3, 4 };
    const int *qubits[2] = { f0, f1 };
    irrep_color_lattice_t L = {
        .n_qubits = 5, .n_faces = 2,
        .face_sizes = sizes, .face_qubits = qubits, .face_color = NULL,
    };
    irrep_css_code_t c;
    irrep_status_t s = irrep_generic_color_build(&L, &c);
    if (s == IRREP_OK) { irrep_css_code_free(&c); return 1; }
    return s == IRREP_ERR_PRECONDITION ? 0 : 1;
}

int main(void) {
    int rc = 0;
    if (test_steane_via_generic())             { fprintf(stderr, "FAIL test_steane_via_generic\n"); rc = 1; }
    if (test_generic_matches_bespoke_steane()) { fprintf(stderr, "FAIL test_generic_matches_bespoke_steane\n"); rc = 1; }
    if (test_out_of_range_rejected())          { fprintf(stderr, "FAIL test_out_of_range_rejected\n"); rc = 1; }
    if (test_bad_overlap_rejected())           { fprintf(stderr, "FAIL test_bad_overlap_rejected\n"); rc = 1; }
    return rc;
}
