/* SPDX-License-Identifier: MIT */
/* Tests for `<irrep/color_codes_2d.h>` — `[[19, 1, 5]]` hex color code.
 *
 * Verifies:
 *   - Construction returns IRREP_OK and produces a 19-qubit code with 9
 *     X-stabilizers and 9 Z-stabilizers (i.e. H_X = H_Z, the color-code
 *     property).
 *   - Each face has the expected weight (6 weight-4 boundary + 3 weight-6
 *     interior).
 *   - CSS orthogonality H_X · H_Zᵀ = 0 (mod 2).
 *   - Full pairwise stabilizer commutativity.
 *   - Stabilizer-group materialisation yields 18 generators on 19 qubits.
 *   - Code distance is 5 (= IRREP_COLOR_HEX_19_1_5_DISTANCE), verified by
 *     brute-force search over all Pauli errors up to weight 4 finding no
 *     undetectable non-stabilizer (the lowest-weight logical sits at
 *     weight 5).
 */
#include "harness.h"
#include <irrep/color_codes_2d.h>
#include <irrep/css_code.h>
#include <irrep/qec_distance.h>
#include <irrep/stabilizer_group.h>

#include <stdio.h>

static int test_hex_19_1_5_structure(void) {
    IRREP_TEST_START("color_hex_19_1_5_structure");

    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_hex_19_1_5(&cs) == IRREP_OK);
    IRREP_ASSERT(cs.n == IRREP_COLOR_HEX_19_1_5_N);
    IRREP_ASSERT(cs.H_X.n_rows == IRREP_COLOR_HEX_19_1_5_N_FACES);
    IRREP_ASSERT(cs.H_Z.n_rows == IRREP_COLOR_HEX_19_1_5_N_FACES);

    /* Face weights: 6 weight-4 boundary + 3 weight-6 interior.
     * Total weight-4 count and weight-6 count, summed against expected. */
    int n_w4 = 0, n_w6 = 0, n_other = 0;
    for (int f = 0; f < IRREP_COLOR_HEX_19_1_5_N_FACES; ++f) {
        int w = 0;
        for (int q = 0; q < IRREP_COLOR_HEX_19_1_5_N; ++q) {
            if (irrep_parity_matrix_get(&cs.H_X, f, q)) ++w;
        }
        if      (w == 4) ++n_w4;
        else if (w == 6) ++n_w6;
        else             ++n_other;
    }
    IRREP_ASSERT(n_w4 == 6);
    IRREP_ASSERT(n_w6 == 3);
    IRREP_ASSERT(n_other == 0);

    /* Color-code defining property: H_X = H_Z bit-for-bit. */
    int diff = 0;
    for (int f = 0; f < IRREP_COLOR_HEX_19_1_5_N_FACES; ++f) {
        for (int q = 0; q < IRREP_COLOR_HEX_19_1_5_N; ++q) {
            if (irrep_parity_matrix_get(&cs.H_X, f, q)
                != irrep_parity_matrix_get(&cs.H_Z, f, q)) ++diff;
        }
    }
    IRREP_ASSERT(diff == 0);

    /* CSS orthogonality and full stabilizer commutativity. */
    IRREP_ASSERT(irrep_css_code_verify(&cs) == IRREP_OK);

    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);
    IRREP_ASSERT(g.n_generators == 18);
    IRREP_ASSERT(g.n == 19);
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);

    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

static int test_hex_19_1_5_distance(void) {
    IRREP_TEST_START("color_hex_19_1_5_distance");

    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_hex_19_1_5(&cs) == IRREP_OK);

    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);

    /* irrep_qec_distance_brute: minimum weight of a Pauli operator that
     * commutes with all stabilizers but is NOT in the stabilizer group.
     * For [[19, 1, 5]] this must be 5. The brute search visits all Pauli
     * supports of weight 1..max_weight and returns the first hit. */
    int d = irrep_qec_distance_brute(&g, /*max_weight*/ 5);
    IRREP_ASSERT(d == IRREP_COLOR_HEX_19_1_5_DISTANCE);

    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

static int test_hex_19_1_5_logicals(void) {
    IRREP_TEST_START("color_hex_19_1_5_logicals");

    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_hex_19_1_5(&cs) == IRREP_OK);
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);

    irrep_pauli_t Lx, Lz;
    IRREP_ASSERT(irrep_color_hex_19_1_5_logical_X(&Lx) == IRREP_OK);
    IRREP_ASSERT(irrep_color_hex_19_1_5_logical_Z(&Lz) == IRREP_OK);

    /* Each canonical logical has weight equal to the distance. */
    IRREP_ASSERT(irrep_pauli_weight(&Lx) == IRREP_COLOR_HEX_19_1_5_DISTANCE);
    IRREP_ASSERT(irrep_pauli_weight(&Lz) == IRREP_COLOR_HEX_19_1_5_DISTANCE);

    /* Anti-commute on the shared corner D0. */
    IRREP_ASSERT(irrep_pauli_symp_inner(&Lx, &Lz) == 1);

    /* Commute with every stabilizer. */
    int anti_Lx = 0, anti_Lz = 0;
    for (int i = 0; i < g.n_generators; ++i) {
        if (!irrep_pauli_commute(&g.gens[i], &Lx)) ++anti_Lx;
        if (!irrep_pauli_commute(&g.gens[i], &Lz)) ++anti_Lz;
    }
    IRREP_ASSERT(anti_Lx == 0);
    IRREP_ASSERT(anti_Lz == 0);

    irrep_pauli_free(&Lx);
    irrep_pauli_free(&Lz);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

static int test_hex_19_1_5_invalid_args(void) {
    IRREP_TEST_START("color_hex_19_1_5_invalid_args");
    IRREP_ASSERT(irrep_color_hex_19_1_5(NULL) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_color_hex_19_1_5_logical_X(NULL) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_color_hex_19_1_5_logical_Z(NULL) == IRREP_ERR_INVALID_ARG);
    return IRREP_TEST_END();
}

/* ====================================================================
 * `[[17, 1, 5]]` 4.8.8 (square-octagon) color code tests.
 * ==================================================================== */

static int test_488_17_1_5_structure(void) {
    IRREP_TEST_START("color_488_17_1_5_structure");

    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_488_17_1_5(&cs) == IRREP_OK);
    IRREP_ASSERT(cs.n == IRREP_COLOR_488_17_1_5_N);
    IRREP_ASSERT(cs.H_X.n_rows == IRREP_COLOR_488_17_1_5_N_FACES);
    IRREP_ASSERT(cs.H_Z.n_rows == IRREP_COLOR_488_17_1_5_N_FACES);

    /* Face weights: 1 weight-8 central octagon + 7 weight-4 faces. */
    int n_w4 = 0, n_w8 = 0, n_other = 0;
    for (int f = 0; f < IRREP_COLOR_488_17_1_5_N_FACES; ++f) {
        int w = 0;
        for (int q = 0; q < IRREP_COLOR_488_17_1_5_N; ++q) {
            if (irrep_parity_matrix_get(&cs.H_X, f, q)) ++w;
        }
        if      (w == 4) ++n_w4;
        else if (w == 8) ++n_w8;
        else             ++n_other;
    }
    IRREP_ASSERT(n_w4 == 7);
    IRREP_ASSERT(n_w8 == 1);
    IRREP_ASSERT(n_other == 0);

    int diff = 0;
    for (int f = 0; f < IRREP_COLOR_488_17_1_5_N_FACES; ++f) {
        for (int q = 0; q < IRREP_COLOR_488_17_1_5_N; ++q) {
            if (irrep_parity_matrix_get(&cs.H_X, f, q)
                != irrep_parity_matrix_get(&cs.H_Z, f, q)) ++diff;
        }
    }
    IRREP_ASSERT(diff == 0);

    IRREP_ASSERT(irrep_css_code_verify(&cs) == IRREP_OK);

    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);
    IRREP_ASSERT(g.n_generators == 16);
    IRREP_ASSERT(g.n == 17);
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);

    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

static int test_488_17_1_5_distance(void) {
    IRREP_TEST_START("color_488_17_1_5_distance");
    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_488_17_1_5(&cs) == IRREP_OK);
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);
    int d = irrep_qec_distance_brute(&g, /*max_weight*/ 5);
    IRREP_ASSERT(d == IRREP_COLOR_488_17_1_5_DISTANCE);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

static int test_488_17_1_5_logicals(void) {
    IRREP_TEST_START("color_488_17_1_5_logicals");

    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_488_17_1_5(&cs) == IRREP_OK);
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);

    irrep_pauli_t Lx, Lz;
    IRREP_ASSERT(irrep_color_488_17_1_5_logical_X(&Lx) == IRREP_OK);
    IRREP_ASSERT(irrep_color_488_17_1_5_logical_Z(&Lz) == IRREP_OK);

    IRREP_ASSERT(irrep_pauli_weight(&Lx) == IRREP_COLOR_488_17_1_5_DISTANCE);
    IRREP_ASSERT(irrep_pauli_weight(&Lz) == IRREP_COLOR_488_17_1_5_DISTANCE);
    IRREP_ASSERT(irrep_pauli_symp_inner(&Lx, &Lz) == 1);

    int anti_Lx = 0, anti_Lz = 0;
    for (int i = 0; i < g.n_generators; ++i) {
        if (!irrep_pauli_commute(&g.gens[i], &Lx)) ++anti_Lx;
        if (!irrep_pauli_commute(&g.gens[i], &Lz)) ++anti_Lz;
    }
    IRREP_ASSERT(anti_Lx == 0);
    IRREP_ASSERT(anti_Lz == 0);

    irrep_pauli_free(&Lx);
    irrep_pauli_free(&Lz);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

static int test_488_17_1_5_invalid_args(void) {
    IRREP_TEST_START("color_488_17_1_5_invalid_args");
    IRREP_ASSERT(irrep_color_488_17_1_5(NULL) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_color_488_17_1_5_logical_X(NULL) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_color_488_17_1_5_logical_Z(NULL) == IRREP_ERR_INVALID_ARG);
    return IRREP_TEST_END();
}

int main(void) {
    int rc = 0;
    rc |= test_hex_19_1_5_structure();
    rc |= test_hex_19_1_5_distance();
    rc |= test_hex_19_1_5_logicals();
    rc |= test_hex_19_1_5_invalid_args();
    rc |= test_488_17_1_5_structure();
    rc |= test_488_17_1_5_distance();
    rc |= test_488_17_1_5_logicals();
    rc |= test_488_17_1_5_invalid_args();
    return rc;
}
