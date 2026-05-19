/* SPDX-License-Identifier: MIT */
/* Tests for `<irrep/color_code_3d.h>` — `[[8, 3, 2]]` cubic 3D color code. */
#include "harness.h"
#include <irrep/color_code_3d.h>
#include <irrep/css_code.h>
#include <irrep/qec_distance.h>
#include <irrep/stabilizer_group.h>

#include <stdio.h>

static int test_cube_8_3_2_structure(void) {
    IRREP_TEST_START("color_3d_cube_8_3_2_structure");

    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_3d_cube_8_3_2(&cs) == IRREP_OK);
    IRREP_ASSERT(cs.n == IRREP_COLOR_3D_CUBE_N);
    IRREP_ASSERT(cs.H_X.n_rows == IRREP_COLOR_3D_CUBE_N_X_STABS);
    IRREP_ASSERT(cs.H_Z.n_rows == IRREP_COLOR_3D_CUBE_N_Z_STABS);

    /* CSS orthogonality: H_X · H_Zᵀ = 0 (mod 2). */
    IRREP_ASSERT(irrep_css_code_verify(&cs) == IRREP_OK);

    /* k = 3 by the [[8, 3, 2]] name (3 logical qubits in a single cube). */
    IRREP_ASSERT(irrep_css_code_logical_qubits(&cs) == 3);

    /* Full pairwise commutativity. */
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);
    IRREP_ASSERT(g.n_generators == IRREP_COLOR_3D_CUBE_N_X_STABS
                                 + IRREP_COLOR_3D_CUBE_N_Z_STABS);
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);

    /* Weights: each face stab is weight 4; the volume stabs are weight 8. */
    int weights_X[4] = {0};
    for (int r = 0; r < IRREP_COLOR_3D_CUBE_N_X_STABS; ++r) {
        for (int q = 0; q < IRREP_COLOR_3D_CUBE_N; ++q) {
            if (irrep_parity_matrix_get(&cs.H_X, r, q)) ++weights_X[r];
        }
    }
    IRREP_ASSERT(weights_X[0] == 4);
    IRREP_ASSERT(weights_X[1] == 4);
    IRREP_ASSERT(weights_X[2] == 4);
    IRREP_ASSERT(weights_X[3] == 8);

    int weights_Z = 0;
    for (int q = 0; q < IRREP_COLOR_3D_CUBE_N; ++q) {
        if (irrep_parity_matrix_get(&cs.H_Z, 0, q)) ++weights_Z;
    }
    IRREP_ASSERT(weights_Z == 8);

    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

static int test_cube_8_3_2_distance(void) {
    IRREP_TEST_START("color_3d_cube_8_3_2_distance");
    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_3d_cube_8_3_2(&cs) == IRREP_OK);
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);

    /* Brute-force enumerate Pauli supports up to weight 2.
     * The expected distance is 2 (= IRREP_COLOR_3D_CUBE_DISTANCE). */
    int d = irrep_qec_distance_brute(&g, /*max_weight*/ 2);
    IRREP_ASSERT(d == IRREP_COLOR_3D_CUBE_DISTANCE);

    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

static int test_cube_8_3_2_invalid_args(void) {
    IRREP_TEST_START("color_3d_cube_8_3_2_invalid_args");
    IRREP_ASSERT(irrep_color_3d_cube_8_3_2(NULL) == IRREP_ERR_INVALID_ARG);
    return IRREP_TEST_END();
}

/* [[15, 1, 3]] Reed-Muller / tetrahedral 3D color code. */
static int test_rm_15_1_3_structure(void) {
    IRREP_TEST_START("color_3d_rm_15_1_3_structure");

    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_3d_rm_15_1_3(&cs) == IRREP_OK);
    IRREP_ASSERT(cs.n == IRREP_COLOR_3D_RM15_N);
    IRREP_ASSERT(cs.H_X.n_rows == IRREP_COLOR_3D_RM15_N_X_STABS);
    IRREP_ASSERT(cs.H_Z.n_rows == IRREP_COLOR_3D_RM15_N_Z_STABS);

    /* CSS orthogonality. */
    IRREP_ASSERT(irrep_css_code_verify(&cs) == IRREP_OK);

    /* k = 1 by the [[15, 1, 3]] name. */
    IRREP_ASSERT(irrep_css_code_logical_qubits(&cs) == 1);

    /* Materialise + commutativity. */
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);
    IRREP_ASSERT(g.n_generators == IRREP_COLOR_3D_RM15_N_X_STABS
                                 + IRREP_COLOR_3D_RM15_N_Z_STABS);
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);

    /* Weight expectations.
     * H_X: all 4 rows weight 8 (one indicator per bit, 8 of 15 points
     *      have bit b set: those with the high half of the bit). */
    for (int r = 0; r < IRREP_COLOR_3D_RM15_N_X_STABS; ++r) {
        int w = 0;
        for (int q = 0; q < IRREP_COLOR_3D_RM15_N; ++q) {
            if (irrep_parity_matrix_get(&cs.H_X, r, q)) ++w;
        }
        IRREP_ASSERT(w == 8);
    }
    /* H_Z: rows 0..3 weight 8 (single-bit), rows 4..9 weight 4 (pair). */
    for (int r = 0; r < 4; ++r) {
        int w = 0;
        for (int q = 0; q < IRREP_COLOR_3D_RM15_N; ++q) {
            if (irrep_parity_matrix_get(&cs.H_Z, r, q)) ++w;
        }
        IRREP_ASSERT(w == 8);
    }
    for (int r = 4; r < 10; ++r) {
        int w = 0;
        for (int q = 0; q < IRREP_COLOR_3D_RM15_N; ++q) {
            if (irrep_parity_matrix_get(&cs.H_Z, r, q)) ++w;
        }
        IRREP_ASSERT(w == 4);
    }
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

static int test_rm_15_1_3_distance(void) {
    IRREP_TEST_START("color_3d_rm_15_1_3_distance");
    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_3d_rm_15_1_3(&cs) == IRREP_OK);
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);

    /* Brute-force enumerate up to weight 3. The expected distance is 3:
     *   - no weight-1 or weight-2 Pauli is in the centralizer minus
     *     the stabilizer group;
     *   - a weight-3 Pauli exists (the PG(3,2) line {0,1,2}). */
    int d = irrep_qec_distance_brute(&g, /*max_weight*/ 3);
    IRREP_ASSERT(d == IRREP_COLOR_3D_RM15_DISTANCE);

    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

/* PROOF that [[15, 1, 3]] Reed-Muller has asymmetric X/Z distances:
 *   d_X = 7  (auto-discovered X-only logical has weight 7)
 *   d_Z = 3  (auto-discovered Z-only logical has weight 3)
 *   d = min(d_X, d_Z) = 3 (matches the published distance)
 *
 * This asymmetry distinguishes the Reed-Muller [[15,1,3]] from
 * symmetric CSS codes like Steane (d_X = d_Z = 3). The Z-side
 * minimum-weight logical exists at weight 3 (a PG(3,2)-style line
 * support); the X-side requires weight 7 (= the minimum distance
 * of the punctured RM(1,4) = [15,4,8] code's complement coset).
 *
 * This is the structural feature that supports the transversal-T
 * gate: triorthogonality of H_X yields the asymmetric distance. */
static int test_rm_15_1_3_dX_dZ_asymmetry(void) {
    IRREP_TEST_START("color_3d_rm_15_1_3_dX_dZ_asymmetry");
    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_color_3d_rm_15_1_3(&cs) == IRREP_OK);
    irrep_pauli_t Lx, Lz;
    int d_X = irrep_css_code_compute_logical_X(&cs, 8, &Lx);
    int d_Z = irrep_css_code_compute_logical_Z(&cs, 4, &Lz);
    IRREP_ASSERT(d_X == 7);
    IRREP_ASSERT(d_Z == 3);
    IRREP_ASSERT(irrep_pauli_weight(&Lx) == 7);
    IRREP_ASSERT(irrep_pauli_weight(&Lz) == 3);
    /* X and Z logicals must anti-commute (encoded Pauli relation). */
    IRREP_ASSERT(!irrep_pauli_commute(&Lx, &Lz));
    irrep_pauli_free(&Lx);
    irrep_pauli_free(&Lz);
    irrep_css_code_free(&cs);
    return IRREP_TEST_END();
}

static int test_rm_15_1_3_invalid_args(void) {
    IRREP_TEST_START("color_3d_rm_15_1_3_invalid_args");
    IRREP_ASSERT(irrep_color_3d_rm_15_1_3(NULL) == IRREP_ERR_INVALID_ARG);
    return IRREP_TEST_END();
}

int main(void) {
    int rc = 0;
    rc |= test_cube_8_3_2_structure();
    rc |= test_cube_8_3_2_distance();
    rc |= test_cube_8_3_2_invalid_args();
    rc |= test_rm_15_1_3_structure();
    rc |= test_rm_15_1_3_distance();
    rc |= test_rm_15_1_3_invalid_args();
    rc |= test_rm_15_1_3_dX_dZ_asymmetry();
    return rc;
}
