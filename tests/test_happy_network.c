/* SPDX-License-Identifier: MIT */
/* Tests for `<irrep/happy_network.h>`: 6-leg perfect tensor + depth-2
 * uncontracted joined stabilizer group. */
#include "harness.h"
#include <irrep/happy_network.h>
#include <irrep/stabilizer_group.h>

#include <stdio.h>

static int test_perfect_tensor_6leg_shape(void) {
    IRREP_TEST_START("happy_perfect_tensor_6leg_shape");
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_happy_perfect_tensor_6leg(&g) == IRREP_OK);
    IRREP_ASSERT(g.n == 6);
    IRREP_ASSERT(g.n_generators == 6);

    /* Boundary stabilizers (g_1..g_4) have weight 4 — same as the
     * [[5, 1, 3]] code generators, with I on leg 0. */
    for (int i = 0; i < 4; ++i) {
        IRREP_ASSERT(irrep_pauli_weight(&g.gens[i]) == 4);
    }
    /* Cross-stabilizers (g_5, g_6) have weight 6 (all-X and all-Z). */
    IRREP_ASSERT(irrep_pauli_weight(&g.gens[4]) == 6);
    IRREP_ASSERT(irrep_pauli_weight(&g.gens[5]) == 6);

    /* All 6 generators pairwise commute. */
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);

    irrep_stabilizer_group_free(&g);
    return IRREP_TEST_END();
}

static int test_depth2_network_shape(void) {
    IRREP_TEST_START("happy_network_depth2_shape");
    irrep_stabilizer_group_t g;
    int pairs[2 * IRREP_HAPPY_DEPTH2_N_CONTRACTIONS];
    IRREP_ASSERT(irrep_happy_network_depth2(&g, pairs) == IRREP_OK);
    IRREP_ASSERT(g.n == IRREP_HAPPY_DEPTH2_N_QUBITS);
    IRREP_ASSERT(g.n_generators ==
                 IRREP_HAPPY_DEPTH2_N_TILES * 6);

    /* All 36 generators pairwise commute (tiles are disjoint). */
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);

    /* Contraction pairs: central qubit (k+1) ↔ tile-k qubit (6(k+1)+1). */
    for (int k = 0; k < IRREP_HAPPY_DEPTH2_N_CONTRACTIONS; ++k) {
        IRREP_ASSERT(pairs[2 * k + 0] == (k + 1));
        IRREP_ASSERT(pairs[2 * k + 1] == 6 * (k + 1) + 1);
    }

    /* Generator weights per tile: same as the 6-leg perfect tensor —
     * four weight-4 and two weight-6, repeated for each of 6 tiles. */
    int w4 = 0, w6 = 0, other = 0;
    for (int r = 0; r < g.n_generators; ++r) {
        int w = irrep_pauli_weight(&g.gens[r]);
        if      (w == 4) ++w4;
        else if (w == 6) ++w6;
        else             ++other;
    }
    IRREP_ASSERT(w4 == IRREP_HAPPY_DEPTH2_N_TILES * 4);   /* 24 */
    IRREP_ASSERT(w6 == IRREP_HAPPY_DEPTH2_N_TILES * 2);   /* 12 */
    IRREP_ASSERT(other == 0);

    irrep_stabilizer_group_free(&g);
    return IRREP_TEST_END();
}

static int test_invalid_args(void) {
    IRREP_TEST_START("happy_network_invalid_args");
    IRREP_ASSERT(irrep_happy_perfect_tensor_6leg(NULL) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_happy_network_depth2(NULL, NULL) == IRREP_ERR_INVALID_ARG);
    return IRREP_TEST_END();
}

int main(void) {
    int rc = 0;
    rc |= test_perfect_tensor_6leg_shape();
    rc |= test_depth2_network_shape();
    rc |= test_invalid_args();
    return rc;
}
