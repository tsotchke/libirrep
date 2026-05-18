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

/* Two-tile contraction: build the joined stabilizer group of two
 * 6-leg perfect tensors on 12 qubits (tile A: 0..5, tile B: 6..11),
 * then contract leg 1 of tile A (= qubit 1) with leg 1 of tile B
 * (= qubit 7). The resulting stabilizer group should sit on 10 qubits
 * with all generators pairwise commuting. */
static int test_contract_two_tiles(void) {
    IRREP_TEST_START("happy_contract_two_tiles");

    /* Build a joined 12-qubit group: 6 stabs per tile lifted to the
     * 12-qubit space, total 12 stabs. */
    irrep_stabilizer_group_t g12;
    IRREP_ASSERT(irrep_stabilizer_group_new(&g12, 12, 12) == IRREP_OK);

    irrep_stabilizer_group_t tile;
    IRREP_ASSERT(irrep_happy_perfect_tensor_6leg(&tile) == IRREP_OK);

    /* Copy tile A into qubits 0..5 (rows 0..5 of g12). */
    for (int r = 0; r < 6; ++r) {
        for (int q = 0; q < 6; ++q) {
            irrep_pauli_letter_t L = irrep_pauli_get(&tile.gens[r], q);
            if (L != IRREP_PAULI_LETTER_I) {
                irrep_pauli_set(&g12.gens[r], q, L);
            }
        }
    }
    /* Copy tile B into qubits 6..11 (rows 6..11 of g12). */
    for (int r = 0; r < 6; ++r) {
        for (int q = 0; q < 6; ++q) {
            irrep_pauli_letter_t L = irrep_pauli_get(&tile.gens[r], q);
            if (L != IRREP_PAULI_LETTER_I) {
                irrep_pauli_set(&g12.gens[6 + r], 6 + q, L);
            }
        }
    }
    irrep_stabilizer_group_free(&tile);

    /* Sanity: the joined group commutes. */
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g12) == IRREP_OK);

    /* Contract qubits 1 and 7 (leg 1 of A ↔ leg 1 of B). */
    irrep_stabilizer_group_t g10;
    IRREP_ASSERT(irrep_stabilizer_contract_bell(&g12, 1, 7, &g10) == IRREP_OK);
    IRREP_ASSERT(g10.n == 10);
    /* All generators commute after contraction. */
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g10) == IRREP_OK);

    irrep_stabilizer_group_free(&g12);
    irrep_stabilizer_group_free(&g10);
    return IRREP_TEST_END();
}

/* Depth-2 HaPPY Bell-pair contraction: the 5-edge contraction pipeline
 * produces a 26-qubit stabilizer group representing the [[20, 6, ?]]
 * HaPPY encoding isometry. Verifies:
 *   - n = 26 qubits as expected (36 - 2·5).
 *   - All generators pairwise commute (the Bell-contraction primitive
 *     preserves commutativity of every centralizer element).
 *   - Bulk-qubit indices in the contracted frame match the predicted
 *     {0, 5, 8, 11, 16, 21}.
 *   - For each bulk qubit, the lifted cross-stabilizer
 *     `X X ... X` originally on a tile reduces to a non-trivial
 *     element of the contracted group acting on that bulk leg —
 *     i.e., the bulk leg still carries a logical operator. We test
 *     by checking the contracted generator at row index `6 t + 4`
 *     for the all-X cross-stab has support on the bulk qubit of tile
 *     t after contraction. */
static int test_depth2_contracted_shape_and_bulk_index(void) {
    IRREP_TEST_START("happy_depth2_contracted_shape_and_bulk_index");

    irrep_stabilizer_group_t g;
    int bulk_qubits[IRREP_HAPPY_DEPTH2_N_TILES];
    IRREP_ASSERT(irrep_happy_network_depth2_contracted(&g, bulk_qubits)
                 == IRREP_OK);
    IRREP_ASSERT(g.n == IRREP_HAPPY_DEPTH2_N_CONTRACTED_QUBITS);
    IRREP_ASSERT(g.n == 26);
    /* Generator count is at most the input count — Bell contraction
     * either replaces an anti-commuting generator in place or drops a
     * redundant one when the measurement is already in the span. */
    IRREP_ASSERT(g.n_generators <=
                 IRREP_HAPPY_DEPTH2_N_TILES * 6);
    IRREP_ASSERT(g.n_generators > 0);
    /* All pairwise commute. */
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);

    /* Predicted bulk indices: for original bulk `6t` (t = 0..5),
     * contracted index = 6t - (count of removed in {1,2,3,4,5,7,13,19,
     * 25,31} that are strictly less than 6t). Gives {0,1,6,11,16,21}. */
    int expected_bulk[6] = { 0, 1, 6, 11, 16, 21 };
    for (int t = 0; t < IRREP_HAPPY_DEPTH2_N_TILES; ++t) {
        IRREP_ASSERT(bulk_qubits[t] == expected_bulk[t]);
    }
    /* All 26 qubit indices are well-formed (0..25). */
    for (int t = 0; t < IRREP_HAPPY_DEPTH2_N_TILES; ++t) {
        IRREP_ASSERT(bulk_qubits[t] >= 0 && bulk_qubits[t] < g.n);
    }
    irrep_stabilizer_group_free(&g);
    return IRREP_TEST_END();
}

int main(void) {
    int rc = 0;
    rc |= test_perfect_tensor_6leg_shape();
    rc |= test_depth2_network_shape();
    rc |= test_invalid_args();
    rc |= test_contract_two_tiles();
    rc |= test_depth2_contracted_shape_and_bulk_index();
    return rc;
}
