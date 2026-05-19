/* SPDX-License-Identifier: MIT */
/* Tests for the HaPPY-code [[5, 1, 3]] perfect-tensor primitive.
 *
 * Verifies:
 *  - 5 qubits, 4 stabilizer generators.
 *  - Each generator has weight 4 (one I per row, four non-I letters).
 *  - All 4 generators pairwise commute (the defining stabilizer property).
 *  - The encoded logical operators X̄ = X X X X X and Z̄ = Z Z Z Z Z both
 *    lie in the centraliser of the stabilizer group (commute with every
 *    generator) and anti-commute with each other on each qubit, giving
 *    a total of 5 anti-commuting positions → odd → anti-commute. (5 is
 *    odd so X̄ and Z̄ anti-commute, defining the encoded qubit's algebra.)
 *  - Every individual qubit error is detected by at least one stabilizer
 *    (this is what makes the code distance ≥ 2). For the [[5, 1, 3]]
 *    code, all weight-1 X, Y, Z errors must produce at least one
 *    syndrome bit.
 */
#include "harness.h"
#include <irrep/happy_code.h>
#include <irrep/qec_distance.h>
#include <irrep/stabilizer_group.h>
#include <stdio.h>

static int test_shape(void) {
    irrep_stabilizer_group_t g;
    if (irrep_happy_perfect_tensor_5_1_3(&g) != IRREP_OK) return 1;
    int rc = (g.n == 5 && g.n_generators == 4) ? 0 : 1;
    /* Each generator weight = 4. */
    for (int i = 0; i < 4 && rc == 0; ++i) {
        if (irrep_pauli_weight(&g.gens[i]) != 4) rc = 1;
    }
    irrep_stabilizer_group_free(&g);
    return rc;
}

static int test_pairwise_commute(void) {
    irrep_stabilizer_group_t g;
    if (irrep_happy_perfect_tensor_5_1_3(&g) != IRREP_OK) return 1;
    int rc = irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK ? 0 : 1;
    irrep_stabilizer_group_free(&g);
    return rc;
}

/* Logical operators X̄ = X⊗5 and Z̄ = Z⊗5 must commute with every stabilizer. */
static int test_logical_in_centraliser(void) {
    irrep_stabilizer_group_t g;
    if (irrep_happy_perfect_tensor_5_1_3(&g) != IRREP_OK) return 1;
    irrep_pauli_t Xbar, Zbar;
    irrep_pauli_new(&Xbar, 5);
    irrep_pauli_new(&Zbar, 5);
    for (int q = 0; q < 5; ++q) {
        irrep_pauli_set(&Xbar, q, IRREP_PAULI_LETTER_X);
        irrep_pauli_set(&Zbar, q, IRREP_PAULI_LETTER_Z);
    }
    int rc = 0;
    for (int i = 0; i < 4; ++i) {
        if (!irrep_pauli_commute(&g.gens[i], &Xbar)) rc = 1;
        if (!irrep_pauli_commute(&g.gens[i], &Zbar)) rc = 1;
    }
    /* X̄ and Z̄ must anti-commute (5 anti-commute positions, all qubits). */
    if (irrep_pauli_commute(&Xbar, &Zbar)) rc = 1;
    irrep_pauli_free(&Xbar);
    irrep_pauli_free(&Zbar);
    irrep_stabilizer_group_free(&g);
    return rc;
}

/* Single-qubit error detection: every weight-1 X, Y, or Z error on any
 * qubit must trigger at least one stabilizer (non-zero syndrome). This
 * proves the code has distance ≥ 2 (and combined with the cyclic
 * structure, distance = 3). */
static int test_single_qubit_errors_detected(void) {
    irrep_stabilizer_group_t g;
    if (irrep_happy_perfect_tensor_5_1_3(&g) != IRREP_OK) return 1;
    int rc = 0;
    for (int q = 0; q < 5; ++q) {
        for (int letter = 1; letter <= 3; ++letter) { /* X, Y, Z */
            irrep_pauli_t err;
            irrep_pauli_new(&err, 5);
            irrep_pauli_set(&err, q, (irrep_pauli_letter_t)letter);
            uint64_t syndrome = 0;
            irrep_stabilizer_syndrome(&g, &err, &syndrome);
            if (syndrome == 0) rc = 1; /* uncaught error → distance < 2 */
            irrep_pauli_free(&err);
        }
    }
    irrep_stabilizer_group_free(&g);
    return rc;
}

/* The [[5, 1, 3]] perfect-tensor code has k = 1 logical qubit. The new
 * symplectic-rank primitive derives this directly from the 4 stabilizer
 * generators (n = 5, rank = 4 in the 2n = 10 dim symplectic space). */
static int test_logical_qubits(void) {
    irrep_stabilizer_group_t g;
    if (irrep_happy_perfect_tensor_5_1_3(&g) != IRREP_OK) return 1;
    int k = irrep_stabilizer_group_n_logical_qubits(&g);
    irrep_stabilizer_group_free(&g);
    if (k != 1) {
        fprintf(stderr, "  [[5,1,3]] k = %d (expected 1)\n", k);
        return 1;
    }
    return 0;
}

/* Brute-force distance verification: [[5, 1, 3]] has d = 3. Weight-3
 * enumeration is trivial (n=5). */
static int test_brute_distance(void) {
    irrep_stabilizer_group_t g;
    if (irrep_happy_perfect_tensor_5_1_3(&g) != IRREP_OK) return 1;
    int d = irrep_qec_distance_brute(&g, /*max_weight*/ 3);
    irrep_stabilizer_group_free(&g);
    if (d != 3) {
        fprintf(stderr, "  [[5,1,3]] brute-distance = %d (expected 3)\n", d);
        return 1;
    }
    return 0;
}

int main(void) {
    int rc = 0;
    if (test_shape())                       { fprintf(stderr, "FAIL test_shape\n"); rc = 1; }
    if (test_pairwise_commute())            { fprintf(stderr, "FAIL test_pairwise_commute\n"); rc = 1; }
    if (test_logical_in_centraliser())      { fprintf(stderr, "FAIL test_logical_in_centraliser\n"); rc = 1; }
    if (test_single_qubit_errors_detected()){ fprintf(stderr, "FAIL test_single_qubit_errors_detected\n"); rc = 1; }
    if (test_logical_qubits())              { fprintf(stderr, "FAIL test_logical_qubits\n"); rc = 1; }
    if (test_brute_distance())              { fprintf(stderr, "FAIL test_brute_distance\n"); rc = 1; }
    return rc;
}
