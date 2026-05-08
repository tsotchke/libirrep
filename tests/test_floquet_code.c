/* SPDX-License-Identifier: MIT */
/* Tests for the Floquet-code measurement-schedule abstraction.
 *
 * Verifies:
 *  - Round allocation and per-round Pauli set/get round-trips.
 *  - Intra-round commutativity check passes when measurements commute,
 *    fails when they anti-commute.
 *  - A 4-qubit, 3-round Floquet schedule mimicking the honeycomb code
 *    (X-bonds, Y-bonds, Z-bonds on the 4-cycle) satisfies:
 *      * Each round's measurements pairwise commute (intra-round).
 *      * Measurements from different rounds anti-commute pairwise on
 *        the qubits they share (the defining Floquet dynamics).
 *  - Materialising a single round as a stabilizer group reproduces the
 *    measurement supports.
 */
#include "harness.h"
#include <irrep/floquet_code.h>
#include <irrep/stabilizer_group.h>
#include <stdio.h>

static int test_floquet_alloc(void) {
    irrep_floquet_code_t f;
    if (irrep_floquet_code_new(&f, 4, 3) != IRREP_OK) return 1;
    int rc = (f.n_qubits == 4 && f.n_rounds == 3) ? 0 : 1;
    /* Allocate 2 measurements in round 1. */
    if (irrep_floquet_round_alloc(&f, 1, 2) != IRREP_OK) rc = 1;
    if (f.rounds[1].n_meas != 2) rc = 1;
    /* Set the Pauli letters. */
    irrep_pauli_set(&f.rounds[1].meas[0], 0, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&f.rounds[1].meas[0], 1, IRREP_PAULI_LETTER_X);
    if (irrep_pauli_get(&f.rounds[1].meas[0], 0) != IRREP_PAULI_LETTER_X) rc = 1;
    if (irrep_pauli_get(&f.rounds[1].meas[0], 1) != IRREP_PAULI_LETTER_X) rc = 1;
    irrep_floquet_code_free(&f);
    return rc;
}

/* Build the 4-qubit honeycomb-style Floquet schedule:
 *   round 0 (red):    X_0 X_1, X_2 X_3
 *   round 1 (green):  Y_1 Y_2, Y_0 Y_3
 *   round 2 (blue):   Z_0 Z_2, Z_1 Z_3
 */
static void build_honeycomb_4q(irrep_floquet_code_t *f) {
    irrep_floquet_code_new(f, 4, 3);
    /* Round 0. */
    irrep_floquet_round_alloc(f, 0, 2);
    irrep_pauli_set(&f->rounds[0].meas[0], 0, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&f->rounds[0].meas[0], 1, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&f->rounds[0].meas[1], 2, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&f->rounds[0].meas[1], 3, IRREP_PAULI_LETTER_X);
    /* Round 1. */
    irrep_floquet_round_alloc(f, 1, 2);
    irrep_pauli_set(&f->rounds[1].meas[0], 1, IRREP_PAULI_LETTER_Y);
    irrep_pauli_set(&f->rounds[1].meas[0], 2, IRREP_PAULI_LETTER_Y);
    irrep_pauli_set(&f->rounds[1].meas[1], 0, IRREP_PAULI_LETTER_Y);
    irrep_pauli_set(&f->rounds[1].meas[1], 3, IRREP_PAULI_LETTER_Y);
    /* Round 2. */
    irrep_floquet_round_alloc(f, 2, 2);
    irrep_pauli_set(&f->rounds[2].meas[0], 0, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&f->rounds[2].meas[0], 2, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&f->rounds[2].meas[1], 1, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&f->rounds[2].meas[1], 3, IRREP_PAULI_LETTER_Z);
}

static int test_floquet_intra_round(void) {
    irrep_floquet_code_t f;
    build_honeycomb_4q(&f);
    int rc = irrep_floquet_code_check(&f) == IRREP_OK ? 0 : 1;
    irrep_floquet_code_free(&f);
    return rc;
}

/* Verify the Floquet-defining property: every measurement in round t
 * anti-commutes with every measurement in round (t+1) mod 3 that shares
 * a qubit. */
static int test_floquet_inter_round_anticomm(void) {
    irrep_floquet_code_t f;
    build_honeycomb_4q(&f);
    int rc = 0;
    /* For each pair (round 0, round 1), (round 1, round 2), (round 2, round 0):
     * count anti-commuting pairs and require each pair to anti-commute. */
    int pair_counts[3][2][2]; /* [round_pair][i][j] = symp inner */
    for (int rp = 0; rp < 3; ++rp) {
        int t1 = rp;
        int t2 = (rp + 1) % 3;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                pair_counts[rp][i][j] = irrep_pauli_symp_inner(
                    &f.rounds[t1].meas[i], &f.rounds[t2].meas[j]);
            }
        }
    }
    /* In our 4-qubit honeycomb-style schedule, every cross-round measurement
     * pair shares exactly one qubit and anti-commutes (XY, YZ, ZX). */
    for (int rp = 0; rp < 3; ++rp) {
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                if (pair_counts[rp][i][j] != 1) {
                    rc = 1;
                }
            }
        }
    }
    irrep_floquet_code_free(&f);
    return rc;
}

/* Materialise round 0 as a stabilizer group, check it passes
 * commutativity (the X⊗X measurements all commute pairwise). */
static int test_floquet_round_to_group(void) {
    irrep_floquet_code_t f;
    build_honeycomb_4q(&f);
    irrep_stabilizer_group_t g;
    int rc = 1;
    if (irrep_floquet_round_to_stabilizer_group(&f, 0, &g) == IRREP_OK) {
        if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK) rc = 0;
        irrep_stabilizer_group_free(&g);
    }
    irrep_floquet_code_free(&f);
    return rc;
}

/* Negative test: a malformed round with anti-commuting Paulis should be
 * rejected. */
static int test_floquet_bad_round_rejected(void) {
    irrep_floquet_code_t f;
    irrep_floquet_code_new(&f, 2, 1);
    irrep_floquet_round_alloc(&f, 0, 2);
    irrep_pauli_set(&f.rounds[0].meas[0], 0, IRREP_PAULI_LETTER_X); /* X_0 */
    irrep_pauli_set(&f.rounds[0].meas[1], 0, IRREP_PAULI_LETTER_Z); /* Z_0 */
    int rc = irrep_floquet_code_check(&f) == IRREP_ERR_PRECONDITION ? 0 : 1;
    irrep_floquet_code_free(&f);
    return rc;
}

int main(void) {
    int rc = 0;
    if (test_floquet_alloc())                  { fprintf(stderr, "FAIL test_floquet_alloc\n"); rc = 1; }
    if (test_floquet_intra_round())            { fprintf(stderr, "FAIL test_floquet_intra_round\n"); rc = 1; }
    if (test_floquet_inter_round_anticomm())   { fprintf(stderr, "FAIL test_floquet_inter_round_anticomm\n"); rc = 1; }
    if (test_floquet_round_to_group())         { fprintf(stderr, "FAIL test_floquet_round_to_group\n"); rc = 1; }
    if (test_floquet_bad_round_rejected())     { fprintf(stderr, "FAIL test_floquet_bad_round_rejected\n"); rc = 1; }
    return rc;
}
