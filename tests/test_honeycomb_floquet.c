/* SPDX-License-Identifier: MIT */
/* Tests for the Hastings-Haah honeycomb Floquet code (PRX Quantum 2021).
 *
 * Verifies:
 *  - Construction shape: 2·Lx·Ly qubits, 3 rounds, Lx·Ly measurements per round.
 *  - Each measurement has weight exactly 2 (a 2-body Pauli on a single edge).
 *  - Each vertex (A or B) appears in exactly one edge per round (each
 *    vertex has 3 incident edges, one of each color).
 *  - Intra-round commutativity: edges of one color form a perfect matching,
 *    so all measurements within a round commute (disjoint support).
 *  - Inter-round anti-commutation count: every B-vertex receives one edge
 *    in each round; the three operators X·X, Y·Y, Z·Z restricted to that
 *    B-vertex are X, Y, Z respectively, which anti-commute pairwise.
 *    Therefore the total count of anti-commuting (round t, round t+1)
 *    measurement pairs equals the number of vertices = 2·Lx·Ly.
 */
#include "harness.h"
#include <irrep/floquet_code.h>
#include <irrep/honeycomb_floquet.h>
#include <irrep/stabilizer_group.h>
#include <stdio.h>

static int test_honeycomb_shape(int Lx, int Ly) {
    irrep_floquet_code_t f;
    if (irrep_honeycomb_floquet_build(Lx, Ly, &f) != IRREP_OK) return 1;
    int rc = (f.n_qubits == 2 * Lx * Ly && f.n_rounds == 3) ? 0 : 1;
    for (int t = 0; t < 3; ++t) {
        if (f.rounds[t].n_meas != Lx * Ly) rc = 1;
        for (int i = 0; i < f.rounds[t].n_meas; ++i) {
            if (irrep_pauli_weight(&f.rounds[t].meas[i]) != 2) rc = 1;
        }
    }
    irrep_floquet_code_free(&f);
    return rc;
}

static int test_honeycomb_intra_round(int Lx, int Ly) {
    irrep_floquet_code_t f;
    if (irrep_honeycomb_floquet_build(Lx, Ly, &f) != IRREP_OK) return 1;
    irrep_status_t s = irrep_floquet_code_check(&f);
    irrep_floquet_code_free(&f);
    return s == IRREP_OK ? 0 : 1;
}

/* For each round, check every vertex appears in exactly one measurement
 * (each vertex has one edge of each color). */
static int test_honeycomb_perfect_matching(int Lx, int Ly) {
    irrep_floquet_code_t f;
    if (irrep_honeycomb_floquet_build(Lx, Ly, &f) != IRREP_OK) return 1;
    int rc = 0;
    for (int t = 0; t < 3 && rc == 0; ++t) {
        for (int q = 0; q < f.n_qubits && rc == 0; ++q) {
            int hits = 0;
            for (int i = 0; i < f.rounds[t].n_meas; ++i) {
                if (irrep_pauli_get(&f.rounds[t].meas[i], q) != IRREP_PAULI_LETTER_I) {
                    ++hits;
                }
            }
            if (hits != 1) rc = 1;
        }
    }
    irrep_floquet_code_free(&f);
    return rc;
}

/* Count anti-commuting pairs across consecutive rounds. Each B-vertex
 * receives exactly one edge of each color, so for each B-vertex there is
 * exactly one anti-commuting pair (one_R, one_G), one (one_G, one_B),
 * one (one_B, one_R). Same for A-vertices.
 * Total anti-commuting pairs across (round 0, round 1) = n_qubits.  */
static int test_honeycomb_inter_round_anticomm(int Lx, int Ly) {
    irrep_floquet_code_t f;
    if (irrep_honeycomb_floquet_build(Lx, Ly, &f) != IRREP_OK) return 1;
    int rc = 0;
    for (int rp = 0; rp < 3; ++rp) {
        int t1 = rp, t2 = (rp + 1) % 3;
        int anticomm = 0;
        for (int i = 0; i < f.rounds[t1].n_meas; ++i) {
            for (int j = 0; j < f.rounds[t2].n_meas; ++j) {
                if (irrep_pauli_symp_inner(&f.rounds[t1].meas[i],
                                            &f.rounds[t2].meas[j]) == 1) {
                    ++anticomm;
                }
            }
        }
        if (anticomm != f.n_qubits) rc = 1;
    }
    irrep_floquet_code_free(&f);
    return rc;
}

/* Round 0 materialised as a stabilizer group must pass commutativity. */
static int test_honeycomb_round0_group(int Lx, int Ly) {
    irrep_floquet_code_t f;
    if (irrep_honeycomb_floquet_build(Lx, Ly, &f) != IRREP_OK) return 1;
    irrep_stabilizer_group_t g;
    int rc = 1;
    if (irrep_floquet_round_to_stabilizer_group(&f, 0, &g) == IRREP_OK) {
        if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK) rc = 0;
        irrep_stabilizer_group_free(&g);
    }
    irrep_floquet_code_free(&f);
    return rc;
}

int main(void) {
    int rc = 0;
    int sizes[][2] = { {2, 2}, {3, 2}, {3, 3}, {4, 4} };
    for (size_t k = 0; k < sizeof(sizes)/sizeof(sizes[0]); ++k) {
        int Lx = sizes[k][0], Ly = sizes[k][1];
        if (test_honeycomb_shape(Lx, Ly)) {
            fprintf(stderr, "FAIL test_honeycomb_shape(%d, %d)\n", Lx, Ly); rc = 1;
        }
        if (test_honeycomb_intra_round(Lx, Ly)) {
            fprintf(stderr, "FAIL test_honeycomb_intra_round(%d, %d)\n", Lx, Ly); rc = 1;
        }
        if (test_honeycomb_perfect_matching(Lx, Ly)) {
            fprintf(stderr, "FAIL test_honeycomb_perfect_matching(%d, %d)\n", Lx, Ly); rc = 1;
        }
        if (test_honeycomb_inter_round_anticomm(Lx, Ly)) {
            fprintf(stderr, "FAIL test_honeycomb_inter_round_anticomm(%d, %d)\n", Lx, Ly); rc = 1;
        }
        if (test_honeycomb_round0_group(Lx, Ly)) {
            fprintf(stderr, "FAIL test_honeycomb_round0_group(%d, %d)\n", Lx, Ly); rc = 1;
        }
    }
    return rc;
}
