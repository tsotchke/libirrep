/* SPDX-License-Identifier: MIT */
/* Tests for CSS Floquet code on square lattice (Davydova et al. 2023).
 *
 * Verifies:
 *  - Construction shape: 2·Lx·Ly qubits, 4 rounds, Lx·Ly weight-2 measurements
 *    per round.
 *  - Intra-round commutativity: each round measures on disjoint pairs.
 *  - Round 0 (XX horizontal) and round 1 (ZZ horizontal) anti-commute on
 *    the same edge: X·Z anti-commutes on each shared qubit, but the edge
 *    has 2 shared qubits → 2 anti-commutations → COMMUTE on a single edge!
 *
 *    Wait — that's wrong. XX on edge {a, b} vs ZZ on same edge {a, b}:
 *    overlap is qubits a and b, both with X·Z anti-commutation.
 *    Symplectic inner = 2 mod 2 = 0, so XX and ZZ on the SAME edge COMMUTE.
 *
 *  Actually the Floquet dynamics comes from anti-commutation BETWEEN
 *  edges sharing a single vertex. So we test:
 *    - Within a round: all measurements pairwise commute (disjoint edges).
 *    - Between consecutive rounds with DIFFERENT ORIENTATION (e.g. round 1
 *      ZZ horizontal vs round 2 XX vertical): pairs sharing one qubit
 *      anti-commute (X·Z anti-commute on the shared qubit).
 *    - Total inter-round anti-commutation between round 1 → round 2
 *      = number of shared (qubit, edge_orientation) pairs = some count
 *      we can compute exactly.
 *  - Each round's measurements form a valid stabilizer subgroup (commute
 *    pairwise) when materialised.
 */
#include "harness.h"
#include <irrep/css_floquet.h>
#include <irrep/floquet_code.h>
#include <irrep/stabilizer_group.h>
#include <stdio.h>

static int test_shape(int Lx, int Ly) {
    irrep_floquet_code_t f;
    if (irrep_css_floquet_square_build(Lx, Ly, &f) != IRREP_OK) return 1;
    int rc = (f.n_qubits == 2*Lx*Ly && f.n_rounds == 4) ? 0 : 1;
    for (int r = 0; r < 4; ++r) {
        if (f.rounds[r].n_meas != Lx*Ly) rc = 1;
        for (int i = 0; i < f.rounds[r].n_meas; ++i) {
            if (irrep_pauli_weight(&f.rounds[r].meas[i]) != 2) rc = 1;
        }
    }
    irrep_floquet_code_free(&f);
    return rc;
}

static int test_intra_round_commute(int Lx, int Ly) {
    irrep_floquet_code_t f;
    if (irrep_css_floquet_square_build(Lx, Ly, &f) != IRREP_OK) return 1;
    int rc = irrep_floquet_code_check(&f) == IRREP_OK ? 0 : 1;
    irrep_floquet_code_free(&f);
    return rc;
}

/* Round 0 (XX on horizontal edges) and Round 1 (ZZ on the same horizontal
 * edges) should pairwise COMMUTE: on a single edge {a,b}, X_a X_b and Z_a Z_b
 * share both qubits, contributing 2 anti-commutation events → 0 mod 2 → commute.
 *
 * For different edges the support is disjoint, so commute trivially. */
static int test_round0_round1_commute(int Lx, int Ly) {
    irrep_floquet_code_t f;
    if (irrep_css_floquet_square_build(Lx, Ly, &f) != IRREP_OK) return 1;
    int rc = 0;
    for (int i = 0; i < f.rounds[0].n_meas; ++i) {
        for (int j = 0; j < f.rounds[1].n_meas; ++j) {
            if (irrep_pauli_symp_inner(&f.rounds[0].meas[i],
                                        &f.rounds[1].meas[j]) != 0) {
                rc = 1;
            }
        }
    }
    irrep_floquet_code_free(&f);
    return rc;
}

/* Round 1 (ZZ horizontal edges) vs Round 2 (XX vertical edges):
 * a horizontal edge ZZ on {A(x,y), B(x,y)} and a vertical edge XX on
 * {B(x,y), A(x,y+1)} share exactly the qubit B(x,y) → X·Z anti-commute on
 * one shared qubit → anti-commute (symp = 1).
 *
 * Total count of (i, j) pairs with anti-commutation = number of (h-edge,
 * v-edge) pairs sharing a B-vertex = n_qubits/2 = Lx·Ly (one per B-vertex,
 * since each B is shared by exactly one h-edge and one v-edge).
 *
 * Plus the symmetric case: h-edge ZZ ∈ {A(x,y), B(x,y)} and v-edge XX ∈
 * {B(x, y-1), A(x, y)} share qubit A(x,y) → also anti-commute.
 *
 * So total anti-commute pairs = 2·Lx·Ly = n_qubits. */
static int test_round1_round2_anticomm_count(int Lx, int Ly) {
    irrep_floquet_code_t f;
    if (irrep_css_floquet_square_build(Lx, Ly, &f) != IRREP_OK) return 1;
    int anti = 0;
    for (int i = 0; i < f.rounds[1].n_meas; ++i) {
        for (int j = 0; j < f.rounds[2].n_meas; ++j) {
            if (irrep_pauli_symp_inner(&f.rounds[1].meas[i],
                                        &f.rounds[2].meas[j]) == 1) {
                ++anti;
            }
        }
    }
    irrep_floquet_code_free(&f);
    return anti == 2 * Lx * Ly ? 0 : 1;
}

static int test_round_to_stabilizer_group(int Lx, int Ly) {
    irrep_floquet_code_t f;
    if (irrep_css_floquet_square_build(Lx, Ly, &f) != IRREP_OK) return 1;
    int rc = 0;
    for (int r = 0; r < 4 && rc == 0; ++r) {
        irrep_stabilizer_group_t g;
        if (irrep_floquet_round_to_stabilizer_group(&f, r, &g) == IRREP_OK) {
            if (irrep_stabilizer_group_check_commutativity(&g) != IRREP_OK) rc = 1;
            irrep_stabilizer_group_free(&g);
        } else {
            rc = 1;
        }
    }
    irrep_floquet_code_free(&f);
    return rc;
}

int main(void) {
    int rc = 0;
    int sizes[][2] = { {2, 2}, {3, 2}, {3, 3}, {4, 4} };
    for (size_t k = 0; k < sizeof(sizes)/sizeof(sizes[0]); ++k) {
        int Lx = sizes[k][0], Ly = sizes[k][1];
        if (test_shape(Lx, Ly))
            { fprintf(stderr, "FAIL test_shape(%d,%d)\n", Lx, Ly); rc = 1; }
        if (test_intra_round_commute(Lx, Ly))
            { fprintf(stderr, "FAIL test_intra_round_commute(%d,%d)\n", Lx, Ly); rc = 1; }
        if (test_round0_round1_commute(Lx, Ly))
            { fprintf(stderr, "FAIL test_round0_round1_commute(%d,%d)\n", Lx, Ly); rc = 1; }
        if (test_round1_round2_anticomm_count(Lx, Ly))
            { fprintf(stderr, "FAIL test_round1_round2_anticomm_count(%d,%d)\n", Lx, Ly); rc = 1; }
        if (test_round_to_stabilizer_group(Lx, Ly))
            { fprintf(stderr, "FAIL test_round_to_stabilizer_group(%d,%d)\n", Lx, Ly); rc = 1; }
    }
    return rc;
}
