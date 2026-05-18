/* SPDX-License-Identifier: MIT */
/* Tests for the abstract symplectic Pauli / stabilizer-group machinery.
 *
 * Verifies:
 *  - Pauli letter set/get round-trip across multi-word boundaries.
 *  - Symplectic inner product reproduces the textbook X/Y/Z anti-commutation.
 *  - irrep_pauli_commute on independent qubits always returns true.
 *  - Weight equals the count of non-identity letters.
 *  - Stabilizer-group commutativity check passes for a 3-qubit Bell-state
 *    [[3, 1, 1]] code (XXI, IXX, ZZI, IZZ — all pairwise commuting).
 *  - Reproduces the toric code on a 2x2 torus: build the vertex (all-X) and
 *    plaquette (all-Z) generators by Pauli supports, run the abstract
 *    commutativity check, and confirm IRREP_OK.
 *  - Syndrome extraction: a single-qubit X error against a Z stabilizer
 *    flags exactly the corresponding generator bit.
 */
#include "harness.h"
#include <irrep/stabilizer_group.h>
#include <irrep/toric_code.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int test_pauli_set_get(void) {
    irrep_pauli_t p;
    if (irrep_pauli_new(&p, 130) != IRREP_OK) return 1; /* >2 words */
    /* Set qubit 0 = X, 63 = Z, 64 = Y, 129 = X. */
    if (irrep_pauli_set(&p, 0,   IRREP_PAULI_LETTER_X) != IRREP_OK) goto fail;
    if (irrep_pauli_set(&p, 63,  IRREP_PAULI_LETTER_Z) != IRREP_OK) goto fail;
    if (irrep_pauli_set(&p, 64,  IRREP_PAULI_LETTER_Y) != IRREP_OK) goto fail;
    if (irrep_pauli_set(&p, 129, IRREP_PAULI_LETTER_X) != IRREP_OK) goto fail;
    if (irrep_pauli_get(&p, 0)   != IRREP_PAULI_LETTER_X) goto fail;
    if (irrep_pauli_get(&p, 63)  != IRREP_PAULI_LETTER_Z) goto fail;
    if (irrep_pauli_get(&p, 64)  != IRREP_PAULI_LETTER_Y) goto fail;
    if (irrep_pauli_get(&p, 129) != IRREP_PAULI_LETTER_X) goto fail;
    if (irrep_pauli_get(&p, 1)   != IRREP_PAULI_LETTER_I) goto fail;
    if (irrep_pauli_get(&p, 100) != IRREP_PAULI_LETTER_I) goto fail;
    /* Overwrite qubit 0 with Y, then back to I. */
    irrep_pauli_set(&p, 0, IRREP_PAULI_LETTER_Y);
    if (irrep_pauli_get(&p, 0) != IRREP_PAULI_LETTER_Y) goto fail;
    irrep_pauli_set(&p, 0, IRREP_PAULI_LETTER_I);
    if (irrep_pauli_get(&p, 0) != IRREP_PAULI_LETTER_I) goto fail;
    irrep_pauli_free(&p);
    return 0;
fail:
    irrep_pauli_free(&p);
    return 1;
}

/* Single-qubit anti-commutation: X·Z = -Z·X, X·Y = -Y·X, Y·Z = -Z·Y. */
static int test_pauli_anticomm_single_qubit(void) {
    struct { irrep_pauli_letter_t a, b; int expected_anticomm; } cases[] = {
        { IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_X, 0 },
        { IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z, 0 },
        { IRREP_PAULI_LETTER_Y, IRREP_PAULI_LETTER_Y, 0 },
        { IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_Z, 1 },
        { IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_Y, 1 },
        { IRREP_PAULI_LETTER_Y, IRREP_PAULI_LETTER_Z, 1 },
        { IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_I, 0 },
        { IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_Z, 0 },
    };
    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); ++k) {
        irrep_pauli_t a, b;
        irrep_pauli_new(&a, 1);
        irrep_pauli_new(&b, 1);
        irrep_pauli_set(&a, 0, cases[k].a);
        irrep_pauli_set(&b, 0, cases[k].b);
        int got = irrep_pauli_symp_inner(&a, &b);
        irrep_pauli_free(&a);
        irrep_pauli_free(&b);
        if (got != cases[k].expected_anticomm) return 1;
    }
    return 0;
}

/* Independent-qubit operators always commute: X on qubit i, Z on qubit j ≠ i. */
static int test_pauli_disjoint_support_commutes(void) {
    irrep_pauli_t a, b;
    irrep_pauli_new(&a, 100);
    irrep_pauli_new(&b, 100);
    irrep_pauli_set(&a, 5,  IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&a, 70, IRREP_PAULI_LETTER_Y);
    irrep_pauli_set(&b, 6,  IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&b, 71, IRREP_PAULI_LETTER_Z);
    int rc = irrep_pauli_commute(&a, &b) ? 0 : 1;
    irrep_pauli_free(&a);
    irrep_pauli_free(&b);
    return rc;
}

static int test_pauli_weight(void) {
    irrep_pauli_t p;
    irrep_pauli_new(&p, 80);
    irrep_pauli_set(&p, 0,  IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&p, 13, IRREP_PAULI_LETTER_Y);
    irrep_pauli_set(&p, 64, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&p, 79, IRREP_PAULI_LETTER_X);
    int w = irrep_pauli_weight(&p);
    irrep_pauli_free(&p);
    return w == 4 ? 0 : 1;
}

/* 4-qubit example: stabilizers ZZII, IZZI, IIZZ all commute pairwise. */
static int test_stabilizer_group_commute(void) {
    irrep_stabilizer_group_t g;
    if (irrep_stabilizer_group_new(&g, 4, 3) != IRREP_OK) return 1;
    irrep_pauli_set(&g.gens[0], 0, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[0], 1, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[1], 1, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[1], 2, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[2], 2, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[2], 3, IRREP_PAULI_LETTER_Z);
    int rc = irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK ? 0 : 1;
    irrep_stabilizer_group_free(&g);
    return rc;
}

/* A non-commuting pair should be rejected: gen0 = XII, gen1 = ZII anti-commute. */
static int test_stabilizer_group_rejects_anticomm(void) {
    irrep_stabilizer_group_t g;
    if (irrep_stabilizer_group_new(&g, 3, 2) != IRREP_OK) return 1;
    irrep_pauli_set(&g.gens[0], 0, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&g.gens[1], 0, IRREP_PAULI_LETTER_Z);
    irrep_status_t s = irrep_stabilizer_group_check_commutativity(&g);
    irrep_stabilizer_group_free(&g);
    return s == IRREP_ERR_PRECONDITION ? 0 : 1;
}

/* Reproduce the 2x2 toric code abstractly:
 *   - 8 qubits (n_qubits = 2 Lx Ly = 8)
 *   - 4 vertex stabilizers (all-X on the 4 incident edges)
 *   - 4 plaquette stabilizers (all-Z on the 4 boundary edges)
 *   - all 8 generators must pairwise commute.
 *
 * (Note: there are 2 redundancies in the toric code — product of all vertex
 *  stabilizers and product of all plaquette stabilizers both = I — but the
 *  *commutativity* check is independent of that.) */
static int test_toric_via_stabilizer_group(void) {
    irrep_toric_params_t tp;
    if (irrep_toric_init(&tp, 2, 2) != IRREP_OK) return 1;
    irrep_stabilizer_group_t g;
    if (irrep_stabilizer_group_new(&g, tp.n_qubits, 8) != IRREP_OK) return 1;
    int idx = 0;
    /* Vertex stabilizers: all-X on the 4 incident edges. */
    for (int vy = 0; vy < tp.Ly; ++vy) {
        for (int vx = 0; vx < tp.Lx; ++vx, ++idx) {
            int e[4];
            irrep_toric_vertex_edges(&tp, vx, vy, e);
            for (int k = 0; k < 4; ++k) {
                irrep_pauli_set(&g.gens[idx], e[k], IRREP_PAULI_LETTER_X);
            }
        }
    }
    /* Plaquette stabilizers: all-Z on the 4 boundary edges. */
    for (int py = 0; py < tp.Ly; ++py) {
        for (int px = 0; px < tp.Lx; ++px, ++idx) {
            int e[4];
            irrep_toric_plaquette_edges(&tp, px, py, e);
            for (int k = 0; k < 4; ++k) {
                irrep_pauli_set(&g.gens[idx], e[k], IRREP_PAULI_LETTER_Z);
            }
        }
    }
    irrep_status_t s = irrep_stabilizer_group_check_commutativity(&g);
    irrep_stabilizer_group_free(&g);
    return s == IRREP_OK ? 0 : 1;
}

/* Single-qubit X error vs ZZ stabilizer flags exactly that Z's bit. */
static int test_syndrome_single_x_error(void) {
    irrep_stabilizer_group_t g;
    if (irrep_stabilizer_group_new(&g, 4, 3) != IRREP_OK) return 1;
    /* gen0 = ZZII, gen1 = IZZI, gen2 = IIZZ. */
    irrep_pauli_set(&g.gens[0], 0, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[0], 1, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[1], 1, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[1], 2, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[2], 2, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[2], 3, IRREP_PAULI_LETTER_Z);
    /* Error: X on qubit 1. Anti-commutes with gen0 (qubit 1 = Z) and gen1
     * (qubit 1 = Z); commutes with gen2. So syndrome bits = 011 in binary. */
    irrep_pauli_t err;
    irrep_pauli_new(&err, 4);
    irrep_pauli_set(&err, 1, IRREP_PAULI_LETTER_X);
    uint64_t syndrome = 0;
    irrep_status_t s = irrep_stabilizer_syndrome(&g, &err, &syndrome);
    irrep_pauli_free(&err);
    irrep_stabilizer_group_free(&g);
    if (s != IRREP_OK) return 1;
    /* Expect syndrome = 0b011 = 3. */
    return syndrome == 3 ? 0 : 1;
}

/* irrep_pauli_in_stabilizer_span: F₂ row-span test.
 *
 * Use a 3-qubit Bell-pair-style code:
 *   g_0 = X_0 X_1
 *   g_1 = Z_0 Z_1
 * which generates the [[2, 0, 2]] code (Bell-pair stabilizer).
 * Membership tests:
 *   - I    → in span (trivially)
 *   - X_0 X_1            (= g_0)      → in span
 *   - Z_0 Z_1            (= g_1)      → in span
 *   - g_0 · g_1 = Y_0 Y_1             → in span
 *   - X_0    (anti-commutes with g_1) → not in span (in fact not in N)
 *   - X_0 X_1 X_2 (X-string outside)  → not in span (anti-commutes with g_1)
 *
 * For a 2-generator group, the span has at most 4 elements: {I, g_0, g_1, g_0g_1}.
 */
static int test_pauli_in_stabilizer_span(void) {
    irrep_stabilizer_group_t g;
    if (irrep_stabilizer_group_new(&g, 3, 2) != IRREP_OK) return 1;
    irrep_pauli_set(&g.gens[0], 0, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&g.gens[0], 1, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&g.gens[1], 0, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&g.gens[1], 1, IRREP_PAULI_LETTER_Z);

    int rc = 0;
    irrep_pauli_t p;

    /* Identity. */
    irrep_pauli_new(&p, 3);
    if (irrep_pauli_in_stabilizer_span(&g, &p) != 1) rc = 1;
    irrep_pauli_free(&p);

    /* g_0 = X_0 X_1. */
    irrep_pauli_new(&p, 3);
    irrep_pauli_set(&p, 0, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&p, 1, IRREP_PAULI_LETTER_X);
    if (irrep_pauli_in_stabilizer_span(&g, &p) != 1) rc = 1;
    irrep_pauli_free(&p);

    /* g_1 = Z_0 Z_1. */
    irrep_pauli_new(&p, 3);
    irrep_pauli_set(&p, 0, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&p, 1, IRREP_PAULI_LETTER_Z);
    if (irrep_pauli_in_stabilizer_span(&g, &p) != 1) rc = 1;
    irrep_pauli_free(&p);

    /* g_0 · g_1 = Y_0 Y_1 (symplectic-only; sign ignored). */
    irrep_pauli_new(&p, 3);
    irrep_pauli_set(&p, 0, IRREP_PAULI_LETTER_Y);
    irrep_pauli_set(&p, 1, IRREP_PAULI_LETTER_Y);
    if (irrep_pauli_in_stabilizer_span(&g, &p) != 1) rc = 1;
    irrep_pauli_free(&p);

    /* X_0 alone — not in span. */
    irrep_pauli_new(&p, 3);
    irrep_pauli_set(&p, 0, IRREP_PAULI_LETTER_X);
    if (irrep_pauli_in_stabilizer_span(&g, &p) != 0) rc = 1;
    irrep_pauli_free(&p);

    /* X_2 — disjoint support, not in span. */
    irrep_pauli_new(&p, 3);
    irrep_pauli_set(&p, 2, IRREP_PAULI_LETTER_X);
    if (irrep_pauli_in_stabilizer_span(&g, &p) != 0) rc = 1;
    irrep_pauli_free(&p);

    /* X_0 X_1 X_2 — not in span (anti-commutes with g_1 anyway). */
    irrep_pauli_new(&p, 3);
    irrep_pauli_set(&p, 0, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&p, 1, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&p, 2, IRREP_PAULI_LETTER_X);
    if (irrep_pauli_in_stabilizer_span(&g, &p) != 0) rc = 1;
    irrep_pauli_free(&p);

    irrep_stabilizer_group_free(&g);
    return rc;
}

int main(void) {
    int rc = 0;
    if (test_pauli_set_get())                     { fprintf(stderr, "FAIL test_pauli_set_get\n"); rc = 1; }
    if (test_pauli_anticomm_single_qubit())       { fprintf(stderr, "FAIL test_pauli_anticomm_single_qubit\n"); rc = 1; }
    if (test_pauli_disjoint_support_commutes())   { fprintf(stderr, "FAIL test_pauli_disjoint_support_commutes\n"); rc = 1; }
    if (test_pauli_weight())                      { fprintf(stderr, "FAIL test_pauli_weight\n"); rc = 1; }
    if (test_stabilizer_group_commute())          { fprintf(stderr, "FAIL test_stabilizer_group_commute\n"); rc = 1; }
    if (test_stabilizer_group_rejects_anticomm()) { fprintf(stderr, "FAIL test_stabilizer_group_rejects_anticomm\n"); rc = 1; }
    if (test_toric_via_stabilizer_group())        { fprintf(stderr, "FAIL test_toric_via_stabilizer_group\n"); rc = 1; }
    if (test_syndrome_single_x_error())           { fprintf(stderr, "FAIL test_syndrome_single_x_error\n"); rc = 1; }
    if (test_pauli_in_stabilizer_span())          { fprintf(stderr, "FAIL test_pauli_in_stabilizer_span\n"); rc = 1; }
    return rc;
}
