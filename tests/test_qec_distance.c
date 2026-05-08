/* SPDX-License-Identifier: MIT */
/* Tests for Pauli arithmetic + brute-force code distance.
 *
 * Verifies:
 *  - Pauli multiplication: textbook examples on small cases.
 *    X·Z = -iY corresponds to sign carry; Z·X has no carry.
 *    P·P = I (squared Pauli is identity, sign tracked).
 *  - Stabilizer-group membership: any single generator and any product
 *    of generators is in the group; the identity is in the group;
 *    a Pauli with disjoint support is NOT in the group.
 *  - **Distance verification on known small codes**:
 *      Steane [[7, 1, 3]]   → d = 3
 *      [[5, 1, 3]] (HaPPY)  → d = 3
 *      Surface d=2 [[4,1,2]] → d = 2
 *      Surface d=3 [[9,1,3]] → d = 3
 *      2x2 toric code        → d = 2
 */
#include "harness.h"
#include <irrep/color_code.h>
#include <irrep/css_code.h>
#include <irrep/happy_code.h>
#include <irrep/qec_distance.h>
#include <irrep/stabilizer_group.h>
#include <irrep/surface_code.h>
#include <irrep/toric3d.h>
#include <irrep/toric_code.h>
#include <stdio.h>
#include <string.h>

static int test_pauli_mul_xz(void) {
    /* X·Z on qubit 0: result = X·Z = -iY (in physical convention).
     * In our symplectic-only convention, support is x=z=1 (no Y mixing
     * at the bit level). The sign carry is z_P · x_Q = 0·1 = 0 for X·Z
     * (because X has z=0; the X part is x=1, z=0).
     * Wait: X·Z means P=X, Q=Z. P has x=1, z=0; Q has x=0, z=1.
     * carry = popcount(z_P AND x_Q) = popcount(0 AND 0) = 0.
     * Result: out.x = 1, out.z = 1, sign = +1 + +1 + 0 = +1 (POS).
     *
     * Z·X means P=Z, Q=X. P has x=0, z=1; Q has x=1, z=0.
     * carry = popcount(z_P AND x_Q) = popcount(1 AND 1) = 1.
     * Result: out.x = 1, out.z = 1, sign = +1 + +1 + 2 = -1 (NEG).
     *
     * So X·Z and Z·X differ by a sign — exactly the canonical anti-comm. */
    irrep_pauli_t X, Z, out;
    irrep_pauli_new(&X, 1);
    irrep_pauli_new(&Z, 1);
    irrep_pauli_new(&out, 1);
    irrep_pauli_set(&X, 0, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&Z, 0, IRREP_PAULI_LETTER_Z);
    int rc = 0;
    /* X · Z: */
    irrep_pauli_multiply(&out, &X, &Z);
    if (out.sign != IRREP_PAULI_SIGN_POS) rc = 1;
    if (irrep_pauli_get(&out, 0) != IRREP_PAULI_LETTER_Y) rc = 1; /* x=z=1 */
    /* Z · X: */
    irrep_pauli_multiply(&out, &Z, &X);
    if (out.sign != IRREP_PAULI_SIGN_NEG) rc = 1;
    if (irrep_pauli_get(&out, 0) != IRREP_PAULI_LETTER_Y) rc = 1;
    irrep_pauli_free(&X);
    irrep_pauli_free(&Z);
    irrep_pauli_free(&out);
    return rc;
}

static int test_pauli_mul_squared(void) {
    /* X·X = I, Z·Z = I, (XZ)(XZ) = ?
     * P = X·Z, Q = X·Z. P · Q has support x_P + x_Q = 0, z_P + z_Q = 0,
     * carry = z_P · x_Q = 1·1 = 1, so sign = +1 + +1 + 2 = NEG.
     * So (XZ)·(XZ) = -I in our convention (which is correct: (XZ)² = -1
     * since XZ·XZ = X·(ZX)·Z = X·(-XZ)·Z = -X²·Z² = -I).
     */
    irrep_pauli_t XZ, out;
    irrep_pauli_new(&XZ, 1);
    irrep_pauli_new(&out, 1);
    irrep_pauli_set(&XZ, 0, IRREP_PAULI_LETTER_Y); /* x=z=1 in our convention */
    irrep_pauli_multiply(&out, &XZ, &XZ);
    int rc = 0;
    if (out.sign != IRREP_PAULI_SIGN_NEG) rc = 1;
    if (out.x[0] != 0 || out.z[0] != 0) rc = 1; /* identity support */
    irrep_pauli_free(&XZ);
    irrep_pauli_free(&out);
    return rc;
}

static int test_group_contains_generator(void) {
    /* Steane: any generator is in the group. */
    irrep_css_code_t cs;
    if (irrep_color_steane(&cs) != IRREP_OK) return 1;
    irrep_stabilizer_group_t g;
    if (irrep_css_code_to_stabilizer_group(&cs, &g) != IRREP_OK) {
        irrep_css_code_free(&cs);
        return 1;
    }
    int rc = 0;
    for (int i = 0; i < g.n_generators; ++i) {
        if (irrep_stabilizer_group_contains(&g, &g.gens[i]) != 1) rc = 1;
    }
    /* Identity should also be in the group. */
    irrep_pauli_t id;
    irrep_pauli_new(&id, 7);
    if (irrep_stabilizer_group_contains(&g, &id) != 1) rc = 1;
    /* X_0 alone (weight-1 X) should NOT be in the group (not a Steane stab). */
    irrep_pauli_set(&id, 0, IRREP_PAULI_LETTER_X);
    if (irrep_stabilizer_group_contains(&g, &id) != 0) rc = 1;
    irrep_pauli_free(&id);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return rc;
}

static int test_group_contains_product(void) {
    /* Product of two generators is in the group. */
    irrep_css_code_t cs;
    irrep_color_steane(&cs);
    irrep_stabilizer_group_t g;
    irrep_css_code_to_stabilizer_group(&cs, &g);
    irrep_pauli_t prod;
    irrep_pauli_new(&prod, 7);
    irrep_pauli_multiply(&prod, &g.gens[0], &g.gens[1]);
    int rc = irrep_stabilizer_group_contains(&g, &prod) == 1 ? 0 : 1;
    irrep_pauli_free(&prod);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return rc;
}

/* Distance of Steane [[7, 1, 3]] = 3. */
static int test_distance_steane(void) {
    irrep_css_code_t cs;
    irrep_color_steane(&cs);
    irrep_stabilizer_group_t g;
    irrep_css_code_to_stabilizer_group(&cs, &g);
    int d = irrep_qec_distance_brute(&g, 4);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return d == 3 ? 0 : 1;
}

/* Distance of [[5, 1, 3]] (HaPPY perfect tensor) = 3. */
static int test_distance_5qubit(void) {
    irrep_stabilizer_group_t g;
    irrep_happy_perfect_tensor_5_1_3(&g);
    int d = irrep_qec_distance_brute(&g, 4);
    irrep_stabilizer_group_free(&g);
    return d == 3 ? 0 : 1;
}

/* Distance of [[d², 1, d]] surface code at d=2 should be 2,
 * at d=3 should be 3. */
static int test_distance_surface(int d_param, int expected) {
    irrep_surface_params_t p;
    irrep_surface_init(&p, d_param);
    irrep_css_code_t cs;
    irrep_surface_build(&p, &cs);
    irrep_stabilizer_group_t g;
    irrep_css_code_to_stabilizer_group(&cs, &g);
    int d = irrep_qec_distance_brute(&g, expected + 1);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return d == expected ? 0 : 1;
}

/* 2x2 toric code: 8 qubits, d=2. */
static int test_distance_toric_2x2(void) {
    irrep_toric_params_t tp;
    irrep_toric_init(&tp, 2, 2);
    irrep_stabilizer_group_t g;
    if (irrep_stabilizer_group_new(&g, tp.n_qubits, 8) != IRREP_OK) return 1;
    int idx = 0;
    for (int vy = 0; vy < tp.Ly; ++vy) {
        for (int vx = 0; vx < tp.Lx; ++vx, ++idx) {
            int e[4];
            irrep_toric_vertex_edges(&tp, vx, vy, e);
            for (int k = 0; k < 4; ++k)
                irrep_pauli_set(&g.gens[idx], e[k], IRREP_PAULI_LETTER_X);
        }
    }
    for (int py = 0; py < tp.Ly; ++py) {
        for (int px = 0; px < tp.Lx; ++px, ++idx) {
            int e[4];
            irrep_toric_plaquette_edges(&tp, px, py, e);
            for (int k = 0; k < 4; ++k)
                irrep_pauli_set(&g.gens[idx], e[k], IRREP_PAULI_LETTER_Z);
        }
    }
    int d = irrep_qec_distance_brute(&g, 3);
    irrep_stabilizer_group_free(&g);
    return d == 2 ? 0 : 1;
}

/* 3D toric L=2: 24 qubits, expected distance = min(Lx, Ly, Lz) = 2.
 * Logical operators are Wilson lines along non-contractible cycles. */
static int test_distance_toric3d_222(void) {
    irrep_toric3d_params_t p;
    irrep_toric3d_init(&p, 2, 2, 2);
    irrep_css_code_t cs;
    if (irrep_toric3d_build(&p, &cs) != IRREP_OK) return 1;
    irrep_stabilizer_group_t g;
    irrep_css_code_to_stabilizer_group(&cs, &g);
    int d = irrep_qec_distance_brute(&g, 3);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return d == 2 ? 0 : 1;
}

/* Surface d=4 [[16, 1, 4]]: brute force at weight up to 4 must find d=4. */
static int test_distance_surface_4(void) {
    irrep_surface_params_t p;
    irrep_surface_init(&p, 4);
    irrep_css_code_t cs;
    irrep_surface_build(&p, &cs);
    irrep_stabilizer_group_t g;
    irrep_css_code_to_stabilizer_group(&cs, &g);
    int d = irrep_qec_distance_brute(&g, 5);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return d == 4 ? 0 : 1;
}

int main(void) {
    int rc = 0;
    if (test_pauli_mul_xz())                  { fprintf(stderr, "FAIL test_pauli_mul_xz\n"); rc = 1; }
    if (test_pauli_mul_squared())             { fprintf(stderr, "FAIL test_pauli_mul_squared\n"); rc = 1; }
    if (test_group_contains_generator())      { fprintf(stderr, "FAIL test_group_contains_generator\n"); rc = 1; }
    if (test_group_contains_product())        { fprintf(stderr, "FAIL test_group_contains_product\n"); rc = 1; }
    if (test_distance_steane())               { fprintf(stderr, "FAIL test_distance_steane\n"); rc = 1; }
    if (test_distance_5qubit())               { fprintf(stderr, "FAIL test_distance_5qubit\n"); rc = 1; }
    if (test_distance_surface(2, 2))          { fprintf(stderr, "FAIL test_distance_surface(d=2)\n"); rc = 1; }
    if (test_distance_surface(3, 3))          { fprintf(stderr, "FAIL test_distance_surface(d=3)\n"); rc = 1; }
    if (test_distance_surface_4())            { fprintf(stderr, "FAIL test_distance_surface_4\n"); rc = 1; }
    if (test_distance_toric_2x2())            { fprintf(stderr, "FAIL test_distance_toric_2x2\n"); rc = 1; }
    if (test_distance_toric3d_222())          { fprintf(stderr, "FAIL test_distance_toric3d_222\n"); rc = 1; }
    return rc;
}
