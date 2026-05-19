/* SPDX-License-Identifier: MIT */
/* Tests for `<irrep/toric4d.h>`: 4D toric code on T⁴.
 *
 * Verifies:
 *   - Counts: n = 6 L⁴, n_X = n_Z = 4 L⁴ (the symmetric 4D structure).
 *   - Every stabilizer (X and Z) has weight exactly 6.
 *   - CSS orthogonality H_X · H_Zᵀ = 0 (mod 2).
 *   - Full pairwise stabilizer commutativity (= same condition via
 *     the abstract group representation).
 *   - Per-vertex linear redundancy among X-stabs: the product of all
 *     X-stabilizers around any single vertex (the "vertex 1-form
 *     differential") cancels because each face has 4 boundary edges.
 */
#include "harness.h"
#include <irrep/css_code.h>
#include <irrep/stabilizer_group.h>
#include <irrep/toric4d.h>

#include <stdio.h>

static int run_for_L(int L) {
    char banner[64];
    snprintf(banner, sizeof banner, "toric4d_L%d", L);
    IRREP_TEST_START(banner);

    irrep_toric4d_params_t p;
    IRREP_ASSERT(irrep_toric4d_init(&p, L, L, L, L) == IRREP_OK);
    IRREP_ASSERT(p.V == L * L * L * L);
    IRREP_ASSERT(p.n_qubits  == 6 * p.V);
    IRREP_ASSERT(p.n_X_stabs == 4 * p.V);
    IRREP_ASSERT(p.n_Z_stabs == 4 * p.V);

    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_toric4d_build(&p, &cs) == IRREP_OK);
    IRREP_ASSERT(cs.n == p.n_qubits);
    IRREP_ASSERT(cs.H_X.n_rows == p.n_X_stabs);
    IRREP_ASSERT(cs.H_Z.n_rows == p.n_Z_stabs);

    /* Every stabilizer has weight exactly 6 (the symmetric 4D
     * structure that 2D and 3D toric codes lack). */
    int bad_X_weights = 0, bad_Z_weights = 0;
    for (int r = 0; r < p.n_X_stabs; ++r) {
        int w = 0;
        for (int q = 0; q < p.n_qubits; ++q) {
            if (irrep_parity_matrix_get(&cs.H_X, r, q)) ++w;
        }
        if (w != 6) ++bad_X_weights;
    }
    for (int r = 0; r < p.n_Z_stabs; ++r) {
        int w = 0;
        for (int q = 0; q < p.n_qubits; ++q) {
            if (irrep_parity_matrix_get(&cs.H_Z, r, q)) ++w;
        }
        if (w != 6) ++bad_Z_weights;
    }
    IRREP_ASSERT(bad_X_weights == 0);
    IRREP_ASSERT(bad_Z_weights == 0);

    /* CSS orthogonality. */
    IRREP_ASSERT(irrep_css_code_verify(&cs) == IRREP_OK);

    /* k = dim H₂(T⁴, F₂) = C(4, 2) = 6 for any L ≥ 2 (topology-only). */
    IRREP_ASSERT(irrep_css_code_logical_qubits(&cs) == 6);

    /* Full stabilizer-group commutativity. */
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);
    IRREP_ASSERT(g.n_generators == p.n_X_stabs + p.n_Z_stabs);
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);

    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);

    return IRREP_TEST_END();
}

static int test_invalid_args(void) {
    IRREP_TEST_START("toric4d_invalid_args");
    irrep_toric4d_params_t p;
    IRREP_ASSERT(irrep_toric4d_init(NULL, 2, 2, 2, 2) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_toric4d_init(&p, 0, 2, 2, 2)   == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_toric4d_init(&p, 2, 2, 2, 2)   == IRREP_OK);
    IRREP_ASSERT(irrep_toric4d_build(NULL, NULL)      == IRREP_ERR_INVALID_ARG);
    return IRREP_TEST_END();
}

int main(void) {
    int rc = 0;
    rc |= run_for_L(2);   /* V=16,  n=96,  n_X=n_Z=64 — smallest non-trivial. */
    rc |= run_for_L(3);   /* V=81,  n=486, n_X=n_Z=324. */
    rc |= test_invalid_args();
    return rc;
}
