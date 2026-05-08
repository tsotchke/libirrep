/* SPDX-License-Identifier: MIT */
/* Tests for the CSS-code abstraction.
 *
 * Verifies:
 *  - Parity-matrix bit set/get correctness (single + multi-word).
 *  - Row-row F2 inner product reproduces textbook overlaps.
 *  - 5-qubit repetition code (n=5, m_X=4 phase, m_Z=4 bit, k=1)
 *    fails CSS verify because XX...X and Z...Z do anticommute (this
 *    is the Steane structure, not pure repetition — let's use Steane).
 *  - **Steane [[7,1,3]] code**: builds H_X = H_Z from the [7,4,3]
 *    Hamming dual, runs CSS-orthogonality verify, materialises as a
 *    stabilizer group, and confirms group commutativity.
 *  - Toric code on a 2x2 torus reconstructed from H_X (vertex incidence)
 *    and H_Z (plaquette incidence) — must pass CSS-verify because every
 *    vertex × plaquette pair shares 0 or 2 edges (even-parity overlap).
 *  - Materialises the toric CSS code as a stabilizer group and re-runs
 *    full pairwise commutativity (consistency check).
 */
#include "harness.h"
#include <irrep/css_code.h>
#include <irrep/stabilizer_group.h>
#include <irrep/toric_code.h>
#include <stdio.h>

static int test_parity_matrix_set_get(void) {
    irrep_parity_matrix_t m;
    if (irrep_parity_matrix_new(&m, 5, 130) != IRREP_OK) return 1;
    if (irrep_parity_matrix_set(&m, 0, 0)   != IRREP_OK) goto fail;
    if (irrep_parity_matrix_set(&m, 0, 63)  != IRREP_OK) goto fail;
    if (irrep_parity_matrix_set(&m, 0, 64)  != IRREP_OK) goto fail;
    if (irrep_parity_matrix_set(&m, 0, 129) != IRREP_OK) goto fail;
    if (irrep_parity_matrix_set(&m, 4, 100) != IRREP_OK) goto fail;
    if (irrep_parity_matrix_get(&m, 0, 0)   != 1) goto fail;
    if (irrep_parity_matrix_get(&m, 0, 63)  != 1) goto fail;
    if (irrep_parity_matrix_get(&m, 0, 64)  != 1) goto fail;
    if (irrep_parity_matrix_get(&m, 0, 129) != 1) goto fail;
    if (irrep_parity_matrix_get(&m, 0, 1)   != 0) goto fail;
    if (irrep_parity_matrix_get(&m, 4, 100) != 1) goto fail;
    if (irrep_parity_matrix_get(&m, 1, 100) != 0) goto fail;
    irrep_parity_matrix_free(&m);
    return 0;
fail:
    irrep_parity_matrix_free(&m);
    return 1;
}

static int test_parity_matrix_inner(void) {
    irrep_parity_matrix_t a, b;
    irrep_parity_matrix_new(&a, 2, 8);
    irrep_parity_matrix_new(&b, 2, 8);
    /* Row 0 of a: 11000000.  Row 0 of b: 11000000. → inner = 2 mod 2 = 0. */
    irrep_parity_matrix_set(&a, 0, 0);
    irrep_parity_matrix_set(&a, 0, 1);
    irrep_parity_matrix_set(&b, 0, 0);
    irrep_parity_matrix_set(&b, 0, 1);
    /* Row 1 of a: 11000000.  Row 1 of b: 10000000. → inner = 1 mod 2 = 1. */
    irrep_parity_matrix_set(&a, 1, 0);
    irrep_parity_matrix_set(&a, 1, 1);
    irrep_parity_matrix_set(&b, 1, 0);
    int rc = 0;
    if (irrep_parity_matrix_row_inner(&a, 0, &b, 0) != 0) rc = 1;
    if (irrep_parity_matrix_row_inner(&a, 1, &b, 1) != 1) rc = 1;
    irrep_parity_matrix_free(&a);
    irrep_parity_matrix_free(&b);
    return rc;
}

/* Steane [[7, 1, 3]] CSS code.
 *
 * H_X = H_Z = parity-check matrix of the [7, 4, 3] Hamming code:
 *
 *   1 0 1 0 1 0 1
 *   0 1 1 0 0 1 1
 *   0 0 0 1 1 1 1
 *
 * The Steane code has H_X = H_Z by construction; CSS orthogonality
 * H_X · H_Zᵀ = 0 holds because the Hamming code is self-dual under
 * the F2 inner product (each pair of parity rows overlaps in an even
 * number of columns: 2 for any two distinct rows).
 *
 * Verify:
 *   row 0 · row 1 = (cols 3, 7) → 2 → 0 mod 2 ✓
 *   row 0 · row 2 = (cols 5, 7) → 2 → 0 mod 2 ✓
 *   row 1 · row 2 = (cols 6, 7) → 2 → 0 mod 2 ✓
 */
static int build_steane(irrep_css_code_t *c) {
    if (irrep_css_code_new(c, 7, 3, 3) != IRREP_OK) return 1;
    /* H_X row 0: cols 0, 2, 4, 6 (1-indexed cols 1, 3, 5, 7). */
    irrep_parity_matrix_set(&c->H_X, 0, 0);
    irrep_parity_matrix_set(&c->H_X, 0, 2);
    irrep_parity_matrix_set(&c->H_X, 0, 4);
    irrep_parity_matrix_set(&c->H_X, 0, 6);
    /* H_X row 1: cols 1, 2, 5, 6. */
    irrep_parity_matrix_set(&c->H_X, 1, 1);
    irrep_parity_matrix_set(&c->H_X, 1, 2);
    irrep_parity_matrix_set(&c->H_X, 1, 5);
    irrep_parity_matrix_set(&c->H_X, 1, 6);
    /* H_X row 2: cols 3, 4, 5, 6. */
    irrep_parity_matrix_set(&c->H_X, 2, 3);
    irrep_parity_matrix_set(&c->H_X, 2, 4);
    irrep_parity_matrix_set(&c->H_X, 2, 5);
    irrep_parity_matrix_set(&c->H_X, 2, 6);
    /* H_Z = H_X for Steane. */
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 7; ++j) {
            if (irrep_parity_matrix_get(&c->H_X, i, j)) {
                irrep_parity_matrix_set(&c->H_Z, i, j);
            }
        }
    }
    return 0;
}

static int test_steane_css_orthogonal(void) {
    irrep_css_code_t c;
    if (build_steane(&c)) return 1;
    irrep_status_t s = irrep_css_code_verify(&c);
    irrep_css_code_free(&c);
    return s == IRREP_OK ? 0 : 1;
}

static int test_steane_to_stabilizer_group(void) {
    irrep_css_code_t c;
    if (build_steane(&c)) return 1;
    irrep_stabilizer_group_t g;
    if (irrep_css_code_to_stabilizer_group(&c, &g) != IRREP_OK) {
        irrep_css_code_free(&c);
        return 1;
    }
    irrep_status_t s = irrep_stabilizer_group_check_commutativity(&g);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&c);
    return s == IRREP_OK ? 0 : 1;
}

/* Build the 2x2 toric code as a CSS code and verify the orthogonality
 * condition (every vertex stabilizer × plaquette stabilizer share an even
 * number of edges). */
static int test_toric_as_css(void) {
    irrep_toric_params_t tp;
    if (irrep_toric_init(&tp, 2, 2) != IRREP_OK) return 1;
    int n_qubits = tp.n_qubits;             /* 8 */
    int m_v = tp.n_vertices;                /* 4 */
    int m_p = tp.n_plaquettes;              /* 4 */
    irrep_css_code_t c;
    if (irrep_css_code_new(&c, n_qubits, m_v, m_p) != IRREP_OK) return 1;
    /* H_X: vertex stabilizers (one row per vertex, 4 incident edges each). */
    int idx = 0;
    for (int vy = 0; vy < tp.Ly; ++vy) {
        for (int vx = 0; vx < tp.Lx; ++vx, ++idx) {
            int e[4];
            irrep_toric_vertex_edges(&tp, vx, vy, e);
            for (int k = 0; k < 4; ++k) {
                irrep_parity_matrix_set(&c.H_X, idx, e[k]);
            }
        }
    }
    idx = 0;
    for (int py = 0; py < tp.Ly; ++py) {
        for (int px = 0; px < tp.Lx; ++px, ++idx) {
            int e[4];
            irrep_toric_plaquette_edges(&tp, px, py, e);
            for (int k = 0; k < 4; ++k) {
                irrep_parity_matrix_set(&c.H_Z, idx, e[k]);
            }
        }
    }
    int rc = irrep_css_code_verify(&c) == IRREP_OK ? 0 : 1;
    /* Also: round-trip through stabilizer group and verify there. */
    if (rc == 0) {
        irrep_stabilizer_group_t g;
        if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
            if (irrep_stabilizer_group_check_commutativity(&g) != IRREP_OK) {
                rc = 1;
            }
            irrep_stabilizer_group_free(&g);
        } else {
            rc = 1;
        }
    }
    irrep_css_code_free(&c);
    return rc;
}

/* Negative test: a bad pair (X-stab on qubit 0, Z-stab on qubit 0) must
 * fail CSS-verify because they anti-commute. */
static int test_css_rejects_anticomm(void) {
    irrep_css_code_t c;
    if (irrep_css_code_new(&c, 4, 1, 1) != IRREP_OK) return 1;
    irrep_parity_matrix_set(&c.H_X, 0, 0); /* X on qubit 0 */
    irrep_parity_matrix_set(&c.H_Z, 0, 0); /* Z on qubit 0 */
    irrep_status_t s = irrep_css_code_verify(&c);
    irrep_css_code_free(&c);
    return s == IRREP_ERR_PRECONDITION ? 0 : 1;
}

int main(void) {
    int rc = 0;
    if (test_parity_matrix_set_get())     { fprintf(stderr, "FAIL test_parity_matrix_set_get\n"); rc = 1; }
    if (test_parity_matrix_inner())       { fprintf(stderr, "FAIL test_parity_matrix_inner\n"); rc = 1; }
    if (test_steane_css_orthogonal())     { fprintf(stderr, "FAIL test_steane_css_orthogonal\n"); rc = 1; }
    if (test_steane_to_stabilizer_group()){ fprintf(stderr, "FAIL test_steane_to_stabilizer_group\n"); rc = 1; }
    if (test_toric_as_css())              { fprintf(stderr, "FAIL test_toric_as_css\n"); rc = 1; }
    if (test_css_rejects_anticomm())      { fprintf(stderr, "FAIL test_css_rejects_anticomm\n"); rc = 1; }
    return rc;
}
