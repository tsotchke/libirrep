/* SPDX-License-Identifier: MIT */
/* Tests for the Tillich-Zemor hypergraph-product CSS qLDPC code.
 *
 * Verifies:
 *  - Trivial case: H_a = H_b = (1) the 1x1 zero-rank check yields a code on
 *    1·1 + 1·1 = 2 qubits with 1 X-stab + 1 Z-stab and CSS-orthogonality.
 *
 *  - Repetition × repetition: classical [3, 1, 3] repetition code parity-
 *    check H_a = H_b = [[1, 1, 0], [0, 1, 1]] (m=2, n=3), giving a quantum
 *    [[3·3 + 2·2, ?, ?]] = [[13, ?, ?]] hypergraph-product code. CSS
 *    orthogonality must hold. (This is the planar surface-code analogue
 *    in the boundary-not-PBC version.)
 *
 *  - Larger repetition×repetition: [4, 1, 4] × [4, 1, 4] giving
 *    [[4·4 + 3·3, ?, ?]] = [[25, ?, ?]] code. CSS-orthogonality and
 *    materialised stabilizer-group commutativity must hold.
 *
 *  - Negative orthogonality test (impossible by construction so just sanity).
 */
#include "harness.h"
#include <irrep/css_code.h>
#include <irrep/hypergraph_product.h>
#include <irrep/qec_distance.h>
#include <irrep/stabilizer_group.h>
#include <stdio.h>

static int build_repetition_check(irrep_parity_matrix_t *H, int n) {
    /* Repetition code parity-check H[i, j] = δ_{j, i} XOR δ_{j, i+1}, m = n-1. */
    if (n < 2) return 1;
    int m = n - 1;
    if (irrep_parity_matrix_new(H, m, n) != IRREP_OK) return 1;
    for (int i = 0; i < m; ++i) {
        irrep_parity_matrix_set(H, i, i);
        irrep_parity_matrix_set(H, i, i + 1);
    }
    return 0;
}

static int test_hgp_trivial(void) {
    /* H_a = H_b = (1) — single 1×1 check. */
    irrep_parity_matrix_t H_a, H_b;
    irrep_parity_matrix_new(&H_a, 1, 1);
    irrep_parity_matrix_new(&H_b, 1, 1);
    irrep_parity_matrix_set(&H_a, 0, 0);
    irrep_parity_matrix_set(&H_b, 0, 0);
    irrep_css_code_t c;
    int rc = 1;
    if (irrep_hypergraph_product_build(&H_a, &H_b, &c) == IRREP_OK) {
        /* n_q = 1·1 + 1·1 = 2, m_X = 1·1 = 1, m_Z = 1·1 = 1. */
        if (c.n == 2 && c.H_X.n_rows == 1 && c.H_Z.n_rows == 1) {
            if (irrep_css_code_verify(&c) == IRREP_OK) {
                rc = 0;
            }
        }
        irrep_css_code_free(&c);
    }
    irrep_parity_matrix_free(&H_a);
    irrep_parity_matrix_free(&H_b);
    return rc;
}

static int test_hgp_repetition_3x3(void) {
    /* [3, 1, 3] × [3, 1, 3]: n_a = n_b = 3, m_a = m_b = 2.
     *   n_q = 9 + 4 = 13, m_X = 2·3 = 6, m_Z = 3·2 = 6. */
    irrep_parity_matrix_t H_a, H_b;
    if (build_repetition_check(&H_a, 3)) return 1;
    if (build_repetition_check(&H_b, 3)) {
        irrep_parity_matrix_free(&H_a);
        return 1;
    }
    irrep_css_code_t c;
    int rc = 1;
    if (irrep_hypergraph_product_build(&H_a, &H_b, &c) == IRREP_OK) {
        if (c.n == 13 && c.H_X.n_rows == 6 && c.H_Z.n_rows == 6) {
            if (irrep_css_code_verify(&c) == IRREP_OK) rc = 0;
        }
        irrep_css_code_free(&c);
    }
    irrep_parity_matrix_free(&H_a);
    irrep_parity_matrix_free(&H_b);
    return rc;
}

static int test_hgp_repetition_4x4(void) {
    /* [4, 1, 4] × [4, 1, 4]: n_a = n_b = 4, m_a = m_b = 3.
     *   n_q = 16 + 9 = 25, m_X = 3·4 = 12, m_Z = 4·3 = 12. */
    irrep_parity_matrix_t H_a, H_b;
    if (build_repetition_check(&H_a, 4)) return 1;
    if (build_repetition_check(&H_b, 4)) {
        irrep_parity_matrix_free(&H_a);
        return 1;
    }
    irrep_css_code_t c;
    int rc = 1;
    if (irrep_hypergraph_product_build(&H_a, &H_b, &c) == IRREP_OK) {
        if (c.n == 25 && c.H_X.n_rows == 12 && c.H_Z.n_rows == 12) {
            if (irrep_css_code_verify(&c) == IRREP_OK) {
                /* Round-trip through the abstract stabilizer group. */
                irrep_stabilizer_group_t g;
                if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
                    if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK) {
                        rc = 0;
                    }
                    irrep_stabilizer_group_free(&g);
                }
            }
        }
        irrep_css_code_free(&c);
    }
    irrep_parity_matrix_free(&H_a);
    irrep_parity_matrix_free(&H_b);
    return rc;
}

/* Cross-check: hand-verify CSS-orthogonality for [3,1,3] x [3,1,3]:
 * For each X-stab row (6 of them) and each Z-stab row (6 of them), the F2
 * inner product must equal 0. We've already done this via irrep_css_code_verify,
 * but this re-verifies with explicit row counts. */
static int test_hgp_explicit_orthogonality(void) {
    irrep_parity_matrix_t H_a, H_b;
    if (build_repetition_check(&H_a, 3)) return 1;
    if (build_repetition_check(&H_b, 3)) {
        irrep_parity_matrix_free(&H_a);
        return 1;
    }
    irrep_css_code_t c;
    int rc = 1;
    if (irrep_hypergraph_product_build(&H_a, &H_b, &c) == IRREP_OK) {
        rc = 0;
        for (int i = 0; i < c.H_X.n_rows; ++i) {
            for (int k = 0; k < c.H_Z.n_rows; ++k) {
                int inner = 0;
                for (int j = 0; j < c.n; ++j) {
                    inner ^= (irrep_parity_matrix_get(&c.H_X, i, j) &
                              irrep_parity_matrix_get(&c.H_Z, k, j));
                }
                if (inner != 0) { rc = 1; goto done; }
            }
        }
done:
        irrep_css_code_free(&c);
    }
    irrep_parity_matrix_free(&H_a);
    irrep_parity_matrix_free(&H_b);
    return rc;
}

/* Named instance: [[13, 1, 3]] HGP of repetition-3 × repetition-3.
 *
 * Verifies: structural counts, CSS orthogonality, k = 1 via the
 * F₂-rank logical-qubit primitive, AND d = 3 via brute-force
 * enumeration of weight ≤ 3 Paulis. */
static int test_hgp_named_13_1_3(void) {
    irrep_css_code_t cs;
    if (irrep_hgp_repetition_3_13_1_3(&cs) != IRREP_OK) return 1;
    int rc = 0;
    if (cs.n != 13) rc = 1;
    if (cs.H_X.n_rows != 6) rc = 1;
    if (cs.H_Z.n_rows != 6) rc = 1;
    if (irrep_css_code_verify(&cs) != IRREP_OK) rc = 1;
    if (irrep_css_code_logical_qubits(&cs) != 1) {
        fprintf(stderr, "  [[13,1,3]] HGP k = %d (expected 1)\n",
                irrep_css_code_logical_qubits(&cs));
        rc = 1;
    }
    /* Distance verification via the brute-force enumerator. */
    irrep_stabilizer_group_t g;
    if (irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK) {
        int d = irrep_qec_distance_brute(&g, /*max_weight*/ 3);
        if (d != 3) {
            fprintf(stderr, "  [[13,1,3]] HGP d = %d (expected 3)\n", d);
            rc = 1;
        }
        irrep_stabilizer_group_free(&g);
    } else {
        rc = 1;
    }
    irrep_css_code_free(&cs);
    return rc;
}

int main(void) {
    int rc = 0;
    if (test_hgp_trivial())               { fprintf(stderr, "FAIL test_hgp_trivial\n"); rc = 1; }
    if (test_hgp_repetition_3x3())        { fprintf(stderr, "FAIL test_hgp_repetition_3x3\n"); rc = 1; }
    if (test_hgp_repetition_4x4())        { fprintf(stderr, "FAIL test_hgp_repetition_4x4\n"); rc = 1; }
    if (test_hgp_explicit_orthogonality()){ fprintf(stderr, "FAIL test_hgp_explicit_orthogonality\n"); rc = 1; }
    if (test_hgp_named_13_1_3())          { fprintf(stderr, "FAIL test_hgp_named_13_1_3\n"); rc = 1; }

    /* Named [[25, 1, 4]] HGP(rep-4, rep-4): structural counts +
     * CSS-orth + k = 1 + brute-force distance = 4. */
    {
        irrep_css_code_t cs;
        if (irrep_hgp_repetition_4_25_1_4(&cs) != IRREP_OK) {
            fprintf(stderr, "FAIL HGP[[25,1,4]] build\n"); rc = 1;
        } else {
            int sub = 0;
            if (cs.n != 25) sub = 1;
            if (cs.H_X.n_rows != 12) sub = 1;
            if (cs.H_Z.n_rows != 12) sub = 1;
            if (irrep_css_code_verify(&cs) != IRREP_OK) sub = 1;
            if (irrep_css_code_logical_qubits(&cs) != 1) sub = 1;
            irrep_stabilizer_group_t g;
            if (irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK) {
                int d = irrep_qec_distance_brute(&g, /*max_weight*/ 4);
                if (d != 4) {
                    fprintf(stderr, "  [[25,1,4]] HGP d = %d (expected 4)\n", d);
                    sub = 1;
                }
                irrep_stabilizer_group_free(&g);
            } else { sub = 1; }
            if (sub) { fprintf(stderr, "FAIL HGP[[25,1,4]]\n"); rc = 1; }
            irrep_css_code_free(&cs);
        }
    }
    return rc;
}
