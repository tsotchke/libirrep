/* SPDX-License-Identifier: MIT */
/* Tests for the Panteleev-Kalachev lifted-product CSS qLDPC code over
 * the cyclic group ring F₂[ℤ/ℓ].
 *
 * Verifies:
 *  - Polynomial-matrix monomial set/get round-trip + XOR cancellation.
 *  - ℓ=1 reduces to the hypergraph product (CSS-orthogonal for any A, B).
 *  - ℓ=4 with 1×1 polynomial matrices A, B (a univariate-bicycle case)
 *    builds a CSS code with n = 2ℓ qubits and CSS-orthogonality holds.
 *  - ℓ=3 with 1×2 polynomial matrices (a small qLDPC instance) builds
 *    and CSS-orthogonality holds.
 *  - ℓ=4 with 2×2 polynomial matrices: a non-trivial lifted product
 *    where lift(A), lift(B) are 8×8 F₂ matrices; CSS verify must pass.
 *  - Materialised stabilizer group has full pairwise commutativity.
 */
#include "harness.h"
#include <irrep/css_code.h>
#include <irrep/lifted_product.h>
#include <irrep/stabilizer_group.h>
#include <stdio.h>

static int test_poly_matrix_set_get(void) {
    irrep_poly_matrix_t M;
    if (irrep_poly_matrix_new(&M, 3, 2, 5) != IRREP_OK) return 1;
    irrep_poly_matrix_add_monomial(&M, 1, 1, 3);
    irrep_poly_matrix_add_monomial(&M, 0, 0, 0);
    irrep_poly_matrix_add_monomial(&M, 0, 0, 4);
    int ok =
        irrep_poly_matrix_get(&M, 1, 1, 3) == 1 &&
        irrep_poly_matrix_get(&M, 0, 0, 0) == 1 &&
        irrep_poly_matrix_get(&M, 0, 0, 4) == 1 &&
        irrep_poly_matrix_get(&M, 0, 0, 1) == 0 &&
        irrep_poly_matrix_get(&M, 2, 1, 0) == 0;
    /* XOR cancel. */
    irrep_poly_matrix_add_monomial(&M, 1, 1, 3);
    if (irrep_poly_matrix_get(&M, 1, 1, 3) != 0) ok = 0;
    irrep_poly_matrix_free(&M);
    return ok ? 0 : 1;
}

/* ℓ=1 case: lifted product over F₂[1] = F₂. The result must be
 * CSS-orthogonal for any choice of A, B. */
static int test_ell1_reduces_to_hypergraph(void) {
    irrep_poly_matrix_t A, B;
    irrep_poly_matrix_new(&A, 2, 3, 1);
    irrep_poly_matrix_new(&B, 2, 3, 1);
    /* Set A = some 2x3 binary matrix, B = some 2x3 binary matrix. */
    irrep_poly_matrix_add_monomial(&A, 0, 0, 0);
    irrep_poly_matrix_add_monomial(&A, 0, 1, 0);
    irrep_poly_matrix_add_monomial(&A, 1, 1, 0);
    irrep_poly_matrix_add_monomial(&A, 1, 2, 0);
    irrep_poly_matrix_add_monomial(&B, 0, 0, 0);
    irrep_poly_matrix_add_monomial(&B, 0, 2, 0);
    irrep_poly_matrix_add_monomial(&B, 1, 1, 0);
    irrep_poly_matrix_add_monomial(&B, 1, 2, 0);

    irrep_css_code_t c;
    int rc = 1;
    if (irrep_lifted_product_build(&A, &B, &c) == IRREP_OK) {
        if (c.n == 1 * (3 * 3 + 2 * 2) && c.H_X.n_rows == 2 * 3 && c.H_Z.n_rows == 3 * 2) {
            if (irrep_css_code_verify(&c) == IRREP_OK) rc = 0;
        }
        irrep_css_code_free(&c);
    }
    irrep_poly_matrix_free(&A);
    irrep_poly_matrix_free(&B);
    return rc;
}

/* ℓ=4, 1×1 polynomial matrices A = 1 + x, B = 1 + x²
 * → n = 4·(1·1 + 1·1) = 8, m_X = m_Z = 4. */
static int test_ell4_1x1(void) {
    irrep_poly_matrix_t A, B;
    irrep_poly_matrix_new(&A, 1, 1, 4);
    irrep_poly_matrix_new(&B, 1, 1, 4);
    irrep_poly_matrix_add_monomial(&A, 0, 0, 0);
    irrep_poly_matrix_add_monomial(&A, 0, 0, 1);
    irrep_poly_matrix_add_monomial(&B, 0, 0, 0);
    irrep_poly_matrix_add_monomial(&B, 0, 0, 2);

    irrep_css_code_t c;
    int rc = 1;
    if (irrep_lifted_product_build(&A, &B, &c) == IRREP_OK) {
        if (c.n == 8 && c.H_X.n_rows == 4 && c.H_Z.n_rows == 4) {
            if (irrep_css_code_verify(&c) == IRREP_OK) {
                irrep_stabilizer_group_t g;
                if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
                    /* k = 2 for this LP instance (verified via F₂-rank). */
                    if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK
                        && irrep_css_code_logical_qubits(&c) == 2) rc = 0;
                    irrep_stabilizer_group_free(&g);
                }
            }
        }
        irrep_css_code_free(&c);
    }
    irrep_poly_matrix_free(&A);
    irrep_poly_matrix_free(&B);
    return rc;
}

/* ℓ=3, A = 1×2, B = 1×2 polynomial matrices
 * Lift sizes: lift(A) is 3×6, lift(B) is 3×6.
 * Lifted product: n = 3·(2·2 + 1·1) = 15, m_X = 3·1·2 = 6, m_Z = 3·2·1 = 6.
 * k = 3 (verified via F₂-rank). */
static int test_ell3_1x2(void) {
    irrep_poly_matrix_t A, B;
    irrep_poly_matrix_new(&A, 1, 2, 3);
    irrep_poly_matrix_new(&B, 1, 2, 3);
    irrep_poly_matrix_add_monomial(&A, 0, 0, 0);
    irrep_poly_matrix_add_monomial(&A, 0, 1, 1);
    irrep_poly_matrix_add_monomial(&B, 0, 0, 0);
    irrep_poly_matrix_add_monomial(&B, 0, 0, 2);
    irrep_poly_matrix_add_monomial(&B, 0, 1, 1);

    irrep_css_code_t c;
    int rc = 1;
    if (irrep_lifted_product_build(&A, &B, &c) == IRREP_OK) {
        if (c.n == 15 && c.H_X.n_rows == 6 && c.H_Z.n_rows == 6) {
            if (irrep_css_code_verify(&c) == IRREP_OK
                && irrep_css_code_logical_qubits(&c) == 3) rc = 0;
        }
        irrep_css_code_free(&c);
    }
    irrep_poly_matrix_free(&A);
    irrep_poly_matrix_free(&B);
    return rc;
}

/* ℓ=4, A = 2×2, B = 2×2: a non-trivial lifted product.
 *   n = 4·(2·2 + 2·2) = 32 qubits.
 *   m_X = 4·2·2 = 16, m_Z = 4·2·2 = 16. */
static int test_ell4_2x2(void) {
    irrep_poly_matrix_t A, B;
    irrep_poly_matrix_new(&A, 2, 2, 4);
    irrep_poly_matrix_new(&B, 2, 2, 4);
    /* A = [[1+x, x²],
     *      [x³, 1]]   */
    irrep_poly_matrix_add_monomial(&A, 0, 0, 0);
    irrep_poly_matrix_add_monomial(&A, 0, 0, 1);
    irrep_poly_matrix_add_monomial(&A, 0, 1, 2);
    irrep_poly_matrix_add_monomial(&A, 1, 0, 3);
    irrep_poly_matrix_add_monomial(&A, 1, 1, 0);
    /* B = [[x, 1+x³],
     *      [1, x²+x]]  */
    irrep_poly_matrix_add_monomial(&B, 0, 0, 1);
    irrep_poly_matrix_add_monomial(&B, 0, 1, 0);
    irrep_poly_matrix_add_monomial(&B, 0, 1, 3);
    irrep_poly_matrix_add_monomial(&B, 1, 0, 0);
    irrep_poly_matrix_add_monomial(&B, 1, 1, 2);
    irrep_poly_matrix_add_monomial(&B, 1, 1, 1);

    irrep_css_code_t c;
    int rc = 1;
    if (irrep_lifted_product_build(&A, &B, &c) == IRREP_OK) {
        if (c.n == 32 && c.H_X.n_rows == 16 && c.H_Z.n_rows == 16) {
            if (irrep_css_code_verify(&c) == IRREP_OK) {
                irrep_stabilizer_group_t g;
                if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
                    if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK) rc = 0;
                    irrep_stabilizer_group_free(&g);
                }
            }
        }
        irrep_css_code_free(&c);
    }
    irrep_poly_matrix_free(&A);
    irrep_poly_matrix_free(&B);
    return rc;
}

int main(void) {
    int rc = 0;
    if (test_poly_matrix_set_get())          { fprintf(stderr, "FAIL test_poly_matrix_set_get\n"); rc = 1; }
    if (test_ell1_reduces_to_hypergraph())   { fprintf(stderr, "FAIL test_ell1_reduces_to_hypergraph\n"); rc = 1; }
    if (test_ell4_1x1())                     { fprintf(stderr, "FAIL test_ell4_1x1\n"); rc = 1; }
    if (test_ell3_1x2())                     { fprintf(stderr, "FAIL test_ell3_1x2\n"); rc = 1; }
    if (test_ell4_2x2())                     { fprintf(stderr, "FAIL test_ell4_2x2\n"); rc = 1; }
    return rc;
}
