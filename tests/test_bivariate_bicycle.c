/* SPDX-License-Identifier: MIT */
/* Tests for bivariate-bicycle CSS qLDPC codes (Bravyi et al. Nature 2024).
 *
 * Verifies:
 *  - Polynomial monomial set/get round-trip in F₂[x,y]/(xˡ-1, yᵐ-1).
 *  - XOR semantics: adding the same monomial twice cancels in F₂.
 *  - Toy [[18, ?, ?]] BB code (ℓ=m=3, A=1+x+y, B=1+x²+y²) builds and
 *    passes CSS orthogonality verify.
 *  - Headline [[72, 12, 6]] BB code (ℓ=m=6, A=x³+y+y², B=y³+x+x²)
 *    builds, passes CSS orthogonality verify, materialises into a
 *    72-qubit / 72-generator stabilizer group, and the abstract
 *    pairwise-commutativity check passes.
 *  - Each row of H_X for [[72,12,6]] has weight 6 (LDPC sanity).
 */
#include "harness.h"
#include <irrep/bivariate_bicycle.h>
#include <irrep/css_code.h>
#include <irrep/stabilizer_group.h>
#include <stdio.h>

static int test_bb_poly_set_get(void) {
    irrep_bb_poly_t p;
    if (irrep_bb_poly_new(&p, 6, 6) != IRREP_OK) return 1;
    irrep_bb_poly_add_monomial(&p, 3, 0); /* x³ */
    irrep_bb_poly_add_monomial(&p, 0, 1); /* y */
    irrep_bb_poly_add_monomial(&p, 0, 2); /* y² */
    int ok =
        irrep_bb_poly_get(&p, 3, 0) == 1 &&
        irrep_bb_poly_get(&p, 0, 1) == 1 &&
        irrep_bb_poly_get(&p, 0, 2) == 1 &&
        irrep_bb_poly_get(&p, 0, 0) == 0 &&
        irrep_bb_poly_get(&p, 5, 5) == 0;
    /* XOR: add x³ again, should cancel. */
    irrep_bb_poly_add_monomial(&p, 3, 0);
    if (irrep_bb_poly_get(&p, 3, 0) != 0) ok = 0;
    irrep_bb_poly_free(&p);
    return ok ? 0 : 1;
}

/* Toy BB code with ℓ=m=3:
 *   A(x,y) = 1 + x + y
 *   B(x,y) = 1 + x² + y²
 * Verify CSS orthogonality (which holds for any A, B in commutative
 * F₂[x,y]/(xˡ-1, yᵐ-1) — the textbook BB construction). */
static int test_bb_toy_18(void) {
    irrep_bb_poly_t A, B;
    irrep_bb_poly_new(&A, 3, 3);
    irrep_bb_poly_new(&B, 3, 3);
    irrep_bb_poly_add_monomial(&A, 0, 0);
    irrep_bb_poly_add_monomial(&A, 1, 0);
    irrep_bb_poly_add_monomial(&A, 0, 1);
    irrep_bb_poly_add_monomial(&B, 0, 0);
    irrep_bb_poly_add_monomial(&B, 2, 0);
    irrep_bb_poly_add_monomial(&B, 0, 2);

    irrep_css_code_t c;
    int rc = 1;
    if (irrep_bb_code_build(&A, &B, &c) == IRREP_OK) {
        if (irrep_css_code_verify(&c) == IRREP_OK) {
            /* n_qubits = 2 * 3 * 3 = 18, m_X = m_Z = 9. */
            if (c.n == 18 && c.H_X.n_rows == 9 && c.H_Z.n_rows == 9) {
                rc = 0;
            }
        }
        irrep_css_code_free(&c);
    }
    irrep_bb_poly_free(&A);
    irrep_bb_poly_free(&B);
    return rc;
}

/* [[72, 12, 6]] BB code from Bravyi et al. Nature 2024:
 *   ℓ = 6, m = 6
 *   A(x, y) = x³ + y + y²
 *   B(x, y) = y³ + x + x²
 * (One of the smallest instances in the Nature paper.) */
static int test_bb_72_12_6(void) {
    irrep_bb_poly_t A, B;
    irrep_bb_poly_new(&A, 6, 6);
    irrep_bb_poly_new(&B, 6, 6);
    irrep_bb_poly_add_monomial(&A, 3, 0); /* x³ */
    irrep_bb_poly_add_monomial(&A, 0, 1); /* y */
    irrep_bb_poly_add_monomial(&A, 0, 2); /* y² */
    irrep_bb_poly_add_monomial(&B, 0, 3); /* y³ */
    irrep_bb_poly_add_monomial(&B, 1, 0); /* x */
    irrep_bb_poly_add_monomial(&B, 2, 0); /* x² */

    irrep_css_code_t c;
    int rc = 1;
    if (irrep_bb_code_build(&A, &B, &c) == IRREP_OK) {
        if (c.n == 72 && c.H_X.n_rows == 36 && c.H_Z.n_rows == 36) {
            if (irrep_css_code_verify(&c) == IRREP_OK) {
                /* Each row of H_X should have weight 6 (3 ones from
                 * A side + 3 from B side). Check row 0. */
                int weight_row0 = 0;
                for (int col = 0; col < 72; ++col) {
                    if (irrep_parity_matrix_get(&c.H_X, 0, col)) weight_row0++;
                }
                if (weight_row0 == 6) {
                    /* Materialise as stabilizer group + verify commutativity. */
                    irrep_stabilizer_group_t g;
                    if (irrep_css_code_to_stabilizer_group(&c, &g) == IRREP_OK) {
                        if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK) {
                            rc = 0;
                        }
                        irrep_stabilizer_group_free(&g);
                    }
                }
            }
        }
        irrep_css_code_free(&c);
    }
    irrep_bb_poly_free(&A);
    irrep_bb_poly_free(&B);
    return rc;
}

/* Named IBM instances: build, verify CSS orthogonality, check structural
 * counts. Distance verification is impractical via brute-force at these
 * sizes; we trust the polynomial-pair construction (Bravyi 2024 Table 3). */
static int test_bb_ibm_72_12_6(void) {
    irrep_css_code_t c;
    if (irrep_bb_code_ibm_72_12_6(&c) != IRREP_OK) return 1;
    int rc = 0;
    if (c.n != 72) rc = 1;
    if (c.H_X.n_rows != 36) rc = 1;
    if (c.H_Z.n_rows != 36) rc = 1;
    if (irrep_css_code_verify(&c) != IRREP_OK) rc = 1;
    if (irrep_css_code_logical_qubits(&c) != 12) {
        fprintf(stderr, "  [[72,12,6]] k = %d (expected 12)\n",
                irrep_css_code_logical_qubits(&c));
        rc = 1;
    }
    irrep_css_code_free(&c);
    return rc;
}

static int test_bb_ibm_144_12_12(void) {
    irrep_css_code_t c;
    if (irrep_bb_code_ibm_144_12_12(&c) != IRREP_OK) return 1;
    int rc = 0;
    if (c.n != 144) rc = 1;
    if (c.H_X.n_rows != 72) rc = 1;
    if (c.H_Z.n_rows != 72) rc = 1;
    if (irrep_css_code_verify(&c) != IRREP_OK) rc = 1;
    if (irrep_css_code_logical_qubits(&c) != 12) {
        fprintf(stderr, "  [[144,12,12]] k = %d (expected 12)\n",
                irrep_css_code_logical_qubits(&c));
        rc = 1;
    }
    irrep_css_code_free(&c);
    return rc;
}

static int test_bb_ibm_288_12_18(void) {
    irrep_css_code_t c;
    if (irrep_bb_code_ibm_288_12_18(&c) != IRREP_OK) return 1;
    int rc = 0;
    if (c.n != 288) rc = 1;
    if (c.H_X.n_rows != 144) rc = 1;
    if (c.H_Z.n_rows != 144) rc = 1;
    if (irrep_css_code_verify(&c) != IRREP_OK) rc = 1;
    if (irrep_css_code_logical_qubits(&c) != 12) {
        fprintf(stderr, "  [[288,12,18]] k = %d (expected 12)\n",
                irrep_css_code_logical_qubits(&c));
        rc = 1;
    }
    irrep_css_code_free(&c);
    return rc;
}

int main(void) {
    int rc = 0;
    if (test_bb_poly_set_get()) { fprintf(stderr, "FAIL test_bb_poly_set_get\n"); rc = 1; }
    if (test_bb_toy_18())       { fprintf(stderr, "FAIL test_bb_toy_18\n"); rc = 1; }
    if (test_bb_72_12_6())      { fprintf(stderr, "FAIL test_bb_72_12_6\n"); rc = 1; }
    if (test_bb_ibm_72_12_6())   { fprintf(stderr, "FAIL test_bb_ibm_72_12_6\n"); rc = 1; }
    if (test_bb_ibm_144_12_12()) { fprintf(stderr, "FAIL test_bb_ibm_144_12_12\n"); rc = 1; }
    if (test_bb_ibm_288_12_18()) { fprintf(stderr, "FAIL test_bb_ibm_288_12_18\n"); rc = 1; }
    return rc;
}
