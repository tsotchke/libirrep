/* SPDX-License-Identifier: MIT */
/* Tests for concatenated CSS codes.
 *
 * Verifies:
 *  - Steane ⊗ Steane = [[49, 1, ≥9]] CSS code:
 *      n_qubits = 49, m_X = 7·3 + 3 = 24, m_Z = 7·3 + 3 = 24, k = 1.
 *      CSS-orthogonality holds at all 3 levels (inner-inner, inner-outer
 *      via X̄/Z̄ centraliser, outer-outer via outer CSS).
 *  - Pure-X / pure-Z precondition rejection.
 */
#include "harness.h"
#include <irrep/color_code.h>
#include <irrep/concatenated_code.h>
#include <irrep/css_code.h>
#include <irrep/stabilizer_group.h>
#include <stdio.h>

static int test_steane_steane_concat(void) {
    irrep_css_code_t inner, outer;
    if (irrep_color_steane(&inner) != IRREP_OK) return 1;
    if (irrep_color_steane(&outer) != IRREP_OK) {
        irrep_css_code_free(&inner);
        return 1;
    }
    /* Steane logical X̄ = X⊗7, Z̄ = Z⊗7. */
    irrep_pauli_t Xbar, Zbar;
    irrep_pauli_new(&Xbar, 7);
    irrep_pauli_new(&Zbar, 7);
    for (int q = 0; q < 7; ++q) {
        irrep_pauli_set(&Xbar, q, IRREP_PAULI_LETTER_X);
        irrep_pauli_set(&Zbar, q, IRREP_PAULI_LETTER_Z);
    }

    irrep_css_code_t cat;
    int rc = 1;
    if (irrep_css_concatenate(&inner, &Xbar, &Zbar, &outer, &cat) == IRREP_OK) {
        if (cat.n == 49 && cat.H_X.n_rows == 24 && cat.H_Z.n_rows == 24) {
            if (irrep_css_code_verify(&cat) == IRREP_OK) {
                /* Stabilizer-group materialisation passes commutativity. */
                irrep_stabilizer_group_t g;
                if (irrep_css_code_to_stabilizer_group(&cat, &g) == IRREP_OK) {
                    if (irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK) {
                        rc = 0;
                    }
                    irrep_stabilizer_group_free(&g);
                }
            }
        }
        irrep_css_code_free(&cat);
    }
    irrep_pauli_free(&Xbar);
    irrep_pauli_free(&Zbar);
    irrep_css_code_free(&inner);
    irrep_css_code_free(&outer);
    return rc;
}

/* Reject if X̄ has a Z bit set or Z̄ has an X bit set. */
static int test_concat_rejects_mixed_logicals(void) {
    irrep_css_code_t inner, outer;
    irrep_color_steane(&inner);
    irrep_color_steane(&outer);
    irrep_pauli_t Xbar, Zbar;
    irrep_pauli_new(&Xbar, 7);
    irrep_pauli_new(&Zbar, 7);
    /* Make Xbar mixed (set Y on qubit 0). */
    irrep_pauli_set(&Xbar, 0, IRREP_PAULI_LETTER_Y);
    for (int q = 0; q < 7; ++q) irrep_pauli_set(&Zbar, q, IRREP_PAULI_LETTER_Z);

    irrep_css_code_t cat;
    irrep_status_t s = irrep_css_concatenate(&inner, &Xbar, &Zbar, &outer, &cat);
    int rc = (s == IRREP_ERR_PRECONDITION) ? 0 : 1;
    if (s == IRREP_OK) irrep_css_code_free(&cat);
    irrep_pauli_free(&Xbar);
    irrep_pauli_free(&Zbar);
    irrep_css_code_free(&inner);
    irrep_css_code_free(&outer);
    return rc;
}

int main(void) {
    int rc = 0;
    if (test_steane_steane_concat())      { fprintf(stderr, "FAIL test_steane_steane_concat\n"); rc = 1; }
    if (test_concat_rejects_mixed_logicals()) { fprintf(stderr, "FAIL test_concat_rejects_mixed_logicals\n"); rc = 1; }
    return rc;
}
