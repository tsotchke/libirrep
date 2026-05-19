/* SPDX-License-Identifier: MIT */
/* Tests for `<irrep/toric_mtc.h>` — the Z₂ × Z₂ toric-code modular
 * tensor category. Mirrors the layout of `test_crane_yetter_ising.c`
 * for Ising. */
#include "harness.h"
#include <irrep/toric_mtc.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>

static int test_basics(void) {
    IRREP_TEST_START("toric_mtc_basics");
    /* All quantum dimensions = 1. */
    for (int a = 0; a < IRREP_TORIC_MTC_N_OBJECTS; ++a) {
        IRREP_ASSERT_NEAR(irrep_toric_mtc_quantum_dim((irrep_toric_mtc_object_t)a),
                          1.0, 1e-12);
    }
    /* Global dim D = 2. */
    IRREP_ASSERT_NEAR(irrep_toric_mtc_global_dim(), 2.0, 1e-12);
    /* Central charge c = 0 (non-chiral). */
    IRREP_ASSERT_NEAR(irrep_toric_mtc_central_charge(), 0.0, 1e-12);
    return IRREP_TEST_END();
}

static int test_fusion_rules(void) {
    IRREP_TEST_START("toric_mtc_fusion_rules");
    irrep_toric_mtc_object_t I_ = IRREP_TORIC_MTC_OBJ_1;
    irrep_toric_mtc_object_t E_ = IRREP_TORIC_MTC_OBJ_E;
    irrep_toric_mtc_object_t M_ = IRREP_TORIC_MTC_OBJ_M;
    irrep_toric_mtc_object_t P_ = IRREP_TORIC_MTC_OBJ_PSI;

    /* Z₂ × Z₂ fusion. */
    IRREP_ASSERT(irrep_toric_mtc_fusion(E_, E_, I_) == 1);
    IRREP_ASSERT(irrep_toric_mtc_fusion(M_, M_, I_) == 1);
    IRREP_ASSERT(irrep_toric_mtc_fusion(P_, P_, I_) == 1);
    IRREP_ASSERT(irrep_toric_mtc_fusion(E_, M_, P_) == 1);
    IRREP_ASSERT(irrep_toric_mtc_fusion(E_, P_, M_) == 1);
    IRREP_ASSERT(irrep_toric_mtc_fusion(M_, P_, E_) == 1);
    /* Forbidden fusions. */
    IRREP_ASSERT(irrep_toric_mtc_fusion(E_, E_, E_) == 0);
    IRREP_ASSERT(irrep_toric_mtc_fusion(E_, E_, M_) == 0);
    IRREP_ASSERT(irrep_toric_mtc_fusion(E_, M_, I_) == 0);
    return IRREP_TEST_END();
}

static int test_S_matrix(void) {
    IRREP_TEST_START("toric_mtc_S_matrix");
    /* S = (1/2) [[1, 1, 1, 1], [1, 1, -1, -1],
     *            [1, -1, 1, -1], [1, -1, -1, 1]]
     * Order: 1, e, m, ψ. */
    const double expected[4][4] = {
        { 1,  1,  1,  1 },
        { 1,  1, -1, -1 },
        { 1, -1,  1, -1 },
        { 1, -1, -1,  1 },
    };
    for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
            double _Complex s = irrep_toric_mtc_S_matrix(
                (irrep_toric_mtc_object_t)a, (irrep_toric_mtc_object_t)b);
            double s_expected = 0.5 * expected[a][b];
            IRREP_ASSERT_NEAR(creal(s), s_expected, 1e-12);
            IRREP_ASSERT_NEAR(cimag(s), 0.0, 1e-12);
        }
    }
    return IRREP_TEST_END();
}

static int test_T_eigenvalues(void) {
    IRREP_TEST_START("toric_mtc_T_eigenvalues");
    IRREP_ASSERT_NEAR(creal(irrep_toric_mtc_T_eigenvalue(IRREP_TORIC_MTC_OBJ_1)),    1.0, 1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_toric_mtc_T_eigenvalue(IRREP_TORIC_MTC_OBJ_E)),    1.0, 1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_toric_mtc_T_eigenvalue(IRREP_TORIC_MTC_OBJ_M)),    1.0, 1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_toric_mtc_T_eigenvalue(IRREP_TORIC_MTC_OBJ_PSI)), -1.0, 1e-12);
    return IRREP_TEST_END();
}

static int test_R_symbols(void) {
    IRREP_TEST_START("toric_mtc_R_symbols");
    irrep_toric_mtc_object_t I_ = IRREP_TORIC_MTC_OBJ_1;
    irrep_toric_mtc_object_t E_ = IRREP_TORIC_MTC_OBJ_E;
    irrep_toric_mtc_object_t M_ = IRREP_TORIC_MTC_OBJ_M;
    irrep_toric_mtc_object_t P_ = IRREP_TORIC_MTC_OBJ_PSI;
    /* Self-braids. */
    IRREP_ASSERT_NEAR(creal(irrep_toric_mtc_R_symbol(E_, E_, I_)), 1.0, 1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_toric_mtc_R_symbol(M_, M_, I_)), 1.0, 1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_toric_mtc_R_symbol(P_, P_, I_)), -1.0, 1e-12);
    /* Mutual e ↔ m: R = i. */
    IRREP_ASSERT_NEAR(cimag(irrep_toric_mtc_R_symbol(E_, M_, P_)), 1.0, 1e-12);
    IRREP_ASSERT_NEAR(cimag(irrep_toric_mtc_R_symbol(M_, E_, P_)), 1.0, 1e-12);
    /* Mutual-statistics phase: R^{em} · R^{me} = (i)(i) = -1. */
    double _Complex prod = irrep_toric_mtc_R_symbol(E_, M_, P_)
                         * irrep_toric_mtc_R_symbol(M_, E_, P_);
    IRREP_ASSERT_NEAR(creal(prod), -1.0, 1e-12);
    IRREP_ASSERT_NEAR(cimag(prod),  0.0, 1e-12);
    /* Forbidden R. */
    IRREP_ASSERT_NEAR(cabs(irrep_toric_mtc_R_symbol(E_, E_, M_)), 0.0, 1e-12);
    return IRREP_TEST_END();
}

/* Consistency proofs: modular S² = I, Verlinde formula, and the
 * topological twist θ_a from R-symbols. */
static int test_consistency_proofs(void) {
    IRREP_TEST_START("toric_mtc_consistency");
    IRREP_ASSERT(irrep_toric_mtc_S_squared_residual()      < 1e-12);
    IRREP_ASSERT(irrep_toric_mtc_verlinde_residual()       < 1e-12);
    IRREP_ASSERT(irrep_toric_mtc_twist_from_R_residual()   < 1e-12);
    return IRREP_TEST_END();
}

/* Walker-Wang admissibility on the 3-simplex for Z₂×Z₂. Counts edge
 * labelings of the tetrahedron with all 4 vertex + 4 face constraints
 * satisfied.
 *
 * Cross-MTC observation: both Z₂×Z₂ and Ising give 16 on the 3-simplex
 * — even though one is abelian and one non-abelian, they share global
 * dimension D = 2, and the WW ground-state-dim on this small geometry
 * depends on D rather than the detailed anyon-model algebra. */
static int test_walker_wang_simplex3(void) {
    IRREP_TEST_START("toric_mtc_walker_wang_simplex3");
    long long c = irrep_toric_mtc_walker_wang_simplex3_full_count();
    IRREP_ASSERT(c == 16);
    return IRREP_TEST_END();
}

static int test_admissibility_basic(void) {
    IRREP_TEST_START("toric_mtc_admissibility");
    irrep_toric_mtc_object_t E_ = IRREP_TORIC_MTC_OBJ_E;
    irrep_toric_mtc_object_t M_ = IRREP_TORIC_MTC_OBJ_M;
    irrep_toric_mtc_object_t P_ = IRREP_TORIC_MTC_OBJ_PSI;
    irrep_toric_mtc_object_t I_ = IRREP_TORIC_MTC_OBJ_1;

    irrep_toric_mtc_object_t emp[3] = { E_, M_, P_ };
    IRREP_ASSERT(irrep_toric_mtc_admissible(emp, 3) == 1);   /* e × m × ψ = 1 ✓ */

    irrep_toric_mtc_object_t eee[3] = { E_, E_, E_ };
    IRREP_ASSERT(irrep_toric_mtc_admissible(eee, 3) == 0);   /* e × e × e = e ≠ 1 */

    irrep_toric_mtc_object_t eemm[4] = { E_, E_, M_, M_ };
    IRREP_ASSERT(irrep_toric_mtc_admissible(eemm, 4) == 1);  /* (e²)(m²) = 1 ✓ */

    irrep_toric_mtc_object_t all_one[3] = { I_, I_, I_ };
    IRREP_ASSERT(irrep_toric_mtc_admissible(all_one, 3) == 1);
    return IRREP_TEST_END();
}

int main(void) {
    int rc = 0;
    rc |= test_basics();
    rc |= test_fusion_rules();
    rc |= test_S_matrix();
    rc |= test_T_eigenvalues();
    rc |= test_R_symbols();
    rc |= test_consistency_proofs();
    rc |= test_admissibility_basic();
    rc |= test_walker_wang_simplex3();
    return rc;
}
