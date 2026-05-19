/* SPDX-License-Identifier: MIT */
/* Tests for `<irrep/fibonacci_mtc.h>` — the Fibonacci anyon MTC. */
#include "harness.h"
#include <irrep/fibonacci_mtc.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>

static int test_basics(void) {
    IRREP_TEST_START("fibonacci_mtc_basics");
    const double PHI = 1.6180339887498948482045868343656;
    IRREP_ASSERT_NEAR(irrep_fib_golden_ratio(), PHI, 1e-12);
    IRREP_ASSERT_NEAR(irrep_fib_quantum_dim(IRREP_FIB_OBJ_1),   1.0, 1e-12);
    IRREP_ASSERT_NEAR(irrep_fib_quantum_dim(IRREP_FIB_OBJ_TAU), PHI, 1e-12);
    /* D² = 1 + φ² = 2 + φ. */
    IRREP_ASSERT_NEAR(irrep_fib_global_dim() * irrep_fib_global_dim(),
                      2.0 + PHI, 1e-12);
    IRREP_ASSERT_NEAR(irrep_fib_central_charge(), 14.0 / 5.0, 1e-12);
    return IRREP_TEST_END();
}

static int test_fusion(void) {
    IRREP_TEST_START("fibonacci_mtc_fusion");
    IRREP_ASSERT(irrep_fib_fusion(IRREP_FIB_OBJ_1, IRREP_FIB_OBJ_1, IRREP_FIB_OBJ_1) == 1);
    IRREP_ASSERT(irrep_fib_fusion(IRREP_FIB_OBJ_1, IRREP_FIB_OBJ_TAU, IRREP_FIB_OBJ_TAU) == 1);
    IRREP_ASSERT(irrep_fib_fusion(IRREP_FIB_OBJ_TAU, IRREP_FIB_OBJ_1, IRREP_FIB_OBJ_TAU) == 1);
    /* τ × τ = 1 + τ */
    IRREP_ASSERT(irrep_fib_fusion(IRREP_FIB_OBJ_TAU, IRREP_FIB_OBJ_TAU, IRREP_FIB_OBJ_1)   == 1);
    IRREP_ASSERT(irrep_fib_fusion(IRREP_FIB_OBJ_TAU, IRREP_FIB_OBJ_TAU, IRREP_FIB_OBJ_TAU) == 1);
    /* Forbidden. */
    IRREP_ASSERT(irrep_fib_fusion(IRREP_FIB_OBJ_1, IRREP_FIB_OBJ_TAU, IRREP_FIB_OBJ_1) == 0);
    return IRREP_TEST_END();
}

static int test_S_matrix(void) {
    IRREP_TEST_START("fibonacci_mtc_S_matrix");
    const double PHI = 1.6180339887498948482045868343656;
    double D = irrep_fib_global_dim();
    /* S = (1/D) [[1, φ], [φ, -1]] */
    IRREP_ASSERT_NEAR(creal(irrep_fib_S_matrix(IRREP_FIB_OBJ_1,   IRREP_FIB_OBJ_1)),   1.0/D, 1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_fib_S_matrix(IRREP_FIB_OBJ_1,   IRREP_FIB_OBJ_TAU)), PHI/D, 1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_fib_S_matrix(IRREP_FIB_OBJ_TAU, IRREP_FIB_OBJ_1)),   PHI/D, 1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_fib_S_matrix(IRREP_FIB_OBJ_TAU, IRREP_FIB_OBJ_TAU)), -1.0/D, 1e-12);
    return IRREP_TEST_END();
}

static int test_T_eigenvalues(void) {
    IRREP_TEST_START("fibonacci_mtc_T_eigenvalues");
    /* T_τ = exp(4πi/5). */
    double _Complex T_tau = irrep_fib_T_eigenvalue(IRREP_FIB_OBJ_TAU);
    IRREP_ASSERT_NEAR(cabs(T_tau - cexp(I * 4.0 * M_PI / 5.0)), 0.0, 1e-12);
    IRREP_ASSERT_NEAR(cabs(irrep_fib_T_eigenvalue(IRREP_FIB_OBJ_1) - 1.0), 0.0, 1e-12);
    return IRREP_TEST_END();
}

static int test_F_symbol(void) {
    IRREP_TEST_START("fibonacci_mtc_F_symbol");
    const double PHI = 1.6180339887498948482045868343656;
    double inv_phi = 1.0 / PHI;
    double inv_sqrt_phi = 1.0 / sqrt(PHI);
    irrep_fib_object_t s = IRREP_FIB_OBJ_TAU;
    irrep_fib_object_t I_ = IRREP_FIB_OBJ_1;
    /* F^{τττ}_τ entries. */
    IRREP_ASSERT_NEAR(creal(irrep_fib_F_symbol(s, s, s, s, I_, I_)),  inv_phi,      1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_fib_F_symbol(s, s, s, s, I_, s )),  inv_sqrt_phi, 1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_fib_F_symbol(s, s, s, s, s,  I_)),  inv_sqrt_phi, 1e-12);
    IRREP_ASSERT_NEAR(creal(irrep_fib_F_symbol(s, s, s, s, s,  s )), -inv_phi,      1e-12);
    return IRREP_TEST_END();
}

static int test_R_symbols(void) {
    IRREP_TEST_START("fibonacci_mtc_R_symbols");
    irrep_fib_object_t s = IRREP_FIB_OBJ_TAU;
    irrep_fib_object_t I_ = IRREP_FIB_OBJ_1;
    /* R^{ττ}_1 = exp(-4πi/5). */
    double _Complex R1 = irrep_fib_R_symbol(s, s, I_);
    IRREP_ASSERT_NEAR(cabs(R1 - cexp(-I * 4.0 * M_PI / 5.0)), 0.0, 1e-12);
    /* R^{ττ}_τ = exp(3πi/5). */
    double _Complex Rt = irrep_fib_R_symbol(s, s, s);
    IRREP_ASSERT_NEAR(cabs(Rt - cexp(I * 3.0 * M_PI / 5.0)), 0.0, 1e-12);
    /* R^{1,a}_a = R^{a,1}_a = 1. */
    IRREP_ASSERT_NEAR(cabs(irrep_fib_R_symbol(I_, s, s) - 1.0), 0.0, 1e-12);
    return IRREP_TEST_END();
}

/* The Fibonacci MTC's full consistency suite, mirroring Ising and Z₂×Z₂. */
static int test_consistency_proofs(void) {
    IRREP_TEST_START("fibonacci_mtc_consistency");
    IRREP_ASSERT(irrep_fib_F_unitarity_residual() < 1e-12);
    IRREP_ASSERT(irrep_fib_S_squared_residual()    < 1e-12);
    IRREP_ASSERT(irrep_fib_verlinde_residual()     < 1e-12);
    IRREP_ASSERT(irrep_fib_twist_from_R_residual() < 1e-12);
    return IRREP_TEST_END();
}

int main(void) {
    int rc = 0;
    rc |= test_basics();
    rc |= test_fusion();
    rc |= test_S_matrix();
    rc |= test_T_eigenvalues();
    rc |= test_F_symbol();
    rc |= test_R_symbols();
    rc |= test_consistency_proofs();
    return rc;
}
