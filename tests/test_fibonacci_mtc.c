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

/* Walker-Wang admissibility-restricted ground-state-dim counts on the
 * five polyhedra (simplex, cube, octahedron, bipyramid, prism) for the
 * Fibonacci MTC. Distinct from Ising and Z₂×Z₂ counts — proves the
 * MTC dependence visible across all three modules.
 *
 * Polyhedral duality (cube=octahedron, bipyramid=prism) holds for
 * Fibonacci just as for Ising and Z₂×Z₂. */
static int test_walker_wang_polyhedra(void) {
    IRREP_TEST_START("fibonacci_mtc_walker_wang_polyhedra");
    long long c_simplex    = irrep_fib_walker_wang_simplex3_full_count();
    long long c_cube       = irrep_fib_walker_wang_cube_full_count();
    long long c_octahedron = irrep_fib_walker_wang_octahedron_full_count();
    long long c_bipyramid  = irrep_fib_walker_wang_tri_bipyramid_full_count();
    long long c_prism      = irrep_fib_walker_wang_tri_prism_full_count();

    IRREP_ASSERT(c_simplex    == 11);
    IRREP_ASSERT(c_cube       == 145);
    IRREP_ASSERT(c_octahedron == 145);
    IRREP_ASSERT(c_bipyramid  == 33);
    IRREP_ASSERT(c_prism      == 33);

    /* Polyhedral duality. */
    IRREP_ASSERT(c_cube == c_octahedron);
    IRREP_ASSERT(c_bipyramid == c_prism);
    return IRREP_TEST_END();
}

/* Admissibility primitive: only #τ == 1 is non-admissible. */
static int test_admissibility(void) {
    IRREP_TEST_START("fibonacci_mtc_admissibility");
    irrep_fib_object_t I_ = IRREP_FIB_OBJ_1;
    irrep_fib_object_t s = IRREP_FIB_OBJ_TAU;
    /* #τ = 0 (all 1's): admissible. */
    {
        irrep_fib_object_t v[3] = { I_, I_, I_ };
        IRREP_ASSERT(irrep_fib_admissible(v, 3) == 1);
    }
    /* #τ = 1: NOT admissible. */
    {
        irrep_fib_object_t v[3] = { s, I_, I_ };
        IRREP_ASSERT(irrep_fib_admissible(v, 3) == 0);
    }
    /* #τ = 2: admissible (τ² = 1 + τ). */
    {
        irrep_fib_object_t v[3] = { s, s, I_ };
        IRREP_ASSERT(irrep_fib_admissible(v, 3) == 1);
    }
    /* #τ = 3: admissible (τ³ = 1 + 2τ). */
    {
        irrep_fib_object_t v[3] = { s, s, s };
        IRREP_ASSERT(irrep_fib_admissible(v, 3) == 1);
    }
    /* #τ = 4: admissible (τ⁴ = 2 + 3τ). */
    {
        irrep_fib_object_t v[4] = { s, s, s, s };
        IRREP_ASSERT(irrep_fib_admissible(v, 4) == 1);
    }
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
    rc |= test_admissibility();
    rc |= test_walker_wang_polyhedra();

    /* Anyonic braiding: σ_1, σ_2 are unitaries and satisfy Yang-Baxter.
     * These two facts together define a faithful representation of the
     * braid group B_3 on the 2-dim fusion space of 4 Fibonacci anyons
     * — the algebraic content of "universal topological quantum
     * computation via Fibonacci braiding". */
    IRREP_TEST_START("fibonacci_mtc_braiding");
    IRREP_ASSERT(irrep_fib_braid_unitarity_residual()   < 1e-12);
    IRREP_ASSERT(irrep_fib_braid_yang_baxter_residual() < 1e-12);
    rc |= IRREP_TEST_END();
    return rc;
}
