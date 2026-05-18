/* SPDX-License-Identifier: MIT */
/* Tests for the Crane-Yetter / Walker-Wang Ising-MTC runtime.
 *
 * Verifies:
 *  - Ising quantum dimensions: d_1 = 1, d_σ = √2, d_ψ = 1.
 *  - Global dimension D = 2; central charge c = 1/2.
 *  - Fusion rules: σ×σ=1+ψ, σ×ψ=σ, ψ×ψ=1.
 *  - Modular S-matrix unitarity (S†S = I) and symmetry.
 *  - Modular T-matrix entries: θ_1=1, θ_σ=exp(iπ/8), θ_ψ=-1.
 *  - CY invariants on canonical 4-manifolds:
 *      S⁴    (χ=2, σ=0)  → 1/4
 *      CP²   (χ=3, σ=1)  → (1/8)·exp(iπ/8)
 *      ‾CP²  (χ=3, σ=-1) → (1/8)·exp(-iπ/8)
 *      S²×S² (χ=4, σ=0)  → 1/16
 *      T⁴    (χ=0, σ=0)  → 1
 */
#include "harness.h"
#include <irrep/crane_yetter_ising.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>

static int approx_eq(double _Complex a, double _Complex b, double tol) {
    return cabs(a - b) < tol;
}

static int test_quantum_dimensions(void) {
    if (fabs(irrep_ising_quantum_dim(IRREP_ISING_OBJ_1)   - 1.0)      > 1e-12) return 1;
    if (fabs(irrep_ising_quantum_dim(IRREP_ISING_OBJ_SIGMA) - M_SQRT2) > 1e-12) return 1;
    if (fabs(irrep_ising_quantum_dim(IRREP_ISING_OBJ_PSI) - 1.0)      > 1e-12) return 1;
    if (fabs(irrep_ising_global_dim() - 2.0) > 1e-12) return 1;
    if (fabs(irrep_ising_central_charge() - 0.5) > 1e-12) return 1;
    return 0;
}

static int test_fusion_rules(void) {
    /* σ×σ = 1 + ψ */
    if (irrep_ising_fusion(IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_1)   != 1) return 1;
    if (irrep_ising_fusion(IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_SIGMA) != 0) return 1;
    if (irrep_ising_fusion(IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_PSI) != 1) return 1;
    /* σ×ψ = σ */
    if (irrep_ising_fusion(IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_PSI, IRREP_ISING_OBJ_SIGMA) != 1) return 1;
    if (irrep_ising_fusion(IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_PSI, IRREP_ISING_OBJ_1)   != 0) return 1;
    /* ψ×ψ = 1 */
    if (irrep_ising_fusion(IRREP_ISING_OBJ_PSI, IRREP_ISING_OBJ_PSI, IRREP_ISING_OBJ_1)   != 1) return 1;
    if (irrep_ising_fusion(IRREP_ISING_OBJ_PSI, IRREP_ISING_OBJ_PSI, IRREP_ISING_OBJ_PSI) != 0) return 1;
    /* Symmetry: N^{ab}_c = N^{ba}_c. */
    if (irrep_ising_fusion(IRREP_ISING_OBJ_PSI, IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_SIGMA) != 1) return 1;
    /* 1 × x = x. */
    if (irrep_ising_fusion(IRREP_ISING_OBJ_1, IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_SIGMA) != 1) return 1;
    if (irrep_ising_fusion(IRREP_ISING_OBJ_1, IRREP_ISING_OBJ_PSI, IRREP_ISING_OBJ_PSI)     != 1) return 1;
    return 0;
}

static int test_S_matrix_unitarity(void) {
    /* S† S should equal identity. Since S is real-symmetric here, this
     * is just S^2 = I (treating S as real). */
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            double _Complex sum = 0;
            for (int k = 0; k < 3; ++k) {
                sum += conj(irrep_ising_S_matrix((irrep_ising_object_t)a, (irrep_ising_object_t)k))
                       * irrep_ising_S_matrix((irrep_ising_object_t)k, (irrep_ising_object_t)b);
            }
            double _Complex expected = (a == b) ? 1.0 + 0.0*I : 0.0 + 0.0*I;
            if (!approx_eq(sum, expected, 1e-12)) return 1;
        }
    }
    /* Symmetry: S_ab = S_ba. */
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            if (!approx_eq(irrep_ising_S_matrix((irrep_ising_object_t)a, (irrep_ising_object_t)b),
                           irrep_ising_S_matrix((irrep_ising_object_t)b, (irrep_ising_object_t)a), 1e-12))
                return 1;
        }
    }
    return 0;
}

static int test_T_eigenvalues(void) {
    if (!approx_eq(irrep_ising_T_eigenvalue(IRREP_ISING_OBJ_1), 1.0, 1e-12)) return 1;
    /* θ_σ = exp(iπ/8) */
    double _Complex expected_sigma = cexp(I * M_PI / 8.0);
    if (!approx_eq(irrep_ising_T_eigenvalue(IRREP_ISING_OBJ_SIGMA), expected_sigma, 1e-12)) return 1;
    /* θ_ψ = -1 */
    if (!approx_eq(irrep_ising_T_eigenvalue(IRREP_ISING_OBJ_PSI), -1.0, 1e-12)) return 1;
    return 0;
}

static int test_CY_invariant(void) {
    /* S⁴: χ=2, σ=0. Z = 1/4. */
    if (!approx_eq(irrep_crane_yetter_ising_invariant(2, 0), 0.25, 1e-12)) return 1;
    /* T⁴: χ=0, σ=0. Z = 1. */
    if (!approx_eq(irrep_crane_yetter_ising_invariant(0, 0), 1.0, 1e-12)) return 1;
    /* CP²: χ=3, σ=1. Z = (1/8) · exp(iπ/8). */
    double _Complex z_cp2 = irrep_crane_yetter_ising_invariant(3, 1);
    double _Complex expected_cp2 = (1.0 / 8.0) * cexp(I * M_PI / 8.0);
    if (!approx_eq(z_cp2, expected_cp2, 1e-12)) return 1;
    /* ‾CP² (orientation-reversed): χ=3, σ=-1. Z = (1/8) · exp(-iπ/8). */
    double _Complex z_cp2_bar = irrep_crane_yetter_ising_invariant(3, -1);
    double _Complex expected_cp2_bar = (1.0 / 8.0) * cexp(-I * M_PI / 8.0);
    if (!approx_eq(z_cp2_bar, expected_cp2_bar, 1e-12)) return 1;
    /* CP² and ‾CP² are complex conjugates (orientation flips σ). */
    if (!approx_eq(z_cp2, conj(z_cp2_bar), 1e-12)) return 1;
    /* S²×S²: χ=4, σ=0. Z = 1/16. */
    if (!approx_eq(irrep_crane_yetter_ising_invariant(4, 0), 1.0 / 16.0, 1e-12)) return 1;
    return 0;
}

int main(void) {
    int rc = 0;
    if (test_quantum_dimensions())  { fprintf(stderr, "FAIL test_quantum_dimensions\n"); rc = 1; }
    if (test_fusion_rules())        { fprintf(stderr, "FAIL test_fusion_rules\n"); rc = 1; }
    if (test_S_matrix_unitarity())  { fprintf(stderr, "FAIL test_S_matrix_unitarity\n"); rc = 1; }
    if (test_T_eigenvalues())       { fprintf(stderr, "FAIL test_T_eigenvalues\n"); rc = 1; }
    if (test_CY_invariant())        { fprintf(stderr, "FAIL test_CY_invariant\n"); rc = 1; }
    return rc;
}
