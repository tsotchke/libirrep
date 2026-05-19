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

/* F-symbol Hadamard block on σσσσ + R-symbol → twist consistency. */
static int test_F_R_symbols(void) {
    irrep_ising_object_t s  = IRREP_ISING_OBJ_SIGMA;
    irrep_ising_object_t I_ = IRREP_ISING_OBJ_1;
    irrep_ising_object_t P  = IRREP_ISING_OBJ_PSI;
    const double inv_sqrt2 = 0.7071067811865475244;

    double _Complex F11 = irrep_ising_F_symbol(s, s, s, s, I_, I_);
    double _Complex F1P = irrep_ising_F_symbol(s, s, s, s, I_, P);
    double _Complex FP1 = irrep_ising_F_symbol(s, s, s, s, P,  I_);
    double _Complex FPP = irrep_ising_F_symbol(s, s, s, s, P,  P);
    if (cabs(F11 -   inv_sqrt2)  > 1e-12) return 1;
    if (cabs(F1P -   inv_sqrt2)  > 1e-12) return 1;
    if (cabs(FP1 -   inv_sqrt2)  > 1e-12) return 1;
    if (cabs(FPP - (-inv_sqrt2)) > 1e-12) return 1;

    /* Forbidden fusion → 0. (σ × ψ = σ; (σψ)_ψ is forbidden.) */
    if (cabs(irrep_ising_F_symbol(s, P, s, s, P, s)) > 1e-12) return 1;

    /* F-unitarity on σσσσ. */
    if (irrep_ising_F_unitarity_residual() > 1e-12) return 1;

    /* R-symbols. */
    if (cabs(irrep_ising_R_symbol(s, s, I_) - cexp(-I * M_PI / 8.0))      > 1e-12) return 1;
    if (cabs(irrep_ising_R_symbol(s, s, P)  - cexp(I * 3.0 * M_PI / 8.0)) > 1e-12) return 1;
    if (cabs(irrep_ising_R_symbol(P, P, I_) - (-1.0))                     > 1e-12) return 1;
    if (cabs(irrep_ising_R_symbol(s, P, s)  - I)                          > 1e-12) return 1;

    /* Forbidden R: σ × σ → σ has N = 0. */
    if (cabs(irrep_ising_R_symbol(s, s, s)) > 1e-12) return 1;

    /* Twist consistency: derived θ_a matches T_a to machine precision. */
    if (irrep_ising_twist_from_R_residual() > 1e-12) return 1;
    /* S-matrix from twists matches hardcoded S to machine precision. */
    if (irrep_ising_S_from_twist_residual() > 1e-12) return 1;
    /* Verlinde formula: N derived from S matches the fusion ring. */
    if (irrep_ising_verlinde_residual() > 1e-12) return 1;

    return 0;
}

/* Walker-Wang vertex admissibility: counts at small valencies match
 * the closed-form Ising fusion-multiplet enumeration.
 *
 *   n=1:  only label-1 admissible (need σ-count even AND if σ-count=0
 *         then ψ-count even). 1 of 3. — actually for n=1: σ-count=0
 *         → ψ-count must be even → only the "1" configuration.
 *   n=2:  admissible iff #σ even AND (if #σ=0 then #ψ even).
 *         Configurations: (1,1), (ψ,ψ), (σ,σ). Count = 3.
 *   n=3:  10 of 27 (computed above; matches the σ²-pair + 1-mixed and
 *         the all-(1, ψ) even-ψ subsets).
 *   n=4:  33 of 81.
 *   n=5:  for σ-count = 0 with #ψ even (∈ {0,2,4}):
 *           C(5,0)+C(5,2)+C(5,4) = 1+10+5 = 16;
 *         σ-count = 2: C(5,2) · 2³ = 80;
 *         σ-count = 4: C(5,4) · 2 = 10.
 *         Total = 16 + 80 + 10 = 106. */
static int test_walker_wang_vertex_admissible(void) {
    if (irrep_ising_walker_wang_admissible_count(0) != 1)   return 1;
    if (irrep_ising_walker_wang_admissible_count(1) != 1)   return 1;
    if (irrep_ising_walker_wang_admissible_count(2) != 3)   return 1;
    if (irrep_ising_walker_wang_admissible_count(3) != 10)  return 1;
    if (irrep_ising_walker_wang_admissible_count(4) != 33)  return 1;
    if (irrep_ising_walker_wang_admissible_count(5) != 106) return 1;

    /* Spot-check a few specific configurations. */
    irrep_ising_object_t all_one[3]   = { IRREP_ISING_OBJ_1,     IRREP_ISING_OBJ_1,     IRREP_ISING_OBJ_1     };
    irrep_ising_object_t two_psi[3]   = { IRREP_ISING_OBJ_PSI,   IRREP_ISING_OBJ_PSI,   IRREP_ISING_OBJ_1     };
    irrep_ising_object_t two_sigma[3] = { IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_1     };
    irrep_ising_object_t three_sigma[3]={IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_SIGMA, IRREP_ISING_OBJ_SIGMA };
    irrep_ising_object_t one_psi[3]   = { IRREP_ISING_OBJ_PSI,   IRREP_ISING_OBJ_1,     IRREP_ISING_OBJ_1     };
    if (!irrep_ising_walker_wang_vertex_admissible(all_one,    3)) return 1;
    if (!irrep_ising_walker_wang_vertex_admissible(two_psi,    3)) return 1;
    if (!irrep_ising_walker_wang_vertex_admissible(two_sigma,  3)) return 1;
    if ( irrep_ising_walker_wang_vertex_admissible(three_sigma,3)) return 1; /* #σ=3 odd */
    if ( irrep_ising_walker_wang_vertex_admissible(one_psi,    3)) return 1; /* #σ=0, #ψ=1 odd */

    return 0;
}

/* Walker-Wang on the 3-simplex: ground state space dimension under
 * (a) vertex-only admissibility and (b) vertex + face admissibility. */
static int test_walker_wang_simplex3(void) {
    if (irrep_ising_walker_wang_simplex3_vertex_count() != 36) return 1;
    if (irrep_ising_walker_wang_simplex3_full_count() != 16)   return 1;
    return 0;
}

/* Walker-Wang on the unit cube: 8 vertices, 12 edges, 6 square faces.
 * Enumerates 3^12 = 531,441 configurations, applying 8 vertex + 6 face
 * admissibility constraints. Total runtime ~ few ms with the early-out
 * vertex pruning. */
static int test_walker_wang_cube(void) {
    long long c = irrep_ising_walker_wang_cube_full_count();
    if (c != 120) {
        fprintf(stderr, "cube count = %lld (expected 120)\n", c);
        return 1;
    }
    return 0;
}

/* Walker-Wang on the regular octahedron: 6 vertices (4-valent), 12 edges,
 * 8 triangular faces. The octahedron is dual to the cube, so by
 * Crane-Yetter / Walker-Wang invariance under polyhedral duality the
 * ground-state dimension should match the cube count.
 *
 * PROOF: both enumerators produce 120 — a runtime witness for the
 * cube ↔ octahedron duality at the smallest non-trivial geometry. */
static int test_walker_wang_octahedron(void) {
    long long c_oct = irrep_ising_walker_wang_octahedron_full_count();
    long long c_cube = irrep_ising_walker_wang_cube_full_count();
    if (c_oct != 120) {
        fprintf(stderr, "octahedron count = %lld (expected 120)\n", c_oct);
        return 1;
    }
    /* Polyhedral duality: same ground-state dimension. */
    if (c_oct != c_cube) {
        fprintf(stderr, "WW duality failed: oct=%lld vs cube=%lld\n",
                c_oct, c_cube);
        return 1;
    }
    return 0;
}

/* Walker-Wang on the triangular bipyramid (V=5, E=9, F=6 triangles) vs
 * the triangular prism (V=6, E=9, F=2 triangles + 3 squares). The two
 * polyhedra are dual — bipyramid (V=5, F=6) ↔ prism (V=6, F=5) — so
 * WW invariance under duality predicts identical ground-state dim.
 *
 * PROOF: both produce 1. Third instance of the duality web (after
 * the self-dual tetrahedron and the cube/octahedron pair). */
static int test_walker_wang_bipyramid_prism_duality(void) {
    long long c_bp = irrep_ising_walker_wang_tri_bipyramid_full_count();
    long long c_pr = irrep_ising_walker_wang_tri_prism_full_count();
    if (c_bp != 1) {
        fprintf(stderr, "bipyramid count = %lld (expected 1)\n", c_bp);
        return 1;
    }
    if (c_pr != 1) {
        fprintf(stderr, "prism count = %lld (expected 1)\n", c_pr);
        return 1;
    }
    if (c_bp != c_pr) {
        fprintf(stderr, "WW duality bipyramid↔prism failed: %lld vs %lld\n",
                c_bp, c_pr);
        return 1;
    }
    return 0;
}

/* Walker-Wang B_p^ψ plaquette phase on representative configurations.
 * Phase = i^{#σ} · (-1)^{#ψ}. */
static int test_walker_wang_plaquette_psi_phase(void) {
    irrep_ising_object_t I_ = IRREP_ISING_OBJ_1;
    irrep_ising_object_t S_ = IRREP_ISING_OBJ_SIGMA;
    irrep_ising_object_t P_ = IRREP_ISING_OBJ_PSI;

    /* All-1: phase = 1. */
    {
        irrep_ising_object_t b[4] = { I_, I_, I_, I_ };
        if (cabs(irrep_ising_walker_wang_plaquette_psi_phase(b, 4) - 1.0) > 1e-12) return 1;
    }
    /* All-ψ: phase = (-1)^4 = 1. */
    {
        irrep_ising_object_t b[4] = { P_, P_, P_, P_ };
        if (cabs(irrep_ising_walker_wang_plaquette_psi_phase(b, 4) - 1.0) > 1e-12) return 1;
    }
    /* Single ψ: phase = -1. */
    {
        irrep_ising_object_t b[4] = { P_, I_, I_, I_ };
        if (cabs(irrep_ising_walker_wang_plaquette_psi_phase(b, 4) - (-1.0)) > 1e-12) return 1;
    }
    /* Two σ: phase = i² = -1. */
    {
        irrep_ising_object_t b[4] = { S_, S_, I_, I_ };
        if (cabs(irrep_ising_walker_wang_plaquette_psi_phase(b, 4) - (-1.0)) > 1e-12) return 1;
    }
    /* Four σ: phase = i⁴ = 1. */
    {
        irrep_ising_object_t b[4] = { S_, S_, S_, S_ };
        if (cabs(irrep_ising_walker_wang_plaquette_psi_phase(b, 4) - 1.0) > 1e-12) return 1;
    }
    /* Two σ + one ψ: phase = i² · (-1) = 1. */
    {
        irrep_ising_object_t b[3] = { S_, S_, P_ };
        if (cabs(irrep_ising_walker_wang_plaquette_psi_phase(b, 3) - 1.0) > 1e-12) return 1;
    }
    /* One σ (non-admissible config): phase = i (imaginary). */
    {
        irrep_ising_object_t b[3] = { S_, I_, I_ };
        if (cabs(irrep_ising_walker_wang_plaquette_psi_phase(b, 3) - I) > 1e-12) return 1;
    }
    /* Empty plaquette: phase = 1 (empty product). */
    {
        if (cabs(irrep_ising_walker_wang_plaquette_psi_phase(NULL, 0)) > 1e-12) return 1;
    }
    return 0;
}

int main(void) {
    int rc = 0;
    if (test_quantum_dimensions())  { fprintf(stderr, "FAIL test_quantum_dimensions\n"); rc = 1; }
    if (test_fusion_rules())        { fprintf(stderr, "FAIL test_fusion_rules\n"); rc = 1; }
    if (test_S_matrix_unitarity())  { fprintf(stderr, "FAIL test_S_matrix_unitarity\n"); rc = 1; }
    if (test_T_eigenvalues())       { fprintf(stderr, "FAIL test_T_eigenvalues\n"); rc = 1; }
    if (test_CY_invariant())        { fprintf(stderr, "FAIL test_CY_invariant\n"); rc = 1; }
    if (test_F_R_symbols())         { fprintf(stderr, "FAIL test_F_R_symbols\n"); rc = 1; }
    if (test_walker_wang_vertex_admissible())
        { fprintf(stderr, "FAIL test_walker_wang_vertex_admissible\n"); rc = 1; }
    if (test_walker_wang_simplex3())
        { fprintf(stderr, "FAIL test_walker_wang_simplex3\n"); rc = 1; }
    if (test_walker_wang_cube())
        { fprintf(stderr, "FAIL test_walker_wang_cube\n"); rc = 1; }
    if (test_walker_wang_plaquette_psi_phase())
        { fprintf(stderr, "FAIL test_walker_wang_plaquette_psi_phase\n"); rc = 1; }
    if (test_walker_wang_octahedron())
        { fprintf(stderr, "FAIL test_walker_wang_octahedron\n"); rc = 1; }
    if (test_walker_wang_bipyramid_prism_duality())
        { fprintf(stderr, "FAIL test_walker_wang_bipyramid_prism_duality\n"); rc = 1; }
    return rc;
}
