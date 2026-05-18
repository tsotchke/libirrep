/* SPDX-License-Identifier: MIT */
/** @file crane_yetter_ising.c
 *  @brief Crane-Yetter / Walker-Wang runtime for the Ising MTC. */
#include <irrep/crane_yetter_ising.h>

#include <complex.h>
#include <math.h>

/* ====================================================================
 * Ising MTC modular data
 *
 * Objects: 0 = 1 (vacuum), 1 = σ (non-abelian anyon), 2 = ψ (fermion).
 * Quantum dimensions: d = (1, √2, 1).
 * Global dimension: D = √(d_1² + d_σ² + d_ψ²) = √(1 + 2 + 1) = 2.
 * Central charge: c = 1/2.
 * Conformal weights: h_1 = 0, h_σ = 1/16, h_ψ = 1/2.
 * Topological twists: θ_a = exp(2πi h_a).
 * ==================================================================== */

double
irrep_ising_quantum_dim(irrep_ising_object_t a)
{
    switch (a) {
        case IRREP_ISING_OBJ_1:     return 1.0;
        case IRREP_ISING_OBJ_SIGMA: return M_SQRT2;
        case IRREP_ISING_OBJ_PSI:   return 1.0;
        default:                    return 0.0;
    }
}

double
irrep_ising_global_dim(void)
{
    /* D = sqrt(1 + 2 + 1) = 2. */
    return 2.0;
}

double
irrep_ising_central_charge(void)
{
    return 0.5;
}

int
irrep_ising_fusion(irrep_ising_object_t a,
                   irrep_ising_object_t b,
                   irrep_ising_object_t c)
{
    /* Symmetric in (a, b). Sort so a ≤ b. */
    if (a > b) { irrep_ising_object_t t = a; a = b; b = t; }

    /* 1 × x = x. */
    if (a == IRREP_ISING_OBJ_1) {
        return (b == c) ? 1 : 0;
    }
    /* σ × σ = 1 + ψ. */
    if (a == IRREP_ISING_OBJ_SIGMA && b == IRREP_ISING_OBJ_SIGMA) {
        return (c == IRREP_ISING_OBJ_1 || c == IRREP_ISING_OBJ_PSI) ? 1 : 0;
    }
    /* σ × ψ = σ. */
    if (a == IRREP_ISING_OBJ_SIGMA && b == IRREP_ISING_OBJ_PSI) {
        return (c == IRREP_ISING_OBJ_SIGMA) ? 1 : 0;
    }
    /* ψ × ψ = 1. */
    if (a == IRREP_ISING_OBJ_PSI && b == IRREP_ISING_OBJ_PSI) {
        return (c == IRREP_ISING_OBJ_1) ? 1 : 0;
    }
    return 0;
}

/* Modular S-matrix entries.
 *
 * S = (1/D) [[d_a · d_b · (sum)]] — Verlinde formula. For Ising:
 *
 *     S = (1/2) [[ 1,  √2,  1 ],
 *                [ √2,  0, -√2],
 *                [ 1, -√2,  1 ]]
 *
 * Properties: S is symmetric, real, unitary (S^2 = C where C is charge
 * conjugation; for self-conjugate Ising, S^2 = I).
 */
double _Complex
irrep_ising_S_matrix(irrep_ising_object_t a, irrep_ising_object_t b)
{
    static const double S_real[3][3] = {
        { 0.5,            M_SQRT2 / 2.0,  0.5            },
        { M_SQRT2 / 2.0,  0.0,            -M_SQRT2 / 2.0 },
        { 0.5,            -M_SQRT2 / 2.0, 0.5            },
    };
    if (a < 0 || a >= 3 || b < 0 || b >= 3) return 0.0 + 0.0 * I;
    return (double _Complex)S_real[a][b];
}

double _Complex
irrep_ising_T_eigenvalue(irrep_ising_object_t a)
{
    /* T_a = exp(2πi h_a) modulo a central-charge phase exp(-2πi c/24).
     * Here we expose the bare topological twists θ_a = exp(2πi h_a):
     *   h_1 = 0    → 1
     *   h_σ = 1/16 → exp(iπ/8)
     *   h_ψ = 1/2  → -1
     */
    switch (a) {
        case IRREP_ISING_OBJ_1:
            return 1.0 + 0.0 * I;
        case IRREP_ISING_OBJ_SIGMA:
            return cexp(I * M_PI / 8.0);
        case IRREP_ISING_OBJ_PSI:
            return -1.0 + 0.0 * I;
        default:
            return 0.0 + 0.0 * I;
    }
}

/* ====================================================================
 * Crane-Yetter invariant on closed orientable 4-manifolds.
 *
 * Z_CY(M) = D^{-χ(M)} · exp(2πi c σ(M) / 8)
 *         = 2^{-χ} · exp(iπ σ / 8)   [Ising case: D=2, c=1/2]
 * ==================================================================== */

double _Complex
irrep_crane_yetter_ising_invariant(int euler_char, int signature)
{
    double D = irrep_ising_global_dim();          /* 2 */
    double c = irrep_ising_central_charge();      /* 1/2 */
    double prefactor = pow(D, -(double)euler_char);
    double phase = M_PI * c * (double)signature / 4.0;
    /* exp(2πi c σ / 8) = exp(iπcσ/4). For c = 1/2: phase = π σ / 8. */
    return prefactor * (cos(phase) + I * sin(phase));
}
