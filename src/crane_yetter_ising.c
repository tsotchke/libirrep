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

/* ====================================================================
 * F-symbols and R-symbols of the Ising MTC.
 *
 * Convention: Kitaev's Anyons-paper Appendix E. Frobenius-Schur
 * indicator κ_σ = +1, so F^{σbσ}_1 = +1 for all b.
 * ==================================================================== */

double _Complex
irrep_ising_F_symbol(irrep_ising_object_t a, irrep_ising_object_t b,
                     irrep_ising_object_t c, irrep_ising_object_t d,
                     irrep_ising_object_t e, irrep_ising_object_t f)
{
    if (irrep_ising_fusion(a, b, e) == 0) return 0.0;
    if (irrep_ising_fusion(e, c, d) == 0) return 0.0;
    if (irrep_ising_fusion(b, c, f) == 0) return 0.0;
    if (irrep_ising_fusion(a, f, d) == 0) return 0.0;

    /* Non-trivial Hadamard block: F^{σσσ}_σ with (e, f) ∈ {1, ψ}². */
    if (a == IRREP_ISING_OBJ_SIGMA && b == IRREP_ISING_OBJ_SIGMA
        && c == IRREP_ISING_OBJ_SIGMA && d == IRREP_ISING_OBJ_SIGMA) {
        const double inv_sqrt2 = 0.7071067811865475244;
        int sign = (e == IRREP_ISING_OBJ_PSI && f == IRREP_ISING_OBJ_PSI) ? -1 : +1;
        return (double)sign * inv_sqrt2;
    }
    return 1.0;
}

double _Complex
irrep_ising_R_symbol(irrep_ising_object_t a, irrep_ising_object_t b,
                     irrep_ising_object_t c)
{
    if (irrep_ising_fusion(a, b, c) == 0) return 0.0;

    /* Trivial: any factor is identity. */
    if (a == IRREP_ISING_OBJ_1 || b == IRREP_ISING_OBJ_1) {
        return 1.0;
    }
    /* σ × σ → 1 or ψ. */
    if (a == IRREP_ISING_OBJ_SIGMA && b == IRREP_ISING_OBJ_SIGMA) {
        if (c == IRREP_ISING_OBJ_1)   return cexp(-I * M_PI / 8.0);
        if (c == IRREP_ISING_OBJ_PSI) return cexp(I * 3.0 * M_PI / 8.0);
    }
    /* σ × ψ → σ and ψ × σ → σ. */
    if ((a == IRREP_ISING_OBJ_SIGMA && b == IRREP_ISING_OBJ_PSI)
        || (a == IRREP_ISING_OBJ_PSI && b == IRREP_ISING_OBJ_SIGMA)) {
        return I;
    }
    /* ψ × ψ → 1. */
    if (a == IRREP_ISING_OBJ_PSI && b == IRREP_ISING_OBJ_PSI) {
        return -1.0;
    }
    return 0.0;
}

double
irrep_ising_twist_from_R_residual(void)
{
    double max_err = 0.0;
    for (int a = 0; a < IRREP_ISING_N_OBJECTS; ++a) {
        irrep_ising_object_t A = (irrep_ising_object_t)a;
        double d_a = irrep_ising_quantum_dim(A);
        double _Complex theta = 0.0;
        for (int c = 0; c < IRREP_ISING_N_OBJECTS; ++c) {
            irrep_ising_object_t C = (irrep_ising_object_t)c;
            int N = irrep_ising_fusion(A, A, C);
            if (N == 0) continue;
            double d_c = irrep_ising_quantum_dim(C);
            double _Complex R = irrep_ising_R_symbol(A, A, C);
            theta += (double)N * d_c * R;
        }
        theta /= d_a;
        double _Complex T_hardcoded = irrep_ising_T_eigenvalue(A);
        double err = cabs(theta - T_hardcoded);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

double
irrep_ising_S_from_twist_residual(void)
{
    double max_err = 0.0;
    double D = irrep_ising_global_dim();
    for (int a = 0; a < IRREP_ISING_N_OBJECTS; ++a) {
        for (int b = 0; b < IRREP_ISING_N_OBJECTS; ++b) {
            irrep_ising_object_t A = (irrep_ising_object_t)a;
            irrep_ising_object_t B = (irrep_ising_object_t)b;
            double _Complex theta_a = irrep_ising_T_eigenvalue(A);
            double _Complex theta_b = irrep_ising_T_eigenvalue(B);
            double _Complex S_derived = 0.0;
            for (int c = 0; c < IRREP_ISING_N_OBJECTS; ++c) {
                irrep_ising_object_t C = (irrep_ising_object_t)c;
                int N = irrep_ising_fusion(A, B, C);
                if (N == 0) continue;
                double _Complex theta_c = irrep_ising_T_eigenvalue(C);
                double d_c = irrep_ising_quantum_dim(C);
                S_derived += (double)N * theta_c / (theta_a * theta_b) * d_c;
            }
            S_derived /= D;
            double _Complex S_hard = irrep_ising_S_matrix(A, B);
            double err = cabs(S_derived - S_hard);
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

double
irrep_ising_F_unitarity_residual(void)
{
    double max_err = 0.0;
    irrep_ising_object_t s = IRREP_ISING_OBJ_SIGMA;
    irrep_ising_object_t ch[2] = {
        IRREP_ISING_OBJ_1, IRREP_ISING_OBJ_PSI,
    };
    for (int fi = 0; fi < 2; ++fi) {
        for (int fj = 0; fj < 2; ++fj) {
            double _Complex acc = 0.0;
            for (int ei = 0; ei < 2; ++ei) {
                double _Complex Fi = irrep_ising_F_symbol(s, s, s, s, ch[ei], ch[fi]);
                double _Complex Fj = irrep_ising_F_symbol(s, s, s, s, ch[ei], ch[fj]);
                acc += Fi * conj(Fj);
            }
            double expected = (fi == fj) ? 1.0 : 0.0;
            double err = cabs(acc - expected);
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

/* ====================================================================
 * Walker-Wang 3+1D Hamiltonian: vertex admissibility
 * ==================================================================== */

int
irrep_ising_walker_wang_vertex_admissible(
    const irrep_ising_object_t *edge_labels, int n)
{
    if (edge_labels == NULL || n < 0) return 0;
    int sigma = 0, psi = 0;
    for (int i = 0; i < n; ++i) {
        switch (edge_labels[i]) {
            case IRREP_ISING_OBJ_SIGMA: ++sigma; break;
            case IRREP_ISING_OBJ_PSI:   ++psi;   break;
            default: break;
        }
    }
    if (sigma % 2 != 0) return 0;
    if (sigma == 0 && (psi % 2) != 0) return 0;
    return 1;
}

long long
irrep_ising_walker_wang_admissible_count(int n)
{
    if (n < 0) return 0;
    if (n == 0) return 1; /* Empty configuration trivially admissible. */
    /* Enumerate all 3^n label tuples. */
    long long total = 1;
    for (int i = 0; i < n; ++i) total *= 3;
    irrep_ising_object_t labels[64];
    if (n > 64) return -1; /* Sanity bound. */
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        long long x = c;
        for (int i = 0; i < n; ++i) {
            labels[i] = (irrep_ising_object_t)(x % 3);
            x /= 3;
        }
        if (irrep_ising_walker_wang_vertex_admissible(labels, n)) ++count;
    }
    return count;
}

/* 3-simplex / tetrahedron edge indexing:
 *   edge index 0 = (0, 1)
 *   edge index 1 = (0, 2)
 *   edge index 2 = (0, 3)
 *   edge index 3 = (1, 2)
 *   edge index 4 = (1, 3)
 *   edge index 5 = (2, 3)
 *
 * Vertex incidences (3 edges per vertex):
 *   v0: {0, 1, 2}   v1: {0, 3, 4}   v2: {1, 3, 5}   v3: {2, 4, 5}
 *
 * Face incidences (3 edges per face):
 *   f(0,1,2): {0, 1, 3}   f(0,1,3): {0, 2, 4}
 *   f(0,2,3): {1, 2, 5}   f(1,2,3): {3, 4, 5}
 */
static const int kSimplex3Vertices[4][3] = {
    {0, 1, 2}, {0, 3, 4}, {1, 3, 5}, {2, 4, 5},
};
static const int kSimplex3Faces[4][3] = {
    {0, 1, 3}, {0, 2, 4}, {1, 2, 5}, {3, 4, 5},
};

static long long
simplex3_count(int with_faces)
{
    irrep_ising_object_t labels[6];
    long long total = 1;
    for (int i = 0; i < 6; ++i) total *= 3;
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        long long x = c;
        for (int i = 0; i < 6; ++i) {
            labels[i] = (irrep_ising_object_t)(x % 3);
            x /= 3;
        }
        int ok = 1;
        for (int v = 0; v < 4 && ok; ++v) {
            irrep_ising_object_t local[3] = {
                labels[kSimplex3Vertices[v][0]],
                labels[kSimplex3Vertices[v][1]],
                labels[kSimplex3Vertices[v][2]],
            };
            if (!irrep_ising_walker_wang_vertex_admissible(local, 3)) ok = 0;
        }
        if (ok && with_faces) {
            for (int f = 0; f < 4 && ok; ++f) {
                irrep_ising_object_t local[3] = {
                    labels[kSimplex3Faces[f][0]],
                    labels[kSimplex3Faces[f][1]],
                    labels[kSimplex3Faces[f][2]],
                };
                if (!irrep_ising_walker_wang_vertex_admissible(local, 3)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
}

long long
irrep_ising_walker_wang_simplex3_vertex_count(void)
{
    return simplex3_count(/*with_faces=*/0);
}

long long
irrep_ising_walker_wang_simplex3_full_count(void)
{
    return simplex3_count(/*with_faces=*/1);
}

double _Complex
irrep_ising_walker_wang_plaquette_psi_phase(
    const irrep_ising_object_t *boundary_labels, int n)
{
    if (boundary_labels == NULL || n < 0) return 0.0;
    double _Complex phase = 1.0 + 0.0 * I;
    for (int i = 0; i < n; ++i) {
        switch (boundary_labels[i]) {
            case IRREP_ISING_OBJ_1:                 break;
            case IRREP_ISING_OBJ_SIGMA: phase *= I; break;
            case IRREP_ISING_OBJ_PSI:   phase *= -1.0; break;
            default: return 0.0;
        }
    }
    return phase;
}
