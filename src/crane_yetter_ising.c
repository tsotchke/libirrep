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
irrep_crane_yetter_connected_sum_residual(void)
{
    /* Test pairs (M, N) and the resulting connected sum M#N. */
    struct { int chi_m, sig_m, chi_n, sig_n; } pairs[] = {
        { 3,  1, 3,  1 },  /* CP² # CP² */
        { 3, -1, 3, -1 },  /* ‾CP² # ‾CP² */
        { 3,  1, 3, -1 },  /* CP² # ‾CP² */
        { 4,  0, 4,  0 },  /* (S² × S²) # (S² × S²) */
        { 2,  0, 3,  1 },  /* S⁴ # CP² = CP² (test) */
    };
    int n_pairs = sizeof(pairs) / sizeof(pairs[0]);
    double max_err = 0.0;
    double _Complex Z_S4 = irrep_crane_yetter_ising_invariant(2, 0);
    for (int i = 0; i < n_pairs; ++i) {
        int chi_sum = pairs[i].chi_m + pairs[i].chi_n - 2;
        int sig_sum = pairs[i].sig_m + pairs[i].sig_n;
        double _Complex Z_MN = irrep_crane_yetter_ising_invariant(chi_sum, sig_sum);
        double _Complex Z_M  = irrep_crane_yetter_ising_invariant(pairs[i].chi_m, pairs[i].sig_m);
        double _Complex Z_N  = irrep_crane_yetter_ising_invariant(pairs[i].chi_n, pairs[i].sig_n);
        double _Complex RHS  = Z_M * Z_N / Z_S4;
        double err = cabs(Z_MN - RHS);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

double
irrep_ising_S_squared_residual(void)
{
    double max_err = 0.0;
    for (int a = 0; a < IRREP_ISING_N_OBJECTS; ++a) {
        for (int b = 0; b < IRREP_ISING_N_OBJECTS; ++b) {
            irrep_ising_object_t A = (irrep_ising_object_t)a;
            irrep_ising_object_t B = (irrep_ising_object_t)b;
            /* (S²)_{ab} = Σ_x S_{ax} · S_{xb} */
            double _Complex acc = 0.0;
            for (int x = 0; x < IRREP_ISING_N_OBJECTS; ++x) {
                irrep_ising_object_t X = (irrep_ising_object_t)x;
                acc += irrep_ising_S_matrix(A, X) * irrep_ising_S_matrix(X, B);
            }
            double expected = (a == b) ? 1.0 : 0.0;
            double err = cabs(acc - expected);
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

double
irrep_ising_verlinde_residual(void)
{
    double max_err = 0.0;
    for (int a = 0; a < IRREP_ISING_N_OBJECTS; ++a) {
        for (int b = 0; b < IRREP_ISING_N_OBJECTS; ++b) {
            for (int c = 0; c < IRREP_ISING_N_OBJECTS; ++c) {
                irrep_ising_object_t A = (irrep_ising_object_t)a;
                irrep_ising_object_t B = (irrep_ising_object_t)b;
                irrep_ising_object_t C = (irrep_ising_object_t)c;
                double _Complex acc = 0.0;
                for (int x = 0; x < IRREP_ISING_N_OBJECTS; ++x) {
                    irrep_ising_object_t X = (irrep_ising_object_t)x;
                    double _Complex S_ax = irrep_ising_S_matrix(A, X);
                    double _Complex S_bx = irrep_ising_S_matrix(B, X);
                    double _Complex S_cx = irrep_ising_S_matrix(C, X);
                    double _Complex S_0x = irrep_ising_S_matrix(IRREP_ISING_OBJ_1, X);
                    if (cabs(S_0x) < 1e-15) continue;
                    acc += S_ax * S_bx * conj(S_cx) / S_0x;
                }
                int N_actual = irrep_ising_fusion(A, B, C);
                double err = cabs(acc - (double)N_actual);
                if (err > max_err) max_err = err;
            }
        }
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

/* Unit cube edge indexing:
 *
 *   Vertices (x, y, z) ∈ {0, 1}³, with binary index v = 4z + 2y + x.
 *   12 edges, grouped by direction:
 *     x-direction (4 edges): edge[0..3] = (0-1), (2-3), (4-5), (6-7)
 *     y-direction (4 edges): edge[4..7] = (0-2), (1-3), (4-6), (5-7)
 *     z-direction (4 edges): edge[8..11] = (0-4), (1-5), (2-6), (3-7)
 *
 * Each vertex has 3 incident edges (one per direction). Each face has
 * 4 boundary edges (a closed cycle around the face).
 *
 *   Vertex 0 (0,0,0): x-edge 0 (0-1), y-edge 4 (0-2), z-edge  8 (0-4)
 *   Vertex 1 (1,0,0): x-edge 0 (0-1), y-edge 5 (1-3), z-edge  9 (1-5)
 *   Vertex 2 (0,1,0): x-edge 1 (2-3), y-edge 4 (0-2), z-edge 10 (2-6)
 *   Vertex 3 (1,1,0): x-edge 1 (2-3), y-edge 5 (1-3), z-edge 11 (3-7)
 *   Vertex 4 (0,0,1): x-edge 2 (4-5), y-edge 6 (4-6), z-edge  8 (0-4)
 *   Vertex 5 (1,0,1): x-edge 2 (4-5), y-edge 7 (5-7), z-edge  9 (1-5)
 *   Vertex 6 (0,1,1): x-edge 3 (6-7), y-edge 6 (4-6), z-edge 10 (2-6)
 *   Vertex 7 (1,1,1): x-edge 3 (6-7), y-edge 7 (5-7), z-edge 11 (3-7)
 *
 * Faces (6 total):
 *   z=0 (bottom): boundary edges {0, 4, 1, 5}  (square 0-1-3-2)
 *   z=1 (top):    {2, 6, 3, 7}                (square 4-5-7-6)
 *   y=0 (front):  {0, 8, 2, 9}                (square 0-1-5-4)
 *   y=1 (back):   {1, 10, 3, 11}              (square 2-3-7-6)
 *   x=0 (left):   {4, 8, 6, 10}               (square 0-2-6-4)
 *   x=1 (right):  {5, 9, 7, 11}               (square 1-3-7-5)
 */
static const int kCubeVertices[8][3] = {
    {0, 4,  8}, {0, 5,  9}, {1, 4, 10}, {1, 5, 11},
    {2, 6,  8}, {2, 7,  9}, {3, 6, 10}, {3, 7, 11},
};
static const int kCubeFaces[6][4] = {
    {0, 4, 1, 5}, {2, 6, 3, 7},
    {0, 8, 2, 9}, {1, 10, 3, 11},
    {4, 8, 6, 10}, {5, 9, 7, 11},
};

/* Octahedron incidences:
 *   Vertices: T (top), M0..M3 (middle ring), B (bottom). 6 vertices total.
 *   Edges (12):
 *     0: T-M0   1: T-M1   2: T-M2   3: T-M3
 *     4: B-M0   5: B-M1   6: B-M2   7: B-M3
 *     8: M0-M1  9: M1-M2 10: M2-M3 11: M3-M0
 *   Each vertex has 4 incident edges (octahedron is 4-valent).
 *   Each triangular face has 3 boundary edges. */
static const int kOctahedronVertices[6][4] = {
    {  0, 1, 2,  3 },     /* T:   TM0, TM1, TM2, TM3 */
    {  4, 5, 6,  7 },     /* B:   BM0, BM1, BM2, BM3 */
    {  0, 4, 8, 11 },     /* M0:  TM0, BM0, M0M1, M3M0 */
    {  1, 5, 8,  9 },     /* M1:  TM1, BM1, M0M1, M1M2 */
    {  2, 6, 9, 10 },     /* M2:  TM2, BM2, M1M2, M2M3 */
    {  3, 7, 10, 11 },    /* M3:  TM3, BM3, M2M3, M3M0 */
};
static const int kOctahedronFaces[8][3] = {
    { 0, 1,  8 },  /* T-M0-M1 */
    { 1, 2,  9 },  /* T-M1-M2 */
    { 2, 3, 10 },  /* T-M2-M3 */
    { 3, 0, 11 },  /* T-M3-M0 */
    { 4, 5,  8 },  /* B-M0-M1 */
    { 5, 6,  9 },  /* B-M1-M2 */
    { 6, 7, 10 },  /* B-M2-M3 */
    { 7, 4, 11 },  /* B-M3-M0 */
};

/* Triangular bipyramid incidences:
 *   Vertices: T (top apex), B (bottom apex), R0, R1, R2 (ring). 5 total.
 *   Edges (9):
 *     0: T-R0  1: T-R1  2: T-R2
 *     3: B-R0  4: B-R1  5: B-R2
 *     6: R0-R1 7: R1-R2 8: R2-R0
 *   T, B are 3-valent; R0, R1, R2 are 4-valent.
 *   6 triangular faces (3 top + 3 bottom). */
static const int kBipyramidVertexSizes[5] = { 3, 3, 4, 4, 4 };
static const int kBipyramidVertices[5][4] = {
    {  0,  1,  2, -1 },     /* T: TR0, TR1, TR2 */
    {  3,  4,  5, -1 },     /* B: BR0, BR1, BR2 */
    {  0,  3,  6,  8 },     /* R0 */
    {  1,  4,  6,  7 },     /* R1 */
    {  2,  5,  7,  8 },     /* R2 */
};
static const int kBipyramidFaces[6][3] = {
    { 0, 1, 6 },  /* T-R0-R1 */
    { 1, 2, 7 },  /* T-R1-R2 */
    { 2, 0, 8 },  /* T-R2-R0 */
    { 3, 4, 6 },  /* B-R0-R1 */
    { 4, 5, 7 },  /* B-R1-R2 */
    { 5, 3, 8 },  /* B-R2-R0 */
};

long long
irrep_ising_walker_wang_tri_bipyramid_full_count(void)
{
    irrep_ising_object_t labels[9];
    long long total = 1;
    for (int i = 0; i < 9; ++i) total *= 3;
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        long long x = c;
        for (int i = 0; i < 9; ++i) {
            labels[i] = (irrep_ising_object_t)(x % 3);
            x /= 3;
        }
        int ok = 1;
        for (int v = 0; v < 5 && ok; ++v) {
            int nv = kBipyramidVertexSizes[v];
            irrep_ising_object_t local[4];
            for (int k = 0; k < nv; ++k) {
                local[k] = labels[kBipyramidVertices[v][k]];
            }
            if (!irrep_ising_walker_wang_vertex_admissible(local, nv)) ok = 0;
        }
        if (ok) {
            for (int f = 0; f < 6 && ok; ++f) {
                irrep_ising_object_t local[3] = {
                    labels[kBipyramidFaces[f][0]],
                    labels[kBipyramidFaces[f][1]],
                    labels[kBipyramidFaces[f][2]],
                };
                if (!irrep_ising_walker_wang_vertex_admissible(local, 3)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
}

/* Triangular prism incidences:
 *   Vertices: T0, T1, T2 (top), B0, B1, B2 (bottom). 6 total, all 3-valent.
 *   Edges (9):
 *     0: T0-T1  1: T1-T2  2: T2-T0
 *     3: B0-B1  4: B1-B2  5: B2-B0
 *     6: T0-B0  7: T1-B1  8: T2-B2
 *   Faces (5): 2 triangles (top, bottom) + 3 squares (sides). */
static const int kPrismVertices[6][3] = {
    { 0, 2, 6 },  /* T0 */
    { 0, 1, 7 },  /* T1 */
    { 1, 2, 8 },  /* T2 */
    { 3, 5, 6 },  /* B0 */
    { 3, 4, 7 },  /* B1 */
    { 4, 5, 8 },  /* B2 */
};
static const int kPrismFaceSizes[5] = { 3, 3, 4, 4, 4 };
static const int kPrismFaces[5][4] = {
    { 0, 1, 2, -1 },  /* Top triangle T0-T1-T2 */
    { 3, 4, 5, -1 },  /* Bottom triangle B0-B1-B2 */
    { 0, 7, 3,  6 },  /* Side square T0-T1-B1-B0 */
    { 1, 8, 4,  7 },  /* Side square T1-T2-B2-B1 */
    { 2, 6, 5,  8 },  /* Side square T2-T0-B0-B2 */
};

long long
irrep_ising_walker_wang_tri_prism_full_count(void)
{
    irrep_ising_object_t labels[9];
    long long total = 1;
    for (int i = 0; i < 9; ++i) total *= 3;
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        long long x = c;
        for (int i = 0; i < 9; ++i) {
            labels[i] = (irrep_ising_object_t)(x % 3);
            x /= 3;
        }
        int ok = 1;
        for (int v = 0; v < 6 && ok; ++v) {
            irrep_ising_object_t local[3] = {
                labels[kPrismVertices[v][0]],
                labels[kPrismVertices[v][1]],
                labels[kPrismVertices[v][2]],
            };
            if (!irrep_ising_walker_wang_vertex_admissible(local, 3)) ok = 0;
        }
        if (ok) {
            for (int f = 0; f < 5 && ok; ++f) {
                int nf = kPrismFaceSizes[f];
                irrep_ising_object_t local[4];
                for (int k = 0; k < nf; ++k) {
                    local[k] = labels[kPrismFaces[f][k]];
                }
                if (!irrep_ising_walker_wang_vertex_admissible(local, nf)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
}

long long
irrep_ising_walker_wang_octahedron_full_count(void)
{
    irrep_ising_object_t labels[12];
    long long total = 1;
    for (int i = 0; i < 12; ++i) total *= 3;
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        long long x = c;
        for (int i = 0; i < 12; ++i) {
            labels[i] = (irrep_ising_object_t)(x % 3);
            x /= 3;
        }
        int ok = 1;
        for (int v = 0; v < 6 && ok; ++v) {
            irrep_ising_object_t local[4] = {
                labels[kOctahedronVertices[v][0]],
                labels[kOctahedronVertices[v][1]],
                labels[kOctahedronVertices[v][2]],
                labels[kOctahedronVertices[v][3]],
            };
            if (!irrep_ising_walker_wang_vertex_admissible(local, 4)) ok = 0;
        }
        if (ok) {
            for (int f = 0; f < 8 && ok; ++f) {
                irrep_ising_object_t local[3] = {
                    labels[kOctahedronFaces[f][0]],
                    labels[kOctahedronFaces[f][1]],
                    labels[kOctahedronFaces[f][2]],
                };
                if (!irrep_ising_walker_wang_vertex_admissible(local, 3)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
}

long long
irrep_ising_walker_wang_cube_full_count(void)
{
    /* 3^12 enumeration; pruning on vertices is applied early to cut
     * the inner constraint loop. */
    irrep_ising_object_t labels[12];
    long long total = 1;
    for (int i = 0; i < 12; ++i) total *= 3;
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        long long x = c;
        for (int i = 0; i < 12; ++i) {
            labels[i] = (irrep_ising_object_t)(x % 3);
            x /= 3;
        }
        int ok = 1;
        /* Vertex constraints first (cheaper; cut more aggressively). */
        for (int v = 0; v < 8 && ok; ++v) {
            irrep_ising_object_t local[3] = {
                labels[kCubeVertices[v][0]],
                labels[kCubeVertices[v][1]],
                labels[kCubeVertices[v][2]],
            };
            if (!irrep_ising_walker_wang_vertex_admissible(local, 3)) ok = 0;
        }
        if (ok) {
            for (int f = 0; f < 6 && ok; ++f) {
                irrep_ising_object_t local[4] = {
                    labels[kCubeFaces[f][0]],
                    labels[kCubeFaces[f][1]],
                    labels[kCubeFaces[f][2]],
                    labels[kCubeFaces[f][3]],
                };
                if (!irrep_ising_walker_wang_vertex_admissible(local, 4)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
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
