/* SPDX-License-Identifier: MIT */
/** @file fibonacci_mtc.c
 *  @brief Implementation of the Fibonacci anyon MTC. */
#include <irrep/fibonacci_mtc.h>

#include <complex.h>
#include <math.h>

static const double PHI = 1.6180339887498948482045868343656;  /* (1 + √5)/2 */

double
irrep_fib_golden_ratio(void)
{
    return PHI;
}

double
irrep_fib_quantum_dim(irrep_fib_object_t a)
{
    return (a == IRREP_FIB_OBJ_TAU) ? PHI : 1.0;
}

double
irrep_fib_global_dim(void)
{
    /* D = √(1 + φ²) = √(2 + φ).  1 + φ² = 1 + φ + 1 = 2 + φ. */
    return sqrt(2.0 + PHI);
}

double
irrep_fib_central_charge(void)
{
    return 14.0 / 5.0;
}

int
irrep_fib_fusion(irrep_fib_object_t a, irrep_fib_object_t b,
                 irrep_fib_object_t c)
{
    /* 1×1 = 1; 1×τ = τ; τ×1 = τ; τ×τ = 1 + τ. */
    if (a == IRREP_FIB_OBJ_1 && b == IRREP_FIB_OBJ_1) return (c == IRREP_FIB_OBJ_1) ? 1 : 0;
    if (a == IRREP_FIB_OBJ_1 && b == IRREP_FIB_OBJ_TAU) return (c == IRREP_FIB_OBJ_TAU) ? 1 : 0;
    if (a == IRREP_FIB_OBJ_TAU && b == IRREP_FIB_OBJ_1) return (c == IRREP_FIB_OBJ_TAU) ? 1 : 0;
    if (a == IRREP_FIB_OBJ_TAU && b == IRREP_FIB_OBJ_TAU) {
        return (c == IRREP_FIB_OBJ_1 || c == IRREP_FIB_OBJ_TAU) ? 1 : 0;
    }
    return 0;
}

double _Complex
irrep_fib_S_matrix(irrep_fib_object_t a, irrep_fib_object_t b)
{
    /* S = (1/D) [[1, φ], [φ, -1]]. */
    double D = irrep_fib_global_dim();
    double entry;
    if (a == IRREP_FIB_OBJ_1 && b == IRREP_FIB_OBJ_1)         entry = 1.0;
    else if (a == IRREP_FIB_OBJ_1 && b == IRREP_FIB_OBJ_TAU)  entry = PHI;
    else if (a == IRREP_FIB_OBJ_TAU && b == IRREP_FIB_OBJ_1)  entry = PHI;
    else                                                       entry = -1.0;
    return entry / D;
}

double _Complex
irrep_fib_T_eigenvalue(irrep_fib_object_t a)
{
    /* T_τ = exp(2πi · h_τ) with h_τ = 2/5 ⇒ T_τ = exp(4πi/5). */
    if (a == IRREP_FIB_OBJ_TAU) return cexp(I * 4.0 * M_PI / 5.0);
    return 1.0;
}

double _Complex
irrep_fib_F_symbol(irrep_fib_object_t a, irrep_fib_object_t b,
                   irrep_fib_object_t c, irrep_fib_object_t d,
                   irrep_fib_object_t e, irrep_fib_object_t f)
{
    /* All sub-fusions must be allowed. */
    if (irrep_fib_fusion(a, b, e) == 0) return 0.0;
    if (irrep_fib_fusion(e, c, d) == 0) return 0.0;
    if (irrep_fib_fusion(b, c, f) == 0) return 0.0;
    if (irrep_fib_fusion(a, f, d) == 0) return 0.0;

    /* Non-trivial Hadamard-like 2×2: F^{τττ}_τ with (e, f) ∈ {1, τ}².
     *
     *   F^{τττ}_τ = [[ 1/φ,    1/√φ ],
     *                [ 1/√φ,  -1/φ  ]]
     */
    if (a == IRREP_FIB_OBJ_TAU && b == IRREP_FIB_OBJ_TAU &&
        c == IRREP_FIB_OBJ_TAU && d == IRREP_FIB_OBJ_TAU) {
        double inv_phi = 1.0 / PHI;
        double inv_sqrt_phi = 1.0 / sqrt(PHI);
        if (e == IRREP_FIB_OBJ_1   && f == IRREP_FIB_OBJ_1)   return  inv_phi;
        if (e == IRREP_FIB_OBJ_1   && f == IRREP_FIB_OBJ_TAU) return  inv_sqrt_phi;
        if (e == IRREP_FIB_OBJ_TAU && f == IRREP_FIB_OBJ_1)   return  inv_sqrt_phi;
        if (e == IRREP_FIB_OBJ_TAU && f == IRREP_FIB_OBJ_TAU) return -inv_phi;
    }
    /* All other allowed F-symbols are 1. */
    return 1.0;
}

double _Complex
irrep_fib_R_symbol(irrep_fib_object_t a, irrep_fib_object_t b,
                   irrep_fib_object_t c)
{
    if (irrep_fib_fusion(a, b, c) == 0) return 0.0;
    /* Trivial cases. */
    if (a == IRREP_FIB_OBJ_1 || b == IRREP_FIB_OBJ_1) return 1.0;
    /* τ × τ → 1 or τ. */
    if (a == IRREP_FIB_OBJ_TAU && b == IRREP_FIB_OBJ_TAU) {
        if (c == IRREP_FIB_OBJ_1)   return cexp(-I * 4.0 * M_PI / 5.0);
        if (c == IRREP_FIB_OBJ_TAU) return cexp( I * 3.0 * M_PI / 5.0);
    }
    return 0.0;
}

/* ====================================================================
 * Consistency proofs
 * ==================================================================== */

double
irrep_fib_F_unitarity_residual(void)
{
    /* F·F* unitary on the τττ block. */
    double max_err = 0.0;
    irrep_fib_object_t s = IRREP_FIB_OBJ_TAU;
    irrep_fib_object_t ch[2] = { IRREP_FIB_OBJ_1, IRREP_FIB_OBJ_TAU };
    for (int fi = 0; fi < 2; ++fi) {
        for (int fj = 0; fj < 2; ++fj) {
            double _Complex acc = 0.0;
            for (int ei = 0; ei < 2; ++ei) {
                double _Complex Fi = irrep_fib_F_symbol(s, s, s, s, ch[ei], ch[fi]);
                double _Complex Fj = irrep_fib_F_symbol(s, s, s, s, ch[ei], ch[fj]);
                acc += Fi * conj(Fj);
            }
            double expected = (fi == fj) ? 1.0 : 0.0;
            double err = cabs(acc - expected);
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

double
irrep_fib_S_squared_residual(void)
{
    double max_err = 0.0;
    for (int a = 0; a < IRREP_FIB_N_OBJECTS; ++a) {
        for (int b = 0; b < IRREP_FIB_N_OBJECTS; ++b) {
            double _Complex acc = 0.0;
            for (int x = 0; x < IRREP_FIB_N_OBJECTS; ++x) {
                acc += irrep_fib_S_matrix((irrep_fib_object_t)a, (irrep_fib_object_t)x)
                     * irrep_fib_S_matrix((irrep_fib_object_t)x, (irrep_fib_object_t)b);
            }
            double expected = (a == b) ? 1.0 : 0.0;
            double err = cabs(acc - expected);
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

double
irrep_fib_verlinde_residual(void)
{
    double max_err = 0.0;
    for (int a = 0; a < IRREP_FIB_N_OBJECTS; ++a) {
        for (int b = 0; b < IRREP_FIB_N_OBJECTS; ++b) {
            for (int c = 0; c < IRREP_FIB_N_OBJECTS; ++c) {
                irrep_fib_object_t A = (irrep_fib_object_t)a;
                irrep_fib_object_t B = (irrep_fib_object_t)b;
                irrep_fib_object_t C = (irrep_fib_object_t)c;
                double _Complex acc = 0.0;
                for (int x = 0; x < IRREP_FIB_N_OBJECTS; ++x) {
                    irrep_fib_object_t X = (irrep_fib_object_t)x;
                    double _Complex S_ax = irrep_fib_S_matrix(A, X);
                    double _Complex S_bx = irrep_fib_S_matrix(B, X);
                    double _Complex S_cx = irrep_fib_S_matrix(C, X);
                    double _Complex S_0x = irrep_fib_S_matrix(IRREP_FIB_OBJ_1, X);
                    if (cabs(S_0x) < 1e-15) continue;
                    acc += S_ax * S_bx * conj(S_cx) / S_0x;
                }
                int N_actual = irrep_fib_fusion(A, B, C);
                double err = cabs(acc - (double)N_actual);
                if (err > max_err) max_err = err;
            }
        }
    }
    return max_err;
}

double
irrep_fib_twist_from_R_residual(void)
{
    /* θ_a = (1/d_a) Σ_c N^{aa}_c d_c R^{aa}_c (self-dual MTC). */
    double max_err = 0.0;
    for (int a = 0; a < IRREP_FIB_N_OBJECTS; ++a) {
        irrep_fib_object_t A = (irrep_fib_object_t)a;
        double d_a = irrep_fib_quantum_dim(A);
        double _Complex theta = 0.0;
        for (int c = 0; c < IRREP_FIB_N_OBJECTS; ++c) {
            irrep_fib_object_t C = (irrep_fib_object_t)c;
            int N = irrep_fib_fusion(A, A, C);
            if (N == 0) continue;
            double d_c = irrep_fib_quantum_dim(C);
            double _Complex R = irrep_fib_R_symbol(A, A, C);
            theta += (double)N * d_c * R;
        }
        theta /= d_a;
        double _Complex T_hardcoded = irrep_fib_T_eigenvalue(A);
        double err = cabs(theta - T_hardcoded);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

/* ====================================================================
 * Walker-Wang admissibility for Fibonacci
 * ==================================================================== */

int
irrep_fib_admissible(const irrep_fib_object_t *labels, int n)
{
    if (labels == NULL || n < 0) return 0;
    int n_tau = 0;
    for (int i = 0; i < n; ++i) {
        if (labels[i] == IRREP_FIB_OBJ_TAU) ++n_tau;
    }
    return (n_tau != 1) ? 1 : 0;
}

/* 3-simplex / cube / octahedron / bipyramid / prism incidences —
 * identical geometry to the Ising/Z₂×Z₂ modules. */
static const int kFibSimplex3Vtx[4][3] = {
    {0, 1, 2}, {0, 3, 4}, {1, 3, 5}, {2, 4, 5},
};
static const int kFibSimplex3Faces[4][3] = {
    {0, 1, 3}, {0, 2, 4}, {1, 2, 5}, {3, 4, 5},
};

long long
irrep_fib_walker_wang_simplex3_full_count(void)
{
    irrep_fib_object_t labels[6];
    long long total = 1LL << 6;  /* 2^6 = 64 */
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        for (int i = 0; i < 6; ++i) {
            labels[i] = ((c >> i) & 1) ? IRREP_FIB_OBJ_TAU : IRREP_FIB_OBJ_1;
        }
        int ok = 1;
        for (int v = 0; v < 4 && ok; ++v) {
            irrep_fib_object_t local[3] = {
                labels[kFibSimplex3Vtx[v][0]],
                labels[kFibSimplex3Vtx[v][1]],
                labels[kFibSimplex3Vtx[v][2]],
            };
            if (!irrep_fib_admissible(local, 3)) ok = 0;
        }
        if (ok) {
            for (int f = 0; f < 4 && ok; ++f) {
                irrep_fib_object_t local[3] = {
                    labels[kFibSimplex3Faces[f][0]],
                    labels[kFibSimplex3Faces[f][1]],
                    labels[kFibSimplex3Faces[f][2]],
                };
                if (!irrep_fib_admissible(local, 3)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
}

static const int kFibCubeVtx[8][3] = {
    {0, 4,  8}, {0, 5,  9}, {1, 4, 10}, {1, 5, 11},
    {2, 6,  8}, {2, 7,  9}, {3, 6, 10}, {3, 7, 11},
};
static const int kFibCubeFaces[6][4] = {
    {0, 4, 1, 5}, {2, 6, 3, 7},
    {0, 8, 2, 9}, {1, 10, 3, 11},
    {4, 8, 6, 10}, {5, 9, 7, 11},
};

long long
irrep_fib_walker_wang_cube_full_count(void)
{
    irrep_fib_object_t labels[12];
    long long total = 1LL << 12;
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        for (int i = 0; i < 12; ++i) {
            labels[i] = ((c >> i) & 1) ? IRREP_FIB_OBJ_TAU : IRREP_FIB_OBJ_1;
        }
        int ok = 1;
        for (int v = 0; v < 8 && ok; ++v) {
            irrep_fib_object_t local[3] = {
                labels[kFibCubeVtx[v][0]],
                labels[kFibCubeVtx[v][1]],
                labels[kFibCubeVtx[v][2]],
            };
            if (!irrep_fib_admissible(local, 3)) ok = 0;
        }
        if (ok) {
            for (int f = 0; f < 6 && ok; ++f) {
                irrep_fib_object_t local[4] = {
                    labels[kFibCubeFaces[f][0]],
                    labels[kFibCubeFaces[f][1]],
                    labels[kFibCubeFaces[f][2]],
                    labels[kFibCubeFaces[f][3]],
                };
                if (!irrep_fib_admissible(local, 4)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
}

static const int kFibOctVtx[6][4] = {
    {  0, 1, 2,  3 }, {  4, 5, 6,  7 },
    {  0, 4, 8, 11 }, {  1, 5, 8,  9 },
    {  2, 6, 9, 10 }, {  3, 7, 10, 11 },
};
static const int kFibOctFaces[8][3] = {
    { 0, 1,  8 }, { 1, 2,  9 }, { 2, 3, 10 }, { 3, 0, 11 },
    { 4, 5,  8 }, { 5, 6,  9 }, { 6, 7, 10 }, { 7, 4, 11 },
};

long long
irrep_fib_walker_wang_octahedron_full_count(void)
{
    irrep_fib_object_t labels[12];
    long long total = 1LL << 12;
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        for (int i = 0; i < 12; ++i) {
            labels[i] = ((c >> i) & 1) ? IRREP_FIB_OBJ_TAU : IRREP_FIB_OBJ_1;
        }
        int ok = 1;
        for (int v = 0; v < 6 && ok; ++v) {
            irrep_fib_object_t local[4] = {
                labels[kFibOctVtx[v][0]],
                labels[kFibOctVtx[v][1]],
                labels[kFibOctVtx[v][2]],
                labels[kFibOctVtx[v][3]],
            };
            if (!irrep_fib_admissible(local, 4)) ok = 0;
        }
        if (ok) {
            for (int f = 0; f < 8 && ok; ++f) {
                irrep_fib_object_t local[3] = {
                    labels[kFibOctFaces[f][0]],
                    labels[kFibOctFaces[f][1]],
                    labels[kFibOctFaces[f][2]],
                };
                if (!irrep_fib_admissible(local, 3)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
}

static const int kFibBipyramidVtxSizes[5] = { 3, 3, 4, 4, 4 };
static const int kFibBipyramidVtx[5][4] = {
    {  0,  1,  2, -1 }, {  3,  4,  5, -1 },
    {  0,  3,  6,  8 }, {  1,  4,  6,  7 }, {  2,  5,  7,  8 },
};
static const int kFibBipyramidFaces[6][3] = {
    { 0, 1, 6 }, { 1, 2, 7 }, { 2, 0, 8 },
    { 3, 4, 6 }, { 4, 5, 7 }, { 5, 3, 8 },
};

long long
irrep_fib_walker_wang_tri_bipyramid_full_count(void)
{
    irrep_fib_object_t labels[9];
    long long total = 1LL << 9;
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        for (int i = 0; i < 9; ++i) {
            labels[i] = ((c >> i) & 1) ? IRREP_FIB_OBJ_TAU : IRREP_FIB_OBJ_1;
        }
        int ok = 1;
        for (int v = 0; v < 5 && ok; ++v) {
            int nv = kFibBipyramidVtxSizes[v];
            irrep_fib_object_t local[4];
            for (int k = 0; k < nv; ++k) local[k] = labels[kFibBipyramidVtx[v][k]];
            if (!irrep_fib_admissible(local, nv)) ok = 0;
        }
        if (ok) {
            for (int f = 0; f < 6 && ok; ++f) {
                irrep_fib_object_t local[3] = {
                    labels[kFibBipyramidFaces[f][0]],
                    labels[kFibBipyramidFaces[f][1]],
                    labels[kFibBipyramidFaces[f][2]],
                };
                if (!irrep_fib_admissible(local, 3)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
}

static const int kFibPrismVtx[6][3] = {
    { 0, 2, 6 }, { 0, 1, 7 }, { 1, 2, 8 },
    { 3, 5, 6 }, { 3, 4, 7 }, { 4, 5, 8 },
};
static const int kFibPrismFaceSizes[5] = { 3, 3, 4, 4, 4 };
static const int kFibPrismFaces[5][4] = {
    { 0, 1, 2, -1 }, { 3, 4, 5, -1 },
    { 0, 7, 3,  6 }, { 1, 8, 4,  7 }, { 2, 6, 5,  8 },
};

long long
irrep_fib_walker_wang_tri_prism_full_count(void)
{
    irrep_fib_object_t labels[9];
    long long total = 1LL << 9;
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        for (int i = 0; i < 9; ++i) {
            labels[i] = ((c >> i) & 1) ? IRREP_FIB_OBJ_TAU : IRREP_FIB_OBJ_1;
        }
        int ok = 1;
        for (int v = 0; v < 6 && ok; ++v) {
            irrep_fib_object_t local[3] = {
                labels[kFibPrismVtx[v][0]],
                labels[kFibPrismVtx[v][1]],
                labels[kFibPrismVtx[v][2]],
            };
            if (!irrep_fib_admissible(local, 3)) ok = 0;
        }
        if (ok) {
            for (int f = 0; f < 5 && ok; ++f) {
                int nf = kFibPrismFaceSizes[f];
                irrep_fib_object_t local[4];
                for (int k = 0; k < nf; ++k) local[k] = labels[kFibPrismFaces[f][k]];
                if (!irrep_fib_admissible(local, nf)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
}

/* ====================================================================
 * Anyonic braiding on 4 Fibonacci anyons.
 *
 * Fusion space basis: |1⟩, |τ⟩ — the two channels for τ × τ → ?
 * fusion of the first pair (1, 2). After fixing the total fusion to
 * vacuum, the last pair's channel is determined.
 *
 * Braid generators σ_1 and σ_3 are diagonal in this basis; σ_2 acts as
 * F · σ_1 · F⁻¹ where F is the τττ F-matrix.
 * ==================================================================== */

void
irrep_fib_braid_sigma_1(double _Complex state[2], int inverse)
{
    double _Complex R1 = cexp(-I * 4.0 * M_PI / 5.0);
    double _Complex Rt = cexp( I * 3.0 * M_PI / 5.0);
    if (inverse) {
        R1 = conj(R1);
        Rt = conj(Rt);
    }
    state[0] *= R1;
    state[1] *= Rt;
}

void
irrep_fib_braid_sigma_2(double _Complex state[2], int inverse)
{
    /* F^{τττ}_τ matrix (symmetric, real, unitary). */
    double inv_phi = 1.0 / PHI;
    double inv_sqrt_phi = 1.0 / sqrt(PHI);
    double F00 =  inv_phi, F01 =  inv_sqrt_phi;
    double F10 =  inv_sqrt_phi, F11 = -inv_phi;
    /* F⁻¹ = F^T = F (symmetric). */

    /* (1) F⁻¹ · state. */
    double _Complex s0 = F00 * state[0] + F10 * state[1];
    double _Complex s1 = F01 * state[0] + F11 * state[1];
    /* (2) diag(R1, Rt) · F⁻¹ · state. */
    double _Complex R1 = cexp(-I * 4.0 * M_PI / 5.0);
    double _Complex Rt = cexp( I * 3.0 * M_PI / 5.0);
    if (inverse) { R1 = conj(R1); Rt = conj(Rt); }
    s0 *= R1;
    s1 *= Rt;
    /* (3) F · (diag · F⁻¹ · state). */
    state[0] = F00 * s0 + F01 * s1;
    state[1] = F10 * s0 + F11 * s1;
}

/* Build the 2×2 matrix form of a generator by applying it to e_0, e_1. */
static void
build_matrix(void (*gen)(double _Complex[2], int), int inv,
             double _Complex M[2][2])
{
    double _Complex e0[2] = { 1.0, 0.0 };
    double _Complex e1[2] = { 0.0, 1.0 };
    gen(e0, inv); gen(e1, inv);
    M[0][0] = e0[0]; M[1][0] = e0[1];
    M[0][1] = e1[0]; M[1][1] = e1[1];
}

/* 2×2 matrix multiply: out = A · B. */
static void
mat2_mul(double _Complex A[2][2], double _Complex B[2][2],
         double _Complex out[2][2])
{
    out[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
    out[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
    out[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
    out[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
}

double
irrep_fib_braid_unitarity_residual(void)
{
    double max_err = 0.0;
    double _Complex M[2][2], M_dag[2][2], prod[2][2];
    /* σ_1 unitarity. */
    build_matrix(irrep_fib_braid_sigma_1, 0, M);
    M_dag[0][0] = conj(M[0][0]); M_dag[0][1] = conj(M[1][0]);
    M_dag[1][0] = conj(M[0][1]); M_dag[1][1] = conj(M[1][1]);
    mat2_mul(M, M_dag, prod);
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
        double expected = (i == j) ? 1.0 : 0.0;
        double err = cabs(prod[i][j] - expected);
        if (err > max_err) max_err = err;
    }
    /* σ_2 unitarity. */
    build_matrix(irrep_fib_braid_sigma_2, 0, M);
    M_dag[0][0] = conj(M[0][0]); M_dag[0][1] = conj(M[1][0]);
    M_dag[1][0] = conj(M[0][1]); M_dag[1][1] = conj(M[1][1]);
    mat2_mul(M, M_dag, prod);
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
        double expected = (i == j) ? 1.0 : 0.0;
        double err = cabs(prod[i][j] - expected);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

double
irrep_fib_braid_yang_baxter_residual(void)
{
    /* σ_1 · σ_2 · σ_1 = σ_2 · σ_1 · σ_2. */
    double _Complex S1[2][2], S2[2][2];
    double _Complex tmp1[2][2], LHS[2][2], tmp2[2][2], RHS[2][2];
    build_matrix(irrep_fib_braid_sigma_1, 0, S1);
    build_matrix(irrep_fib_braid_sigma_2, 0, S2);
    mat2_mul(S1, S2, tmp1);
    mat2_mul(tmp1, S1, LHS);
    mat2_mul(S2, S1, tmp2);
    mat2_mul(tmp2, S2, RHS);
    double max_err = 0.0;
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
        double err = cabs(LHS[i][j] - RHS[i][j]);
        if (err > max_err) max_err = err;
    }
    return max_err;
}
