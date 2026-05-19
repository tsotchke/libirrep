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
