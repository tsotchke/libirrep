/* SPDX-License-Identifier: MIT */
/** @file toric_mtc.c
 *  @brief Implementation of the Z₂ × Z₂ toric-code MTC. */
#include <irrep/toric_mtc.h>

#include <complex.h>
#include <math.h>

/* Z₂ × Z₂ fusion table:
 *
 *           1       e       m       ψ
 *      1 |  1       e       m       ψ
 *      e |  e       1       ψ       m
 *      m |  m       ψ       1       e
 *      ψ |  ψ       m       e       1
 *
 * Indexed as fusion_table[a][b] = a × b. */
static const irrep_toric_mtc_object_t kFusionTable[4][4] = {
    { IRREP_TORIC_MTC_OBJ_1,   IRREP_TORIC_MTC_OBJ_E,
      IRREP_TORIC_MTC_OBJ_M,   IRREP_TORIC_MTC_OBJ_PSI },
    { IRREP_TORIC_MTC_OBJ_E,   IRREP_TORIC_MTC_OBJ_1,
      IRREP_TORIC_MTC_OBJ_PSI, IRREP_TORIC_MTC_OBJ_M   },
    { IRREP_TORIC_MTC_OBJ_M,   IRREP_TORIC_MTC_OBJ_PSI,
      IRREP_TORIC_MTC_OBJ_1,   IRREP_TORIC_MTC_OBJ_E   },
    { IRREP_TORIC_MTC_OBJ_PSI, IRREP_TORIC_MTC_OBJ_M,
      IRREP_TORIC_MTC_OBJ_E,   IRREP_TORIC_MTC_OBJ_1   },
};

double
irrep_toric_mtc_quantum_dim(irrep_toric_mtc_object_t a)
{
    (void)a;
    return 1.0; /* All simples of an abelian MTC have d = 1. */
}

double
irrep_toric_mtc_global_dim(void)
{
    return 2.0; /* D = √(4 · 1²) = 2. */
}

double
irrep_toric_mtc_central_charge(void)
{
    return 0.0; /* Non-chiral / doubled topological order. */
}

int
irrep_toric_mtc_fusion(irrep_toric_mtc_object_t a,
                       irrep_toric_mtc_object_t b,
                       irrep_toric_mtc_object_t c)
{
    if (a < 0 || a >= IRREP_TORIC_MTC_N_OBJECTS) return 0;
    if (b < 0 || b >= IRREP_TORIC_MTC_N_OBJECTS) return 0;
    if (c < 0 || c >= IRREP_TORIC_MTC_N_OBJECTS) return 0;
    return (kFusionTable[a][b] == c) ? 1 : 0;
}

double _Complex
irrep_toric_mtc_S_matrix(irrep_toric_mtc_object_t a,
                         irrep_toric_mtc_object_t b)
{
    /* S = (1/2) × character table of Z₂×Z₂. Entry sign computed from
     * the symplectic bilinear form on the Z₂×Z₂ lattice:
     *
     *   Identify e = (1, 0), m = (0, 1), ψ = (1, 1) in (Z_2)². The
     *   "symplectic" form is the standard mutual-statistics form:
     *   <(a_e, a_m), (b_e, b_m)> = a_e · b_m + a_m · b_e (mod 2).
     *
     *   S_{ab} = (1/D) · (-1)^{<a, b>} = (1/2) · (-1)^{a_e b_m + a_m b_e}.
     *
     * Lookups: a_e = (a == E or a == PSI), a_m = (a == M or a == PSI). */
    int a_e = (a == IRREP_TORIC_MTC_OBJ_E   || a == IRREP_TORIC_MTC_OBJ_PSI) ? 1 : 0;
    int a_m = (a == IRREP_TORIC_MTC_OBJ_M   || a == IRREP_TORIC_MTC_OBJ_PSI) ? 1 : 0;
    int b_e = (b == IRREP_TORIC_MTC_OBJ_E   || b == IRREP_TORIC_MTC_OBJ_PSI) ? 1 : 0;
    int b_m = (b == IRREP_TORIC_MTC_OBJ_M   || b == IRREP_TORIC_MTC_OBJ_PSI) ? 1 : 0;
    int sign_pow = (a_e * b_m + a_m * b_e) & 1;
    double sign = sign_pow ? -1.0 : +1.0;
    return 0.5 * sign;
}

double _Complex
irrep_toric_mtc_T_eigenvalue(irrep_toric_mtc_object_t a)
{
    return (a == IRREP_TORIC_MTC_OBJ_PSI) ? -1.0 : 1.0;
}

double _Complex
irrep_toric_mtc_F_symbol(irrep_toric_mtc_object_t a, irrep_toric_mtc_object_t b,
                         irrep_toric_mtc_object_t c, irrep_toric_mtc_object_t d,
                         irrep_toric_mtc_object_t e, irrep_toric_mtc_object_t f)
{
    /* Abelian MTC: F = 1 if all four sub-fusions are allowed, 0 else.
     * (For the trivial 2-cocycle representative; non-trivial cocycle
     * would give signs but is gauge-equivalent here.) */
    if (irrep_toric_mtc_fusion(a, b, e) == 0) return 0.0;
    if (irrep_toric_mtc_fusion(e, c, d) == 0) return 0.0;
    if (irrep_toric_mtc_fusion(b, c, f) == 0) return 0.0;
    if (irrep_toric_mtc_fusion(a, f, d) == 0) return 0.0;
    return 1.0;
}

double _Complex
irrep_toric_mtc_R_symbol(irrep_toric_mtc_object_t a,
                         irrep_toric_mtc_object_t b,
                         irrep_toric_mtc_object_t c)
{
    if (irrep_toric_mtc_fusion(a, b, c) == 0) return 0.0;

    /* Trivial: either factor is identity. */
    if (a == IRREP_TORIC_MTC_OBJ_1 || b == IRREP_TORIC_MTC_OBJ_1) {
        return 1.0;
    }
    /* Self-braids. */
    if (a == b) {
        if (a == IRREP_TORIC_MTC_OBJ_E)   return 1.0;
        if (a == IRREP_TORIC_MTC_OBJ_M)   return 1.0;
        if (a == IRREP_TORIC_MTC_OBJ_PSI) return -1.0;
    }
    /* Mutual e↔m (and derived ψ braids).
     *
     * Convention: half-monodromy R^{em}_ψ = R^{me}_ψ = i. Then the
     * mutual-statistics phase is R^{em}_ψ · R^{me}_ψ = i² = -1. */
    if ((a == IRREP_TORIC_MTC_OBJ_E && b == IRREP_TORIC_MTC_OBJ_M)   ||
        (a == IRREP_TORIC_MTC_OBJ_M && b == IRREP_TORIC_MTC_OBJ_E)   ||
        (a == IRREP_TORIC_MTC_OBJ_E && b == IRREP_TORIC_MTC_OBJ_PSI) ||
        (a == IRREP_TORIC_MTC_OBJ_PSI && b == IRREP_TORIC_MTC_OBJ_E) ||
        (a == IRREP_TORIC_MTC_OBJ_M && b == IRREP_TORIC_MTC_OBJ_PSI) ||
        (a == IRREP_TORIC_MTC_OBJ_PSI && b == IRREP_TORIC_MTC_OBJ_M)) {
        return I;
    }
    return 0.0; /* unreachable */
}

/* ====================================================================
 * Consistency proofs
 * ==================================================================== */

double
irrep_toric_mtc_S_squared_residual(void)
{
    double max_err = 0.0;
    for (int a = 0; a < IRREP_TORIC_MTC_N_OBJECTS; ++a) {
        for (int b = 0; b < IRREP_TORIC_MTC_N_OBJECTS; ++b) {
            double _Complex acc = 0.0;
            for (int x = 0; x < IRREP_TORIC_MTC_N_OBJECTS; ++x) {
                acc += irrep_toric_mtc_S_matrix((irrep_toric_mtc_object_t)a,
                                                (irrep_toric_mtc_object_t)x)
                     * irrep_toric_mtc_S_matrix((irrep_toric_mtc_object_t)x,
                                                (irrep_toric_mtc_object_t)b);
            }
            double expected = (a == b) ? 1.0 : 0.0;
            double err = cabs(acc - expected);
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

double
irrep_toric_mtc_verlinde_residual(void)
{
    double max_err = 0.0;
    for (int a = 0; a < IRREP_TORIC_MTC_N_OBJECTS; ++a) {
        for (int b = 0; b < IRREP_TORIC_MTC_N_OBJECTS; ++b) {
            for (int c = 0; c < IRREP_TORIC_MTC_N_OBJECTS; ++c) {
                double _Complex acc = 0.0;
                for (int x = 0; x < IRREP_TORIC_MTC_N_OBJECTS; ++x) {
                    double _Complex S_ax = irrep_toric_mtc_S_matrix(
                        (irrep_toric_mtc_object_t)a, (irrep_toric_mtc_object_t)x);
                    double _Complex S_bx = irrep_toric_mtc_S_matrix(
                        (irrep_toric_mtc_object_t)b, (irrep_toric_mtc_object_t)x);
                    double _Complex S_cx = irrep_toric_mtc_S_matrix(
                        (irrep_toric_mtc_object_t)c, (irrep_toric_mtc_object_t)x);
                    double _Complex S_0x = irrep_toric_mtc_S_matrix(
                        IRREP_TORIC_MTC_OBJ_1, (irrep_toric_mtc_object_t)x);
                    if (cabs(S_0x) < 1e-15) continue;
                    acc += S_ax * S_bx * conj(S_cx) / S_0x;
                }
                int N_actual = irrep_toric_mtc_fusion(
                    (irrep_toric_mtc_object_t)a, (irrep_toric_mtc_object_t)b,
                    (irrep_toric_mtc_object_t)c);
                double err = cabs(acc - (double)N_actual);
                if (err > max_err) max_err = err;
            }
        }
    }
    return max_err;
}

double
irrep_toric_mtc_twist_from_R_residual(void)
{
    double max_err = 0.0;
    for (int a = 0; a < IRREP_TORIC_MTC_N_OBJECTS; ++a) {
        irrep_toric_mtc_object_t A = (irrep_toric_mtc_object_t)a;
        /* θ_a = R^{aa}_1 for abelian self-dual MTC (single channel). */
        double _Complex theta = irrep_toric_mtc_R_symbol(A, A, IRREP_TORIC_MTC_OBJ_1);
        double _Complex T_hardcoded = irrep_toric_mtc_T_eigenvalue(A);
        double err = cabs(theta - T_hardcoded);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

/* ====================================================================
 * Walker-Wang on the 3-simplex (toric MTC version).
 * ==================================================================== */

int
irrep_toric_mtc_admissible(const irrep_toric_mtc_object_t *labels, int n)
{
    if (labels == NULL || n < 0) return 0;
    if (n == 0) return 1;
    irrep_toric_mtc_object_t prod = labels[0];
    for (int i = 1; i < n; ++i) {
        prod = kFusionTable[prod][labels[i]];
    }
    return prod == IRREP_TORIC_MTC_OBJ_1 ? 1 : 0;
}

/* Re-use the tetrahedron edge incidences from the Ising module — same
 * 3-simplex geometry, different MTC. */
static const int kSimplex3VtxToric[4][3] = {
    {0, 1, 2}, {0, 3, 4}, {1, 3, 5}, {2, 4, 5},
};
static const int kSimplex3FaceToric[4][3] = {
    {0, 1, 3}, {0, 2, 4}, {1, 2, 5}, {3, 4, 5},
};

long long
irrep_toric_mtc_walker_wang_simplex3_full_count(void)
{
    irrep_toric_mtc_object_t labels[6];
    long long total = 1;
    for (int i = 0; i < 6; ++i) total *= 4;  /* 4^6 = 4096 */
    long long count = 0;
    for (long long c = 0; c < total; ++c) {
        long long x = c;
        for (int i = 0; i < 6; ++i) {
            labels[i] = (irrep_toric_mtc_object_t)(x % 4);
            x /= 4;
        }
        int ok = 1;
        for (int v = 0; v < 4 && ok; ++v) {
            irrep_toric_mtc_object_t local[3] = {
                labels[kSimplex3VtxToric[v][0]],
                labels[kSimplex3VtxToric[v][1]],
                labels[kSimplex3VtxToric[v][2]],
            };
            if (!irrep_toric_mtc_admissible(local, 3)) ok = 0;
        }
        if (ok) {
            for (int f = 0; f < 4 && ok; ++f) {
                irrep_toric_mtc_object_t local[3] = {
                    labels[kSimplex3FaceToric[f][0]],
                    labels[kSimplex3FaceToric[f][1]],
                    labels[kSimplex3FaceToric[f][2]],
                };
                if (!irrep_toric_mtc_admissible(local, 3)) ok = 0;
            }
        }
        if (ok) ++count;
    }
    return count;
}
