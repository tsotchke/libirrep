/* SPDX-License-Identifier: MIT */
/* M4: Clebsch-Gordan coefficients and Wigner 3j symbols.
 *
 * Primary kernel: Schulten–Gordon three-term recurrence for the Wigner 3j
 * symbol over varying j₁ at fixed (j₂, j₃, m₂, m₃) (Schulten & Gordon, J.
 * Math. Phys. 16, 1961, 1975; Luscombe & Luban, Phys. Rev. E 57, 7274,
 * 1998). Forward recurrence from j_min = max(|j₂−j₃|, |m₁|) to
 * j_max = j₂+j₃ with arbitrary initial value, then normalise via the sum
 * rule Σ_j (2j+1) · 3j(j, j₂, j₃; m₁, m₂, m₃)² = 1 and fix the overall
 * sign from the closed-form at j_max, (−1)^{j₂−j₃−m₁}.
 *
 * This replaces the earlier Racah-log-gamma single-sum form, which lost
 * precision to catastrophic cancellation in the alternating-sign sum past
 * j ≈ 20 and produced NaN past j ≈ 60 near the triangle edge.
 *
 * Clebsch-Gordan is derived from the 3j symbol via
 *   ⟨j₁ m₁; j₂ m₂ | J M⟩ = (−1)^{j₁−j₂+M} · √(2J+1) · (j₁ j₂ J; m₁ m₂ −M).
 *
 * Doubled-integer arguments (`_2j` variants) are primary; integer wrappers
 * multiply by 2. See also Varshalovich §8; Edmonds §3.
 */

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include <irrep/clebsch_gordan.h>

#define UNUSED(x) ((void)(x))

static inline int iabs_(int x) { return x < 0 ? -x : x; }

/* -------------------------------------------------------------------------- *
 * Schulten-Gordon forward recurrence for the Wigner-3j series.               *
 *                                                                            *
 * For fixed (j₂, j₃, m₁, m₂, m₃) with m₁ + m₂ + m₃ = 0, the 3j symbol        *
 * T(j) = (j, j₂, j₃; m₁, m₂, m₃) satisfies                                   *
 *                                                                            *
 *   j · E(j+1) · T(j+1) + F(j) · T(j) + (j+1) · E(j) · T(j−1) = 0            *
 *                                                                            *
 * where                                                                      *
 *   E(j)² = [j² − (j₂−j₃)²] · [(j₂+j₃+1)² − j²] · [j² − m₁²]                 *
 *   F(j)  = −(2j+1) · { [j₂(j₂+1) − j₃(j₃+1)] · m₁ + j(j+1) · (m₂ − m₃) }    *
 *                                                                            *
 * Forward recurrence (starting with T(j_min) arbitrary, T(j_min − 1) = 0)    *
 * is stable in the classical region; the arbitrary overall scale and sign   *
 * are fixed by the sum rule Σ_j (2j+1) T(j)² = 1 and the sign of T(j_max)   *
 * from the closed-form (−1)^{j₂−j₃−m₁}.                                     *
 * -------------------------------------------------------------------------- */

#define IRREP_3J_MAX_SERIES 512

static double three_j_E_squared_(double j, double j2mj3, double j2pj3p1, double m1) {
    return (j*j - j2mj3*j2mj3)
         * (j2pj3p1*j2pj3p1 - j*j)
         * (j*j - m1*m1);
}

static int wigner_3j_series_(double j2, double j3,
                             double m1, double m2, double m3,
                             double j_min, double j_max,
                             double *T_out, int N) {
    if (N <= 0 || N > IRREP_3J_MAX_SERIES) return -1;

    double j2mj3   = j2 - j3;
    double j2pj3p1 = j2 + j3 + 1.0;

    /* Backward recurrence from j_max down to j_min.  The coefficient
     * (j+1)·E(j) in the T(j-1) position is nonzero for all j > j_min
     * (E(j) > 0 for j ∈ (j_min, j_max]); forward recurrence would hit a
     * 0/0 division when j_min = 0 and m₁ = 0 because the coefficient of
     * T(j+1) is j·E(j+1), zero at j = 0.
     *
     *   T(j-1) = -[ j·E(j+1)·T(j+1) + F(j)·T(j) ] / [ (j+1)·E(j) ]
     *
     * Initialised with T(j_max+1) = 0 (outside support) and T(j_max) = 1
     * (arbitrary; normalised below).  T(j_max - 1) then follows from the
     * recurrence at j = j_max. */
    double T_next = 0.0;                  /* T(j_max + 1) */
    double T_curr = 1.0;                  /* T(j_max); arbitrary */
    T_out[N - 1] = T_curr;

    for (int idx = N - 2; idx >= 0; --idx) {
        double j   = j_min + (double)(idx + 1);   /* recurrence index yielding T(j-1) */
        double jp1 = j + 1.0;

        double Ej_sq   = three_j_E_squared_(j,   j2mj3, j2pj3p1, m1);
        double Ejp1_sq = three_j_E_squared_(jp1, j2mj3, j2pj3p1, m1);
        double Ej   = (Ej_sq   > 0.0) ? sqrt(Ej_sq)   : 0.0;
        double Ejp1 = (Ejp1_sq > 0.0) ? sqrt(Ejp1_sq) : 0.0;

        double Fj = -(2.0*j + 1.0)
                  * ( (j2*(j2 + 1.0) - j3*(j3 + 1.0)) * m1
                    + j*(j + 1.0) * (m2 - m3) );

        double T_prev;
        if (Ej > 0.0) {
            T_prev = -(j * Ejp1 * T_next + Fj * T_curr) / (jp1 * Ej);
        } else {
            /* Shouldn't happen inside the loop: Ej > 0 for j ∈ (j_min, j_max]. */
            T_prev = 0.0;
        }

        T_out[idx] = T_prev;
        T_next     = T_curr;
        T_curr     = T_prev;
    }

    /* Sum-rule normalisation: Σ_j (2j+1) · T(j)² = 1. */
    double norm_sq = 0.0;
    for (int idx = 0; idx < N; ++idx) {
        double j = j_min + (double)idx;
        norm_sq += (2.0*j + 1.0) * T_out[idx] * T_out[idx];
    }
    if (norm_sq <= 0.0) {
        for (int idx = 0; idx < N; ++idx) T_out[idx] = 0.0;
        return 0;
    }
    double norm = sqrt(norm_sq);
    for (int idx = 0; idx < N; ++idx) T_out[idx] /= norm;

    /* Sign convention: at j_max = j₂+j₃, the 3j symbol equals
     *   (−1)^{j₂−j₃−m₁} · √(positive factor).
     * j₂−j₃−m₁ is always integer for valid 3j (parity constraints).
     * If recurrence came out with the wrong overall sign, flip. */
    int    exp_int = (int)lround(j2 - j3 - m1);
    double expected_sign = (exp_int & 1) ? -1.0 : 1.0;
    if (T_out[N-1] * expected_sign < 0.0) {
        for (int idx = 0; idx < N; ++idx) T_out[idx] = -T_out[idx];
    }
    return 0;
}

/* -------------------------------------------------------------------------- *
 * Clebsch-Gordan — Racah formula                                             *
 * -------------------------------------------------------------------------- */

double irrep_cg_2j(int two_j1, int two_m1,
                   int two_j2, int two_m2,
                   int two_J,  int two_M) {
    /* CG ↔ 3j relation:
     *   ⟨j₁ m₁ j₂ m₂ | J M⟩ = (−1)^{j₁−j₂+M} · √(2J+1) · (j₁ j₂ J; m₁ m₂ −M)
     *
     * m-conservation check here; the remaining selection rules are enforced
     * inside irrep_wigner_3j_2j. */
    if (two_m1 + two_m2 != two_M) return 0.0;
    double three_j = irrep_wigner_3j_2j(two_j1, two_m1,
                                        two_j2, two_m2,
                                        two_J,  -two_M);
    if (three_j == 0.0) return 0.0;
    int phase_half = (two_j1 - two_j2 + two_M) / 2;
    double phase   = (phase_half & 1) ? -1.0 : 1.0;
    return phase * sqrt((double)two_J + 1.0) * three_j;
}

double irrep_cg(int j1, int m1, int j2, int m2, int J, int M) {
    return irrep_cg_2j(2*j1, 2*m1, 2*j2, 2*m2, 2*J, 2*M);
}

/* -------------------------------------------------------------------------- *
 * Wigner 3j via Schulten-Gordon recurrence.                                  *
 * -------------------------------------------------------------------------- */

double irrep_wigner_3j_2j(int two_j1, int two_m1,
                          int two_j2, int two_m2,
                          int two_j3, int two_m3) {
    /* selection rules */
    if (two_j1 < 0 || two_j2 < 0 || two_j3 < 0) return 0.0;
    if (two_m1 + two_m2 + two_m3 != 0) return 0.0;
    if (two_m1 < -two_j1 || two_m1 > two_j1) return 0.0;
    if (two_m2 < -two_j2 || two_m2 > two_j2) return 0.0;
    if (two_m3 < -two_j3 || two_m3 > two_j3) return 0.0;
    if ((two_j1 + two_m1) & 1) return 0.0;
    if ((two_j2 + two_m2) & 1) return 0.0;
    if ((two_j3 + two_m3) & 1) return 0.0;
    if ((two_j1 + two_j2 + two_j3) & 1) return 0.0;
    if (two_j1 < iabs_(two_j2 - two_j3) || two_j1 > two_j2 + two_j3) return 0.0;

    int two_j1_min = iabs_(two_j2 - two_j3);
    if (iabs_(two_m1) > two_j1_min) two_j1_min = iabs_(two_m1);
    int two_j1_max = two_j2 + two_j3;

    /* The series steps by 1 in j₁, i.e. by 2 in two_j₁; parity must agree. */
    if ((two_j1 - two_j1_min) & 1) return 0.0;
    int N   = (two_j1_max - two_j1_min) / 2 + 1;
    int idx = (two_j1     - two_j1_min) / 2;
    if (idx < 0 || idx >= N) return 0.0;

    double j2 = 0.5 * (double)two_j2;
    double j3 = 0.5 * (double)two_j3;
    double m1 = 0.5 * (double)two_m1;
    double m2 = 0.5 * (double)two_m2;
    double m3 = 0.5 * (double)two_m3;
    double j_min = 0.5 * (double)two_j1_min;
    double j_max = 0.5 * (double)two_j1_max;

    double T[IRREP_3J_MAX_SERIES];
    if (wigner_3j_series_(j2, j3, m1, m2, m3, j_min, j_max, T, N) != 0) return 0.0;
    return T[idx];
}

double irrep_wigner_3j(int j1, int m1, int j2, int m2, int j3, int m3) {
    return irrep_wigner_3j_2j(2*j1, 2*m1, 2*j2, 2*m2, 2*j3, 2*m3);
}

/* -------------------------------------------------------------------------- *
 * Dense CG cache                                                             *
 *                                                                            *
 * Five-dimensional flat table, indexed by                                    *
 *   (two_j1, two_m1 + two_j1_max, two_j2, two_m2 + two_j2_max, two_J).       *
 * Storage is (two_j1_max + 1) × (2·two_j1_max + 1) × (two_j2_max + 1)        *
 *          × (2·two_j2_max + 1) × (two_j1_max + two_j2_max + 1) doubles;     *
 * invalid combos stored as 0. For two_j1_max = two_j2_max = 8 (int j ≤ 4)    *
 * this is ~3 MB; for 16/16 (int j ≤ 8) it is ~80 MB and callers should use   *
 * a smaller-scope table.                                                     *
 * -------------------------------------------------------------------------- */

struct cg_table {
    int    two_j1_max;
    int    two_j2_max;
    size_t stride_j1;
    size_t stride_m1;
    size_t stride_j2;
    size_t stride_m2;
    size_t stride_J;
    size_t data_size;
    double *data;
};

cg_table_t *irrep_cg_table_build_2j(int two_j1_max, int two_j2_max) {
    if (two_j1_max < 0 || two_j2_max < 0) return NULL;

    cg_table_t *t = calloc(1, sizeof(*t));
    if (!t) return NULL;
    t->two_j1_max = two_j1_max;
    t->two_j2_max = two_j2_max;

    size_t nm1 = (size_t)(2 * two_j1_max + 1);
    size_t nj2 = (size_t)(two_j2_max + 1);
    size_t nm2 = (size_t)(2 * two_j2_max + 1);
    size_t nJ  = (size_t)(two_j1_max + two_j2_max + 1);

    t->stride_J  = 1;
    t->stride_m2 = nJ;
    t->stride_j2 = nm2 * nJ;
    t->stride_m1 = nj2 * nm2 * nJ;
    t->stride_j1 = nm1 * nj2 * nm2 * nJ;

    size_t total = (size_t)(two_j1_max + 1) * t->stride_j1;
    if (total == 0) total = 1;
    t->data_size = total;
    t->data = calloc(total, sizeof(double));
    if (!t->data) { free(t); return NULL; }

    for (int two_j1 = 0; two_j1 <= two_j1_max; ++two_j1) {
        for (int two_m1 = -two_j1; two_m1 <= two_j1; ++two_m1) {
            for (int two_j2 = 0; two_j2 <= two_j2_max; ++two_j2) {
                for (int two_m2 = -two_j2; two_m2 <= two_j2; ++two_m2) {
                    int two_M = two_m1 + two_m2;
                    int two_J_min = iabs_(two_j1 - two_j2);
                    int two_J_top = two_j1 + two_j2;
                    for (int two_J = two_J_min; two_J <= two_J_top; two_J += 2) {
                        if (iabs_(two_M) > two_J) continue;
                        double v = irrep_cg_2j(two_j1, two_m1,
                                               two_j2, two_m2,
                                               two_J,  two_M);
                        size_t idx = (size_t)two_j1 * t->stride_j1
                                   + (size_t)(two_m1 + two_j1_max) * t->stride_m1
                                   + (size_t)two_j2 * t->stride_j2
                                   + (size_t)(two_m2 + two_j2_max) * t->stride_m2
                                   + (size_t)two_J * t->stride_J;
                        if (idx < t->data_size) t->data[idx] = v;
                    }
                }
            }
        }
    }
    return t;
}

cg_table_t *irrep_cg_table_build(int j1_max, int j2_max) {
    return irrep_cg_table_build_2j(2 * j1_max, 2 * j2_max);
}

void irrep_cg_table_free(cg_table_t *t) {
    if (!t) return;
    free(t->data);
    free(t);
}

double irrep_cg_lookup_2j(const cg_table_t *t,
                          int two_j1, int two_m1,
                          int two_j2, int two_m2,
                          int two_J,  int two_M) {
    if (!t || !t->data) return 0.0;
    if (two_j1 < 0 || two_j1 > t->two_j1_max) return 0.0;
    if (two_j2 < 0 || two_j2 > t->two_j2_max) return 0.0;
    if (two_m1 < -two_j1 || two_m1 > two_j1) return 0.0;
    if (two_m2 < -two_j2 || two_m2 > two_j2) return 0.0;
    if (two_m1 + two_m2 != two_M) return 0.0;
    if (two_J < 0 || two_J > t->two_j1_max + t->two_j2_max) return 0.0;
    if (two_J < iabs_(two_j1 - two_j2) || two_J > two_j1 + two_j2) return 0.0;
    if (iabs_(two_M) > two_J) return 0.0;

    size_t idx = (size_t)two_j1 * t->stride_j1
               + (size_t)(two_m1 + t->two_j1_max) * t->stride_m1
               + (size_t)two_j2 * t->stride_j2
               + (size_t)(two_m2 + t->two_j2_max) * t->stride_m2
               + (size_t)two_J * t->stride_J;
    if (idx >= t->data_size) return 0.0;
    return t->data[idx];
}

double irrep_cg_lookup(const cg_table_t *t,
                       int j1, int m1, int j2, int m2, int J, int M) {
    return irrep_cg_lookup_2j(t, 2*j1, 2*m1, 2*j2, 2*m2, 2*J, 2*M);
}

/* -------------------------------------------------------------------------- *
 * Dense block                                                                *
 *                                                                            *
 * Layout:                                                                    *
 *   out[((m1 + j1) · (2 j2 + 1) + (m2 + j2)) · (2 J + 1) + (M + J)]          *
 *     = ⟨j1 m1; j2 m2 | J M⟩.                                                *
 * -------------------------------------------------------------------------- */

void irrep_cg_block(int j1, int j2, int J, double *out) {
    int d1 = 2 * j1 + 1;
    int d2 = 2 * j2 + 1;
    int dJ = 2 * J  + 1;
    for (int i1 = 0; i1 < d1; ++i1) {
        int m1 = i1 - j1;
        for (int i2 = 0; i2 < d2; ++i2) {
            int m2 = i2 - j2;
            for (int iJ = 0; iJ < dJ; ++iJ) {
                int M = iJ - J;
                out[(i1 * d2 + i2) * dJ + iJ] = irrep_cg(j1, m1, j2, m2, J, M);
            }
        }
    }
}
