/* SPDX-License-Identifier: MIT */
/* M4: Clebsch-Gordan coefficients and Wigner 3j symbols.
 *
 * Algorithm: Racah single-sum formula in log-gamma form. Stable well past
 * j = 50 in double precision; `long double` is not needed. Doubled-integer
 * arguments (`_2j` variants) are primary; integer wrappers multiply by 2.
 *
 * See Racah, Phys. Rev. 62, 438 (1942); Varshalovich §8; Edmonds §3.
 */

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include <irrep/clebsch_gordan.h>

#define UNUSED(x) ((void)(x))

static inline int iabs_(int x) { return x < 0 ? -x : x; }

/* -------------------------------------------------------------------------- *
 * Clebsch-Gordan — Racah formula                                             *
 * -------------------------------------------------------------------------- */

double irrep_cg_2j(int two_j1, int two_m1,
                   int two_j2, int two_m2,
                   int two_J,  int two_M) {
    /* ---- selection rules ---- */
    if (two_j1 < 0 || two_j2 < 0 || two_J < 0) return 0.0;
    if (two_m1 + two_m2 != two_M) return 0.0;
    if (two_m1 < -two_j1 || two_m1 > two_j1) return 0.0;
    if (two_m2 < -two_j2 || two_m2 > two_j2) return 0.0;
    if (two_M  < -two_J  || two_M  > two_J ) return 0.0;
    if (two_J < iabs_(two_j1 - two_j2) || two_J > two_j1 + two_j2) return 0.0;
    /* parity: m and j must share parity; triangle sum must be integer */
    if ((two_j1 + two_m1) & 1) return 0.0;
    if ((two_j2 + two_m2) & 1) return 0.0;
    if ((two_J  + two_M ) & 1) return 0.0;
    if ((two_j1 + two_j2 + two_J) & 1) return 0.0;

    /* All of the following are guaranteed non-negative integers. */
    int a = (two_j1 + two_j2 - two_J) / 2;          /* j1 + j2 − J        */
    int b = (two_j1 - two_m1) / 2;                  /* j1 − m1            */
    int c = (two_j2 + two_m2) / 2;                  /* j2 + m2            */
    int d = (two_J  - two_j2 + two_m1) / 2;         /* J − j2 + m1        */
    int e = (two_J  - two_j1 - two_m2) / 2;         /* J − j1 − m2        */
    int f = (two_j1 - two_j2 + two_J) / 2;          /* j1 − j2 + J        */
    int g = (-two_j1 + two_j2 + two_J) / 2;         /* −j1 + j2 + J       */
    int h = (two_j1 + two_j2 + two_J) / 2;          /* j1 + j2 + J        */
    int j1pm1 = (two_j1 + two_m1) / 2;              /* j1 + m1            */
    int j2mm2 = (two_j2 - two_m2) / 2;              /* j2 − m2            */
    int JpM   = (two_J  + two_M ) / 2;              /* J  + M             */
    int JmM   = (two_J  - two_M ) / 2;              /* J  − M             */

    /* ln Δ² = ln(a!) + ln(f!) + ln(g!) − ln((h + 1)!) */
    double log_delta_sq = lgamma(a + 1) + lgamma(f + 1)
                        + lgamma(g + 1) - lgamma(h + 2);
    /* ln[(j1+m1)!(j1−m1)!(j2+m2)!(j2−m2)!(J+M)!(J−M)!] */
    double log_pre = lgamma(j1pm1 + 1) + lgamma(b + 1)
                   + lgamma(c + 1)     + lgamma(j2mm2 + 1)
                   + lgamma(JpM + 1)   + lgamma(JmM + 1);

    /* k-range: every factorial argument in the Racah sum must be non-negative. */
    int k_min = 0;
    if (-d > k_min) k_min = -d;
    if (-e > k_min) k_min = -e;
    int k_max = a;
    if (b < k_max) k_max = b;
    if (c < k_max) k_max = c;
    if (k_min > k_max) return 0.0;

    /* Σ_k (−1)^k / [k! (a−k)! (b−k)! (c−k)! (d+k)! (e+k)!] */
    double sum = 0.0;
    for (int k = k_min; k <= k_max; ++k) {
        double log_term = -(lgamma(k + 1)
                          + lgamma(a - k + 1)
                          + lgamma(b - k + 1)
                          + lgamma(c - k + 1)
                          + lgamma(d + k + 1)
                          + lgamma(e + k + 1));
        double term = exp(log_term);
        if (k & 1) term = -term;
        sum += term;
    }

    double log_overall = 0.5 * log((double)two_J + 1.0)
                       + 0.5 * log_delta_sq
                       + 0.5 * log_pre;
    return exp(log_overall) * sum;
}

double irrep_cg(int j1, int m1, int j2, int m2, int J, int M) {
    return irrep_cg_2j(2*j1, 2*m1, 2*j2, 2*m2, 2*J, 2*M);
}

/* -------------------------------------------------------------------------- *
 * Wigner 3j from CG                                                          *
 *                                                                            *
 *   (j1 j2 j3; m1 m2 m3)                                                     *
 *     = (−1)^{j1 − j2 − m3} / √(2 j3 + 1) · ⟨j1 m1; j2 m2 | j3, −m3⟩.        *
 * -------------------------------------------------------------------------- */

double irrep_wigner_3j_2j(int two_j1, int two_m1,
                          int two_j2, int two_m2,
                          int two_j3, int two_m3) {
    if (two_m1 + two_m2 + two_m3 != 0) return 0.0;
    double cg = irrep_cg_2j(two_j1, two_m1, two_j2, two_m2, two_j3, -two_m3);
    if (cg == 0.0) return 0.0;
    int phase_half = (two_j1 - two_j2 - two_m3) / 2;
    double phase   = (phase_half & 1) ? -1.0 : 1.0;
    return phase * cg / sqrt((double)two_j3 + 1.0);
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
