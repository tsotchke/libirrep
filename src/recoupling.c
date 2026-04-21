/* SPDX-License-Identifier: MIT */
/* M8: Wigner 6j, 9j, Racah W — single-sum formulas in log-gamma form. */

#include <math.h>
#include <stdlib.h>

#include <irrep/recoupling.h>

#define UNUSED(x) ((void)(x))

static inline int iabs_(int x) {
    return x < 0 ? -x : x;
}

/* Triangle check: |j1 - j2| <= j3 <= j1 + j2 and j1+j2+j3 integer. */
static int triangle_2j_(int two_j1, int two_j2, int two_j3) {
    if (two_j1 < 0 || two_j2 < 0 || two_j3 < 0)
        return 0;
    if ((two_j1 + two_j2 + two_j3) & 1)
        return 0;
    if (two_j3 < iabs_(two_j1 - two_j2))
        return 0;
    if (two_j3 > two_j1 + two_j2)
        return 0;
    return 1;
}

/* Logarithm of Δ² = (j1+j2-j3)! (j1-j2+j3)! (-j1+j2+j3)! / (j1+j2+j3+1)! */
static double log_delta_sq_2j_(int two_j1, int two_j2, int two_j3) {
    int a = (two_j1 + two_j2 - two_j3) / 2;
    int b = (two_j1 - two_j2 + two_j3) / 2;
    int c = (-two_j1 + two_j2 + two_j3) / 2;
    int h = (two_j1 + two_j2 + two_j3) / 2;
    return lgamma(a + 1) + lgamma(b + 1) + lgamma(c + 1) - lgamma(h + 2);
}

/* Racah single-sum form:
 *   {j1 j2 j3; j4 j5 j6} = Δ(j1,j2,j3) Δ(j4,j5,j3) Δ(j4,j2,j6) Δ(j1,j5,j6)
 *                        · Σ_t (−1)^t (t+1)! / [(t−a)!(t−b)!(t−c)!(t−d)!(e−t)!(f−t)!(g−t)!]
 * with a = j1+j2+j3 etc., t ∈ [max(a,b,c,d), min(e,f,g)]. */
double irrep_wigner_6j_2j(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6) {
    if (!triangle_2j_(two_j1, two_j2, two_j3))
        return 0.0;
    if (!triangle_2j_(two_j4, two_j5, two_j3))
        return 0.0;
    if (!triangle_2j_(two_j4, two_j2, two_j6))
        return 0.0;
    if (!triangle_2j_(two_j1, two_j5, two_j6))
        return 0.0;

    int a = (two_j1 + two_j2 + two_j3) / 2;
    int b = (two_j4 + two_j5 + two_j3) / 2;
    int c = (two_j4 + two_j2 + two_j6) / 2;
    int d = (two_j1 + two_j5 + two_j6) / 2;
    int e = (two_j1 + two_j2 + two_j4 + two_j5) / 2;
    int f = (two_j2 + two_j3 + two_j5 + two_j6) / 2;
    int g = (two_j1 + two_j3 + two_j4 + two_j6) / 2;

    int t_min = a;
    if (b > t_min)
        t_min = b;
    if (c > t_min)
        t_min = c;
    if (d > t_min)
        t_min = d;
    int t_max = e;
    if (f < t_max)
        t_max = f;
    if (g < t_max)
        t_max = g;
    if (t_min > t_max)
        return 0.0;

    double log_deltas =
        0.5 * (log_delta_sq_2j_(two_j1, two_j2, two_j3) + log_delta_sq_2j_(two_j4, two_j5, two_j3) +
               log_delta_sq_2j_(two_j4, two_j2, two_j6) + log_delta_sq_2j_(two_j1, two_j5, two_j6));

    double sum = 0.0;
    for (int t = t_min; t <= t_max; ++t) {
        double log_term = lgamma(t + 2) - lgamma(t - a + 1) - lgamma(t - b + 1) -
                          lgamma(t - c + 1) - lgamma(t - d + 1) - lgamma(e - t + 1) -
                          lgamma(f - t + 1) - lgamma(g - t + 1);
        double term = exp(log_term);
        if (t & 1)
            term = -term;
        sum += term;
    }
    return exp(log_deltas) * sum;
}

double irrep_wigner_6j(int j1, int j2, int j3, int j4, int j5, int j6) {
    return irrep_wigner_6j_2j(2 * j1, 2 * j2, 2 * j3, 2 * j4, 2 * j5, 2 * j6);
}

/* 9j as single sum over k of three 6j products (Edmonds 6.4.3):
 *   {j1 j2 j3; j4 j5 j6; j7 j8 j9}
 *     = Σ_k (2k+1)(−1)^{2k} {j1 j4 j7; j8 j9 k} {j2 j5 j8; j4 k j6} {j3 j6 j9; k j1 j2}.
 */
double irrep_wigner_9j_2j(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6,
                          int two_j7, int two_j8, int two_j9) {
    int two_k_max = two_j1 + two_j9;
    if (two_j4 + two_j8 < two_k_max)
        two_k_max = two_j4 + two_j8;
    if (two_j2 + two_j6 < two_k_max)
        two_k_max = two_j2 + two_j6;

    int two_k_min = 0;
    int a = iabs_(two_j1 - two_j9);
    int b = iabs_(two_j4 - two_j8);
    int c = iabs_(two_j2 - two_j6);
    if (a > two_k_min)
        two_k_min = a;
    if (b > two_k_min)
        two_k_min = b;
    if (c > two_k_min)
        two_k_min = c;
    if (two_k_min > two_k_max)
        return 0.0;

    double sum = 0.0;
    for (int two_k = two_k_min; two_k <= two_k_max; ++two_k) {
        double s1 = irrep_wigner_6j_2j(two_j1, two_j4, two_j7, two_j8, two_j9, two_k);
        if (s1 == 0.0)
            continue;
        double s2 = irrep_wigner_6j_2j(two_j2, two_j5, two_j8, two_j4, two_k, two_j6);
        if (s2 == 0.0)
            continue;
        double s3 = irrep_wigner_6j_2j(two_j3, two_j6, two_j9, two_k, two_j1, two_j2);
        if (s3 == 0.0)
            continue;
        double phase = (two_k & 1) ? -1.0 : 1.0;
        sum += (double)(two_k + 1) * phase * s1 * s2 * s3;
    }
    return sum;
}

double irrep_wigner_9j(int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8, int j9) {
    return irrep_wigner_9j_2j(2 * j1, 2 * j2, 2 * j3, 2 * j4, 2 * j5, 2 * j6, 2 * j7, 2 * j8,
                              2 * j9);
}

/* Racah W(j1, j2, J, j3; j12, j23) = (−1)^{j1+j2+j3+J} {j1 j2 j12; j3 J j23}. */
double irrep_racah_w(int j1, int j2, int J, int j3, int j12, int j23) {
    double s6j = irrep_wigner_6j(j1, j2, j12, j3, J, j23);
    return ((j1 + j2 + j3 + J) & 1) ? -s6j : s6j;
}
