/* SPDX-License-Identifier: MIT */
/* M5: Wigner-D and small-d matrices.
 *
 * Small-d via Sakurai (3.8.33) direct sum in log-gamma form — stable well
 * past j = 50 in double precision. Full D factors as
 *
 *     D^j_{m'm}(α, β, γ) = e^{−i m' α} · d^j_{m'm}(β) · e^{−i m γ}.
 *
 * Analytic ∂d/∂β from term-wise differentiation of the direct formula.
 * Multiset block-diagonal builder for applying D(R) to a feature vector.
 */

#include <complex.h>
#include <math.h>
#include <stddef.h>

#include <irrep/types.h>
#include <irrep/wigner_d.h>

#define UNUSED(x) ((void)(x))
#define IRREP_MAX_BLOCK_DIM (2 * IRREP_L_MAX + 1)

/* -------------------------------------------------------------------------- *
 * small-d: Sakurai (3.8.33) direct sum, evaluated in log-gamma form          *
 * -------------------------------------------------------------------------- */

double irrep_wigner_d_small_2j(int two_j, int two_mp, int two_m, double beta) {
    if (two_j < 0) return 0.0;
    if (two_m  < -two_j || two_m  > two_j ) return 0.0;
    if (two_mp < -two_j || two_mp > two_j ) return 0.0;
    if ((two_j + two_m ) & 1) return 0.0;
    if ((two_j + two_mp) & 1) return 0.0;

    int p   = (two_j + two_m ) / 2;  /* j + m       */
    int q   = (two_j - two_m ) / 2;  /* j − m       */
    int r   = (two_j + two_mp) / 2;  /* j + m'      */
    int s   = (two_j - two_mp) / 2;  /* j − m'      */
    int dmm = (two_m - two_mp) / 2;  /* m − m'      */
    int dmp = (two_m - two_mp) / 2;  /* m − m' (alias; keeps naming local) */
    UNUSED(dmp);

    int k_min = 0;
    if (dmm > k_min) k_min = dmm;
    int k_max = p;
    if (s < k_max) k_max = s;
    if (k_min > k_max) return 0.0;

    double cb2 = cos(0.5 * beta);
    double sb2 = sin(0.5 * beta);

    double half_log_fact = 0.5 * (lgamma(p + 1) + lgamma(q + 1)
                                + lgamma(r + 1) + lgamma(s + 1));

    double sum = 0.0;
    for (int k = k_min; k <= k_max; ++k) {
        int p1 = two_j - 2 * k + (two_m  - two_mp) / 2;   /* 2j − 2k + m − m' */
        int p2 =         2 * k + (two_mp - two_m ) / 2;   /* 2k − m + m'      */

        double log_denom = lgamma(p - k + 1) + lgamma(k + 1)
                         + lgamma(s - k + 1) + lgamma(k - dmm + 1);
        double log_A = half_log_fact - log_denom;
        double sign  = ((k - dmm) & 1) ? -1.0 : 1.0;

        double term = sign * exp(log_A) * pow(cb2, (double)p1) * pow(sb2, (double)p2);
        sum += term;
    }
    return sum;
}

double irrep_wigner_d_small(int j, int mp, int m, double beta) {
    return irrep_wigner_d_small_2j(2 * j, 2 * mp, 2 * m, beta);
}

/* -------------------------------------------------------------------------- *
 * Full D = e^{−i m' α} · d^j_{m'm}(β) · e^{−i m γ}                           *
 * -------------------------------------------------------------------------- */

double _Complex irrep_wigner_D_2j(int two_j, int two_mp, int two_m,
                                  double alpha, double beta, double gamma) {
    double d = irrep_wigner_d_small_2j(two_j, two_mp, two_m, beta);
    if (d == 0.0) return 0.0;
    double phase = -(0.5 * two_mp * alpha + 0.5 * two_m * gamma);
    return d * (cos(phase) + I * sin(phase));
}

double _Complex irrep_wigner_D(int j, int mp, int m,
                               double alpha, double beta, double gamma) {
    return irrep_wigner_D_2j(2 * j, 2 * mp, 2 * m, alpha, beta, gamma);
}

void irrep_wigner_D_matrix(int j, double _Complex *out,
                           double alpha, double beta, double gamma) {
    if (j < 0) return;
    int d = 2 * j + 1;
    for (int imp = 0; imp < d; ++imp) {
        int mp = imp - j;
        for (int im = 0; im < d; ++im) {
            int m = im - j;
            out[imp * d + im] = irrep_wigner_D(j, mp, m, alpha, beta, gamma);
        }
    }
}

void irrep_wigner_d_matrix(int j, double *out, double beta) {
    if (j < 0) return;
    int d = 2 * j + 1;
    for (int imp = 0; imp < d; ++imp) {
        int mp = imp - j;
        for (int im = 0; im < d; ++im) {
            int m = im - j;
            out[imp * d + im] = irrep_wigner_d_small(j, mp, m, beta);
        }
    }
}

/* -------------------------------------------------------------------------- *
 * Multiset block-diagonal D                                                  *
 *                                                                            *
 * For an irrep_multiset_t m = ⊕ mult_i (l_i, parity_i), writes a             *
 * total_dim × total_dim complex matrix where each (l_i)-block appears        *
 * mult_i times along the diagonal. Parity plays no role under SO(3) —        *
 * rotations preserve parity — so both parities use the same D^{l}(R).        *
 * -------------------------------------------------------------------------- */

void irrep_wigner_D_multiset(const irrep_multiset_t *m, double _Complex *out,
                             double alpha, double beta, double gamma) {
    if (!m) return;
    int total = m->total_dim;
    for (int i = 0; i < total * total; ++i) out[i] = 0.0;

    double _Complex block[IRREP_MAX_BLOCK_DIM * IRREP_MAX_BLOCK_DIM];
    int offset = 0;

    for (int t = 0; t < m->num_terms; ++t) {
        int l        = m->labels[t].l;
        int mult     = m->multiplicities[t];
        int bdim     = 2 * l + 1;
        irrep_wigner_D_matrix(l, block, alpha, beta, gamma);

        for (int c = 0; c < mult; ++c) {
            for (int i = 0; i < bdim; ++i) {
                for (int j = 0; j < bdim; ++j) {
                    out[(offset + i) * total + (offset + j)] = block[i * bdim + j];
                }
            }
            offset += bdim;
        }
    }
}

/* -------------------------------------------------------------------------- *
 * Analytic ∂d^j_{m'm}/∂β                                                      *
 *                                                                            *
 * Differentiate each term A · cos(β/2)^{p1} · sin(β/2)^{p2} with respect to   *
 * β. The chain rule gives                                                    *
 *                                                                            *
 *   −(p1/2) · sin(β/2) · cos(β/2)^{p1−1} · sin(β/2)^{p2}                     *
 *   + (p2/2) · cos(β/2) · cos(β/2)^{p1}   · sin(β/2)^{p2−1}                  *
 *                                                                            *
 * Guard pow(x, −1) at the poles: when p1 = 0 or p2 = 0 the corresponding     *
 * contribution is multiplied by zero anyway, so skip the whole term.         *
 * -------------------------------------------------------------------------- */

double irrep_wigner_d_small_dbeta(int j, int mp, int m, double beta) {
    int two_j = 2 * j, two_mp = 2 * mp, two_m = 2 * m;

    if (two_j < 0) return 0.0;
    if (two_m  < -two_j || two_m  > two_j ) return 0.0;
    if (two_mp < -two_j || two_mp > two_j ) return 0.0;
    if ((two_j + two_m ) & 1) return 0.0;
    if ((two_j + two_mp) & 1) return 0.0;

    int pv  = (two_j + two_m ) / 2;
    int qv  = (two_j - two_m ) / 2;
    int rv  = (two_j + two_mp) / 2;
    int sv  = (two_j - two_mp) / 2;
    int dmm = (two_m - two_mp) / 2;

    int k_min = 0;
    if (dmm > k_min) k_min = dmm;
    int k_max = pv;
    if (sv < k_max) k_max = sv;
    if (k_min > k_max) return 0.0;

    double cb2 = cos(0.5 * beta);
    double sb2 = sin(0.5 * beta);

    double half_log_fact = 0.5 * (lgamma(pv + 1) + lgamma(qv + 1)
                                + lgamma(rv + 1) + lgamma(sv + 1));

    double sum = 0.0;
    for (int k = k_min; k <= k_max; ++k) {
        int p1 = two_j - 2 * k + (two_m  - two_mp) / 2;
        int p2 =         2 * k + (two_mp - two_m ) / 2;

        double log_denom = lgamma(pv - k + 1) + lgamma(k + 1)
                         + lgamma(sv - k + 1) + lgamma(k - dmm + 1);
        double log_A = half_log_fact - log_denom;
        double A     = exp(log_A);
        double sign  = ((k - dmm) & 1) ? -1.0 : 1.0;

        double contrib = 0.0;
        if (p1 > 0) {
            contrib -= 0.5 * (double)p1 * sb2 * pow(cb2, (double)(p1 - 1))
                                              * pow(sb2, (double)p2);
        }
        if (p2 > 0) {
            contrib += 0.5 * (double)p2 * cb2 * pow(cb2, (double)p1)
                                              * pow(sb2, (double)(p2 - 1));
        }
        sum += sign * A * contrib;
    }
    return sum;
}
