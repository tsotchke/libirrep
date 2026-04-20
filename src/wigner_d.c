/* SPDX-License-Identifier: MIT */
/* M5: Wigner-D and small-d matrices.
 *
 * Small-d via Jacobi-polynomial form (Edmonds 4.1.23 / NIST DLMF §34.2):
 *
 *     d^j_{m'm}(β) = ξ · sqrt( (j+m)!(j-m)! / (j+m')!(j-m')! )
 *                    · (cos β/2)^{m+m'} · (sin β/2)^{m−m'}
 *                    · P_{j−m}^{(m−m', m+m')}(cos β)
 *
 * valid in the canonical region m ≥ |m'|; other regions reduced by the
 * Varshalovich §4.4.1 symmetries (swap m ↔ m' with a (−1)^{m−m'} phase,
 * and/or (m, m') → (−m', −m)). The Jacobi polynomial itself is evaluated
 * by the NIST DLMF §18.9.1 three-term recurrence in n — stable forward for
 * non-negative integer (α, β) at x ∈ [−1, 1]. Measured unitarity at
 * (α,β,γ) = (0.3, 0.9, 1.5) is better than 1e-12 up to j ≈ 80 and bounded
 * by the lgamma precision limit (j ≈ 170) past that.
 *
 * The earlier Sakurai (3.8.33) direct sum was replaced here: it loses
 * precision to catastrophic cancellation past j ≈ 20 (j = 50 unitarity
 * err ≈ 2e-3, j = 60 unitarity err O(10), j = 80 divergent).
 *
 * Full D factors as
 *     D^j_{m'm}(α, β, γ) = e^{−i m' α} · d^j_{m'm}(β) · e^{−i m γ}.
 *
 * Analytic ∂d/∂β from product rule + Jacobi derivative identity
 *     (d/dx) P_n^{(α,β)}(x) = (n + α + β + 1)/2 · P_{n−1}^{(α+1, β+1)}(x).
 * Multiset block-diagonal builder for applying D(R) to a feature vector.
 */

#include <complex.h>
#include <math.h>
#include <stddef.h>

#include <irrep/types.h>
#include <irrep/wigner_d.h>

#define UNUSED(x) ((void)(x))
#define IRREP_MAX_BLOCK_DIM (2 * IRREP_L_MAX + 1)

/* x^n for non-negative integer n via repeated squaring — faster than
 * libm pow() on runtime integer exponents and yields identical results
 * across platforms (libm pow has sub-ULP variation across implementations). */
static inline double ipow_(double x, int n) {
    double r = 1.0;
    double base = x;
    while (n > 0) {
        if (n & 1) r *= base;
        n >>= 1;
        if (n) base *= base;
    }
    return r;
}

/* -------------------------------------------------------------------------- *
 * Jacobi polynomial P_n^{(α,β)}(x) via DLMF §18.9.1 forward recurrence.      *
 * Stable for (α, β) non-negative and x ∈ [−1, 1]; no cancellation across n.  *
 * -------------------------------------------------------------------------- */

static double jacobi_P_(int n, int alpha, int beta, double x) {
    if (n < 0)  return 0.0;
    if (n == 0) return 1.0;
    double a  = (double)alpha;
    double b  = (double)beta;
    double p0 = 1.0;
    double p1 = 0.5 * ((a - b) + (a + b + 2.0) * x);
    if (n == 1) return p1;
    for (int k = 1; k < n; ++k) {
        double dk = (double)k;
        double two_k_ab = 2.0 * dk + a + b;
        double denom = 2.0 * (dk + 1.0) * (dk + a + b + 1.0) * two_k_ab;
        double c1    = (two_k_ab + 1.0) * (two_k_ab * (two_k_ab + 2.0) * x + (a * a - b * b));
        double c2    = 2.0 * (dk + a) * (dk + b) * (two_k_ab + 2.0);
        double p2    = (c1 * p1 - c2 * p0) / denom;
        p0 = p1;
        p1 = p2;
    }
    return p1;
}

/* (d/dx) P_n^{(α,β)}(x) = (n + α + β + 1)/2 · P_{n-1}^{(α+1, β+1)}(x) */
static double jacobi_dP_dx_(int n, int alpha, int beta, double x) {
    if (n <= 0) return 0.0;
    double factor = 0.5 * (double)(n + alpha + beta + 1);
    return factor * jacobi_P_(n - 1, alpha + 1, beta + 1, x);
}

/* -------------------------------------------------------------------------- *
 * small-d: Jacobi-polynomial form with symmetry reduction to m ≥ |m'|.       *
 * -------------------------------------------------------------------------- */

/* Reduce (two_j, two_mp, two_m) to the canonical region two_m ≥ |two_mp| via
 * the Varshalovich §4.4.1 symmetries. Returns the accumulated sign. */
static double canonicalise_(int two_j, int *two_mp, int *two_m) {
    double sign = 1.0;
    int    mp   = *two_mp;
    int    m    = *two_m;
    if (m < mp) {
        /* d^j_{m',m} = (−1)^{m'−m} · d^j_{m,m'}  (parity guarantees 2m−2m' is even) */
        int delta_half = (mp - m) / 2;
        if (delta_half & 1) sign = -sign;
        int tmp = mp; mp = m; m = tmp;
    }
    if (m < -mp) {
        /* d^j_{m',m} = d^j_{-m,-m'}  (no sign change) */
        int tmp = -mp; mp = -m; m = tmp;
    }
    (void)two_j;
    *two_mp = mp;
    *two_m  = m;
    return sign;
}

double irrep_wigner_d_small_2j(int two_j, int two_mp, int two_m, double beta) {
    if (two_j < 0) return 0.0;
    if (two_m  < -two_j || two_m  > two_j ) return 0.0;
    if (two_mp < -two_j || two_mp > two_j ) return 0.0;
    if ((two_j + two_m ) & 1) return 0.0;
    if ((two_j + two_mp) & 1) return 0.0;

    double sign = canonicalise_(two_j, &two_mp, &two_m);
    /* Canonical region: two_m ≥ |two_mp|. Integer labels: */
    int mu_int = (two_m - two_mp) / 2;   /* m − m' ≥ 0 */
    int nu_int = (two_m + two_mp) / 2;   /* m + m' ≥ 0 */
    int n      = (two_j - two_m ) / 2;   /* j − m     ≥ 0 */

    int jpm  = (two_j + two_m ) / 2;
    int jmm  = (two_j - two_m ) / 2;
    int jpmp = (two_j + two_mp) / 2;
    int jmmp = (two_j - two_mp) / 2;

    double lg_norm = 0.5 * (lgamma(jpm + 1.0) + lgamma(jmm + 1.0)
                          - lgamma(jpmp + 1.0) - lgamma(jmmp + 1.0));
    double norm    = exp(lg_norm);

    double cb2 = cos(0.5 * beta);
    double sb2 = sin(0.5 * beta);
    double x   = cos(beta);

    double cb2_pow = ipow_(cb2, nu_int);
    double sb2_pow = ipow_(sb2, mu_int);
    double Pn      = jacobi_P_(n, mu_int, nu_int, x);

    return sign * norm * cb2_pow * sb2_pow * Pn;
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
    /* Exploit two Varshalovich §4.4.1 symmetries to cut the work ≈4× on the
     * full d×d matrix:
     *   d^j_{m', m}(β) = d^j_{−m, −m'}(β)                (antipodal swap)
     *   d^j_{m', m}(β) = (−1)^{m'−m} · d^j_{m, m'}(β)    (transpose with sign)
     * Compute only cells where m ≥ |m'| AND (m', m) lex ≤ (−m, −m'); fill
     * the three symmetric images from each. */
    for (int imp = 0; imp < d; ++imp) out[imp * d + (d - 1 - imp)] = 0.0;  /* zero scratch on anti-diagonal — refilled below */

    for (int imp = 0; imp < d; ++imp) {
        int mp = imp - j;
        int im_lo = (mp < 0) ? -mp + j : mp + j;   /* smallest im with m ≥ |m'| */
        for (int im = im_lo; im < d; ++im) {
            int m = im - j;
            double v = irrep_wigner_d_small(j, mp, m, beta);
            out[imp * d + im] = v;                                   /* (m', m) */
            /* Image 1: (−m, −m') → same value */
            int imp_a = -m  + j;
            int im_a  = -mp + j;
            out[imp_a * d + im_a] = v;
            /* Image 2: (m, m') → (-1)^{m-m'} · v  (skip if equal to (m', m)) */
            if (mp != m) {
                double sign = ((m - mp) & 1) ? -v : v;
                out[im * d + imp] = sign;
                /* Image 3: (−m', −m) → (-1)^{m-m'} · v  (mirror of image 2) */
                int imp_b = -mp + j;
                int im_b  = -m  + j;
                if (imp_b != im || im_b != imp) {
                    out[imp_b * d + im_b] = sign;
                }
            }
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

    double sign = canonicalise_(two_j, &two_mp, &two_m);
    int mu_int = (two_m - two_mp) / 2;
    int nu_int = (two_m + two_mp) / 2;
    int n      = (two_j - two_m ) / 2;

    int jpm  = (two_j + two_m ) / 2;
    int jmm  = (two_j - two_m ) / 2;
    int jpmp = (two_j + two_mp) / 2;
    int jmmp = (two_j - two_mp) / 2;

    double lg_norm = 0.5 * (lgamma(jpm + 1.0) + lgamma(jmm + 1.0)
                          - lgamma(jpmp + 1.0) - lgamma(jmmp + 1.0));
    double norm    = exp(lg_norm);

    double cb2 = cos(0.5 * beta);
    double sb2 = sin(0.5 * beta);
    double x   = cos(beta);
    double sbeta = sin(beta);

    double Pn   = jacobi_P_   (n, mu_int, nu_int, x);
    double dPdx = jacobi_dP_dx_(n, mu_int, nu_int, x);

    /* f(β) = sign · norm · cb2^ν · sb2^μ · P_n(cos β)
     * f'(β) = sign · norm · [ d/dβ(cb2^ν · sb2^μ) · P + cb2^ν · sb2^μ · (−sin β) · dP/dx ]
     * d/dβ(cb2^ν · sb2^μ) = (μ/2) cb2^{ν+1} sb2^{μ−1} − (ν/2) cb2^{ν−1} sb2^{μ+1}
     * (Using d(cos β/2)/dβ = −sb2/2 and d(sin β/2)/dβ = cb2/2.)
     */
    double env_d = 0.0;
    if (mu_int > 0) {
        env_d += 0.5 * (double)mu_int * ipow_(cb2, nu_int + 1)
                                      * ipow_(sb2, mu_int - 1);
    }
    if (nu_int > 0) {
        env_d -= 0.5 * (double)nu_int * ipow_(cb2, nu_int - 1)
                                      * ipow_(sb2, mu_int + 1);
    }
    double env   = ipow_(cb2, nu_int) * ipow_(sb2, mu_int);
    double total = env_d * Pn + env * (-sbeta) * dPdx;

    return sign * norm * total;
}
