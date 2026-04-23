/* SPDX-License-Identifier: MIT */
/* AVX2 + FMA (x86_64) batch kernel for `irrep_wigner_d_matrix_batch`.
 *
 * Four β values per `__m256d`: lane 0 = β_b, lane 1 = β_{b+1}, lane 2 = β_{b+2},
 * lane 3 = β_{b+3}. Same vectorisation shape as the NEON kernel — the Jacobi
 * recurrence coefficients are scalar-broadcast across lanes; only `x = cos β`
 * differs. Tail handling: remainder of `n_betas % 4` falls to the scalar
 * `irrep_wigner_d_matrix` path.
 */

#if (defined(__x86_64__) || defined(_M_X64)) && defined(__AVX2__) && defined(__FMA__)

#include <immintrin.h>
#include <math.h>
#include <stddef.h>

#include <irrep/wigner_d.h>

#pragma STDC FP_CONTRACT OFF

extern double irrep_wigner_d_small(int j, int mp, int m, double beta);
extern void   irrep_wigner_d_matrix(int j, double *out, double beta);

static inline __m256d ipow_v4_(__m256d x, int n) {
    __m256d r = _mm256_set1_pd(1.0);
    __m256d base = x;
    while (n > 0) {
        if (n & 1)
            r = _mm256_mul_pd(r, base);
        n >>= 1;
        if (n)
            base = _mm256_mul_pd(base, base);
    }
    return r;
}

/* Jacobi P_n^{(α,β)} recurrence over four x-lanes. α, β, n shared. */
static __m256d jacobi_P_v4_(int n, int alpha, int beta, __m256d x_v) {
    if (n < 0)
        return _mm256_setzero_pd();
    if (n == 0)
        return _mm256_set1_pd(1.0);
    double  a = (double)alpha;
    double  b = (double)beta;
    __m256d p0 = _mm256_set1_pd(1.0);
    double  p1_const = 0.5 * (a - b);
    double  p1_slope = 0.5 * (a + b + 2.0);
    /* p1 = p1_slope * x + p1_const */
    __m256d p1 = _mm256_fmadd_pd(x_v, _mm256_set1_pd(p1_slope), _mm256_set1_pd(p1_const));
    if (n == 1)
        return p1;
    for (int k = 1; k < n; ++k) {
        double dk = (double)k;
        double two_k_ab = 2.0 * dk + a + b;
        double denom = 2.0 * (dk + 1.0) * (dk + a + b + 1.0) * two_k_ab;
        double slope = (two_k_ab + 1.0) * two_k_ab * (two_k_ab + 2.0);
        double offset = (two_k_ab + 1.0) * (a * a - b * b);
        double c2 = 2.0 * (dk + a) * (dk + b) * (two_k_ab + 2.0);
        /* c1_v = slope * x + offset */
        __m256d c1_v = _mm256_fmadd_pd(x_v, _mm256_set1_pd(slope), _mm256_set1_pd(offset));
        /* num = c1_v * p1 - c2 * p0 */
        __m256d num = _mm256_fmsub_pd(c1_v, p1, _mm256_mul_pd(p0, _mm256_set1_pd(c2)));
        __m256d p2 = _mm256_mul_pd(num, _mm256_set1_pd(1.0 / denom));
        p0 = p1;
        p1 = p2;
    }
    return p1;
}

static double canonicalise_int_(int *mp_io, int *m_io) {
    double sign = 1.0;
    int    mp = *mp_io;
    int    m = *m_io;
    if (m < mp) {
        int delta = mp - m;
        if (delta & 1)
            sign = -sign;
        int tmp = mp;
        mp = m;
        m = tmp;
    }
    if (m < -mp) {
        int tmp = -mp;
        mp = -m;
        m = tmp;
    }
    *mp_io = mp;
    *m_io = m;
    return sign;
}

static __m256d d_small_quad_(int j, int mp_in, int m_in, __m256d cb2_v, __m256d sb2_v,
                             __m256d x_v) {
    int mp = mp_in, m = m_in;
    if (m > j || m < -j || mp > j || mp < -j)
        return _mm256_setzero_pd();
    double sign = canonicalise_int_(&mp, &m);

    int mu = m - mp;
    int nu = m + mp;
    int n = j - m;

    int jpm = j + m, jmm = j - m, jpmp = j + mp, jmmp = j - mp;
    double lg_norm =
        0.5 * (lgamma((double)jpm + 1.0) + lgamma((double)jmm + 1.0) -
               lgamma((double)jpmp + 1.0) - lgamma((double)jmmp + 1.0));
    double norm = exp(lg_norm);

    __m256d cb2_pow = ipow_v4_(cb2_v, nu);
    __m256d sb2_pow = ipow_v4_(sb2_v, mu);
    __m256d Pn = jacobi_P_v4_(n, mu, nu, x_v);

    __m256d prod = _mm256_mul_pd(cb2_pow, sb2_pow);
    prod = _mm256_mul_pd(prod, Pn);
    prod = _mm256_mul_pd(prod, _mm256_set1_pd(sign * norm));
    return prod;
}

static void fill_matrix_quad_(int j, const double *betas, double **outs) {
    int d = 2 * j + 1;
    for (int k = 0; k < 4; ++k)
        for (int imp = 0; imp < d; ++imp)
            outs[k][imp * d + (d - 1 - imp)] = 0.0;

    __m256d beta_v = _mm256_loadu_pd(betas);
    __m256d half = _mm256_set1_pd(0.5);
    __m256d beta_half = _mm256_mul_pd(beta_v, half);
    double  bh[4];
    _mm256_storeu_pd(bh, beta_half);
    double cb2[4] = {cos(bh[0]), cos(bh[1]), cos(bh[2]), cos(bh[3])};
    double sb2[4] = {sin(bh[0]), sin(bh[1]), sin(bh[2]), sin(bh[3])};
    double x[4]   = {cos(betas[0]), cos(betas[1]), cos(betas[2]), cos(betas[3])};
    __m256d cb2_v = _mm256_loadu_pd(cb2);
    __m256d sb2_v = _mm256_loadu_pd(sb2);
    __m256d x_v = _mm256_loadu_pd(x);

    for (int imp = 0; imp < d; ++imp) {
        int mp = imp - j;
        int im_lo = (mp < 0) ? -mp + j : mp + j;
        for (int im = im_lo; im < d; ++im) {
            int     m = im - j;
            __m256d v = d_small_quad_(j, mp, m, cb2_v, sb2_v, x_v);
            double  lanes[4];
            _mm256_storeu_pd(lanes, v);
            for (int k = 0; k < 4; ++k) {
                double *o = outs[k];
                o[imp * d + im] = lanes[k];
                int imp_a = -m + j, im_a = -mp + j;
                o[imp_a * d + im_a] = lanes[k];
                if (mp != m) {
                    double s = ((m - mp) & 1) ? -lanes[k] : lanes[k];
                    o[im * d + imp] = s;
                    int imp_b = -mp + j, im_b = -m + j;
                    if (imp_b != im || im_b != imp)
                        o[imp_b * d + im_b] = s;
                }
            }
        }
    }
}

void irrep_wigner_d_matrix_batch_avx2(int j, size_t n_betas, const double *betas, double *out) {
    int    d = 2 * j + 1;
    size_t stride = (size_t)d * (size_t)d;
    size_t b = 0;
    for (; b + 4 <= n_betas; b += 4) {
        double *outs[4] = {out + (b + 0) * stride, out + (b + 1) * stride,
                           out + (b + 2) * stride, out + (b + 3) * stride};
        fill_matrix_quad_(j, betas + b, outs);
    }
    for (; b < n_betas; ++b)
        irrep_wigner_d_matrix(j, out + b * stride, betas[b]);
}

#else
/* Placeholder on non-AVX2 targets: keep the translation unit valid ISO C. */
typedef int irrep_wigner_d_avx2_disabled_t;
#endif /* x86_64 && AVX2 && FMA */
