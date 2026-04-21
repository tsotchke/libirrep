/* SPDX-License-Identifier: MIT */
/* AVX2 (x86_64) kernels for batched radial / cutoff functions.
 *
 * Mirrors src/radial_neon.c but four elements per __m256d. Same
 * arithmetic order as the scalar reference; the AVX2 FMA intrinsics
 * (_mm256_fnmadd_pd, _mm256_fmadd_pd) align with the fmsub/fmadd
 * patterns that clang emits on the scalar side when FP_CONTRACT is
 * on, giving bit-exact output. FP_CONTRACT OFF is therefore NOT used
 * here — the goal is to match the scalar path, which IS contracted
 * on x86_64 optimized builds.
 *
 * The NEON equivalent keeps FP_CONTRACT on for the same reason.
 *
 * Compiled into the library only on x86_64 hosts with AVX2 + FMA
 * available. Dispatch in src/simd_runtime.c selects these kernels
 * when irrep_cpu_has_avx2() && irrep_cpu_has_fma().
 */

#if (defined(__x86_64__) || defined(_M_X64)) && defined(__AVX2__) && defined(__FMA__)

#include <immintrin.h>
#include <stddef.h>

#include "internal/dispatch.h"

/* NequIP polynomial cutoff:
 *   f_p(u) = 1 − (p+1)(p+2)/2 · u^p
 *              + p (p+2)      · u^{p+1}
 *              − p (p+1)/2    · u^{p+2},        u = r / r_cut.
 */
void irrep_cutoff_polynomial_batch_avx2(size_t N, const double *r, double r_cut, int p,
                                        double *out) {
    if (r_cut <= 0.0 || p < 1) {
        for (size_t i = 0; i < N; ++i)
            out[i] = 0.0;
        return;
    }
    const double  c1 = 0.5 * (double)(p + 1) * (double)(p + 2);
    const double  c2 = (double)p * (double)(p + 2);
    const double  c3 = 0.5 * (double)p * (double)(p + 1);
    const double  inv_rcut = 1.0 / r_cut;

    const __m256d vc1 = _mm256_set1_pd(c1);
    const __m256d vc2 = _mm256_set1_pd(c2);
    const __m256d vc3 = _mm256_set1_pd(c3);
    const __m256d vone = _mm256_set1_pd(1.0);
    const __m256d vzero = _mm256_setzero_pd();
    const __m256d vr_cut = _mm256_set1_pd(r_cut);
    const __m256d vinv_rcut = _mm256_set1_pd(inv_rcut);

    size_t        i = 0;
    for (; i + 4 <= N; i += 4) {
        __m256d vr = _mm256_loadu_pd(r + i);

        /* Validity mask: r ∈ [0, r_cut). */
        __m256d ge_zero = _mm256_cmp_pd(vr, vzero, _CMP_GE_OS);
        __m256d lt_cut = _mm256_cmp_pd(vr, vr_cut, _CMP_LT_OS);
        __m256d valid = _mm256_and_pd(ge_zero, lt_cut);

        __m256d vu = _mm256_mul_pd(vr, vinv_rcut);

        /* u^p by repeated multiplies. */
        __m256d vup = vu;
        for (int k = 1; k < p; ++k)
            vup = _mm256_mul_pd(vup, vu);
        __m256d vup1 = _mm256_mul_pd(vup, vu);
        __m256d vup2 = _mm256_mul_pd(vup1, vu);

        /* t = 1 − c1·u^p + c2·u^{p+1} − c3·u^{p+2} */
        __m256d t = _mm256_fnmadd_pd(vc1, vup, vone); /* 1 − c1·u^p */
        t = _mm256_fmadd_pd(vc2, vup1, t);            /* + c2·u^{p+1} */
        t = _mm256_fnmadd_pd(vc3, vup2, t);           /* − c3·u^{p+2} */

        t = _mm256_and_pd(t, valid);
        _mm256_storeu_pd(out + i, t);
    }
    for (; i < N; ++i) {
        double ri = r[i];
        if (ri < 0.0 || ri >= r_cut) {
            out[i] = 0.0;
            continue;
        }
        double u = ri * inv_rcut;
        double up = u;
        for (int k = 1; k < p; ++k)
            up *= u;
        double up1 = up * u;
        double up2 = up1 * u;
        out[i] = 1.0 - c1 * up + c2 * up1 - c3 * up2;
    }
}

/* f_p'(r) = −p(p+1)(p+2) / (2 r_cut) · u^{p−1} · (1 − u)². */
void irrep_cutoff_polynomial_d_batch_avx2(size_t N, const double *r, double r_cut, int p,
                                          double *out) {
    if (r_cut <= 0.0 || p < 1) {
        for (size_t i = 0; i < N; ++i)
            out[i] = 0.0;
        return;
    }
    const double  pre = -(double)p * (double)(p + 1) * (double)(p + 2) / (2.0 * r_cut);
    const double  inv_rcut = 1.0 / r_cut;

    const __m256d vpre = _mm256_set1_pd(pre);
    const __m256d vone = _mm256_set1_pd(1.0);
    const __m256d vzero = _mm256_setzero_pd();
    const __m256d vr_cut = _mm256_set1_pd(r_cut);
    const __m256d vinv_rcut = _mm256_set1_pd(inv_rcut);

    size_t        i = 0;
    for (; i + 4 <= N; i += 4) {
        __m256d vr = _mm256_loadu_pd(r + i);
        __m256d valid = _mm256_and_pd(_mm256_cmp_pd(vr, vzero, _CMP_GE_OS),
                                      _mm256_cmp_pd(vr, vr_cut, _CMP_LT_OS));

        __m256d vu = _mm256_mul_pd(vr, vinv_rcut);

        __m256d vupm1 = vone;
        for (int k = 0; k < p - 1; ++k)
            vupm1 = _mm256_mul_pd(vupm1, vu);

        __m256d vom = _mm256_sub_pd(vone, vu);
        /* pre·u^{p-1}·(1-u)² — left-to-right to match scalar associativity. */
        __m256d t = _mm256_mul_pd(vpre, vupm1);
        t = _mm256_mul_pd(t, vom);
        t = _mm256_mul_pd(t, vom);

        t = _mm256_and_pd(t, valid);
        _mm256_storeu_pd(out + i, t);
    }
    for (; i < N; ++i) {
        double ri = r[i];
        if (ri < 0.0 || ri >= r_cut) {
            out[i] = 0.0;
            continue;
        }
        double u = ri * inv_rcut;
        double upm1 = 1.0;
        for (int k = 0; k < p - 1; ++k)
            upm1 *= u;
        double om = 1.0 - u;
        out[i] = pre * upm1 * om * om;
    }
}

#else
typedef int irrep_radial_avx2_stub_t;
#endif
