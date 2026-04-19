/* SPDX-License-Identifier: MIT */
/* NEON (aarch64) kernels for batched radial/cutoff functions.
 *
 * Bit-exact against the scalar reference on representative inputs: the
 * implementation uses the same scalar fused multiply-add tree, expressed
 * over float64x2_t lanes. No transcendentals — the polynomial cutoffs
 * reduce to repeated multiplies.
 *
 * Compiled into the library only on aarch64 hosts. The dispatch table in
 * src/simd_runtime.c selects these kernels when irrep_cpu_has_neon() is true.
 */

#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)

#include <arm_neon.h>
#include <stddef.h>

#include "internal/dispatch.h"

/* NequIP polynomial cutoff:
 *   f_p(u) = 1 − (p+1)(p+2)/2 · u^p
 *              + p (p+2)      · u^{p+1}
 *              − p (p+1)/2    · u^{p+2},        u = r / r_cut,
 *
 * supported on u ∈ [0, 1); zero for u ≥ 1. */
void irrep_cutoff_polynomial_batch_neon(size_t N, const double *r,
                                         double r_cut, int p, double *out) {
    if (r_cut <= 0.0 || p < 1) {
        for (size_t i = 0; i < N; ++i) out[i] = 0.0;
        return;
    }
    const double c1 = 0.5 * (double)(p + 1) * (double)(p + 2);
    const double c2 =       (double)p       * (double)(p + 2);
    const double c3 = 0.5 * (double)p       * (double)(p + 1);
    const double inv_rcut = 1.0 / r_cut;

    const float64x2_t vc1       = vdupq_n_f64(c1);
    const float64x2_t vc2       = vdupq_n_f64(c2);
    const float64x2_t vc3       = vdupq_n_f64(c3);
    const float64x2_t vone      = vdupq_n_f64(1.0);
    const float64x2_t vzero     = vdupq_n_f64(0.0);
    const float64x2_t vinv_rcut = vdupq_n_f64(inv_rcut);

    size_t i = 0;
    for (; i + 2 <= N; i += 2) {
        float64x2_t vr = vld1q_f64(r + i);

        /* mask: lanes where r ∈ [0, r_cut); elsewhere the result is 0. */
        uint64x2_t valid = vandq_u64(
            vcgeq_f64(vr, vzero),
            vcltq_f64(vr, vdupq_n_f64(r_cut)));

        float64x2_t vu = vmulq_f64(vr, vinv_rcut);

        /* u^p by repeated squaring / multiply (p ≥ 1 runtime value). */
        float64x2_t vup = vu;
        for (int k = 1; k < p; ++k) vup = vmulq_f64(vup, vu);
        float64x2_t vup1 = vmulq_f64(vup, vu);
        float64x2_t vup2 = vmulq_f64(vup1, vu);

        /* Match the scalar accumulation order and fused-mul-add exactly —
         * clang contracts the scalar `1 - c1*a + c2*b - c3*c` to fmsub/fmadd
         * on Apple Silicon; we use the corresponding NEON FMA intrinsics so
         * the lanes are bit-identical to the scalar reference. */
        float64x2_t t = vfmsq_f64(vone, vc1, vup);          /* 1 - c1·u^p */
        t = vfmaq_f64(t, vc2, vup1);                        /* +  c2·u^{p+1} */
        t = vfmsq_f64(t, vc3, vup2);                        /* -  c3·u^{p+2} */

        /* Zero out invalid lanes. */
        t = vreinterpretq_f64_u64(
                vandq_u64(vreinterpretq_u64_f64(t), valid));

        vst1q_f64(out + i, t);
    }
    /* Tail: one element. */
    for (; i < N; ++i) {
        double ri = r[i];
        if (ri < 0.0 || ri >= r_cut) { out[i] = 0.0; continue; }
        double u  = ri * inv_rcut;
        double up = u;
        for (int k = 1; k < p; ++k) up *= u;
        double up1 = up * u;
        double up2 = up1 * u;
        out[i] = 1.0 - c1 * up + c2 * up1 - c3 * up2;
    }
}

/* Derivative:
 *   f_p'(r) = −p(p+1)(p+2) / (2 r_cut) · u^{p−1} · (1 − u)². */
void irrep_cutoff_polynomial_d_batch_neon(size_t N, const double *r,
                                           double r_cut, int p, double *out) {
    if (r_cut <= 0.0 || p < 1) {
        for (size_t i = 0; i < N; ++i) out[i] = 0.0;
        return;
    }
    const double pre =
        -(double)p * (double)(p + 1) * (double)(p + 2) / (2.0 * r_cut);
    const double inv_rcut = 1.0 / r_cut;

    const float64x2_t vpre      = vdupq_n_f64(pre);
    const float64x2_t vone      = vdupq_n_f64(1.0);
    const float64x2_t vzero     = vdupq_n_f64(0.0);
    const float64x2_t vinv_rcut = vdupq_n_f64(inv_rcut);

    size_t i = 0;
    for (; i + 2 <= N; i += 2) {
        float64x2_t vr = vld1q_f64(r + i);

        uint64x2_t valid = vandq_u64(
            vcgeq_f64(vr, vzero),
            vcltq_f64(vr, vdupq_n_f64(r_cut)));

        float64x2_t vu = vmulq_f64(vr, vinv_rcut);

        /* u^{p-1}: for p=1 this is 1. */
        float64x2_t vupm1 = vone;
        for (int k = 0; k < p - 1; ++k) vupm1 = vmulq_f64(vupm1, vu);

        float64x2_t vom = vsubq_f64(vone, vu);
        /* Match scalar left-to-right associativity: pre*upm1*om*om */
        float64x2_t t = vmulq_f64(vpre, vupm1);
        t = vmulq_f64(t, vom);
        t = vmulq_f64(t, vom);

        t = vreinterpretq_f64_u64(
                vandq_u64(vreinterpretq_u64_f64(t), valid));

        vst1q_f64(out + i, t);
    }
    for (; i < N; ++i) {
        double ri = r[i];
        if (ri < 0.0 || ri >= r_cut) { out[i] = 0.0; continue; }
        double u = ri * inv_rcut;
        double upm1 = 1.0;
        for (int k = 0; k < p - 1; ++k) upm1 *= u;
        double om = 1.0 - u;
        out[i] = pre * upm1 * om * om;
    }
}

#else
/* Keep the translation unit non-empty on non-aarch64 hosts to avoid
 * "file contained no source" warnings. */
typedef int irrep_radial_neon_stub_t;
#endif
