/* SPDX-License-Identifier: MIT */
/* NEON (aarch64) batch-flat kernel for `irrep_tp_apply_weighted_batch_flat`.
 *
 * Dim-first layout: for a buffer `a_in` of shape `[a_dim][batch]`, the
 * element at irrep-slot `i` of sample `bi` lives at `a_in[i * batch + bi]`.
 * Consecutive batch lanes are contiguous, so two lanes per `float64x2_t`
 * give strict contiguous vector loads/stores without gather/scatter.
 *
 * Inner loop, per CG sparse entry (ima, imb, imc, cg), per u-copy, per
 * path:
 *
 *     c_row[bi..bi+1] += (w · cg) · a_row[bi..bi+1] · b_row[bi..bi+1]
 *
 * The (w · cg) prefactor is scalar-broadcast; a_row / b_row / c_row are
 * contiguous over bi. Tail (`batch % 2`) falls to scalar.
 */

#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)

#include <arm_neon.h>
#include <stddef.h>
#include <string.h>

#include <irrep/tensor_product.h>

#pragma STDC FP_CONTRACT OFF

/* The internal tp_path / tp_nz_entry / tp_descriptor types are declared in
 * src/tensor_product.c. Mirror the minimum shape needed to read the nz
 * tables. (The canonical declaration can't live in a public header because
 * it carries compiler-dependent struct packing decisions.) */
struct tp_nz_entry_view_ {
    short  ima;
    short  imb;
    short  imc;
    short  pad_;
    double cg;
};

struct tp_path_view_ {
    int                            i_a, i_b, i_c;
    int                            l_a, l_b, l_c;
    int                            d_a, d_b, d_c;
    int                            mult_u, mult_v, mult_w;
    int                            offset_a, offset_b, offset_c;
    int                            weight_offset;
    double                        *cg;
    struct tp_nz_entry_view_      *nz;
    int                            n_nz;
};

struct tp_descriptor_view_ {
    int                       mode;
    int                       num_paths;
    int                       num_weights;
    int                       a_dim, b_dim, c_dim;
    struct tp_path_view_     *paths;
};

void irrep_tp_apply_weighted_batch_flat_neon(const tp_descriptor_t *desc_opaque, size_t batch,
                                             const double *weights, const double *a_in,
                                             const double *b_in, double *c_out) {
    const struct tp_descriptor_view_ *desc = (const struct tp_descriptor_view_ *)desc_opaque;

    /* Zero c_out (bulk memset — faster than a per-element SIMD zero). */
    memset(c_out, 0, (size_t)desc->c_dim * batch * sizeof(double));

    for (int k = 0; k < desc->num_paths; ++k) {
        const struct tp_path_view_     *p = &desc->paths[k];
        double                          w = weights ? weights[k] : 1.0;
        int                             d_a = p->d_a, d_b = p->d_b, d_c = p->d_c;
        const struct tp_nz_entry_view_ *nz = p->nz;
        int                             n_nz = p->n_nz;

        for (int u = 0; u < p->mult_u; ++u) {
            const double *a_block = a_in + (size_t)(p->offset_a + u * d_a) * batch;
            const double *b_block = b_in + (size_t)(p->offset_b + u * d_b) * batch;
            /* */ double *c_block = c_out + (size_t)(p->offset_c + u * d_c) * batch;

            for (int e = 0; e < n_nz; ++e) {
                const double *a_row = a_block + (size_t)nz[e].ima * batch;
                const double *b_row = b_block + (size_t)nz[e].imb * batch;
                /* */ double *c_row = c_block + (size_t)nz[e].imc * batch;
                double        coef = w * nz[e].cg;
                float64x2_t   coef_v = vdupq_n_f64(coef);

                size_t bi = 0;
                for (; bi + 2 <= batch; bi += 2) {
                    float64x2_t a_v = vld1q_f64(a_row + bi);
                    float64x2_t b_v = vld1q_f64(b_row + bi);
                    float64x2_t c_v = vld1q_f64(c_row + bi);
                    /* c += coef · a · b  (two fmas over the lanes) */
                    c_v = vfmaq_f64(c_v, vmulq_f64(coef_v, a_v), b_v);
                    vst1q_f64(c_row + bi, c_v);
                }
                for (; bi < batch; ++bi)
                    c_row[bi] += coef * a_row[bi] * b_row[bi];
            }
        }
    }
}

#else
typedef int irrep_tp_apply_weighted_batch_flat_neon_disabled_t;
#endif
