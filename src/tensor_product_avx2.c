/* SPDX-License-Identifier: MIT */
/* AVX2+FMA (x86_64) batch-flat kernel for `irrep_tp_apply_weighted_batch_flat`.
 * Four contiguous batch lanes per `__m256d`. Mirror of the NEON kernel. */

#if (defined(__x86_64__) || defined(_M_X64)) && defined(__AVX2__) && defined(__FMA__)

#include <immintrin.h>
#include <stddef.h>
#include <string.h>

#include <irrep/tensor_product.h>

#pragma STDC FP_CONTRACT OFF

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

void irrep_tp_apply_weighted_batch_flat_avx2(const tp_descriptor_t *desc_opaque, size_t batch,
                                             const double *weights, const double *a_in,
                                             const double *b_in, double *c_out) {
    const struct tp_descriptor_view_ *desc = (const struct tp_descriptor_view_ *)desc_opaque;
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
                __m256d       coef_v = _mm256_set1_pd(coef);

                size_t bi = 0;
                for (; bi + 4 <= batch; bi += 4) {
                    __m256d a_v = _mm256_loadu_pd(a_row + bi);
                    __m256d b_v = _mm256_loadu_pd(b_row + bi);
                    __m256d c_v = _mm256_loadu_pd(c_row + bi);
                    /* c += coef · a · b  (fmadd: c = a·b + c with coef folded) */
                    __m256d prod = _mm256_mul_pd(coef_v, a_v);
                    c_v = _mm256_fmadd_pd(prod, b_v, c_v);
                    _mm256_storeu_pd(c_row + bi, c_v);
                }
                for (; bi < batch; ++bi)
                    c_row[bi] += coef * a_row[bi] * b_row[bi];
            }
        }
    }
}

#else
typedef int irrep_tp_apply_weighted_batch_flat_avx2_disabled_t;
#endif
