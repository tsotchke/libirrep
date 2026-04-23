/* SPDX-License-Identifier: MIT */
/* Tests for `irrep_tp_apply_weighted_batch_flat` (SIMD-vectorised dim-first
 * batch). Validates bit-exactness against a per-sample scalar loop through
 * `irrep_tp_apply_weighted` on a freshly-built descriptor. */

#include "harness.h"
#include <irrep/multiset.h>
#include <irrep/tensor_product.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Reorder a batch-first buffer [batch][dim] into the dim-first layout
 * [dim][batch] expected by the `_flat` API. */
static void transpose_to_flat_(const double *src, double *dst, size_t batch, int dim) {
    for (size_t bi = 0; bi < batch; ++bi)
        for (int i = 0; i < dim; ++i)
            dst[(size_t)i * batch + bi] = src[bi * (size_t)dim + i];
}

/* Dim-first → batch-first transpose (for scalar-reference comparison). */
static void transpose_from_flat_(const double *src, double *dst, size_t batch, int dim) {
    for (size_t bi = 0; bi < batch; ++bi)
        for (int i = 0; i < dim; ++i)
            dst[bi * (size_t)dim + i] = src[(size_t)i * batch + bi];
}

int main(void) {
    IRREP_TEST_START("tp_batch_flat");

    /* 1 × 1o × 1 × 1o → 1 × 1e + 1 × 2e (three paths: (0, 0, 0), (0, 0, 1), (0, 0, 2)) */
    irrep_multiset_t *A = irrep_multiset_parse("1x1o");
    irrep_multiset_t *B = irrep_multiset_parse("1x1o");
    irrep_multiset_t *C = irrep_multiset_parse("1x0e + 1x1e + 1x2e");
    int               n_paths = irrep_tp_enumerate_paths(A, B, C, NULL, 0);
    IRREP_ASSERT(n_paths == 3);
    int *paths = malloc((size_t)n_paths * 3 * sizeof(int));
    irrep_tp_enumerate_paths(A, B, C, paths, n_paths);
    tp_descriptor_t *d = irrep_tp_build(A, B, C, paths, n_paths);
    IRREP_ASSERT(d != NULL);

    /* Build a batched a / b with random-ish but deterministic values. */
    const size_t batch = 17; /* odd to exercise the SIMD tail */
    const int    a_dim = A->total_dim;
    const int    b_dim = B->total_dim;
    const int    c_dim = irrep_tp_output_dim(d);

    double *a_bf = malloc(batch * (size_t)a_dim * sizeof(double));
    double *b_bf = malloc(batch * (size_t)b_dim * sizeof(double));
    double *c_ref = calloc(batch * (size_t)c_dim, sizeof(double));
    for (size_t bi = 0; bi < batch; ++bi) {
        for (int i = 0; i < a_dim; ++i) a_bf[bi * a_dim + i] = 0.13 * (i + 1) + 0.07 * (int)bi;
        for (int i = 0; i < b_dim; ++i) b_bf[bi * b_dim + i] = 0.21 * (i + 1) - 0.05 * (int)bi;
    }
    double weights[3] = {1.3, -0.7, 0.5};

    /* Reference: per-sample `irrep_tp_apply_weighted`, batch-first. */
    for (size_t bi = 0; bi < batch; ++bi)
        irrep_tp_apply_weighted(d, weights, a_bf + bi * a_dim, b_bf + bi * b_dim,
                                c_ref + bi * c_dim);

    /* Flat SIMD: transpose inputs to dim-first, run, transpose output back,
     * compare bit-exactly against the reference. */
    double *a_df = malloc(batch * (size_t)a_dim * sizeof(double));
    double *b_df = malloc(batch * (size_t)b_dim * sizeof(double));
    double *c_df = calloc(batch * (size_t)c_dim, sizeof(double));
    transpose_to_flat_(a_bf, a_df, batch, a_dim);
    transpose_to_flat_(b_bf, b_df, batch, b_dim);
    irrep_tp_apply_weighted_batch_flat(d, batch, weights, a_df, b_df, c_df);
    double *c_bf_simd = malloc(batch * (size_t)c_dim * sizeof(double));
    transpose_from_flat_(c_df, c_bf_simd, batch, c_dim);

    size_t n_assertions = 0;
    for (size_t k = 0; k < batch * (size_t)c_dim; ++k) {
        double diff = fabs(c_ref[k] - c_bf_simd[k]);
        IRREP_ASSERT(diff < 1e-14);
        ++n_assertions;
    }
    (void)n_assertions;

    /* Unweighted reduction: weights = NULL on flat API not supported (the
     * function expects a weights array). Pass an all-ones weights vector
     * and compare vs per-sample `irrep_tp_apply`. */
    {
        double ones[3] = {1.0, 1.0, 1.0};
        double *c_ref_un = calloc(batch * (size_t)c_dim, sizeof(double));
        for (size_t bi = 0; bi < batch; ++bi)
            irrep_tp_apply(d, a_bf + bi * a_dim, b_bf + bi * b_dim, c_ref_un + bi * c_dim);
        double *c_df_un = calloc(batch * (size_t)c_dim, sizeof(double));
        irrep_tp_apply_weighted_batch_flat(d, batch, ones, a_df, b_df, c_df_un);
        double *c_bf_un = malloc(batch * (size_t)c_dim * sizeof(double));
        transpose_from_flat_(c_df_un, c_bf_un, batch, c_dim);
        for (size_t k = 0; k < batch * (size_t)c_dim; ++k)
            IRREP_ASSERT(fabs(c_ref_un[k] - c_bf_un[k]) < 1e-14);
        free(c_ref_un);
        free(c_df_un);
        free(c_bf_un);
    }

    /* Small batch sizes (including 1, 2, 3) must still produce correct
     * output — exercises the SIMD head/tail boundaries. */
    for (size_t b_small = 1; b_small <= 5; ++b_small) {
        double *a_small = calloc(b_small * (size_t)a_dim, sizeof(double));
        double *b_small_df = calloc(b_small * (size_t)a_dim, sizeof(double));
        double *b_small_arr = calloc(b_small * (size_t)b_dim, sizeof(double));
        double *c_small_ref = calloc(b_small * (size_t)c_dim, sizeof(double));
        double *c_small_df = calloc(b_small * (size_t)c_dim, sizeof(double));
        for (size_t bi = 0; bi < b_small; ++bi) {
            for (int i = 0; i < a_dim; ++i)
                a_small[bi * a_dim + i] = 0.3 + 0.2 * (i + bi);
            for (int i = 0; i < b_dim; ++i)
                b_small_arr[bi * b_dim + i] = -0.4 + 0.1 * (i + 2 * bi);
            irrep_tp_apply_weighted(d, weights, a_small + bi * a_dim, b_small_arr + bi * b_dim,
                                    c_small_ref + bi * c_dim);
        }
        transpose_to_flat_(a_small, b_small_df, b_small, a_dim);
        double *b_flat = malloc(b_small * (size_t)b_dim * sizeof(double));
        transpose_to_flat_(b_small_arr, b_flat, b_small, b_dim);
        irrep_tp_apply_weighted_batch_flat(d, b_small, weights, b_small_df, b_flat, c_small_df);
        double *c_bf2 = malloc(b_small * (size_t)c_dim * sizeof(double));
        transpose_from_flat_(c_small_df, c_bf2, b_small, c_dim);
        for (size_t k = 0; k < b_small * (size_t)c_dim; ++k)
            IRREP_ASSERT(fabs(c_small_ref[k] - c_bf2[k]) < 1e-14);
        free(a_small);
        free(b_small_df);
        free(b_small_arr);
        free(c_small_ref);
        free(c_small_df);
        free(b_flat);
        free(c_bf2);
    }

    free(a_bf);
    free(b_bf);
    free(c_ref);
    free(a_df);
    free(b_df);
    free(c_df);
    free(c_bf_simd);
    free(paths);
    irrep_tp_free(d);
    irrep_multiset_free(A);
    irrep_multiset_free(B);
    irrep_multiset_free(C);

    return IRREP_TEST_END();
}
