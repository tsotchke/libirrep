/* SPDX-License-Identifier: MIT */
/* Tensor-product benchmarks. Two workloads:
 *   tp_apply_empty   — dispatch overhead on a zero-path descriptor.
 *   tp_apply_nequip  — NequIP-like l_max = 3 descriptor, the shape that
 *                      dominates equivariant-NN inference; the sparse
 *                      CG inner kernel is the regression target.
 *   tp_apply_nequip_weighted — same with per-path scalar weights.
 */
#include "harness.h"
#include <irrep/tensor_product.h>
#include <irrep/multiset.h>

int main(void) {
    /* -------- empty descriptor: measures pure call overhead -------- */
    {
        irrep_multiset_t *a = irrep_multiset_new(0);
        irrep_multiset_t *b = irrep_multiset_new(0);
        irrep_multiset_t *c = irrep_multiset_new(0);
        tp_descriptor_t  *d = irrep_tp_build(a, b, c, NULL, 0);

        const long iterations = 100000;
        double a_in[1] = { 0 }, b_in[1] = { 0 }, c_out[1] = { 0 };
        double t0 = irrep_bench_now_ns();
        for (long i = 0; i < iterations; ++i) {
            irrep_tp_apply(d, a_in, b_in, c_out);
        }
        double t1 = irrep_bench_now_ns();
        irrep_bench_report("tp_apply_empty", iterations, t1 - t0, 1);

        irrep_tp_free(d);
        irrep_multiset_free(a);
        irrep_multiset_free(b);
        irrep_multiset_free(c);
    }

    /* -------- NequIP-like descriptor ------------------------------ *
     * a: 1x0e + 1x1o + 1x2e   (dim 9)
     * b: 1x0e + 1x1o + 1x2e + 1x3o   (dim 16)
     * c: 1x0e + 1x1o + 1x2e   (dim 9)
     * 13 auto-enumerated paths — the shape NequIP / MACE batch
     * inference spends most of its time in.                            */
    {
        irrep_multiset_t *a = irrep_multiset_parse("1x0e + 1x1o + 1x2e");
        irrep_multiset_t *b = irrep_multiset_parse("1x0e + 1x1o + 1x2e + 1x3o");
        irrep_multiset_t *c = irrep_multiset_parse("1x0e + 1x1o + 1x2e");
        int paths[64];
        int np = irrep_tp_enumerate_paths(a, b, c, paths, 64);
        tp_descriptor_t *d = irrep_tp_build(a, b, c, paths, np);

        double A[9]  = { 1, 0.5, 0.3, -0.2, 0.1, 0.4, -0.1, 0.6, -0.3 };
        double B[16] = { 1, 0.2, -0.3, 0.4, 0.5, -0.1, 0.2, 0.3, -0.4,
                         0.1, 0.2, -0.3, 0.5, -0.6, 0.7, -0.1 };
        double C[9]  = { 0 };
        double w[32] = { 0 };
        for (int i = 0; i < np; ++i) w[i] = 0.1 * (i + 1);

        const long iters = 500000;

        double t0 = irrep_bench_now_ns();
        for (long i = 0; i < iters; ++i) irrep_tp_apply(d, A, B, C);
        double t1 = irrep_bench_now_ns();
        irrep_bench_report("tp_apply_nequip", iters, t1 - t0, 1);

        double t2 = irrep_bench_now_ns();
        for (long i = 0; i < iters; ++i) irrep_tp_apply_weighted(d, w, A, B, C);
        double t3 = irrep_bench_now_ns();
        irrep_bench_report("tp_apply_nequip_weighted", iters, t3 - t2, 1);

        irrep_tp_free(d);
        irrep_multiset_free(a);
        irrep_multiset_free(b);
        irrep_multiset_free(c);
    }
    return 0;
}
