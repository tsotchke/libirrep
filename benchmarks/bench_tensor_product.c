/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/tensor_product.h>
#include <irrep/multiset.h>

int main(void) {
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
    return 0;
}
