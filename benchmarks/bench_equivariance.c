/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/so3.h>
#include <irrep/tensor_product.h>
#include <irrep/multiset.h>

int main(void) {
    /* M1 stub — real equivariance benchmark composes random rotation,
     * rotation of input feature, tp_apply, rotation of output, then measures
     * agreement cost. Wire up in M7/M14. */
    const long iterations = 10000;

    double t0 = irrep_bench_now_ns();
    for (long i = 0; i < iterations; ++i) {
        (void)irrep_quat_identity();
    }
    double t1 = irrep_bench_now_ns();

    irrep_bench_report("equivariance_stub", iterations, t1 - t0, 1);
    return 0;
}
