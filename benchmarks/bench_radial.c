/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/radial.h>

int main(void) {
    const long iterations = 1000000;
    double acc = 0.0;

    double t0 = irrep_bench_now_ns();
    for (long i = 0; i < iterations; ++i) {
        acc += irrep_rbf_bessel(3, 0.75, 1.0);
    }
    double t1 = irrep_bench_now_ns();
    (void)acc;

    irrep_bench_report("rbf_bessel", iterations, t1 - t0, 1);
    return 0;
}
