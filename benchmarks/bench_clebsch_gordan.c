/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/clebsch_gordan.h>

int main(void) {
    const long iterations = 1000000;

    double t0 = irrep_bench_now_ns();
    double acc = 0.0;
    for (long i = 0; i < iterations; ++i) {
        acc += irrep_cg(1, 0, 1, 0, 1, 0);
    }
    double t1 = irrep_bench_now_ns();
    (void)acc;

    irrep_bench_report("cg_110_110_110", iterations, t1 - t0, 1);

    cg_table_t *t = irrep_cg_table_build(4, 4);
    irrep_cg_table_free(t);
    return 0;
}
