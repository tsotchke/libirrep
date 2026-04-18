/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/wigner_d.h>

int main(void) {
    const int  j          = 4;
    const long iterations = 100000;
    double _Complex D[(2 * 4 + 1) * (2 * 4 + 1)];

    double t0 = irrep_bench_now_ns();
    for (long i = 0; i < iterations; ++i) {
        irrep_wigner_D_matrix(j, D, 0.1, 0.2, 0.3);
    }
    double t1 = irrep_bench_now_ns();

    irrep_bench_report("wigner_D_matrix_j4", iterations, t1 - t0,
                       (2 * j + 1) * (2 * j + 1));
    return 0;
}
