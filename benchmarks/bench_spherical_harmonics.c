/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/spherical_harmonics.h>

int main(void) {
    const int   l_max      = 4;
    const long  iterations = 100000;
    const int   out_dim    = (l_max + 1) * (l_max + 1);

    double out[(4 + 1) * (4 + 1)];
    double r_hat[3] = { 0.3, 0.4, 0.5 /* not unit; M1 stub doesn't care */ };

    double t0 = irrep_bench_now_ns();
    for (long i = 0; i < iterations; ++i) {
        irrep_sph_harm_cart_all(l_max, out, r_hat);
    }
    double t1 = irrep_bench_now_ns();

    irrep_bench_report("sph_harm_cart_all_l4", iterations, t1 - t0, out_dim);
    return 0;
}
