/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/spherical_harmonics.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
    const int   l_max   = 4;
    const int   out_dim = (l_max + 1) * (l_max + 1);

    /* Single-call throughput. */
    {
        const long iterations = 100000;
        double out[(4 + 1) * (4 + 1)];
        double r_hat[3] = { 0.3, 0.4, 0.5 };

        double t0 = irrep_bench_now_ns();
        for (long i = 0; i < iterations; ++i) {
            irrep_sph_harm_cart_all(l_max, out, r_hat);
        }
        double t1 = irrep_bench_now_ns();
        irrep_bench_report("sph_harm_cart_all_l4", iterations, t1 - t0, out_dim);
    }

    /* Batched throughput across a NequIP-representative edge count. Drives
     * the NEON batch kernel when available; the 1-ns min-abs gate in
     * perf_compare ignores the noise on very short single-edge calls. */
    {
        const size_t N = 4096;
        const long   iterations = 200;
        double *rhats = malloc(N * 3 * sizeof(double));
        double *outs  = malloc(N * (size_t)out_dim * sizeof(double));
        for (size_t i = 0; i < N; ++i) {
            double t = 0.17 * (double)i + 0.05;
            double x = sin(t), y = cos(0.7 * t), z = sin(1.3 * t);
            double n = sqrt(x*x + y*y + z*z);
            rhats[3*i+0] = x / n;
            rhats[3*i+1] = y / n;
            rhats[3*i+2] = z / n;
        }

        double t0 = irrep_bench_now_ns();
        for (long i = 0; i < iterations; ++i) {
            irrep_sph_harm_cart_all_batch(l_max, N, rhats, outs);
        }
        double t1 = irrep_bench_now_ns();
        irrep_bench_report("sph_harm_cart_all_batch_l4_4096E",
                           iterations, t1 - t0, (long)N);

        free(rhats); free(outs);
    }
    return 0;
}
