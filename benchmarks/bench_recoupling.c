/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/recoupling.h>

int main(void) {
    /* Wigner 6j evaluations across a range of small-j configurations
     * representative of atomic / NN-coupling use cases. */
    const long iterations = 200000;

    /* j=2 spin-orbit recoupling triple */
    {
        double t0 = irrep_bench_now_ns();
        double acc = 0.0;
        for (long i = 0; i < iterations; ++i)
            acc += irrep_wigner_6j(2, 2, 2, 2, 2, 2);
        double t1 = irrep_bench_now_ns();
        (void)acc;
        irrep_bench_report("wigner_6j_222222", iterations, t1 - t0, 1);
    }

    /* mixed j=4, j=2 — typical d-electron LS recoupling */
    {
        double t0 = irrep_bench_now_ns();
        double acc = 0.0;
        for (long i = 0; i < iterations; ++i)
            acc += irrep_wigner_6j(4, 4, 2, 2, 2, 4);
        double t1 = irrep_bench_now_ns();
        (void)acc;
        irrep_bench_report("wigner_6j_442224", iterations, t1 - t0, 1);
    }

    /* larger j=6 to exercise a longer t-loop */
    {
        double t0 = irrep_bench_now_ns();
        double acc = 0.0;
        for (long i = 0; i < iterations; ++i)
            acc += irrep_wigner_6j(6, 6, 4, 4, 4, 4);
        double t1 = irrep_bench_now_ns();
        (void)acc;
        irrep_bench_report("wigner_6j_664444", iterations, t1 - t0, 1);
    }

    /* 9j composed of three 6j calls — representative 4-spin recoupling */
    {
        double t0 = irrep_bench_now_ns();
        double acc = 0.0;
        for (long i = 0; i < iterations; ++i)
            acc += irrep_wigner_9j(2, 2, 2, 2, 2, 2, 2, 2, 0);
        double t1 = irrep_bench_now_ns();
        (void)acc;
        irrep_bench_report("wigner_9j_222222220", iterations, t1 - t0, 1);
    }

    return 0;
}
