/* SPDX-License-Identifier: MIT */
#ifndef LIBIRREP_BENCH_HARNESS_H
#define LIBIRREP_BENCH_HARNESS_H

#include <stdio.h>
#include <time.h>

static inline double irrep_bench_now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

/* Emit one JSON record to stdout. Caller provides ops-per-call scaling. */
static inline void irrep_bench_report(const char *name,
                                       long iterations,
                                       double total_ns,
                                       long ops_per_iter) {
    double ns_per_op   = total_ns / ((double)iterations * (double)ops_per_iter);
    double ops_per_sec = 1e9 / ns_per_op;
    printf("{\"name\":\"%s\",\"iterations\":%ld,\"ops_per_iter\":%ld,"
           "\"ns_per_op\":%.3f,\"ops_per_sec\":%.3e}\n",
           name, iterations, ops_per_iter, ns_per_op, ops_per_sec);
}

#endif /* LIBIRREP_BENCH_HARNESS_H */
