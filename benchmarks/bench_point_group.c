/* SPDX-License-Identifier: MIT */
/* Projection throughput across all four supported point groups, on a
 * representative multiset (2x0e + 1x1o + 1x2e — matches the small NequIP
 * hidden-space shape).
 *
 * Per-call cost scales roughly as `|G| · Σ_l (2l+1)²` since each element
 * invokes a (2l+1)×(2l+1) real-D multiply per block. The dominant cost
 * for small multisets is the `irrep_euler_zyz_from_rot + wigner_D_matrix`
 * construction per (element, l), run once per call — cacheable in a
 * future revision. */
#include "harness.h"
#include <irrep/multiset.h>
#include <irrep/point_group.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void bench_one_(const char *name_prefix, irrep_point_group_t g,
                       long iterations) {
    irrep_pg_table_t *t = irrep_pg_table_build(g);
    if (!t) return;
    irrep_multiset_t *m = irrep_multiset_parse("2x0e + 1x1o + 1x2e");
    if (!m) { irrep_pg_table_free(t); return; }
    int N = m->total_dim;

    double *in  = calloc((size_t)N, sizeof(double));
    double *out = calloc((size_t)N, sizeof(double));
    for (int i = 0; i < N; ++i) in[i] = 0.137 * (double)i - 0.5;

    /* Project onto the trivial irrep (μ = 0). */
    double t0 = irrep_bench_now_ns();
    for (long i = 0; i < iterations; ++i) {
        irrep_pg_project(t, /*mu=*/0, m, in, out);
    }
    double t1 = irrep_bench_now_ns();
    char n[96];
    snprintf(n, sizeof(n), "pg_project_%s_trivial", name_prefix);
    irrep_bench_report(n, iterations, t1 - t0, 1);

    free(in); free(out);
    irrep_multiset_free(m);
    irrep_pg_table_free(t);
}

int main(void) {
    bench_one_("c4v", IRREP_PG_C4V, 20000);
    bench_one_("d6",  IRREP_PG_D6,  20000);
    bench_one_("c3v", IRREP_PG_C3V, 20000);
    bench_one_("d3",  IRREP_PG_D3,  20000);
    return 0;
}
