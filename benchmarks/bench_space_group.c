/* SPDX-License-Identifier: MIT */
/* Benchmarks for site-permutation application.
 *
 * The 1.3 Target: < 1 µs per group action on the
 * 6×6 × 3 = 108-site kagome cluster so the full-orbit sum over the
 * p6mm order-432 group costs < 0.5 ms per MCMC step — subordinate to
 * the wavefunction evaluation. The bench measures three surfaces:
 *
 *   1. Raw permutation fetch: `irrep_space_group_apply(g, s)`.
 *   2. Full permutation write: `irrep_space_group_permutation(g, out)`.
 *   3. Configuration pullback: `irrep_space_group_apply_config` on a
 *      length-num_sites double buffer — the surface most callers use
 *      from inside the orbit-sum inner loop.
 *
 * All three are (ultimately) a `memcpy` by design — the permutation
 * table is pre-computed at build. The benchmark exists to catch any
 * future regression from added runtime bookkeeping.
 */

#include "harness.h"
#include <irrep/lattice.h>
#include <irrep/space_group.h>
#include <stdlib.h>
#include <string.h>

static void bench_full_orbit_(const char *name,
                              irrep_lattice_kind_t kind,
                              irrep_wallpaper_t wp,
                              int L,
                              long iterations) {
    irrep_lattice_t *lat = irrep_lattice_build(kind, L, L);
    if (!lat) return;
    irrep_space_group_t *G = irrep_space_group_build(lat, wp);
    if (!G) { irrep_lattice_free(lat); return; }

    int n     = irrep_space_group_num_sites(G);
    int order = irrep_space_group_order(G);

    double *cfg_in  = calloc((size_t)n, sizeof(double));
    double *cfg_out = calloc((size_t)n, sizeof(double));
    for (int s = 0; s < n; ++s) cfg_in[s] = 0.137 * (double)s - 0.5;

    /* Sum-over-group cost: the inner loop a symmetric-projection caller runs. */
    double t0 = irrep_bench_now_ns();
    for (long it = 0; it < iterations; ++it) {
        for (int g = 0; g < order; ++g) {
            irrep_space_group_apply_config(G, g, cfg_in, cfg_out);
        }
    }
    double t1 = irrep_bench_now_ns();
    char n_buf[96];
    snprintf(n_buf, sizeof(n_buf), "sg_apply_config_%s_L%d_per_element", name, L);
    /* ops_per_iter = order (one apply per group element, per iteration). */
    irrep_bench_report(n_buf, iterations, t1 - t0, order);

    /* Single-element permutation copy. */
    t0 = irrep_bench_now_ns();
    int *perm = malloc((size_t)n * sizeof(int));
    for (long it = 0; it < iterations * order; ++it) {
        int g = (int)(it % order);
        irrep_space_group_permutation(G, g, perm);
    }
    t1 = irrep_bench_now_ns();
    snprintf(n_buf, sizeof(n_buf), "sg_permutation_copy_%s_L%d", name, L);
    irrep_bench_report(n_buf, iterations * order, t1 - t0, 1);
    free(perm);

    /* Single-site apply. */
    long site_its = iterations * (long)order * (long)n;
    t0 = irrep_bench_now_ns();
    long probe = 0;
    for (long it = 0; it < site_its; ++it) {
        int g = (int)((it / n) % order);
        int s = (int)(it % n);
        probe += irrep_space_group_apply(G, g, s);
    }
    t1 = irrep_bench_now_ns();
    snprintf(n_buf, sizeof(n_buf), "sg_apply_site_%s_L%d", name, L);
    irrep_bench_report(n_buf, site_its, t1 - t0, 1);
    /* anti-DCE */
    if (probe == 0xBADFEE1L) fprintf(stderr, "probe=%ld\n", probe);

    free(cfg_in); free(cfg_out);
    irrep_space_group_free(G);
    irrep_lattice_free(lat);
}

int main(void) {
    /* Kagome 6×6 is the 6×6 target. 432 elements × 108 sites = 46,656
     * permutation entries per direction; full-orbit apply should take
     * << 1 ms at < 1 µs per element. */
    bench_full_orbit_("kagome_p6mm",    IRREP_LATTICE_KAGOME,    IRREP_WALLPAPER_P6MM, 6,  100);
    /* Triangular and square for comparison. */
    bench_full_orbit_("triangular_p6mm",IRREP_LATTICE_TRIANGULAR,IRREP_WALLPAPER_P6MM, 6,  500);
    bench_full_orbit_("square_p4mm",    IRREP_LATTICE_SQUARE,    IRREP_WALLPAPER_P4MM, 6, 1000);
    return 0;
}
