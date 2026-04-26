/* SPDX-License-Identifier: MIT */
/* Demonstrates the cache-aware workflow for iterative symmetry-resolved
 * ED research. Rep-table and sector-binding persistence turn the 13 s
 * N=27 rebuild into a 1 ms file read.
 *
 * Typical research iteration:
 *   1. First invocation: compute rep table + sector binding, save both.
 *   2. Subsequent invocations (change seed, parameter sweep, debug run):
 *      load from disk, run Lanczos.
 *
 *   make examples
 *   ./build/bin/ed_cache_demo         # first run: build + save
 *   ./build/bin/ed_cache_demo         # second run: load + apply
 *   rm -f /tmp/libirrep_kagome27_*    # reset cache
 *
 * Output shows timings per phase for both the cold-start and the
 * cached invocation — the difference is the user-visible win. */
#include <irrep/config_project.h>
#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>
#include <irrep/space_group.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* On-disk cache paths. In a real workflow these should go in a
 * per-project directory; /tmp is fine for demonstration. */
#define REP_CACHE_PATH    "/tmp/libirrep_kagome27_reps.bin"
#define SECTOR_CACHE_PATH "/tmp/libirrep_kagome27_gamma_a1.bin"

int main(void) {
    printf("=== libirrep cache-workflow demo — N=27 kagome ===\n\n");

    double t_total = now_sec();

    /* Lattice + Hamiltonian (cheap, built every run). */
    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
    int N  = irrep_space_group_num_sites(G);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N, nb, bi, bj, 1.0);

    /* Step 1: rep table. Try load; on miss, build + save. */
    double t = now_sec();
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_load(G, 13, REP_CACHE_PATH);
    int rep_from_cache = (T != NULL);
    if (!T) {
        T = irrep_sg_rep_table_build(G, 13);
        if (irrep_sg_rep_table_save(T, REP_CACHE_PATH) != 0)
            fprintf(stderr, "  (warning: rep-table save failed)\n");
    }
    double t_rep = now_sec() - t;
    printf("rep table %s in %.3f s (%lld reps)\n",
           rep_from_cache ? "LOADED from cache" : "BUILT + cached",
           t_rep, irrep_sg_rep_table_count(T));

    /* Step 2: Γ-A_1 sector binding. Try load; on miss, build + save.
     * The binding depends on (H, T) so the cache is keyed on the
     * Hamiltonian coupling + lattice; if the user changes J values
     * between runs they should remove the cache file. */
    t = now_sec();
    long long expected_dim = irrep_sg_rep_table_count(T);
    irrep_sg_heisenberg_sector_t *S =
        irrep_sg_heisenberg_sector_load(SECTOR_CACHE_PATH, expected_dim);
    int sec_from_cache = (S != NULL);
    if (!S) {
        S = irrep_sg_heisenberg_sector_build(H, T);
        if (irrep_sg_heisenberg_sector_save(S, SECTOR_CACHE_PATH) != 0)
            fprintf(stderr, "  (warning: sector save failed)\n");
    }
    double t_sec = now_sec() - t;
    printf("sector     %s in %.3f s (dim = %lld)\n",
           sec_from_cache ? "LOADED from cache" : "BUILT + cached",
           t_sec, irrep_sg_heisenberg_sector_dim(S));

    /* Step 3: Lanczos in the Γ-A_1 sector — the actual physics. */
    long long n = irrep_sg_heisenberg_sector_dim(S);
    double _Complex *seed = malloc((size_t)n * sizeof(double _Complex));
    for (long long i = 0; i < n; ++i)
        seed[i] = 0.13 * sin(0.41 * (double)i) + I * 0.07 * cos(0.23 * (double)i);
    double eigs[4];
    t = now_sec();
    irrep_status_t rc = irrep_lanczos_eigvals_reorth(
        irrep_sg_heisenberg_sector_apply, S, n, 4, 60, seed, eigs);
    double t_lanc = now_sec() - t;
    if (rc != IRREP_OK) {
        fprintf(stderr, "Lanczos failed (%d)\n", rc);
        return 1;
    }
    printf("lanczos    (60 iter, 4 eigs)     %.3f s\n", t_lanc);

    printf("\n  Γ-A_1 ground state:  E_0 / J = %+.8f    E_0 / N = %+.8f\n",
           eigs[0], eigs[0] / (double)N);

    double total = now_sec() - t_total;
    printf("\n  Total wall-clock: %.3f s\n", total);
    if (rep_from_cache && sec_from_cache) {
        printf("  (Both cached — subsequent runs will be this fast.)\n");
    } else {
        printf("  (Caches populated — next run will load instead of rebuild.)\n");
        printf("  Try running again: `./build/bin/ed_cache_demo`.\n");
    }

    free(seed);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
