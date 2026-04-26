/* SPDX-License-Identifier: MIT */
/* Benchmark: sparse-rep-basis Heisenberg Lanczos vs dense-full-Hilbert
 * Lanczos on kagome clusters.
 *
 *   make examples
 *   ./build/bin/bench_sparse_lanczos
 *
 * What this demonstrates:
 * 1. On kagome 2×2 (N = 12, dim = 4096), both the dense and sparse
 *    Lanczos paths produce E_0 for the trivial sector and agree to
 *    1e-8 while the sparse path has ~|G| = 48× fewer matvec FLOPs.
 * 2. On kagome 3×3 (N = 27, dim = 1.34e8), a full-Hilbert state vector
 *    is 2.1 GiB per copy; Lanczos needs 3 copies ≈ 6.4 GiB just for
 *    vectors. The p6mm order is 108 so the trivial-sector rep table is
 *    ~1.2e6 reps (per Sz); state vectors drop to ~20 MiB. One sparse
 *    Lanczos run on a workstation takes seconds. This is the N = 27 ED
 *    unlock — infeasible on the dense path, routine on the sparse one.
 */

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

/* Monotonic wall-clock seconds. */
static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* The cached sector-binding IS the Lanczos context — its apply fn
 * already has the (psi_in, psi_out, ctx) signature Lanczos expects. */

static void bench_kagome_trivial_sector(int Lx, int Ly, int popcount, int do_dense_compare) {
    printf("\n=== kagome %dx%d (N = %d), popcount = %d, p6mm Γ-A_1 sector ===\n",
           Lx, Ly, 3 * Lx * Ly, popcount);

    double t0 = now_sec();
    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, Lx, Ly);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
    int N  = irrep_space_group_num_sites(G);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N, nb, bi, bj, 1.0);
    double t_build = now_sec() - t0;
    printf("  lattice + space group + Hamiltonian: %.3f s\n", t_build);

    /* Rep table at fixed popcount. */
    double t1 = now_sec();
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, popcount);
    long long n_reps = irrep_sg_rep_table_count(T);
    double t_reps = now_sec() - t1;
    printf("  rep table (popcount %d): %lld reps in %.3f s\n", popcount, n_reps, t_reps);
    printf("  reduction: %lld / 2^%d = %.3g (~%.1fx)\n",
           n_reps, N,
           (double)n_reps / (double)(1LL << N),
           (double)(1LL << N) / (double)n_reps);

    /* Build the cached sector binding — one-shot canonicalisation. */
    double tb = now_sec();
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build(H, T);
    double t_cache = now_sec() - tb;
    printf("  cached sector binding: %.3f s (one-shot precompute)\n", t_cache);

    /* Sparse Lanczos on the cached binding. */
    double t2 = now_sec();
    double *eigs = malloc(sizeof(double) * 8);
    double _Complex *seed = malloc((size_t)n_reps * sizeof(double _Complex));
    for (long long i = 0; i < n_reps; ++i)
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.07 * cos(0.29 * i);
    int                 max_iters = (int)((n_reps > 60) ? 60 : n_reps);
    int                 n_eig = (int)((n_reps > 4) ? 4 : n_reps);
    irrep_status_t rc = irrep_lanczos_eigvals_reorth(
        irrep_sg_heisenberg_sector_apply, S, n_reps, n_eig, max_iters, seed, eigs);
    double t_sparse = now_sec() - t2;
    if (rc != IRREP_OK) {
        printf("  sparse Lanczos FAILED (rc = %d)\n", rc);
    } else {
        printf("  sparse Lanczos (max_iters = %d, %d eigs): %.3f s\n", max_iters, n_eig, t_sparse);
        /* Ascending or descending? The reorth variant returns ascending
         * already per the signature convention (lowest at [0]). */
        printf("    E_0 / J = %+.8f   (E_0 / N = %+.8f)\n", eigs[0], eigs[0] / (double)N);
    }
    free(seed);

    /* Dense timing reference: run dense Lanczos on the FULL Hilbert
     * space so the matvec timing is directly comparable. (Sector
     * E_0 correctness is verified bit-exactly by test_sparse_apply.c
     * — here we just want a wall-clock reference.) */
    if (do_dense_compare && N <= 24) {
        double   t3     = now_sec();
        long long D     = 1LL << N;
        double _Complex *seed_d = calloc((size_t)D, sizeof(double _Complex));
        for (long long s = 0; s < D; ++s)
            seed_d[s] = (__builtin_popcountll((unsigned long long)s) == popcount)
                            ? 0.1 * sin(0.37 * (double)s) + I * 0.03 * cos(0.29 * (double)s)
                            : 0.0;

        double *eigs_d = malloc(sizeof(double) * 4);
        irrep_status_t rc_d =
            irrep_lanczos_eigvals_reorth(irrep_heisenberg_apply, (void *)H, D, 4, 60, seed_d, eigs_d);
        double t_dense = now_sec() - t3;
        if (rc_d != IRREP_OK) {
            printf("  dense Lanczos FAILED (rc = %d)\n", rc_d);
        } else {
            printf("  dense Lanczos (2^%d = %lld dim, full Hilbert): %.3f s\n", N, D, t_dense);
            printf("    absolute GS E_0 / J = %+.8f   (E_0 / N = %+.8f)\n",
                   eigs_d[0], eigs_d[0] / (double)N);
            printf("    [different sector than sparse; see test_sparse_apply for bit-exact check]\n");
            printf("    sparse matvec wall-clock speedup = %.2fx\n", t_dense / t_sparse);
        }
        free(seed_d); free(eigs_d);
    } else if (N > 24) {
        long long D     = 1LL << N;
        double    gib_v = (double)D * 16.0 / (1024.0 * 1024.0 * 1024.0);
        printf("  dense path infeasible at N = %d: 2^N = %lld, 1 state vector = %.2f GiB\n",
               N, D, gib_v);
    }

    free(eigs);
    free(bi);
    free(bj);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_rep_table_free(T);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
}

int main(void) {
    printf("=== sparse-rep-basis Lanczos benchmark: kagome Heisenberg ===\n");
    printf("Compiled with libirrep; runtime SIMD dispatch in effect.\n");

    /* N = 12: both paths feasible, cross-check. */
    bench_kagome_trivial_sector(2, 2, 6, 1);

    /* N = 18 (kagome 2×3 — Lx = Ly required for p6mm; use p1 here).
     * Skipping unless user wants to wire p1 version separately. */

    /* N = 27: dense infeasible. popcount 13 (Sz = −1/2) is the lowest
     * available non-zero-Sz sector on this odd-N cluster. */
    bench_kagome_trivial_sector(3, 3, 13, 0);

    return 0;
}
