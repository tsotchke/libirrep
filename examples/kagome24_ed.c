/* SPDX-License-Identifier: MIT */
/* glibc feature-test: CLOCK_MONOTONIC is POSIX.1b (1993), not in strict
 * C11. Define before <time.h> is pulled in by any header. */
#define _POSIX_C_SOURCE 199309L
/* 24-site kagome Heisenberg ED on a 2×4 torus, via the sparse Lanczos
 * solver in irrep/rdm.h. Hilbert space 2^24 = 16.78 M complex amplitudes
 * (≈ 256 MB per state vector). Demonstrates that the libirrep 1.3
 * substrate scales to the next cluster beyond 18-site without requiring
 * dense block-ED: Lanczos converges to the ground state in ~50 iterations
 * on an on-the-fly H-apply (never materialising H as a dense matrix).
 *
 * Memory footprint: 3 Lanczos state vectors (~768 MB) plus the H-apply
 * callback working buffer. Runtime on M2 Ultra: a few minutes depending
 * on iteration count and convergence.
 *
 *   make examples
 *   ./build/bin/kagome24_ed
 */

#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N   24
#define DIM (1LL << N)        /* 16 777 216 */

static double now_s_(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(void) {
    printf("===  24-site kagome Heisenberg ED via Lanczos (2×4 torus)  ===\n\n");

    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    int Nsites = irrep_lattice_num_sites(L);
    int nb     = irrep_lattice_num_bonds_nn(L);
    printf("cluster:  %d sites,  %d NN bonds\n", Nsites, nb);
    printf("Hilbert space: 2^%d = %lld  (%.1f MB per state vector)\n\n",
           N, DIM, (double)DIM * sizeof(double _Complex) / (1024.0*1024.0));

    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);

    /* Promoted to a library primitive; see include/irrep/hamiltonian.h. */
    irrep_heisenberg_t *H = irrep_heisenberg_new(Nsites, nb, bi, bj, /*J=*/1.0);
    if (!H) { fprintf(stderr, "heisenberg_new failed\n"); return 1; }

    /* Seed in the S_z = 0 sector with pseudo-random amplitudes. */
    printf("allocating & seeding 3 state vectors (%lld MB) ...\n",
           3LL * DIM * (long long)sizeof(double _Complex) / (1LL << 20));
    double _Complex *seed = malloc(sizeof(double _Complex) * DIM);
    if (!seed) { fprintf(stderr, "OOM\n"); return 1; }
    uint64_t rng = 0x13579bdf13579ULL;
    long long num_sz0 = 0;
    for (long long s = 0; s < DIM; ++s) {
        if (__builtin_popcountll(s) == N / 2) {
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            double val = (double)(rng >> 32) / (double)0xFFFFFFFFULL - 0.5;
            seed[s] = val;
            ++num_sz0;
        } else {
            seed[s] = 0.0;
        }
    }
    printf("S_z = 0 subspace dimension: %lld   (C(24, 12) = 2 704 156)\n\n",
           num_sz0);

    /* Run Lanczos. */
    const int k_wanted = 2;           /* ground + first excited in Sz=0 */
    const int max_iters = 80;         /* usually converges in < 50 */
    double ev[2] = {0, 0};

    printf("running Lanczos (max %d iterations, 2 eigenvalues wanted) ...\n",
           max_iters);
    double t0 = now_s_();
    irrep_status_t rc = irrep_lanczos_eigvals(
        irrep_heisenberg_apply, H, DIM, k_wanted, max_iters, seed, ev);
    double t1 = now_s_();
    if (rc != IRREP_OK) { fprintf(stderr, "Lanczos failed: %d\n", rc); return 1; }
    printf("  Lanczos: %.2f s  (%.2f s per iteration)\n",
           t1 - t0, (t1 - t0) / max_iters);

    double E0 = ev[0];
    double E1 = ev[1];
    printf("\nE_0              = %+.8f J\n", E0);
    printf("E_0 / N_site     = %+.8f J\n", E0 / Nsites);
    printf("    (24-site kagome ED literature: roughly −0.438 to −0.443 J/site\n"
           "     depending on cluster geometry; 2×4 torus published as ≈ −0.441)\n");
    printf("\nE_1              = %+.8f J\n", E1);
    printf("Δ (E_1 − E_0) in Sz=0 sector = %+.8f J\n", E1 - E0);

    /* Lowest triplet in Sz=1 sector. */
    printf("\nrunning Lanczos in S_z = 1 sector for lowest triplet ...\n");
    memset(seed, 0, sizeof(double _Complex) * DIM);
    rng = 0x24680ace24680ULL;
    for (long long s = 0; s < DIM; ++s) {
        if (__builtin_popcountll(s) == N / 2 - 1) {
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            seed[s] = (double)(rng >> 32) / (double)0xFFFFFFFFULL - 0.5;
        }
    }
    double ev_trip[1];
    t0 = now_s_();
    rc = irrep_lanczos_eigvals(irrep_heisenberg_apply, H, DIM,
                               /*k_wanted=*/1, max_iters, seed, ev_trip);
    t1 = now_s_();
    if (rc != IRREP_OK) { fprintf(stderr, "Lanczos Sz=1 failed: %d\n", rc); return 1; }
    double spin_gap = ev_trip[0] - E0;
    printf("  Lanczos: %.2f s\n", t1 - t0);
    printf("E_triplet         = %+.8f J\n", ev_trip[0]);
    printf("spin gap Δ_S      = %+.8f J\n", spin_gap);

    /* ------------------------------------------------------------------ *
     * Finite-size scaling summary.                                        *
     *                                                                    *
     *   Δ_S(1/N) = Δ_∞ + a · (1/N) + b · (1/N²)                            *
     *                                                                    *
     * Fit linear and quadratic models to the three libirrep ED values    *
     * at 12, 18, and this run's 24-site cluster. With only three points  *
     * we cannot distinguish a slow gapless flow from a finite-gap        *
     * saturation, but the linear extrapolation gives the cleanest        *
     * numerical estimate of the thermodynamic limit.                     *
     * ------------------------------------------------------------------ */
    {
        int    Ns[3]  = { 12, 18, 24 };
        double Ds[3]  = { 0.3827, 0.2835, spin_gap };   /* use this run's Δ */
        double Es[3]  = { -0.4537, -0.4471, E0 / Nsites };

        double invN[3], invN2[3];
        for (int k = 0; k < 3; ++k) {
            invN [k] = 1.0 / (double)Ns[k];
            invN2[k] = invN[k] * invN[k];
        }

        /* Least-squares linear fit  y = y∞ + a/N */
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int k = 0; k < 3; ++k) {
            sx  += invN[k];
            sy  += Ds[k];
            sxx += invN2[k];
            sxy += invN[k] * Ds[k];
        }
        double n = 3.0;
        double a_lin = (n * sxy - sx * sy) / (n * sxx - sx * sx);
        double Dinf_lin = (sy - a_lin * sx) / n;

        printf("\n==========================================================\n");
        printf("  finite-size scaling summary:  Δ_S(N) for kagome Heisenberg\n");
        printf("==========================================================\n");
        printf("  N_site | 1/N       | Δ_S (J)      | E_0/N (J)\n");
        printf("  -------+-----------+--------------+-----------\n");
        for (int k = 0; k < 3; ++k) {
            printf("  %6d | %.6f  | %+9.6f    | %+9.6f\n",
                   Ns[k], invN[k], Ds[k], Es[k]);
        }
        printf("  -------+-----------+--------------+-----------\n");
        printf("  linear 1/N extrapolation:\n");
        printf("    Δ_S(N→∞) ≈ %+.6f J    (slope a = %+.4f J)\n",
               Dinf_lin, a_lin);
        printf("    For reference: YHW/White DMRG 2011 (cylinder) → 0.13 J;\n");
        printf("    gapless Dirac spin liquid hypothesis → 0.\n");

        /* Ground-state energy per site extrapolation */
        sx = sy = sxx = sxy = 0;
        for (int k = 0; k < 3; ++k) {
            sx  += invN[k];
            sy  += Es[k];
            sxx += invN2[k];
            sxy += invN[k] * Es[k];
        }
        a_lin = (n * sxy - sx * sy) / (n * sxx - sx * sx);
        double Einf_lin = (sy - a_lin * sx) / n;
        printf("    (E_0/N)(N→∞) ≈ %+.6f J   (published thermodynamic limit ≈ −0.437)\n",
               Einf_lin);
    }

    free(seed);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_lattice_free(L);
    return 0;
}
