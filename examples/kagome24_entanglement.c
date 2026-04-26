/* SPDX-License-Identifier: MIT */
/* Fourth N data point in the kagome-Heisenberg γ scaling series.
 * Kagome 2×4 rectangular (N=24, p1 only).
 *
 * Sector dim ≈ 338k (C(24,12)/|G|=8). Lanczos + 10-eigvec recovery +
 * singlet filter + unfold to 2^24 = 16.8M dense state + 7 partial
 * traces = ~3-4 minutes wall-clock on M2 Ultra with OpenMP.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_entanglement
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

#define N_SITES  24
#define D_FULL   (1LL << N_SITES)
#define POPCOUNT 12

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void unfold_sector_to_dense(const irrep_space_group_t *G,
                                   const irrep_sg_rep_table_t *T,
                                   int order,
                                   const double _Complex *w,
                                   const double _Complex *psi_sector,
                                   long long sector_dim,
                                   double _Complex *psi_full) {
    memset(psi_full, 0, (size_t)D_FULL * sizeof(double _Complex));
    uint64_t        *opos = malloc((size_t)order * sizeof(uint64_t));
    double _Complex *oamp = malloc((size_t)order * sizeof(double _Complex));
    long long n_reps = irrep_sg_rep_table_count(T);
    long long written = 0;
    for (long long k = 0; k < n_reps && written < sector_dim; ++k) {
        uint64_t u = irrep_sg_rep_table_get(T, k);
        int n_ent = 0;
        for (int g = 0; g < order; ++g) {
            uint64_t gx = irrep_space_group_apply_bits(G, g, u);
            int found = -1;
            for (int e = 0; e < n_ent; ++e)
                if (opos[e] == gx) { found = e; break; }
            if (found >= 0) oamp[found] += w[g];
            else { opos[n_ent] = gx; oamp[n_ent] = w[g]; ++n_ent; }
        }
        double nn = 0.0;
        for (int e = 0; e < n_ent; ++e)
            nn += creal(oamp[e] * conj(oamp[e]));
        if (nn < 1e-20) continue;
        double inv = 1.0 / sqrt(nn);
        double _Complex coef = psi_sector[written];
        for (int e = 0; e < n_ent; ++e)
            psi_full[opos[e]] += coef * inv * oamp[e];
        ++written;
    }
    free(opos); free(oamp);
}

int main(void) {
    printf("=== Kagome N=24 (2×4, p1) singlet-sector γ extraction ===\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);

    double t = now_sec();
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);
    printf("  rep table: %lld reps in %.2f s\n", irrep_sg_rep_table_count(T), now_sec() - t);

    t = now_sec();
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
    double _Complex chi_ones[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu =
        irrep_sg_little_group_irrep_new(lg, chi_ones, 1);
    irrep_sg_heisenberg_sector_t *S =
        irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);
    printf("  sector build: dim = %lld in %.2f s\n", sector_dim, now_sec() - t);

    /* Lanczos, 10 eigvecs. */
    t = now_sec();
    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.1 * sin(0.37 * (double)i) + I * 0.05 * cos(0.23 * (double)i);
    const int K = 10;
    double          *eigs_all = malloc((size_t)K * sizeof(double));
    double _Complex *psi_all  = malloc((size_t)K * sector_dim * sizeof(double _Complex));
    irrep_status_t rc = irrep_lanczos_eigvecs_reorth(
        irrep_sg_heisenberg_sector_apply, S, sector_dim, K, 200, seed, eigs_all, psi_all);
    if (rc != IRREP_OK) { fprintf(stderr, "Lanczos failed\n"); return 1; }
    printf("  Lanczos (200 iters, %d eigvecs): %.2f s\n", K, now_sec() - t);

    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);

    /* Find lowest F=+1 (singlet). */
    double _Complex *psi_full = malloc((size_t)D_FULL * sizeof(double _Complex));
    uint64_t mask = (1ULL << N_SITES) - 1;
    printf("\n  F-scan of lowest %d eigenvalues:\n", K);
    printf("  %-3s  %13s  %13s\n", "k", "E_k", "⟨ψ|F|ψ⟩");
    int singlet_k = -1;
    for (int k = 0; k < K; ++k) {
        unfold_sector_to_dense(G, T, order, w,
                               psi_all + (size_t)k * sector_dim, sector_dim, psi_full);
        double _Complex Fk = 0.0;
        for (long long s = 0; s < D_FULL; ++s)
            Fk += conj(psi_full[s]) * psi_full[s ^ mask];
        double Fk_re = creal(Fk);
        printf("  %-3d  %+.8f  %+.8f\n", k, eigs_all[k], Fk_re);
        if (singlet_k < 0 && Fk_re > 0.5) singlet_k = k;
    }
    if (singlet_k < 0) { fprintf(stderr, "  no F=+1 state found\n"); return 1; }
    printf("\n  → singlet GS at k=%d: E = %+.8f   E/N = %+.8f\n",
           singlet_k, eigs_all[singlet_k], eigs_all[singlet_k] / (double)N_SITES);

    /* Unfold the singlet state. */
    unfold_sector_to_dense(G, T, order, w,
                           psi_all + (size_t)singlet_k * sector_dim,
                           sector_dim, psi_full);

    /* Single-site entropy (should be exactly ln 2 for singlet). */
    {
        int site0[1] = {0};
        double _Complex rho1[4];
        irrep_partial_trace(N_SITES, 2, psi_full, site0, 1, rho1);
        double S1 = irrep_entropy_vonneumann(rho1, 2);
        printf("  S_{|A|=1} = %+.10f   (ln 2 = %+.10f, Δ = %.2e)\n",
               S1, log(2.0), fabs(S1 - log(2.0)));
    }

    /* Kitaev-Preskill γ (|A|=|B|=|C|=3). */
    printf("\n  Three-region γ (|A|=|B|=|C|=3):\n");
    int A[3] = {0, 1, 2};
    int B[3] = {8, 9, 10};
    int C[3] = {16, 17, 18};
    int regions[7][24];
    int sizes[7];
    for (int i = 0; i < 3; ++i) {
        regions[0][i] = A[i]; regions[1][i] = B[i]; regions[2][i] = C[i];
    }
    sizes[0] = sizes[1] = sizes[2] = 3;
    for (int i = 0; i < 3; ++i) { regions[3][i]=A[i]; regions[3][i+3]=B[i]; } sizes[3]=6;
    for (int i = 0; i < 3; ++i) { regions[4][i]=A[i]; regions[4][i+3]=C[i]; } sizes[4]=6;
    for (int i = 0; i < 3; ++i) { regions[5][i]=B[i]; regions[5][i+3]=C[i]; } sizes[5]=6;
    for (int i = 0; i < 3; ++i) {
        regions[6][i]=A[i]; regions[6][i+3]=B[i]; regions[6][i+6]=C[i];
    }
    sizes[6] = 9;
    const char *labels[7] = {"S_A","S_B","S_C","S_{AB}","S_{AC}","S_{BC}","S_{ABC}"};
    double S_R[7];
    printf("  %-8s  %12s  %s\n", "region", "entropy", "[t]");
    printf("  --------  ------------  -----\n");
    for (int r = 0; r < 7; ++r) {
        long long dim_R = 1LL << sizes[r];
        double _Complex *rho_R = malloc((size_t)dim_R * dim_R * sizeof(double _Complex));
        double tr = now_sec();
        irrep_partial_trace(N_SITES, 2, psi_full, regions[r], sizes[r], rho_R);
        S_R[r] = irrep_entropy_vonneumann(rho_R, (int)dim_R);
        printf("  %-8s  %+12.8f  [%.2f s]\n", labels[r], S_R[r], now_sec() - tr);
        free(rho_R);
    }
    double gamma = S_R[0] + S_R[1] + S_R[2] - S_R[3] - S_R[4] - S_R[5] + S_R[6];
    printf("\n  γ (Kitaev-Preskill) = %+.8f\n", gamma);
    printf("  (log 2 = %+.8f expected for gapped Z_2)\n", log(2.0));

    printf("\n  Total wall-clock: %.2f s\n", now_sec() - t0);

    free(psi_full); free(w); free(seed); free(eigs_all); free(psi_all);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
