/* SPDX-License-Identifier: MIT */
/* Third N data point in the kagome-Heisenberg γ extraction series.
 *
 * Cluster: kagome 2×3 rectangular (N=18, 6 translations in p1 only —
 * the rectangular aspect ratio breaks p6mm point-group symmetry, but
 * translations still act). Sz=0 sector (popcount=9).
 *
 * This complements N=12 (p6mm 2×2 square cluster, high symmetry) and
 * N=27 (p6mm 3×3 square cluster) with an intermediate N at a DIFFERENT
 * cluster shape. For KH γ-scaling the shape-dependence matters — the
 * methodology paper should report γ across several cluster geometries.
 *
 *   make examples
 *   ./build/bin/kagome18_entanglement
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

#define N_SITES 18
#define D_FULL  (1LL << N_SITES)
#define POPCOUNT 9

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(void) {
    printf("=== Kagome N=18 (2×3, p1) Γ-point GS entanglement ===\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 3);
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
    /* p1 has only the trivial irrep. */
    double _Complex chi_ones[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu =
        irrep_sg_little_group_irrep_new(lg, chi_ones, 1);
    irrep_sg_heisenberg_sector_t *S =
        irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);
    printf("  sector (Γ, trivial): dim = %lld in %.2f s\n", sector_dim, now_sec() - t);

    t = now_sec();
    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.13 * sin(0.41 * (double)i) + I * 0.07 * cos(0.23 * (double)i);
    /* Request 10 eigvecs to identify the lowest F=+1 state — the p1
     * trivial sector on this cluster contains both singlets and triplets
     * (p1 has no symmetry to separate them). */
    int K = 10;
    double *eigs_all = malloc((size_t)K * sizeof(double));
    double _Complex *psi_all = malloc((size_t)K * sector_dim * sizeof(double _Complex));
    irrep_status_t rc = irrep_lanczos_eigvecs_reorth(
        irrep_sg_heisenberg_sector_apply, S, sector_dim, K, 200, seed, eigs_all, psi_all);
    if (rc != IRREP_OK) { fprintf(stderr, "Lanczos failed\n"); return 1; }
    printf("  Lanczos (200 iters, %d eigvecs): %.2f s\n", K, now_sec() - t);
    printf("    10 lowest eigenvalues:\n");
    for (int k = 0; k < K; ++k)
        printf("      E[%d] = %+.8f\n", k, eigs_all[k]);

    /* Identify the lowest F=+1 (singlet) state via the Lanczos eigenvectors
     * first (unfold+measure F for each), so the downstream γ computation
     * uses the SINGLET GS — directly comparable to N=12 p6mm Γ-A_1. */
    int order = irrep_space_group_order(G);
    double _Complex *pw_early = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, pw_early);
    uint64_t        *opos_early = malloc((size_t)order * sizeof(uint64_t));
    double _Complex *oamp_early = malloc((size_t)order * sizeof(double _Complex));
    double _Complex *psi_scan = malloc((size_t)D_FULL * sizeof(double _Complex));
    uint64_t mask_early = (1ULL << N_SITES) - 1;
    int singlet_k = -1;
    double singlet_F = 0;
    for (int kk = 0; kk < K; ++kk) {
        memset(psi_scan, 0, (size_t)D_FULL * sizeof(double _Complex));
        double _Complex *src = psi_all + (size_t)kk * sector_dim;
        long long n_reps_early = irrep_sg_rep_table_count(T);
        long long written_e = 0;
        for (long long kr = 0; kr < n_reps_early && written_e < sector_dim; ++kr) {
            uint64_t u = irrep_sg_rep_table_get(T, kr);
            int n_ent = 0;
            for (int g = 0; g < order; ++g) {
                uint64_t gx = irrep_space_group_apply_bits(G, g, u);
                int found = -1;
                for (int e = 0; e < n_ent; ++e)
                    if (opos_early[e] == gx) { found = e; break; }
                if (found >= 0) oamp_early[found] += pw_early[g];
                else { opos_early[n_ent] = gx; oamp_early[n_ent] = pw_early[g]; ++n_ent; }
            }
            double nn = 0.0;
            for (int e = 0; e < n_ent; ++e)
                nn += creal(oamp_early[e] * conj(oamp_early[e]));
            if (nn < 1e-20) continue;
            double inv = 1.0 / sqrt(nn);
            for (int e = 0; e < n_ent; ++e)
                psi_scan[opos_early[e]] += src[written_e] * inv * oamp_early[e];
            ++written_e;
        }
        double _Complex Fk = 0.0;
        for (long long s = 0; s < D_FULL; ++s)
            Fk += conj(psi_scan[s]) * psi_scan[s ^ mask_early];
        if (singlet_k < 0 && creal(Fk) > 0.5) { singlet_k = kk; singlet_F = creal(Fk); }
    }
    free(psi_scan); free(opos_early); free(oamp_early); free(pw_early);

    if (singlet_k < 0) {
        fprintf(stderr, "  no F=+1 state found among 10 lowest — increase K\n");
        return 1;
    }
    printf("  → selecting singlet GS at k=%d (E=%+.8f, F=%+.4f)\n",
           singlet_k, eigs_all[singlet_k], singlet_F);

    double eigs[4];
    double _Complex *psi_sector = malloc((size_t)sector_dim * sizeof(double _Complex));
    eigs[0] = eigs_all[singlet_k];
    memcpy(psi_sector, psi_all + (size_t)singlet_k * sector_dim,
           (size_t)sector_dim * sizeof(double _Complex));
    printf("    E_singlet / J = %+.8f   (E_singlet / N = %+.8f)\n",
           eigs[0], eigs[0] / (double)N_SITES);

    /* Residual check. */
    double _Complex *Hpsi = malloc((size_t)sector_dim * sizeof(double _Complex));
    irrep_sg_heisenberg_sector_apply(psi_sector, Hpsi, S);
    double res2 = 0.0;
    for (long long i = 0; i < sector_dim; ++i) {
        double _Complex r = Hpsi[i] - eigs[0] * psi_sector[i];
        res2 += creal(r * conj(r));
    }
    printf("    ||H·ψ − E·ψ|| = %.3e\n", sqrt(res2));
    free(Hpsi); free(seed);

    /* Unfold to dense 2^18 = 262144. */
    double _Complex *psi_full = calloc((size_t)D_FULL, sizeof(double _Complex));
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);

    t = now_sec();
    uint64_t        *orbit_pos = malloc((size_t)order * sizeof(uint64_t));
    double _Complex *orbit_amp = malloc((size_t)order * sizeof(double _Complex));
    long long n_reps = irrep_sg_rep_table_count(T);
    long long written = 0;
    for (long long k = 0; k < n_reps && written < sector_dim; ++k) {
        uint64_t u = irrep_sg_rep_table_get(T, k);
        int n_entries = 0;
        for (int g = 0; g < order; ++g) {
            uint64_t gx = irrep_space_group_apply_bits(G, g, u);
            int found = -1;
            for (int e = 0; e < n_entries; ++e)
                if (orbit_pos[e] == gx) { found = e; break; }
            if (found >= 0) orbit_amp[found] += w[g];
            else { orbit_pos[n_entries] = gx; orbit_amp[n_entries] = w[g]; ++n_entries; }
        }
        double nn = 0.0;
        for (int e = 0; e < n_entries; ++e)
            nn += creal(orbit_amp[e] * conj(orbit_amp[e]));
        if (nn < 1e-20) continue;
        double inv = 1.0 / sqrt(nn);
        double _Complex coef = psi_sector[written];
        for (int e = 0; e < n_entries; ++e)
            psi_full[orbit_pos[e]] += coef * inv * orbit_amp[e];
        ++written;
    }
    printf("  unfold to dense 2^%d: %.2f s\n", N_SITES, now_sec() - t);
    free(orbit_pos); free(orbit_amp); free(w); free(psi_sector);

    double full_norm = 0.0;
    for (long long s = 0; s < D_FULL; ++s)
        full_norm += creal(psi_full[s] * conj(psi_full[s]));
    printf("    ||ψ_full||² = %.10f\n", full_norm);

    /* Entanglement entropies. */
    printf("\n  S_A for contiguous bipartitions:\n");
    printf("  %3s   %12s\n", "|A|", "S_A");
    printf("  ---   ------------\n");
    for (int nA = 1; nA <= 9; ++nA) {
        int sites_A[18];
        for (int i = 0; i < nA; ++i) sites_A[i] = i;
        long long dim_A = 1LL << nA;
        double _Complex *rho_A = malloc((size_t)dim_A * dim_A * sizeof(double _Complex));
        irrep_partial_trace(N_SITES, 2, psi_full, sites_A, nA, rho_A);
        double S = irrep_entropy_vonneumann(rho_A, (int)dim_A);
        printf("  %3d   %+12.8f\n", nA, S);
        free(rho_A);
    }

    /* Three-region Kitaev-Preskill γ. */
    printf("\n  Three-region γ (|A|=|B|=|C|=3):\n");
    int A[3] = {0, 1, 2};
    int B[3] = {6, 7, 8};
    int C[3] = {12, 13, 14};
    int regions[7][18];
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
    printf("  %-8s  %12s\n", "region", "entropy");
    printf("  --------  ------------\n");
    for (int r = 0; r < 7; ++r) {
        long long dim_R = 1LL << sizes[r];
        double _Complex *rho_R = malloc((size_t)dim_R * dim_R * sizeof(double _Complex));
        irrep_partial_trace(N_SITES, 2, psi_full, regions[r], sizes[r], rho_R);
        S_R[r] = irrep_entropy_vonneumann(rho_R, (int)dim_R);
        printf("  %-8s  %+12.8f\n", labels[r], S_R[r]);
        free(rho_R);
    }
    double gamma = S_R[0] + S_R[1] + S_R[2] - S_R[3] - S_R[4] - S_R[5] + S_R[6];
    printf("\n  γ (Kitaev-Preskill) = %+.8f\n", gamma);
    printf("  (log 2 = %+.8f expected for gapped Z_2)\n", log(2.0));

    /* Spin-flip diagnostic on the lowest-energy state. */
    double _Complex Fexp = 0.0;
    uint64_t mask = (1ULL << N_SITES) - 1;
    for (long long s = 0; s < D_FULL; ++s)
        Fexp += conj(psi_full[s]) * psi_full[s ^ mask];
    printf("\n  Spin-flip ⟨ψ_0|F|ψ_0⟩ = %+.8f %+.8fi\n", creal(Fexp), cimag(Fexp));

    /* Scan all 10 eigenstates, identify F-sector of each. Report lowest
     * F=+1 state's energy — this is the singlet-sector ground state,
     * directly comparable to N=12 p6mm Γ-A_1 (which was F=+1). */
    printf("\n  F-eigenvalue per Lanczos state (filter for singlet sector):\n");
    printf("  %-3s  %13s  %13s\n", "k", "E_k", "⟨ψ_k|F|ψ_k⟩");
    printf("  ---  -------------  ------------\n");
    uint64_t        *opos = malloc((size_t)order * sizeof(uint64_t));
    double _Complex *oamp = malloc((size_t)order * sizeof(double _Complex));
    double _Complex *psi_k_full = malloc((size_t)D_FULL * sizeof(double _Complex));
    double _Complex *pw = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, pw);
    int lowest_f_plus = -1;
    for (int k = 0; k < K; ++k) {
        memset(psi_k_full, 0, (size_t)D_FULL * sizeof(double _Complex));
        long long written_k = 0;
        double _Complex *psi_k_sector = psi_all + (size_t)k * sector_dim;
        for (long long kr = 0; kr < n_reps && written_k < sector_dim; ++kr) {
            uint64_t u = irrep_sg_rep_table_get(T, kr);
            int n_entries = 0;
            for (int g = 0; g < order; ++g) {
                uint64_t gx = irrep_space_group_apply_bits(G, g, u);
                int found = -1;
                for (int e = 0; e < n_entries; ++e)
                    if (opos[e] == gx) { found = e; break; }
                if (found >= 0) oamp[found] += pw[g];
                else { opos[n_entries] = gx; oamp[n_entries] = pw[g]; ++n_entries; }
            }
            double nn = 0.0;
            for (int e = 0; e < n_entries; ++e)
                nn += creal(oamp[e] * conj(oamp[e]));
            if (nn < 1e-20) continue;
            double inv = 1.0 / sqrt(nn);
            double _Complex coef = psi_k_sector[written_k];
            for (int e = 0; e < n_entries; ++e)
                psi_k_full[opos[e]] += coef * inv * oamp[e];
            ++written_k;
        }
        double _Complex Fk = 0.0;
        for (long long s = 0; s < D_FULL; ++s)
            Fk += conj(psi_k_full[s]) * psi_k_full[s ^ mask];
        double Fk_re = creal(Fk);
        printf("  %-3d  %+.8f  %+.8f\n", k, eigs_all[k], Fk_re);
        if (lowest_f_plus < 0 && Fk_re > 0.5) lowest_f_plus = k;
    }
    if (lowest_f_plus >= 0) {
        printf("\n  → lowest F=+1 (singlet-sector) state: k=%d, E=%+.8f (E/N=%+.8f)\n",
               lowest_f_plus, eigs_all[lowest_f_plus],
               eigs_all[lowest_f_plus] / (double)N_SITES);
    } else {
        printf("\n  → no F=+1 state found in the first %d eigenvalues\n", K);
    }
    free(opos); free(oamp); free(psi_k_full); free(pw); free(eigs_all); free(psi_all);

    printf("\n  Total wall-clock: %.2f s\n", now_sec() - t0);

    free(psi_full);
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
