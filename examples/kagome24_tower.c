/* SPDX-License-Identifier: MIT */
/* N=24 kagome 2×4 topological GS-tower check.
 *
 * A Z₂ spin liquid on a torus has 4 exponentially-near-degenerate
 * ground states (one per topological sector = (e-charge, m-flux)
 * combination). This is a sector-free, region-free, γ-free signature
 * of topological order. On the 2×4 kagome torus we enumerate the
 * lowest 20 F=+1 singlet states across all 8 Bloch momenta and look
 * for a tower of 4 closely-spaced states well-separated from the
 * rest.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_tower
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

#define N_SITES 24
#define D_FULL  (1LL << N_SITES)
#define POPCOUNT 12

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void unfold(const irrep_space_group_t *G,
                   const irrep_sg_rep_table_t *T,
                   int order, const double _Complex *w,
                   const double _Complex *psi_sector, long long sector_dim,
                   double _Complex *psi_full) {
    memset(psi_full, 0, (size_t)D_FULL * sizeof(double _Complex));
    uint64_t        *op = malloc((size_t)order * sizeof(uint64_t));
    double _Complex *oa = malloc((size_t)order * sizeof(double _Complex));
    long long n_reps = irrep_sg_rep_table_count(T);
    long long written = 0;
    for (long long k = 0; k < n_reps && written < sector_dim; ++k) {
        uint64_t u = irrep_sg_rep_table_get(T, k);
        int ne = 0;
        for (int g = 0; g < order; ++g) {
            uint64_t gx = irrep_space_group_apply_bits(G, g, u);
            int f = -1;
            for (int e = 0; e < ne; ++e) if (op[e] == gx) { f = e; break; }
            if (f >= 0) oa[f] += w[g];
            else { op[ne] = gx; oa[ne] = w[g]; ++ne; }
        }
        double nn = 0.0;
        for (int e = 0; e < ne; ++e) nn += creal(oa[e] * conj(oa[e]));
        if (nn < 1e-20) continue;
        double inv = 1.0 / sqrt(nn);
        double _Complex c = psi_sector[written];
        for (int e = 0; e < ne; ++e) psi_full[op[e]] += c * inv * oa[e];
        ++written;
    }
    free(op); free(oa);
}

typedef struct {
    double E;
    double F;
    int    kx, ky;
    int    level;   /* Lanczos index within its sector */
} state_t;

static int cmp_state(const void *a, const void *b) {
    double ea = ((const state_t *)a)->E;
    double eb = ((const state_t *)b)->E;
    return (ea < eb) ? -1 : (ea > eb) ? 1 : 0;
}

int main(void) {
    printf("=== N=24 kagome 2×4 — topological GS-tower scan ===\n");
    printf("    Looking for 4 nearly-degenerate F=+1 singlet states\n");
    printf("    (Z₂ topological signature on torus)\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);

    const int K_PER_SECTOR = 20;
    const int Lx = 2, Ly = 4;
    const int N_SECTORS = Lx * Ly;

    state_t *all_states = malloc(sizeof(state_t) * K_PER_SECTOR * N_SECTORS);
    int n_states = 0;

    printf("  Scanning %d k-sectors × %d states = %d spectrum samples …\n\n",
           N_SECTORS, K_PER_SECTOR, N_SECTORS * K_PER_SECTOR);

    double _Complex *psi_full = malloc((size_t)D_FULL * sizeof(double _Complex));
    uint64_t mask = (1ULL << N_SITES) - 1;

    for (int my = 0; my < Ly; ++my)
        for (int mx = 0; mx < Lx; ++mx) {
            double tk = now_sec();
            irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, mx, my);
            double _Complex chi1[1] = {1.0 + 0.0 * I};
            irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
            irrep_sg_heisenberg_sector_t *S =
                irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
            long long sdim = irrep_sg_heisenberg_sector_dim(S);

            double _Complex *seed = malloc((size_t)sdim * sizeof(double _Complex));
            for (long long i = 0; i < sdim; ++i)
                seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
            int K = K_PER_SECTOR > sdim ? (int)sdim : K_PER_SECTOR;
            double *eigs = malloc((size_t)K * sizeof(double));
            double _Complex *psi_all =
                malloc((size_t)K * sdim * sizeof(double _Complex));
            int it = 300 > sdim ? (int)sdim : 300;
            irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                         sdim, K, it, seed, eigs, psi_all);

            int order = irrep_space_group_order(G);
            double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
            irrep_sg_projector_weights(lg, mu, w);

            int singlets_found = 0;
            for (int k = 0; k < K; ++k) {
                unfold(G, T, order, w, psi_all + (size_t)k * sdim,
                       sdim, psi_full);
                double _Complex Fk = 0;
                for (long long s = 0; s < D_FULL; ++s)
                    Fk += conj(psi_full[s]) * psi_full[s ^ mask];
                if (creal(Fk) > 0.5) {
                    all_states[n_states].E  = eigs[k];
                    all_states[n_states].F  = creal(Fk);
                    all_states[n_states].kx = mx;
                    all_states[n_states].ky = my;
                    all_states[n_states].level = k;
                    ++n_states;
                    ++singlets_found;
                }
            }
            printf("  k=(%d,%d) dim=%-6lld  (%d singlets in lowest %d, %.1fs)\n",
                   mx, my, sdim, singlets_found, K, now_sec() - tk);

            free(seed); free(eigs); free(psi_all); free(w);
            irrep_sg_heisenberg_sector_free(S);
            irrep_sg_little_group_irrep_free(mu);
            irrep_sg_little_group_free(lg);
        }
    free(psi_full);

    qsort(all_states, n_states, sizeof(state_t), cmp_state);

    printf("\n━━━ Combined F=+1 singlet spectrum (lowest 15) ━━━\n");
    printf("  %-3s %-14s %-12s %-7s %-6s\n", "#", "E", "ΔE from E₀", "(kx,ky)", "level");
    printf("  --- -------------- ------------ ------- ------\n");
    int show = n_states > 15 ? 15 : n_states;
    for (int i = 0; i < show; ++i) {
        double dE = all_states[i].E - all_states[0].E;
        printf("  %-3d %+14.8f %+12.2e (%d,%d)   %d\n",
               i, all_states[i].E, dE,
               all_states[i].kx, all_states[i].ky,
               all_states[i].level);
    }

    /* Check if the top-4 singlet states form a tower:
     * Criterion: ΔE(top 4) ≪ ΔE(4th gap) i.e. E_3 - E_0 < (E_4 - E_3).
     * A stronger Z₂ signal: ΔE_01, ΔE_02, ΔE_03 all < 10% of E_0's bandwidth.
     */
    if (n_states >= 5) {
        double E0 = all_states[0].E;
        double dE_tower_top = all_states[3].E - all_states[0].E;
        double gap_after    = all_states[4].E - all_states[3].E;
        double bandwidth_guess = 0.1 * fabs(E0);
        printf("\n━━━ Tower analysis ━━━\n");
        printf("  4-state tower span:     ΔE_top = %+.4f  (%.3f %% of |E₀|)\n",
               dE_tower_top, 100.0 * dE_tower_top / fabs(E0));
        printf("  Gap after top 4:        ΔE_gap = %+.4f\n", gap_after);
        printf("  Tower-to-gap ratio:     ΔE_top / ΔE_gap = %.4f\n",
               dE_tower_top / gap_after);
        if (dE_tower_top < bandwidth_guess && gap_after > dE_tower_top)
            printf("  ✓ Four-state tower CONSISTENT with Z₂ topological degeneracy.\n");
        else if (gap_after > dE_tower_top)
            printf("  ~ Four-state grouping present but spread is larger than 10%% of E₀;\n"
                   "    could be Z₂ with finite-size splitting, or non-topological low-lying states.\n");
        else
            printf("  ✗ No clean four-state tower; the near-degenerate structure does not\n"
                   "    match the topological prediction on this cluster.\n");

        /* Report the 4 candidates' momenta — for Z₂ we expect 4 different k's. */
        printf("\n  Tower momenta:\n");
        int k_unique = 0;
        int seen[8][8] = {{0}};
        for (int i = 0; i < 4; ++i) {
            int mx = all_states[i].kx, my = all_states[i].ky;
            if (!seen[mx][my]) { ++k_unique; seen[mx][my] = 1; }
            printf("    state %d: k = (%d, %d)\n", i, mx, my);
        }
        printf("  Distinct k-points among top 4: %d / 4\n", k_unique);
        if (k_unique == 4)
            printf("  ✓ 4 distinct momenta — topologically-distinct sectors, consistent with Z₂ tower.\n");
        else if (k_unique >= 2)
            printf("  ~ Some degeneracy within k-sector — may indicate level crossings.\n");
        else
            printf("  ✗ All top 4 states at same k — not a topological tower.\n");
    }

    free(all_states);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
