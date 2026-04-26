/* SPDX-License-Identifier: MIT */
/* Unfiltered tower scan on kagome 2×4 N=24, reporting F=±1 eigenvalues
 * for the lowest 8 states in every k-sector. Looks for the 2 missing
 * topological sectors that may be F = −1 (spin-flip-odd) instead of F=+1.
 *
 * On Z₂ spin liquids, the 4-fold torus GS manifold typically splits as
 * 2 F=+1 + 2 F=−1 under the global spin-flip symmetry Z_2^{spin}. The
 * earlier F-filtered tower scan missed F=−1 states by design.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_tower_F
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

int main(void) {
    printf("=== Kagome 2×4 N=24 unfiltered tower scan (F=±1 labelled) ===\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);

    const int K_PER = 8;
    const int Lx = 2, Ly = 4;
    uint64_t mask = (1ULL << N_SITES) - 1;

    double _Complex *psi_full = malloc((size_t)D_FULL * sizeof(double _Complex));

    printf("  %-5s %-3s %-14s %-8s %-6s\n", "(k)", "lv", "E", "F", "type");
    printf("  --------------------------------------------\n");

    /* Collect all states globally. */
    typedef struct { double E; double F; int mx, my, lv; } entry_t;
    entry_t *all = malloc(sizeof(entry_t) * Lx * Ly * K_PER);
    int ne = 0;

    for (int my = 0; my < Ly; ++my)
        for (int mx = 0; mx < Lx; ++mx) {
            irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, mx, my);
            double _Complex chi1[1] = {1.0 + 0.0 * I};
            irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
            irrep_sg_heisenberg_sector_t *S =
                irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
            long long sdim = irrep_sg_heisenberg_sector_dim(S);
            double _Complex *seed = malloc((size_t)sdim * sizeof(double _Complex));
            for (long long i = 0; i < sdim; ++i)
                seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
            int K = K_PER > sdim ? (int)sdim : K_PER;
            double *eigs = malloc((size_t)K * sizeof(double));
            double _Complex *psi_all =
                malloc((size_t)K * sdim * sizeof(double _Complex));
            int it = 200 > sdim ? (int)sdim : 200;
            irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                         sdim, K, it, seed, eigs, psi_all);
            int order = irrep_space_group_order(G);
            double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
            irrep_sg_projector_weights(lg, mu, w);

            for (int k = 0; k < K; ++k) {
                unfold(G, T, order, w, psi_all + (size_t)k * sdim,
                       sdim, psi_full);
                double _Complex Fk = 0;
                for (long long s = 0; s < D_FULL; ++s)
                    Fk += conj(psi_full[s]) * psi_full[s ^ mask];
                const char *tl = (creal(Fk) > 0.5) ? "F=+1" :
                                  (creal(Fk) < -0.5) ? "F=-1" : "mixed";
                printf("  (%d,%d)  %-3d %+14.8f %+8.4f %s\n",
                       mx, my, k, eigs[k], creal(Fk), tl);
                all[ne].E = eigs[k]; all[ne].F = creal(Fk);
                all[ne].mx = mx; all[ne].my = my; all[ne].lv = k;
                ++ne;
            }
            free(seed); free(eigs); free(psi_all); free(w);
            irrep_sg_heisenberg_sector_free(S);
            irrep_sg_little_group_irrep_free(mu);
            irrep_sg_little_group_free(lg);
            fflush(stdout);
        }

    /* Sort globally. */
    for (int i = 0; i < ne-1; ++i)
        for (int j = i+1; j < ne; ++j)
            if (all[j].E < all[i].E) {
                entry_t tmp = all[i]; all[i] = all[j]; all[j] = tmp;
            }

    printf("\n━━━ Global combined spectrum (lowest 20) ━━━\n");
    printf("  %-3s %-14s %-8s %-8s %-9s\n", "#", "E", "ΔE_0", "F", "(kx,ky)");
    for (int i = 0; i < 20 && i < ne; ++i)
        printf("  %-3d %+14.8f %+8.4f %+8.4f (%d,%d,lv%d)\n",
               i, all[i].E, all[i].E - all[0].E, all[i].F,
               all[i].mx, all[i].my, all[i].lv);

    /* Separate F=+1 and F=−1 spectra. */
    printf("\n━━━ Lowest F=+1 and F=−1 states ━━━\n");
    int cnt_pos = 0, cnt_neg = 0;
    printf("  F=+1 lowest 6:\n");
    for (int i = 0; i < ne && cnt_pos < 6; ++i) {
        if (all[i].F > 0.5) {
            printf("    E=%+14.8f  ΔE_0=%+.4f  (%d,%d,lv%d)\n",
                   all[i].E, all[i].E - all[0].E,
                   all[i].mx, all[i].my, all[i].lv);
            ++cnt_pos;
        }
    }
    printf("  F=-1 lowest 6:\n");
    for (int i = 0; i < ne && cnt_neg < 6; ++i) {
        if (all[i].F < -0.5) {
            printf("    E=%+14.8f  ΔE_0=%+.4f  (%d,%d,lv%d)\n",
                   all[i].E, all[i].E - all[0].E,
                   all[i].mx, all[i].my, all[i].lv);
            ++cnt_neg;
        }
    }

    free(all); free(psi_full);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
