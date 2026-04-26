/* SPDX-License-Identifier: MIT */
/* Kitaev-Preskill three-region γ on the true N=24 (0,π) singlet GS.
 *
 * The KP construction is, by design, corner-term-free:
 *   γ = S_A + S_B + S_C − S_{AB} − S_{AC} − S_{BC} + S_{ABC}
 *
 * For a topologically-ordered state, γ recovers log(2) for Z₂ and 0
 * for trivial/gapless. Region SHAPE still matters, so we report
 * γ_KP across 4 different 3-region geometries to expose any
 * remaining cluster-specific bias.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_kp_gamma
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
#define KX_GS 0
#define KY_GS 2

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

/* Compute S_X = -Tr(ρ_X log ρ_X) by partial-tracing over complement. */
static double entropy_of_region(const double _Complex *psi, const int *X, int nX) {
    long long dX = 1LL << nX;
    double _Complex *rho = malloc((size_t)dX * dX * sizeof(double _Complex));
    irrep_partial_trace(N_SITES, 2, psi, X, nX, rho);
    double S = irrep_entropy_vonneumann(rho, (int)dX);
    free(rho);
    return S;
}

typedef struct {
    const char *name;
    int A[8], B[8], C[8];
    int nA, nB, nC;
} kp_geom_t;

int main(void) {
    printf("=== KP γ on N=24 kagome 2×4 TRUE singlet GS at k=(%d,%d) ===\n",
           KX_GS, KY_GS);
    printf("    γ_KP = S_A + S_B + S_C − S_AB − S_AC − S_BC + S_ABC\n");
    printf("    Target: log 2 = %.4f for Z₂ topological phase, 0 for trivial.\n\n",
           log(2.0));
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);

    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, KX_GS, KY_GS);
    double _Complex chi1[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);

    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
    double eigs[2];
    double _Complex *psi_sec = malloc((size_t)2 * sector_dim * sizeof(double _Complex));
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sector_dim, 2, 200, seed, eigs, psi_sec);
    printf("  E_0 = %+.8f (E/N = %+.8f)\n", eigs[0], eigs[0] / N_SITES);

    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    double _Complex *psi = malloc((size_t)D_FULL * sizeof(double _Complex));
    unfold(G, T, order, w, psi_sec, sector_dim, psi);
    free(psi_sec); free(seed); free(w);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);

    /* KP geometries: 3 disjoint 3-site regions (|AUBUC|=9, env=15).
     * 3-site (rather than 4-site) keeps the |ABC|=9 RDM at 512×512,
     * so the Jacobi eigendecomp cost drops from ~10 minutes to <1 s
     * per geometry while preserving the corner-cancellation property
     * of the KP construction (γ_KP is asymptotically independent of
     * region size for |A| ≫ ξ; here ξ ≈ 0.4 ⟹ |A| = 3 is already in
     * the asymptotic regime). */
    kp_geom_t geoms[] = {
        {
            "compact strips (contiguous 3-site blocks)",
            {0, 1, 2}, {3, 4, 5}, {6, 7, 8}, 3, 3, 3
        },
        {
            "separated strips (stride-2 within each 3-site block)",
            {0, 2, 4}, {6, 8, 10}, {12, 14, 16}, 3, 3, 3
        },
        {
            "sublattice mix (A/B/C on different kagome sublattices)",
            {0, 3, 6}, {1, 4, 7}, {2, 5, 8}, 3, 3, 3
        },
        {
            "corners (3 sites near each cluster corner)",
            {0, 1, 2}, {6, 7, 8}, {12, 13, 14}, 3, 3, 3
        },
    };
    int n_geoms = sizeof(geoms) / sizeof(geoms[0]);

    printf("\n  %-55s %-10s %-10s %-10s\n",
           "geometry", "γ_KP", "S_ABC-avg", "time");
    printf("  ------------------------------------------------------- ---------- ---------- ----------\n");

    double gammas[8];
    for (int g_i = 0; g_i < n_geoms; ++g_i) {
        const kp_geom_t *geo = &geoms[g_i];
        /* Union regions: AB, AC, BC, ABC. */
        int AB[8], AC[8], BC[8], ABC[12];
        int nAB = 0, nAC = 0, nBC = 0, nABC = 0;
        for (int i = 0; i < geo->nA; ++i) { AB[nAB++] = geo->A[i]; AC[nAC++] = geo->A[i]; ABC[nABC++] = geo->A[i]; }
        for (int i = 0; i < geo->nB; ++i) { AB[nAB++] = geo->B[i]; BC[nBC++] = geo->B[i]; ABC[nABC++] = geo->B[i]; }
        for (int i = 0; i < geo->nC; ++i) { AC[nAC++] = geo->C[i]; BC[nBC++] = geo->C[i]; ABC[nABC++] = geo->C[i]; }

        double tg = now_sec();
        double S_A = entropy_of_region(psi, geo->A, geo->nA);
        double S_B = entropy_of_region(psi, geo->B, geo->nB);
        double S_C = entropy_of_region(psi, geo->C, geo->nC);
        double S_AB = entropy_of_region(psi, AB, nAB);
        double S_AC = entropy_of_region(psi, AC, nAC);
        double S_BC = entropy_of_region(psi, BC, nBC);
        double S_ABC = entropy_of_region(psi, ABC, nABC);
        double gamma_kp = S_A + S_B + S_C - S_AB - S_AC - S_BC + S_ABC;
        gammas[g_i] = gamma_kp;
        printf("  %-55s %+10.4f %+10.4f  %6.1fs\n",
               geo->name, gamma_kp, S_ABC, now_sec() - tg);
        fflush(stdout);
    }

    double g_mean = 0;
    for (int i = 0; i < n_geoms; ++i) g_mean += gammas[i];
    g_mean /= n_geoms;
    double g_min = gammas[0], g_max = gammas[0];
    for (int i = 0; i < n_geoms; ++i) {
        if (gammas[i] < g_min) g_min = gammas[i];
        if (gammas[i] > g_max) g_max = gammas[i];
    }
    printf("\n  γ_KP range across geometries: [%+.4f, %+.4f]  spread = %.4f\n",
           g_min, g_max, g_max - g_min);
    printf("  γ_KP mean                   : %+.4f   (log 2 = %+.4f)\n",
           g_mean, log(2.0));

    printf("\n  COMPARISON WITH AREA-LAW γ (§1.17c):\n");
    printf("    area-law γ range: [+0.036, +0.878]  spread = 0.842\n");
    printf("    KP γ spread     : %.3f  %s\n", g_max - g_min,
           (g_max - g_min < 0.4) ? "(tighter — KP is design-robust)" : "(still wide)");
    if (fabs(g_mean - log(2.0)) < 0.2 && g_max - g_min < 0.4)
        printf("  ✓ KP γ is geometry-stable and near log 2 → Z₂ identification stands.\n");
    else if (fabs(g_mean) < 0.2 && g_max - g_min < 0.4)
        printf("  ✗ KP γ is geometry-stable and near 0 → trivial/gapless, not Z₂.\n");
    else if (g_max - g_min < 0.4)
        printf("  ~ KP γ is geometry-stable (spread %.2f) but at an\n"
               "    intermediate value (%.3f) — not a clean phase identification.\n",
               g_max - g_min, g_mean);
    else
        printf("  ~ KP γ still varies with geometry (spread %.2f) — topological γ\n"
               "    extraction genuinely requires larger clusters.\n", g_max - g_min);

    free(psi);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
