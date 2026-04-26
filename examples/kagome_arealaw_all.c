/* SPDX-License-Identifier: MIT */
/* Area-law γ extrapolation across all three accessible kagome singlet
 * clusters (N=12, 18, 24). Unified analysis:
 *
 *   S_A = α(N) · |∂A| + β(N)
 *   γ(N) = −β(N)
 *
 * For topologically-ordered phases, γ(N) → γ_topological as N → ∞.
 * Reports γ at each N + variance of the linear fit (R²).
 *
 * Produces a publishable three-point γ scaling table.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome_arealaw_all
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

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void unfold_to_dense(const irrep_space_group_t *G,
                            const irrep_sg_rep_table_t *T,
                            int order, const double _Complex *w,
                            const double _Complex *psi_sector, long long sector_dim,
                            int num_sites,
                            double _Complex *psi_full) {
    long long D = 1LL << num_sites;
    memset(psi_full, 0, (size_t)D * sizeof(double _Complex));
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

static int boundary_bonds(const int *A, int nA, const int *bi, const int *bj, int nb) {
    char in_A[64] = {0};
    for (int k = 0; k < nA; ++k) in_A[A[k]] = 1;
    int boundary = 0;
    for (int b = 0; b < nb; ++b)
        if (in_A[bi[b]] != in_A[bj[b]]) ++boundary;
    return boundary;
}

/* Find the lowest F=+1 eigenvalue-eigenvector among k_scan Lanczos states. */
static int find_singlet(const irrep_space_group_t *G, const irrep_sg_rep_table_t *T,
                        const irrep_sg_little_group_t *lg,
                        const irrep_sg_little_group_irrep_t *mu,
                        const double _Complex *psi_all, long long sector_dim,
                        int k_scan, int num_sites) {
    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    long long D = 1LL << num_sites;
    double _Complex *psi_full = malloc((size_t)D * sizeof(double _Complex));
    uint64_t mask = (1ULL << num_sites) - 1;
    int result = -1;
    for (int k = 0; k < k_scan; ++k) {
        unfold_to_dense(G, T, order, w,
                        psi_all + (size_t)k * sector_dim, sector_dim, num_sites, psi_full);
        double _Complex Fk = 0.0;
        for (long long s = 0; s < D; ++s) Fk += conj(psi_full[s]) * psi_full[s ^ mask];
        if (creal(Fk) > 0.5) { result = k; break; }
    }
    free(psi_full); free(w);
    return result;
}

static void run_cluster(int Lx, int Ly, irrep_wallpaper_t wp, const char *label) {
    printf("\n=== %s (kagome %d×%d, %s) ===\n", label, Lx, Ly,
           wp == IRREP_WALLPAPER_P6MM ? "p6mm" : "p1");

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, Lx, Ly);
    irrep_space_group_t *G = irrep_space_group_build(L, wp);
    int num_sites = irrep_space_group_num_sites(G);
    int popcount = num_sites / 2;
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(num_sites, nb, bi, bj, 1.0);

    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, popcount);
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
    /* For p6mm use the A_1 named irrep; for p1, a custom trivial. */
    irrep_sg_little_group_irrep_t *mu;
    if (wp == IRREP_WALLPAPER_P6MM)
        mu = irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_A1);
    else {
        double _Complex chi1[1] = {1.0 + 0.0 * I};
        mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
    }
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);

    /* 10 eigvecs, then F-scan for singlet. */
    int K = 10;
    double          *eigs = malloc(sizeof(double) * K);
    double _Complex *psi_all = malloc((size_t)K * sector_dim * sizeof(double _Complex));
    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
    int max_it = 200 > sector_dim ? (int)sector_dim : 200;
    int n_eig = K > sector_dim ? (int)sector_dim : K;
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sector_dim, n_eig, max_it, seed, eigs, psi_all);
    free(seed);

    int singlet_k = find_singlet(G, T, lg, mu, psi_all, sector_dim, n_eig, num_sites);
    if (singlet_k < 0) {
        printf("  no F=+1 state found; skipping\n");
        goto cleanup;
    }
    printf("  singlet at k=%d, E=%+.6f (E/N=%+.6f)\n",
           singlet_k, eigs[singlet_k], eigs[singlet_k] / num_sites);

    /* Unfold singlet GS. */
    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    long long D = 1LL << num_sites;
    double _Complex *psi_full = malloc((size_t)D * sizeof(double _Complex));
    unfold_to_dense(G, T, order, w,
                    psi_all + (size_t)singlet_k * sector_dim, sector_dim,
                    num_sites, psi_full);

    /* Area-law fit: S_A for |A|=1..min(N/2, 9). */
    int max_A = (num_sites / 2 < 9) ? num_sites / 2 : 9;
    double xs[10], ys[10];
    int valid = 0;
    printf("  %3s %5s %12s\n", "|A|", "|∂A|", "S_A");
    for (int nA = 1; nA <= max_A; ++nA) {
        int A[24];
        for (int j = 0; j < nA; ++j) A[j] = j;
        int bdry = boundary_bonds(A, nA, bi, bj, nb);
        long long dA = 1LL << nA;
        double _Complex *rho = malloc((size_t)dA * dA * sizeof(double _Complex));
        irrep_partial_trace(num_sites, 2, psi_full, A, nA, rho);
        double S_A = irrep_entropy_vonneumann(rho, (int)dA);
        printf("  %3d %5d %+12.8f\n", nA, bdry, S_A);
        xs[valid] = bdry; ys[valid] = S_A; ++valid;
        free(rho);
    }

    double mx=0, my=0, mxx=0, mxy=0;
    for (int i = 0; i < valid; ++i) { mx+=xs[i]; my+=ys[i]; mxx+=xs[i]*xs[i]; mxy+=xs[i]*ys[i]; }
    mx /= valid; my /= valid; mxx /= valid; mxy /= valid;
    double alpha = (mxy - mx*my) / (mxx - mx*mx);
    double beta  = my - alpha * mx;
    double ss_res=0, ss_tot=0;
    for (int i = 0; i < valid; ++i) {
        double pred = alpha * xs[i] + beta;
        ss_res += (ys[i]-pred)*(ys[i]-pred);
        ss_tot += (ys[i]-my)*(ys[i]-my);
    }
    double R2 = 1.0 - ss_res / ss_tot;
    printf("  → α = %+.4f   γ_ext = %+.4f   R² = %.4f\n", alpha, -beta, R2);

    free(psi_full); free(w);
cleanup:
    free(eigs); free(psi_all);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
}

int main(void) {
    printf("=== Kagome-Heisenberg γ extraction via area-law fit ===\n");
    printf("    S_A = α·|∂A| − γ + O(1)\n");
    printf("    γ predictions: log 2 = %+.6f (gapped Z_2), 0 (gapless)\n",
           log(2.0));

    double t0 = now_sec();
    run_cluster(2, 2, IRREP_WALLPAPER_P6MM, "N=12");
    run_cluster(2, 3, IRREP_WALLPAPER_P1,   "N=18");
    run_cluster(2, 4, IRREP_WALLPAPER_P1,   "N=24");
    printf("\n  Total wall-clock: %.2f s\n", now_sec() - t0);
    return 0;
}
