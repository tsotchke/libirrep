/* SPDX-License-Identifier: MIT */
/* γ(area-law) bootstrap on the true N=24 k=(0,π) singlet GS.
 *
 * The audit (§1.17b, D) showed that contiguous, stride-2 and
 * mixed-sublattice region families give γ = 0.88 / 0.62 / 0.04
 * respectively — a 0.84 spread on the SAME ground state.
 *
 * Generate N random region choices of size |A| ∈ [1, 9] for each
 * sample; fit S_A = α|∂A| + β; collect the γ = −β distribution.
 * If γ ≈ log 2 is physical, the mode is near 0.69. If it's an
 * artifact of the contiguous region family, the mode is elsewhere
 * and log 2 is a tail event.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_gamma_bootstrap
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

#define N_SAMPLES 20   /* random-region samples; each fits 9 |A| values */

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

static int boundary_bonds(const int *A, int nA, const int *bi, const int *bj, int nb) {
    char in_A[64] = {0};
    for (int k = 0; k < nA; ++k) in_A[A[k]] = 1;
    int boundary = 0;
    for (int b = 0; b < nb; ++b)
        if (in_A[bi[b]] != in_A[bj[b]]) ++boundary;
    return boundary;
}

/* Random permutation of [0, N_SITES) via Fisher-Yates.  Caller reads
 * the first |A| entries for a random region. */
static void random_permutation(int *perm, int n, unsigned *state) {
    for (int i = 0; i < n; ++i) perm[i] = i;
    for (int i = n - 1; i > 0; --i) {
        *state = (*state) * 1103515245u + 12345u;
        int j = (int)(((*state) >> 8) % (unsigned)(i + 1));
        int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
    }
}

static int cmp_dbl(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

int main(void) {
    printf("=== γ bootstrap on N=24 kagome TRUE singlet GS, k=(%d,%d) ===\n",
           KX_GS, KY_GS);
    printf("    %d random region samples × 9 region sizes each\n\n", N_SAMPLES);
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

    /* Bootstrap loop. */
    double *gammas = malloc(sizeof(double) * N_SAMPLES);
    double *R2s    = malloc(sizeof(double) * N_SAMPLES);
    double *alphas = malloc(sizeof(double) * N_SAMPLES);
    unsigned rng_state = 12345u;

    /* For each sample, pick 9 nested regions: pick a permutation, then
     * use its first nA sites for nA=1..9. The 9 points on the S_A vs
     * |∂A| plane come from that single random ordering. This is
     * faithful to the "nested region family" methodology but randomises
     * WHICH family. */
    const int nA_max = 9;
    for (int r = 0; r < N_SAMPLES; ++r) {
        int perm[32];
        random_permutation(perm, N_SITES, &rng_state);

        double xs[16], ys[16];
        double ts = now_sec();
        for (int nA = 1; nA <= nA_max; ++nA) {
            int A[16];
            for (int i = 0; i < nA; ++i) A[i] = perm[i];
            int bdry = boundary_bonds(A, nA, bi, bj, nb);
            long long dA = 1LL << nA;
            double _Complex *rho = malloc((size_t)dA * dA * sizeof(double _Complex));
            irrep_partial_trace(N_SITES, 2, psi, A, nA, rho);
            double S_A = irrep_entropy_vonneumann(rho, (int)dA);
            xs[nA - 1] = bdry; ys[nA - 1] = S_A;
            free(rho);
        }
        int npts = nA_max;
        double mxv=0, myv=0, mxxv=0, mxyv=0;
        for (int i = 0; i < npts; ++i) {
            mxv+=xs[i]; myv+=ys[i]; mxxv+=xs[i]*xs[i]; mxyv+=xs[i]*ys[i];
        }
        mxv/=npts; myv/=npts; mxxv/=npts; mxyv/=npts;
        double alpha = (mxyv - mxv*myv)/(mxxv - mxv*mxv);
        double beta  = myv - alpha*mxv;
        double rs=0, ts_tot=0;
        for (int i = 0; i < npts; ++i) {
            double pred = alpha*xs[i] + beta;
            rs += (ys[i]-pred)*(ys[i]-pred);
            ts_tot += (ys[i]-myv)*(ys[i]-myv);
        }
        double R2 = 1.0 - rs/ts_tot;
        gammas[r] = -beta;
        R2s[r]    = R2;
        alphas[r] = alpha;
        printf("  sample %2d: α=%+.4f γ=%+.4f R²=%.4f (%.1fs)\n",
               r, alpha, -beta, R2, now_sec() - ts);
        fflush(stdout);
    }

    /* Distribution stats. */
    double g_mean = 0, g_var = 0;
    for (int r = 0; r < N_SAMPLES; ++r) g_mean += gammas[r];
    g_mean /= N_SAMPLES;
    for (int r = 0; r < N_SAMPLES; ++r) g_var += (gammas[r]-g_mean)*(gammas[r]-g_mean);
    g_var /= (N_SAMPLES - 1);
    double g_std = sqrt(g_var);

    double *g_sorted = malloc(sizeof(double) * N_SAMPLES);
    memcpy(g_sorted, gammas, sizeof(double) * N_SAMPLES);
    qsort(g_sorted, N_SAMPLES, sizeof(double), cmp_dbl);
    double g_median = g_sorted[N_SAMPLES / 2];
    double g_q25 = g_sorted[N_SAMPLES / 4];
    double g_q75 = g_sorted[(3 * N_SAMPLES) / 4];
    double g_min = g_sorted[0], g_max = g_sorted[N_SAMPLES - 1];

    printf("\n━━━ γ distribution over %d random regions ━━━\n", N_SAMPLES);
    printf("  mean     = %+.4f\n", g_mean);
    printf("  std      = %.4f\n", g_std);
    printf("  median   = %+.4f\n", g_median);
    printf("  IQR      = [%+.4f, %+.4f]\n", g_q25, g_q75);
    printf("  full     = [%+.4f, %+.4f]\n", g_min, g_max);
    printf("  log 2    = %+.4f   (for comparison)\n", log(2.0));

    /* Histogram, 20 bins across [g_min, g_max]. */
    int nbins = 20;
    int hist[32] = {0};
    double lo = g_min - 1e-6, hi = g_max + 1e-6;
    double bw = (hi - lo) / nbins;
    for (int r = 0; r < N_SAMPLES; ++r) {
        int b = (int)((gammas[r] - lo) / bw);
        if (b < 0) b = 0; if (b >= nbins) b = nbins - 1;
        hist[b]++;
    }
    printf("\n  histogram:\n");
    for (int b = 0; b < nbins; ++b) {
        double lo_b = lo + b * bw, hi_b = lo_b + bw;
        printf("  [%+6.3f, %+6.3f]  %3d  ", lo_b, hi_b, hist[b]);
        for (int k = 0; k < hist[b]; ++k) printf("█");
        printf("\n");
    }

    /* Assessment. */
    double contiguous_gamma = +0.878;   /* from §1.17c */
    double log2 = log(2.0);
    printf("\n  SPECIAL POINTS:\n");
    /* how many samples fall within ±0.1 of log 2? */
    int near_log2 = 0, near_zero = 0;
    for (int r = 0; r < N_SAMPLES; ++r) {
        if (fabs(gammas[r] - log2) < 0.1) ++near_log2;
        if (fabs(gammas[r]) < 0.1) ++near_zero;
    }
    printf("  samples within ±0.1 of log 2 : %d / %d  (%.0f%%)\n",
           near_log2, N_SAMPLES, 100.0*near_log2/N_SAMPLES);
    printf("  samples within ±0.1 of 0     : %d / %d  (%.0f%%)\n",
           near_zero, N_SAMPLES, 100.0*near_zero/N_SAMPLES);
    printf("  contiguous-family γ (§1.17c): %+.4f  — %s\n",
           contiguous_gamma,
           (fabs(contiguous_gamma - g_mean) < g_std) ? "within 1σ of mean" :
           (fabs(contiguous_gamma - g_mean) < 2*g_std) ? "within 2σ of mean" :
           "tail-region of distribution");

    printf("\n  CONCLUSION:\n");
    if (fabs(g_median - log2) < 0.1 && g_std < 0.15)
        printf("    ✓ Random-region γ distribution concentrates near log 2\n"
               "      with tight width — γ identification survives.\n");
    else if (g_std > 0.3)
        printf("    ✗ γ distribution is too wide (σ=%.2f); finite-N γ\n"
               "      extraction is region-dominated on this cluster.\n", g_std);
    else
        printf("    ~ γ distribution (median=%+.3f, σ=%.2f) does not clearly\n"
               "      align with any single phase prediction.\n", g_median, g_std);

    free(gammas); free(R2s); free(alphas); free(g_sorted); free(psi);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
