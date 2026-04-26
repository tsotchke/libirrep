/* SPDX-License-Identifier: MIT */
/* Methodology null-control: area-law γ extraction on the square-lattice
 * Heisenberg antiferromagnet. The published thermodynamic-limit state
 * is Néel-ordered (spontaneously broken SU(2)); topological γ = 0.
 *
 * If our area-law fit gives γ_square ≈ 0 while γ_kagome ≈ log 2 at
 * matched cluster size, the methodology cleanly distinguishes gapped
 * topological phases from spontaneously-ordered phases — a smoking-gun
 * validation of the library's γ extraction pipeline.
 *
 * Cluster: square 4×4 (N=16), p4mm wallpaper group. Γ-A_1 sector
 * (singlet GS known to sit here).
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/square_arealaw_null
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

#define N_SITES 16
#define D_FULL  (1LL << N_SITES)
#define POPCOUNT 8

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

int main(void) {
    printf("=== Square-lattice Heisenberg N=16 — γ NULL CONTROL ===\n");
    printf("    Published prediction: Néel-ordered, γ = 0\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P4MM);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);

    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);
    printf("  rep table: %lld reps (p4mm order %d)\n",
           irrep_sg_rep_table_count(T), irrep_space_group_order(G));

    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
    irrep_sg_little_group_irrep_t *A1 =
        irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_A1);
    irrep_sg_heisenberg_sector_t *S =
        irrep_sg_heisenberg_sector_build_at_k(H, T, lg, A1);
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);
    printf("  sector (Γ, A_1): dim = %lld\n", sector_dim);

    /* Lanczos with eigvec. */
    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.1 * sin(0.41 * i) + I * 0.05 * cos(0.29 * i);
    double          *eigs = malloc(sizeof(double) * 4);
    double _Complex *psi_all = malloc((size_t)4 * sector_dim * sizeof(double _Complex));
    int max_it = 150 > sector_dim ? (int)sector_dim : 150;
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sector_dim, 4, max_it, seed, eigs, psi_all);
    printf("  GS: E = %+.6f, E/N = %+.6f (published ≈ -0.7015)\n",
           eigs[0], eigs[0] / N_SITES);
    free(seed);

    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, A1, w);

    double _Complex *psi_full = malloc((size_t)D_FULL * sizeof(double _Complex));
    unfold(G, T, order, w, psi_all, sector_dim, psi_full);

    /* F-diagnostic. */
    double _Complex Fexp = 0.0;
    uint64_t mask = (1ULL << N_SITES) - 1;
    for (long long s = 0; s < D_FULL; ++s)
        Fexp += conj(psi_full[s]) * psi_full[s ^ mask];
    printf("  ⟨ψ|F|ψ⟩ = %+.8f (should be +1 for singlet Néel GS)\n", creal(Fexp));

    /* Area-law fit on nested contiguous regions. */
    printf("\n  Area-law fit:\n");
    printf("  %3s %5s %12s\n", "|A|", "|∂A|", "S_A");
    int sizes[] = {1, 2, 3, 4, 5, 6, 7, 8};
    int npts = sizeof(sizes) / sizeof(sizes[0]);
    double xs[16], ys[16];
    for (int i = 0; i < npts; ++i) {
        int nA = sizes[i];
        int A[16];
        for (int j = 0; j < nA; ++j) A[j] = j;
        int bdry = boundary_bonds(A, nA, bi, bj, nb);
        long long dA = 1LL << nA;
        double _Complex *rho = malloc((size_t)dA * dA * sizeof(double _Complex));
        irrep_partial_trace(N_SITES, 2, psi_full, A, nA, rho);
        double S_A = irrep_entropy_vonneumann(rho, (int)dA);
        printf("  %3d %5d %+12.8f\n", nA, bdry, S_A);
        xs[i] = bdry; ys[i] = S_A;
        free(rho);
    }

    double mx=0, my=0, mxx=0, mxy=0;
    for (int i = 0; i < npts; ++i) { mx+=xs[i]; my+=ys[i]; mxx+=xs[i]*xs[i]; mxy+=xs[i]*ys[i]; }
    mx /= npts; my /= npts; mxx /= npts; mxy /= npts;
    double alpha = (mxy - mx*my) / (mxx - mx*mx);
    double beta  = my - alpha * mx;
    double ss_res=0, ss_tot=0;
    for (int i = 0; i < npts; ++i) {
        double pred = alpha * xs[i] + beta;
        ss_res += (ys[i]-pred)*(ys[i]-pred);
        ss_tot += (ys[i]-my)*(ys[i]-my);
    }
    double R2 = 1.0 - ss_res / ss_tot;

    printf("\n  α = %+.6f   β = %+.6f\n", alpha, beta);
    printf("  γ_extrapolated = %+.6f   R² = %.4f\n", -beta, R2);
    printf("\n  NULL-CONTROL VERDICT:\n");
    if (fabs(-beta) < 0.2)
        printf("  ✓ γ ≈ 0 → methodology correctly identifies Néel (non-topological) phase\n");
    else if (-beta > 0.5)
        printf("  ✗ γ ≫ 0 → methodology has systematic bias overestimating γ at finite N\n");
    else
        printf("  ~ γ intermediate — finite-size correction present but not smoking-gun\n");

    printf("\n  Total wall-clock: %.2f s\n", now_sec() - t0);

    free(psi_full); free(psi_all); free(eigs); free(w);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(A1);
    irrep_sg_little_group_free(lg);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
