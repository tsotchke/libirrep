/* SPDX-License-Identifier: MIT */
/* Real-space spin-spin correlation ⟨S_0 · S_r⟩ as a function of separation
 * on the N=24 kagome Heisenberg singlet ground state. Two independent fits
 * discriminate the Z₂-gapped vs Dirac-gapless spin-liquid scenarios:
 *
 *   Gapped (Z₂): ⟨S·S⟩(r) ~ A · exp(−r/ξ)      — finite correlation length
 *   Gapless (Dirac): ⟨S·S⟩(r) ~ A · r^(−η)    — algebraic decay, no scale
 *
 * This complements the two diagnostics already on the N=24 GS:
 *   § 1.9–1.11  γ (area-law) ≈ log 2   → topological phase
 *   § 1.14      max S(k)/N = 0.031      → no magnetic order
 * Adding a finite ξ here tips the balance towards the Z₂ picture; finding
 * pure algebraic decay would tip it towards the Dirac picture.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome_correlation_length
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

static double corr_si_sj(const double _Complex *psi, int i, int j) {
    if (i == j) return 0.75;
    long long mi = 1LL << i, mj = 1LL << j;
    double diag = 0.0, offd = 0.0;
    for (long long s = 0; s < D_FULL; ++s) {
        double a2 = creal(psi[s] * conj(psi[s]));
        int zi = (int)((s >> i) & 1), zj = (int)((s >> j) & 1);
        diag += ((zi == zj) ? 0.25 : -0.25) * a2;
        if (zi != zj) {
            long long sp = s ^ mi ^ mj;
            offd += 0.5 * creal(conj(psi[sp]) * psi[s]);
        }
    }
    return diag + offd;
}

/* Minimum-image distance on a 2D torus parameterised by its primitive
 * lattice vectors (a1, a2). Loops the three neighbouring images to pick
 * the shortest representative. */
static double mi_distance(double rx, double ry,
                          const double a1[2], const double a2[2],
                          int Lx, int Ly) {
    double best = 1e300;
    for (int nx = -1; nx <= 1; ++nx)
        for (int ny = -1; ny <= 1; ++ny) {
            double dx = rx + nx * Lx * a1[0] + ny * Ly * a2[0];
            double dy = ry + nx * Lx * a1[1] + ny * Ly * a2[1];
            double r = sqrt(dx*dx + dy*dy);
            if (r < best) best = r;
        }
    return best;
}

int main(void) {
    printf("=== Kagome N=24 real-space correlation length ξ ===\n");
    printf("    ⟨S_0·S_r⟩ — exponential fit ↔ Z₂ gapped; power-law ↔ Dirac gapless\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);

    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
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
                                 sector_dim, 2, 150, seed, eigs, psi_sec);
    printf("  E_0/N = %+.8f\n", eigs[0] / N_SITES);

    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    double _Complex *psi_full = malloc((size_t)D_FULL * sizeof(double _Complex));
    unfold(G, T, order, w, psi_sec, sector_dim, psi_full);

    double nn = 0.0;
    for (long long s = 0; s < D_FULL; ++s) nn += creal(psi_full[s] * conj(psi_full[s]));
    printf("  ||ψ||² = %.10f\n\n", nn);

    /* Site positions + primitive lattice vectors for minimum-image. */
    double site_x[N_SITES], site_y[N_SITES];
    for (int i = 0; i < N_SITES; ++i) {
        double xy[2];
        irrep_lattice_site_position(L, i, xy);
        site_x[i] = xy[0]; site_y[i] = xy[1];
    }
    double a1[2], a2[2];
    irrep_lattice_primitive_vectors(L, a1, a2);
    int Lx = irrep_lattice_Lx(L), Ly = irrep_lattice_Ly(L);

    /* Compute ⟨S_0·S_r⟩ for all j ≥ 1, translate-average over all i to
     * average out boundary-specific fluctuations (with PBC this recovers
     * translation invariance in expectation). */
    printf("  Computing translation-averaged ⟨S_0·S_r⟩ …\n");
    double tC = now_sec();
    int npairs = 0;
    typedef struct { double r; double c; } pt_t;
    pt_t pts[N_SITES * N_SITES];
    for (int i = 0; i < N_SITES; ++i) {
        for (int j = 0; j < N_SITES; ++j) {
            if (i == j) continue;
            double dx = site_x[j] - site_x[i];
            double dy = site_y[j] - site_y[i];
            double r = mi_distance(dx, dy, a1, a2, Lx, Ly);
            pts[npairs].r = r;
            pts[npairs].c = corr_si_sj(psi_full, i, j);
            ++npairs;
        }
    }
    printf("  correlations done in %.2f s over %d pairs\n\n",
           now_sec() - tC, npairs);

    /* Bin by distance (tolerance 1e-3) and report |avg ⟨S·S⟩|. */
    double uniq_r[64]; int uniq_cnt[64]; double uniq_sum[64];
    int nu = 0;
    for (int p = 0; p < npairs; ++p) {
        int found = -1;
        for (int u = 0; u < nu; ++u)
            if (fabs(uniq_r[u] - pts[p].r) < 1e-3) { found = u; break; }
        if (found >= 0) { uniq_sum[found] += pts[p].c; ++uniq_cnt[found]; }
        else {
            uniq_r[nu] = pts[p].r;
            uniq_sum[nu] = pts[p].c;
            uniq_cnt[nu] = 1;
            ++nu;
        }
    }
    /* Sort by r. */
    for (int a = 0; a < nu - 1; ++a)
        for (int b = a + 1; b < nu; ++b)
            if (uniq_r[b] < uniq_r[a]) {
                double tr = uniq_r[a]; uniq_r[a] = uniq_r[b]; uniq_r[b] = tr;
                double ts = uniq_sum[a]; uniq_sum[a] = uniq_sum[b]; uniq_sum[b] = ts;
                int   tc = uniq_cnt[a]; uniq_cnt[a] = uniq_cnt[b]; uniq_cnt[b] = tc;
            }
    printf("  Translation-averaged C(r) = (1/n_r) Σ ⟨S_i·S_j⟩ at |r_i−r_j|=r:\n");
    printf("  %-6s %5s %12s %12s\n", "r", "n_r", "C(r)", "|C(r)|");
    printf("  ------ ----- ------------ ------------\n");
    for (int u = 0; u < nu; ++u) {
        double avg = uniq_sum[u] / uniq_cnt[u];
        printf("  %6.4f %5d %+12.6e %12.6e\n",
               uniq_r[u], uniq_cnt[u], avg, fabs(avg));
    }

    /* Exponential fit: log|C(r)| = log A − r/ξ   on r ≥ r_min.
     * Power-law fit:  log|C(r)| = log A − η·log r   on same r ≥ r_min.
     * Skip r=0 (on-site) and any C(r) below 1e-6 (noise floor). Restrict
     * to r up to half the torus diameter to avoid wrap-around saturation. */
    /* Widen the fit window: the torus is 2×4 kagome, so min-image
     * distances beyond the shortest half-diameter are still real
     * minimum-image separations; they stop being bulk-representative
     * only when they approach the large half-diameter. Use the longer
     * half-diameter as the cap. */
    double a1m = sqrt(a1[0]*a1[0] + a1[1]*a1[1]);
    double a2m = sqrt(a2[0]*a2[0] + a2[1]*a2[1]);
    double rmax = 0.5 * ( (Lx*a1m > Ly*a2m) ? Lx*a1m : Ly*a2m );
    double rmin = 0.99;
    (void)Lx; (void)Ly;

    int    nfit = 0;
    double xr[64], xlr[64], yl[64];
    for (int u = 0; u < nu; ++u) {
        double avg = uniq_sum[u] / uniq_cnt[u];
        double aC = fabs(avg);
        if (aC < 1e-6) continue;
        if (uniq_r[u] < rmin) continue;
        if (uniq_r[u] > rmax) continue;
        xr[nfit]  = uniq_r[u];
        xlr[nfit] = log(uniq_r[u]);
        yl[nfit]  = log(aC);
        ++nfit;
    }
    printf("\n  Fit window: %.3f ≤ r ≤ %.3f  (%d points)\n", rmin, rmax, nfit);

    /* Two linear regressions. */
    double mx=0, my=0, mxx=0, mxy=0;
    for (int i = 0; i < nfit; ++i) { mx+=xr[i]; my+=yl[i]; mxx+=xr[i]*xr[i]; mxy+=xr[i]*yl[i]; }
    mx/=nfit; my/=nfit; mxx/=nfit; mxy/=nfit;
    double slope_exp = (mxy - mx*my) / (mxx - mx*mx);
    double icpt_exp  = my - slope_exp*mx;
    double xi        = -1.0 / slope_exp;
    double res_exp=0, tot_exp=0;
    for (int i = 0; i < nfit; ++i) {
        double p = slope_exp*xr[i] + icpt_exp;
        res_exp += (yl[i]-p)*(yl[i]-p);
        tot_exp += (yl[i]-my)*(yl[i]-my);
    }
    double R2_exp = 1.0 - res_exp/tot_exp;

    mx=my=mxx=mxy=0;
    for (int i = 0; i < nfit; ++i) { mx+=xlr[i]; my+=yl[i]; mxx+=xlr[i]*xlr[i]; mxy+=xlr[i]*yl[i]; }
    mx/=nfit; my/=nfit; mxx/=nfit; mxy/=nfit;
    double slope_pow = (mxy - mx*my) / (mxx - mx*mx);
    double icpt_pow  = my - slope_pow*mx;
    double eta       = -slope_pow;
    double res_pow=0, tot_pow=0;
    for (int i = 0; i < nfit; ++i) {
        double p = slope_pow*xlr[i] + icpt_pow;
        res_pow += (yl[i]-p)*(yl[i]-p);
        tot_pow += (yl[i]-my)*(yl[i]-my);
    }
    double R2_pow = 1.0 - res_pow/tot_pow;

    printf("\n  Exponential fit:  |C(r)| = A·exp(−r/ξ)\n");
    printf("    ξ (correlation length)  = %+.4f (units of NN bond length)\n", xi);
    printf("    log A                   = %+.4f   A = %.4e\n", icpt_exp, exp(icpt_exp));
    printf("    R²                      = %.4f\n", R2_exp);

    printf("\n  Power-law fit:    |C(r)| = A·r^(−η)\n");
    printf("    η (decay exponent)      = %+.4f\n", eta);
    printf("    log A                   = %+.4f   A = %.4e\n", icpt_pow, exp(icpt_pow));
    printf("    R²                      = %.4f\n", R2_pow);

    printf("\n  VERDICT:\n");
    if (R2_exp > R2_pow + 0.05) {
        printf("    ✓ Exponential wins by ΔR² = %+.3f  →  GAPPED Z₂ spin liquid\n",
               R2_exp - R2_pow);
        printf("      ξ = %.3f bond-lengths (finite)\n", xi);
    } else if (R2_pow > R2_exp + 0.05) {
        printf("    ✓ Power-law wins by ΔR² = %+.3f  →  GAPLESS Dirac spin liquid\n",
               R2_pow - R2_exp);
        printf("      η = %.3f (algebraic decay)\n", eta);
    } else {
        printf("    ~ Inconclusive at N=24 (|ΔR²| = %.3f).\n",
               fabs(R2_exp - R2_pow));
        printf("      Finite N lacks enough r-range to discriminate.\n");
    }

    printf("\n  Total wall-clock: %.2f s\n", now_sec() - t0);

    free(psi_full); free(psi_sec); free(seed); free(w);
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
