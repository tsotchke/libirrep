/* SPDX-License-Identifier: MIT */
/* Finite-size scaling of the real-space correlation length ξ across three
 * kagome-Heisenberg singlet clusters (N = 12, 18, 24). Key observable
 * for discriminating gapped (ξ stays short as N → ∞) from gapless
 * (ξ diverges with N) phases.
 *
 * Selection of singlet GS uses the spin-flip symmetry ⟨ψ|F|ψ⟩ = +1
 * criterion — necessary at N=18 p1 where the lowest Lanczos eigenvalue
 * is a triplet, and applied uniformly for cross-validation at all N.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome_xi_scaling
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

static double corr_si_sj(int num_sites, const double _Complex *psi, int i, int j) {
    if (i == j) return 0.75;
    long long D = 1LL << num_sites;
    long long mi = 1LL << i, mj = 1LL << j;
    double diag = 0.0, offd = 0.0;
    for (long long s = 0; s < D; ++s) {
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

typedef struct {
    double E_singlet;
    double E_per_site;
    double F_eigenvalue;
    int    singlet_k;
    double xi;
    double R2_exp;
    double eta;
    double R2_pow;
    double rmax_used;
    int    nfit;
} xi_result_t;

static xi_result_t compute_xi(int Lx_cells, int Ly_cells, int popcount, int num_sites,
                              int K_eigs, int max_iters) {
    xi_result_t R = {0};
    R.singlet_k = -1;
    long long D = 1LL << num_sites;

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, Lx_cells, Ly_cells);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(num_sites, nb, bi, bj, 1.0);

    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, popcount);
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
    double _Complex chi1[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);

    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
    double *eigs_all = malloc((size_t)K_eigs * sizeof(double));
    double _Complex *psi_all = malloc((size_t)K_eigs * sector_dim * sizeof(double _Complex));
    int max_it = max_iters > sector_dim ? (int)sector_dim : max_iters;
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sector_dim, K_eigs, max_it, seed, eigs_all, psi_all);

    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    double _Complex *psi_full = malloc((size_t)D * sizeof(double _Complex));

    /* Scan the K lowest to find the first F = +1 (singlet) state. */
    uint64_t mask = (num_sites < 64) ? ((1ULL << num_sites) - 1) : ~0ULL;
    for (int kk = 0; kk < K_eigs; ++kk) {
        unfold_to_dense(G, T, order, w, psi_all + (size_t)kk*sector_dim,
                        sector_dim, num_sites, psi_full);
        double _Complex Fk = 0.0;
        for (long long s = 0; s < D; ++s)
            Fk += conj(psi_full[s]) * psi_full[s ^ mask];
        if (creal(Fk) > 0.5) {
            R.singlet_k = kk;
            R.F_eigenvalue = creal(Fk);
            R.E_singlet = eigs_all[kk];
            R.E_per_site = eigs_all[kk] / num_sites;
            break; /* psi_full now holds the singlet GS */
        }
    }
    if (R.singlet_k < 0) {
        fprintf(stderr, "  WARN: no F=+1 state found in %d lowest — report first state\n", K_eigs);
        R.singlet_k = 0;
        R.E_singlet = eigs_all[0];
        R.E_per_site = eigs_all[0] / num_sites;
        unfold_to_dense(G, T, order, w, psi_all, sector_dim, num_sites, psi_full);
    }

    /* Geometry. */
    double a1[2], a2[2];
    irrep_lattice_primitive_vectors(L, a1, a2);
    int Lx = irrep_lattice_Lx(L), Ly = irrep_lattice_Ly(L);

    double *sx = malloc(sizeof(double) * num_sites);
    double *sy = malloc(sizeof(double) * num_sites);
    for (int i = 0; i < num_sites; ++i) {
        double xy[2];
        irrep_lattice_site_position(L, i, xy);
        sx[i] = xy[0]; sy[i] = xy[1];
    }

    /* Translation-averaged C(r). */
    int npairs = 0;
    int cap_pairs = num_sites * num_sites;
    double *rr = malloc(sizeof(double) * cap_pairs);
    double *cc = malloc(sizeof(double) * cap_pairs);
    for (int i = 0; i < num_sites; ++i)
        for (int j = 0; j < num_sites; ++j) {
            if (i == j) continue;
            double dx = sx[j] - sx[i], dy = sy[j] - sy[i];
            double r = mi_distance(dx, dy, a1, a2, Lx, Ly);
            rr[npairs] = r;
            cc[npairs] = corr_si_sj(num_sites, psi_full, i, j);
            ++npairs;
        }

    /* Bin by r. */
    enum { NMAX = 128 };
    double uniq_r[NMAX]; double uniq_sum[NMAX]; int uniq_cnt[NMAX];
    int nu = 0;
    for (int p = 0; p < npairs; ++p) {
        int f = -1;
        for (int u = 0; u < nu; ++u)
            if (fabs(uniq_r[u] - rr[p]) < 1e-3) { f = u; break; }
        if (f >= 0) { uniq_sum[f] += cc[p]; ++uniq_cnt[f]; }
        else if (nu < NMAX) {
            uniq_r[nu] = rr[p]; uniq_sum[nu] = cc[p]; uniq_cnt[nu] = 1; ++nu;
        }
    }
    for (int a = 0; a < nu - 1; ++a)
        for (int b = a + 1; b < nu; ++b)
            if (uniq_r[b] < uniq_r[a]) {
                double tr = uniq_r[a]; uniq_r[a] = uniq_r[b]; uniq_r[b] = tr;
                double ts = uniq_sum[a]; uniq_sum[a] = uniq_sum[b]; uniq_sum[b] = ts;
                int   tc = uniq_cnt[a]; uniq_cnt[a] = uniq_cnt[b]; uniq_cnt[b] = tc;
            }

    double a1m = sqrt(a1[0]*a1[0] + a1[1]*a1[1]);
    double a2m = sqrt(a2[0]*a2[0] + a2[1]*a2[1]);
    R.rmax_used = 0.5 * ( (Lx*a1m > Ly*a2m) ? Lx*a1m : Ly*a2m );
    double rmin = 0.99;

    int nf = 0;
    double xr[NMAX], xlr[NMAX], yl[NMAX];
    for (int u = 0; u < nu; ++u) {
        double avg = uniq_sum[u] / uniq_cnt[u];
        double aC = fabs(avg);
        if (aC < 1e-6) continue;
        if (uniq_r[u] < rmin) continue;
        if (uniq_r[u] > R.rmax_used) continue;
        xr[nf] = uniq_r[u];
        xlr[nf] = log(uniq_r[u]);
        yl[nf] = log(aC);
        ++nf;
    }
    R.nfit = nf;

    if (nf >= 2) {
        double mx=0, my=0, mxx=0, mxy=0;
        for (int i = 0; i < nf; ++i) { mx+=xr[i]; my+=yl[i]; mxx+=xr[i]*xr[i]; mxy+=xr[i]*yl[i]; }
        mx/=nf; my/=nf; mxx/=nf; mxy/=nf;
        double slope = (mxy - mx*my) / (mxx - mx*mx);
        double icpt  = my - slope*mx;
        R.xi = -1.0 / slope;
        double rs=0, ts=0;
        for (int i = 0; i < nf; ++i) {
            double p = slope*xr[i] + icpt;
            rs += (yl[i]-p)*(yl[i]-p);
            ts += (yl[i]-my)*(yl[i]-my);
        }
        R.R2_exp = 1.0 - rs/ts;

        mx=my=mxx=mxy=0;
        for (int i = 0; i < nf; ++i) { mx+=xlr[i]; my+=yl[i]; mxx+=xlr[i]*xlr[i]; mxy+=xlr[i]*yl[i]; }
        mx/=nf; my/=nf; mxx/=nf; mxy/=nf;
        double slope2 = (mxy - mx*my) / (mxx - mx*mx);
        double icpt2  = my - slope2*mx;
        R.eta = -slope2;
        rs=ts=0;
        for (int i = 0; i < nf; ++i) {
            double p = slope2*xlr[i] + icpt2;
            rs += (yl[i]-p)*(yl[i]-p);
            ts += (yl[i]-my)*(yl[i]-my);
        }
        R.R2_pow = 1.0 - rs/ts;
    }

    free(sx); free(sy); free(rr); free(cc);
    free(psi_full); free(psi_all); free(eigs_all); free(seed); free(w);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return R;
}

int main(void) {
    printf("=== Kagome correlation-length finite-size scaling (N=12, 18, 24) ===\n");
    printf("    ξ(N) across three singlet-GS clusters: does ξ saturate?\n\n");
    double t0 = now_sec();

    printf("  Computing N=12 (2×2 kagome)…\n");
    double t12 = now_sec();
    xi_result_t R12 = compute_xi(2, 2, 6, 12, 4, 150);
    printf("    singlet at k=%d (F=%+.4f)  E_0/N=%+.6f  ξ=%.4f  η=%.4f  R²_exp=%.3f  R²_pow=%.3f  (nfit=%d, rmax=%.2f, %.1fs)\n",
           R12.singlet_k, R12.F_eigenvalue, R12.E_per_site, R12.xi, R12.eta,
           R12.R2_exp, R12.R2_pow, R12.nfit, R12.rmax_used, now_sec() - t12);

    printf("\n  Computing N=18 (2×3 kagome, F-filter)…\n");
    double t18 = now_sec();
    xi_result_t R18 = compute_xi(2, 3, 9, 18, 10, 200);
    printf("    singlet at k=%d (F=%+.4f)  E_0/N=%+.6f  ξ=%.4f  η=%.4f  R²_exp=%.3f  R²_pow=%.3f  (nfit=%d, rmax=%.2f, %.1fs)\n",
           R18.singlet_k, R18.F_eigenvalue, R18.E_per_site, R18.xi, R18.eta,
           R18.R2_exp, R18.R2_pow, R18.nfit, R18.rmax_used, now_sec() - t18);

    printf("\n  Computing N=24 (2×4 kagome)…\n");
    double t24 = now_sec();
    xi_result_t R24 = compute_xi(2, 4, 12, 24, 4, 200);
    printf("    singlet at k=%d (F=%+.4f)  E_0/N=%+.6f  ξ=%.4f  η=%.4f  R²_exp=%.3f  R²_pow=%.3f  (nfit=%d, rmax=%.2f, %.1fs)\n",
           R24.singlet_k, R24.F_eigenvalue, R24.E_per_site, R24.xi, R24.eta,
           R24.R2_exp, R24.R2_pow, R24.nfit, R24.rmax_used, now_sec() - t24);

    printf("\n  Finite-size scaling ξ(1/N) and η(1/N):\n");
    printf("  %-6s %-10s %-10s %-10s %-10s %-10s\n",
           "N", "1/N", "ξ", "η", "E_0/N", "nfit");
    printf("  ------ ---------- ---------- ---------- ---------- ----------\n");
    printf("  %-6d %-10.5f %-10.4f %-10.4f %-10.6f %-10d\n",
           12, 1.0/12, R12.xi, R12.eta, R12.E_per_site, R12.nfit);
    printf("  %-6d %-10.5f %-10.4f %-10.4f %-10.6f %-10d\n",
           18, 1.0/18, R18.xi, R18.eta, R18.E_per_site, R18.nfit);
    printf("  %-6d %-10.5f %-10.4f %-10.4f %-10.6f %-10d\n",
           24, 1.0/24, R24.xi, R24.eta, R24.E_per_site, R24.nfit);

    /* Three-point linear extrapolation: ξ(1/N) = ξ(∞) + α · (1/N). */
    double x1 = 1.0/12, x2 = 1.0/18, x3 = 1.0/24;
    double y1 = R12.xi, y2 = R18.xi, y3 = R24.xi;
    double mx = (x1+x2+x3)/3.0;
    double my = (y1+y2+y3)/3.0;
    double sxy = (x1-mx)*(y1-my) + (x2-mx)*(y2-my) + (x3-mx)*(y3-my);
    double sxx = (x1-mx)*(x1-mx) + (x2-mx)*(x2-mx) + (x3-mx)*(x3-mx);
    double slope = sxy / sxx;
    double xi_inf = my - slope * mx;
    printf("\n  Three-point 1/N → 0 extrapolation of ξ:\n");
    printf("    ξ(∞)   = %+.4f  (linear fit in 1/N, slope=%+.4f)\n", xi_inf, slope);

    printf("\n  INTERPRETATION:\n");
    double xi_max_N = R12.xi > R18.xi ? (R12.xi > R24.xi ? R12.xi : R24.xi)
                                       : (R18.xi > R24.xi ? R18.xi : R24.xi);
    if (xi_max_N < 2.0) {
        printf("    ✓ ξ(N) stays below 2 NN-bond-lengths for all N ∈ {12,18,24}.\n");
        printf("      This is incompatible with a gapless Dirac spin liquid,\n");
        printf("      which would require ξ(N) → ∞ as N grows. Consistent\n");
        printf("      with a gapped Z₂ spin liquid where ξ is bounded by\n");
        printf("      the (inverse) spin gap.\n");
    } else {
        printf("    ~ ξ(N) reaches %.2f bond-lengths at some N; harder to rule\n", xi_max_N);
        printf("      out algebraic decay without more cluster sizes.\n");
    }

    printf("\n  Total wall-clock: %.2f s\n", now_sec() - t0);
    return 0;
}
