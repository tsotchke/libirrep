/* SPDX-License-Identifier: MIT */
/* Full diagnostic suite on the TRUE N=24 kagome 2×4 singlet GS.
 *
 * The adversarial audit (`kagome24_audit.c`) found that the global
 * singlet GS of the 2×4 kagome torus lives at k = (0, 2), i.e. ky = π
 * along the long cluster axis — NOT at Γ. This script recomputes
 * every diagnostic from §§1.14–1.18 on the correct GS, and repeats
 * the region-family area-law γ stability check.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_true_gs
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

static void apply_bond(int i, int j,
                       const double _Complex *psi_in, double _Complex *phi_out) {
    long long mi = 1LL << i, mj = 1LL << j;
    for (long long s = 0; s < D_FULL; ++s) {
        int zi = (int)((s >> i) & 1), zj = (int)((s >> j) & 1);
        double diag = (zi == zj) ? +0.25 : -0.25;
        phi_out[s] += diag * psi_in[s];
        if (zi != zj) {
            long long sp = s ^ mi ^ mj;
            phi_out[sp] += 0.5 * psi_in[s];
        }
    }
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
    printf("=== N=24 kagome 2×4 — full diagnostics on TRUE singlet GS at k=(%d,%d) ===\n\n",
           KX_GS, KY_GS);
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);

    /* Build GS at k=(0,2). */
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, KX_GS, KY_GS);
    double _Complex chi1[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);
    printf("  sector dim at (0,2): %lld\n", sector_dim);

    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
    double eigs[4];
    double _Complex *psi_sec = malloc((size_t)4 * sector_dim * sizeof(double _Complex));
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sector_dim, 4, 200, seed, eigs, psi_sec);
    printf("  Lanczos: E_0=%+.8f, E_1=%+.8f, E_2=%+.8f, E_3=%+.8f\n",
           eigs[0], eigs[1], eigs[2], eigs[3]);
    printf("  E_0/N = %+.8f  (vs Γ-sector singlet %+.8f, diff %+.2e)\n",
           eigs[0] / N_SITES, -10.70614979 / N_SITES,
           eigs[0] - (-10.70614979));

    /* Lanczos residual. */
    double _Complex *Hpsi = malloc((size_t)sector_dim * sizeof(double _Complex));
    irrep_sg_heisenberg_sector_apply(psi_sec, Hpsi, S);
    double res2 = 0;
    for (long long i = 0; i < sector_dim; ++i) {
        double _Complex r = Hpsi[i] - eigs[0] * psi_sec[i];
        res2 += creal(r * conj(r));
    }
    printf("  ||H·ψ − E·ψ|| = %.3e\n", sqrt(res2));
    free(Hpsi);

    /* Unfold. */
    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    double _Complex *psi = malloc((size_t)D_FULL * sizeof(double _Complex));
    unfold(G, T, order, w, psi_sec, sector_dim, psi);
    double pnorm = 0;
    for (long long s = 0; s < D_FULL; ++s) pnorm += creal(psi[s] * conj(psi[s]));
    printf("  ||ψ_full||² = %.10f\n", pnorm);

    /* F-check: spin-flip parity. */
    uint64_t mask = (1ULL << N_SITES) - 1;
    double _Complex Fval = 0;
    for (long long s = 0; s < D_FULL; ++s) Fval += conj(psi[s]) * psi[s ^ mask];
    printf("  ⟨ψ|F|ψ⟩ = %+.6f %+.6ei (should be +1 for even S)\n",
           creal(Fval), cimag(Fval));

    /* ⟨S²⟩ = N·3/4 + 2·Σ_{i<j} ⟨S_i·S_j⟩. */
    printf("\n━ ⟨S²⟩ check ━\n");
    double tC = now_sec();
    double S2 = 0.75 * N_SITES;
    double *Cmat = malloc(sizeof(double) * N_SITES * N_SITES);
    for (int i = 0; i < N_SITES; ++i) Cmat[i * N_SITES + i] = 0.75;
    for (int i = 0; i < N_SITES; ++i)
        for (int j = i + 1; j < N_SITES; ++j) {
            double c = corr_si_sj(psi, i, j);
            Cmat[i * N_SITES + j] = c;
            Cmat[j * N_SITES + i] = c;
            S2 += 2.0 * c;
        }
    printf("  ⟨S²⟩ = %+.6e (computed in %.1fs)\n", S2, now_sec() - tC);

    /* ═══ DIAGNOSTIC 1: area-law γ with three region families ═══ */
    printf("\n━ Area-law γ across three region families ━\n");
    const int fam[3][12] = {
        {0,1,2,3,4,5,6,7,8,9,10,11},
        {0,2,4,6,8,10,12,14,16,18,20,22},
        {0,3,6,9,12,15,18,21, 1,4,7,10},
    };
    const char *fam_name[3] = {"contiguous", "stride-2", "mixed-sublatt"};

    double gamma_by_family[3];
    for (int f = 0; f < 3; ++f) {
        double xs[12], ys[12];
        int nA_max = 9;
        for (int nA = 1; nA <= nA_max; ++nA) {
            int A[16];
            for (int i = 0; i < nA; ++i) A[i] = fam[f][i];
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
        double R2=0, ss_res=0, ss_tot=0;
        for (int i = 0; i < npts; ++i) {
            double pred = alpha*xs[i] + beta;
            ss_res += (ys[i]-pred)*(ys[i]-pred);
            ss_tot += (ys[i]-myv)*(ys[i]-myv);
        }
        R2 = 1.0 - ss_res/ss_tot;
        gamma_by_family[f] = -beta;
        printf("  %-14s α=%+.4f β=%+.4f γ=%+.4f (R²=%.4f)\n",
               fam_name[f], alpha, beta, -beta, R2);
    }
    double g_min = gamma_by_family[0], g_max = gamma_by_family[0];
    double g_mean = 0;
    for (int f = 0; f < 3; ++f) {
        if (gamma_by_family[f] < g_min) g_min = gamma_by_family[f];
        if (gamma_by_family[f] > g_max) g_max = gamma_by_family[f];
        g_mean += gamma_by_family[f];
    }
    g_mean /= 3.0;
    printf("  γ range: [%+.4f, %+.4f]  spread=%+.4f  mean=%+.4f  (vs log 2 = %+.4f)\n",
           g_min, g_max, g_max - g_min, g_mean, log(2.0));

    /* ═══ DIAGNOSTIC 2: correlation length ξ ═══ */
    printf("\n━ Correlation length ξ ━\n");
    double site_x[N_SITES], site_y[N_SITES];
    for (int i = 0; i < N_SITES; ++i) {
        double xy[2];
        irrep_lattice_site_position(L, i, xy);
        site_x[i] = xy[0]; site_y[i] = xy[1];
    }
    double a1[2], a2[2];
    irrep_lattice_primitive_vectors(L, a1, a2);
    int Lx = irrep_lattice_Lx(L), Ly = irrep_lattice_Ly(L);

    int cap_pairs = N_SITES * N_SITES;
    double *rr = malloc(sizeof(double) * cap_pairs);
    double *cc = malloc(sizeof(double) * cap_pairs);
    int npairs = 0;
    for (int i = 0; i < N_SITES; ++i)
        for (int j = 0; j < N_SITES; ++j) {
            if (i == j) continue;
            double dx = site_x[j] - site_x[i], dy = site_y[j] - site_y[i];
            double best = 1e300;
            for (int nx = -1; nx <= 1; ++nx)
                for (int ny = -1; ny <= 1; ++ny) {
                    double ddx = dx + nx*Lx*a1[0] + ny*Ly*a2[0];
                    double ddy = dy + nx*Lx*a1[1] + ny*Ly*a2[1];
                    double r = sqrt(ddx*ddx + ddy*ddy);
                    if (r < best) best = r;
                }
            rr[npairs] = best;
            cc[npairs] = Cmat[i*N_SITES + j];
            ++npairs;
        }
    enum { NMAX = 128 };
    double uniq_r[NMAX], uniq_sum[NMAX]; int uniq_cnt[NMAX]; int nu = 0;
    for (int p = 0; p < npairs; ++p) {
        int f = -1;
        for (int u = 0; u < nu; ++u)
            if (fabs(uniq_r[u] - rr[p]) < 1e-3) { f = u; break; }
        if (f >= 0) { uniq_sum[f] += cc[p]; ++uniq_cnt[f]; }
        else { uniq_r[nu] = rr[p]; uniq_sum[nu] = cc[p]; uniq_cnt[nu] = 1; ++nu; }
    }
    for (int a = 0; a < nu - 1; ++a)
        for (int b = a + 1; b < nu; ++b)
            if (uniq_r[b] < uniq_r[a]) {
                double tr = uniq_r[a]; uniq_r[a] = uniq_r[b]; uniq_r[b] = tr;
                double ts = uniq_sum[a]; uniq_sum[a] = uniq_sum[b]; uniq_sum[b] = ts;
                int   tc = uniq_cnt[a]; uniq_cnt[a] = uniq_cnt[b]; uniq_cnt[b] = tc;
            }
    printf("  r        n_r         C(r)           |C(r)|\n");
    for (int u = 0; u < nu; ++u)
        printf("  %6.4f %5d   %+12.6e %12.6e\n",
               uniq_r[u], uniq_cnt[u],
               uniq_sum[u]/uniq_cnt[u], fabs(uniq_sum[u]/uniq_cnt[u]));

    double a1m = sqrt(a1[0]*a1[0] + a1[1]*a1[1]);
    double a2m = sqrt(a2[0]*a2[0] + a2[1]*a2[1]);
    double rmax = 0.5 * ((Lx*a1m > Ly*a2m) ? Lx*a1m : Ly*a2m);
    int nf = 0;
    double xr[NMAX], xlr[NMAX], yl[NMAX];
    for (int u = 0; u < nu; ++u) {
        double avg = uniq_sum[u]/uniq_cnt[u];
        double aC = fabs(avg);
        if (aC < 1e-6 || uniq_r[u] < 0.99 || uniq_r[u] > rmax) continue;
        xr[nf] = uniq_r[u]; xlr[nf] = log(uniq_r[u]); yl[nf] = log(aC); ++nf;
    }
    double mxv=0, myv=0, mxxv=0, mxyv=0;
    for (int i = 0; i < nf; ++i) { mxv+=xr[i]; myv+=yl[i]; mxxv+=xr[i]*xr[i]; mxyv+=xr[i]*yl[i]; }
    mxv/=nf; myv/=nf; mxxv/=nf; mxyv/=nf;
    double slope = (mxyv-mxv*myv)/(mxxv-mxv*mxv);
    double xi = -1.0/slope;
    double icpt = myv - slope*mxv;
    double rs=0,ts=0;
    for (int i = 0; i < nf; ++i) {
        double p = slope*xr[i] + icpt;
        rs += (yl[i]-p)*(yl[i]-p); ts += (yl[i]-myv)*(yl[i]-myv);
    }
    double R2e = 1.0 - rs/ts;
    mxv=myv=mxxv=mxyv=0;
    for (int i = 0; i < nf; ++i) { mxv+=xlr[i]; myv+=yl[i]; mxxv+=xlr[i]*xlr[i]; mxyv+=xlr[i]*yl[i]; }
    mxv/=nf; myv/=nf; mxxv/=nf; mxyv/=nf;
    double slope2 = (mxyv-mxv*myv)/(mxxv-mxv*mxv);
    double icpt2 = myv - slope2*mxv;
    double eta = -slope2;
    rs=ts=0;
    for (int i = 0; i < nf; ++i) {
        double p = slope2*xlr[i] + icpt2;
        rs += (yl[i]-p)*(yl[i]-p); ts += (yl[i]-myv)*(yl[i]-myv);
    }
    double R2p = 1.0 - rs/ts;
    printf("  ξ fit (nf=%d, rmax=%.2f): ξ=%.4f (R²_exp=%.4f)  η=%.4f (R²_pow=%.4f)\n",
           nf, rmax, xi, R2e, eta, R2p);

    /* ═══ DIAGNOSTIC 3: spin structure factor S(k) ═══ */
    printf("\n━ Spin structure factor S(k) ━\n");
    double b1[2], b2[2];
    irrep_lattice_reciprocal_vectors(L, b1, b2);
    double Smax = 0; int Smx = 0, Smy = 0;
    for (int my = 0; my < Ly; ++my)
        for (int mx = 0; mx < Lx; ++mx) {
            double kx = (double)mx/Lx*b1[0] + (double)my/Ly*b2[0];
            double ky = (double)mx/Lx*b1[1] + (double)my/Ly*b2[1];
            double Sk = 0;
            for (int i = 0; i < N_SITES; ++i)
                for (int j = 0; j < N_SITES; ++j) {
                    double dx = site_x[i]-site_x[j], dy = site_y[i]-site_y[j];
                    Sk += Cmat[i*N_SITES+j] * cos(kx*dx + ky*dy);
                }
            Sk /= N_SITES;
            printf("  (%d,%d)  S(k) = %+.4f\n", mx, my, Sk);
            if (fabs(Sk) > fabs(Smax)) { Smax = Sk; Smx = mx; Smy = my; }
        }
    printf("  max |S(k)|/N = %.4f at (%d,%d)\n", fabs(Smax)/N_SITES, Smx, Smy);

    /* ═══ DIAGNOSTIC 4: dimer-dimer S_D and T_disc ═══ */
    printf("\n━ Dimer diagnostics ━\n");
    double _Complex **phi = malloc((size_t)nb * sizeof(double _Complex *));
    for (int b = 0; b < nb; ++b) {
        phi[b] = calloc((size_t)D_FULL, sizeof(double _Complex));
        apply_bond(bi[b], bj[b], psi, phi[b]);
    }
    double *Bmean = malloc(sizeof(double)*nb);
    for (int b = 0; b < nb; ++b) {
        double _Complex x = 0;
        for (long long s = 0; s < D_FULL; ++s) x += conj(psi[s])*phi[b][s];
        Bmean[b] = creal(x);
    }
    double Bavg = 0;
    for (int b = 0; b < nb; ++b) Bavg += Bmean[b];
    Bavg /= nb;
    double Bvar = 0;
    for (int b = 0; b < nb; ++b) Bvar += (Bmean[b]-Bavg)*(Bmean[b]-Bavg);
    double Bstd = sqrt(Bvar/nb);
    printf("  ⟨B_b⟩ avg = %+.6f   std = %.3e   (ratio std/|avg| = %.2f%%)\n",
           Bavg, Bstd, 100.0*Bstd/fabs(Bavg));

    double *Dmat = malloc((size_t)nb*nb*sizeof(double));
    double tm = now_sec();
    for (int a = 0; a < nb; ++a)
        for (int c = a; c < nb; ++c) {
            double _Complex x = 0;
            for (long long s = 0; s < D_FULL; ++s) x += conj(phi[a][s])*phi[c][s];
            double d = creal(x) - Bmean[a]*Bmean[c];
            Dmat[a*nb+c] = d; Dmat[c*nb+a] = d;
        }
    printf("  D matrix computed in %.1fs\n", now_sec()-tm);

    double *bx = malloc(sizeof(double)*nb);
    double *by = malloc(sizeof(double)*nb);
    for (int b = 0; b < nb; ++b) {
        double xi2[2], xj[2];
        irrep_lattice_site_position(L, bi[b], xi2);
        irrep_lattice_site_position(L, bj[b], xj);
        bx[b] = 0.5*(xi2[0]+xj[0]); by[b] = 0.5*(xi2[1]+xj[1]);
    }
    double Td_g = 0, Td_max_ng = 0, Sd_max = 0;
    for (int my = 0; my < Ly; ++my)
        for (int mx = 0; mx < Lx; ++mx) {
            double kx = (double)mx/Lx*b1[0] + (double)my/Ly*b2[0];
            double ky = (double)mx/Lx*b1[1] + (double)my/Ly*b2[1];
            double _Complex acc = 0;
            for (int b = 0; b < nb; ++b)
                acc += Bmean[b] * (cos(kx*bx[b]+ky*by[b]) + I*sin(kx*bx[b]+ky*by[b]));
            double Td = creal(acc*conj(acc))/nb;
            double Sd = 0;
            for (int a = 0; a < nb; ++a)
                for (int c = 0; c < nb; ++c) {
                    double dx = bx[a]-bx[c], dy = by[a]-by[c];
                    Sd += Dmat[a*nb+c]*cos(kx*dx + ky*dy);
                }
            Sd /= nb;
            if (mx == 0 && my == 0) Td_g = Td;
            else if (Td > Td_max_ng) Td_max_ng = Td;
            if (fabs(Sd) > fabs(Sd_max)) Sd_max = Sd;
        }
    printf("  T_disc(Γ) = %.4f   max T_disc(k≠Γ) = %.4f   ratio = %.4f\n",
           Td_g, Td_max_ng, Td_max_ng/Td_g);
    printf("  max |S_D(k)|/N_bonds = %.4f\n", fabs(Sd_max)/nb);

    /* Cleanup. */
    for (int b = 0; b < nb; ++b) free(phi[b]);
    free(phi); free(Bmean); free(Dmat); free(bx); free(by);
    free(rr); free(cc); free(Cmat);
    free(psi); free(psi_sec); free(seed); free(w);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
