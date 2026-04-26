/* SPDX-License-Identifier: MIT */
/* Adversarial audit of the N=24 kagome-Heisenberg spin-liquid claim.
 * Independent checks of the ansatz assumptions and numerical quality.
 *
 *   A — ⟨S²⟩ eigenvalue check: F=+1 singlet identification
 *   B — Lanczos residual ||H·ψ − E·ψ||²
 *   C — E(k) sweep across the full 2×4 Brillouin zone with F-filter
 *         → is Γ the global GS?
 *   D — Area-law γ with three different contiguous region families
 *         → region-choice stability of γ
 *   F — 5-seed Lanczos bootstrap → GS basin-uniqueness
 *   H — ξ extraction across three (rmin, rmax) windows
 *         → fit-window sensitivity of the correlation length
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_audit
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

static int boundary_bonds(const int *A, int nA, const int *bi, const int *bj, int nb) {
    char in_A[64] = {0};
    for (int k = 0; k < nA; ++k) in_A[A[k]] = 1;
    int boundary = 0;
    for (int b = 0; b < nb; ++b)
        if (in_A[bi[b]] != in_A[bj[b]]) ++boundary;
    return boundary;
}

/* Build GS at momentum (kx, ky), F-filter. Returns first F=+1 state. */
typedef struct {
    int kx, ky;
    double E;
    double F;
    int    kk; /* which Lanczos state was picked */
} k_gs_t;

static k_gs_t lanczos_gs_at_k_with_F(const irrep_heisenberg_t *H,
                                     const irrep_space_group_t *G,
                                     const irrep_sg_rep_table_t *T,
                                     int kx, int ky, int K_eigs, int max_it,
                                     double seed_seed,
                                     double _Complex *psi_full_out /* 2^N or NULL */,
                                     double *E_out, int *sector_dim_out) {
    k_gs_t R = { kx, ky, +INFINITY, 0.0, -1 };

    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, kx, ky);
    double _Complex chi1[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    if (!S) {
        irrep_sg_little_group_irrep_free(mu);
        irrep_sg_little_group_free(lg);
        return R;
    }
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);
    if (sector_dim_out) *sector_dim_out = (int)sector_dim;

    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.1 * sin(seed_seed * i) + I * 0.05 * cos((seed_seed + 1.3) * i);
    double *eigs_all = malloc((size_t)K_eigs * sizeof(double));
    double _Complex *psi_all =
        malloc((size_t)K_eigs * sector_dim * sizeof(double _Complex));
    int it = max_it > sector_dim ? (int)sector_dim : max_it;
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sector_dim, K_eigs, it, seed, eigs_all, psi_all);

    /* F-filter: unfold each, measure F, pick first F=+1. */
    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    double _Complex *psi_full = malloc((size_t)D_FULL * sizeof(double _Complex));
    uint64_t mask = (1ULL << N_SITES) - 1;

    for (int kk = 0; kk < K_eigs; ++kk) {
        unfold(G, T, order, w, psi_all + (size_t)kk * sector_dim,
               sector_dim, psi_full);
        double _Complex Fk = 0;
        for (long long s = 0; s < D_FULL; ++s) Fk += conj(psi_full[s]) * psi_full[s ^ mask];
        if (creal(Fk) > 0.5) {
            R.kk = kk;
            R.E  = eigs_all[kk];
            R.F  = creal(Fk);
            if (psi_full_out) memcpy(psi_full_out, psi_full, (size_t)D_FULL * sizeof(double _Complex));
            if (E_out) *E_out = eigs_all[kk];
            break;
        }
    }

    free(psi_full); free(psi_all); free(eigs_all); free(seed); free(w);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
    return R;
}

int main(void) {
    printf("=== Kagome N=24 ADVERSARIAL AUDIT ===\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);
    int Lx = irrep_lattice_Lx(L), Ly = irrep_lattice_Ly(L);

    /* --- Check C: E(k) sweep across all (kx, ky) in 2×4 BZ. -------- */
    printf("━━━ Check C: E₀(k) sweep across 8 BZ points (F-filter singlet) ━━━\n");
    printf("  %-8s %-10s %-5s %-14s %-10s %-6s\n",
           "(kx,ky)", "dim", "kk", "E₀", "E₀/N", "F");

    double E_gamma = +INFINITY, E_Gamma_global = 0;
    double _Complex *psi_gamma = malloc((size_t)D_FULL * sizeof(double _Complex));
    double *E_by_k = malloc(sizeof(double) * Lx * Ly);

    for (int my = 0; my < Ly; ++my)
        for (int mx = 0; mx < Lx; ++mx) {
            int gamma = (mx == 0 && my == 0);
            double _Complex *pfull = gamma ? psi_gamma : NULL;
            double E;
            int sdim;
            double tk = now_sec();
            k_gs_t R = lanczos_gs_at_k_with_F(H, G, T, mx, my, 6, 150, 0.37,
                                              pfull, &E, &sdim);
            double dt = now_sec() - tk;
            E_by_k[my * Lx + mx] = R.E;
            if (gamma) { E_gamma = R.E; E_Gamma_global = R.E; }
            printf("  (%d,%d)    %-10d %-5d %+14.8f %+10.6f %+6.4f  [%.1fs]\n",
                   mx, my, sdim, R.kk, R.E, R.E / N_SITES, R.F, dt);
        }

    double E_min = E_gamma;
    int    kmin_mx = 0, kmin_my = 0;
    for (int my = 0; my < Ly; ++my)
        for (int mx = 0; mx < Lx; ++mx)
            if (E_by_k[my * Lx + mx] < E_min) {
                E_min = E_by_k[my * Lx + mx]; kmin_mx = mx; kmin_my = my;
            }
    printf("\n  → global minimum singlet E₀ = %+.8f at k = (%d, %d)\n",
           E_min, kmin_mx, kmin_my);
    printf("  → Γ point singlet E₀          = %+.8f\n", E_gamma);
    printf("  → Γ − global = %+.2e\n", E_gamma - E_min);
    if (fabs(E_gamma - E_min) < 1e-6)
        printf("  ✓ Γ IS the global GS; ansatz's Γ-restriction is justified.\n");
    else if (E_gamma > E_min)
        printf("  ✗ Γ is NOT the global GS — true singlet at (%d,%d) is %+.2e lower.\n",
               kmin_mx, kmin_my, E_gamma - E_min);
    else
        printf("  (non-standard) Γ lies above the other-k singlet by %.2e\n", E_min - E_gamma);

    /* From here on work with the Γ GS psi_gamma. */
    double _Complex *psi = psi_gamma;
    double E0 = E_gamma;
    double pnorm = 0;
    for (long long s = 0; s < D_FULL; ++s) pnorm += creal(psi[s] * conj(psi[s]));
    printf("  ||ψ_Γ||² = %.10f  (should be 1.0)\n", pnorm);

    /* --- Check A: ⟨S²⟩ eigenvalue. ---------------------------------- */
    printf("\n━━━ Check A: total spin ⟨S²⟩ for F=+1 GS ━━━\n");
    printf("  Definition: S² = (Σ_i S_i)² = Σ_i ⟨S_i²⟩ + Σ_{i≠j} ⟨S_i·S_j⟩.\n");
    printf("  Exact values: singlet S=0 ⇒ 0; triplet S=1 ⇒ 2; S=2 ⇒ 6.\n");
    double tA = now_sec();
    double S2 = 0.75 * N_SITES;
    for (int i = 0; i < N_SITES; ++i)
        for (int j = i + 1; j < N_SITES; ++j) {
            double c = corr_si_sj(psi, i, j);
            S2 += 2.0 * c;  /* symmetric pair */
        }
    printf("  ⟨S²⟩ = %+.6e   (computed in %.1fs)\n", S2, now_sec() - tA);
    if (fabs(S2) < 1e-4)
        printf("  ✓ S² ≈ 0 → genuine singlet, F=+1 filter is confirmed correct.\n");
    else if (fabs(S2 - 2.0) < 1e-3)
        printf("  ✗ S² ≈ 2 → TRIPLET disguised as F=+1! Ansatz is wrong.\n");
    else
        printf("  ~ S² = %.4f → intermediate/unclean. Inspect.\n", S2);

    /* --- Check B: Lanczos residual ||H·ψ_sector − E·ψ_sector||². ---- */
    printf("\n━━━ Check B: Lanczos residual at Γ ━━━\n");
    irrep_sg_little_group_t *lg_g = irrep_sg_little_group_build(G, 0, 0);
    double _Complex chi1[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu_g = irrep_sg_little_group_irrep_new(lg_g, chi1, 1);
    irrep_sg_heisenberg_sector_t *S_g = irrep_sg_heisenberg_sector_build_at_k(H, T, lg_g, mu_g);
    long long sdim_g = irrep_sg_heisenberg_sector_dim(S_g);
    double _Complex *seed = malloc((size_t)sdim_g * sizeof(double _Complex));
    for (long long i = 0; i < sdim_g; ++i)
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
    double eigs_B[4];
    double _Complex *psi_B_all = malloc((size_t)4 * sdim_g * sizeof(double _Complex));
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S_g,
                                 sdim_g, 4, 200, seed, eigs_B, psi_B_all);
    double _Complex *Hpsi = malloc((size_t)sdim_g * sizeof(double _Complex));
    irrep_sg_heisenberg_sector_apply(psi_B_all, Hpsi, S_g);
    double res2 = 0;
    for (long long i = 0; i < sdim_g; ++i) {
        double _Complex r = Hpsi[i] - eigs_B[0] * psi_B_all[i];
        res2 += creal(r * conj(r));
    }
    printf("  ||H·ψ − E·ψ|| = %.3e  (sector-basis, max_it=200)\n", sqrt(res2));
    if (sqrt(res2) < 1e-6) printf("  ✓ Lanczos converged to machine-precision GS.\n");
    else if (sqrt(res2) < 1e-3) printf("  ~ Loose convergence; increase max_it.\n");
    else printf("  ✗ Bad convergence — GS not reliably identified.\n");
    free(Hpsi);

    /* --- Check F: 5-seed Lanczos bootstrap. ------------------------ */
    printf("\n━━━ Check F: 5-seed Lanczos bootstrap at Γ ━━━\n");
    printf("  %-4s %-14s %-14s %-10s\n", "seed", "E₀", "ΔE", "|⟨ψ_seed|ψ_B[0]⟩|");
    printf("  ---- -------------- -------------- ----------\n");
    double seed_seeds[5] = {0.37, 0.13, 0.71, 1.01, 2.23};
    for (int s_i = 0; s_i < 5; ++s_i) {
        double ss = seed_seeds[s_i];
        for (long long i = 0; i < sdim_g; ++i)
            seed[i] = 0.1 * sin(ss * i) + I * 0.05 * cos((ss + 1.3) * i);
        double eigs_s[2];
        double _Complex *psi_s = malloc((size_t)2 * sdim_g * sizeof(double _Complex));
        irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S_g,
                                     sdim_g, 2, 200, seed, eigs_s, psi_s);
        double _Complex ov = 0;
        for (long long i = 0; i < sdim_g; ++i) ov += conj(psi_B_all[i]) * psi_s[i];
        printf("  %-4d %+14.8f %+14.2e %+10.6f\n",
               s_i, eigs_s[0], eigs_s[0] - eigs_B[0], cabs(ov));
        free(psi_s);
    }
    printf("  (|overlap| ≈ 1 ⇒ unique GS; |overlap| < 1 ⇒ degeneracy)\n");

    free(seed); free(psi_B_all);
    irrep_sg_heisenberg_sector_free(S_g);
    irrep_sg_little_group_irrep_free(mu_g);
    irrep_sg_little_group_free(lg_g);

    /* --- Check D: area-law γ with 3 region families. --------------- */
    printf("\n━━━ Check D: area-law γ stability under region-family choice ━━━\n");
    /* Family 1: contiguous [0, |A|)  (this is our default).
     * Family 2: every-other-site (stride 2): 0, 2, 4, …
     * Family 3: first half of each translational orbit: 0, 3, 6, … (kagome has 3 sublattices, stride 3)
     */
    const int fam[3][12] = {
        {0,1,2,3,4,5,6,7,8,9,10,11},        /* contiguous */
        {0,2,4,6,8,10,12,14,16,18,20,22},   /* stride 2 */
        {0,3,6,9,12,15,18,21, 1,4,7,10},    /* mixed sublattice */
    };
    const char *fam_name[3] = {"contiguous   ", "stride-2     ", "mixed-sublatt"};
    int fam_len[3] = {12, 12, 12};

    for (int f = 0; f < 3; ++f) {
        printf("\n  family: %s\n", fam_name[f]);
        printf("  %-5s %-6s %12s\n", "|A|", "|∂A|", "S_A");
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
            printf("    %-3d %-6d %+12.6f\n", nA, bdry, S_A);
            xs[nA - 1] = bdry; ys[nA - 1] = S_A;
            free(rho);
        }
        int npts = nA_max;
        double mx=0, my=0, mxx=0, mxy=0;
        for (int i = 0; i < npts; ++i) { mx+=xs[i]; my+=ys[i]; mxx+=xs[i]*xs[i]; mxy+=xs[i]*ys[i]; }
        mx /= npts; my /= npts; mxx /= npts; mxy /= npts;
        double alpha = (mxy - mx*my) / (mxx - mx*mx);
        double beta  = my - alpha * mx;
        double gamma_ext = -beta;
        double ss_res=0, ss_tot=0;
        for (int i = 0; i < npts; ++i) {
            double pred = alpha * xs[i] + beta;
            ss_res += (ys[i]-pred)*(ys[i]-pred);
            ss_tot += (ys[i]-my)*(ys[i]-my);
        }
        double R2 = 1.0 - ss_res/ss_tot;
        (void)fam_len;
        printf("  ⇒ α = %+.4f,  β = %+.4f,  γ = %+.4f  (R² = %.4f)\n",
               alpha, beta, gamma_ext, R2);
    }
    printf("\n  Consistency interpretation:\n");
    printf("    If γ spread ≫ log 2 = 0.69 across families,\n");
    printf("    the γ extraction is region-choice dominated, not physical.\n");

    /* --- Check H: ξ fit-window sensitivity. ------------------------ */
    printf("\n━━━ Check H: ξ fit-window sensitivity ━━━\n");
    /* Re-use C_{ij} on psi, then fit with different (rmin, rmax). */
    double site_x[N_SITES], site_y[N_SITES];
    for (int i = 0; i < N_SITES; ++i) {
        double xy[2];
        irrep_lattice_site_position(L, i, xy);
        site_x[i] = xy[0]; site_y[i] = xy[1];
    }
    double a1[2], a2[2];
    irrep_lattice_primitive_vectors(L, a1, a2);
    /* Precomputed all pairwise C_ij. */
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
                    double ddx = dx + nx * Lx * a1[0] + ny * Ly * a2[0];
                    double ddy = dy + nx * Lx * a1[1] + ny * Ly * a2[1];
                    double r = sqrt(ddx*ddx + ddy*ddy);
                    if (r < best) best = r;
                }
            rr[npairs] = best;
            cc[npairs] = corr_si_sj(psi, i, j);
            ++npairs;
        }
    enum { NMAX = 128 };
    double uniq_r[NMAX], uniq_sum[NMAX]; int uniq_cnt[NMAX];
    int nu = 0;
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
    double windows[4][2] = { {0.99, 2.00}, {0.99, 3.00}, {0.99, 4.00}, {1.50, 4.00} };
    const char *window_name[4] = { "0.99–2.00  (short-range only)",
                                    "0.99–3.00  (intermediate)",
                                    "0.99–4.00  (full clean range)",
                                    "1.50–4.00  (skip r=1, tail only)" };
    printf("  %-36s %-8s %-8s %-8s %-8s\n", "window", "ξ", "R²_exp", "η", "R²_pow");
    for (int k = 0; k < 4; ++k) {
        double rmin = windows[k][0], rmax = windows[k][1];
        int nf = 0;
        double xr[NMAX], xlr[NMAX], yl[NMAX];
        for (int u = 0; u < nu; ++u) {
            double avg = uniq_sum[u] / uniq_cnt[u];
            double aC = fabs(avg);
            if (aC < 1e-6) continue;
            if (uniq_r[u] < rmin || uniq_r[u] > rmax) continue;
            xr[nf] = uniq_r[u];
            xlr[nf] = log(uniq_r[u]);
            yl[nf] = log(aC);
            ++nf;
        }
        double mxx=0, myy=0, mxxm=0, mxym=0;
        for (int i = 0; i < nf; ++i) { mxx+=xr[i]; myy+=yl[i]; mxxm+=xr[i]*xr[i]; mxym+=xr[i]*yl[i]; }
        mxx/=nf; myy/=nf; mxxm/=nf; mxym/=nf;
        double slope = (mxym - mxx*myy) / (mxxm - mxx*mxx);
        double icpt  = myy - slope*mxx;
        double xi = -1.0/slope;
        double rs=0, ts=0;
        for (int i = 0; i < nf; ++i) {
            double p = slope*xr[i] + icpt;
            rs += (yl[i]-p)*(yl[i]-p);
            ts += (yl[i]-myy)*(yl[i]-myy);
        }
        double R2e = 1.0 - rs/ts;
        mxx=myy=mxxm=mxym=0;
        for (int i = 0; i < nf; ++i) { mxx+=xlr[i]; myy+=yl[i]; mxxm+=xlr[i]*xlr[i]; mxym+=xlr[i]*yl[i]; }
        mxx/=nf; myy/=nf; mxxm/=nf; mxym/=nf;
        double slope2 = (mxym - mxx*myy) / (mxxm - mxx*mxx);
        double icpt2  = myy - slope2*mxx;
        double eta = -slope2;
        rs=ts=0;
        for (int i = 0; i < nf; ++i) {
            double p = slope2*xlr[i] + icpt2;
            rs += (yl[i]-p)*(yl[i]-p);
            ts += (yl[i]-myy)*(yl[i]-myy);
        }
        double R2p = 1.0 - rs/ts;
        printf("  %-36s %-8.4f %-8.4f %-8.4f %-8.4f  (nf=%d)\n",
               window_name[k], xi, R2e, eta, R2p, nf);
    }
    printf("  (ξ spread ≫ 0.2 ⇒ fit-window dominates; qualitative gap conclusion survives if R²_exp > R²_pow consistently.)\n");

    free(rr); free(cc);
    free(psi_gamma); free(E_by_k);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);

    printf("\n━━━ AUDIT total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
