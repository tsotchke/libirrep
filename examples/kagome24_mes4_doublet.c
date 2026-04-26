/* SPDX-License-Identifier: MIT */
/* Kagome 2×4 N=24 4-STATE MES γ on the DOUBLET-AUGMENTED 4-state
 * candidate manifold:
 *
 *   |ψ_0⟩ = singlet GS at k = (0, 0),                 F = +1
 *   |ψ_1⟩ = singlet GS at k = (0, π),                 F = +1
 *   |ψ_2⟩ = (|ψ_(0,1)⟩ + |ψ_(0,3)⟩)/√2, F-symmetrised
 *   |ψ_3⟩ = (|ψ_(0,1)⟩ − |ψ_(0,3)⟩)/√2, F-symmetrised
 *
 * The unfiltered tower scan found (0,1) and (0,3) as a degenerate
 * doublet pair under ky → −ky, both at E = −10.65970 with
 * ⟨F⟩ = −0.30. Symmetric / antisymmetric combinations give
 * F = ±1 eigenstates spanning this 2-dim doublet sector.
 *
 * User-suggested investigation: the Z₂ topological 4-state manifold
 * typically decomposes as 2 singlets + 1 doublet-pair, not 4 singlets
 * at BZ vertices. This is the correct 4-state manifold candidate.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_mes4_doublet
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
#define NS 4

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

/* Build lowest F=+1 state (singlet GS) at given k. */
static void build_gs_at_k_singlet(const irrep_heisenberg_t *H,
                                  const irrep_space_group_t *G,
                                  const irrep_sg_rep_table_t *T,
                                  int kx, int ky, int K_eigs,
                                  double _Complex *psi_full_out,
                                  double *E_out, double *F_out) {
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, kx, ky);
    double _Complex chi1[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sdim = irrep_sg_heisenberg_sector_dim(S);

    double _Complex *seed = malloc((size_t)sdim * sizeof(double _Complex));
    for (long long i = 0; i < sdim; ++i)
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
    double *eigs = malloc((size_t)K_eigs * sizeof(double));
    double _Complex *psi_all = malloc((size_t)K_eigs * sdim * sizeof(double _Complex));
    int it = 200 > sdim ? (int)sdim : 200;
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sdim, K_eigs, it, seed, eigs, psi_all);
    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    double _Complex *psi_scan = malloc((size_t)D_FULL * sizeof(double _Complex));
    uint64_t mask = (1ULL << N_SITES) - 1;
    int found = -1;
    for (int k = 0; k < K_eigs; ++k) {
        unfold(G, T, order, w, psi_all + (size_t)k*sdim, sdim, psi_scan);
        double _Complex Fk = 0;
        for (long long s = 0; s < D_FULL; ++s) Fk += conj(psi_scan[s]) * psi_scan[s ^ mask];
        if (creal(Fk) > 0.5) {
            memcpy(psi_full_out, psi_scan, (size_t)D_FULL * sizeof(double _Complex));
            *E_out = eigs[k]; *F_out = creal(Fk);
            found = k; break;
        }
    }
    if (found < 0) { *E_out = eigs[0]; *F_out = 0; memcpy(psi_full_out, psi_scan, (size_t)D_FULL * sizeof(double _Complex)); }
    free(seed); free(eigs); free(psi_all); free(w); free(psi_scan);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
}

/* Build lowest state at given k (no F filter). */
static void build_gs_at_k_lowest(const irrep_heisenberg_t *H,
                                 const irrep_space_group_t *G,
                                 const irrep_sg_rep_table_t *T,
                                 int kx, int ky,
                                 double _Complex *psi_full_out,
                                 double *E_out, double *F_out) {
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, kx, ky);
    double _Complex chi1[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sdim = irrep_sg_heisenberg_sector_dim(S);

    double _Complex *seed = malloc((size_t)sdim * sizeof(double _Complex));
    for (long long i = 0; i < sdim; ++i)
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
    double eigs[2];
    double _Complex *psi_all = malloc((size_t)2 * sdim * sizeof(double _Complex));
    int it = 200 > sdim ? (int)sdim : 200;
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sdim, 2, it, seed, eigs, psi_all);
    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    unfold(G, T, order, w, psi_all, sdim, psi_full_out);
    uint64_t mask = (1ULL << N_SITES) - 1;
    double _Complex Fk = 0;
    for (long long s = 0; s < D_FULL; ++s)
        Fk += conj(psi_full_out[s]) * psi_full_out[s ^ mask];
    *E_out = eigs[0]; *F_out = creal(Fk);
    free(seed); free(psi_all); free(w);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
}

static void cross_rdm(const double _Complex *psi_i,
                      const double _Complex *psi_j,
                      const int *A, int nA,
                      double _Complex *rho_out) {
    long long dA = 1LL << nA;
    long long dE = D_FULL >> nA;
    memset(rho_out, 0, (size_t)dA * dA * sizeof(double _Complex));
    int in_A[64] = {0}, pos_A[64] = {0}, pos_E[64] = {0};
    int ia = 0, ie = 0;
    for (int k = 0; k < nA; ++k) { in_A[A[k]] = 1; pos_A[A[k]] = ia++; }
    for (int s = 0; s < N_SITES; ++s) if (!in_A[s]) pos_E[s] = ie++;
    long long *s_from_ae = malloc((size_t)dA * dE * sizeof(long long));
    for (long long s = 0; s < D_FULL; ++s) {
        long long a = 0, e = 0;
        for (int k = 0; k < N_SITES; ++k) {
            int bit = (int)((s >> k) & 1);
            if (in_A[k]) a |= ((long long)bit) << pos_A[k];
            else          e |= ((long long)bit) << pos_E[k];
        }
        s_from_ae[a * dE + e] = s;
    }
    for (long long a = 0; a < dA; ++a)
        for (long long b = 0; b < dA; ++b) {
            double _Complex acc = 0;
            for (long long e = 0; e < dE; ++e)
                acc += conj(psi_i[s_from_ae[a * dE + e]]) *
                       psi_j[s_from_ae[b * dE + e]];
            rho_out[a * dA + b] = acc;
        }
    free(s_from_ae);
}

static double entropy_ns(double _Complex *const *r, int N, int nA,
                         const double _Complex *c) {
    long long dA = 1LL << nA;
    double _Complex *rho = malloc((size_t)dA * dA * sizeof(double _Complex));
    memset(rho, 0, (size_t)dA * dA * sizeof(double _Complex));
    for (int a = 0; a < N; ++a)
        for (int b = 0; b < N; ++b) {
            double _Complex w = conj(c[a]) * c[b];
            for (long long k = 0; k < dA * dA; ++k)
                rho[k] += w * r[a * N + b][k];
        }
    double S = irrep_entropy_vonneumann(rho, (int)dA);
    free(rho);
    return S;
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
    printf("=== Kagome 2×4 N=24: 4-state MES γ on DOUBLET-AUGMENTED manifold ===\n");
    printf("    |ψ_0⟩ = singlet at (0,0)\n");
    printf("    |ψ_1⟩ = singlet at (0,π)\n");
    printf("    |ψ_2⟩, |ψ_3⟩ = sym/anti combinations of (0,1), (0,3) doublet\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);

    /* Build the 2 singlet GSes and 2 doublet states. */
    double _Complex *psi = malloc((size_t)NS * D_FULL * sizeof(double _Complex));
    double _Complex *psi_01 = malloc((size_t)D_FULL * sizeof(double _Complex));
    double _Complex *psi_03 = malloc((size_t)D_FULL * sizeof(double _Complex));
    double Es[NS], Fs[NS];

    printf("  Building singlet GS at (0,0)…\n"); fflush(stdout);
    double tG = now_sec();
    build_gs_at_k_singlet(H, G, T, 0, 0, 6, psi + 0*D_FULL, &Es[0], &Fs[0]);
    printf("    E=%+.8f, F=%+.4f (%.1fs)\n", Es[0], Fs[0], now_sec() - tG); fflush(stdout);

    printf("  Building singlet GS at (0,π)…\n"); fflush(stdout);
    tG = now_sec();
    build_gs_at_k_singlet(H, G, T, 0, 2, 6, psi + 1*D_FULL, &Es[1], &Fs[1]);
    printf("    E=%+.8f, F=%+.4f (%.1fs)\n", Es[1], Fs[1], now_sec() - tG); fflush(stdout);

    printf("  Building lowest state at (0,1)…\n"); fflush(stdout);
    tG = now_sec();
    double E_01, F_01;
    build_gs_at_k_lowest(H, G, T, 0, 1, psi_01, &E_01, &F_01);
    printf("    E=%+.8f, F=%+.4f (%.1fs)\n", E_01, F_01, now_sec() - tG); fflush(stdout);

    printf("  Building lowest state at (0,3)…\n"); fflush(stdout);
    tG = now_sec();
    double E_03, F_03;
    build_gs_at_k_lowest(H, G, T, 0, 3, psi_03, &E_03, &F_03);
    printf("    E=%+.8f, F=%+.4f (%.1fs)\n", E_03, F_03, now_sec() - tG); fflush(stdout);

    /* Symmetric / antisymmetric combinations of the doublet. */
    double _Complex *psi_sym = psi + 2*D_FULL;
    double _Complex *psi_ant = psi + 3*D_FULL;
    double sqrt2 = sqrt(2.0);
    for (long long s = 0; s < D_FULL; ++s) {
        psi_sym[s] = (psi_01[s] + psi_03[s]) / sqrt2;
        psi_ant[s] = (psi_01[s] - psi_03[s]) / sqrt2;
    }

    /* Normalisations of sym/anti. */
    double n_sym = 0, n_ant = 0;
    for (long long s = 0; s < D_FULL; ++s) {
        n_sym += creal(psi_sym[s] * conj(psi_sym[s]));
        n_ant += creal(psi_ant[s] * conj(psi_ant[s]));
    }
    printf("  ||ψ_sym||² = %.6f,  ||ψ_ant||² = %.6f\n", n_sym, n_ant);
    /* Re-normalize. */
    double sf_sym = 1.0/sqrt(n_sym);
    double sf_ant = 1.0/sqrt(n_ant);
    for (long long s = 0; s < D_FULL; ++s) {
        psi_sym[s] *= sf_sym;
        psi_ant[s] *= sf_ant;
    }

    /* F-eigenvalue measurements on sym/anti. */
    uint64_t mask = (1ULL << N_SITES) - 1;
    double _Complex Fsym = 0, Fant = 0;
    for (long long s = 0; s < D_FULL; ++s) {
        Fsym += conj(psi_sym[s]) * psi_sym[s ^ mask];
        Fant += conj(psi_ant[s]) * psi_ant[s ^ mask];
    }
    Es[2] = E_01;  /* both doublet states are at E_01 = E_03 */
    Es[3] = E_03;
    Fs[2] = creal(Fsym);
    Fs[3] = creal(Fant);
    printf("  ⟨ψ_sym|F|ψ_sym⟩ = %+.4f\n", creal(Fsym));
    printf("  ⟨ψ_ant|F|ψ_ant⟩ = %+.4f\n", creal(Fant));
    fflush(stdout);

    free(psi_01); free(psi_03);

    /* Orthogonality check. */
    printf("\n  Orthogonality of the 4-state basis:\n");
    for (int a = 0; a < NS; ++a) {
        for (int b = 0; b < NS; ++b) {
            double _Complex ov = 0;
            for (long long s = 0; s < D_FULL; ++s)
                ov += conj(psi[(size_t)a*D_FULL+s]) * psi[(size_t)b*D_FULL+s];
            printf("  %.3e", cabs(ov));
        }
        printf("\n");
    }

    /* Region series |A| ∈ {2, 3, 4, 5}. */
    int As[4] = {2, 3, 4, 5};
    int nsizes = 4;
    int bdry[4];
    double S_mes_arr[4], S_ind_arr[4];

    for (int si = 0; si < nsizes; ++si) {
        int nA = As[si];
        int reg[8];
        for (int i = 0; i < nA; ++i) reg[i] = i;
        bdry[si] = boundary_bonds(reg, nA, bi, bj, nb);
        printf("\n━ |A|=%d, |∂A|=%d ━\n", nA, bdry[si]); fflush(stdout);
        long long dA = 1LL << nA;

        double _Complex **r = malloc((size_t)NS*NS * sizeof(double _Complex *));
        double tf = now_sec();
        for (int a = 0; a < NS; ++a)
            for (int b = 0; b < NS; ++b) {
                r[a * NS + b] = malloc((size_t)dA * dA * sizeof(double _Complex));
                cross_rdm(psi + (size_t)a * D_FULL,
                          psi + (size_t)b * D_FULL, reg, nA, r[a * NS + b]);
            }
        printf("  cross-RDMs (%.1fs)\n", now_sec() - tf);
        fflush(stdout);

        double S_i[NS];
        for (int a = 0; a < NS; ++a) {
            double _Complex c[NS] = {0};
            c[a] = 1;
            S_i[a] = entropy_ns(r, NS, nA, c);
        }
        double S_ind_mean = 0;
        for (int a = 0; a < NS; ++a) S_ind_mean += S_i[a];
        S_ind_mean /= NS;
        S_ind_arr[si] = S_ind_mean;
        printf("  S_i:");
        for (int a = 0; a < NS; ++a) printf(" %.4f", S_i[a]);
        printf("  mean %.4f\n", S_ind_mean);

        double tmes = now_sec();
        int Nth = 6, Nph = 6;
        double S_min = +INFINITY, S_max = -INFINITY;
        double c_min[NS] = {0};
        for (int it1 = 0; it1 <= Nth; ++it1) {
            double th1 = it1 * 0.5 * M_PI / Nth;
            for (int it2 = 0; it2 <= Nth; ++it2) {
                double th2 = it2 * 0.5 * M_PI / Nth;
                for (int it3 = 0; it3 <= Nth; ++it3) {
                    double th3 = it3 * 0.5 * M_PI / Nth;
                    for (int ip1 = 0; ip1 < Nph; ++ip1) {
                        double ph1 = ip1 * 2.0 * M_PI / Nph;
                        for (int ip2 = 0; ip2 < Nph; ++ip2) {
                            double ph2 = ip2 * 2.0 * M_PI / Nph;
                            for (int ip3 = 0; ip3 < Nph; ++ip3) {
                                double ph3 = ip3 * 2.0 * M_PI / Nph;
                                double _Complex c[NS];
                                c[0] = cos(th1);
                                c[1] = sin(th1) * cos(th2) * cexp(I * ph1);
                                c[2] = sin(th1) * sin(th2) * cos(th3) * cexp(I * ph2);
                                c[3] = sin(th1) * sin(th2) * sin(th3) * cexp(I * ph3);
                                double S = entropy_ns(r, NS, nA, c);
                                if (S < S_min) {
                                    S_min = S;
                                    for (int k = 0; k < NS; ++k) c_min[k] = cabs(c[k]);
                                }
                                if (S > S_max) S_max = S;
                            }
                        }
                    }
                }
            }
        }
        S_mes_arr[si] = S_min;
        printf("  MES: S_min=%.4f  S_max=%.4f  range=%.4f  (%.1fs)\n",
               S_min, S_max, S_max - S_min, now_sec() - tmes);
        printf("  S_min − S_ind_mean = %+.4f   (log 2 = %.4f)\n",
               S_min - S_ind_mean, log(2.0));
        printf("  c_MES ≈ (");
        for (int k = 0; k < NS; ++k)
            printf("%+.3f%s", c_min[k], k < NS-1 ? " " : "");
        printf(")\n");
        fflush(stdout);

        for (int a = 0; a < NS * NS; ++a) free(r[a]);
        free(r);
    }

    /* Fit. */
    printf("\n━ Linear fit S vs |∂A| ━\n");
    printf("  |A|  |∂A|  S_individual_mean  S_MES\n");
    for (int si = 0; si < nsizes; ++si)
        printf("  %-3d  %-4d  %+.4f             %+.4f\n",
               As[si], bdry[si], S_ind_arr[si], S_mes_arr[si]);

    double mx=0, myi=0, mym=0, mxx=0, mxyi=0, mxym=0;
    for (int si = 0; si < nsizes; ++si) {
        mx   += bdry[si]; myi  += S_ind_arr[si]; mym  += S_mes_arr[si];
        mxx  += (double)bdry[si] * bdry[si];
        mxyi += (double)bdry[si] * S_ind_arr[si];
        mxym += (double)bdry[si] * S_mes_arr[si];
    }
    mx /= nsizes; myi /= nsizes; mym /= nsizes;
    mxx /= nsizes; mxyi /= nsizes; mxym /= nsizes;
    double alpha_i = (mxyi - mx*myi) / (mxx - mx*mx);
    double beta_i  = myi - alpha_i * mx;
    double alpha_m = (mxym - mx*mym) / (mxx - mx*mx);
    double beta_m  = mym - alpha_m * mx;
    double ss_res_i=0, ss_tot_i=0, ss_res_m=0, ss_tot_m=0;
    for (int si = 0; si < nsizes; ++si) {
        double p_i = alpha_i * bdry[si] + beta_i;
        double p_m = alpha_m * bdry[si] + beta_m;
        ss_res_i += (S_ind_arr[si] - p_i) * (S_ind_arr[si] - p_i);
        ss_tot_i += (S_ind_arr[si] - myi) * (S_ind_arr[si] - myi);
        ss_res_m += (S_mes_arr[si] - p_m) * (S_mes_arr[si] - p_m);
        ss_tot_m += (S_mes_arr[si] - mym) * (S_mes_arr[si] - mym);
    }
    printf("\n  Individual-state fit:\n");
    printf("    S_ind = %+.4f |∂A| + %+.4f   γ_ind = %+.4f   (R²=%.4f)\n",
           alpha_i, beta_i, -beta_i, 1 - ss_res_i/ss_tot_i);
    printf("  MES fit:\n");
    printf("    S_MES = %+.4f |∂A| + %+.4f   γ_MES = %+.4f   (R²=%.4f)\n",
           alpha_m, beta_m, -beta_m, 1 - ss_res_m/ss_tot_m);
    printf("\n  γ_ind − γ_MES = %+.4f\n", -beta_i - (-beta_m));
    printf("\n  CALIBRATION CONTROL (Kitaev A-phase, §1.20):\n");
    printf("    4-state Z₂ manifold: γ_ind=+1.353, γ_MES=+1.257, Δ=0.097 (14%% of log 2)\n");
    printf("    γ_MES recovers 2·log 2 (= +1.386) within 9.3%%.\n");
    printf("\n  Kagome 4-state (BZ vertices, earlier): γ_ind=0.81, γ_MES=0.79, Δ=0.02\n");
    printf("  Kagome 4-state (doublet-augmented): γ_ind=%.3f, γ_MES=%.3f, Δ=%.3f\n",
           -beta_i, -beta_m, -beta_i-(-beta_m));

    free(psi);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
