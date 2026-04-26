/* SPDX-License-Identifier: MIT */
/* Modular S-matrix extraction on the kagome 2×4 N=24 4-state
 * doublet-augmented candidate manifold:
 *
 *   |ψ_0⟩ = singlet GS at k = (0, 0)              F = +1
 *   |ψ_1⟩ = singlet GS at k = (0, π)              F = +1
 *   |ψ_2⟩, |ψ_3⟩ = sym/anti combinations of (0, 1) + (0, 3) doublet
 *
 * The Kitaev calibration (§1.24) showed that on a known-Z₂
 * topological state, wrapping-region MES produces an overlap matrix
 * with mean magnitude 0.47 ± 0.17 (predicted 0.5 for Z₂), while
 * local regions produce a near-identity matrix.
 *
 * If kagome is Z₂ topological with this 4-state manifold, the same
 * structural shift should appear here. If not, the wrapping-region
 * matrix should look like the local-region one (near-identity), and
 * Z₂ identification fails.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_modular_s
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

static void project_orthogonal(double _Complex *c, int N,
                               const double _Complex *subspace_v, int nv) {
    for (int k = 0; k < nv; ++k) {
        const double _Complex *v = subspace_v + (size_t)k * N;
        double _Complex ov = 0;
        for (int a = 0; a < N; ++a) ov += conj(v[a]) * c[a];
        for (int a = 0; a < N; ++a) c[a] -= ov * v[a];
    }
    double nrm = 0;
    for (int a = 0; a < N; ++a) nrm += creal(c[a] * conj(c[a]));
    nrm = sqrt(nrm);
    if (nrm > 1e-12) for (int a = 0; a < N; ++a) c[a] /= nrm;
}

static double find_mes_in_complement(
        double _Complex *const *r, int N, int nA,
        const double _Complex *subspace_v, int nv,
        double _Complex *c_out) {
    int Nth = 5, Nph = 5;
    double S_min = +INFINITY;
    double _Complex c_best[NS] = {0};
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
                            if (nv > 0) project_orthogonal(c, N, subspace_v, nv);
                            double nrm = 0;
                            for (int a = 0; a < N; ++a) nrm += creal(c[a]*conj(c[a]));
                            if (nrm < 0.9) continue;
                            double S = entropy_ns(r, N, nA, c);
                            if (S < S_min) {
                                S_min = S;
                                for (int k = 0; k < N; ++k) c_best[k] = c[k];
                            }
                        }
                    }
                }
            }
        }
    }
    for (int k = 0; k < N; ++k) c_out[k] = c_best[k];
    return S_min;
}

static void find_mes_basis(double _Complex *const *r, int N, int nA,
                           double _Complex *mes_basis,
                           double *S_mes_arr) {
    for (int m = 0; m < N; ++m) {
        double _Complex c[NS];
        double S = find_mes_in_complement(r, N, nA, mes_basis, m, c);
        double nrm = 0;
        for (int a = 0; a < N; ++a) nrm += creal(c[a] * conj(c[a]));
        nrm = sqrt(nrm);
        for (int a = 0; a < N; ++a) c[a] /= nrm;
        for (int a = 0; a < N; ++a) mes_basis[(size_t)m * N + a] = c[a];
        S_mes_arr[m] = S;
    }
}

int main(void) {
    printf("=== Kagome 2×4 N=24 modular S-matrix extraction ===\n");
    printf("    4-state manifold: 2 singlets {(0,0), (0,π)} + (0,1)/(0,3) doublet\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);

    /* Build 4-state manifold. */
    double _Complex *psi = malloc((size_t)NS * D_FULL * sizeof(double _Complex));
    double _Complex *psi_01 = malloc((size_t)D_FULL * sizeof(double _Complex));
    double _Complex *psi_03 = malloc((size_t)D_FULL * sizeof(double _Complex));
    double Es[NS], Fs[NS];

    printf("  Building |ψ_0⟩ singlet at (0,0)…\n"); fflush(stdout);
    double tG = now_sec();
    build_gs_at_k_singlet(H, G, T, 0, 0, 6, psi + 0*D_FULL, &Es[0], &Fs[0]);
    printf("    E=%+.6f F=%+.4f (%.1fs)\n", Es[0], Fs[0], now_sec() - tG); fflush(stdout);

    printf("  Building |ψ_1⟩ singlet at (0,π)…\n"); fflush(stdout);
    tG = now_sec();
    build_gs_at_k_singlet(H, G, T, 0, 2, 6, psi + 1*D_FULL, &Es[1], &Fs[1]);
    printf("    E=%+.6f F=%+.4f (%.1fs)\n", Es[1], Fs[1], now_sec() - tG); fflush(stdout);

    printf("  Building lowest at (0,1)…\n"); fflush(stdout);
    tG = now_sec();
    double E_01, F_01;
    build_gs_at_k_lowest(H, G, T, 0, 1, psi_01, &E_01, &F_01);
    printf("    E=%+.6f F=%+.4f (%.1fs)\n", E_01, F_01, now_sec() - tG); fflush(stdout);

    printf("  Building lowest at (0,3)…\n"); fflush(stdout);
    tG = now_sec();
    double E_03, F_03;
    build_gs_at_k_lowest(H, G, T, 0, 3, psi_03, &E_03, &F_03);
    printf("    E=%+.6f F=%+.4f (%.1fs)\n", E_03, F_03, now_sec() - tG); fflush(stdout);

    /* Sym/anti combinations. */
    double _Complex *psi_sym = psi + 2*D_FULL;
    double _Complex *psi_ant = psi + 3*D_FULL;
    double sqrt2 = sqrt(2.0);
    for (long long s = 0; s < D_FULL; ++s) {
        psi_sym[s] = (psi_01[s] + psi_03[s]) / sqrt2;
        psi_ant[s] = (psi_01[s] - psi_03[s]) / sqrt2;
    }
    double n_sym = 0, n_ant = 0;
    for (long long s = 0; s < D_FULL; ++s) {
        n_sym += creal(psi_sym[s] * conj(psi_sym[s]));
        n_ant += creal(psi_ant[s] * conj(psi_ant[s]));
    }
    double sf_sym = 1.0/sqrt(n_sym);
    double sf_ant = 1.0/sqrt(n_ant);
    for (long long s = 0; s < D_FULL; ++s) {
        psi_sym[s] *= sf_sym;
        psi_ant[s] *= sf_ant;
    }
    free(psi_01); free(psi_03);
    Es[2] = E_01; Es[3] = E_03;

    printf("\n  Tower: E=[%+.4f, %+.4f, %+.4f, %+.4f]  span=%.4f\n",
           Es[0], Es[1], Es[2], Es[3],
           fmax(fmax(Es[0],Es[1]),fmax(Es[2],Es[3]))
           - fmin(fmin(Es[0],Es[1]),fmin(Es[2],Es[3])));
    fflush(stdout);

    /* Region A: y=0 row, sites {0,1,2,3,4,5} (wraps y-cycle on kagome 2×4). */
    int region_A[6] = {0, 1, 2, 3, 4, 5};
    /* Region B: local 6-site block elsewhere, sites {6,7,8,9,10,11}. */
    int region_B[6] = {6, 7, 8, 9, 10, 11};
    int nA = 6;
    long long dA = 1LL << nA;

    /* Region A: wrapping. */
    printf("\n━ Region A: y=0 row [0,1,2,3,4,5] (wraps y-cycle) ━\n");
    fflush(stdout);
    double _Complex **rA = malloc((size_t)NS*NS * sizeof(double _Complex *));
    double tf = now_sec();
    for (int a = 0; a < NS; ++a)
        for (int b = 0; b < NS; ++b) {
            rA[a*NS+b] = malloc((size_t)dA*dA * sizeof(double _Complex));
            cross_rdm(psi + (size_t)a*D_FULL, psi + (size_t)b*D_FULL,
                      region_A, nA, rA[a*NS+b]);
        }
    printf("  cross-RDMs precomputed (%.1fs)\n", now_sec() - tf);
    fflush(stdout);

    double _Complex mes_A[NS*NS];
    double S_mes_A[NS];
    double tmesA = now_sec();
    find_mes_basis(rA, NS, nA, mes_A, S_mes_A);
    printf("  MES_A entropies: %.4f  %.4f  %.4f  %.4f  (spread %.4f, %.1fs)\n",
           S_mes_A[0], S_mes_A[1], S_mes_A[2], S_mes_A[3],
           S_mes_A[3]-S_mes_A[0], now_sec() - tmesA);
    fflush(stdout);
    printf("  MES_A basis (|c_α|):\n");
    for (int m = 0; m < NS; ++m) {
        printf("    ");
        for (int a = 0; a < NS; ++a) printf(" %.3f", cabs(mes_A[m*NS+a]));
        printf("\n");
    }

    /* Region B: local. */
    printf("\n━ Region B: y=1 row [6,7,8,9,10,11] (also wraps y-cycle) ━\n");
    fflush(stdout);
    double _Complex **rB = malloc((size_t)NS*NS * sizeof(double _Complex *));
    tf = now_sec();
    for (int a = 0; a < NS; ++a)
        for (int b = 0; b < NS; ++b) {
            rB[a*NS+b] = malloc((size_t)dA*dA * sizeof(double _Complex));
            cross_rdm(psi + (size_t)a*D_FULL, psi + (size_t)b*D_FULL,
                      region_B, nA, rB[a*NS+b]);
        }
    printf("  cross-RDMs precomputed (%.1fs)\n", now_sec() - tf);
    fflush(stdout);

    double _Complex mes_B[NS*NS];
    double S_mes_B[NS];
    double tmesB = now_sec();
    find_mes_basis(rB, NS, nA, mes_B, S_mes_B);
    printf("  MES_B entropies: %.4f  %.4f  %.4f  %.4f  (spread %.4f, %.1fs)\n",
           S_mes_B[0], S_mes_B[1], S_mes_B[2], S_mes_B[3],
           S_mes_B[3]-S_mes_B[0], now_sec() - tmesB);
    fflush(stdout);

    /* Overlap matrix. */
    printf("\n━ |⟨MES_A | MES_B⟩| matrix ━\n");
    double S_matrix[NS*NS];
    for (int a = 0; a < NS; ++a)
        for (int b = 0; b < NS; ++b) {
            double _Complex ov = 0;
            for (int k = 0; k < NS; ++k)
                ov += conj(mes_A[a*NS+k]) * mes_B[b*NS+k];
            S_matrix[a*NS+b] = cabs(ov);
        }
    for (int a = 0; a < NS; ++a) {
        printf("    ");
        for (int b = 0; b < NS; ++b) printf("%.4f  ", S_matrix[a*NS+b]);
        printf("\n");
    }

    double mean = 0, var = 0;
    for (int k = 0; k < NS*NS; ++k) mean += S_matrix[k];
    mean /= NS*NS;
    for (int k = 0; k < NS*NS; ++k) var += (S_matrix[k]-mean)*(S_matrix[k]-mean);
    double rms_from_z2 = 0;
    for (int k = 0; k < NS*NS; ++k) {
        double d = S_matrix[k] - 0.5;
        rms_from_z2 += d*d;
    }
    rms_from_z2 = sqrt(rms_from_z2/(NS*NS));
    printf("\n  mean|S| = %.4f   std = %.4f   RMS from |S_Z₂| = 0.5: %.4f\n",
           mean, sqrt(var/(NS*NS)), rms_from_z2);
    printf("\n  CALIBRATION (Kitaev N=24 §1.24, both regions wrapping):\n");
    printf("    mean|S|=0.470, std=0.171, RMS from S_toric = 0.173\n");
    printf("    Kitaev local pair: mean|S|=0.291, near-identity\n");

    if (fabs(mean - 0.5) < 0.1 && rms_from_z2 < 0.25)
        printf("\n  → Kagome modular S structure CONSISTENT with Z₂ — overlap matrix\n"
               "    delocalised with mean ≈ 0.5, similar to Kitaev. Strong Z₂ signal.\n");
    else if (mean < 0.4 && var < 0.05*0.05)
        printf("\n  → Kagome modular S CLOSE TO IDENTITY — wrapping regions don't\n"
               "    resolve topological structure. Z₂ identification disfavoured.\n");
    else
        printf("\n  → Kagome modular S in intermediate regime; needs comparison\n"
               "    against control of similar manifold truncation.\n");

    for (int k = 0; k < NS*NS; ++k) { free(rA[k]); free(rB[k]); }
    free(rA); free(rB);
    free(psi);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
