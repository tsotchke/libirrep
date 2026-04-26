/* SPDX-License-Identifier: MIT */
/* Kitaev honeycomb A-phase MES γ scaling: S_MES vs |∂A| regression.
 *
 * The preceding Kitaev MES test at |A|=4 found Δ = S_min − mean(S_i)
 * of 0.07 nats (stride-2 family), about 10× below log 2. That size
 * is dominated by subleading finite-size corrections. The cleanest
 * γ extraction in the MES literature (Zhang-Grover-Turner 2012)
 * fits S_MES = α|∂A| − γ + O(exp(−L/ξ)) across multiple |A|.
 *
 * This script:
 *   1. Lanczos to 120 iters → try to catch the 4th Z₂ topological
 *      sector state (we previously found 3 / 4 expected states).
 *   2. For |A| ∈ {2, 3, 4, 5, 6}, compute 9 cross-RDMs and the MES
 *      on a downsampled CP³ grid.
 *   3. Fit S_MES = α|∂A| − γ_MES.
 *   4. Independently fit S_individual = α|∂A| − γ_ind.
 *   5. Difference γ_ind − γ_MES = how much the individual-state
 *      γ over-estimates the true topological invariant.
 *
 * Positive-control target: γ_MES should be within ~log 2 of the
 * actual topological γ = log 2 (exact for Kitaev A-phase toric code).
 * Deviation quantifies the methodology's finite-size resolution
 * ceiling at N = 24.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kitaev_mes_scaling
 */

#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define LX 3
#define LY 4
#define N_SITES (2 * LX * LY)
#define D_FULL  (1LL << N_SITES)
#define N_STATES 4

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

typedef struct {
    int     nb;
    int    *bi, *bj, *bt;
    double  Kx, Ky, Kz;
} kitaev_ctx_t;

static void kitaev_apply(const double _Complex *psi,
                         double _Complex *out, void *opaque) {
    kitaev_ctx_t *K = (kitaev_ctx_t *)opaque;
    memset(out, 0, (size_t)D_FULL * sizeof(double _Complex));
    for (int b = 0; b < K->nb; ++b) {
        int i = K->bi[b], j = K->bj[b], t = K->bt[b];
        long long mi = 1LL << i, mj = 1LL << j;
        double coup = (t == 0) ? K->Kx : (t == 1) ? K->Ky : K->Kz;
        for (long long s = 0; s < D_FULL; ++s) {
            int zi = (int)((s >> i) & 1), zj = (int)((s >> j) & 1);
            if (t == 2) {
                double sgn = (zi == zj) ? +1.0 : -1.0;
                out[s] += -coup * sgn * psi[s];
            } else if (t == 0) {
                long long sp = s ^ mi ^ mj;
                out[sp] += -coup * psi[s];
            } else {
                long long sp = s ^ mi ^ mj;
                double sgn = (zi == zj) ? -1.0 : +1.0;
                out[sp] += -coup * sgn * psi[s];
            }
        }
    }
}

static int classify_bond(double dx, double dy) {
    if (fabs(dy) < 0.1) return 0;
    if (dy > 0)         return 1;
    return 2;
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

/* Assemble ρ_A(c) from NxN cross-RDM tensor and coefficient vector c. */
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
    printf("=== Kitaev A-phase N=24: MES γ from S_MES vs |∂A| scaling ===\n\n");
    double t0 = now_sec();

    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_HONEYCOMB, LX, LY);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    int *bt = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    double a1[2], a2[2];
    irrep_lattice_primitive_vectors(L, a1, a2);
    for (int b = 0; b < nb; ++b) {
        double pi_[2], pj_[2];
        irrep_lattice_site_position(L, bi[b], pi_);
        irrep_lattice_site_position(L, bj[b], pj_);
        double dx0 = pj_[0] - pi_[0], dy0 = pj_[1] - pi_[1];
        double best_r = 1e300, bdx = 0, bdy = 0;
        for (int nx = -1; nx <= 1; ++nx)
            for (int ny = -1; ny <= 1; ++ny) {
                double dx = dx0 + nx * LX * a1[0] + ny * LY * a2[0];
                double dy = dy0 + nx * LX * a1[1] + ny * LY * a2[1];
                double r = sqrt(dx*dx + dy*dy);
                if (r < best_r) { best_r = r; bdx = dx; bdy = dy; }
            }
        bt[b] = classify_bond(bdx, bdy);
    }
    kitaev_ctx_t K = { nb, bi, bj, bt, 0.1, 0.1, 1.0 };

    double _Complex *seed = malloc((size_t)D_FULL * sizeof(double _Complex));
    for (long long s = 0; s < D_FULL; ++s)
        seed[s] = 0.1 * sin(0.37 * s) + I * 0.05 * cos(0.23 * s);
    double sn = 0;
    for (long long s = 0; s < D_FULL; ++s) sn += creal(seed[s] * conj(seed[s]));
    sn = sqrt(sn);
    for (long long s = 0; s < D_FULL; ++s) seed[s] /= sn;

    double eigs[N_STATES + 2];
    double _Complex *psi = malloc((size_t)(N_STATES + 2) * D_FULL * sizeof(double _Complex));
    double tL = now_sec();
    printf("  Lanczos (120 iters, dim=%lld)…\n", D_FULL); fflush(stdout);
    irrep_lanczos_eigvecs_reorth(kitaev_apply, &K, D_FULL, N_STATES + 2, 120,
                                 seed, eigs, psi);
    printf("    done (%.1fs). E:", now_sec() - tL);
    for (int k = 0; k < N_STATES + 2; ++k) printf(" %+.6f", eigs[k]);
    printf("\n    tower span E_{N-1}-E_0=%+.4e, gap E_N-E_{N-1}=%+.4e\n",
           eigs[N_STATES-1] - eigs[0], eigs[N_STATES] - eigs[N_STATES-1]);
    fflush(stdout);
    free(seed);

    /* Region series: nested contiguous [0, |A|) for |A| in {2, 3, 4, 5}.
     * Could extend to |A|=6 but 64×64 RDM + grid = longer. */
    int As[4] = {2, 3, 4, 5};
    int nsizes = 4;
    int *region[4];
    int bdry[4];
    double S_mes_arr[4], S_ind_arr[4];

    for (int si = 0; si < nsizes; ++si) {
        int nA = As[si];
        region[si] = malloc(sizeof(int) * nA);
        for (int i = 0; i < nA; ++i) region[si][i] = i;
        bdry[si] = boundary_bonds(region[si], nA, bi, bj, nb);

        printf("\n━ |A|=%d, |∂A|=%d ━\n", nA, bdry[si]); fflush(stdout);
        long long dA = 1LL << nA;

        /* 16 cross-RDMs for 4×4 case. */
        int Nst = N_STATES;
        double _Complex **r = malloc((size_t)Nst * Nst * sizeof(double _Complex *));
        double tf = now_sec();
        for (int a = 0; a < Nst; ++a)
            for (int b = 0; b < Nst; ++b) {
                r[a * Nst + b] = malloc((size_t)dA * dA * sizeof(double _Complex));
                cross_rdm(psi + (size_t)a * D_FULL,
                          psi + (size_t)b * D_FULL, region[si], nA, r[a * Nst + b]);
            }
        printf("  %d×%d cross-RDMs (%.1fs)\n", Nst, Nst, now_sec() - tf);
        fflush(stdout);

        /* Individual entropies. */
        double S_i[N_STATES];
        for (int a = 0; a < Nst; ++a) {
            double _Complex c[N_STATES] = {0};
            c[a] = 1;
            S_i[a] = entropy_ns(r, Nst, nA, c);
        }
        double S_ind_mean = 0;
        for (int a = 0; a < Nst; ++a) S_ind_mean += S_i[a];
        S_ind_mean /= Nst;
        S_ind_arr[si] = S_ind_mean;
        printf("  S_i:");
        for (int a = 0; a < Nst; ++a) printf(" %.4f", S_i[a]);
        printf("   mean %.4f\n", S_ind_mean);

        /* MES on CP^{N-1}: for N=4, coords are (θ1, θ2, θ3, φ1, φ2, φ3).
         * Parameterise as:
         *   c_0 = cos θ1
         *   c_1 = sin θ1 cos θ2        e^{iφ1}
         *   c_2 = sin θ1 sin θ2 cos θ3 e^{iφ2}
         *   c_3 = sin θ1 sin θ2 sin θ3 e^{iφ3}
         * θ_k ∈ [0, π/2], φ_k ∈ [0, 2π). A coarse 6⁶ = 46656 grid. */
        double tmes = now_sec();
        int Nth = 6, Nph = 6;
        double S_min = +INFINITY, S_max = -INFINITY;
        double c_min[N_STATES] = {0};
        long long grid_cnt = 0;
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
                                double _Complex c[N_STATES];
                                c[0] = cos(th1);
                                c[1] = sin(th1) * cos(th2) * cexp(I * ph1);
                                c[2] = sin(th1) * sin(th2) * cos(th3) * cexp(I * ph2);
                                c[3] = sin(th1) * sin(th2) * sin(th3) * cexp(I * ph3);
                                double S = entropy_ns(r, Nst, nA, c);
                                if (S < S_min) {
                                    S_min = S;
                                    for (int k = 0; k < Nst; ++k)
                                        c_min[k] = cabs(c[k]);
                                }
                                if (S > S_max) S_max = S;
                                ++grid_cnt;
                            }
                        }
                    }
                }
            }
        }
        S_mes_arr[si] = S_min;
        printf("  MES (grid=%lld): S_min=%.4f  S_max=%.4f  range=%.4f  (%.1fs)\n",
               grid_cnt, S_min, S_max, S_max - S_min, now_sec() - tmes);
        printf("  S_min − S_ind_mean = %.4f   (log 2 = %.4f)\n",
               S_min - S_ind_mean, log(2.0));
        printf("  c_MES ≈ (");
        for (int k = 0; k < Nst; ++k)
            printf("%+.3f%s", c_min[k], k < Nst-1 ? " " : "");
        printf(")\n");
        fflush(stdout);

        for (int a = 0; a < Nst * Nst; ++a) free(r[a]);
        free(r);
    }

    /* Linear fit S = α|∂A| − γ on both series. */
    printf("\n━ Linear fit S vs |∂A| ━\n");
    printf("  |A|   |∂A|   S_individual_mean   S_MES\n");
    for (int si = 0; si < nsizes; ++si)
        printf("  %-3d   %-4d   %+.4f              %+.4f\n",
               As[si], bdry[si], S_ind_arr[si], S_mes_arr[si]);

    double mx=0, myi=0, mym=0, mxx=0, mxyi=0, mxym=0;
    for (int si = 0; si < nsizes; ++si) {
        mx   += bdry[si];
        myi  += S_ind_arr[si];
        mym  += S_mes_arr[si];
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
    /* Compute R² for both. */
    double ss_res_i=0, ss_tot_i=0, ss_res_m=0, ss_tot_m=0;
    for (int si = 0; si < nsizes; ++si) {
        double p_i = alpha_i * bdry[si] + beta_i;
        double p_m = alpha_m * bdry[si] + beta_m;
        ss_res_i += (S_ind_arr[si] - p_i) * (S_ind_arr[si] - p_i);
        ss_tot_i += (S_ind_arr[si] - myi) * (S_ind_arr[si] - myi);
        ss_res_m += (S_mes_arr[si] - p_m) * (S_mes_arr[si] - p_m);
        ss_tot_m += (S_mes_arr[si] - mym) * (S_mes_arr[si] - mym);
    }
    double R2_i = 1 - ss_res_i/ss_tot_i;
    double R2_m = 1 - ss_res_m/ss_tot_m;
    printf("\n  Individual-state fit:\n");
    printf("    S_ind = %+.4f |∂A| + %+.4f   γ_ind = %+.4f   (R²=%.4f)\n",
           alpha_i, beta_i, -beta_i, R2_i);
    printf("  MES fit:\n");
    printf("    S_MES = %+.4f |∂A| + %+.4f   γ_MES = %+.4f   (R²=%.4f)\n",
           alpha_m, beta_m, -beta_m, R2_m);
    printf("\n  γ_ind − γ_MES = %+.4f     (Kitaev prediction: γ = log 2 = +%.4f)\n",
           -beta_i - (-beta_m), log(2.0));
    printf("  γ_MES vs log 2: deviation = %+.4f  (%.1f%% of log 2)\n",
           -beta_m - log(2.0), 100.0 * (-beta_m - log(2.0)) / log(2.0));

    for (int si = 0; si < nsizes; ++si) free(region[si]);
    free(psi); free(bi); free(bj); free(bt);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
