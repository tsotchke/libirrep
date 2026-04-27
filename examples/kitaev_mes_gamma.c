/* SPDX-License-Identifier: MIT */
/* MES γ extraction on the Kitaev honeycomb A-phase 3-state near-
 * degenerate tower (N = 24, 3×4 torus, K_z=1, K_x=K_y=0.1).
 *
 * The audit of kagome MES (§…) showed that a 2-state manifold from
 * (Γ, (0,π)) spans at most 1 topological sector and MES reduces S_A
 * by only <0.03 nats. On exactly-solvable Kitaev A-phase we know
 * γ_exact = log 2 and the near-degenerate manifold we found has
 * 3 states (within 3.6e-4 of each other) — truncated from the
 * expected 4 but likely spans multiple Z₂ sectors.
 *
 * This script:
 *   1. Builds the 3 near-degenerate Kitaev GSes via Lanczos.
 *   2. Precomputes the 3×3 cross-RDM tensor on 3 region families.
 *   3. Scans the CP² manifold of normalised 3-coefficient
 *      superpositions on a (θ₁, θ₂, φ₁, φ₂) grid.
 *   4. Reports the minimum S_A; checks whether it approaches
 *      S_A − log 2 relative to the individual-state S_A.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kitaev_mes_gamma
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define LX 3
#define LY 4
#define N_SITES (2 * LX * LY)
#define D_FULL  (1LL << N_SITES)

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
    if (fabs(dy) < 0.1)          return 0;
    if (dy > 0)                  return 1;
    return 2;
}

/* Cross RDM ρ_A^{ij}[a, b] = Σ_env ⟨a, env|ψ_i⟩* ⟨ψ_j|b, env⟩*
 * As in kagome24_mes_gamma.c, enumerate (a, e) decomposition. */
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

/* Assemble ρ_A(c) from 3×3 cross-RDM tensor and coefficient vector c.
 * c is a length-3 complex vector with |c|²=1. */
static double entropy_3state(double _Complex *const r[3][3], int nA,
                             const double _Complex *c) {
    long long dA = 1LL << nA;
    double _Complex *rho = malloc((size_t)dA * dA * sizeof(double _Complex));
    memset(rho, 0, (size_t)dA * dA * sizeof(double _Complex));
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            double _Complex w = conj(c[a]) * c[b];
            for (long long k = 0; k < dA * dA; ++k)
                rho[k] += w * r[a][b][k];
        }
    double S = irrep_entropy_vonneumann(rho, (int)dA);
    free(rho);
    return S;
}

int main(void) {
    printf("=== Kitaev honeycomb A-phase 3-state MES γ ===\n");
    printf("    3×4 torus, K_z=1, K_x=K_y=0.1, tower span 3.6e-4 J\n\n");
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

    double eigs[4];
    double _Complex *psi = malloc((size_t)4 * D_FULL * sizeof(double _Complex));
    double tL = now_sec();
    printf("  Lanczos (80 iters, dim=%lld)…\n", D_FULL); fflush(stdout);
    irrep_lanczos_eigvecs_reorth(kitaev_apply, &K, D_FULL, 4, 80, seed, eigs, psi);
    printf("    done (%.1fs). E: %+.8f  %+.8f  %+.8f  %+.8f\n",
           now_sec() - tL, eigs[0], eigs[1], eigs[2], eigs[3]);
    printf("    tower span E_2 - E_0 = %+.4e,  gap E_3 - E_2 = %+.4e\n",
           eigs[2] - eigs[0], eigs[3] - eigs[2]);
    fflush(stdout);
    free(seed);

    /* Three region families, |A|=4. */
    int fam[3][8] = {
        {0,1,2,3},
        {0,2,4,6},
        {0,3,6,9},
    };
    const char *fam_name[3] = {"contiguous", "stride-2  ", "sublatt   "};
    int nA = 4;
    long long dA = 1LL << nA;

    for (int f = 0; f < 3; ++f) {
        double tf = now_sec();
        printf("\n━ |A|=%d, family %s ━\n", nA, fam_name[f]);

        /* Precompute 9 cross-RDMs. */
        double _Complex *r[3][3];
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                r[a][b] = malloc((size_t)dA * dA * sizeof(double _Complex));
                cross_rdm(psi + (size_t)a * D_FULL,
                          psi + (size_t)b * D_FULL, fam[f], nA, r[a][b]);
            }
        printf("  cross-RDMs precomputed (%.1fs)\n", now_sec() - tf);
        fflush(stdout);

        /* Individual-state entropies. */
        double S_i[3];
        for (int a = 0; a < 3; ++a) {
            double _Complex c[3] = {0, 0, 0};
            c[a] = 1;
            S_i[a] = entropy_3state(r, nA, c);
        }
        printf("  S(ψ_0)=%.4f  S(ψ_1)=%.4f  S(ψ_2)=%.4f\n",
               S_i[0], S_i[1], S_i[2]);

        /* CP² parameterization: c = (cos θ₁, sin θ₁ cos θ₂ e^{iφ₁}, sin θ₁ sin θ₂ e^{iφ₂}).
         * θ₁ ∈ [0, π/2], θ₂ ∈ [0, π/2], φ₁, φ₂ ∈ [0, 2π). Scan 16⁴ = 65 K points. */
        int Nth = 16, Nph = 16;
        double S_min = +INFINITY, S_max = -INFINITY;
        double th1_min=0, th2_min=0, ph1_min=0, ph2_min=0;
        double c_min[3] = {0, 0, 0};
        long long tot_grid = (long long)(Nth + 1) * (Nth + 1) * Nph * Nph;
        long long cnt = 0;
        for (int it1 = 0; it1 <= Nth; ++it1) {
            double th1 = it1 * 0.5 * M_PI / Nth;
            for (int it2 = 0; it2 <= Nth; ++it2) {
                double th2 = it2 * 0.5 * M_PI / Nth;
                for (int ip1 = 0; ip1 < Nph; ++ip1) {
                    double ph1 = ip1 * 2.0 * M_PI / Nph;
                    for (int ip2 = 0; ip2 < Nph; ++ip2) {
                        double ph2 = ip2 * 2.0 * M_PI / Nph;
                        double _Complex c[3];
                        c[0] = cos(th1);
                        c[1] = sin(th1) * cos(th2) * cexp(I * ph1);
                        c[2] = sin(th1) * sin(th2) * cexp(I * ph2);
                        double S = entropy_3state(r, nA, c);
                        if (S < S_min) {
                            S_min = S;
                            th1_min = th1; th2_min = th2;
                            ph1_min = ph1; ph2_min = ph2;
                            c_min[0] = creal(c[0]);
                            c_min[1] = cabs(c[1]);
                            c_min[2] = cabs(c[2]);
                        }
                        if (S > S_max) S_max = S;
                        ++cnt;
                    }
                }
            }
        }
        double S_mean_i = (S_i[0] + S_i[1] + S_i[2]) / 3.0;
        printf("  MES: S_min=%.4f  S_max=%.4f  range=%.4f  (%lld grid pts, %.1fs)\n",
               S_min, S_max, S_max - S_min, tot_grid, now_sec() - tf);
        printf("  S_min − mean(S_i) = %.4f  (vs −log 2 = %.4f)\n",
               S_min - S_mean_i, -log(2.0));
        printf("  c_MES ≈ (%+.3f, %+.3f, %+.3f)  (|c|₁ amplitudes)\n",
               c_min[0], c_min[1], c_min[2]);
        fflush(stdout);

        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) free(r[a][b]);
    }

    free(psi); free(bi); free(bj); free(bt);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
