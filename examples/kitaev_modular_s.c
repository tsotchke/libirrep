/* SPDX-License-Identifier: MIT */
/* Modular S-matrix extraction on the Kitaev honeycomb A-phase.
 *
 * Zhang-Grover-Turner (PRB 85, 235151, 2012) prescription:
 *
 *   Given D² near-degenerate GSes spanning the topological manifold,
 *   the Minimum-Entanglement-State (MES) basis for a given region
 *   shape is the basis of 4 orthogonal states that each locally
 *   minimise S_A(c). Two inequivalent region shapes give two MES
 *   bases; the 4×4 overlap matrix between them is the modular
 *   S-matrix (up to a global phase / permutation of anyon labels).
 *
 * For Z₂ toric code the expected S matrix is
 *
 *           1  [ 1  1  1  1 ]
 *   S  =  --- [ 1  1 -1 -1 ]
 *           2  [ 1 -1  1 -1 ]
 *              [ 1 -1 -1  1 ]
 *
 * which is the character table of Z₂ × Z₂.
 *
 * Procedure:
 *   1. Lanczos 4 tower states (120 iters).
 *   2. For region_1 (horizontal strip): precompute 4×4 cross-RDMs.
 *   3. Find 4 orthogonal CP³ MES states via multi-start gradient
 *      descent (50 random starts; cluster by orthogonality).
 *   4. Repeat for region_2 (vertical strip).
 *   5. Compute 4×4 overlap matrix between MES bases.
 *   6. Compare to S_toric.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kitaev_modular_s
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
#define NS 4

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

/* Assemble ρ_A(c) and compute S_A. */
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

/* Project c onto the orthogonal complement of subspace_v (nv vectors). */
static void project_orthogonal(double _Complex *c, int N,
                               const double _Complex *subspace_v, int nv) {
    for (int k = 0; k < nv; ++k) {
        const double _Complex *v = subspace_v + (size_t)k * N;
        double _Complex ov = 0;
        for (int a = 0; a < N; ++a) ov += conj(v[a]) * c[a];
        for (int a = 0; a < N; ++a) c[a] -= ov * v[a];
    }
    /* Renormalise. */
    double nrm = 0;
    for (int a = 0; a < N; ++a) nrm += creal(c[a] * conj(c[a]));
    nrm = sqrt(nrm);
    if (nrm > 1e-12) {
        for (int a = 0; a < N; ++a) c[a] /= nrm;
    }
}

/* Find MES in orthogonal complement of subspace_v via grid scan.
 * Grid density auto-scales with dA cost so big |A| stays tractable. */
static double find_mes_in_complement(
        double _Complex *const *r, int N, int nA,
        const double _Complex *subspace_v, int nv,
        double _Complex *c_out) {
    int Nth, Nph;
    if (nA <= 6) { Nth = 5; Nph = 5; }        /* 27k grid */
    else if (nA == 7) { Nth = 3; Nph = 3; }   /* 1728 grid */
    else { Nth = 2; Nph = 2; }                 /* 216 grid — coarse for |A|≥8 */
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

/* Find 4 orthogonal MES states. */
static void find_mes_basis(double _Complex *const *r, int N, int nA,
                           double _Complex *mes_basis /* NS × NS, row-major */,
                           double *S_mes_arr) {
    for (int m = 0; m < N; ++m) {
        double _Complex c[NS];
        double S = find_mes_in_complement(r, N, nA, mes_basis, m, c);
        /* Renormalise. */
        double nrm = 0;
        for (int a = 0; a < N; ++a) nrm += creal(c[a] * conj(c[a]));
        nrm = sqrt(nrm);
        for (int a = 0; a < N; ++a) c[a] /= nrm;
        for (int a = 0; a < N; ++a) mes_basis[(size_t)m * N + a] = c[a];
        S_mes_arr[m] = S;
    }
}

int main(void) {
    printf("=== Kitaev honeycomb A-phase modular S-matrix extraction ===\n");
    printf("    3×4 torus, K_z=1, K_x=K_y=0.1\n\n");
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

    double eigs[NS+1];
    double _Complex *psi = malloc((size_t)(NS+1) * D_FULL * sizeof(double _Complex));
    double tL = now_sec();
    printf("  Lanczos (120 iters)…\n"); fflush(stdout);
    irrep_lanczos_eigvecs_reorth(kitaev_apply, &K, D_FULL, NS+1, 120,
                                 seed, eigs, psi);
    printf("    done (%.1fs). E:", now_sec() - tL);
    for (int k = 0; k < NS+1; ++k) printf(" %+.6f", eigs[k]);
    printf("\n    tower span E_%d − E_0 = %+.4e   gap E_%d − E_%d = %+.4e\n",
           NS-1, eigs[NS-1] - eigs[0], NS, NS-1, eigs[NS] - eigs[NS-1]);
    fflush(stdout);
    free(seed);

    /* TOPOLOGICALLY INEQUIVALENT wrapping regions, both cutting the
     * torus with non-contractible boundaries:
     *   Region A wraps x-cycle: entire y=0 row, 6 sites
     *   Region B wraps y-cycle: entire x=0 column, 8 sites
     * The modular S matrix appears as the unitary change-of-basis
     * between the two resulting MES bases. A local-region pair
     * (both non-wrapping) would give identity overlap — which is
     * what kitaev_modular_s originally found with the wrong regions. */
    int region_A[6] = {0, 1, 2, 3, 4, 5};         /* y=0 row, wraps x */
    int region_B[8] = {0, 1, 6, 7, 12, 13, 18, 19}; /* x=0 column, wraps y */
    int nA_A = 6, nA_B = 8;
    long long dA_A = 1LL << nA_A;
    long long dA_B = 1LL << nA_B;
    long long dA = dA_A;  /* alias for legacy code below */
    int nA = nA_A;
    (void)dA;  (void)nA;

    /* Region A: y=0 row (wraps x-cycle), 6 sites. */
    printf("\n━ Region A: y=0 row [0,1,2,3,4,5] (wraps x-cycle) ━\n");
    fflush(stdout);
    double _Complex **rA = malloc((size_t)NS*NS * sizeof(double _Complex *));
    double tf = now_sec();
    for (int a = 0; a < NS; ++a)
        for (int b = 0; b < NS; ++b) {
            rA[a * NS + b] = malloc((size_t)dA_A * dA_A * sizeof(double _Complex));
            cross_rdm(psi + (size_t)a * D_FULL, psi + (size_t)b * D_FULL,
                      region_A, nA_A, rA[a * NS + b]);
        }
    printf("  cross-RDMs precomputed (%.1fs)\n", now_sec() - tf);
    fflush(stdout);

    double _Complex mes_A[NS*NS];
    double S_mes_A[NS];
    double tmesA = now_sec();
    find_mes_basis(rA, NS, nA_A, mes_A, S_mes_A);
    printf("  4 orthogonal MES entropies: %.4f  %.4f  %.4f  %.4f  (%.1fs)\n",
           S_mes_A[0], S_mes_A[1], S_mes_A[2], S_mes_A[3], now_sec() - tmesA);
    fflush(stdout);
    /* Print MES basis (4×4 matrix). */
    printf("  MES_A basis (rows = MES states, columns = Lanczos states):\n");
    for (int m = 0; m < NS; ++m) {
        printf("    ");
        for (int a = 0; a < NS; ++a) {
            double _Complex c = mes_A[m * NS + a];
            printf(" %+.3f%+.3fi  ", creal(c), cimag(c));
        }
        printf("  (|c|² = ");
        double nrm = 0;
        for (int a = 0; a < NS; ++a) nrm += creal(mes_A[m*NS+a]*conj(mes_A[m*NS+a]));
        printf("%.3f)\n", nrm);
    }

    /* Region B: x=0 column (wraps y-cycle), 8 sites. */
    printf("\n━ Region B: x=0 column [0,1,6,7,12,13,18,19] (wraps y-cycle) ━\n");
    fflush(stdout);
    double _Complex **rB = malloc((size_t)NS*NS * sizeof(double _Complex *));
    tf = now_sec();
    for (int a = 0; a < NS; ++a)
        for (int b = 0; b < NS; ++b) {
            rB[a * NS + b] = malloc((size_t)dA_B * dA_B * sizeof(double _Complex));
            cross_rdm(psi + (size_t)a * D_FULL, psi + (size_t)b * D_FULL,
                      region_B, nA_B, rB[a * NS + b]);
        }
    printf("  cross-RDMs precomputed (%.1fs)\n", now_sec() - tf);
    fflush(stdout);

    double _Complex mes_B[NS*NS];
    double S_mes_B[NS];
    double tmesB = now_sec();
    find_mes_basis(rB, NS, nA_B, mes_B, S_mes_B);
    printf("  4 orthogonal MES entropies: %.4f  %.4f  %.4f  %.4f  (%.1fs)\n",
           S_mes_B[0], S_mes_B[1], S_mes_B[2], S_mes_B[3], now_sec() - tmesB);
    fflush(stdout);
    printf("  MES_B basis:\n");
    for (int m = 0; m < NS; ++m) {
        printf("    ");
        for (int a = 0; a < NS; ++a) {
            double _Complex c = mes_B[m * NS + a];
            printf(" %+.3f%+.3fi  ", creal(c), cimag(c));
        }
        printf("\n");
    }

    /* Overlap matrix S_αβ = ⟨MES_A^α | MES_B^β⟩ = Σ_k conj(mes_A[α,k]) · mes_B[β,k]. */
    printf("\n━ Modular S-matrix (overlaps |⟨MES_A | MES_B⟩|) ━\n");
    double S_matrix[NS*NS];
    for (int a = 0; a < NS; ++a)
        for (int b = 0; b < NS; ++b) {
            double _Complex ov = 0;
            for (int k = 0; k < NS; ++k)
                ov += conj(mes_A[a * NS + k]) * mes_B[b * NS + k];
            S_matrix[a * NS + b] = cabs(ov);
        }
    printf("  row=MES_A, col=MES_B:\n");
    for (int a = 0; a < NS; ++a) {
        printf("    ");
        for (int b = 0; b < NS; ++b) printf("%.4f  ", S_matrix[a * NS + b]);
        printf("\n");
    }

    /* Compare to S_toric = (1/2) [[1,1,1,1], [1,1,-1,-1], [1,-1,1,-1], [1,-1,-1,1]]. */
    const double S_toric[NS*NS] = {
        0.5, 0.5, 0.5, 0.5,
        0.5, 0.5, 0.5, 0.5,
        0.5, 0.5, 0.5, 0.5,
        0.5, 0.5, 0.5, 0.5,
    };
    printf("\n  Expected |S_toric| = 0.5 for every entry (Z₂ character table magnitudes).\n");
    double mean = 0, var = 0;
    for (int k = 0; k < NS*NS; ++k) mean += S_matrix[k];
    mean /= (NS*NS);
    for (int k = 0; k < NS*NS; ++k) var += (S_matrix[k]-mean)*(S_matrix[k]-mean);
    double rms_from_ideal = 0;
    for (int k = 0; k < NS*NS; ++k) {
        double d = S_matrix[k] - S_toric[k];
        rms_from_ideal += d * d;
    }
    rms_from_ideal = sqrt(rms_from_ideal / (NS*NS));
    printf("  Observed mean |S| = %.4f,  std = %.4f\n", mean, sqrt(var/(NS*NS)));
    printf("  RMS deviation from S_toric: %.4f\n", rms_from_ideal);
    if (rms_from_ideal < 0.15 && fabs(mean - 0.5) < 0.1)
        printf("  ✓ Structure consistent with Z₂ toric-code modular S-matrix.\n");
    else
        printf("  ~ Structure does not cleanly match S_toric; refinement needed.\n");

    /* Cleanup. */
    for (int k = 0; k < NS*NS; ++k) { free(rA[k]); free(rB[k]); }
    free(rA); free(rB);
    free(psi); free(bi); free(bj); free(bt);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
