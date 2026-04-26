/* SPDX-License-Identifier: MIT */
/* Li-Haldane entanglement spectrum of the Kitaev honeycomb A-phase.
 *
 * Calibration reference for kagome24_es.c.  Computes the reduced-density-
 * matrix eigenvalue spectrum ξ_i = −log λ_i for the A-phase ground state,
 * resolved by Sz_A sector.  The spectrum is the central Li-Haldane (2008)
 * diagnostic: a Z₂ topological phase shows an entanglement gap separating
 * a small set of "topological" low-ξ levels from a bulk of high-ξ levels,
 * with counting that matches the edge-mode theory.
 *
 * Geometry: Kitaev honeycomb 3×4 (N=24 sites), K_z=1, K_x=K_y=0.1 (A-phase).
 * Region A: y=0 row, sites {0..5}, nA=6, wraps x-cycle.
 * Region B: x=0 column, sites {0,1,6,7,12,13,18,19}, nA=8, wraps y-cycle.
 *
 *   make examples
 *   ./build/bin/kitaev_es
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

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

typedef struct { int nb, *bi, *bj, *bt; double Kx, Ky, Kz; } kitaev_ctx_t;

static void kitaev_apply(const double _Complex *in, double _Complex *out, void *op) {
    kitaev_ctx_t *K = op;
    memset(out, 0, (size_t)D_FULL * sizeof *out);
    for (int b = 0; b < K->nb; ++b) {
        int i = K->bi[b], j = K->bj[b], t = K->bt[b];
        long long mi = 1LL << i, mj = 1LL << j;
        double coup = t == 0 ? K->Kx : t == 1 ? K->Ky : K->Kz;
        for (long long s = 0; s < D_FULL; ++s) {
            int zi = (s >> i) & 1, zj = (s >> j) & 1;
            if (t == 2) {
                out[s] += -coup * ((zi == zj) ? 1.0 : -1.0) * in[s];
            } else if (t == 0) {
                out[s ^ mi ^ mj] += -coup * in[s];
            } else {
                out[s ^ mi ^ mj] += -coup * ((zi == zj) ? -1.0 : 1.0) * in[s];
            }
        }
    }
}

static int classify_bond(double dx, double dy) {
    if (fabs(dy) < 0.1) return 0;
    return dy > 0 ? 1 : 2;
}

static void cross_rdm(const double _Complex *pi, const double _Complex *pj,
                      const int *A, int nA, double _Complex *rho) {
    long long dA = 1LL << nA, dE = D_FULL >> nA;
    memset(rho, 0, (size_t)dA * dA * sizeof *rho);
    int in_A[64] = {0}, pos_A[64] = {0}, pos_E[64] = {0};
    int ia = 0, ie = 0;
    for (int k = 0; k < nA; ++k) { in_A[A[k]] = 1; pos_A[A[k]] = ia++; }
    for (int s = 0; s < N_SITES; ++s) if (!in_A[s]) pos_E[s] = ie++;
    long long *sf = malloc((size_t)dA * dE * sizeof(long long));
    for (long long s = 0; s < D_FULL; ++s) {
        long long a = 0, e = 0;
        for (int k = 0; k < N_SITES; ++k) {
            int bit = (s >> k) & 1;
            if (in_A[k]) a |= (long long)bit << pos_A[k];
            else          e |= (long long)bit << pos_E[k];
        }
        sf[a * dE + e] = s;
    }
    for (long long a = 0; a < dA; ++a)
        for (long long b = 0; b < dA; ++b) {
            double _Complex acc = 0;
            for (long long e = 0; e < dE; ++e)
                acc += conj(pi[sf[a*dE+e]]) * pj[sf[b*dE+e]];
            rho[a*dA+b] = acc;
        }
    free(sf);
}

/* Block-diagonalise ρ_A (dA×dA, row-major) by Sz_A = popcount(a) − nA/2.
 * Returns 2^nA entries sorted by ξ = −log λ ascending. */
static int entanglement_spectrum(const double _Complex *rho, int nA,
                                  int *Sz_A_out, double *xi_out) {
    long long dA = 1LL << nA;
    /* Max block size = C(nA, nA/2) — allocate conservatively. */
    double _Complex *sub = malloc((size_t)dA * dA * sizeof *sub);
    int *idxs = malloc((size_t)dA * sizeof *idxs);
    double *ev = malloc((size_t)dA * sizeof *ev);
    int n_total = 0;

    for (int nu = 0; nu <= nA; ++nu) {
        int bsz = 0;
        for (long long a = 0; a < dA; ++a)
            if (__builtin_popcountll((unsigned long long)a) == nu)
                idxs[bsz++] = (int)a;
        for (int i = 0; i < bsz; ++i)
            for (int j = 0; j < bsz; ++j)
                sub[i*bsz+j] = rho[(long long)idxs[i]*dA + idxs[j]];
        irrep_hermitian_eigvals(bsz, sub, ev);
        int Sz = nu - nA / 2;
        for (int k = 0; k < bsz; ++k) {
            double lam = ev[k] < 1e-300 ? 1e-300 : ev[k];
            xi_out[n_total]   = -log(lam);
            Sz_A_out[n_total] = Sz;
            ++n_total;
        }
    }
    free(sub); free(idxs); free(ev);

    /* Sort by ξ ascending (insertion sort; n ≤ 256 for nA ≤ 8). */
    for (int i = 1; i < n_total; ++i) {
        double xv = xi_out[i]; int Sv = Sz_A_out[i];
        int j = i - 1;
        while (j >= 0 && xi_out[j] > xv) {
            xi_out[j+1] = xi_out[j]; Sz_A_out[j+1] = Sz_A_out[j]; --j;
        }
        xi_out[j+1] = xv; Sz_A_out[j+1] = Sv;
    }
    return n_total;
}

static void print_es(const int *Sz, const double *xi, int n_total,
                     int n_print, int nA, const char *label) {
    int lim = n_print < n_total ? n_print : n_total;
    printf("  %s  (%d of %d eigenvalues)\n", label, lim, n_total);
    printf("  %5s  %6s  %12s\n", "rank", "Sz_A", "ξ = -log λ");
    for (int i = 0; i < lim; ++i)
        printf("  %5d  %+6d  %12.6f\n", i+1, Sz[i], xi[i]);

    /* Largest gap between consecutive ξ values in the printed window. */
    int gap_idx = 1; double gap_max = 0;
    for (int i = 1; i < lim; ++i) {
        double g = xi[i] - xi[i-1];
        if (g > gap_max) { gap_max = g; gap_idx = i; }
    }
    printf("  → Entanglement gap: Δξ = %.4f  (between rank %d and %d)\n",
           gap_max, gap_idx, gap_idx+1);

    /* Count levels per Sz_A below the gap. */
    int cnt[25] = {0};
    for (int i = 0; i < gap_idx; ++i) cnt[Sz[i] + 12]++;
    printf("    Levels below gap by Sz_A:");
    for (int s = -nA/2; s <= nA/2; ++s)
        if (cnt[s+12]) printf("  Sz=%+d:%d", s, cnt[s+12]);
    printf("\n");
}

int main(void) {
    printf("=== Kitaev honeycomb A-phase: Li-Haldane entanglement spectrum ===\n");
    printf("    3×4 torus, K_z=1, K_x=K_y=0.1  (N=24)\n\n");
    double t0 = now_sec();

    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_HONEYCOMB, LX, LY);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(nb * sizeof *bi);
    int *bj = malloc(nb * sizeof *bj);
    int *bt = malloc(nb * sizeof *bt);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    double a1[2], a2[2];
    irrep_lattice_primitive_vectors(L, a1, a2);
    for (int b = 0; b < nb; ++b) {
        double pi[2], pj[2];
        irrep_lattice_site_position(L, bi[b], pi);
        irrep_lattice_site_position(L, bj[b], pj);
        double dx0 = pj[0]-pi[0], dy0 = pj[1]-pi[1];
        double best_r = 1e300, bdx = 0, bdy = 0;
        for (int nx = -1; nx <= 1; ++nx)
            for (int ny = -1; ny <= 1; ++ny) {
                double dx = dx0 + nx*LX*a1[0] + ny*LY*a2[0];
                double dy = dy0 + nx*LX*a1[1] + ny*LY*a2[1];
                double r = sqrt(dx*dx+dy*dy);
                if (r < best_r) { best_r = r; bdx = dx; bdy = dy; }
            }
        bt[b] = classify_bond(bdx, bdy);
    }
    kitaev_ctx_t K = { nb, bi, bj, bt, 0.1, 0.1, 1.0 };

    double _Complex *seed = malloc((size_t)D_FULL * sizeof *seed);
    for (long long s = 0; s < D_FULL; ++s)
        seed[s] = 0.1*sin(0.37*s) + I*0.05*cos(0.23*s);
    double sn = 0;
    for (long long s = 0; s < D_FULL; ++s) sn += creal(seed[s]*conj(seed[s]));
    sn = sqrt(sn);
    for (long long s = 0; s < D_FULL; ++s) seed[s] /= sn;

    double eigs[2];
    double _Complex *psi = malloc((size_t)2 * D_FULL * sizeof *psi);
    printf("  Lanczos GS (120 iters)…"); fflush(stdout);
    double tL = now_sec();
    irrep_lanczos_eigvecs_reorth(kitaev_apply, &K, D_FULL, 2, 120, seed, eigs, psi);
    printf(" E_0 = %+.6f  (%.1fs)\n\n", eigs[0], now_sec()-tL);
    free(seed);

    /* ---- Region A: y=0 row {0..5}, nA=6, wraps x-cycle ---- */
    int region_A[] = {0,1,2,3,4,5};
    int nA_A = 6;
    long long dA_A = 1LL << nA_A;
    printf("━ Region A: y=0 row {0..5}  (nA=6, wraps x-cycle) ━\n");
    double _Complex *rhoA = malloc((size_t)dA_A * dA_A * sizeof *rhoA);
    double trd = now_sec();
    cross_rdm(psi, psi, region_A, nA_A, rhoA);
    printf("  RDM computed (%.1fs)\n", now_sec()-trd);
    double _Complex *rhoA_copy = malloc((size_t)dA_A * dA_A * sizeof *rhoA_copy);
    memcpy(rhoA_copy, rhoA, (size_t)dA_A * dA_A * sizeof *rhoA_copy);
    printf("  S₁ = %.6f\n", irrep_entropy_vonneumann(rhoA_copy, (int)dA_A));
    free(rhoA_copy);
    int   *Sz_A  = malloc((size_t)dA_A * sizeof *Sz_A);
    double *xi_A = malloc((size_t)dA_A * sizeof *xi_A);
    int nES_A = entanglement_spectrum(rhoA, nA_A, Sz_A, xi_A);
    print_es(Sz_A, xi_A, nES_A, 32, nA_A, "ES_A");
    free(rhoA); free(Sz_A); free(xi_A);

    /* ---- Region B: x=0 column {0,1,6,7,12,13,18,19}, nA=8, wraps y-cycle ---- */
    int region_B[] = {0,1,6,7,12,13,18,19};
    int nA_B = 8;
    long long dA_B = 1LL << nA_B;
    printf("\n━ Region B: x=0 column {0,1,6,7,12,13,18,19}  (nA=8, wraps y-cycle) ━\n");
    double _Complex *rhoB = malloc((size_t)dA_B * dA_B * sizeof *rhoB);
    trd = now_sec();
    cross_rdm(psi, psi, region_B, nA_B, rhoB);
    printf("  RDM computed (%.1fs)\n", now_sec()-trd);
    double _Complex *rhoB_copy = malloc((size_t)dA_B * dA_B * sizeof *rhoB_copy);
    memcpy(rhoB_copy, rhoB, (size_t)dA_B * dA_B * sizeof *rhoB_copy);
    printf("  S₁ = %.6f\n", irrep_entropy_vonneumann(rhoB_copy, (int)dA_B));
    free(rhoB_copy);
    int   *Sz_B  = malloc((size_t)dA_B * sizeof *Sz_B);
    double *xi_B = malloc((size_t)dA_B * sizeof *xi_B);
    int nES_B = entanglement_spectrum(rhoB, nA_B, Sz_B, xi_B);
    print_es(Sz_B, xi_B, nES_B, 48, nA_B, "ES_B");
    free(rhoB); free(Sz_B); free(xi_B);

    free(psi); free(bi); free(bj); free(bt);
    irrep_lattice_free(L);
    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec()-t0);
    return 0;
}
