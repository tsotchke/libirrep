/* SPDX-License-Identifier: MIT */
/* Li-Haldane entanglement spectrum of the kagome Heisenberg S=1/2 ground state.
 *
 * Applies the Li-Haldane (2008) diagnostic to the kagome 2×4 (N=24) system.
 * For a Z₂ spin liquid the entanglement Hamiltonian H_E = −log ρ_A of a
 * topological bipartition should exhibit an entanglement gap: a clear
 * separation between low-ξ "topological edge" levels and bulk high-ξ levels,
 * with counting matching the edge-mode theory.
 *
 * The two regions follow the ZGT protocol (§1.26):
 *   Region A: y=0 row {0..5},  nA=6, wraps x-cycle  (von Neumann S₁)
 *   Region B: x=0 column {0,1,2,6,7,8,12,13,14,18,19,20}, nA=12, wraps y-cycle (S₁)
 *
 * Reference: compare to kitaev_es.c (Z₂ A-phase positive control).
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_es 2>&1 | tee /tmp/kagome24_es.log
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

#ifdef _OPENMP
#include <omp.h>
#endif

#define N_SITES 24
#define D_FULL  (1LL << N_SITES)
#define POPCOUNT 12

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ---- unfold momentum sector → dense 2^N vector (from kagome24_true_zgt.c) ---- */
static void unfold(const irrep_space_group_t *G,
                   const irrep_sg_rep_table_t *T,
                   int order, const double _Complex *w,
                   const double _Complex *psi_sector, long long sector_dim,
                   double _Complex *psi_full) {
    memset(psi_full, 0, (size_t)D_FULL * sizeof *psi_full);
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
            if (f >= 0) oa[f] += w[g]; else { op[ne] = gx; oa[ne] = w[g]; ++ne; }
        }
        double nn = 0.0;
        for (int e = 0; e < ne; ++e) nn += creal(oa[e]*conj(oa[e]));
        if (nn < 1e-20) continue;
        double inv = 1.0 / sqrt(nn);
        double _Complex c = psi_sector[written];
        for (int e = 0; e < ne; ++e) psi_full[op[e]] += c * inv * oa[e];
        ++written;
    }
    free(op); free(oa);
}

/* ---- build F=+1 ground state at k=(kx,ky) ---- */
static int build_singlet_at_k(const irrep_heisenberg_t *H,
                               const irrep_space_group_t *G,
                               const irrep_sg_rep_table_t *T,
                               int kx, int ky,
                               double _Complex *psi_full_out,
                               double *E_out, double *F_out) {
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, kx, ky);
    double _Complex chi1[1] = {1.0};
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sdim = irrep_sg_heisenberg_sector_dim(S);

    double _Complex *seed = malloc((size_t)sdim * sizeof(double _Complex));
    for (long long i = 0; i < sdim; ++i)
        seed[i] = 0.13*sin(0.41*(double)i) + I*0.07*cos(0.23*(double)i);
    int K_eigs = 3;
    double *eigs = malloc((size_t)K_eigs * sizeof(double));
    double _Complex *psi_all = malloc((size_t)K_eigs * sdim * sizeof(double _Complex));
    int it = K_eigs * 10 < (int)sdim ? K_eigs * 10 : (int)sdim;
    if (it < 200) it = 200 < (int)sdim ? 200 : (int)sdim;
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sdim, K_eigs, it, seed, eigs, psi_all);

    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    double _Complex *psi_scan = malloc((size_t)D_FULL * sizeof(double _Complex));
    uint64_t mask = (1ULL << N_SITES) - 1;
    int found = -1;
    for (int k = 0; k < K_eigs; ++k) {
        unfold(G, T, order, w, psi_all + (size_t)k * sdim, sdim, psi_scan);
        double _Complex Fk = 0;
        for (long long s = 0; s < D_FULL; ++s)
            Fk += conj(psi_scan[s]) * psi_scan[s ^ mask];
        if (creal(Fk) > 0.5) {
            memcpy(psi_full_out, psi_scan, (size_t)D_FULL * sizeof(double _Complex));
            *E_out = eigs[k]; *F_out = creal(Fk);
            found = k; break;
        }
    }
    if (found < 0) {
        unfold(G, T, order, w, psi_all, sdim, psi_full_out);
        double _Complex Fk = 0;
        for (long long s = 0; s < D_FULL; ++s)
            Fk += conj(psi_full_out[s]) * psi_full_out[s ^ mask];
        *E_out = eigs[0]; *F_out = creal(Fk);
    }
    free(seed); free(eigs); free(psi_all); free(w); free(psi_scan);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
    return found;
}

/* ---- cross-RDM ρ^{ij}[a,b] = Σ_e ψ_i*[a,e] ψ_j[b,e] (OpenMP for large |A|) ---- */
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
#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic, 64)
#endif
    for (long long a = 0; a < dA; ++a)
        for (long long b = 0; b < dA; ++b) {
            double _Complex acc = 0;
            const long long *sa = sf + a*dE;
            const long long *sb = sf + b*dE;
            for (long long e = 0; e < dE; ++e)
                acc += conj(pi[sa[e]]) * pj[sb[e]];
            rho[a*dA+b] = acc;
        }
    free(sf);
}

/* ---- entanglement spectrum, block-diagonalised by Sz_A ---- */
static int entanglement_spectrum(const double _Complex *rho, int nA,
                                  int *Sz_A_out, double *xi_out) {
    long long dA = 1LL << nA;
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

    /* Sort by ξ ascending. */
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

/* von Neumann entropy from entanglement spectrum ξ_i = -log λ_i. */
static double entropy_from_es(const double *xi, int n) {
    double S = 0;
    for (int i = 0; i < n; ++i) {
        double lam = exp(-xi[i]);
        if (lam > 1e-300) S += lam * xi[i];  /* S₁ = Σ λ_i ξ_i = -Σ λ_i log λ_i */
    }
    return S;
}

static void print_es(const int *Sz, const double *xi, int n_total,
                     int n_print, int nA, const char *label) {
    int lim = n_print < n_total ? n_print : n_total;
    printf("  %s  (%d of %d eigenvalues)\n", label, lim, n_total);
    printf("  %5s  %6s  %12s\n", "rank", "Sz_A", "ξ = -log λ");
    for (int i = 0; i < lim; ++i)
        printf("  %5d  %+6d  %12.6f\n", i+1, Sz[i], xi[i]);

    int gap_idx = 1; double gap_max = 0;
    for (int i = 1; i < lim; ++i) {
        double g = xi[i] - xi[i-1];
        if (g > gap_max) { gap_max = g; gap_idx = i; }
    }
    printf("  → Entanglement gap: Δξ = %.4f  (between rank %d and %d)\n",
           gap_max, gap_idx, gap_idx+1);

    int cnt[25] = {0};
    for (int i = 0; i < gap_idx; ++i) cnt[Sz[i]+12]++;
    printf("    Levels below gap by Sz_A:");
    for (int s = -nA/2; s <= nA/2; ++s)
        if (cnt[s+12]) printf("  Sz=%+d:%d", s, cnt[s+12]);
    printf("\n");
}

int main(void) {
    printf("=== Kagome 2×4 N=24: Li-Haldane entanglement spectrum ===\n");
    printf("    GS at k=(kx=0, ky=π)  [true lowest state, E≈−10.760]\n");
#ifdef _OPENMP
    printf("    OpenMP: %d threads\n", omp_get_max_threads());
#endif
    printf("\n");
    double t0 = now_sec();

    /* Build kagome Heisenberg + P1 space group. */
    irrep_lattice_t *lat = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(lat, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(lat);
    int *bi = malloc(nb * sizeof *bi);
    int *bj = malloc(nb * sizeof *bj);
    irrep_lattice_fill_bonds_nn(lat, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);

    /* GS at k=(kx=0, ky=2) = (kx=0, ky=π) on 2×4 grid. */
    double _Complex *psi = malloc((size_t)D_FULL * sizeof(double _Complex));
    double E_gs, F_gs;
    double tgs = now_sec();
    printf("  Building GS at k=(kx=0, ky=π)…"); fflush(stdout);
    build_singlet_at_k(H, G, T, 0, 2, psi, &E_gs, &F_gs);
    printf(" E=%+.6f  F=%+.4f  (%.1fs)\n\n", E_gs, F_gs, now_sec()-tgs);

    /* ---- Region A: y=0 row {0..5}, nA=6, wraps x-cycle ---- */
    int region_A[] = {0,1,2,3,4,5};
    int nA_A = 6;
    long long dA_A = 1LL << nA_A;
    printf("━ Region A: y=0 row {0..5}  (nA=6, wraps x-cycle) ━\n");
    double _Complex *rhoA = malloc((size_t)dA_A * dA_A * sizeof *rhoA);
    double trd = now_sec();
    cross_rdm(psi, psi, region_A, nA_A, rhoA);
    printf("  RDM computed (%.1fs)\n", now_sec()-trd);
    int   *Sz_A  = malloc((size_t)dA_A * sizeof *Sz_A);
    double *xi_A = malloc((size_t)dA_A * sizeof *xi_A);
    int nES_A = entanglement_spectrum(rhoA, nA_A, Sz_A, xi_A);
    printf("  S₁ = %.6f  (from ES eigenvalues)\n", entropy_from_es(xi_A, nES_A));
    print_es(Sz_A, xi_A, nES_A, 32, nA_A, "ES_A");
    free(rhoA); free(Sz_A); free(xi_A);

    /* ---- Region B: x=0 column, nA=12, wraps y-cycle ---- */
    int region_B[] = {0,1,2,6,7,8,12,13,14,18,19,20};
    int nA_B = 12;
    long long dA_B = 1LL << nA_B;
    printf("\n━ Region B: x=0 column {0,1,2,6,7,8,12,13,14,18,19,20}  (nA=12, wraps y-cycle) ━\n");
    printf("  (cross-RDM uses OpenMP)\n");
    double _Complex *rhoB = malloc((size_t)dA_B * dA_B * sizeof *rhoB);
    trd = now_sec();
    cross_rdm(psi, psi, region_B, nA_B, rhoB);
    printf("  RDM computed (%.1fs)\n", now_sec()-trd);
    int   *Sz_B  = malloc((size_t)dA_B * sizeof *Sz_B);
    double *xi_B = malloc((size_t)dA_B * sizeof *xi_B);
    int nES_B = entanglement_spectrum(rhoB, nA_B, Sz_B, xi_B);
    printf("  S₁ = %.6f  (from ES eigenvalues)\n", entropy_from_es(xi_B, nES_B));
    print_es(Sz_B, xi_B, nES_B, 64, nA_B, "ES_B");
    free(rhoB); free(Sz_B); free(xi_B);

    free(psi); free(bi); free(bj);
    irrep_sg_rep_table_free(T);
    irrep_space_group_free(G);
    irrep_heisenberg_free(H);
    irrep_lattice_free(lat);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec()-t0);
    return 0;
}
