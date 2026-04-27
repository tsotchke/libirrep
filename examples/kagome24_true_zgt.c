/* SPDX-License-Identifier: MIT */
/* True Zhang-Grover-Turner modular S-matrix on kagome 2×4 N=24.
 *
 * Fixes two critical errors in kagome24_modular_s.c:
 *
 *   1. MANIFOLD (fixed): uses {(kx=0,ky=0), (kx=0,ky=π), (kx=π,ky=0),
 *      (kx=π,ky=π)} — the 4 BZ corners where Z₂ topological sectors should
 *      live, each with F=+1 (genuine Sz=0 singlets). Previous code used
 *      (0,1)/(0,3) doublet pair with F≈−0.30, which are spin-wave excitations
 *      outside any Z₂ topological tower.
 *
 *   2. REGIONS (fixed): Region A = y=0 row {0..5}, 6 sites, wraps x-cycle.
 *      Region B = x=0 unit-cell column {0,1,2,6,7,8,12,13,14,18,19,20},
 *      12 sites, wraps y-cycle. These are topologically INEQUIVALENT cuts.
 *      Previous code used two rows both wrapping the same x-cycle.
 *
 * Computational techniques for the 12-site Region B:
 *   - OpenMP-parallelised cross_rdm (20 cores → ~10 min instead of 2.5 h)
 *   - G-tensor precomputation for Renyi-2 entropy: avoids repeated Jacobi
 *     eigendecomposition of 4096×4096 matrices during the MES grid scan.
 *     Region B entropies reported as Renyi-2 (S₂); Region A as von Neumann (S₁).
 *
 * Expected result if Z₂ spin liquid: all |⟨MES_A|MES_B⟩| entries ≈ 0.5
 * (uniform 4×4 matrix = Z₂ character table magnitudes).
 *
 * Kitaev calibration (§1.24): Kitaev A-phase with same wrapping protocol
 * gives mean|S|=0.470 (Z₂ prediction 0.5, 6% finite-size error).
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_true_zgt 2>&1 | tee /tmp/kagome24_true_zgt.log
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

#define N_SITES 24
#define D_FULL  (1LL << N_SITES)
#define POPCOUNT 12
#define NS 4

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* ---- unfold symmetry sector → dense 2^N vector ---- */
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
            if (f >= 0) oa[f] += w[g]; else { op[ne] = gx; oa[ne] = w[g]; ++ne; }
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

/* ---- build F=+1 ground state at k=(kx,ky) ---- */
static int build_singlet_at_k(const irrep_heisenberg_t *H,
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
        seed[i] = 0.13 * sin(0.41 * (double)i) + I * 0.07 * cos(0.23 * (double)i);
    double *eigs = malloc((size_t)K_eigs * sizeof(double));
    double _Complex *psi_all = malloc((size_t)K_eigs * sdim * sizeof(double _Complex));
    int it = (K_eigs * 10 < (int)sdim) ? K_eigs * 10 : (int)sdim;
    if (it < 200) it = (200 < (int)sdim) ? 200 : (int)sdim;
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
        /* No F=+1 state found; return lowest eigenstate with its F value. */
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

/* ---- cross-reduced-density-matrix ρ^{ij}_A[a,b] = Σ_e ψ_i*[a,e] ψ_j[b,e] ---- */
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
#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic, 64)
#endif
    for (long long a = 0; a < dA; ++a) {
        for (long long b = 0; b < dA; ++b) {
            double _Complex acc = 0;
            const long long *sa = s_from_ae + a * dE;
            const long long *sb = s_from_ae + b * dE;
            for (long long e = 0; e < dE; ++e)
                acc += conj(psi_i[sa[e]]) * psi_j[sb[e]];
            rho_out[a * dA + b] = acc;
        }
    }
    free(s_from_ae);
}

/* ---- assemble ρ_A(c) and compute S₁ (von Neumann) ---- */
static double entropy_s1(double _Complex *const *r, int N, int nA,
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

/* ---- Renyi-2 G-tensor approach for large regions ----
 * Precompute G[a,b,c,d] = Tr(ρ^{ab} × ρ^{cd}) = Σ_{ij} ρ^{ab}[i,j] ρ^{cd}[j,i]
 * Then S₂(c) = -log Re(Σ_{abcd} c_a* c_b c_c* c_d G[a,b,c,d])
 * This avoids eigendecomposition during the grid scan.
 */
static void precompute_g_tensor(double _Complex *const *r, int N, long long dA,
                                double _Complex *G) {
    /* Precompute transposed cross-RDMs to make inner product cache-friendly. */
    double _Complex **rT = malloc((size_t)N * N * sizeof(double _Complex *));
    for (int a = 0; a < N; ++a)
    for (int b = 0; b < N; ++b) {
        rT[a * N + b] = malloc((size_t)dA * dA * sizeof(double _Complex));
        const double _Complex *src = r[a * N + b];
        double _Complex *dst = rT[a * N + b];
        for (long long i = 0; i < dA; ++i)
        for (long long j = 0; j < dA; ++j)
            dst[i * dA + j] = src[j * dA + i];  /* transpose */
    }

    /* G[abcd] = ⟨ρ^{ab}, ρ^{cd}_T⟩_F = dot(vec(ρ^{ab}), vec(ρ^{cd}_T)) */
    int N4 = N * N * N * N;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int idx = 0; idx < N4; ++idx) {
        int a = idx / (N * N * N);
        int b = (idx / (N * N)) % N;
        int c = (idx / N) % N;
        int d = idx % N;
        const double _Complex *rab  = r[a * N + b];
        const double _Complex *rcdT = rT[c * N + d];
        double _Complex g = 0;
        long long sz = dA * dA;
        for (long long k = 0; k < sz; ++k)
            g += rab[k] * rcdT[k];
        G[idx] = g;
    }

    for (int k = 0; k < N * N; ++k) free(rT[k]);
    free(rT);
}

static double entropy_s2_fast(const double _Complex *G, int N,
                               const double _Complex *c) {
    double _Complex tr2 = 0;
    for (int a = 0; a < N; ++a)
    for (int b = 0; b < N; ++b)
    for (int c2 = 0; c2 < N; ++c2)
    for (int d = 0; d < N; ++d) {
        double _Complex w = conj(c[a]) * c[b] * conj(c[c2]) * c[d];
        tr2 += w * G[(size_t)((a * N * N * N) + (b * N * N) + (c2 * N) + d)];
    }
    double t = creal(tr2);
    return (t > 1e-300) ? -log(t) : 100.0;
}

/* ---- orthogonal complement projection ---- */
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

/* ---- find one MES in the complement of subspace_v, using the provided
 *     entropy function.  entropy_fn receives either s1_fn or s2_fn dispatch.
 * ---- */
typedef double (*entropy_fn_t)(void *ctx, int N, int nA, const double _Complex *c);

/* S₁ wrapper */
typedef struct { double _Complex *const *r; } s1_ctx_t;
static double s1_dispatch(void *ctx, int N, int nA, const double _Complex *c) {
    s1_ctx_t *s = (s1_ctx_t *)ctx;
    return entropy_s1(s->r, N, nA, c);
}

/* S₂ wrapper */
typedef struct { const double _Complex *G; } s2_ctx_t;
static double s2_dispatch(void *ctx, int N, int nA, const double _Complex *c) {
    s2_ctx_t *s = (s2_ctx_t *)ctx;
    (void)nA;
    return entropy_s2_fast(s->G, N, c);
}

static double find_mes_in_complement(
        entropy_fn_t efn, void *ectx, int N, int nA,
        const double _Complex *subspace_v, int nv,
        double _Complex *c_out) {
    int Nth, Nph;
    if (nA <= 6) { Nth = 5; Nph = 5; }
    else if (nA <= 8) { Nth = 4; Nph = 4; }
    else { Nth = 3; Nph = 3; }   /* denser grid for S₂ (cheap) */
    double S_min = +1e300;
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
                            double nr = 0;
                            for (int a = 0; a < N; ++a) nr += creal(c[a] * conj(c[a]));
                            if (nr < 0.9) continue;
                            double S = efn(ectx, N, nA, c);
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

static void find_mes_basis(entropy_fn_t efn, void *ectx, int N, int nA,
                           double _Complex *mes_basis, double *S_arr) {
    for (int m = 0; m < N; ++m) {
        double _Complex c[NS];
        double S = find_mes_in_complement(efn, ectx, N, nA, mes_basis, m, c);
        double nrm = 0;
        for (int a = 0; a < N; ++a) nrm += creal(c[a] * conj(c[a]));
        nrm = sqrt(nrm);
        for (int a = 0; a < N; ++a) c[a] /= nrm;
        for (int a = 0; a < N; ++a) mes_basis[(size_t)m * N + a] = c[a];
        S_arr[m] = S;
    }
}

/* ======================================================= */
int main(void) {
    printf("=== Kagome 2×4 N=24: True ZGT modular S-matrix ===\n");
    printf("    Manifold: {(kx=0,ky=0), (kx=0,ky=π), (kx=π,ky=0), (kx=π,ky=π)}\n");
    printf("    Regions:  A=row(y=0, |A|=6, x-wrap)  "
           "B=col(x=0, |A|=12, y-wrap)\n\n");
    double t0 = now_sec();

#ifdef _OPENMP
    printf("  OpenMP: %d threads\n\n", omp_get_max_threads());
#endif

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);

    /* ---- 1. Build 4-state manifold at all 4 BZ corners ---- */
    printf("  Building 4-state Z₂ candidate manifold at all BZ corners:\n");
    double _Complex *psi = malloc((size_t)NS * D_FULL * sizeof(double _Complex));
    double Es[NS], Fs[NS];
    /* On 2×4 torus: kx ∈ {0,1}, ky ∈ {0,1,2,3}.
     * Z₂ corners: (0,0), (0,ky=2=π), (kx=1=π,0), (kx=1=π,ky=2=π). */
    int kx_tab[NS] = {0, 0, 1, 1};
    int ky_tab[NS] = {0, 2, 0, 2};
    const char *k_label[NS] = {"(kx=0,ky=0)", "(kx=0,ky=π)", "(kx=π,ky=0)", "(kx=π,ky=π)"};
    for (int s = 0; s < NS; ++s) {
        double ts = now_sec();
        int found = build_singlet_at_k(H, G, T, kx_tab[s], ky_tab[s], 8,
                                       psi + (size_t)s * D_FULL, &Es[s], &Fs[s]);
        printf("    %s: E=%+.6f  F=%+.4f  %s (%.1fs)\n",
               k_label[s], Es[s], Fs[s],
               found >= 0 ? "F=+1 found" : "NO F=+1 (fallback)", now_sec() - ts);
        fflush(stdout);
    }
    double E_min = Es[0], E_max = Es[0];
    for (int s = 1; s < NS; ++s) {
        if (Es[s] < E_min) E_min = Es[s];
        if (Es[s] > E_max) E_max = Es[s];
    }
    printf("\n  Tower span: Δ = %.4f J   (finite-size splitting)\n", E_max - E_min);
    int all_singlets = 1;
    for (int s = 0; s < NS; ++s) if (Fs[s] < 0.5) all_singlets = 0;
    printf("  All 4 states have F=+1: %s\n\n",
           all_singlets ? "YES  (genuine singlets — Z₂ candidate manifold)" :
                          "NO   (mixed F — manifold contains non-singlets)");
    fflush(stdout);

    /* ---- 2. Region A: y=0 row, 6 sites, wraps x-cycle ---- */
    /* On kagome 2×4 with row-major ordering: y=0 row = sites 0..5 */
    int region_A[6] = {0, 1, 2, 3, 4, 5};
    int nA_A = 6;
    long long dA_A = 1LL << nA_A;

    printf("━ Region A: y=0 row {0..5}  (|A|=6, wraps x-cycle)  S₁ entropy ━\n");
    fflush(stdout);
    double _Complex **rA = malloc((size_t)NS * NS * sizeof(double _Complex *));
    double tA = now_sec();
    for (int a = 0; a < NS; ++a)
    for (int b = 0; b < NS; ++b) {
        rA[a * NS + b] = malloc((size_t)dA_A * dA_A * sizeof(double _Complex));
        cross_rdm(psi + (size_t)a * D_FULL, psi + (size_t)b * D_FULL,
                  region_A, nA_A, rA[a * NS + b]);
    }
    printf("  cross-RDMs precomputed (%.1fs)\n", now_sec() - tA); fflush(stdout);

    s1_ctx_t ctxA; ctxA.r = (double _Complex *const *)rA;
    double _Complex mes_A[NS * NS];
    double S_mes_A[NS];
    double tmA = now_sec();
    find_mes_basis(s1_dispatch, &ctxA, NS, nA_A, mes_A, S_mes_A);
    printf("  MES_A S₁ entropies: %.4f  %.4f  %.4f  %.4f  "
           "(spread %.4f, %.1fs)\n",
           S_mes_A[0], S_mes_A[1], S_mes_A[2], S_mes_A[3],
           S_mes_A[3] - S_mes_A[0], now_sec() - tmA);
    fflush(stdout);
    printf("  MES_A basis (|c_α|):\n");
    for (int m = 0; m < NS; ++m) {
        printf("    ");
        for (int a = 0; a < NS; ++a) printf(" %.3f", cabs(mes_A[m * NS + a]));
        printf("\n");
    }

    /* ---- 3. Region B: x=0 column, 12 sites, wraps y-cycle ---- */
    /* Unit cells at col_x=0, all 4 rows: {0,1,2}, {6,7,8}, {12,13,14}, {18,19,20} */
    int region_B[12] = {0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20};
    int nA_B = 12;
    long long dA_B = 1LL << nA_B;

    printf("\n━ Region B: x=0 column {0,1,2,6,7,8,12,13,14,18,19,20}"
           "  (|A|=12, wraps y-cycle)  S₂ entropy ━\n");
    printf("  (cross-RDM computation uses OpenMP; MES uses Renyi-2 fast entropy)\n");
    fflush(stdout);
    double _Complex **rB = malloc((size_t)NS * NS * sizeof(double _Complex *));
    double tB = now_sec();
    for (int a = 0; a < NS; ++a)
    for (int b = 0; b < NS; ++b) {
        rB[a * NS + b] = malloc((size_t)dA_B * dA_B * sizeof(double _Complex));
        cross_rdm(psi + (size_t)a * D_FULL, psi + (size_t)b * D_FULL,
                  region_B, nA_B, rB[a * NS + b]);
    }
    printf("  cross-RDMs precomputed (%.1fs)\n", now_sec() - tB); fflush(stdout);

    /* Precompute G tensor for Renyi-2 fast MES. */
    double _Complex *G_tensor = malloc((size_t)NS * NS * NS * NS * sizeof(double _Complex));
    double tG = now_sec();
    precompute_g_tensor((double _Complex *const *)rB, NS, dA_B, G_tensor);
    printf("  G tensor precomputed (%.1fs)\n", now_sec() - tG); fflush(stdout);

    s2_ctx_t ctxB; ctxB.G = G_tensor;
    double _Complex mes_B[NS * NS];
    double S_mes_B[NS];
    double tmB = now_sec();
    find_mes_basis(s2_dispatch, &ctxB, NS, nA_B, mes_B, S_mes_B);
    printf("  MES_B S₂ entropies: %.4f  %.4f  %.4f  %.4f  "
           "(spread %.4f, %.1fs)\n",
           S_mes_B[0], S_mes_B[1], S_mes_B[2], S_mes_B[3],
           S_mes_B[3] - S_mes_B[0], now_sec() - tmB);
    fflush(stdout);
    printf("  MES_B basis (|c_α|):\n");
    for (int m = 0; m < NS; ++m) {
        printf("    ");
        for (int a = 0; a < NS; ++a) printf(" %.3f", cabs(mes_B[m * NS + a]));
        printf("\n");
    }
    free(G_tensor);

    /* ---- 4. Overlap matrix S_αβ = |⟨MES_A^α | MES_B^β⟩| ---- */
    printf("\n━ |⟨MES_A | MES_B⟩| overlap matrix (modular S candidates) ━\n");
    printf("  row=MES_A (x-wrap), col=MES_B (y-wrap):\n");
    double S_matrix[NS * NS];
    for (int a = 0; a < NS; ++a)
    for (int b = 0; b < NS; ++b) {
        double _Complex ov = 0;
        for (int k = 0; k < NS; ++k)
            ov += conj(mes_A[a * NS + k]) * mes_B[b * NS + k];
        S_matrix[a * NS + b] = cabs(ov);
    }
    for (int a = 0; a < NS; ++a) {
        printf("    ");
        for (int b = 0; b < NS; ++b) printf("%.4f  ", S_matrix[a * NS + b]);
        printf("\n");
    }

    double mean_s = 0, var_s = 0, rms_z2 = 0;
    for (int k = 0; k < NS * NS; ++k) mean_s += S_matrix[k];
    mean_s /= (NS * NS);
    for (int k = 0; k < NS * NS; ++k)
        var_s += (S_matrix[k] - mean_s) * (S_matrix[k] - mean_s);
    for (int k = 0; k < NS * NS; ++k) {
        double d = S_matrix[k] - 0.5;
        rms_z2 += d * d;
    }
    printf("\n  mean|S| = %.4f   std = %.4f   RMS from |S_Z₂|=0.5: %.4f\n",
           mean_s, sqrt(var_s / (NS * NS)), sqrt(rms_z2 / (NS * NS)));

    /* Detect block structure */
    double off01 = 0, off23 = 0;
    for (int a = 0; a < 2; ++a)
    for (int b = 2; b < 4; ++b) off01 += S_matrix[a * NS + b] + S_matrix[b * NS + a];
    off01 /= 8;
    for (int a = 0; a < 2; ++a)
    for (int b = 0; b < 2; ++b) off23 += S_matrix[a * NS + b];
    off23 /= 4;
    printf("  Off-diagonal blocks (0,1)×(2,3): mean = %.4f\n", off01);
    printf("  Upper-left 2×2 block (0,1)×(0,1): mean = %.4f\n", off23);

    printf("\n  CALIBRATION (Kitaev §1.24, same protocol, inequivalent regions):\n");
    printf("    mean|S|=0.470, RMS=0.173 for Z₂ A-phase\n");

    if (fabs(mean_s - 0.5) < 0.12 && sqrt(rms_z2 / (NS * NS)) < 0.25 && off01 > 0.3)
        printf("\n  → CONSISTENT with Z₂ topological order:\n"
               "    uniform overlap matrix, all 4 sectors coupled.\n");
    else if (off01 < 0.1)
        printf("\n  → BLOCK-DIAGONAL structure: singlet 2-block decoupled from\n"
               "    other 2-block. NOT consistent with 4-anyon Z₂ identification.\n");
    else
        printf("\n  → INTERMEDIATE: partial coupling. Finite-size effects or\n"
               "    manifold impurity (non-singlet states) distort signal.\n");

    /* ---- Cleanup ---- */
    for (int k = 0; k < NS * NS; ++k) { free(rA[k]); free(rB[k]); }
    free(rA); free(rB); free(psi);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
