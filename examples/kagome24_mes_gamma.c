/* SPDX-License-Identifier: MIT */
/* MES (Minimum Entanglement State) extraction of γ on the N=24
 * kagome 2×4 near-degenerate singlet manifold.
 *
 *   Zhang-Grover-Turner, Phys. Rev. B 85, 235151 (2012):
 *
 *   For a topologically-ordered GS with quantum dimension D, the
 *   torus hosts D² near-degenerate GS. For any bipartition A, the
 *   minimum-entropy superposition |Ψ_MES⟩ = c₀|ψ₀⟩ + c₁|ψ₁⟩ + …
 *   satisfies S_A(|Ψ_MES⟩) = α|∂A| − log D + O(exp(−L/ξ)).
 *
 * Two-state implementation: tower scan (§1.17b) found two
 * low-lying F=+1 singlets at Γ and (0, π) with ΔE = 0.054 J. We
 * minimise S_A over the projective two-state family
 *
 *     |Ψ(θ, φ)⟩ = cos θ |ψ_Γ⟩ + e^{iφ} sin θ |ψ_(0,π)⟩,
 *
 * precomputing four cross-RDMs
 *
 *     ρ_A^{ij}[a,b] = Σ_env ⟨a,env|ψ_i⟩⟨ψ_j|b,env⟩
 *
 * once, then assembling
 *
 *     ρ_A(θ,φ) = Σ_{ij} c_i c_j^* ρ_A^{ij}
 *
 * on a (θ, φ) grid and diagonalising.  Minimum S_A over (θ, φ)
 * is the two-state-manifold MES. If the FULL manifold has d > 2
 * degenerate states, this is a *lower bound* on the true γ_MES.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome24_mes_gamma
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

/* Build the unfolded singlet GS at a given (kx, ky). */
static void build_gs_at_k(const irrep_heisenberg_t *H,
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
            found = k;
            break;
        }
    }
    if (found < 0) { *E_out = eigs[0]; *F_out = 0; memcpy(psi_full_out, psi_scan, (size_t)D_FULL * sizeof(double _Complex)); }

    free(seed); free(eigs); free(psi_all); free(w); free(psi_scan);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
}

/* Build the cross reduced density matrix
 *   ρ_A^{ij}[a, b] = Σ_env ⟨a, env | ψ_i⟩ ⟨ψ_j | b, env⟩
 * by direct enumeration. Loops over all 2^N states, each contributes
 * to one (a, env) decomposition; accumulates into a 2^|A| × 2^|A|
 * matrix. */
static void cross_rdm(const double _Complex *psi_i,
                      const double _Complex *psi_j,
                      const int *A, int nA,
                      double _Complex *rho_out /* 2^nA × 2^nA */) {
    long long dA = 1LL << nA;
    long long dE = D_FULL >> nA;
    memset(rho_out, 0, (size_t)dA * dA * sizeof(double _Complex));

    /* Build bitmask for "is site in A" and map site → position within A/env. */
    int in_A[64] = {0}, pos_A[64] = {0}, pos_E[64] = {0};
    int ia = 0, ie = 0;
    for (int k = 0; k < nA; ++k) { in_A[A[k]] = 1; pos_A[A[k]] = ia++; }
    for (int s = 0; s < N_SITES; ++s) if (!in_A[s]) pos_E[s] = ie++;

    /* For each full-state index, extract (a, env) indices, accumulate. */
    /* To get ρ[a, b] = Σ_env ψ_i(a, env)^* ψ_j(b, env), we traverse
     * full states once per row. Simpler: enumerate (env, a) and (env, b)
     * together. O(D_FULL × 2^nA) work. For N=24, |A|=3: 16M × 8 = 128M ops. */
    for (long long s = 0; s < D_FULL; ++s) {
        long long a = 0, e = 0;
        for (int k = 0; k < N_SITES; ++k) {
            int bit = (int)((s >> k) & 1);
            if (in_A[k]) a |= ((long long)bit) << pos_A[k];
            else          e |= ((long long)bit) << pos_E[k];
        }
        /* s is the (a, e) combined; we need ψ_i(b, e)^* ψ_j(a, e)
         * accumulated into ρ[a, b]. Simplest: for each env, hold all
         * 2^|A| amplitudes in a strip. But that requires an inner loop
         * over "sister" full states with same e but different a. */
        /* Approach: accumulate row-vectors. v_i[e, a] = ψ_i(a, e).
         * Then ρ[a, b] = Σ_e v_j[e, a] v_i[e, b]^*.
         * Too much memory: dA × dE complex = 2^N complex total. OK,
         * same as psi_full. Actually we can just refer to psi_i and
         * psi_j through their (a, e) decomposition, which is above. */
        (void)a; (void)e;
    }
    /* Direct implementation: for each pair of basis strips (a, b)
     * with shared env mask, sum over env. This is O(dA² × dE). For
     * nA=3: 8 × 8 × 2^21 = 134M ops; for nA=9: 512 × 512 × 2^15 = 8.6G;
     * for nA=12: 4096² × 4096 = 64G. Fine for nA ≤ 9. */
    for (long long e = 0; e < dE; ++e) {
        /* Reconstruct the full-state indices for (a, e) and (b, e). */
        /* Fastest: build a lookup table mapping each full s to (a, e). */
    }
    /* Build the explicit lookup. */
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
    for (long long a = 0; a < dA; ++a) {
        for (long long b = 0; b < dA; ++b) {
            double _Complex acc = 0;
            for (long long e = 0; e < dE; ++e) {
                long long s_ae = s_from_ae[a * dE + e];
                long long s_be = s_from_ae[b * dE + e];
                acc += conj(psi_i[s_ae]) * psi_j[s_be];
            }
            rho_out[a * dA + b] = acc;
        }
    }
    free(s_from_ae);
}

/* For a 2-state manifold, assemble ρ_A(θ, φ) from the four cross-RDMs
 * and return S_A. */
static double entropy_2state(const double _Complex *r00,
                             const double _Complex *r01,
                             const double _Complex *r10,
                             const double _Complex *r11,
                             int nA, double theta, double phi) {
    long long dA = 1LL << nA;
    double c0 = cos(theta), c1 = sin(theta);
    double _Complex p01 = c0 * c1 * cexp(I * phi);
    double _Complex p10 = c0 * c1 * cexp(-I * phi);
    double c00 = c0 * c0, c11 = c1 * c1;
    double _Complex *rho = malloc((size_t)dA * dA * sizeof(double _Complex));
    for (long long i = 0; i < dA * dA; ++i)
        rho[i] = c00 * r00[i] + c11 * r11[i] + p01 * r10[i] + p10 * r01[i];
    double S = irrep_entropy_vonneumann(rho, (int)dA);
    free(rho);
    return S;
}

int main(void) {
    printf("=== N=24 kagome 2×4 MES-γ extraction on 2-state manifold ===\n");
    printf("    Zhang-Grover-Turner 2012 construction\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);

    double _Complex *psi0 = malloc((size_t)D_FULL * sizeof(double _Complex));
    double _Complex *psi1 = malloc((size_t)D_FULL * sizeof(double _Complex));

    printf("  Building GS at (0,0)…\n"); fflush(stdout);
    double E0, F0;
    double tG = now_sec();
    build_gs_at_k(H, G, T, 0, 0, 6, psi0, &E0, &F0);
    printf("    E=%+.8f, F=%+.4f (%.1fs)\n", E0, F0, now_sec() - tG); fflush(stdout);

    printf("  Building GS at (0,π)…\n"); fflush(stdout);
    double E1, F1;
    double tG2 = now_sec();
    build_gs_at_k(H, G, T, 0, 2, 6, psi1, &E1, &F1);
    printf("    E=%+.8f, F=%+.4f (%.1fs)\n", E1, F1, now_sec() - tG2); fflush(stdout);

    /* Verify orthogonality (states at different k are orthogonal on
     * a translation-invariant Hamiltonian, to machine precision). */
    double _Complex ov = 0;
    for (long long s = 0; s < D_FULL; ++s) ov += conj(psi0[s]) * psi1[s];
    printf("  ⟨ψ_Γ|ψ_(0,π)⟩ = %+.2e %+.2ei  (should be ~0)\n",
           creal(ov), cimag(ov));

    /* MES scan on (θ, φ) grid. Three region choices at |A|=3 and
     * |A|=4 each (fast). */
    int region_sizes[2] = {3, 4};
    for (int rs_i = 0; rs_i < 2; ++rs_i) {
        int nA = region_sizes[rs_i];
        int fam[3][8] = {
            {0,1,2,3,4,5,6,7},
            {0,2,4,6,8,10,12,14},
            {0,3,6,9,12,15,18,21},
        };
        const char *fam_name[3] = {"contiguous", "stride-2 ", "sublatt  "};

        printf("\n━ |A|=%d, MES scan on (θ, φ) ━\n", nA);
        for (int f = 0; f < 3; ++f) {
            double tf = now_sec();
            long long dA = 1LL << nA;
            size_t rdm_size = dA * dA;
            double _Complex *r00 = malloc(rdm_size * sizeof(double _Complex));
            double _Complex *r01 = malloc(rdm_size * sizeof(double _Complex));
            double _Complex *r10 = malloc(rdm_size * sizeof(double _Complex));
            double _Complex *r11 = malloc(rdm_size * sizeof(double _Complex));
            cross_rdm(psi0, psi0, fam[f], nA, r00);
            cross_rdm(psi0, psi1, fam[f], nA, r01);
            cross_rdm(psi1, psi0, fam[f], nA, r10);
            cross_rdm(psi1, psi1, fam[f], nA, r11);

            /* (θ, φ) grid: 32 × 32. Minimum and maximum S_A. */
            double S_min = +INFINITY, S_max = -INFINITY;
            double th_min = 0, ph_min = 0;
            int Nth = 32, Nph = 32;
            for (int ith = 0; ith <= Nth; ++ith) {
                double th = ith * 0.5 * M_PI / Nth;
                for (int iph = 0; iph < Nph; ++iph) {
                    double ph = iph * 2.0 * M_PI / Nph;
                    double S = entropy_2state(r00, r01, r10, r11, nA, th, ph);
                    if (S < S_min) { S_min = S; th_min = th; ph_min = ph; }
                    if (S > S_max) S_max = S;
                }
            }
            /* Pure ψ_0 and ψ_1 entropies as references. */
            double S_psi0 = entropy_2state(r00, r01, r10, r11, nA, 0.0, 0.0);
            double S_psi1 = entropy_2state(r00, r01, r10, r11, nA, M_PI/2, 0.0);
            printf("  family %s  S(ψ_Γ)=%.4f  S(ψ_(0,π))=%.4f  MES: S_min=%.4f  S_max=%.4f  (θ*=%.2f φ*=%.2f) (%.1fs)\n",
                   fam_name[f], S_psi0, S_psi1, S_min, S_max,
                   th_min, ph_min, now_sec() - tf);
            fflush(stdout);

            free(r00); free(r01); free(r10); free(r11);
        }
    }

    free(psi0); free(psi1);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
