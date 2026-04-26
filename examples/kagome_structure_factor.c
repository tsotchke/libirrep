/* SPDX-License-Identifier: MIT */
/* Static spin structure factor S(k) = (1/N) Σ_{ij} e^{-ik·(r_i−r_j)} ⟨S_i·S_j⟩
 * on the N=24 kagome Heisenberg singlet ground state. Complements the
 * area-law γ extraction as an independent spin-liquid diagnostic:
 *
 *   Magnetically ordered (Néel etc.): S(k) has sharp Bragg peaks at the
 *     ordering momentum (K-point for 120° Néel on kagome).
 *   Spin liquid:                       S(k) is featureless / broadly
 *     distributed over the BZ; no Bragg peaks.
 *
 * For genuinely spin-liquid kagome, max S(k) should be of order 1,
 * not scaling with N. If max S(k) ≪ N, no long-range order.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome_structure_factor
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

/* ⟨S_i · S_j⟩ for a dense 2^N state. Spin-1/2:
 *   S_i · S_j = S^z_i S^z_j + (1/2)(S^+_i S^-_j + S^-_i S^+_j)
 * For i == j: S_i² = 3/4 (spin-1/2 eigenvalue). */
static double corr_si_sj(const double _Complex *psi, int i, int j) {
    if (i == j) return 0.75;
    long long mi = 1LL << i, mj = 1LL << j;
    double diag = 0.0, offd = 0.0;
    for (long long s = 0; s < D_FULL; ++s) {
        double a2 = creal(psi[s] * conj(psi[s]));
        int zi = (int)((s >> i) & 1), zj = (int)((s >> j) & 1);
        /* S^z S^z = ±1/4 (parallel = +, anti = −). */
        diag += ((zi == zj) ? 0.25 : -0.25) * a2;
        /* (1/2)(S^+ S^- + h.c.) flips anti-aligned pair with coeff 1/2. */
        if (zi != zj) {
            long long sp = s ^ mi ^ mj;
            /* Matrix element ⟨sp|S^+_i S^-_j + h.c.|s⟩ = 1 if (zi, zj)=(1,0)
             * (S^-_j removes j bit, S^+_i adds i bit → wait no),
             * let me reconsider. For bit bit_i=1, bit_j=0 (anti-aligned with
             * i up, j down): S^+_i S^-_j acts as (raise i already max, 0)
             *                and S^-_i S^+_j lowers i (1→0) + raises j (0→1)
             *                → flips both bits with coefficient 1.
             * So the matrix element is 1 when one of the spins flips
             * correctly. After the bit flip operator s ⊕ mi ⊕ mj:
             *   if (zi, zj) was (1, 0): bit i flips to 0, bit j flips to 1
             *     ⟨new|S^-_i S^+_j|s⟩ = 1 (coeff in (1/2)(...) term is 1/2)
             *   if (zi, zj) was (0, 1): bit i flips to 1, bit j flips to 0
             *     ⟨new|S^+_i S^-_j|s⟩ = 1 (coeff 1/2)
             * In both cases: offdiag matrix element = 1/2.
             */
            offd += 0.5 * creal(conj(psi[sp]) * psi[s]);
        }
    }
    return diag + offd;
}

int main(void) {
    printf("=== Kagome N=24 static structure factor S(k) ===\n");
    printf("    Independent spin-liquid diagnostic: featureless S(k)\n");
    printf("    ↔ no long-range order.\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);

    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
    double _Complex chi1[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);

    /* Singlet GS (F=+1 at N=24 is the lowest eigenvalue). */
    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
    double eigs[2];
    double _Complex *psi_sec = malloc((size_t)2 * sector_dim * sizeof(double _Complex));
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sector_dim, 2, 150, seed, eigs, psi_sec);
    printf("  E_0/N = %+.8f\n", eigs[0] / N_SITES);

    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    double _Complex *psi_full = malloc((size_t)D_FULL * sizeof(double _Complex));
    unfold(G, T, order, w, psi_sec, sector_dim, psi_full);

    /* Verify normalisation. */
    double nn = 0.0;
    for (long long s = 0; s < D_FULL; ++s) nn += creal(psi_full[s] * conj(psi_full[s]));
    printf("  ||ψ||² = %.10f\n", nn);

    /* Compute all ⟨S_i · S_j⟩. */
    double *C = malloc(sizeof(double) * N_SITES * N_SITES);
    printf("\n  Computing ⟨S_i·S_j⟩ for all %d pairs...\n",
           (N_SITES * (N_SITES + 1)) / 2);
    double tC = now_sec();
    for (int i = 0; i < N_SITES; ++i) {
        for (int j = i; j < N_SITES; ++j) {
            double c = corr_si_sj(psi_full, i, j);
            C[i * N_SITES + j] = c;
            C[j * N_SITES + i] = c;
        }
    }
    printf("  correlations done in %.2f s\n", now_sec() - tC);

    /* Sanity: on-site ⟨S²⟩ = 3/4. */
    for (int i = 0; i < 3; ++i)
        if (fabs(C[i * N_SITES + i] - 0.75) > 1e-10) {
            fprintf(stderr, "  WARN: ⟨S_%d·S_%d⟩ = %.6f ≠ 3/4\n", i, i,
                    C[i * N_SITES + i]);
        }

    /* Get site positions. */
    double site_x[N_SITES], site_y[N_SITES];
    for (int i = 0; i < N_SITES; ++i) {
        double xy[2];
        irrep_lattice_site_position(L, i, xy);
        site_x[i] = xy[0]; site_y[i] = xy[1];
    }

    /* Enumerate Bloch momenta on 2×4 torus. Reciprocal-lattice basis. */
    double b1[2], b2[2];
    irrep_lattice_reciprocal_vectors(L, b1, b2);
    int Lx = irrep_lattice_Lx(L), Ly = irrep_lattice_Ly(L);

    /* Max-magnitude k-point. */
    double S_max = 0, k_max[2] = {0, 0};
    int    k_max_mx = 0, k_max_my = 0;

    printf("\n  S(k) at allowed Bloch momenta:\n");
    printf("  %-8s %-25s %12s\n", "(mx,my)", "k (1/a)", "S(k)");
    printf("  --------  -------------------------  ------\n");
    for (int my = 0; my < Ly; ++my) {
        for (int mx = 0; mx < Lx; ++mx) {
            double kx = (double)mx / Lx * b1[0] + (double)my / Ly * b2[0];
            double ky = (double)mx / Lx * b1[1] + (double)my / Ly * b2[1];
            double Sk_re = 0;
            for (int i = 0; i < N_SITES; ++i)
                for (int j = 0; j < N_SITES; ++j) {
                    double dx = site_x[i] - site_x[j];
                    double dy = site_y[i] - site_y[j];
                    double phase = kx * dx + ky * dy;
                    Sk_re += C[i * N_SITES + j] * cos(phase);
                }
            Sk_re /= N_SITES;
            printf("  (%d,%d)     (%+.4f, %+.4f)             %+10.4f\n",
                   mx, my, kx, ky, Sk_re);
            if (fabs(Sk_re) > fabs(S_max)) {
                S_max = Sk_re; k_max[0] = kx; k_max[1] = ky;
                k_max_mx = mx; k_max_my = my;
            }
        }
    }

    printf("\n  max |S(k)| = %.4f at (mx,my)=(%d,%d), k=(%.4f, %.4f)\n",
           fabs(S_max), k_max_mx, k_max_my, k_max[0], k_max[1]);
    printf("\n  INTERPRETATION:\n");
    printf("    |S(k)| ≪ N=%d everywhere ⟹ no long-range order (spin liquid).\n", N_SITES);
    printf("    Peak |S(k)| / N = %.3f (dimensionless).\n", fabs(S_max) / N_SITES);
    if (fabs(S_max) / N_SITES < 0.2)
        printf("    ✓ FEATURELESS S(k) — consistent with spin-liquid phase.\n");
    else if (fabs(S_max) / N_SITES < 0.4)
        printf("    ~ mild structure — short-range correlations.\n");
    else
        printf("    ✗ Sharp peak — suggests magnetic ordering.\n");

    printf("\n  Total wall-clock: %.2f s\n", now_sec() - t0);

    free(C); free(psi_full); free(psi_sec); free(seed); free(w);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
