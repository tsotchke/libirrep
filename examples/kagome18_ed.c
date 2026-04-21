/* SPDX-License-Identifier: MIT */
/* 18-site kagome Heisenberg ED on a 2×3 torus — a larger, non-square
 * cluster where the libirrep substrate's performance posture starts to
 * matter. Hilbert space 2^18 = 262144 (≈ 4 MB per state vector); H is
 * assembled on the fly (36 NN bonds × sparse matvec per apply) and never
 * materialised dense. Shifted power iteration converges to E_0 in ~500
 * iterations; runtime is a few seconds on M2 Ultra.
 *
 * Note on symmetry: the 2×3 torus breaks C_6 (Lx ≠ Ly), so we use the
 * p1 wallpaper group (translations only) rather than p6mm. Block ED
 * in the Γ-translation sector remains available via
 * `irrep_sg_adapted_basis`, but each block is still ~45k-dim — too big
 * for our O(n³) Jacobi diagonaliser and better left to sparse Lanczos
 * by the caller. Here we just measure E_0 directly.
 *
 * Expected: E_0 / N_site in the range −0.43 to −0.44 J (18-site kagome
 * published values, Waldtmann et al. 1998, Läuchli 2011, etc.).
 *
 *   make examples
 *   ./build/bin/kagome18_ed
 */

#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>
#include <irrep/spin_project.h>

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N   18
#define DIM (1LL << N)        /* 262144 */

/* H-apply is `irrep_heisenberg_apply` from <irrep/hamiltonian.h>. */

static double power_iterate(double _Complex *psi, double _Complex *scratch,
                            const irrep_heisenberg_t *H,
                            double shift, int iters) {
    double norm = 0.0;
    for (long long s = 0; s < DIM; ++s) norm += creal(psi[s])*creal(psi[s]) + cimag(psi[s])*cimag(psi[s]);
    norm = sqrt(norm);
    if (norm < 1e-300) return NAN;
    for (long long s = 0; s < DIM; ++s) psi[s] /= norm;

    double E_prev = 0.0;
    for (int it = 0; it < iters; ++it) {
        irrep_heisenberg_apply(psi, scratch, (void *)H);
        double E_now = 0.0;
        for (long long s = 0; s < DIM; ++s) E_now += creal(conj(psi[s]) * scratch[s]);

        for (long long s = 0; s < DIM; ++s) scratch[s] = shift * psi[s] - scratch[s];
        norm = 0.0;
        for (long long s = 0; s < DIM; ++s) norm += creal(scratch[s])*creal(scratch[s]) + cimag(scratch[s])*cimag(scratch[s]);
        norm = sqrt(norm);
        if (norm < 1e-300) return NAN;
        for (long long s = 0; s < DIM; ++s) psi[s] = scratch[s] / norm;

        if (it > 10 && fabs(E_now - E_prev) < 1e-10) break;
        E_prev = E_now;
    }

    irrep_heisenberg_apply(psi, scratch, (void *)H);
    double E = 0.0;
    for (long long s = 0; s < DIM; ++s) E += creal(conj(psi[s]) * scratch[s]);
    return E;
}

int main(void) {
    printf("===  18-site kagome Heisenberg ED (2×3 torus)  ===\n\n");

    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 3);
    int Nsites = irrep_lattice_num_sites(L);
    int nb     = irrep_lattice_num_bonds_nn(L);
    printf("cluster:  %d sites,  %d NN bonds\n", Nsites, nb);
    printf("Hilbert space: 2^%d = %lld  (≈ %lld MB per state vector)\n",
           N, DIM, (long long)(DIM * sizeof(double _Complex) / (1LL << 20)));

    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);

    irrep_heisenberg_t *H = irrep_heisenberg_new(Nsites, nb, bi, bj, /*J=*/1.0);
    if (!H) { fprintf(stderr, "irrep_heisenberg_new failed\n"); return 1; }

    /* Seed in the S_z = 0 sector with a pseudo-random amplitude pattern.
     * A uniform seed can have accidental zero overlap with specific
     * momentum / parity sectors of the ground state; pseudorandom mixing
     * guarantees a non-zero overlap. */
    double _Complex *psi     = malloc(sizeof(double _Complex) * DIM);
    double _Complex *scratch = malloc(sizeof(double _Complex) * DIM);
    memset(psi, 0, sizeof(double _Complex) * DIM);
    long long num_sz0 = 0;
    uint64_t rng = 0x123456789abcdefULL;
    for (long long s = 0; s < DIM; ++s) {
        if (__builtin_popcountll(s) == N / 2) {
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            double val = (double)(rng >> 32) / (double)0xFFFFFFFFULL - 0.5;
            psi[s] = val;
            ++num_sz0;
        }
    }
    printf("S_z = 0 subspace dimension: %lld   (C(18, 9) = 48620)\n\n",
           num_sz0);

    printf("running shifted power iteration (shift = +15 J, max 2000 iters) ...\n");
    double E0 = power_iterate(psi, scratch, H,
                              /*shift=*/15.0, /*iters=*/2000);
    printf("E_0           = %+.8f J\n", E0);
    printf("E_0 / N_site  = %+.8f J\n", E0 / Nsites);
    printf("                (published 18-site kagome: −0.43 to −0.44 J)\n");

    /* J-sector check: should be pure singlet */
    double _Complex *pJ = malloc(sizeof(double _Complex) * DIM);
    printf("\nrunning total-J projection (J = 0) ...\n");
    irrep_spin_project_spin_half(/*two_J=*/0, N, 12, 6, 12, psi, pJ);
    double w_J0 = 0.0;
    for (long long s = 0; s < DIM; ++s) w_J0 += creal(pJ[s])*creal(pJ[s]) + cimag(pJ[s])*cimag(pJ[s]);
    printf("‖P_{J=0} |gs⟩‖² = %.8f   (expect ≈ 1 for a singlet ground state)\n", w_J0);
    free(pJ);

    /* Bipartite entanglement at an up-triangle cut (3 sites of cell (0,0)). */
    int sites_A[3] = { 0, 1, 2 };
    double _Complex rho_A[64];
    irrep_partial_trace(N, 2, psi, sites_A, 3, rho_A);
    double S = irrep_entropy_vonneumann(rho_A, 8);
    printf("\nbipartite entanglement:\n");
    printf("  S_VN (3-site up-triangle cut) = %.6f nats  (%.4f bits)\n",
           S, S / log(2.0));
    double tr = 0.0;
    for (int k = 0; k < 8; ++k) tr += creal(rho_A[k*8 + k]);
    printf("  Tr(ρ_A) = %.6f   (expect 1.0)\n", tr);

    /* Correlation matrix C_ij = ⟨gs|S_i · S_j|gs⟩ (reused by NN energy + S(k)). */
    static double C[18][18];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) { C[i][j] = 0.75; continue; }
            long long mi = 1LL << i, mj = 1LL << j;
            double acc = 0.0;
            for (long long s = 0; s < DIM; ++s) {
                int zi = (int)((s >> i) & 1);
                int zj = (int)((s >> j) & 1);
                double ps = creal(conj(psi[s]) * psi[s]);
                acc += 0.25 * (zi == zj ? +1.0 : -1.0) * ps;
                if (zi != zj) {
                    long long s_flip = s ^ mi ^ mj;
                    acc += 0.5 * creal(conj(psi[s]) * psi[s_flip]);
                }
            }
            C[i][j] = acc;
        }
    }

    /* NN bond energy consistency check */
    double avg_nn = 0.0;
    for (int b = 0; b < nb; ++b) avg_nn += C[bi[b]][bj[b]];
    avg_nn /= nb;
    printf("\n⟨S_i·S_j⟩_NN = %+.6f J   (per-bond energy)\n", avg_nn);
    printf("E_0 / nb     = %+.6f J   (consistency check)\n", E0 / nb);

    /* Static structure factor S(k) over the 2×3 BZ (6 k-points). */
    double pos[18][2];
    for (int s = 0; s < N; ++s) irrep_lattice_site_position(L, s, pos[s]);
    double b1[2], b2[2];
    irrep_lattice_reciprocal_vectors(L, b1, b2);
    int Lx = irrep_lattice_Lx(L), Ly = irrep_lattice_Ly(L);

    printf("\nstatic structure factor S(k) on the 2×3 BZ:\n");
    printf("  %-6s  %9s  %9s     %9s\n", "(n1,n2)", "kx·a", "ky·a", "S(k)");
    for (int n2 = 0; n2 < Ly; ++n2) {
        for (int n1 = 0; n1 < Lx; ++n1) {
            double kx = (n1 / (double)Lx) * b1[0] + (n2 / (double)Ly) * b2[0];
            double ky = (n1 / (double)Lx) * b1[1] + (n2 / (double)Ly) * b2[1];
            double _Complex sum = 0.0;
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    double dot = kx * (pos[i][0] - pos[j][0]) + ky * (pos[i][1] - pos[j][1]);
                    sum += cexp(I * dot) * C[i][j];
                }
            }
            double Sk = creal(sum) / (double)N;
            printf("  (%d,%d)   %+9.4f  %+9.4f     %9.6f\n",
                   n1, n2, kx, ky, Sk);
        }
    }

    /* Lowest triplet via S_z = 1 seed: popcount = 9 − 1 = 8. */
    printf("\nrunning S_z = 1 power iteration for lowest triplet ...\n");
    double _Complex *psi_trip = malloc(sizeof(double _Complex) * DIM);
    memset(psi_trip, 0, sizeof(double _Complex) * DIM);
    rng = 0xbeef1234beefULL;
    for (long long s = 0; s < DIM; ++s) {
        if (__builtin_popcountll(s) == N / 2 - 1) {
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            psi_trip[s] = (double)(rng >> 32) / (double)0xFFFFFFFFULL - 0.5;
        }
    }
    double E_trip = power_iterate(psi_trip, scratch, H,
                                  /*shift=*/15.0, /*iters=*/2000);
    double spin_gap = E_trip - E0;
    printf("E_triplet         = %+.8f J\n", E_trip);
    printf("spin gap  Δ_S = E_triplet − E_0 = %+.8f J\n", spin_gap);
    printf("                      (12-site value for reference: 0.3827 J)\n");
    free(psi_trip);

    free(psi); free(scratch);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_lattice_free(L);
    return 0;
}
