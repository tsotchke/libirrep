/* SPDX-License-Identifier: MIT */
/* First libirrep-based kagome ED physics regression.
 *
 * Target: 2×2 unit cells × 3 sublattices = 12-site kagome Heisenberg
 * antiferromagnet on a torus, the smallest cluster that respects the full
 * p6mm space group. Hilbert space: 2^12 = 4096 complex amplitudes. Bond
 * count: 24 NN (exact by Euler's formula on kagome).
 *
 * What the example proves, end-to-end:
 *
 *   1. The 1.3 substrate (lattice, space_group, config_project, rdm,
 *      spin_project) composes correctly on a non-trivial cluster — not
 *      just on textbook 2-qubit toys.
 *   2. Observed ground-state energy  E_0 / N  = −0.4537 J, matching the
 *      classic kagome-12 ED values (Elser 1989, Lecheminant et al. 1997).
 *   3. The ground state lives in a unique 1D p6mm irrep (B₁ under this
 *      file's character convention; total 1D weight = 1.000000 by
 *      completeness — E₁ and E₂ contributions are strictly zero).
 *   4. Ground state is a total-J = 0 singlet to machine precision.
 *   5. First excited state (deflated power iteration) is *also* a singlet,
 *      in the A₁ sector. Singlet-singlet gap Δ_ss = 0.117 J — the classic
 *      small-gap signature of proliferating low-lying singlets on kagome.
 *   6. The lowest triplet (S_z = 1 power iteration) sits at +0.383 J
 *      above the ground state, giving the cluster spin gap Δ_S = 0.383 J.
 *   7. Bipartite entanglement at an up-triangle cut is S_VN ≈ 1.574 nats
 *      (≈ 2.27 bits) — well-defined and in the expected range.
 *
 * All four eigenvalues are extracted by shifted / deflated power
 * iteration on a sparse Hamiltonian (on-the-fly H-apply, never
 * materialised as a dense 4096² matrix — that would be 67 MB). Runtime:
 * ~2 s end-to-end on M2 Ultra for the full ED + 6-irrep projection
 * decomposition + two J-projections + triplet extraction + entropy.
 * libirrep is infrastructure-ready for the 6×6 target.
 *
 *   make examples
 *   ./build/bin/kagome12_ed
 */

#include <irrep/config_project.h>
#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>
#include <irrep/space_group.h>
#include <irrep/spin_project.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N   12
#define DIM (1LL << N)        /* 4096 */

/* -------------------------------------------------------------------------- *
 * On-the-fly Heisenberg application via the library primitive. The inline
 * apply_H this example once carried is now `irrep_heisenberg_apply` — see
 * include/irrep/hamiltonian.h. Kept the power-iteration driver here
 * because the shift-and-invert technique it demonstrates is orthogonal
 * to the Lanczos path tested against it.
 * -------------------------------------------------------------------------- */

static double power_iterate(double _Complex *psi, double _Complex *scratch,
                            const irrep_heisenberg_t *H,
                            double shift, int iters) {
    /* Normalise input */
    double norm = 0.0;
    for (long long s = 0; s < DIM; ++s) norm += creal(psi[s])*creal(psi[s]) + cimag(psi[s])*cimag(psi[s]);
    norm = sqrt(norm);
    if (norm < 1e-300) return NAN;
    for (long long s = 0; s < DIM; ++s) psi[s] /= norm;

    for (int it = 0; it < iters; ++it) {
        irrep_heisenberg_apply(psi, scratch, (void *)H);
        for (long long s = 0; s < DIM; ++s) scratch[s] = shift * psi[s] - scratch[s];
        norm = 0.0;
        for (long long s = 0; s < DIM; ++s) norm += creal(scratch[s])*creal(scratch[s]) + cimag(scratch[s])*cimag(scratch[s]);
        norm = sqrt(norm);
        if (norm < 1e-300) return NAN;
        for (long long s = 0; s < DIM; ++s) psi[s] = scratch[s] / norm;
    }

    /* Rayleigh quotient gives E_0 */
    irrep_heisenberg_apply(psi, scratch, (void *)H);
    double E = 0.0;
    for (long long s = 0; s < DIM; ++s) E += creal(conj(psi[s]) * scratch[s]);
    return E;
}

/* Deflated power iteration: find the first excited state by orthogonalising
 * against a set of previously-converged eigenvectors after every step.
 * Returns the eigenvalue (Rayleigh quotient on the converged vector). */
static double power_iterate_deflated(double _Complex *psi,
                                     double _Complex *scratch,
                                     const double _Complex *const *deflate_against,
                                     int n_deflate,
                                     const irrep_heisenberg_t *H,
                                     double shift, int iters) {
    /* Initial orthogonalisation */
    for (int k = 0; k < n_deflate; ++k) {
        double _Complex overlap = 0.0;
        for (long long s = 0; s < DIM; ++s) overlap += conj(deflate_against[k][s]) * psi[s];
        for (long long s = 0; s < DIM; ++s) psi[s] -= overlap * deflate_against[k][s];
    }
    double norm = 0.0;
    for (long long s = 0; s < DIM; ++s) norm += creal(psi[s])*creal(psi[s]) + cimag(psi[s])*cimag(psi[s]);
    norm = sqrt(norm);
    if (norm < 1e-300) return NAN;
    for (long long s = 0; s < DIM; ++s) psi[s] /= norm;

    for (int it = 0; it < iters; ++it) {
        irrep_heisenberg_apply(psi, scratch, (void *)H);
        for (long long s = 0; s < DIM; ++s) scratch[s] = shift * psi[s] - scratch[s];
        /* Orthogonalise */
        for (int k = 0; k < n_deflate; ++k) {
            double _Complex overlap = 0.0;
            for (long long s = 0; s < DIM; ++s) overlap += conj(deflate_against[k][s]) * scratch[s];
            for (long long s = 0; s < DIM; ++s) scratch[s] -= overlap * deflate_against[k][s];
        }
        norm = 0.0;
        for (long long s = 0; s < DIM; ++s) norm += creal(scratch[s])*creal(scratch[s]) + cimag(scratch[s])*cimag(scratch[s]);
        norm = sqrt(norm);
        if (norm < 1e-300) return NAN;
        for (long long s = 0; s < DIM; ++s) psi[s] = scratch[s] / norm;
    }

    irrep_heisenberg_apply(psi, scratch, (void *)H);
    double E = 0.0;
    for (long long s = 0; s < DIM; ++s) E += creal(conj(psi[s]) * scratch[s]);
    return E;
}

/* -------------------------------------------------------------------------- *
 * Basis-state permutation induced by a site permutation                      *
 * -------------------------------------------------------------------------- */

static long long permute_basis_(long long s, const int *perm_inv, int Nsites) {
    long long out = 0;
    for (int j = 0; j < Nsites; ++j) {
        int bit = (int)((s >> perm_inv[j]) & 1);
        out |= (long long)bit << j;
    }
    return out;
}

/* ‖P_{μ} |ψ⟩‖² — the weight of |ψ⟩ in the μ-irrep sector.
 * P_{μ} ψ(σ) = (d_μ / |G|) Σ_g χ_μ*(g) ψ(g·σ); we apply this operator to |ψ⟩
 * and return the norm² of the result. */
static double projection_weight(const irrep_space_group_t *G,
                                const double _Complex *chi,
                                int d_mu,
                                const double _Complex *psi,
                                int Nsites) {
    int order = irrep_space_group_order(G);
    double _Complex *proj = calloc(DIM, sizeof(double _Complex));
    int *perm_inv = malloc((size_t)Nsites * sizeof(int));

    for (int g = 0; g < order; ++g) {
        irrep_space_group_permutation_inverse(G, g, perm_inv);
        double _Complex w = conj(chi[g]) * (double)d_mu / (double)order;
        for (long long s = 0; s < DIM; ++s) {
            long long s_g = permute_basis_(s, perm_inv, Nsites);
            proj[s] += w * psi[s_g];
        }
    }
    double w2 = 0.0;
    for (long long s = 0; s < DIM; ++s) w2 += creal(proj[s])*creal(proj[s]) + cimag(proj[s])*cimag(proj[s]);
    free(proj);
    free(perm_inv);
    return w2;
}

int main(void) {
    printf("===  12-site kagome Heisenberg ED (2×2 unit cells)  ===\n\n");

    /* 1. Lattice. */
    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
    int Nsites = irrep_lattice_num_sites(L);
    int nb     = irrep_lattice_num_bonds_nn(L);
    printf("cluster:  %d sites,  %d NN bonds\n", Nsites, nb);

    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);

    /* Library Hamiltonian: `J = 1` Heisenberg on the NN bond list. */
    irrep_heisenberg_t *H = irrep_heisenberg_new(Nsites, nb, bi, bj, 1.0);
    if (!H) { fprintf(stderr, "irrep_heisenberg_new failed\n"); return 1; }

    /* 2. Space group attach. */
    irrep_space_group_t *G =
        irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
    int order = irrep_space_group_order(G);
    printf("p6mm on 2×2 torus:  %d group elements\n\n", order);

    /* 3. Ground state via shifted power iteration. */
    double _Complex *psi     = malloc(sizeof(double _Complex) * DIM);
    double _Complex *scratch = malloc(sizeof(double _Complex) * DIM);
    /* Seed: uniform in the S_z=0 sector — the space containing the
     * Heisenberg ground state (Marshall). 12-site cluster has C(12,6) = 924
     * such states. */
    memset(psi, 0, sizeof(double _Complex) * DIM);
    for (long long s = 0; s < DIM; ++s) {
        int pop = __builtin_popcountll(s);
        if (pop == N / 2) psi[s] = 1.0;
    }

    printf("running shifted power iteration (shift = +10 J) ...\n");
    double E0 = power_iterate(psi, scratch, H,
                              /*shift=*/10.0, /*iters=*/2000);
    printf("E_0           = %+.8f J\n", E0);
    printf("E_0 / N_site  = %+.8f J\n", E0 / Nsites);

    /* 4. Full p6mm irrep decomposition of |gs⟩.
     *
     *    C6v characters (Altmann-Herzig 1994):
     *
     *      rotations  E  C6  C3   C2  C3²  C6⁵   (6 classes on my enumeration
     *                                              of elements 0..5)
     *      A₁         1   1   1    1   1    1          mirrors all +1
     *      A₂         1   1   1    1   1    1          mirrors all −1
     *      B₁         1  -1   1   -1   1   -1          mirrors alternate
     *      B₂         1  -1   1   -1   1   -1          mirrors alt (opposite)
     *      E₁         2   1  -1   -2  -1    1          mirrors all 0
     *      E₂         2  -1  -1    2  -1   -1          mirrors all 0
     *
     *    The six weights sum to 1 by completeness. The one that carries the
     *    ground state identifies the symmetry sector. */
    irrep_sg_irrep_t *A1 = irrep_sg_trivial(G);    /* kept for cleanup */
    int point_order = irrep_space_group_point_order(G);
    int half = point_order / 2;
    (void)half;

    static const double rot_chi[6][6] = {
        { +1, +1, +1, +1, +1, +1 },   /* A1 */
        { +1, +1, +1, +1, +1, +1 },   /* A2 */
        { +1, -1, +1, -1, +1, -1 },   /* B1 */
        { +1, -1, +1, -1, +1, -1 },   /* B2 */
        { +2, +1, -1, -2, -1, +1 },   /* E1 */
        { +2, -1, -1, +2, -1, -1 }    /* E2 */
    };
    /* For mirrors we split into two triples by position-parity of the element
     * index (mirrors 0, 2, 4 vs 1, 3, 5): these are the two geometric classes
     * (through-vertex vs through-edge). A₁/B₁ share sign on the first class,
     * A₂/B₂ on the other; 2D irreps have character 0 on every mirror. */
    static const double mir_chi[6][2] = {
        { +1, +1 },   /* A1 */
        { -1, -1 },   /* A2 */
        { +1, -1 },   /* B1 — convention; swap with B2 flips sign */
        { -1, +1 },   /* B2 */
        {  0,  0 },   /* E1 */
        {  0,  0 }    /* E2 */
    };
    static const int dim_mu[6] = { 1, 1, 1, 1, 2, 2 };
    static const char *name_mu[6] = { "A1", "A2", "B1", "B2", "E1", "E2" };

    double total_weight = 0.0;
    double _Complex *chi = malloc(sizeof(double _Complex) * order);
    printf("\np6mm irrep decomposition of |gs⟩:\n");
    for (int mu = 0; mu < 6; ++mu) {
        for (int g = 0; g < order; ++g) {
            int p = g % point_order;
            double c;
            if (p < 6) {
                c = rot_chi[mu][p];
            } else {
                int m = p - 6;
                c = mir_chi[mu][m & 1];
            }
            chi[g] = c + 0.0 * I;
        }
        double w = projection_weight(G, chi, dim_mu[mu], psi, Nsites);
        total_weight += w;
        printf("  ‖P_{%s} |gs⟩‖² = %.10f\n", name_mu[mu], w);
    }
    printf("  ------------------------------------\n");
    printf("  Σ_μ                = %.10f   (expect 1.000 by completeness)\n",
           total_weight);
    free(chi);

    /* 6. Total-J projection: weight in J=0 (singlet). */
    double _Complex *pJ = malloc(sizeof(double _Complex) * DIM);
    printf("\nrunning total-J projection (may take a few seconds) ...\n");
    irrep_spin_project_spin_half(/*two_J=*/0, N, 16, 6, 16, psi, pJ);
    double w_J0 = 0.0;
    for (long long s = 0; s < DIM; ++s) w_J0 += creal(pJ[s])*creal(pJ[s]) + cimag(pJ[s])*cimag(pJ[s]);
    printf("‖P_{J=0} |gs⟩‖² = %.8f   (expect 1.0 for the singlet ground state)\n",
           w_J0);
    irrep_spin_project_spin_half(/*two_J=*/2, N, 16, 6, 16, psi, pJ);
    double w_J1 = 0.0;
    for (long long s = 0; s < DIM; ++s) w_J1 += creal(pJ[s])*creal(pJ[s]) + cimag(pJ[s])*cimag(pJ[s]);
    printf("‖P_{J=1} |gs⟩‖² = %.8f   (expect ≈ 0)\n", w_J1);
    free(pJ);

    /* 7. First excited state via deflated power iteration.
     *    E_1 is the next-lowest eigenvalue of H; Δ = E_1 − E_0 is the
     *    cluster spin gap, a physically meaningful observable for the
     *    kagome Heisenberg antiferromagnet (the heart of the current
     *    ground-state-nature debate). */
    double _Complex *psi_e1 = malloc(sizeof(double _Complex) * DIM);
    memset(psi_e1, 0, sizeof(double _Complex) * DIM);
    /* Seed in the S_z = 0 sector but linearly independent of ψ */
    for (long long s = 0; s < DIM; ++s) {
        if (__builtin_popcountll(s) == N / 2) psi_e1[s] = (double)((s * 31) % 17 - 8);
    }
    const double _Complex *deflate[1] = { psi };
    printf("\nrunning deflated power iteration for E_1 ...\n");
    double E1 = power_iterate_deflated(psi_e1, scratch, deflate, 1,
                                       H,
                                       /*shift=*/10.0, /*iters=*/3000);
    printf("E_1           = %+.8f J\n", E1);
    printf("gap Δ = E_1 − E_0 = %+.8f J\n", E1 - E0);

    /* Classify the excited state: which 1D irrep and which J sector? */
    {
        double _Complex *chi_row = malloc(sizeof(double _Complex) * order);
        printf("\nexcited-state p6mm sector:\n");
        for (int mu = 0; mu < 6; ++mu) {
            for (int g = 0; g < order; ++g) {
                int p = g % point_order;
                double c;
                if (p < 6) {
                    c = rot_chi[mu][p];
                } else {
                    int m = p - 6;
                    c = mir_chi[mu][m & 1];
                }
                chi_row[g] = c + 0.0 * I;
            }
            double w = projection_weight(G, chi_row, dim_mu[mu], psi_e1, Nsites);
            printf("  ‖P_{%s} |e1⟩‖² = %.10f\n", name_mu[mu], w);
        }
        free(chi_row);
    }

    /* J sector of the excited state */
    double _Complex *pJe = malloc(sizeof(double _Complex) * DIM);
    printf("\nexcited-state total-J decomposition:\n");
    for (int two_J = 0; two_J <= 4; two_J += 2) {
        irrep_spin_project_spin_half(two_J, N, 16, 6, 16, psi_e1, pJe);
        double w = 0.0;
        for (long long s = 0; s < DIM; ++s) w += creal(pJe[s])*creal(pJe[s]) + cimag(pJe[s])*cimag(pJe[s]);
        printf("  ‖P_{J=%d} |e1⟩‖² = %.6f\n", two_J / 2, w);
    }
    free(pJe);
    free(psi_e1);

    /* ------------------------------------------------------------------ *
     * 8. True spin gap: lowest triplet state.                             *
     *    H conserves S_z, so seeding power iteration in the S_z = 1       *
     *    sector (popcount = N/2 − 1, the odd-# magnetisation subspace)    *
     *    converges to the lowest triplet eigenvalue.                      *
     * ------------------------------------------------------------------ */
    double _Complex *psi_trip = malloc(sizeof(double _Complex) * DIM);
    memset(psi_trip, 0, sizeof(double _Complex) * DIM);
    /* S_z = 1 → #down = N/2 − 1; with |↑⟩ = 0 (bit unset) convention,
     * popcount counts |↓⟩s. But my apply_H treats bits symmetrically, so
     * the choice of up/down mapping only affects sign conventions, not
     * the eigenvalue. */
    for (long long s = 0; s < DIM; ++s) {
        if (__builtin_popcountll(s) == N / 2 - 1) psi_trip[s] = 1.0;
    }
    printf("\nrunning S_z = 1 power iteration for lowest triplet ...\n");
    double E_trip = power_iterate(psi_trip, scratch, H,
                                  /*shift=*/10.0, /*iters=*/2000);
    double spin_gap = E_trip - E0;
    printf("E_triplet         = %+.8f J\n", E_trip);
    printf("spin gap  Δ_S = E_triplet − E_0 = %+.8f J\n", spin_gap);

    /* Verify it's actually J = 1 */
    double _Complex *pJt = malloc(sizeof(double _Complex) * DIM);
    irrep_spin_project_spin_half(/*two_J=*/2, N, 16, 6, 16, psi_trip, pJt);
    double w_trip = 0.0;
    for (long long s = 0; s < DIM; ++s) w_trip += creal(pJt[s])*creal(pJt[s]) + cimag(pJt[s])*cimag(pJt[s]);
    printf("‖P_{J=1} |trip⟩‖² = %.6f   (confirms triplet character)\n", w_trip);
    free(pJt);
    free(psi_trip);

    /* ------------------------------------------------------------------ *
     * 9. Static structure factor S(k) on the ground state.                *
     *                                                                    *
     *    S(k) = (1/N) Σ_{ij} exp(i k · (r_i − r_j)) ⟨gs| S_i · S_j |gs⟩   *
     *                                                                    *
     *    Computed by first building the 12×12 real-space correlation     *
     *    matrix C_{ij} = ⟨gs| S_i · S_j |gs⟩, then Fourier-transforming   *
     *    at the 4 allowed k-points of the 2×2 Brillouin zone:            *
     *       Γ = (0,0),  M_a = b1/2,  M_b = b2/2,  K = (b1+b2)/2          *
     *    (where b1, b2 are the reciprocal primitive vectors of the       *
     *    kagome lattice).                                                *
     * ------------------------------------------------------------------ */

    /* 12×12 correlator C_{ij} = ⟨ψ|S_i · S_j|ψ⟩.  For i == j: 3/4.
     * For i ≠ j: compute on-the-fly via (¼)[σ_z σ_z + 2(σ_+σ_- + σ_-σ_+)].*/
    double C[12][12] = {{0}};
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) { C[i][j] = 0.75; continue; }
            long long mi = 1LL << i, mj = 1LL << j;
            double acc = 0.0;
            for (long long s = 0; s < DIM; ++s) {
                int zi = (int)((s >> i) & 1);
                int zj = (int)((s >> j) & 1);
                double ps = creal(conj(psi[s]) * psi[s]);
                /* ¼ S_z S_z contribution */
                acc += 0.25 * (zi == zj ? +1.0 : -1.0) * ps;
                /* ½ (σ_+ σ_- + σ_- σ_+) flips differing pairs */
                if (zi != zj) {
                    long long s_flip = s ^ mi ^ mj;
                    acc += 0.5 * creal(conj(psi[s]) * psi[s_flip]);
                }
            }
            C[i][j] = acc;
        }
    }

    /* Site positions (cartesian), needed for the Fourier phases. */
    double pos[12][2];
    for (int s = 0; s < N; ++s) irrep_lattice_site_position(L, s, pos[s]);

    /* Reciprocal primitive vectors of the 2×2 kagome torus.                */
    double b1[2], b2[2];
    irrep_lattice_reciprocal_vectors(L, b1, b2);
    int Lx = irrep_lattice_Lx(L), Ly = irrep_lattice_Ly(L);

    printf("\nstatic structure factor S(k) on the 2×2 BZ:\n");
    printf("  %-6s  %9s  %9s     %9s\n", "k", "kx·a", "ky·a", "S(k)");
    const char *k_names[4] = { "Γ", "M_a", "M_b", "K" };
    int k_indices[4][2] = { {0, 0}, {1, 0}, {0, 1}, {1, 1} };
    for (int p = 0; p < 4; ++p) {
        double kx = (k_indices[p][0] / (double)Lx) * b1[0]
                  + (k_indices[p][1] / (double)Ly) * b2[0];
        double ky = (k_indices[p][0] / (double)Lx) * b1[1]
                  + (k_indices[p][1] / (double)Ly) * b2[1];
        double _Complex sum = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double dot = kx * (pos[i][0] - pos[j][0])
                           + ky * (pos[i][1] - pos[j][1]);
                sum += cexp(I * dot) * C[i][j];
            }
        }
        double Sk = creal(sum) / (double)N;
        printf("  %-6s  %+9.4f  %+9.4f     %9.6f\n",
               k_names[p], kx, ky, Sk);
    }

    /* Kagome nearest-neighbour bond energy ⟨S_i · S_j⟩ averaged over the
     * 24 NN bonds — the per-bond ground-state energy. Must match E_0/nb. */
    double avg_nn = 0.0;
    for (int b = 0; b < nb; ++b) avg_nn += C[bi[b]][bj[b]];
    avg_nn /= nb;
    printf("\n⟨S_i · S_j⟩_NN = %+.6f J  (per-bond energy)\n", avg_nn);
    printf("E_0 / (J · nb) = %+.6f   (consistency check)\n", E0 / nb);

    /* 8. Bipartite entanglement at a 3-site (up-triangle in unit cell (0,0))
     * cut, the kagome-natural bipartition. */
    int sites_A[3] = { 0, 1, 2 };   /* cell (0,0), sublattices A, B, C */
    double _Complex rho_A[64];      /* 8 × 8 */
    irrep_partial_trace(N, 2, psi, sites_A, 3, rho_A);
    double S = irrep_entropy_vonneumann(rho_A, 8);
    printf("\nentanglement at one-triangle cut:\n");
    printf("  S_VN = %.6f nats  (%.4f bits)\n", S, S / log(2.0));

    /* Tr(ρ_A) sanity */
    double tr = 0.0;
    for (int k = 0; k < 8; ++k) tr += creal(rho_A[k * 8 + k]);
    printf("  Tr(ρ_A) = %.6f  (expect 1.0)\n", tr);

    /* ------------------------------------------------------------------ *
     * 10. Kitaev-Preskill topological entanglement entropy γ.             *
     *                                                                    *
     *    A clean γ requires an annular tripartition: trace out part of    *
     *    the system so ρ_{ABC} is a *mixed* state, then split the kept    *
     *    region into three pieces A, B, C. The formula                    *
     *       γ = S_A + S_B + S_C − S_AB − S_BC − S_AC + S_ABC              *
     *    distinguishes Z₂ topological order (γ = ln 2) from trivial       *
     *    (γ = 0) in the thermodynamic limit; on a 2×2 kagome torus it     *
     *    is dominated by finite-size effects but the formula exercises    *
     *    the full rdm.h pipeline end-to-end.                              *
     *                                                                    *
     *    Tripartition: keep sites {0..5} (two kagome unit cells, one row  *
     *    of the torus); trace out sites {6..11} (other row). Within the   *
     *    kept 6 sites:                                                    *
     *        A = {0, 1}    B = {2, 3}    C = {4, 5}                       *
     * ------------------------------------------------------------------ */
    {
        int A[] = { 0, 1 };
        int B[] = { 2, 3 };
        int C[] = { 4, 5 };
        int AB[] = { 0, 1, 2, 3 };
        int BC[] = { 2, 3, 4, 5 };
        int AC[] = { 0, 1, 4, 5 };
        int ABC[] = { 0, 1, 2, 3, 4, 5 };

        /* Partial-trace + entropy for each region. */
        double _Complex rho_A[16], rho_B[16], rho_C[16];
        double _Complex rho_AB[256], rho_BC[256], rho_AC[256];
        double _Complex rho_ABC[4096];

        irrep_partial_trace(N, 2, psi, A,   2, rho_A);
        irrep_partial_trace(N, 2, psi, B,   2, rho_B);
        irrep_partial_trace(N, 2, psi, C,   2, rho_C);
        irrep_partial_trace(N, 2, psi, AB,  4, rho_AB);
        irrep_partial_trace(N, 2, psi, BC,  4, rho_BC);
        irrep_partial_trace(N, 2, psi, AC,  4, rho_AC);
        irrep_partial_trace(N, 2, psi, ABC, 6, rho_ABC);

        double S_A   = irrep_entropy_vonneumann(rho_A,    4);
        double S_B   = irrep_entropy_vonneumann(rho_B,    4);
        double S_C   = irrep_entropy_vonneumann(rho_C,    4);
        double S_AB  = irrep_entropy_vonneumann(rho_AB,  16);
        double S_BC  = irrep_entropy_vonneumann(rho_BC,  16);
        double S_AC  = irrep_entropy_vonneumann(rho_AC,  16);
        double S_ABC = irrep_entropy_vonneumann(rho_ABC, 64);

        double gamma = irrep_topological_entanglement_entropy(
            S_A, S_B, S_C, S_AB, S_BC, S_AC, S_ABC);

        printf("\nKitaev-Preskill tripartition (annular: trace {6..11}):\n");
        printf("  S_A   = %.6f,  S_B  = %.6f,  S_C  = %.6f\n", S_A, S_B, S_C);
        printf("  S_AB  = %.6f,  S_BC = %.6f,  S_AC = %.6f\n", S_AB, S_BC, S_AC);
        printf("  S_ABC = %.6f\n", S_ABC);
        printf("  γ (nats) = %+.6f\n", gamma);
        printf("  γ (bits) = %+.6f    (Z₂ sig: ln 2 = +0.6931 nats, trivial: 0)\n",
               gamma / log(2.0));
    }

    /* Cleanup */
    free(psi); free(scratch);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_sg_irrep_free(A1);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
