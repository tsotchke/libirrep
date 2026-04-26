/* SPDX-License-Identifier: MIT */
/* End-to-end entanglement-entropy calculation on the kagome-Heisenberg
 * ground state. Demonstrates the complete research pipeline:
 *
 *   1. Build lattice + space group + Hamiltonian
 *   2. Build representative table at fixed Sz (popcount = N/2)
 *   3. Build dense Γ-A_1 sector Hamiltonian (small N, fast path)
 *   4. Diagonalise to get E_0 + ground-state eigenvector in sector basis
 *   5. Unfold the sector-basis GS to a full 2^N state vector
 *   6. Compute reduced density matrices ρ_A for several bipartitions
 *   7. Compute von Neumann entropy S_A for each
 *
 * This is the first data point in an N-scaling series for the Kitaev-
 * Preskill topological entanglement entropy γ (S_A + S_B + S_C - S_AB -
 * S_BC - S_AC + S_ABC for three non-overlapping regions). At N=12 the
 * number won't be converged (log 2 for Z_2 spin liquid, 0 for Dirac
 * gapless), but the pipeline produces a well-defined finite-size value.
 *
 *   make examples
 *   ./build/bin/kagome12_entanglement
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

#define N_SITES 12
#define D_FULL  (1LL << N_SITES)

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Get the ground-state (lowest-eigenvalue) eigenvector via the Jacobi
 * eigendecomposition. For small n this is O(n^3) and converges to
 * machine precision in a few sweeps — more robust than power iteration
 * on near-degenerate spectra like kagome's singlet tower. */
static void dense_gs_eigenvector(int n, const double _Complex *H_in,
                                 double _Complex *psi_out, double *E0_out) {
    double _Complex *A = malloc((size_t)n * n * sizeof(double _Complex));
    double _Complex *V = malloc((size_t)n * n * sizeof(double _Complex));
    double          *w = malloc((size_t)n * sizeof(double));
    memcpy(A, H_in, (size_t)n * n * sizeof(double _Complex));
    irrep_hermitian_eigendecomp(n, A, w, V);
    /* Eigenvalues sorted descending; lowest = w[n-1]; vector at column n-1. */
    *E0_out = w[n - 1];
    for (int i = 0; i < n; ++i)
        psi_out[i] = V[(size_t)i * n + (n - 1)];
    free(A); free(V); free(w);
}

/* Unfold a sector-basis coefficient vector into a dense 2^N computational-
 * basis state. For each rep u in the filtered table, construct |ũ⟩ =
 * (1/√σ_u) Σ_g w_g |g·u⟩ explicitly and accumulate c_k · |ũ_k⟩. */
static void unfold_sector_to_dense(const irrep_space_group_t *G,
                                   const irrep_sg_rep_table_t *T,
                                   const irrep_sg_little_group_t *lg,
                                   const irrep_sg_little_group_irrep_t *mu_k,
                                   const double _Complex *psi_sector,
                                   int sector_dim,
                                   double _Complex *psi_full /* length 2^N */) {
    memset(psi_full, 0, (size_t)D_FULL * sizeof(double _Complex));
    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu_k, w);

    double _Complex *basis_vec = malloc((size_t)D_FULL * sizeof(double _Complex));

    long long n_reps = irrep_sg_rep_table_count(T);
    int written = 0;
    for (long long k = 0; k < n_reps && written < sector_dim; ++k) {
        uint64_t u = irrep_sg_rep_table_get(T, k);
        /* Build raw projected vector. */
        memset(basis_vec, 0, (size_t)D_FULL * sizeof(double _Complex));
        for (int g = 0; g < order; ++g) {
            uint64_t gx = irrep_space_group_apply_bits(G, g, u);
            basis_vec[gx] += w[g];
        }
        double nn = 0.0;
        for (long long s = 0; s < D_FULL; ++s)
            nn += creal(basis_vec[s] * conj(basis_vec[s]));
        if (nn < 1e-20) continue; /* annihilated rep — skipped by dense builder */
        double inv = 1.0 / sqrt(nn);
        double _Complex coef = psi_sector[written];
        for (long long s = 0; s < D_FULL; ++s)
            psi_full[s] += coef * inv * basis_vec[s];
        ++written;
    }
    free(basis_vec); free(w);
}

int main(void) {
    printf("=== Kagome-Heisenberg N=12 entanglement-entropy pipeline ===\n\n");

    double t0 = now_sec();

    /* 1. Lattice + Hamiltonian. */
    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);

    /* 2. Rep table at Sz = 0 (popcount = 6). */
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 6);
    long long cap = irrep_sg_rep_table_count(T);
    printf("  rep table: %lld reps at popcount = 6 (Sz = 0)\n", cap);

    /* 3. Dense Γ-A_1 sector Hamiltonian. */
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
    irrep_sg_little_group_irrep_t *A1 =
        irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_A1);

    double _Complex *Hdense_full = malloc((size_t)cap * cap * sizeof(double _Complex));
    int sector_dim = irrep_sg_heisenberg_sector_build_dense(H, T, lg, A1, Hdense_full, (int)cap);
    printf("  sector (Γ, A_1) dim = %d\n", sector_dim);

    /* Extract compact dim × dim block. */
    double _Complex *Hblock = malloc((size_t)sector_dim * sector_dim * sizeof(double _Complex));
    for (int r = 0; r < sector_dim; ++r)
        for (int c = 0; c < sector_dim; ++c)
            Hblock[(size_t)r * sector_dim + c] = Hdense_full[(size_t)r * cap + c];

    /* 4. Eigenvector via power iteration. */
    double _Complex *psi_sector = malloc((size_t)sector_dim * sizeof(double _Complex));
    double E0;
    dense_gs_eigenvector(sector_dim, Hblock, psi_sector, &E0);
    printf("  E_0 = %+.8f J   (E_0 / N = %+.8f)\n", E0, E0 / (double)N_SITES);

    /* Verify against reference (kagome12_k_resolved_ed: Γ-A_1 gives -5.328392). */
    if (fabs(E0 - (-5.328392)) < 1e-4)
        printf("  (matches published Γ-A_1 E_0 = -5.328392 J)\n");
    else
        printf("  (WARNING: expected -5.328392, got %.6f)\n", E0);

    /* 5. Unfold to 2^N = 4096. */
    double _Complex *psi_full = malloc((size_t)D_FULL * sizeof(double _Complex));
    unfold_sector_to_dense(G, T, lg, A1, psi_sector, sector_dim, psi_full);
    /* Verify: should be normalised. */
    double norm2 = 0.0;
    for (long long s = 0; s < D_FULL; ++s)
        norm2 += creal(psi_full[s] * conj(psi_full[s]));
    printf("  unfold: ||ψ_full||² = %.10f (should be 1.0)\n", norm2);

    /* 6. Entanglement entropy vs subsystem size. For each |A| ∈ [1, 6],
     * pick sites [0, |A|) as region A, remaining as B. */
    printf("\n  Entanglement entropy S_A vs subsystem size |A|:\n");
    printf("  %3s   %12s\n", "|A|", "S_A");
    printf("  ---   ------------\n");
    for (int nA = 1; nA <= 6; ++nA) {
        int sites_A[12];
        for (int i = 0; i < nA; ++i) sites_A[i] = i;
        long long dim_A = 1LL << nA;
        double _Complex *rho_A = malloc((size_t)dim_A * dim_A * sizeof(double _Complex));
        irrep_partial_trace(N_SITES, 2, psi_full, sites_A, nA, rho_A);
        double S = irrep_entropy_vonneumann(rho_A, (int)dim_A);
        printf("  %3d   %+12.8f\n", nA, S);
        free(rho_A);
    }

    /* 7. Kitaev-Preskill three-region γ extraction.
     *   γ = S_A + S_B + S_C − S_{A∪B} − S_{A∪C} − S_{B∪C} + S_{A∪B∪C}
     * Three disjoint 3-site regions leave 3 sites as environment.
     * Thermodynamic value: log(2) for gapped Z_2 topological order,
     * 0 for trivial or gapless U(1) phase. */
    printf("\n  Three-region Kitaev-Preskill γ construction:\n");
    int A[3] = {0, 1, 2};
    int B[3] = {3, 4, 5};
    int C[3] = {6, 7, 8};

    /* Compute the 7 entropies.  For region R, build an index array and
     * compute S_R via irrep_partial_trace + von Neumann entropy. */
    double S_R[7]; /* 0:A, 1:B, 2:C, 3:A∪B, 4:A∪C, 5:B∪C, 6:A∪B∪C */
    int regions[7][12];
    int sizes[7];
    for (int i = 0; i < 3; ++i) { regions[0][i] = A[i]; regions[1][i] = B[i]; regions[2][i] = C[i]; }
    sizes[0] = sizes[1] = sizes[2] = 3;
    /* A∪B */
    for (int i = 0; i < 3; ++i) { regions[3][i] = A[i]; regions[3][i+3] = B[i]; }
    sizes[3] = 6;
    /* A∪C */
    for (int i = 0; i < 3; ++i) { regions[4][i] = A[i]; regions[4][i+3] = C[i]; }
    sizes[4] = 6;
    /* B∪C */
    for (int i = 0; i < 3; ++i) { regions[5][i] = B[i]; regions[5][i+3] = C[i]; }
    sizes[5] = 6;
    /* A∪B∪C */
    for (int i = 0; i < 3; ++i) {
        regions[6][i] = A[i]; regions[6][i+3] = B[i]; regions[6][i+6] = C[i];
    }
    sizes[6] = 9;

    const char *labels[7] = {"S_A", "S_B", "S_C", "S_{AB}", "S_{AC}", "S_{BC}", "S_{ABC}"};
    printf("  %-8s  %12s\n", "region", "entropy");
    printf("  --------  ------------\n");
    for (int r = 0; r < 7; ++r) {
        long long dim_R = 1LL << sizes[r];
        double _Complex *rho_R = malloc((size_t)dim_R * dim_R * sizeof(double _Complex));
        irrep_partial_trace(N_SITES, 2, psi_full, regions[r], sizes[r], rho_R);
        S_R[r] = irrep_entropy_vonneumann(rho_R, (int)dim_R);
        printf("  %-8s  %+12.8f\n", labels[r], S_R[r]);
        free(rho_R);
    }

    double gamma = S_R[0] + S_R[1] + S_R[2]
                   - S_R[3] - S_R[4] - S_R[5] + S_R[6];
    printf("\n  γ (Kitaev-Preskill) = %+.8f\n", gamma);
    printf("  (log 2 = %+.8f expected for gapped Z_2; 0 for trivial/gapless)\n", log(2.0));
    printf("  (At N=12 this is dominated by finite-size; publication-grade\n");
    printf("   extraction needs matching larger-N runs + extrapolation.)\n");

    /* Spin-flip Z₂ diagnostic: apply F (bit-complement) to the dense GS
     * and compute ⟨ψ|F|ψ⟩. Since F² = I and [H, F] = 0 at zero field,
     * the GS must be an F-eigenstate with eigenvalue ±1. The value
     * selects the Z₂ sector of the GS — a genuine physics property. */
    double _Complex Fexp = 0.0;
    uint64_t mask = ((1ULL << N_SITES) - 1);
    for (long long s = 0; s < D_FULL; ++s)
        Fexp += conj(psi_full[s]) * psi_full[s ^ mask];
    printf("\n  Spin-flip diagnostic: ⟨ψ|F|ψ⟩ = %+.8f %+.8fi\n",
           creal(Fexp), cimag(Fexp));
    if (fabs(creal(Fexp) - 1.0) < 1e-6)
        printf("  → GS is in the F=+1 (even) sector\n");
    else if (fabs(creal(Fexp) + 1.0) < 1e-6)
        printf("  → GS is in the F=-1 (odd) sector\n");
    else
        printf("  → GS is mixed under F (unexpected for zero-field Heisenberg)\n");

    printf("\n  Total wall-clock: %.2f s\n", now_sec() - t0);
    printf("\n  This is the first N in an extrapolation series. The Kitaev-Preskill\n");
    printf("  γ-extraction requires matching larger-N data with three-region\n");
    printf("  topology (γ = S_A + S_B + S_C − S_AB − S_BC − S_AC + S_ABC).\n");

    free(Hdense_full); free(Hblock); free(psi_sector); free(psi_full);
    free(bi); free(bj);
    irrep_sg_little_group_irrep_free(A1);
    irrep_sg_little_group_free(lg);
    irrep_sg_rep_table_free(T);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
