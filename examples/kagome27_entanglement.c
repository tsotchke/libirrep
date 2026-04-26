/* SPDX-License-Identifier: MIT */
/* Entanglement extraction on the N=27 kagome-Heisenberg Γ-A_1 ground
 * state. Beyond dense-diagonalisation reach (sector dim = 186 616);
 * uses sparse Lanczos with reorth + explicit eigenvector recovery
 * (lanczos_eigvecs_reorth), then unfolds to a 2^27 ≈ 134 million-dim
 * dense state for partial trace.
 *
 * Memory budget: dense state vector = 2.0 GB. Partial trace on |A| ≤ 10
 * sites uses an RDM matrix up to 2^10 × 2^10 = 1024² complex = 16 MB.
 * Workstation-scale but not laptop-scale.
 *
 * Physical prediction (exact): at Sz = p/N (here p = 13, N = 27), the
 * single-site RDM is diagonal with (p/N, 1−p/N) by translation symmetry,
 * so S_{|A|=1} = −(p/N)·ln(p/N) − (1−p/N)·ln(1−p/N). For kagome 3×3
 * popcount=13 that's −(13/27)·ln(13/27) − (14/27)·ln(14/27) ≈ 0.69246.
 * Matching this pins the full pipeline correctness at N = 27. At even
 * N with Sz = 0 (e.g., kagome12_entanglement), the prediction collapses
 * to S_1 = ln 2 exactly via SU(2) symmetry of the singlet GS.
 *
 * Wall-clock on Apple M2 Ultra (USE_OPENMP=1):
 *   rep table      ~4 s
 *   sector build   ~13 s
 *   Lanczos (60 iters, 1 eigenvector)  ~90 s
 *   unfold to 2^27 ~20 s
 *   partial trace per |A|  ~30 s
 * Total: ~3-4 minutes.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome27_entanglement
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

#define N_SITES 27
#define POPCOUNT 13 /* Sz = -1/2 (odd N; closest-to-zero integer Sz) */

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(void) {
    printf("=== Kagome N=27 Γ-A_1 ground-state entanglement ===\n\n");
    double t0 = now_sec();

    /* 1. Setup. */
    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);

    /* 2. Rep table + Γ-A_1 sector binding. */
    double t = now_sec();
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);
    printf("  rep table (popcount %d): %lld reps in %.2f s\n",
           POPCOUNT, irrep_sg_rep_table_count(T), now_sec() - t);

    t = now_sec();
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
    irrep_sg_little_group_irrep_t *A1 =
        irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_A1);
    irrep_sg_heisenberg_sector_t *S =
        irrep_sg_heisenberg_sector_build_at_k(H, T, lg, A1);
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);
    printf("  sector (Γ, A_1) built: dim = %lld in %.2f s\n", sector_dim, now_sec() - t);

    /* 3. Sparse Lanczos with eigenvector recovery. */
    t = now_sec();
    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.13 * sin(0.41 * (double)i) + I * 0.07 * cos(0.23 * (double)i);
    double eigs[4];
    double _Complex *psi_sector = malloc((size_t)sector_dim * sizeof(double _Complex));
    /* 150 iterations — eigenvector accuracy is more demanding than
     * eigenvalue accuracy (Ritz-vector residual scales ~ (ratio_gap)^k).
     * For SU(2) check at 1e-8 precision, 60 is insufficient. */
    irrep_status_t rc = irrep_lanczos_eigvecs_reorth(
        irrep_sg_heisenberg_sector_apply, S, sector_dim, 1, 150, seed, eigs, psi_sector);
    if (rc != IRREP_OK) {
        fprintf(stderr, "Lanczos eigvec extraction failed: rc = %d\n", rc);
        return 1;
    }
    printf("  sparse Lanczos (60 iters, 1 eigvec): %.2f s\n", now_sec() - t);
    printf("    E_0 / J = %+.8f   (E_0 / N = %+.8f)\n", eigs[0], eigs[0] / (double)N_SITES);
    /* Verify residual ||H·ψ − E·ψ||. */
    double _Complex *Hpsi = malloc((size_t)sector_dim * sizeof(double _Complex));
    irrep_sg_heisenberg_sector_apply(psi_sector, Hpsi, S);
    double res2 = 0.0;
    for (long long i = 0; i < sector_dim; ++i) {
        double _Complex r = Hpsi[i] - eigs[0] * psi_sector[i];
        res2 += creal(r * conj(r));
    }
    printf("    ||H·ψ − E·ψ|| = %.3e (should be small)\n", sqrt(res2));
    free(Hpsi); free(seed);

    /* 4. Unfold sector → dense 2^N. Construct |ũ_k⟩ for each rep and
     * accumulate psi_sector[k] · |ũ_k⟩. */
    long long D = 1LL << N_SITES;
    printf("\n  allocating dense state vector: %.2f GB\n",
           (double)D * 16.0 / (1024.0 * 1024.0 * 1024.0));
    double _Complex *psi_full = calloc((size_t)D, sizeof(double _Complex));
    if (!psi_full) {
        fprintf(stderr, "OOM on dense state vector\n");
        return 1;
    }

    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, A1, w);

    /* Unfold efficiently: for each rep, only |G| = 108 orbit sites receive
     * non-zero amplitude. Compact storage (pos, amp) arrays of that size,
     * not a full 2^N scratch. Order-|G| work per rep, not order-2^N. */
    t = now_sec();
    uint64_t        *orbit_pos = malloc((size_t)order * sizeof(uint64_t));
    double _Complex *orbit_amp = malloc((size_t)order * sizeof(double _Complex));
    long long n_reps = irrep_sg_rep_table_count(T);
    long long written = 0;
    for (long long k = 0; k < n_reps && written < sector_dim; ++k) {
        uint64_t u = irrep_sg_rep_table_get(T, k);
        /* Compute distinct (g·u, accumulated-weight) pairs. Since distinct
         * g may map to the same g·u (Stab elements), we aggregate. */
        int n_entries = 0;
        for (int g = 0; g < order; ++g) {
            uint64_t gx = irrep_space_group_apply_bits(G, g, u);
            /* Linear scan for existing entry (order fits in O(100)). */
            int found = -1;
            for (int e = 0; e < n_entries; ++e) {
                if (orbit_pos[e] == gx) { found = e; break; }
            }
            if (found >= 0) orbit_amp[found] += w[g];
            else { orbit_pos[n_entries] = gx; orbit_amp[n_entries] = w[g]; ++n_entries; }
        }
        /* Norm². */
        double nn = 0.0;
        for (int e = 0; e < n_entries; ++e)
            nn += creal(orbit_amp[e] * conj(orbit_amp[e]));
        if (nn < 1e-20) continue;
        double inv = 1.0 / sqrt(nn);
        double _Complex coef = psi_sector[written];
        /* Accumulate into psi_full. */
        for (int e = 0; e < n_entries; ++e)
            psi_full[orbit_pos[e]] += coef * inv * orbit_amp[e];
        ++written;
    }
    printf("  unfold to dense 2^%d: %.2f s\n", N_SITES, now_sec() - t);
    free(orbit_pos); free(orbit_amp); free(w); free(psi_sector);

    /* Normalisation check. */
    double full_norm = 0.0;
    for (long long s = 0; s < D; ++s)
        full_norm += creal(psi_full[s] * conj(psi_full[s]));
    printf("    ||ψ_full||² = %.10f\n", full_norm);

    /* 5. Entanglement entropies at small |A|. */
    printf("\n  Entanglement entropy S_A vs subsystem size:\n");
    printf("  %3s   %12s   %s\n", "|A|", "S_A", "comment");
    printf("  ---   ------------   -------\n");
    for (int nA = 1; nA <= 6; ++nA) {
        t = now_sec();
        int sites_A[27];
        for (int i = 0; i < nA; ++i) sites_A[i] = i;
        long long dim_A = 1LL << nA;
        double _Complex *rho_A = malloc((size_t)dim_A * dim_A * sizeof(double _Complex));
        irrep_partial_trace(N_SITES, 2, psi_full, sites_A, nA, rho_A);
        double S = irrep_entropy_vonneumann(rho_A, (int)dim_A);
        const char *comment = "";
        if (nA == 1) {
            /* Exact Sz-resolved prediction for single-site S_1 on popcount p. */
            double p_up = (double)POPCOUNT / (double)N_SITES;
            double p_dn = 1.0 - p_up;
            double S1_exact = -p_up * log(p_up) - p_dn * log(p_dn);
            double err = fabs(S - S1_exact);
            if (err < 1e-8)      comment = "matches exact Sz-resolved prediction";
            else if (err < 1e-4) comment = "≈ exact Sz prediction";
        }
        printf("  %3d   %+12.8f   %s  [%.1f s]\n", nA, S, comment, now_sec() - t);
        free(rho_A);
    }

    /* Kitaev-Preskill three-region γ. |A|=|B|=|C|=3 matches the N=12
     * setup (same region size, larger environment). S_{ABC} RDM is
     * 2^9 = 512, keeping partial_trace tractable at N=27
     * (dB·dA² = 2^18·2^18 = 68 G ops ≈ ~60 s). */
    printf("\n  Three-region Kitaev-Preskill γ construction (|A|=|B|=|C|=3):\n");
    int A[3] = {0, 1, 2};
    int B[3] = {9, 10, 11};
    int C[3] = {18, 19, 20};

    int regions[7][27];
    int sizes[7];
    for (int i = 0; i < 3; ++i) { regions[0][i] = A[i]; regions[1][i] = B[i]; regions[2][i] = C[i]; }
    sizes[0] = sizes[1] = sizes[2] = 3;
    for (int i = 0; i < 3; ++i) { regions[3][i] = A[i]; regions[3][i+3] = B[i]; }
    sizes[3] = 6;
    for (int i = 0; i < 3; ++i) { regions[4][i] = A[i]; regions[4][i+3] = C[i]; }
    sizes[4] = 6;
    for (int i = 0; i < 3; ++i) { regions[5][i] = B[i]; regions[5][i+3] = C[i]; }
    sizes[5] = 6;
    for (int i = 0; i < 3; ++i) {
        regions[6][i] = A[i]; regions[6][i+3] = B[i]; regions[6][i+6] = C[i];
    }
    sizes[6] = 9;

    const char *labels[7] = {"S_A", "S_B", "S_C", "S_{AB}", "S_{AC}", "S_{BC}", "S_{ABC}"};
    double S_R[7];
    printf("  %-8s  %12s  %s\n", "region", "entropy", "[t]");
    printf("  --------  ------------  -----\n");
    for (int r = 0; r < 7; ++r) {
        long long dim_R = 1LL << sizes[r];
        double _Complex *rho_R = malloc((size_t)dim_R * dim_R * sizeof(double _Complex));
        if (!rho_R) { fprintf(stderr, "OOM on region %d RDM\n", r); return 1; }
        double tr = now_sec();
        irrep_partial_trace(N_SITES, 2, psi_full, regions[r], sizes[r], rho_R);
        S_R[r] = irrep_entropy_vonneumann(rho_R, (int)dim_R);
        printf("  %-8s  %+12.8f  [%.1f s]\n", labels[r], S_R[r], now_sec() - tr);
        free(rho_R);
    }

    double gamma = S_R[0] + S_R[1] + S_R[2]
                   - S_R[3] - S_R[4] - S_R[5] + S_R[6];
    printf("\n  γ (Kitaev-Preskill) = %+.8f\n", gamma);
    printf("  (log 2 = %+.8f expected for gapped Z_2; 0 for trivial/gapless)\n", log(2.0));

    printf("\n  Total wall-clock: %.2f s\n", now_sec() - t0);
    printf("\n  N = 27 provides the second data point in an N-scaling series\n");
    printf("  (prior: N=12 in kagome12_entanglement). Kitaev-Preskill γ\n");
    printf("  needs three-region geometry + extrapolation to ∞.\n");

    free(psi_full);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(A1);
    irrep_sg_little_group_free(lg);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
