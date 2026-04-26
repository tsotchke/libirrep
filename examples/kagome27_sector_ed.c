/* SPDX-License-Identifier: MIT */
/* N = 27 kagome-Heisenberg symmetry-resolved exact diagonalisation.
 *
 * Enumerates every 1D (k, μ_k) sector on the 3×3 kagome cluster (N = 27,
 * popcount = 13, Sz = −1/2) and runs sparse Lanczos in each. The full
 * Hilbert dimension at this popcount is C(27, 13) = 20 058 300; the
 * sparse path reduces this by ~|G| = 108 for an average per-sector
 * dimension around ~186 k.
 *
 * This is the primitive kagome-torus ED benchmark at N = 27:
 *   - Γ-point: C_6v little group — 4 valid 1D irreps (A_1, A_2, B_1, B_2).
 *     E_1 and E_2 are 2D and skipped here.
 *   - K-points (1/3, 2/3) and (2/3, 1/3): C_3v — 2 valid 1D irreps
 *     (A_1, A_2). E is 2D and skipped here.
 *   - 6 generic low-symmetry k-points: little group is order 2 (identity
 *     + one mirror) → 2 irreps each.
 *
 * Dense path is infeasible: one state vector at N = 27 is 2.0 GiB; reorth
 * Lanczos at 60 iters needs ~120 GiB. Sparse path runs in seconds per
 * sector on a workstation.
 *
 *   make examples
 *   ./build/bin/kagome27_sector_ed
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

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static const char *name_of(irrep_lg_named_irrep_t n) {
    switch (n) {
        case IRREP_LG_IRREP_A1: return "A_1";
        case IRREP_LG_IRREP_A2: return "A_2";
        case IRREP_LG_IRREP_B1: return "B_1";
        case IRREP_LG_IRREP_B2: return "B_2";
        case IRREP_LG_IRREP_E:  return "E";
        case IRREP_LG_IRREP_E1: return "E_1";
        case IRREP_LG_IRREP_E2: return "E_2";
        default: return "?";
    }
}

static void do_sector(const irrep_heisenberg_t *H, const irrep_sg_rep_table_t *T,
                      irrep_space_group_t *G, int kx, int ky, const char *k_label,
                      irrep_lg_named_irrep_t name) {
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, kx, ky);
    if (!lg) return;
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_named(lg, name);
    if (!mu) {
        /* Not valid for this little group (e.g. B_1 on C_3v). Silent skip. */
        irrep_sg_little_group_free(lg);
        return;
    }

    double tb = now_sec();
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    double t_build = now_sec() - tb;
    if (!S) {
        printf("  %-10s  %-5s  (build failed — likely 2D irrep)\n", k_label, name_of(name));
        irrep_sg_little_group_irrep_free(mu);
        irrep_sg_little_group_free(lg);
        return;
    }
    long long dim = irrep_sg_heisenberg_sector_dim(S);
    if (dim == 0) {
        printf("  %-10s  %-5s  dim=0 (annihilated)\n", k_label, name_of(name));
        irrep_sg_heisenberg_sector_free(S);
        irrep_sg_little_group_irrep_free(mu);
        irrep_sg_little_group_free(lg);
        return;
    }

    double _Complex *seed = malloc((size_t)dim * sizeof(double _Complex));
    for (long long i = 0; i < dim; ++i)
        seed[i] = 0.13 * sin(0.41 * i + 0.1) + I * 0.07 * cos(0.23 * i);
    double eig[4];
    int max_it = (int)((dim > 60) ? 60 : dim);
    int n_eig  = (int)((dim > 4) ? 4 : dim);

    double tl = now_sec();
    irrep_status_t rc = irrep_lanczos_eigvals_reorth(
        irrep_sg_heisenberg_sector_apply, S, dim, n_eig, max_it, seed, eig);
    double t_lanc = now_sec() - tl;

    if (rc != IRREP_OK)
        printf("  %-10s  %-5s  %9lld  (Lanczos rc=%d)\n", k_label, name_of(name), dim, rc);
    else
        printf("  %-10s  %-5s  %9lld  %+11.8f  %+11.8f  build %6.2fs  lanczos %5.2fs\n",
               k_label, name_of(name), dim, eig[0], eig[0] / 27.0, t_build, t_lanc);

    free(seed);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
}

int main(void) {
    printf("=== N=27 kagome Heisenberg: (k, μ_k)-resolved sparse ED ===\n");
    printf("    3×3 kagome torus, popcount = 13 (Sz = −1/2, odd N)\n\n");

    double t0 = now_sec();
    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
    int N  = irrep_space_group_num_sites(G);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N, nb, bi, bj, 1.0);

    printf("    N = %d sites, %d NN bonds, p6mm order = %d\n",
           N, nb, irrep_space_group_order(G));

    double tr = now_sec();
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 13);
    double t_reps = now_sec() - tr;
    printf("    rep table: %lld reps built in %.2f s\n", irrep_sg_rep_table_count(T), t_reps);
    printf("    (C(27, 13) = 20 058 300; reduction ≈ %.1fx)\n\n",
           20058300.0 / (double)irrep_sg_rep_table_count(T));

    printf("  %-10s  %-5s  %9s  %11s  %11s\n", "k", "μ", "dim", "E_0 (J)", "E_0/N");
    printf("  ----------  -----  ---------  -----------  -----------\n");

    /* Γ-point: C_6v, 1D irreps A_1, A_2, B_1, B_2. */
    do_sector(H, T, G, 0, 0, "Γ", IRREP_LG_IRREP_A1);
    do_sector(H, T, G, 0, 0, "Γ", IRREP_LG_IRREP_A2);
    do_sector(H, T, G, 0, 0, "Γ", IRREP_LG_IRREP_B1);
    do_sector(H, T, G, 0, 0, "Γ", IRREP_LG_IRREP_B2);

    /* K-points: C_3v, A_1 and A_2 are 1D. */
    do_sector(H, T, G, 1, 2, "K",  IRREP_LG_IRREP_A1);
    do_sector(H, T, G, 1, 2, "K",  IRREP_LG_IRREP_A2);
    do_sector(H, T, G, 2, 1, "K'", IRREP_LG_IRREP_A1);
    do_sector(H, T, G, 2, 1, "K'", IRREP_LG_IRREP_A2);

    double t_total = now_sec() - t0;
    printf("\n  Total wall-clock: %.2f s (incl. rep-table build %.2f s)\n",
           t_total, t_reps);

    free(bi); free(bj);
    irrep_sg_rep_table_free(T);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
