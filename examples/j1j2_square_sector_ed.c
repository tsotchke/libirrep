/* SPDX-License-Identifier: MIT */
/* J₁-J₂ square-lattice Heisenberg on 4×4 (N=16), p4mm, (k, μ_k)-resolved ED.
 *
 * Demonstrates the sparse stack on a different space group (p4mm, not
 * p6mm) with dual bond sets (nearest-neighbour + next-nearest diagonal).
 * The proposed deconfined critical point (DQCP) between Néel antiferromagnet
 * and valence-bond solid phases sits around J₂/J₁ ≈ 0.5, making this the
 * canonical "square-lattice frustrated-magnetism" benchmark.
 *
 * At N=16 with p4mm (|G|=128) the Sz=0 Hilbert space (C(16,8)=12870) reduces
 * to ~100 representatives per sector — small enough that this runs in seconds.
 *
 *   make examples
 *   ./build/bin/j1j2_square_sector_ed
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
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static const char *name_of(irrep_lg_named_irrep_t n) {
    switch (n) {
        case IRREP_LG_IRREP_A1: return "A_1";
        case IRREP_LG_IRREP_A2: return "A_2";
        case IRREP_LG_IRREP_B1: return "B_1";
        case IRREP_LG_IRREP_B2: return "B_2";
        default: return "?";
    }
}

static void do_sector(const irrep_heisenberg_t *H, const irrep_sg_rep_table_t *T,
                      irrep_space_group_t *G, int kx, int ky, const char *k_label,
                      irrep_lg_named_irrep_t name) {
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, kx, ky);
    if (!lg) return;
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_named(lg, name);
    if (!mu) { irrep_sg_little_group_free(lg); return; }

    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    if (!S) {
        irrep_sg_little_group_irrep_free(mu);
        irrep_sg_little_group_free(lg);
        return;
    }
    long long dim = irrep_sg_heisenberg_sector_dim(S);
    if (dim == 0) {
        irrep_sg_heisenberg_sector_free(S);
        irrep_sg_little_group_irrep_free(mu);
        irrep_sg_little_group_free(lg);
        return;
    }

    double _Complex *seed = malloc((size_t)dim * sizeof(double _Complex));
    for (long long i = 0; i < dim; ++i)
        seed[i] = 0.1 * sin(0.41 * i) + I * 0.05 * cos(0.23 * i);
    double eig[4];
    int max_it = (int)((dim > 80) ? 80 : dim);
    int n_eig  = (int)((dim > 4) ? 4 : dim);
    irrep_status_t rc = irrep_lanczos_eigvals_reorth(
        irrep_sg_heisenberg_sector_apply, S, dim, n_eig, max_it, seed, eig);
    if (rc == IRREP_OK)
        printf("  %-6s  %-5s  %6lld  %+11.8f  %+11.8f\n",
               k_label, name_of(name), dim, eig[0], eig[0] / 16.0);
    free(seed);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
}

static void run_at(double J2_over_J1) {
    double t0 = now_sec();
    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P4MM);
    int N  = irrep_space_group_num_sites(G);
    int n_nn  = irrep_lattice_num_bonds_nn(L);
    int n_nnn = irrep_lattice_num_bonds_nnn(L);
    int *nn_i  = malloc(sizeof(int) * n_nn);
    int *nn_j  = malloc(sizeof(int) * n_nn);
    int *nnn_i = malloc(sizeof(int) * n_nnn);
    int *nnn_j = malloc(sizeof(int) * n_nnn);
    irrep_lattice_fill_bonds_nn(L, nn_i, nn_j);
    irrep_lattice_fill_bonds_nnn(L, nnn_i, nnn_j);

    irrep_heisenberg_t *H = irrep_heisenberg_j1j2_new(N, n_nn, nn_i, nn_j, 1.0,
                                                     n_nnn, nnn_i, nnn_j, J2_over_J1);

    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 8); /* Sz = 0 */

    printf("=== J_1-J_2 square, 4×4 (N=%d), p4mm, J_2/J_1 = %.2f ===\n", N, J2_over_J1);
    printf("    %d NN bonds, %d NNN bonds, rep table: %lld reps\n\n",
           n_nn, n_nnn, irrep_sg_rep_table_count(T));
    printf("  %-6s  %-5s  %6s  %11s  %11s\n", "k", "μ", "dim", "E_0 (J)", "E_0/N");
    printf("  ------  -----  ------  -----------  -----------\n");

    /* Γ: full C_4v (A_1, A_2, B_1, B_2; E is 2D — skipped). */
    do_sector(H, T, G, 0, 0, "Γ", IRREP_LG_IRREP_A1);
    do_sector(H, T, G, 0, 0, "Γ", IRREP_LG_IRREP_A2);
    do_sector(H, T, G, 0, 0, "Γ", IRREP_LG_IRREP_B1);
    do_sector(H, T, G, 0, 0, "Γ", IRREP_LG_IRREP_B2);
    /* X at (2, 0) = ½ b₁: C_2v. */
    do_sector(H, T, G, 2, 0, "X",  IRREP_LG_IRREP_A1);
    do_sector(H, T, G, 2, 0, "X",  IRREP_LG_IRREP_A2);
    do_sector(H, T, G, 2, 0, "X",  IRREP_LG_IRREP_B1);
    do_sector(H, T, G, 2, 0, "X",  IRREP_LG_IRREP_B2);
    /* M at (2, 2) = ½ b₁ + ½ b₂: C_4v. */
    do_sector(H, T, G, 2, 2, "M",  IRREP_LG_IRREP_A1);
    do_sector(H, T, G, 2, 2, "M",  IRREP_LG_IRREP_A2);
    do_sector(H, T, G, 2, 2, "M",  IRREP_LG_IRREP_B1);
    do_sector(H, T, G, 2, 2, "M",  IRREP_LG_IRREP_B2);

    printf("\n  Total wall-clock: %.2f s\n\n", now_sec() - t0);

    free(nn_i); free(nn_j); free(nnn_i); free(nnn_j);
    irrep_sg_rep_table_free(T);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
}

int main(void) {
    run_at(0.0);  /* Pure NN Heisenberg — Néel AFM */
    run_at(0.5);  /* DQCP regime — where spin-liquid candidates emerge */
    run_at(1.0);  /* Columnar-VBS regime — strongly frustrated */
    return 0;
}
