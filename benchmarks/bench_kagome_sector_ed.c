/* SPDX-License-Identifier: MIT */
/* Benchmark: symmetry-resolved sparse Lanczos on the 3×3 kagome torus
 * (N = 27, p6mm, popcount = 13 = Sz −1/2). Exercises the full sparse
 * stack — rep-table enumeration, (k, μ_k) sector build, Lanczos on the
 * cached CSR — and emits JSON records per phase for regression tracking.
 *
 * The dense path is infeasible at this size (2 GiB per state vector,
 * ~120 GiB for reorth Lanczos at 60 iters), so this harness is the
 * primary performance gate for the sparse-stack hot path. */
#include "harness.h"

#include <irrep/config_project.h>
#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>
#include <irrep/space_group.h>

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void emit_phase(const char *name, double ns) {
    printf("{\"name\":\"%s\",\"iterations\":1,\"ops_per_iter\":1,"
           "\"ns_per_op\":%.0f,\"ops_per_sec\":%.3e}\n",
           name, ns, 1e9 / ns);
}

int main(void) {
    /* Build the lattice + Hamiltonian once. */
    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
    int N  = irrep_space_group_num_sites(G);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N, nb, bi, bj, 1.0);

    /* Phase 1: rep-table enumeration at Sz = −1/2. */
    double t = irrep_bench_now_ns();
    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 13);
    double t_reps = irrep_bench_now_ns() - t;
    emit_phase("kagome27_sz_half/rep_table_build", t_reps);

    /* Phase 2: sector build — trivial (Γ, A_1). */
    t = irrep_bench_now_ns();
    irrep_sg_heisenberg_sector_t *S_triv = irrep_sg_heisenberg_sector_build(H, T);
    double t_triv_build = irrep_bench_now_ns() - t;
    emit_phase("kagome27_sz_half/sector_build_trivial", t_triv_build);

    /* Phase 3: Lanczos on trivial sector. */
    long long dim_triv = irrep_sg_rep_table_count(T);
    double _Complex *seed = malloc((size_t)dim_triv * sizeof(double _Complex));
    for (long long i = 0; i < dim_triv; ++i)
        seed[i] = 0.13 * sin(0.41 * (double)i) + I * 0.07 * cos(0.23 * (double)i);
    double eig[4];

    t = irrep_bench_now_ns();
    irrep_lanczos_eigvals_reorth(irrep_sg_heisenberg_sector_apply, S_triv,
                                 dim_triv, 4, 60, seed, eig);
    double t_triv_lanc = irrep_bench_now_ns() - t;
    emit_phase("kagome27_sz_half/lanczos_trivial_60iter", t_triv_lanc);

    /* Sanity-check the ground state reproduces the example's value. */
    if (fabs(eig[0] - (-11.60977827)) > 1e-3) {
        fprintf(stderr, "bench_kagome_sector_ed: Γ-A_1 E_0 drift: got %.8f\n", eig[0]);
    }
    irrep_sg_heisenberg_sector_free(S_triv);

    /* Phase 4: sector build — (K, A_1), full (k, μ_k) path. */
    irrep_sg_little_group_t *lg_K = irrep_sg_little_group_build(G, 1, 2);
    irrep_sg_little_group_irrep_t *K_A1 =
        irrep_sg_little_group_irrep_named(lg_K, IRREP_LG_IRREP_A1);

    t = irrep_bench_now_ns();
    irrep_sg_heisenberg_sector_t *S_K =
        irrep_sg_heisenberg_sector_build_at_k(H, T, lg_K, K_A1);
    double t_kA1_build = irrep_bench_now_ns() - t;
    emit_phase("kagome27_sz_half/sector_build_K_A1", t_kA1_build);

    /* Phase 5: Lanczos on (K, A_1). */
    long long dim_K = irrep_sg_heisenberg_sector_dim(S_K);
    double _Complex *seed_K = malloc((size_t)dim_K * sizeof(double _Complex));
    for (long long i = 0; i < dim_K; ++i)
        seed_K[i] = 0.11 * sin(0.37 * (double)i) + I * 0.09 * cos(0.19 * (double)i);

    t = irrep_bench_now_ns();
    irrep_lanczos_eigvals_reorth(irrep_sg_heisenberg_sector_apply, S_K,
                                 dim_K, 4, 60, seed_K, eig);
    double t_kA1_lanc = irrep_bench_now_ns() - t;
    emit_phase("kagome27_sz_half/lanczos_K_A1_60iter", t_kA1_lanc);

    /* The K-point sector carries the Dirac-cone physics — E_0 should be
     * substantially below the Γ-A_1 result (~-16.55 J on 3×3 kagome). */
    if (eig[0] > -15.0) {
        fprintf(stderr, "bench_kagome_sector_ed: K-A_1 E_0 unexpected: %.8f\n", eig[0]);
    }

    free(seed_K);
    irrep_sg_heisenberg_sector_free(S_K);
    irrep_sg_little_group_irrep_free(K_A1);
    irrep_sg_little_group_free(lg_K);

    free(seed);
    free(bi);
    free(bj);
    irrep_sg_rep_table_free(T);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
