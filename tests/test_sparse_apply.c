/* SPDX-License-Identifier: MIT */
/* Cross-check the sparse rep-basis Heisenberg apply against a dense
 * construction: for each representative u_k in the rep table, build
 * |ũ_k⟩ = (1/√orbit) · Σ_{x ∈ orbit} |x⟩ as a full-Hilbert vector, apply
 * the existing dense Heisenberg apply, project back onto each rep basis
 * vector. The result must match the sparse apply column-for-column
 * (within floating-point associativity tolerance).
 *
 * Tested on: 4-site Heisenberg ring (p1), kagome 2×2 Heisenberg (p6mm,
 * trivial sector = Γ, A_1) at multiple popcount sectors.
 *
 * This is the correctness gate for the sparse path to N ≥ 27 ED. */
#include "harness.h"
#include <irrep/config_project.h>
#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>
#include <irrep/space_group.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Materialise |ũ_k⟩ into a dense 2^N complex vector. */
static void build_dense_basis_vec(const irrep_space_group_t *G, uint64_t rep,
                                  int orbit_size, int num_sites,
                                  double _Complex *dense_out, long long D) {
    memset(dense_out, 0, (size_t)D * sizeof(double _Complex));
    int order = irrep_space_group_order(G);
    double inv_sqrt = 1.0 / sqrt((double)orbit_size);
    /* Walk every group element; mark each distinct image. Duplicates from
     * stabiliser elements produce the correct weighted sum because we
     * assign 0 → (1/√orbit) on first encounter and leave it untouched
     * (computational basis is already orthonormal). */
    (void)num_sites;
    for (int g = 0; g < order; ++g) {
        uint64_t x = irrep_space_group_apply_bits(G, g, rep);
        if (dense_out[(long long)x] == 0.0)
            dense_out[(long long)x] = inv_sqrt;
    }
}

int main(void) {
    IRREP_TEST_START("sparse_apply");

    /* ---- 4-site Heisenberg ring, p1 (translations only) ---------------- */
    {
        int N = 4;
        int bi[] = {0, 1, 2, 3}, bj[] = {1, 2, 3, 0};
        double J = 1.0;
        irrep_heisenberg_t *H = irrep_heisenberg_new(N, 4, bi, bj, J);

        /* Lattice + p1 rep table at popcount 2 (Sz = 0). */
        irrep_lattice_t     *L  = irrep_lattice_build(IRREP_LATTICE_SQUARE, 2, 2);
        irrep_space_group_t *G  = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
        irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 2);
        IRREP_ASSERT(T != NULL);
        long long n = irrep_sg_rep_table_count(T);
        IRREP_ASSERT(n > 0);

        long long D = 1LL << N;
        double _Complex *dense_vec = malloc((size_t)D * sizeof(double _Complex));
        double _Complex *Hb_dense  = malloc((size_t)D * sizeof(double _Complex));
        double _Complex *sparse_in  = malloc((size_t)n * sizeof(double _Complex));
        double _Complex *sparse_out = malloc((size_t)n * sizeof(double _Complex));

        /* For each column k: sparse y = H |ũ_k⟩; dense y via apply+project. */
        for (long long k = 0; k < n; ++k) {
            /* Sparse column. */
            for (long long j = 0; j < n; ++j)
                sparse_in[j] = (j == k) ? 1.0 : 0.0;
            irrep_heisenberg_apply_in_sector(H, T, sparse_in, sparse_out);

            /* Dense column: build |ũ_k⟩, apply dense H, project to |ũ_j⟩. */
            uint64_t rep_k = irrep_sg_rep_table_get(T, k);
            int      o_k   = irrep_sg_rep_table_orbit_size(T, k);
            build_dense_basis_vec(G, rep_k, o_k, N, dense_vec, D);
            irrep_heisenberg_apply(dense_vec, Hb_dense, H);
            /* Project: (H_dense)_{j,k} = ⟨ũ_j|H|ũ_k⟩. */
            for (long long j = 0; j < n; ++j) {
                uint64_t rep_j = irrep_sg_rep_table_get(T, j);
                int      o_j   = irrep_sg_rep_table_orbit_size(T, j);
                double _Complex *bas_j = malloc((size_t)D * sizeof(double _Complex));
                build_dense_basis_vec(G, rep_j, o_j, N, bas_j, D);
                double _Complex proj = 0.0;
                for (long long s = 0; s < D; ++s)
                    proj += conj(bas_j[s]) * Hb_dense[s];
                double diff = cabs(proj - sparse_out[j]);
                IRREP_ASSERT(diff < 1e-13);
                free(bas_j);
            }
        }

        free(dense_vec); free(Hb_dense); free(sparse_in); free(sparse_out);
        irrep_sg_rep_table_free(T);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
        irrep_heisenberg_free(H);
    }

    /* ---- Kagome 2×2 Heisenberg, p6mm (Γ, A_1) sector, Sz = 0 ----------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int N  = irrep_space_group_num_sites(G); /* 12 */
        int nb = irrep_lattice_num_bonds_nn(L);
        int *bi = malloc(sizeof(int) * nb);
        int *bj = malloc(sizeof(int) * nb);
        irrep_lattice_fill_bonds_nn(L, bi, bj);

        irrep_heisenberg_t   *H = irrep_heisenberg_new(N, nb, bi, bj, 1.0);
        irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 6);
        long long n = irrep_sg_rep_table_count(T);
        IRREP_ASSERT(n == 30); /* matches test_canonicalise */

        long long D = 1LL << N; /* 4096 */
        double _Complex *dense_vec = malloc((size_t)D * sizeof(double _Complex));
        double _Complex *Hb_dense  = malloc((size_t)D * sizeof(double _Complex));
        double _Complex *sparse_in  = malloc((size_t)n * sizeof(double _Complex));
        double _Complex *sparse_out = malloc((size_t)n * sizeof(double _Complex));
        double _Complex *bas_j      = malloc((size_t)D * sizeof(double _Complex));

        double max_err = 0.0;
        for (long long k = 0; k < n; ++k) {
            for (long long j = 0; j < n; ++j)
                sparse_in[j] = (j == k) ? 1.0 : 0.0;
            irrep_heisenberg_apply_in_sector(H, T, sparse_in, sparse_out);

            uint64_t rep_k = irrep_sg_rep_table_get(T, k);
            int      o_k   = irrep_sg_rep_table_orbit_size(T, k);
            build_dense_basis_vec(G, rep_k, o_k, N, dense_vec, D);
            irrep_heisenberg_apply(dense_vec, Hb_dense, H);

            for (long long j = 0; j < n; ++j) {
                uint64_t rep_j = irrep_sg_rep_table_get(T, j);
                int      o_j   = irrep_sg_rep_table_orbit_size(T, j);
                build_dense_basis_vec(G, rep_j, o_j, N, bas_j, D);
                double _Complex proj = 0.0;
                for (long long s = 0; s < D; ++s)
                    proj += conj(bas_j[s]) * Hb_dense[s];
                double diff = cabs(proj - sparse_out[j]);
                if (diff > max_err)
                    max_err = diff;
                IRREP_ASSERT(diff < 1e-12);
            }
        }
        (void)max_err;

        free(bas_j); free(dense_vec); free(Hb_dense);
        free(sparse_in); free(sparse_out);
        free(bi); free(bj);
        irrep_sg_rep_table_free(T);
        irrep_heisenberg_free(H);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Full (k, μ_k) sparse apply: Γ-B_1 on kagome 2×2 ---------------
     * The absolute KH ground state on kagome 2×2 sits at (Γ, B_1) with
     * E_0 = −5.44487522 J. The full-sector sparse apply via
     * irrep_sg_heisenberg_sector_build_at_k must reproduce this. */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int N  = irrep_space_group_num_sites(G);
        int nb = irrep_lattice_num_bonds_nn(L);
        int *bi = malloc(sizeof(int) * nb);
        int *bj = malloc(sizeof(int) * nb);
        irrep_lattice_fill_bonds_nn(L, bi, bj);
        irrep_heisenberg_t *H = irrep_heisenberg_new(N, nb, bi, bj, 1.0);

        irrep_sg_rep_table_t          *T  = irrep_sg_rep_table_build(G, 6);
        irrep_sg_little_group_t       *lg = irrep_sg_little_group_build(G, 0, 0);
        irrep_sg_little_group_irrep_t *B1 =
            irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_B1);
        IRREP_ASSERT(B1 != NULL);

        irrep_sg_heisenberg_sector_t *S_b1 =
            irrep_sg_heisenberg_sector_build_at_k(H, T, lg, B1);
        IRREP_ASSERT(S_b1 != NULL);
        long long dim_b1 = irrep_sg_heisenberg_sector_dim(S_b1);
        IRREP_ASSERT(dim_b1 > 0);
        IRREP_ASSERT(dim_b1 <= irrep_sg_rep_table_count(T));

        /* Lanczos → E_0 in Γ-B_1. */
        double _Complex *seed = malloc((size_t)dim_b1 * sizeof(double _Complex));
        for (long long i = 0; i < dim_b1; ++i)
            seed[i] = 0.13 * sin(0.41 * i) + I * 0.07 * cos(0.23 * i);
        double eig[4];
        int max_it = (int)((dim_b1 > 40) ? 40 : dim_b1);
        int n_eig  = (int)((dim_b1 > 4) ? 4 : dim_b1);
        irrep_status_t rc = irrep_lanczos_eigvals_reorth(
            irrep_sg_heisenberg_sector_apply, S_b1, dim_b1, n_eig, max_it, seed, eig);
        IRREP_ASSERT(rc == IRREP_OK);
        IRREP_ASSERT(fabs(eig[0] - (-5.44487522)) < 1e-5);

        free(seed);
        irrep_sg_heisenberg_sector_free(S_b1);

        /* Sanity: A_1 sector same dim as trivial-sector rep table? Not
         * necessarily — σ_u can be zero when char sum annihilates — but
         * for trivial 1D irrep (all +1s) σ_u = |Stab(u)| > 0 always. So
         * A_1 dim == rep_table count. */
        irrep_sg_little_group_irrep_t *A1 =
            irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_A1);
        irrep_sg_heisenberg_sector_t *S_a1 =
            irrep_sg_heisenberg_sector_build_at_k(H, T, lg, A1);
        IRREP_ASSERT(S_a1 != NULL);
        IRREP_ASSERT(irrep_sg_heisenberg_sector_dim(S_a1) ==
                     irrep_sg_rep_table_count(T));

        /* A_1 Lanczos: E_0 matches the trivial-sector sparse apply at −5.32839. */
        double _Complex *seed_a1 = malloc((size_t)irrep_sg_rep_table_count(T) *
                                          sizeof(double _Complex));
        for (long long i = 0; i < irrep_sg_rep_table_count(T); ++i)
            seed_a1[i] = 0.1 * sin(0.37 * i);
        double eig_a1[4];
        rc = irrep_lanczos_eigvals_reorth(irrep_sg_heisenberg_sector_apply, S_a1,
                                          irrep_sg_rep_table_count(T), 4, 30, seed_a1, eig_a1);
        IRREP_ASSERT(rc == IRREP_OK);
        IRREP_ASSERT(fabs(eig_a1[0] - (-5.32839240)) < 1e-5);

        free(seed_a1);
        irrep_sg_heisenberg_sector_free(S_a1);

        /* 2D irreps: sparse build goes through dense-under-the-hood.
         * Consumer API unchanged; Lanczos on the CSR output must reproduce
         * the known Γ-E_1 ground state (-3.299516 J on kagome 12). */
        irrep_sg_little_group_irrep_t *E1 =
            irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_E1);
        irrep_sg_heisenberg_sector_t *S_E1_sp =
            irrep_sg_heisenberg_sector_build_at_k(H, T, lg, E1);
        IRREP_ASSERT(S_E1_sp != NULL);
        long long dim_E1_sp = irrep_sg_heisenberg_sector_dim(S_E1_sp);
        IRREP_ASSERT(dim_E1_sp == 20);

        double _Complex *seed_E1 = malloc((size_t)dim_E1_sp * sizeof(double _Complex));
        for (long long i = 0; i < dim_E1_sp; ++i)
            seed_E1[i] = 0.11 * sin(0.31 * i) + I * 0.07 * cos(0.19 * i);
        double eig_E1[4];
        int ni_E1 = (int)((dim_E1_sp > 40) ? 40 : dim_E1_sp);
        int ne_E1 = (int)((dim_E1_sp > 4) ? 4 : dim_E1_sp);
        irrep_status_t rc_E1 = irrep_lanczos_eigvals_reorth(
            irrep_sg_heisenberg_sector_apply, S_E1_sp, dim_E1_sp, ne_E1, ni_E1,
            seed_E1, eig_E1);
        IRREP_ASSERT(rc_E1 == IRREP_OK);
        IRREP_ASSERT(fabs(eig_E1[0] - (-3.299516)) < 1e-5);

        free(seed_E1);
        irrep_sg_heisenberg_sector_free(S_E1_sp);
        irrep_sg_little_group_irrep_free(E1);

        /* ---- Dense reference cross-check: sparse Γ-B_1 energies must
         * match the dense-projector construction bit-exactly (to 1e-10
         * Lanczos precision). This validates the dense reference; once
         * validated it becomes the ground-truth oracle for the in-progress
         * 2D sparse formula work. */
        irrep_sg_little_group_irrep_t *B1_ref =
            irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_B1);
        irrep_sg_heisenberg_sector_t *S_B1_sp =
            irrep_sg_heisenberg_sector_build_at_k(H, T, lg, B1_ref);
        long long dim_B1_sp = irrep_sg_heisenberg_sector_dim(S_B1_sp);

        /* Sparse Γ-B_1 E_0. */
        double _Complex *seed_sp = malloc((size_t)dim_B1_sp * sizeof(double _Complex));
        for (long long i = 0; i < dim_B1_sp; ++i)
            seed_sp[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
        double eig_sp[4];
        int ni_sp = (int)((dim_B1_sp > 40) ? 40 : dim_B1_sp);
        int ne_sp = (int)((dim_B1_sp > 4) ? 4 : dim_B1_sp);
        irrep_lanczos_eigvals_reorth(irrep_sg_heisenberg_sector_apply, S_B1_sp,
                                     dim_B1_sp, ne_sp, ni_sp, seed_sp, eig_sp);

        /* Dense Γ-B_1: build matrix, diagonalise. */
        long long D_full = 1LL << N;
        (void)D_full;
        long long cap_ref = irrep_sg_rep_table_count(T);
        double _Complex *Hmat_ref = malloc((size_t)cap_ref * cap_ref * sizeof(double _Complex));
        int dim_B1_dense = irrep_sg_heisenberg_sector_build_dense(H, T, lg, B1_ref,
                                                                   Hmat_ref, (int)cap_ref);
        IRREP_ASSERT(dim_B1_dense > 0);
        IRREP_ASSERT(dim_B1_dense == dim_B1_sp);

        /* Extract the dim_B1_dense × dim_B1_dense block from the top-left
         * of Hmat_ref and diagonalise. (cap_ref is the row stride.) */
        double _Complex *Hblock = malloc((size_t)dim_B1_dense * dim_B1_dense * sizeof(double _Complex));
        for (int r = 0; r < dim_B1_dense; ++r)
            for (int c = 0; c < dim_B1_dense; ++c)
                Hblock[(size_t)r * dim_B1_dense + c] =
                    Hmat_ref[(size_t)r * cap_ref + c];
        double *eigs_dense = malloc(sizeof(double) * dim_B1_dense);
        irrep_hermitian_eigvals(dim_B1_dense, Hblock, eigs_dense);
        double E0_dense = eigs_dense[dim_B1_dense - 1]; /* sorted descending */

        /* Cross-check: sparse Lanczos GS == dense GS. */
        IRREP_ASSERT(fabs(eig_sp[0] - E0_dense) < 1e-10);
        /* And the well-known absolute kagome KH GS value. */
        IRREP_ASSERT(fabs(E0_dense - (-5.44487522)) < 1e-5);

        free(eigs_dense); free(Hblock); free(Hmat_ref);
        free(seed_sp);
        irrep_sg_heisenberg_sector_free(S_B1_sp);
        irrep_sg_little_group_irrep_free(B1_ref);

        /* ---- Dense 2D (Γ, E_1) cross-check against published N=12 values
         * from kagome12_k_resolved_ed:
         *   Γ-E_1 E_0 = -3.299516 J (dim 20 in character-projector basis).
         * Since the sparse path rejects 2D at the API boundary, the dense
         * reference IS the production 2D path for small-N. */
        irrep_sg_little_group_irrep_t *E1_ref =
            irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_E1);
        long long cap_2d = irrep_sg_rep_table_count(T);
        double _Complex *H_E1 = malloc((size_t)cap_2d * cap_2d * sizeof(double _Complex));
        int dim_E1_dense =
            irrep_sg_heisenberg_sector_build_dense(H, T, lg, E1_ref, H_E1, (int)cap_2d);
        IRREP_ASSERT(dim_E1_dense > 0);

        double _Complex *B_E1 = malloc((size_t)dim_E1_dense * dim_E1_dense *
                                       sizeof(double _Complex));
        for (int r = 0; r < dim_E1_dense; ++r)
            for (int c = 0; c < dim_E1_dense; ++c)
                B_E1[(size_t)r * dim_E1_dense + c] = H_E1[(size_t)r * cap_2d + c];

        /* Hermiticity of the dense 2D matrix: must be < 1e-12. */
        double herm_E1 = 0.0;
        for (int r = 0; r < dim_E1_dense; ++r)
            for (int c = 0; c < dim_E1_dense; ++c) {
                double e = cabs(B_E1[(size_t)r * dim_E1_dense + c] -
                                conj(B_E1[(size_t)c * dim_E1_dense + r]));
                if (e > herm_E1) herm_E1 = e;
            }
        IRREP_ASSERT(herm_E1 < 1e-12);

        double *eigs_E1 = malloc(sizeof(double) * dim_E1_dense);
        irrep_hermitian_eigvals(dim_E1_dense, B_E1, eigs_E1);
        double E0_E1 = eigs_E1[dim_E1_dense - 1];
        IRREP_ASSERT(fabs(E0_E1 - (-3.299516)) < 1e-5);

        free(eigs_E1); free(B_E1); free(H_E1);
        irrep_sg_little_group_irrep_free(E1_ref);

        irrep_sg_little_group_irrep_free(A1);
        irrep_sg_little_group_irrep_free(B1);
        irrep_sg_little_group_free(lg);
        irrep_sg_rep_table_free(T);
        free(bi); free(bj);
        irrep_heisenberg_free(H);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Sector-binding save/load round-trip --------------------------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int N  = irrep_space_group_num_sites(G);
        int nb = irrep_lattice_num_bonds_nn(L);
        int *bi = malloc(sizeof(int) * nb);
        int *bj = malloc(sizeof(int) * nb);
        irrep_lattice_fill_bonds_nn(L, bi, bj);
        irrep_heisenberg_t   *H = irrep_heisenberg_new(N, nb, bi, bj, 1.0);
        irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 6);

        /* Trivial-sector (real coefficients). */
        irrep_sg_heisenberg_sector_t *S  = irrep_sg_heisenberg_sector_build(H, T);
        long long n_reps = irrep_sg_heisenberg_sector_dim(S);

        const char *path_triv = "/tmp/libirrep_sector_triv.bin";
        IRREP_ASSERT(irrep_sg_heisenberg_sector_save(S, path_triv) == 0);

        irrep_sg_heisenberg_sector_t *S2 = irrep_sg_heisenberg_sector_load(path_triv, n_reps);
        IRREP_ASSERT(S2 != NULL);
        IRREP_ASSERT(irrep_sg_heisenberg_sector_dim(S2) == n_reps);

        /* Apply both, compare bit-exactly. */
        double _Complex *x  = malloc((size_t)n_reps * sizeof(double _Complex));
        double _Complex *y1 = malloc((size_t)n_reps * sizeof(double _Complex));
        double _Complex *y2 = malloc((size_t)n_reps * sizeof(double _Complex));
        for (long long i = 0; i < n_reps; ++i)
            x[i] = 0.13 * sin(0.27 * i) + I * 0.09 * cos(0.17 * i);
        irrep_sg_heisenberg_sector_apply(x, y1, S);
        irrep_sg_heisenberg_sector_apply(x, y2, S2);
        for (long long i = 0; i < n_reps; ++i)
            IRREP_ASSERT(cabs(y1[i] - y2[i]) < 1e-15);

        /* Bad expected dim rejects. */
        IRREP_ASSERT(irrep_sg_heisenberg_sector_load(path_triv, n_reps + 1) == NULL);
        /* Bad path. */
        IRREP_ASSERT(irrep_sg_heisenberg_sector_load("/nonexistent", n_reps) == NULL);
        /* NULL args. */
        IRREP_ASSERT(irrep_sg_heisenberg_sector_save(NULL, path_triv) == -1);
        IRREP_ASSERT(irrep_sg_heisenberg_sector_save(S, NULL) == -1);
        IRREP_ASSERT(irrep_sg_heisenberg_sector_load(NULL, -1) == NULL);

        irrep_sg_heisenberg_sector_free(S2);

        /* Corrupted-magic rejection. */
        {
            FILE *fp = fopen(path_triv, "r+b");
            IRREP_ASSERT(fp != NULL);
            char bad = 'X';
            fwrite(&bad, 1, 1, fp);
            fclose(fp);
            IRREP_ASSERT(irrep_sg_heisenberg_sector_load(path_triv, -1) == NULL);
        }

        /* Truncated (empty) file. */
        {
            const char *tp = "/tmp/libirrep_sector_trunc.bin";
            FILE *fp = fopen(tp, "wb");
            IRREP_ASSERT(fp != NULL);
            fclose(fp);
            IRREP_ASSERT(irrep_sg_heisenberg_sector_load(tp, -1) == NULL);
            remove(tp);
        }

        irrep_sg_heisenberg_sector_free(S);
        remove(path_triv);

        /* Complex-sector path ((k, μ_k), B_1 at Γ). */
        irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
        irrep_sg_little_group_irrep_t *B1 =
            irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_B1);
        irrep_sg_heisenberg_sector_t *S_B1 =
            irrep_sg_heisenberg_sector_build_at_k(H, T, lg, B1);
        long long dim_B1 = irrep_sg_heisenberg_sector_dim(S_B1);

        const char *path_cplx = "/tmp/libirrep_sector_cplx.bin";
        IRREP_ASSERT(irrep_sg_heisenberg_sector_save(S_B1, path_cplx) == 0);
        irrep_sg_heisenberg_sector_t *S_B1_2 =
            irrep_sg_heisenberg_sector_load(path_cplx, dim_B1);
        IRREP_ASSERT(S_B1_2 != NULL);

        /* Apply both, compare. */
        double _Complex *xb  = malloc((size_t)dim_B1 * sizeof(double _Complex));
        double _Complex *yb1 = malloc((size_t)dim_B1 * sizeof(double _Complex));
        double _Complex *yb2 = malloc((size_t)dim_B1 * sizeof(double _Complex));
        for (long long i = 0; i < dim_B1; ++i)
            xb[i] = 0.11 * sin(0.31 * i) + I * 0.07 * cos(0.19 * i);
        irrep_sg_heisenberg_sector_apply(xb, yb1, S_B1);
        irrep_sg_heisenberg_sector_apply(xb, yb2, S_B1_2);
        for (long long i = 0; i < dim_B1; ++i)
            IRREP_ASSERT(cabs(yb1[i] - yb2[i]) < 1e-15);

        free(xb); free(yb1); free(yb2);
        irrep_sg_heisenberg_sector_free(S_B1_2);
        irrep_sg_heisenberg_sector_free(S_B1);
        irrep_sg_little_group_irrep_free(B1);
        irrep_sg_little_group_free(lg);
        remove(path_cplx);

        free(x); free(y1); free(y2);
        free(bi); free(bj);
        irrep_sg_rep_table_free(T);
        irrep_heisenberg_free(H);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Hermiticity: H_{jk} = conj(H_{kj}) across all reps ------------ */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int N  = irrep_space_group_num_sites(G);
        int nb = irrep_lattice_num_bonds_nn(L);
        int *bi = malloc(sizeof(int) * nb);
        int *bj = malloc(sizeof(int) * nb);
        irrep_lattice_fill_bonds_nn(L, bi, bj);

        irrep_heisenberg_t   *H = irrep_heisenberg_new(N, nb, bi, bj, 1.0);
        irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 6);
        long long n = irrep_sg_rep_table_count(T);

        double _Complex *Hmat   = malloc((size_t)n * n * sizeof(double _Complex));
        double _Complex *sp_in  = malloc((size_t)n * sizeof(double _Complex));
        double _Complex *sp_out = malloc((size_t)n * sizeof(double _Complex));

        for (long long k = 0; k < n; ++k) {
            for (long long j = 0; j < n; ++j)
                sp_in[j] = (j == k) ? 1.0 : 0.0;
            irrep_heisenberg_apply_in_sector(H, T, sp_in, sp_out);
            for (long long j = 0; j < n; ++j)
                Hmat[(size_t)j * n + k] = sp_out[j]; /* column k */
        }
        for (long long j = 0; j < n; ++j) {
            for (long long k = 0; k < n; ++k) {
                double _Complex a = Hmat[(size_t)j * n + k];
                double _Complex b = Hmat[(size_t)k * n + j];
                IRREP_ASSERT(cabs(a - conj(b)) < 1e-13);
            }
        }

        /* Diagonalise and compare GS to the known Γ-A_1 Sz=0 lowest
         * eigenvalue. The A_1 sector's full-Sz GS on kagome 2×2 is
         * E_0 = −5.328392 J at dim 144; restricted to Sz=0 it must
         * coincide (verified by per-popcount enumeration). */
        double *eigvals = malloc(sizeof(double) * n);
        irrep_hermitian_eigvals(n, Hmat, eigvals);
        double E0 = eigvals[n - 1]; /* descending */
        IRREP_ASSERT(fabs(E0 - (-5.328392)) < 1e-5);

        free(eigvals); free(Hmat); free(sp_in); free(sp_out);
        free(bi); free(bj);
        irrep_sg_rep_table_free(T);
        irrep_heisenberg_free(H);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    return IRREP_TEST_END();
}
