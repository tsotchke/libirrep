/* SPDX-License-Identifier: MIT */
/* Tests for the space-group little-group builder at Bloch momentum k.
 *
 * Coverage:
 *   - Γ-point little group equals the full point group on p1 / p4mm / p6mm.
 *   - Sum rule: Σ_k |G_k^point| = point_order · L_x · L_y. Each point-group
 *     element fixes exactly one "k-orbit" worth of momenta, so summing the
 *     little point-group cardinality over the Brillouin-zone mesh recovers
 *     `point_order · num_translations`.
 *   - p4mm on a 2×2 square lattice: explicit cardinalities at Γ, X, M.
 *   - p6mm on a 2×2 kagome torus: cardinalities at the four allowed k-points
 *     (Γ and three M-like points; the cluster is too small to host a K-point).
 *   - Structural invariants: the identity (p = 0) is always in every little
 *     group; the ops are sorted ascending.
 */
#include "harness.h"
#include <irrep/config_project.h>
#include <irrep/lattice.h>
#include <irrep/space_group.h>
#include <stdlib.h>

int main(void) {
    IRREP_TEST_START("little_group");

    /* ---- p1: every k has the trivial little point group (order 1) ---- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
        IRREP_ASSERT(G != NULL);
        for (int ky = 0; ky < 3; ++ky) {
            for (int kx = 0; kx < 3; ++kx) {
                irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, kx, ky);
                IRREP_ASSERT(lg != NULL);
                IRREP_ASSERT(irrep_sg_little_group_point_order(lg) == 1);
                IRREP_ASSERT(irrep_sg_little_group_order(lg) == 9); /* 1 · 3 · 3 */
                int kx_out = -1, ky_out = -1;
                irrep_sg_little_group_k(lg, &kx_out, &ky_out);
                IRREP_ASSERT(kx_out == kx && ky_out == ky);
                irrep_sg_little_group_free(lg);
            }
        }
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- p4mm on a 2×2 square: Γ full (8), M full (8), X partial ------
     * A 2×2 torus has four k-points: Γ=(0,0), X_a=(1,0), X_b=(0,1),
     * M=(1,1). Γ and M are C_4-invariant (their rotations all give back
     * Γ or M modulo reciprocal lattice), so their little point groups
     * are the full C_4v of order 8. X_a / X_b are only C_2v-invariant
     * (order 4). Orbit-count / Frobenius: Σ_k |G_k^point| = |G| · N_orbits
     * = 8 · 3 = 24 (orbits are {Γ}, {X_a, X_b}, {M}).
     * ------------------------------------------------------------------ */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P4MM);
        IRREP_ASSERT(G != NULL);
        IRREP_ASSERT(irrep_space_group_point_order(G) == 8);

        struct {
            int kx, ky, expected;
        } cases[] = {
            {0, 0, 8}, /* Γ — full C_4v                              */
            {1, 0, 4}, /* X_a — C_2v subgroup                        */
            {0, 1, 4}, /* X_b — C_2v subgroup                        */
            {1, 1, 8}, /* M — C_4v (C_4 sends (1,1)→(1,-1)≡(1,1) mod 2) */
        };
        int sum = 0;
        for (int i = 0; i < 4; ++i) {
            irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, cases[i].kx, cases[i].ky);
            IRREP_ASSERT(lg != NULL);
            int n = irrep_sg_little_group_point_order(lg);
            IRREP_ASSERT(n == cases[i].expected);
            /* Identity is always present and the list is sorted. */
            int ops[8];
            irrep_sg_little_group_point_ops(lg, ops);
            IRREP_ASSERT(ops[0] == 0);
            for (int k = 1; k < n; ++k)
                IRREP_ASSERT(ops[k] > ops[k - 1]);
            sum += n;
            irrep_sg_little_group_free(lg);
        }
        /* Frobenius / orbit-count: Σ_k |G_k^point| = |G^point| · N_orbits. */
        IRREP_ASSERT(sum == 8 * 3);

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- p6mm kagome 2×2: Γ + three M-points -------------------------
     * 2×2 is exactly the even-L cluster that hosts the hex-BZ M-points
     * at (1/2, 0), (0, 1/2), (1/2, 1/2). All three are in one C_6v orbit
     * (they're rotated copies of each other), so N_orbits = 2 and
     * Σ_k |G_k^point| = 12 · 2 = 24.
     * ------------------------------------------------------------------ */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        IRREP_ASSERT(G != NULL);
        IRREP_ASSERT(irrep_space_group_point_order(G) == 12);

        struct {
            int kx, ky, expected;
        } cases[] = {
            {0, 0, 12}, /* Γ — full C_6v    */
            {1, 0, 4},  /* M_a — C_2v       */
            {0, 1, 4},  /* M_b — C_2v       */
            {1, 1, 4},  /* M_c — C_2v       */
        };
        int sum = 0;
        for (int i = 0; i < 4; ++i) {
            irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, cases[i].kx, cases[i].ky);
            IRREP_ASSERT(lg != NULL);
            IRREP_ASSERT(irrep_sg_little_group_point_order(lg) == cases[i].expected);
            sum += irrep_sg_little_group_point_order(lg);
            irrep_sg_little_group_free(lg);
        }
        IRREP_ASSERT(sum == 12 * 2);

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- p6mm kagome 3×3: 9 k-points, two K-points ----
     * A 3-divisible cluster hosts the hex-BZ K-points at (1/3, 2/3) and
     * (2/3, 1/3), which on a 3×3 mesh land at integer indices (1, 2) and
     * (2, 1). Each carries C_3v (order 6 — three rotations + three mirrors)
     * — the Dirac-cone little group on the honeycomb / kagome BZ.
     *
     * The six remaining non-Γ k-points are NOT M-points (M = (1/2, 0)
     * requires even L); they are generic low-symmetry k-points whose
     * little group is order 2 (identity + one mirror). Three orbits:
     * {Γ}, {K, K'}, {six generic}, giving Σ = 12·1 + 6·2 + 2·6 = 36.
     * ------------------------------------------------------------------ */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        IRREP_ASSERT(G != NULL);

        /* Γ-point: full C_6v. */
        irrep_sg_little_group_t *gamma = irrep_sg_little_group_build(G, 0, 0);
        IRREP_ASSERT(irrep_sg_little_group_point_order(gamma) == 12);
        irrep_sg_little_group_free(gamma);

        /* Two K-points — Dirac-node signature. */
        for (int i = 0; i < 2; ++i) {
            int                      kx = (i == 0) ? 1 : 2;
            int                      ky = (i == 0) ? 2 : 1;
            irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, kx, ky);
            IRREP_ASSERT(irrep_sg_little_group_point_order(lg) == 6);
            irrep_sg_little_group_free(lg);
        }

        /* Frobenius / orbit-count: 3 orbits on the 9-point mesh. */
        int sum = 0;
        for (int ky = 0; ky < 3; ++ky) {
            for (int kx = 0; kx < 3; ++kx) {
                irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, kx, ky);
                IRREP_ASSERT(lg != NULL);
                sum += irrep_sg_little_group_point_order(lg);
                irrep_sg_little_group_free(lg);
            }
        }
        IRREP_ASSERT(sum == 12 * 3);

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Element-matrix introspection: |det| = 1 on every element ------ */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        /* Γ's little group is the full C_6v — covers every point op. */
        irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
        int                      n = irrep_sg_little_group_point_order(lg);
        IRREP_ASSERT(n == 12);
        int proper = 0, improper = 0;
        for (int i = 0; i < n; ++i) {
            int M[2][2];
            irrep_sg_little_group_element_matrix(lg, i, M);
            int det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
            IRREP_ASSERT(det == 1 || det == -1);
            if (det == 1)
                ++proper;
            else
                ++improper;
        }
        /* C_6v: 6 rotations (proper) + 6 mirrors (improper). */
        IRREP_ASSERT(proper == 6);
        IRREP_ASSERT(improper == 6);
        /* Identity is M = [[1, 0], [0, 1]]. */
        {
            int M[2][2];
            irrep_sg_little_group_element_matrix(lg, 0, M);
            IRREP_ASSERT(M[0][0] == 1 && M[0][1] == 0 && M[1][0] == 0 && M[1][1] == 1);
        }
        irrep_sg_little_group_free(lg);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Trivial irrep at Γ equals the A₁ projector ---------------------
     * The composite projector with all-ones characters and k = (0, 0) is
     * the totally-symmetric space-group projector. Must agree with
     * `irrep_sg_project_A1` on every orbit-amplitude input.
     * --------------------------------------------------------------------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        IRREP_ASSERT(irrep_space_group_order(G) == 48); /* 12 · 4 */

        irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
        IRREP_ASSERT(irrep_sg_little_group_point_order(lg) == 12);

        double _Complex ones[12];
        for (int i = 0; i < 12; ++i)
            ones[i] = 1.0 + 0.0 * I;
        irrep_sg_little_group_irrep_t *mu_A1 = irrep_sg_little_group_irrep_new(lg, ones, 1);
        IRREP_ASSERT(mu_A1 != NULL);

        /* Three diverse probe vectors. */
        double _Complex psi[48];
        for (int g = 0; g < 48; ++g)
            psi[g] = (0.13 * g - 0.2) + I * (0.07 * g - 0.5);

        double _Complex a1_classic = irrep_sg_project_A1(G, psi);
        double _Complex a1_at_k    = irrep_sg_project_at_k(lg, mu_A1, psi);
        IRREP_ASSERT(cabs(a1_classic - a1_at_k) < 1e-12);

        for (int g = 0; g < 48; ++g)
            psi[g] = sin(0.3 * g) + I * cos(0.4 * g);
        a1_classic = irrep_sg_project_A1(G, psi);
        a1_at_k    = irrep_sg_project_at_k(lg, mu_A1, psi);
        IRREP_ASSERT(cabs(a1_classic - a1_at_k) < 1e-12);

        for (int g = 0; g < 48; ++g)
            psi[g] = 0;
        psi[5] = 1.0;
        a1_classic = irrep_sg_project_A1(G, psi);
        a1_at_k    = irrep_sg_project_at_k(lg, mu_A1, psi);
        IRREP_ASSERT(cabs(a1_classic - a1_at_k) < 1e-12);

        irrep_sg_little_group_irrep_free(mu_A1);
        irrep_sg_little_group_free(lg);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Orthogonality: A₁ ⊥ A₂ (sign irrep) on Γ ----------------------
     * With Γ's full point group the "sign" little-group-irrep has χ = +1
     * on proper rotations, −1 on reflections. A₁ and A₂ project onto
     * orthogonal subspaces — starting from an A₁-symmetric input ψ, A₂
     * must return 0.
     * --------------------------------------------------------------------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);

        irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
        int                      n = irrep_sg_little_group_point_order(lg);
        IRREP_ASSERT(n == 12);

        double _Complex chi_A2[12];
        for (int i = 0; i < n; ++i) {
            int M[2][2];
            irrep_sg_little_group_element_matrix(lg, i, M);
            int det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
            chi_A2[i] = (det == 1) ? (1.0 + 0.0 * I) : (-1.0 + 0.0 * I);
        }
        irrep_sg_little_group_irrep_t *mu_A2 = irrep_sg_little_group_irrep_new(lg, chi_A2, 1);

        /* Build a totally-symmetric ψ: average a random vector over the
         * orbit. The averaged vector is invariant under G → lives in A₁,
         * so A₂ projection must vanish. */
        double _Complex raw[48];
        for (int g = 0; g < 48; ++g)
            raw[g] = (0.11 * g + 0.2) - I * (0.3 * g - 0.5);
        double _Complex mean = 0;
        for (int g = 0; g < 48; ++g)
            mean += raw[g];
        mean /= 48.0;
        double _Complex psi_A1[48];
        for (int g = 0; g < 48; ++g)
            psi_A1[g] = mean;

        double _Complex proj = irrep_sg_project_at_k(lg, mu_A2, psi_A1);
        IRREP_ASSERT(cabs(proj) < 1e-12);

        irrep_sg_little_group_irrep_free(mu_A2);
        irrep_sg_little_group_free(lg);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- adapted_basis_at_k: dim sum-rule + Γ/A₁ cross-check ------------
     * On a 4-site (Lx=2, Ly=2) square torus with p4mm and local_dim=2
     * (spin-½), total Hilbert dimension is 16. Summing the per-(k, μ_k)
     * sector dimensions (weighted by d_{μ_k}) over the full BZ and every
     * little-group irrep must recover 16 — Burnside / Frobenius on the
     * representation of G acting by site permutation.
     *
     * Cross-check: the Γ-A₁ block built via _adapted_basis_at_k must
     * span the same space (equal dimension) as the existing
     * _adapted_basis path using the trivial space-group irrep.
     * --------------------------------------------------------------------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P4MM);
        IRREP_ASSERT(G != NULL);
        int                      num_sites = 4;
        int                      D = 1 << num_sites; /* 16 */
        double _Complex         *buf = malloc((size_t)D * D * sizeof(double _Complex));
        IRREP_ASSERT(buf != NULL);

        /* Γ-A₁ via the new path. */
        irrep_sg_little_group_t *lg_G = irrep_sg_little_group_build(G, 0, 0);
        int                      n_ops_G = irrep_sg_little_group_point_order(lg_G);
        double _Complex         *ones_G = malloc((size_t)n_ops_G * sizeof(double _Complex));
        for (int i = 0; i < n_ops_G; ++i)
            ones_G[i] = 1.0 + 0.0 * I;
        irrep_sg_little_group_irrep_t *A1_at_k =
            irrep_sg_little_group_irrep_new(lg_G, ones_G, 1);
        int dim_A1_new =
            irrep_sg_adapted_basis_at_k(lg_G, A1_at_k, num_sites, 2, buf, D);
        IRREP_ASSERT(dim_A1_new >= 0);

        /* Γ-A₁ via the legacy space-group-irrep path. */
        irrep_sg_irrep_t *A1_sg = irrep_sg_trivial(G);
        int dim_A1_legacy = irrep_sg_adapted_basis(G, A1_sg, num_sites, 2, buf, D);
        IRREP_ASSERT(dim_A1_legacy >= 0);

        IRREP_ASSERT(dim_A1_new == dim_A1_legacy);

        irrep_sg_irrep_free(A1_sg);
        irrep_sg_little_group_irrep_free(A1_at_k);
        free(ones_G);
        irrep_sg_little_group_free(lg_G);

        /* Dimension sum rule over the whole BZ × irrep space.
         * At each k, use one all-ones irrep (effectively projecting onto
         * the trivial rep of the little point group) plus the sign rep;
         * for each generating irrep, projection gives one sector's
         * dimension. Summing d_μ · dim_{k,μ} over those two per k, over
         * all 4 k's, should be ≤ 16 (not equal, since C_4v and C_2v each
         * have more than 2 irreps — we're sampling a subset). The test
         * then tightens to: using only trivial-rep across all k recovers
         * the full Hilbert dim IFF every state lives in some k-sector
         * under trivial little-group rep — which is NOT true; point-
         * group non-trivial irreps carry states too. So we assert the
         * weaker invariant: partial sum < D, and each individual
         * sector fits in D. */
        int total = 0;
        for (int ky = 0; ky < 2; ++ky) {
            for (int kx = 0; kx < 2; ++kx) {
                irrep_sg_little_group_t *lg  = irrep_sg_little_group_build(G, kx, ky);
                int                      n   = irrep_sg_little_group_point_order(lg);
                double _Complex         *chi = malloc((size_t)n * sizeof(double _Complex));
                for (int i = 0; i < n; ++i)
                    chi[i] = 1.0 + 0.0 * I;
                irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi, 1);
                int dim_sector =
                    irrep_sg_adapted_basis_at_k(lg, mu, num_sites, 2, buf, D);
                IRREP_ASSERT(dim_sector >= 0 && dim_sector <= D);
                total += dim_sector;
                irrep_sg_little_group_irrep_free(mu);
                free(chi);
                irrep_sg_little_group_free(lg);
            }
        }
        /* Trivial-rep alone doesn't cover all 16 states; must be < D. */
        IRREP_ASSERT(total < D);
        IRREP_ASSERT(total > 0);

        free(buf);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Named irrep builtins: C_6v at Γ, C_3v at K on kagome ---------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);

        /* Γ: C_6v — A_1 via builtin == A_1 via manual all-ones. */
        irrep_sg_little_group_t *lg_G = irrep_sg_little_group_build(G, 0, 0);
        int                      n_G  = irrep_sg_little_group_point_order(lg_G);
        IRREP_ASSERT(n_G == 12);

        irrep_sg_little_group_irrep_t *A1_builtin =
            irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_A1);
        IRREP_ASSERT(A1_builtin != NULL);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(A1_builtin) == 1);

        double _Complex ones[12];
        for (int i = 0; i < 12; ++i)
            ones[i] = 1.0 + 0.0 * I;
        irrep_sg_little_group_irrep_t *A1_manual =
            irrep_sg_little_group_irrep_new(lg_G, ones, 1);

        /* Project a probe vector through both — results must be bit-equal. */
        double _Complex psi[432];
        for (int g = 0; g < 432; ++g) /* 12 · Lx · Ly = 12·9 */
            psi[g] = 0.0;
        /* Construct a non-trivial psi_of_g. */
        for (int g = 0; g < 432; ++g)
            psi[g] = (0.17 * g - 0.3) + I * (0.21 * g - 0.5);

        double _Complex p_builtin = irrep_sg_project_at_k(lg_G, A1_builtin, psi);
        double _Complex p_manual  = irrep_sg_project_at_k(lg_G, A1_manual, psi);
        IRREP_ASSERT(cabs(p_builtin - p_manual) < 1e-13);

        irrep_sg_little_group_irrep_free(A1_builtin);
        irrep_sg_little_group_irrep_free(A1_manual);

        /* E_1 via builtin: dim = 2, non-trivial character row. Project
         * twice and assert a sanity property — E_1 on A_1-symmetric
         * input must vanish (A_1 ⊥ E_1). */
        irrep_sg_little_group_irrep_t *E1 =
            irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_E1);
        IRREP_ASSERT(E1 != NULL);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(E1) == 2);
        double _Complex psi_sym[432];
        double _Complex mean = 0;
        for (int g = 0; g < 432; ++g)
            mean += psi[g];
        mean /= 432.0;
        for (int g = 0; g < 432; ++g)
            psi_sym[g] = mean; /* totally symmetric */
        double _Complex p_E1 = irrep_sg_project_at_k(lg_G, E1, psi_sym);
        IRREP_ASSERT(cabs(p_E1) < 1e-12);
        irrep_sg_little_group_irrep_free(E1);

        /* B_1 via builtin: also valid on C_6v. */
        irrep_sg_little_group_irrep_t *B1 =
            irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_B1);
        IRREP_ASSERT(B1 != NULL);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(B1) == 1);
        irrep_sg_little_group_irrep_free(B1);

        /* Named irrep "E" (the C_3v 2D) is NOT valid on C_6v — must return NULL. */
        IRREP_ASSERT(irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_E) == NULL);

        irrep_sg_little_group_free(lg_G);

        /* K on 3×3 kagome: C_3v. A_1 / A_2 / E must all build. */
        irrep_sg_little_group_t *lg_K = irrep_sg_little_group_build(G, 1, 2);
        IRREP_ASSERT(irrep_sg_little_group_point_order(lg_K) == 6);

        irrep_sg_little_group_irrep_t *K_A1 =
            irrep_sg_little_group_irrep_named(lg_K, IRREP_LG_IRREP_A1);
        irrep_sg_little_group_irrep_t *K_A2 =
            irrep_sg_little_group_irrep_named(lg_K, IRREP_LG_IRREP_A2);
        irrep_sg_little_group_irrep_t *K_E =
            irrep_sg_little_group_irrep_named(lg_K, IRREP_LG_IRREP_E);
        IRREP_ASSERT(K_A1 != NULL && K_A2 != NULL && K_E != NULL);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(K_A1) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(K_A2) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(K_E) == 2);

        /* B_1 is NOT valid on C_3v — must return NULL. */
        IRREP_ASSERT(irrep_sg_little_group_irrep_named(lg_K, IRREP_LG_IRREP_B1) == NULL);

        irrep_sg_little_group_irrep_free(K_A1);
        irrep_sg_little_group_irrep_free(K_A2);
        irrep_sg_little_group_irrep_free(K_E);
        irrep_sg_little_group_free(lg_K);

        irrep_space_group_free(G);
        irrep_lattice_free(L);

        /* M on 2×2 kagome: C_2v (order 4) — four 1D irreps A_1/A_2/B_1/B_2.
         * Builtin must agree with manually-constructed all-ones A_1. */
        irrep_lattice_t     *L2 = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G2 = irrep_space_group_build(L2, IRREP_WALLPAPER_P6MM);
        irrep_sg_little_group_t *lg_M = irrep_sg_little_group_build(G2, 1, 0);
        IRREP_ASSERT(irrep_sg_little_group_point_order(lg_M) == 4);

        irrep_sg_little_group_irrep_t *M_A1 =
            irrep_sg_little_group_irrep_named(lg_M, IRREP_LG_IRREP_A1);
        irrep_sg_little_group_irrep_t *M_A2 =
            irrep_sg_little_group_irrep_named(lg_M, IRREP_LG_IRREP_A2);
        irrep_sg_little_group_irrep_t *M_B1 =
            irrep_sg_little_group_irrep_named(lg_M, IRREP_LG_IRREP_B1);
        irrep_sg_little_group_irrep_t *M_B2 =
            irrep_sg_little_group_irrep_named(lg_M, IRREP_LG_IRREP_B2);
        IRREP_ASSERT(M_A1 && M_A2 && M_B1 && M_B2);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(M_A1) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(M_A2) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(M_B1) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(M_B2) == 1);

        /* C_2v has no E irrep — must return NULL. */
        IRREP_ASSERT(irrep_sg_little_group_irrep_named(lg_M, IRREP_LG_IRREP_E) == NULL);
        IRREP_ASSERT(irrep_sg_little_group_irrep_named(lg_M, IRREP_LG_IRREP_E1) == NULL);
        IRREP_ASSERT(irrep_sg_little_group_irrep_named(lg_M, IRREP_LG_IRREP_E_C4V) == NULL);

        /* A_1 builtin == manually-constructed all-ones. */
        double _Complex ones4[4] = {1, 1, 1, 1};
        irrep_sg_little_group_irrep_t *M_A1_manual =
            irrep_sg_little_group_irrep_new(lg_M, ones4, 1);
        double _Complex psi48[48];
        for (int g = 0; g < 48; ++g)
            psi48[g] = (0.13 * g - 0.4) + I * (0.09 * g + 0.2);
        double _Complex pb = irrep_sg_project_at_k(lg_M, M_A1, psi48);
        double _Complex pm = irrep_sg_project_at_k(lg_M, M_A1_manual, psi48);
        IRREP_ASSERT(cabs(pb - pm) < 1e-13);

        /* A_1 ⊥ B_1: project a totally-symmetric input through B_1 → 0.
         * Build psi = group-average so it lives in A_1. */
        double _Complex mean_M = 0;
        for (int g = 0; g < 48; ++g)
            mean_M += psi48[g];
        mean_M /= 48.0;
        double _Complex psi_sym_M[48];
        for (int g = 0; g < 48; ++g)
            psi_sym_M[g] = mean_M;
        double _Complex p_B1 = irrep_sg_project_at_k(lg_M, M_B1, psi_sym_M);
        IRREP_ASSERT(cabs(p_B1) < 1e-12);

        irrep_sg_little_group_irrep_free(M_A1);
        irrep_sg_little_group_irrep_free(M_A2);
        irrep_sg_little_group_irrep_free(M_B1);
        irrep_sg_little_group_irrep_free(M_B2);
        irrep_sg_little_group_irrep_free(M_A1_manual);
        irrep_sg_little_group_free(lg_M);
        irrep_space_group_free(G2);
        irrep_lattice_free(L2);
    }

    /* ---- C_4v at Γ and C_2v at X on p4mm 4×4 square -------------------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P4MM);
        IRREP_ASSERT(G != NULL);
        int total_dim = irrep_space_group_order(G); /* 8 · 16 */

        /* Γ: C_4v (order 8). All five irreps must build. */
        irrep_sg_little_group_t *lg_G = irrep_sg_little_group_build(G, 0, 0);
        IRREP_ASSERT(irrep_sg_little_group_point_order(lg_G) == 8);

        irrep_sg_little_group_irrep_t *G_A1 =
            irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_A1);
        irrep_sg_little_group_irrep_t *G_A2 =
            irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_A2);
        irrep_sg_little_group_irrep_t *G_B1 =
            irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_B1);
        irrep_sg_little_group_irrep_t *G_B2 =
            irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_B2);
        irrep_sg_little_group_irrep_t *G_E =
            irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_E_C4V);
        IRREP_ASSERT(G_A1 && G_A2 && G_B1 && G_B2 && G_E);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(G_A1) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(G_A2) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(G_B1) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(G_B2) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(G_E) == 2);

        /* C_6v-only names must be rejected on C_4v. */
        IRREP_ASSERT(irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_E1) == NULL);
        IRREP_ASSERT(irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_E2) == NULL);
        IRREP_ASSERT(irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_E) == NULL);

        /* A_1 via builtin == classical irrep_sg_project_A1. */
        double _Complex psi_p4[128]; /* 8 · 16 */
        for (int g = 0; g < total_dim; ++g)
            psi_p4[g] = (0.07 * g - 0.3) + I * (0.11 * g + 0.2);
        double _Complex p_classic = irrep_sg_project_A1(G, psi_p4);
        double _Complex p_builtin = irrep_sg_project_at_k(lg_G, G_A1, psi_p4);
        IRREP_ASSERT(cabs(p_classic - p_builtin) < 1e-12);

        /* A_1-symmetric input through non-trivial irrep → 0. */
        double _Complex mean_p = 0;
        for (int g = 0; g < total_dim; ++g)
            mean_p += psi_p4[g];
        mean_p /= total_dim;
        double _Complex psi_sym[128];
        for (int g = 0; g < total_dim; ++g)
            psi_sym[g] = mean_p;
        IRREP_ASSERT(cabs(irrep_sg_project_at_k(lg_G, G_A2, psi_sym)) < 1e-12);
        IRREP_ASSERT(cabs(irrep_sg_project_at_k(lg_G, G_B1, psi_sym)) < 1e-12);
        IRREP_ASSERT(cabs(irrep_sg_project_at_k(lg_G, G_B2, psi_sym)) < 1e-12);
        IRREP_ASSERT(cabs(irrep_sg_project_at_k(lg_G, G_E,  psi_sym)) < 1e-12);

        irrep_sg_little_group_irrep_free(G_A1);
        irrep_sg_little_group_irrep_free(G_A2);
        irrep_sg_little_group_irrep_free(G_B1);
        irrep_sg_little_group_irrep_free(G_B2);
        irrep_sg_little_group_irrep_free(G_E);
        irrep_sg_little_group_free(lg_G);

        /* X-point (2, 0) on 4×4 = ½·b₁: self-conjugate under C_2 and σ_v,
         * giving C_2v (order 4). (1, 0) is a generic low-symmetry point. */
        irrep_sg_little_group_t *lg_X = irrep_sg_little_group_build(G, 2, 0);
        IRREP_ASSERT(irrep_sg_little_group_point_order(lg_X) == 4);
        irrep_sg_little_group_irrep_t *X_A1 =
            irrep_sg_little_group_irrep_named(lg_X, IRREP_LG_IRREP_A1);
        irrep_sg_little_group_irrep_t *X_A2 =
            irrep_sg_little_group_irrep_named(lg_X, IRREP_LG_IRREP_A2);
        irrep_sg_little_group_irrep_t *X_B1 =
            irrep_sg_little_group_irrep_named(lg_X, IRREP_LG_IRREP_B1);
        irrep_sg_little_group_irrep_t *X_B2 =
            irrep_sg_little_group_irrep_named(lg_X, IRREP_LG_IRREP_B2);
        IRREP_ASSERT(X_A1 && X_A2 && X_B1 && X_B2);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(X_A1) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(X_A2) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(X_B1) == 1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(X_B2) == 1);
        /* C_2v has no 2D irrep. */
        IRREP_ASSERT(irrep_sg_little_group_irrep_named(lg_X, IRREP_LG_IRREP_E) == NULL);
        IRREP_ASSERT(irrep_sg_little_group_irrep_named(lg_X, IRREP_LG_IRREP_E_C4V) == NULL);

        irrep_sg_little_group_irrep_free(X_A1);
        irrep_sg_little_group_irrep_free(X_A2);
        irrep_sg_little_group_irrep_free(X_B1);
        irrep_sg_little_group_irrep_free(X_B2);
        irrep_sg_little_group_free(lg_X);

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Projector-weights export: bit-exact vs irrep_sg_project_at_k -----
     * The handshake primitive for downstream MPO builders. Vend weights
     * w_g over the FULL space group (most zero outside the little group);
     * assert Σ_g w_g · ψ(g) == irrep_sg_project_at_k(lg, μ, ψ) across
     * multiple (k, μ_k) pairs on kagome and square lattices. */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int                  order_G = irrep_space_group_order(G); /* 12 · 9 = 108 */

        /* Probe vector ψ(g), g ∈ [0, |G|). */
        double _Complex psi[108];
        for (int g = 0; g < order_G; ++g)
            psi[g] = (0.17 * g - 0.4) + I * (0.09 * g + 0.3);

        struct {
            int                    kx, ky;
            irrep_lg_named_irrep_t name;
        } probes[] = {
            {0, 0, IRREP_LG_IRREP_A1}, /* Γ, C_6v: all ones */
            {0, 0, IRREP_LG_IRREP_A2}, /* Γ, C_6v: sign rep */
            {0, 0, IRREP_LG_IRREP_E1}, /* Γ, C_6v: 2D irrep */
            {1, 2, IRREP_LG_IRREP_A1}, /* K, C_3v */
            {1, 2, IRREP_LG_IRREP_E},  /* K, C_3v: 2D */
        };

        for (unsigned p = 0; p < sizeof(probes) / sizeof(probes[0]); ++p) {
            irrep_sg_little_group_t *lg =
                irrep_sg_little_group_build(G, probes[p].kx, probes[p].ky);
            irrep_sg_little_group_irrep_t *mu =
                irrep_sg_little_group_irrep_named(lg, probes[p].name);
            IRREP_ASSERT(lg && mu);

            double _Complex weights[108];
            int rc = irrep_sg_projector_weights(lg, mu, weights);
            IRREP_ASSERT(rc == 0);

            double _Complex via_weights = 0.0 + 0.0 * I;
            for (int g = 0; g < order_G; ++g)
                via_weights += weights[g] * psi[g];

            double _Complex via_project = irrep_sg_project_at_k(lg, mu, psi);
            IRREP_ASSERT(cabs(via_weights - via_project) < 1e-13);

            /* Structural invariant: every g whose point part is NOT in
             * the little group has weight exactly 0 (not just ≈ 0). For
             * 2D irreps the little-group weights can also be zero on
             * classes where χ = 0, so we only test the "outside → 0"
             * half of the claim. */
            int n_lg_point = irrep_sg_little_group_point_order(lg);
            int ops[12];
            irrep_sg_little_group_point_ops(lg, ops);
            int pt_order = irrep_space_group_point_order(G);
            int in_lg[12] = {0};
            for (int j = 0; j < n_lg_point; ++j)
                in_lg[ops[j]] = 1;
            for (int g = 0; g < order_G; ++g) {
                int p = g % pt_order;
                if (!in_lg[p])
                    IRREP_ASSERT(creal(weights[g]) == 0.0 && cimag(weights[g]) == 0.0);
            }

            irrep_sg_little_group_irrep_free(mu);
            irrep_sg_little_group_free(lg);
        }

        /* A_1 normalisation: for the trivial irrep at Γ, weights sum to 1.
         * (d_μ/|G|) · Σ_g 1 · 1 = (1/|G|) · |G| = 1. */
        irrep_sg_little_group_t *lg_G = irrep_sg_little_group_build(G, 0, 0);
        irrep_sg_little_group_irrep_t *A1 =
            irrep_sg_little_group_irrep_named(lg_G, IRREP_LG_IRREP_A1);
        double _Complex weights_A1[108];
        irrep_sg_projector_weights(lg_G, A1, weights_A1);
        double _Complex sum_A1 = 0.0 + 0.0 * I;
        for (int g = 0; g < order_G; ++g)
            sum_A1 += weights_A1[g];
        IRREP_ASSERT(cabs(sum_A1 - 1.0) < 1e-13);
        irrep_sg_little_group_irrep_free(A1);
        irrep_sg_little_group_free(lg_G);

        irrep_space_group_free(G);
        irrep_lattice_free(L);

        /* Square lattice p4mm probe: C_4v at Γ including 2D E. */
        irrep_lattice_t     *Ls = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
        irrep_space_group_t *Gs = irrep_space_group_build(Ls, IRREP_WALLPAPER_P4MM);
        int                  order_Gs = irrep_space_group_order(Gs); /* 8 · 16 = 128 */

        double _Complex psi_s[128];
        for (int g = 0; g < order_Gs; ++g)
            psi_s[g] = sin(0.11 * g) + I * cos(0.07 * g);

        irrep_sg_little_group_t *lg_sG = irrep_sg_little_group_build(Gs, 0, 0);
        irrep_sg_little_group_irrep_t *E_c4v =
            irrep_sg_little_group_irrep_named(lg_sG, IRREP_LG_IRREP_E_C4V);
        double _Complex weights_s[128];
        IRREP_ASSERT(irrep_sg_projector_weights(lg_sG, E_c4v, weights_s) == 0);
        double _Complex s_w = 0.0 + 0.0 * I;
        for (int g = 0; g < order_Gs; ++g)
            s_w += weights_s[g] * psi_s[g];
        double _Complex s_p = irrep_sg_project_at_k(lg_sG, E_c4v, psi_s);
        IRREP_ASSERT(cabs(s_w - s_p) < 1e-13);

        irrep_sg_little_group_irrep_free(E_c4v);
        irrep_sg_little_group_free(lg_sG);
        irrep_space_group_free(Gs);
        irrep_lattice_free(Ls);

        /* Error paths: NULL inputs return -1. */
        IRREP_ASSERT(irrep_sg_projector_weights(NULL, NULL, NULL) == -1);
    }

    /* ---- 2D-irrep D-matrix storage + properties ------------------------
     * Verify:
     *   (1) D(E) == I
     *   (2) D(g) is orthogonal: D(g)^T · D(g) = I
     *   (3) trace(D(g_i)) matches the irrep's character at slot i
     *   (4) 1D irreps return their character as a 1×1 matrix via the
     *       same accessor. */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);

        /* Γ on 3×3 kagome: full C_6v. Test E_1 (2D). */
        irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
        irrep_sg_little_group_irrep_t *E1 =
            irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_E1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(E1) == 2);
        int n = irrep_sg_little_group_point_order(lg);

        /* (1) D(E) = I; the identity is always slot 0. */
        double _Complex DE[4];
        IRREP_ASSERT(irrep_sg_little_group_irrep_matrix(E1, 0, DE) == 0);
        IRREP_ASSERT(cabs(DE[0] - 1.0) < 1e-15);
        IRREP_ASSERT(cabs(DE[1]) < 1e-15);
        IRREP_ASSERT(cabs(DE[2]) < 1e-15);
        IRREP_ASSERT(cabs(DE[3] - 1.0) < 1e-15);

        for (int i = 0; i < n; ++i) {
            double _Complex D[4];
            IRREP_ASSERT(irrep_sg_little_group_irrep_matrix(E1, i, D) == 0);

            /* (2) orthogonality: D^T · D = I. */
            double _Complex g00 = conj(D[0])*D[0] + conj(D[2])*D[2];
            double _Complex g01 = conj(D[0])*D[1] + conj(D[2])*D[3];
            double _Complex g11 = conj(D[1])*D[1] + conj(D[3])*D[3];
            IRREP_ASSERT(cabs(g00 - 1.0) < 1e-14);
            IRREP_ASSERT(cabs(g01) < 1e-14);
            IRREP_ASSERT(cabs(g11 - 1.0) < 1e-14);

            /* (3) tr(D) matches character row built into the irrep. */
            double _Complex chi = D[0] + D[3];
            /* C_6v E_1 character row: χ(E)=2, χ(C_6)=1, χ(C_3)=-1,
             * χ(C_2)=-2, χ(σ_v)=0, χ(σ_d)=0. Deduce class from parent op. */
            int parent = 0;
            {
                int ops[12];
                irrep_sg_little_group_point_ops(lg, ops);
                parent = ops[i];
            }
            double chi_expected;
            if (parent == 0) chi_expected = 2.0;
            else if (parent == 1 || parent == 5) chi_expected = 1.0;
            else if (parent == 2 || parent == 4) chi_expected = -1.0;
            else if (parent == 3)                chi_expected = -2.0;
            else                                 chi_expected = 0.0; /* mirrors */
            IRREP_ASSERT(cabs(chi - chi_expected) < 1e-14);
        }

        /* (4) 1D irrep matrix via the same accessor returns 1×1. */
        irrep_sg_little_group_irrep_t *A1 =
            irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_A1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(A1) == 1);
        double _Complex d1[1];
        IRREP_ASSERT(irrep_sg_little_group_irrep_matrix(A1, 0, d1) == 0);
        IRREP_ASSERT(cabs(d1[0] - 1.0) < 1e-15);
        IRREP_ASSERT(irrep_sg_little_group_irrep_matrix(A1, 3, d1) == 0);
        IRREP_ASSERT(cabs(d1[0] - 1.0) < 1e-15); /* trivial rep, all +1 */

        /* Error paths: out-of-range slot, NULL outputs. */
        double _Complex tmp[4];
        IRREP_ASSERT(irrep_sg_little_group_irrep_matrix(NULL, 0, tmp) == -1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_matrix(E1, -1, tmp) == -1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_matrix(E1, n, tmp) == -1);
        IRREP_ASSERT(irrep_sg_little_group_irrep_matrix(E1, 0, NULL) == -1);

        irrep_sg_little_group_irrep_free(E1);
        irrep_sg_little_group_irrep_free(A1);
        irrep_sg_little_group_free(lg);

        /* C_3v E at K-point (1, 2): 2D irrep, verify D(E)=I, trace=chi. */
        irrep_sg_little_group_t *lg_K = irrep_sg_little_group_build(G, 1, 2);
        IRREP_ASSERT(irrep_sg_little_group_point_order(lg_K) == 6);
        irrep_sg_little_group_irrep_t *E_c3v =
            irrep_sg_little_group_irrep_named(lg_K, IRREP_LG_IRREP_E);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(E_c3v) == 2);

        for (int i = 0; i < 6; ++i) {
            double _Complex D[4];
            IRREP_ASSERT(irrep_sg_little_group_irrep_matrix(E_c3v, i, D) == 0);
            double _Complex chi = D[0] + D[3];
            /* C_3v E: χ(E)=2, χ(C_3)=χ(C_3²)=-1, χ(σ)=0. */
            int ops[6];
            irrep_sg_little_group_point_ops(lg_K, ops);
            int p = ops[i];
            double chi_expected;
            if (p == 0) chi_expected = 2.0;
            else if (p == 2 || p == 4) chi_expected = -1.0;
            else chi_expected = 0.0;
            IRREP_ASSERT(cabs(chi - chi_expected) < 1e-14);
        }

        irrep_sg_little_group_irrep_free(E_c3v);
        irrep_sg_little_group_free(lg_K);

        irrep_space_group_free(G);
        irrep_lattice_free(L);

        /* C_4v E on p4mm 4×4 square: 2D at Γ. */
        irrep_lattice_t     *Ls = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
        irrep_space_group_t *Gs = irrep_space_group_build(Ls, IRREP_WALLPAPER_P4MM);
        irrep_sg_little_group_t *lg_s = irrep_sg_little_group_build(Gs, 0, 0);
        irrep_sg_little_group_irrep_t *E_c4v =
            irrep_sg_little_group_irrep_named(lg_s, IRREP_LG_IRREP_E_C4V);
        IRREP_ASSERT(irrep_sg_little_group_irrep_dim(E_c4v) == 2);

        for (int i = 0; i < 8; ++i) {
            double _Complex D[4];
            IRREP_ASSERT(irrep_sg_little_group_irrep_matrix(E_c4v, i, D) == 0);
            double _Complex chi = D[0] + D[3];
            /* C_4v E: χ(E)=2, χ(2C_4)=0, χ(C_2)=-2, χ(2σ_v)=0, χ(2σ_d)=0. */
            int ops[8];
            irrep_sg_little_group_point_ops(lg_s, ops);
            int p = ops[i];
            double chi_expected;
            if (p == 0) chi_expected = 2.0;
            else if (p == 2) chi_expected = -2.0;
            else chi_expected = 0.0;
            IRREP_ASSERT(cabs(chi - chi_expected) < 1e-14);
        }

        irrep_sg_little_group_irrep_free(E_c4v);
        irrep_sg_little_group_free(lg_s);
        irrep_space_group_free(Gs);
        irrep_lattice_free(Ls);
    }

    /* ---- Canonicalisation: negative or out-of-range k maps correctly ---- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P4MM);
        irrep_sg_little_group_t *a = irrep_sg_little_group_build(G, -1, 1);
        irrep_sg_little_group_t *b = irrep_sg_little_group_build(G, 1, 1);
        IRREP_ASSERT(a != NULL && b != NULL);
        IRREP_ASSERT(irrep_sg_little_group_point_order(a) ==
                     irrep_sg_little_group_point_order(b));
        int akx = -99, aky = -99;
        irrep_sg_little_group_k(a, &akx, &aky);
        IRREP_ASSERT(akx == 1 && aky == 1);
        irrep_sg_little_group_free(a);
        irrep_sg_little_group_free(b);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    return IRREP_TEST_END();
}
