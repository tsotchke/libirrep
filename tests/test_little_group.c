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
        int                  order = irrep_space_group_order(G); /* 12 · 4 = 48 */

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
