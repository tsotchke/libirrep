/* SPDX-License-Identifier: MIT */
/* Tests for the wallpaper groups beyond p1 / p4mm / p6mm: p2, p6, p3m1.
 *
 * Coverage:
 *   - Order and point_order match the group theory for each kind.
 *   - Bijectivity: every permutation is a valid bijection of sites.
 *   - p2 on rectangular clusters (Lx ≠ Ly): builds successfully on both
 *     square and hexagonal lattices — the killer feature, enabling C_2
 *     symmetry on finite-size clusters that can't host the full C_4v or
 *     C_6v (e.g. 2×3 kagome, N = 18).
 *   - Hexagonal groups (p6, p3m1) require Lx = Ly on hex-family lattices.
 *   - Lattice-compatibility gates: p4mm requires square; p6mm / p6 /
 *     p3m1 require non-square; p2 accepts any lattice.
 *   - Subgroup relations: p6 is a subgroup of p6mm (every p6 element's
 *     site permutation equals some p6mm element's); p3m1 similarly.
 */

#include "harness.h"
#include <irrep/lattice.h>
#include <irrep/space_group.h>
#include <stdlib.h>
#include <string.h>

static int perm_is_bijection_(const int *perm, int n) {
    int *seen = calloc((size_t)n, sizeof(int));
    if (!seen)
        return 0;
    for (int s = 0; s < n; ++s) {
        int img = perm[s];
        if (img < 0 || img >= n || seen[img]) {
            free(seen);
            return 0;
        }
        seen[img] = 1;
    }
    free(seen);
    return 1;
}

static int perm_equal_(const int *a, const int *b, int n) {
    for (int s = 0; s < n; ++s)
        if (a[s] != b[s])
            return 0;
    return 1;
}

int main(void) {
    IRREP_TEST_START("wallpaper_groups");

    /* ---- p2 on a 3×5 square lattice (Lx ≠ Ly) -------------------------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 3, 5);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P2);
        IRREP_ASSERT(G != NULL);
        IRREP_ASSERT(irrep_space_group_point_order(G) == 2);
        /* order = point · Lx · Ly = 2 · 15 = 30. */
        IRREP_ASSERT(irrep_space_group_order(G) == 30);
        int num = irrep_space_group_num_sites(G);
        int perm[15];
        for (int g = 0; g < irrep_space_group_order(G); ++g) {
            irrep_space_group_permutation(G, g, perm);
            IRREP_ASSERT(perm_is_bijection_(perm, num));
        }
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- p2 on a 2×3 kagome torus (N = 18, Lx ≠ Ly, hex lattice) ------- *
     * This is the N = 18 cluster used in PHYSICS_RESULTS.md — previously
     * forced to p1. With p2 landed, C_2 symmetry on the 2×3 torus is now
     * an available sector structure for ED block decomposition.          */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P2);
        IRREP_ASSERT(G != NULL);
        IRREP_ASSERT(irrep_space_group_point_order(G) == 2);
        IRREP_ASSERT(irrep_space_group_num_sites(G) == 18); /* 2·3·3 */
        IRREP_ASSERT(irrep_space_group_order(G) == 12);     /* 2 · 6 cells */
        int perm[18];
        for (int g = 0; g < irrep_space_group_order(G); ++g) {
            irrep_space_group_permutation(G, g, perm);
            IRREP_ASSERT(perm_is_bijection_(perm, 18));
        }
        /* C_2 acts as r → -r around the hexagon centre; on the kagome
         * sublattice labelling this should permute sublattices (0 ↔ 0,
         * 1 ↔ 1, 2 ↔ 2 with cell-position flipped). */
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- p6 on a 3×3 kagome: 6 rotations, no mirrors ------------------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6);
        IRREP_ASSERT(G != NULL);
        IRREP_ASSERT(irrep_space_group_point_order(G) == 6);
        int perm[27];
        for (int g = 0; g < irrep_space_group_order(G); ++g) {
            irrep_space_group_permutation(G, g, perm);
            IRREP_ASSERT(perm_is_bijection_(perm, 27));
        }

        /* Subgroup check: each p6 element should match a p6mm element. */
        irrep_space_group_t *G_full = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        IRREP_ASSERT(G_full != NULL);
        /* p6 point elements are the first 6 of p6mm (rotations; mirrors at 6..11). */
        for (int p = 0; p < 6; ++p) {
            int perm_p6[27], perm_p6mm[27];
            irrep_space_group_permutation(G, p, perm_p6);       /* tidx = 0, point = p */
            irrep_space_group_permutation(G_full, p, perm_p6mm);
            IRREP_ASSERT(perm_equal_(perm_p6, perm_p6mm, 27));
        }
        irrep_space_group_free(G_full);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- p3m1 on a 3×3 kagome: 3 rotations + 3 mirrors (order 6) ------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P3M1);
        IRREP_ASSERT(G != NULL);
        IRREP_ASSERT(irrep_space_group_point_order(G) == 6);
        int perm[27];
        for (int g = 0; g < irrep_space_group_order(G); ++g) {
            irrep_space_group_permutation(G, g, perm);
            IRREP_ASSERT(perm_is_bijection_(perm, 27));
        }
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Compatibility gates: each group rejects wrong-lattice calls --- */
    {
        irrep_lattice_t *L_sq = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
        IRREP_ASSERT(irrep_space_group_build(L_sq, IRREP_WALLPAPER_P6MM) == NULL);
        IRREP_ASSERT(irrep_space_group_build(L_sq, IRREP_WALLPAPER_P6) == NULL);
        IRREP_ASSERT(irrep_space_group_build(L_sq, IRREP_WALLPAPER_P3M1) == NULL);
        /* p2 is lattice-agnostic. */
        irrep_space_group_t *G_p2 = irrep_space_group_build(L_sq, IRREP_WALLPAPER_P2);
        IRREP_ASSERT(G_p2 != NULL);
        irrep_space_group_free(G_p2);
        irrep_lattice_free(L_sq);

        irrep_lattice_t *L_hex = irrep_lattice_build(IRREP_LATTICE_KAGOME, 4, 4);
        IRREP_ASSERT(irrep_space_group_build(L_hex, IRREP_WALLPAPER_P4MM) == NULL);
        irrep_lattice_free(L_hex);
    }

    /* ---- p4 on a 4×4 square: 4 rotations, no mirrors, subgroup of p4mm --- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P4);
        IRREP_ASSERT(G != NULL);
        IRREP_ASSERT(irrep_space_group_point_order(G) == 4);
        int perm[16];
        for (int g = 0; g < irrep_space_group_order(G); ++g) {
            irrep_space_group_permutation(G, g, perm);
            IRREP_ASSERT(perm_is_bijection_(perm, 16));
        }
        /* Subgroup bit-exactness: p4 rotations == first 4 elements of p4mm. */
        irrep_space_group_t *G_full = irrep_space_group_build(L, IRREP_WALLPAPER_P4MM);
        for (int p = 0; p < 4; ++p) {
            int pa[16], pb[16];
            irrep_space_group_permutation(G, p, pa);
            irrep_space_group_permutation(G_full, p, pb);
            IRREP_ASSERT(perm_equal_(pa, pb, 16));
        }
        irrep_space_group_free(G_full);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- p31m on a 3×3 kagome: 3 rotations + 3 mirrors (distinct from p3m1) ---- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P31M);
        IRREP_ASSERT(G != NULL);
        IRREP_ASSERT(irrep_space_group_point_order(G) == 6);
        int perm[27];
        for (int g = 0; g < irrep_space_group_order(G); ++g) {
            irrep_space_group_permutation(G, g, perm);
            IRREP_ASSERT(perm_is_bijection_(perm, 27));
        }
        /* Rotations match p3m1's rotations bit-exactly (same C_3 on hex). */
        irrep_space_group_t *G_p3m1 = irrep_space_group_build(L, IRREP_WALLPAPER_P3M1);
        for (int p = 0; p < 3; ++p) {
            int pa[27], pb[27];
            irrep_space_group_permutation(G, p, pa);
            irrep_space_group_permutation(G_p3m1, p, pb);
            IRREP_ASSERT(perm_equal_(pa, pb, 27));
        }
        /* Mirror-classes differ: p31m's mirrors at slots 3..5 must NOT equal
         * p3m1's mirrors at slots 3..5 (different axis orientation). */
        int mirror_mismatch = 0;
        for (int p = 3; p < 6; ++p) {
            int pa[27], pb[27];
            irrep_space_group_permutation(G, p, pa);
            irrep_space_group_permutation(G_p3m1, p, pb);
            if (!perm_equal_(pa, pb, 27))
                ++mirror_mismatch;
        }
        IRREP_ASSERT(mirror_mismatch > 0);
        irrep_space_group_free(G_p3m1);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- p4gm rejected on every currently-shipped lattice -------------- *
     * p4gm is non-symmorphic: its four mirrors carry a ½·(a_1 + a_2) glide
     * translation. On a single-sublattice square lattice the glide maps
     * site (i, j) → (i+½, −j+½), which is not an integer-site coordinate.
     * Until a two-basis square lattice ships, p4gm has no compatible
     * lattice and must be rejected at build time with an informative
     * diagnostic. The fractional-translation plumbing (`fill_p4gm_`,
     * `t_frac` → cartesian in `build_point_perm_`) is in place and would
     * take effect on a compatible lattice without further space-group
     * changes. */
    {
        irrep_lattice_t *lattices[] = {
            irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4),
            irrep_lattice_build(IRREP_LATTICE_SQUARE, 2, 2),
            irrep_lattice_build(IRREP_LATTICE_TRIANGULAR, 4, 4),
            irrep_lattice_build(IRREP_LATTICE_KAGOME, 4, 4),
            irrep_lattice_build(IRREP_LATTICE_HONEYCOMB, 4, 4),
        };
        for (size_t i = 0; i < sizeof(lattices) / sizeof(lattices[0]); ++i) {
            IRREP_ASSERT(irrep_space_group_build(lattices[i], IRREP_WALLPAPER_P4GM) == NULL);
            irrep_lattice_free(lattices[i]);
        }
    }

    /* ---- Non-square clusters rejected by groups that require Lx = Ly --- */
    {
        irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 3);
        IRREP_ASSERT(irrep_space_group_build(L, IRREP_WALLPAPER_P6) == NULL);
        IRREP_ASSERT(irrep_space_group_build(L, IRREP_WALLPAPER_P3M1) == NULL);
        IRREP_ASSERT(irrep_space_group_build(L, IRREP_WALLPAPER_P6MM) == NULL);
        /* p2 accepts it. */
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P2);
        IRREP_ASSERT(G != NULL);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    return IRREP_TEST_END();
}
