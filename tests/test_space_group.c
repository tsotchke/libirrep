/* SPDX-License-Identifier: MIT */
/* Tests for 2D wallpaper-group site permutations.
 *
 * Coverage:
 *   - Build error paths (wrong lattice for p4mm / p6mm, Lx != Ly).
 *   - Order matches expected point_order × num_cells for each supported
 *     (lattice, wallpaper) combination.
 *   - Element 0 is the identity permutation.
 *   - Every permutation is a bijection.
 *   - forward · inverse == identity.
 *   - apply_config reproduces the permutation.
 *   - Research-target 6×6 kagome p6mm = 432 elements × 108 sites.
 *   - Group closure: for any two g1, g2 the composition is in the group.
 *   - Specific geometric sanity checks: C2 on square lattice sends site s
 *     to its antipode through the origin.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "harness.h"

#include <irrep/lattice.h>
#include <irrep/space_group.h>

static int is_bijection_(const int *perm, int n) {
    unsigned char *seen = calloc((size_t)n, 1);
    for (int i = 0; i < n; ++i) {
        if (perm[i] < 0 || perm[i] >= n || seen[perm[i]]) { free(seen); return 0; }
        seen[perm[i]] = 1;
    }
    free(seen);
    return 1;
}

int main(void) {
    IRREP_TEST_START("space_group");

    /* ------------------------------------------------------------------ */
    /* Error paths                                                        */
    /* ------------------------------------------------------------------ */
    IRREP_ASSERT(irrep_space_group_build(NULL, IRREP_WALLPAPER_P1) == NULL);

    irrep_lattice_t *sq = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
    IRREP_ASSERT(irrep_space_group_build(sq, IRREP_WALLPAPER_P6MM) == NULL);

    irrep_lattice_t *sq_rect = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 2);
    IRREP_ASSERT(irrep_space_group_build(sq_rect, IRREP_WALLPAPER_P4MM) == NULL);

    irrep_lattice_t *tri = irrep_lattice_build(IRREP_LATTICE_TRIANGULAR, 4, 4);
    IRREP_ASSERT(irrep_space_group_build(tri, IRREP_WALLPAPER_P4MM) == NULL);

    /* free(NULL) is a no-op */
    irrep_space_group_free(NULL);

    /* ------------------------------------------------------------------ */
    /* p1 on square: translations only                                    */
    /* ------------------------------------------------------------------ */
    irrep_space_group_t *g1 = irrep_space_group_build(sq, IRREP_WALLPAPER_P1);
    IRREP_ASSERT(g1 != NULL);
    IRREP_ASSERT(irrep_space_group_point_order(g1) == 1);
    IRREP_ASSERT(irrep_space_group_num_sites(g1)   == 16);
    IRREP_ASSERT(irrep_space_group_order(g1)       == 16);
    IRREP_ASSERT(irrep_space_group_kind(g1)        == IRREP_WALLPAPER_P1);

    /* Element 0 is identity */
    for (int s = 0; s < 16; ++s) {
        IRREP_ASSERT(irrep_space_group_apply(g1, 0, s) == s);
    }

    /* All 16 translations are bijections */
    int *perm = malloc(sizeof(int) * 16);
    for (int g = 0; g < 16; ++g) {
        irrep_space_group_permutation(g1, g, perm);
        IRREP_ASSERT(is_bijection_(perm, 16));
    }
    free(perm);
    irrep_space_group_free(g1);

    /* ------------------------------------------------------------------ */
    /* p4mm on square 4×4 = 16 cells × 8 point ops = 128 elements         */
    /* ------------------------------------------------------------------ */
    irrep_space_group_t *g4 = irrep_space_group_build(sq, IRREP_WALLPAPER_P4MM);
    IRREP_ASSERT(g4 != NULL);
    IRREP_ASSERT(irrep_space_group_point_order(g4) == 8);
    IRREP_ASSERT(irrep_space_group_order(g4)       == 128);

    /* Element 0 is identity */
    for (int s = 0; s < 16; ++s) {
        IRREP_ASSERT(irrep_space_group_apply(g4, 0, s) == s);
    }

    /* Element 2 is C2 (180° rotation). On the square lattice with origin (0,0),
     * C2 sends site (ix, iy) → (-ix, -iy) mod (Lx, Ly).  For Lx=Ly=4: s=0
     * stays at 0, s=1 (ix=1, iy=0) maps to ix=-1 mod 4 = 3, iy=0 → site 3. */
    IRREP_ASSERT(irrep_space_group_apply(g4, 2, 0) == 0);
    IRREP_ASSERT(irrep_space_group_apply(g4, 2, 1) == 3);

    /* Every one of the 128 permutations is a bijection */
    int *q4 = malloc(sizeof(int) * 16);
    for (int g = 0; g < 128; ++g) {
        irrep_space_group_permutation(g4, g, q4);
        IRREP_ASSERT(is_bijection_(q4, 16));
    }

    /* Inverse composes to identity */
    int *q4_inv = malloc(sizeof(int) * 16);
    for (int g = 0; g < 128; g += 17) {
        irrep_space_group_permutation(g4, g, q4);
        irrep_space_group_permutation_inverse(g4, g, q4_inv);
        int ok = 1;
        for (int s = 0; s < 16 && ok; ++s) if (q4_inv[q4[s]] != s) ok = 0;
        IRREP_ASSERT(ok);
    }
    free(q4); free(q4_inv);

    /* apply_config matches the permutation (pullback) */
    double cfg_in[16], cfg_out[16];
    for (int i = 0; i < 16; ++i) cfg_in[i] = (double)(i * 7 % 13);
    int *perm_check = malloc(sizeof(int) * 16);
    for (int g = 0; g < 128; g += 11) {
        irrep_space_group_apply_config(g4, g, cfg_in, cfg_out);
        irrep_space_group_permutation_inverse(g4, g, perm_check);
        int ok = 1;
        for (int s = 0; s < 16 && ok; ++s) {
            if (cfg_out[s] != cfg_in[perm_check[s]]) ok = 0;
        }
        IRREP_ASSERT(ok);
    }
    free(perm_check);
    irrep_space_group_free(g4);

    /* ------------------------------------------------------------------ */
    /* p6mm on triangular 4×4: 4×4 cells × 12 ops = 192 elements          */
    /* ------------------------------------------------------------------ */
    irrep_space_group_t *g6t = irrep_space_group_build(tri, IRREP_WALLPAPER_P6MM);
    IRREP_ASSERT(g6t != NULL);
    IRREP_ASSERT(irrep_space_group_point_order(g6t) == 12);
    IRREP_ASSERT(irrep_space_group_order(g6t)       == 192);

    int *qt = malloc(sizeof(int) * 16);
    for (int g = 0; g < 192; g += 13) {
        irrep_space_group_permutation(g6t, g, qt);
        IRREP_ASSERT(is_bijection_(qt, 16));
    }
    free(qt);
    irrep_space_group_free(g6t);

    /* ------------------------------------------------------------------ */
    /* p6mm on honeycomb 3×3: 18 sites × 12 ops × 9 cells = 108 elements  */
    /* ------------------------------------------------------------------ */
    irrep_lattice_t *hc = irrep_lattice_build(IRREP_LATTICE_HONEYCOMB, 3, 3);
    IRREP_ASSERT(hc != NULL);
    irrep_space_group_t *g6h = irrep_space_group_build(hc, IRREP_WALLPAPER_P6MM);
    IRREP_ASSERT(g6h != NULL);
    IRREP_ASSERT(irrep_space_group_num_sites(g6h) == 18);
    IRREP_ASSERT(irrep_space_group_order(g6h)     == 108);

    int *qh = malloc(sizeof(int) * 18);
    for (int g = 0; g < 108; g += 7) {
        irrep_space_group_permutation(g6h, g, qh);
        IRREP_ASSERT(is_bijection_(qh, 18));
    }
    free(qh);
    irrep_space_group_free(g6h);
    irrep_lattice_free(hc);

    /* ------------------------------------------------------------------ */
    /* Research target: p6mm on 6×6 kagome — 108 sites × 432 elements      */
    /* ------------------------------------------------------------------ */
    irrep_lattice_t *kg = irrep_lattice_build(IRREP_LATTICE_KAGOME, 6, 6);
    IRREP_ASSERT(kg != NULL);
    irrep_space_group_t *g6k = irrep_space_group_build(kg, IRREP_WALLPAPER_P6MM);
    IRREP_ASSERT(g6k != NULL);
    IRREP_ASSERT(irrep_space_group_num_sites(g6k)   == 108);
    IRREP_ASSERT(irrep_space_group_point_order(g6k) ==  12);
    IRREP_ASSERT(irrep_space_group_order(g6k)       == 432);

    /* Every permutation is a bijection (sample 1/8 of them) */
    int *qk = malloc(sizeof(int) * 108);
    for (int g = 0; g < 432; g += 8) {
        irrep_space_group_permutation(g6k, g, qk);
        IRREP_ASSERT(is_bijection_(qk, 108));
    }

    /* Element 0 (E, translation 0) is identity */
    for (int s = 0; s < 108; s += 11) {
        IRREP_ASSERT(irrep_space_group_apply(g6k, 0, s) == s);
    }

    /* Inverse round-trip for a few elements */
    int *qk_inv = malloc(sizeof(int) * 108);
    for (int g = 0; g < 432; g += 37) {
        irrep_space_group_permutation(g6k, g, qk);
        irrep_space_group_permutation_inverse(g6k, g, qk_inv);
        int ok = 1;
        for (int s = 0; s < 108 && ok; ++s) if (qk_inv[qk[s]] != s) ok = 0;
        IRREP_ASSERT(ok);
    }
    free(qk); free(qk_inv);

    /* Group closure check: composition of two elements is itself
     * a permutation found in the table.  We verify that g1 · g2
     * equals *some* group element g3. */
    int *pa = malloc(sizeof(int) * 108);
    int *pb = malloc(sizeof(int) * 108);
    int *pab = malloc(sizeof(int) * 108);
    /* Sample: g1 = 17 (C3 + translation), g2 = 139 */
    irrep_space_group_permutation(g6k,  17, pa);
    irrep_space_group_permutation(g6k, 139, pb);
    for (int s = 0; s < 108; ++s) pab[s] = pa[pb[s]];
    int found = 0;
    int *cand = malloc(sizeof(int) * 108);
    for (int g = 0; g < 432; ++g) {
        irrep_space_group_permutation(g6k, g, cand);
        int match = 1;
        for (int s = 0; s < 108; ++s) if (cand[s] != pab[s]) { match = 0; break; }
        if (match) { found = 1; break; }
    }
    IRREP_ASSERT(found);
    free(pa); free(pb); free(pab); free(cand);

    irrep_space_group_free(g6k);
    irrep_lattice_free(kg);

    irrep_lattice_free(sq);
    irrep_lattice_free(sq_rect);
    irrep_lattice_free(tri);

    return IRREP_TEST_END();
}
