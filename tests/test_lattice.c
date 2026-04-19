/* SPDX-License-Identifier: MIT */
/* Tests for the 2D lattice module.
 *
 * Coverage:
 *   - Builder error paths (bad Lx/Ly, bad kind).
 *   - Site count, cells, sublattices per lattice family.
 *   - Primitive-vector values match the documented conventions.
 *   - Reciprocal-vector duality: a_i · b_j = 2π δ_{ij}.
 *   - NN bond count matches site-count × NN-coordination / 2 on
 *     lattices large enough that PBC does not collapse any bond.
 *   - Every emitted bond pair lies at the expected NN distance (1.0) or
 *     NNN distance (lattice-dependent) up to 1e-12.
 *   - Translation composes correctly and obeys PBC wrap.
 *   - k-grid generates Lx*Ly unique momenta on the square lattice.
 *   - Targeted 6×6×3 = 108-site kagome cluster that the 1.3 scope
 *     pins as the primary NQS target.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "harness.h"

#include <irrep/lattice.h>

static double dist2_(const double a[2], const double b[2]) {
    double dx = a[0] - b[0], dy = a[1] - b[1];
    return dx*dx + dy*dy;
}

static void check_all_bonds_at_(const irrep_lattice_t *L,
                                int nb, const int *ii, const int *jj,
                                double expected_len, double tol) {
    for (int k = 0; k < nb; ++k) {
        double pi[2], pj[2];
        irrep_lattice_site_position(L, ii[k], pi);
        irrep_lattice_site_position(L, jj[k], pj);
        double d = sqrt(dist2_(pi, pj));
        /* One bond may span the PBC wrap; the geometry check should measure
         * the minimum-image distance, not the raw cartesian one. Fix up by
         * taking the lattice-primitive-vector-modulo nearest image. */
        double a1[2], a2[2];
        irrep_lattice_primitive_vectors(L, a1, a2);
        int Lx = irrep_lattice_Lx(L), Ly = irrep_lattice_Ly(L);
        double min_d = d;
        for (int tx = -1; tx <= 1; ++tx) {
            for (int ty = -1; ty <= 1; ++ty) {
                double shift[2] = { tx * Lx * a1[0] + ty * Ly * a2[0],
                                    tx * Lx * a1[1] + ty * Ly * a2[1] };
                double p[2] = { pj[0] + shift[0], pj[1] + shift[1] };
                double trial = sqrt(dist2_(pi, p));
                if (trial < min_d) min_d = trial;
            }
        }
        IRREP_ASSERT_NEAR(min_d, expected_len, tol);
    }
}

static void check_bonds_sorted_unique_(int nb, const int *ii, const int *jj) {
    for (int k = 0; k < nb; ++k) IRREP_ASSERT(ii[k] < jj[k]);
    for (int k = 1; k < nb; ++k) {
        IRREP_ASSERT(ii[k-1] < ii[k] || (ii[k-1] == ii[k] && jj[k-1] < jj[k]));
    }
}

int main(void) {
    IRREP_TEST_START("lattice");

    /* ------------------------------------------------------------------ */
    /* Error paths                                                        */
    /* ------------------------------------------------------------------ */
    IRREP_ASSERT(irrep_lattice_build(IRREP_LATTICE_SQUARE, 1, 4) == NULL);
    IRREP_ASSERT(irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 1) == NULL);
    IRREP_ASSERT(irrep_lattice_build(IRREP_LATTICE_SQUARE, 0, 0) == NULL);
    IRREP_ASSERT(irrep_lattice_build((irrep_lattice_kind_t)99, 4, 4) == NULL);

    /* free(NULL) is a no-op */
    irrep_lattice_free(NULL);

    /* ------------------------------------------------------------------ */
    /* SQUARE 4×4                                                         */
    /* ------------------------------------------------------------------ */
    irrep_lattice_t *sq = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
    IRREP_ASSERT(sq != NULL);
    IRREP_ASSERT(irrep_lattice_num_sites(sq)       == 16);
    IRREP_ASSERT(irrep_lattice_num_cells(sq)       == 16);
    IRREP_ASSERT(irrep_lattice_sites_per_cell(sq)  == 1);
    IRREP_ASSERT(irrep_lattice_kind(sq)            == IRREP_LATTICE_SQUARE);

    /* 4 NN per site × 16 sites / 2 = 32 */
    IRREP_ASSERT(irrep_lattice_num_bonds_nn(sq)    == 32);
    /* 4 NNN per site × 16 sites / 2 = 32 */
    IRREP_ASSERT(irrep_lattice_num_bonds_nnn(sq)   == 32);

    double a1[2], a2[2], b1[2], b2[2];
    irrep_lattice_primitive_vectors(sq, a1, a2);
    IRREP_ASSERT_NEAR(a1[0], 1.0, 1e-15);
    IRREP_ASSERT_NEAR(a1[1], 0.0, 1e-15);
    IRREP_ASSERT_NEAR(a2[0], 0.0, 1e-15);
    IRREP_ASSERT_NEAR(a2[1], 1.0, 1e-15);

    irrep_lattice_reciprocal_vectors(sq, b1, b2);
    /* a_i · b_j = 2π δ_{ij} */
    IRREP_ASSERT_NEAR(a1[0]*b1[0] + a1[1]*b1[1], 2.0*M_PI, 1e-12);
    IRREP_ASSERT_NEAR(a2[0]*b2[0] + a2[1]*b2[1], 2.0*M_PI, 1e-12);
    IRREP_ASSERT_NEAR(a1[0]*b2[0] + a1[1]*b2[1], 0.0,      1e-12);
    IRREP_ASSERT_NEAR(a2[0]*b1[0] + a2[1]*b1[1], 0.0,      1e-12);

    int *ii = malloc(sizeof(int) * 64);
    int *jj = malloc(sizeof(int) * 64);

    irrep_lattice_fill_bonds_nn(sq, ii, jj);
    check_bonds_sorted_unique_(32, ii, jj);
    check_all_bonds_at_(sq, 32, ii, jj, 1.0, 1e-12);

    irrep_lattice_fill_bonds_nnn(sq, ii, jj);
    check_bonds_sorted_unique_(32, ii, jj);
    check_all_bonds_at_(sq, 32, ii, jj, sqrt(2.0), 1e-12);

    /* Translation */
    int s0 = 0;
    int s1 = irrep_lattice_translate(sq, s0, 1, 0);
    IRREP_ASSERT(s1 != -1 && s1 != s0);
    int s_back = irrep_lattice_translate(sq, s1, -1, 0);
    IRREP_ASSERT(s_back == s0);
    /* Full loop around PBC: +Lx, 0 */
    IRREP_ASSERT(irrep_lattice_translate(sq, s0, 4, 0) == s0);

    /* ------------------------------------------------------------------ */
    /* TRIANGULAR 4×4                                                     */
    /* ------------------------------------------------------------------ */
    irrep_lattice_t *tri = irrep_lattice_build(IRREP_LATTICE_TRIANGULAR, 4, 4);
    IRREP_ASSERT(tri != NULL);
    IRREP_ASSERT(irrep_lattice_num_sites(tri) == 16);
    /* 6 NN per site × 16 / 2 = 48 */
    IRREP_ASSERT(irrep_lattice_num_bonds_nn(tri) == 48);
    /* 6 NNN per site × 16 / 2 = 48 */
    IRREP_ASSERT(irrep_lattice_num_bonds_nnn(tri) == 48);

    int *ii_t = malloc(sizeof(int) * 96);
    int *jj_t = malloc(sizeof(int) * 96);
    irrep_lattice_fill_bonds_nn(tri, ii_t, jj_t);
    check_bonds_sorted_unique_(48, ii_t, jj_t);
    check_all_bonds_at_(tri, 48, ii_t, jj_t, 1.0, 1e-12);
    irrep_lattice_fill_bonds_nnn(tri, ii_t, jj_t);
    check_bonds_sorted_unique_(48, ii_t, jj_t);
    check_all_bonds_at_(tri, 48, ii_t, jj_t, sqrt(3.0), 1e-12);
    free(ii_t); free(jj_t);

    /* ------------------------------------------------------------------ */
    /* HONEYCOMB 3×3                                                      */
    /* ------------------------------------------------------------------ */
    irrep_lattice_t *hc = irrep_lattice_build(IRREP_LATTICE_HONEYCOMB, 3, 3);
    IRREP_ASSERT(hc != NULL);
    IRREP_ASSERT(irrep_lattice_num_sites(hc)      == 18);
    IRREP_ASSERT(irrep_lattice_sites_per_cell(hc) ==  2);
    /* 3 NN per site × 18 / 2 = 27 */
    IRREP_ASSERT(irrep_lattice_num_bonds_nn(hc)   == 27);
    /* 6 NNN per site × 18 / 2 = 54 */
    IRREP_ASSERT(irrep_lattice_num_bonds_nnn(hc)  == 54);

    int *ii_h = malloc(sizeof(int) * 108);
    int *jj_h = malloc(sizeof(int) * 108);
    irrep_lattice_fill_bonds_nn(hc, ii_h, jj_h);
    check_bonds_sorted_unique_(27, ii_h, jj_h);
    check_all_bonds_at_(hc, 27, ii_h, jj_h, 1.0, 1e-12);
    irrep_lattice_fill_bonds_nnn(hc, ii_h, jj_h);
    check_bonds_sorted_unique_(54, ii_h, jj_h);
    check_all_bonds_at_(hc, 54, ii_h, jj_h, sqrt(3.0), 1e-12);
    free(ii_h); free(jj_h);

    /* All NN bonds must cross sublattices on honeycomb */
    irrep_lattice_fill_bonds_nn(hc, NULL, NULL);  /* smoke: NULL fill ok */
    int n_nn = irrep_lattice_num_bonds_nn(hc);
    int *ii_hnn = malloc(sizeof(int)*n_nn), *jj_hnn = malloc(sizeof(int)*n_nn);
    irrep_lattice_fill_bonds_nn(hc, ii_hnn, jj_hnn);
    for (int k = 0; k < n_nn; ++k) {
        int si = irrep_lattice_sublattice_of(hc, ii_hnn[k]);
        int sj = irrep_lattice_sublattice_of(hc, jj_hnn[k]);
        IRREP_ASSERT(si != sj);
    }
    free(ii_hnn); free(jj_hnn);

    /* ------------------------------------------------------------------ */
    /* KAGOME — 6×6 cluster.  6×6 = 108 sites.                    */
    /* ------------------------------------------------------------------ */
    irrep_lattice_t *kg = irrep_lattice_build(IRREP_LATTICE_KAGOME, 6, 6);
    IRREP_ASSERT(kg != NULL);
    IRREP_ASSERT(irrep_lattice_num_sites(kg)      == 108);
    IRREP_ASSERT(irrep_lattice_sites_per_cell(kg) ==   3);
    IRREP_ASSERT(irrep_lattice_num_cells(kg)      ==  36);

    /* NN: 4 per site × 108 / 2 = 216 */
    IRREP_ASSERT(irrep_lattice_num_bonds_nn(kg)   == 216);
    /* NNN (across hexagon, mixed sublattice): 4 per site × 108 / 2 = 216 */
    IRREP_ASSERT(irrep_lattice_num_bonds_nnn(kg)  == 216);

    int *ii_k = malloc(sizeof(int) * 432);
    int *jj_k = malloc(sizeof(int) * 432);
    irrep_lattice_fill_bonds_nn(kg, ii_k, jj_k);
    check_bonds_sorted_unique_(216, ii_k, jj_k);
    check_all_bonds_at_(kg, 216, ii_k, jj_k, 1.0, 1e-12);
    irrep_lattice_fill_bonds_nnn(kg, ii_k, jj_k);
    check_bonds_sorted_unique_(216, ii_k, jj_k);
    check_all_bonds_at_(kg, 216, ii_k, jj_k, sqrt(3.0), 1e-12);

    /* Every kagome NN bond is between adjacent corners of a triangle, so
     * the two endpoints must always be on different sublattices. */
    irrep_lattice_fill_bonds_nn(kg, ii_k, jj_k);
    for (int k = 0; k < 216; ++k) {
        int si = irrep_lattice_sublattice_of(kg, ii_k[k]);
        int sj = irrep_lattice_sublattice_of(kg, jj_k[k]);
        IRREP_ASSERT(si != sj);
    }

    /* Every site has exactly 4 NN on 6×6 kagome */
    int *deg = calloc(108, sizeof(int));
    for (int k = 0; k < 216; ++k) { deg[ii_k[k]]++; deg[jj_k[k]]++; }
    int all4 = 1;
    for (int s = 0; s < 108; ++s) if (deg[s] != 4) all4 = 0;
    IRREP_ASSERT(all4);
    free(deg);
    free(ii_k); free(jj_k);

    /* Cell / index round-trip */
    for (int s = 0; s < 108; s += 17) {
        int ix, iy;
        IRREP_ASSERT(irrep_lattice_cell_of(kg, s, &ix, &iy) == IRREP_OK);
        int sub = irrep_lattice_sublattice_of(kg, s);
        IRREP_ASSERT(irrep_lattice_site_index(kg, ix, iy, sub) == s);
    }

    /* Translation by (Lx, 0) or (0, Ly) is identity under PBC */
    for (int s = 0; s < 108; s += 11) {
        IRREP_ASSERT(irrep_lattice_translate(kg, s, 6, 0) == s);
        IRREP_ASSERT(irrep_lattice_translate(kg, s, 0, 6) == s);
        IRREP_ASSERT(irrep_lattice_translate(kg, s, 6, 6) == s);
    }

    /* k-grid: 36 unique k-points on 6×6 kagome */
    double *kx = malloc(sizeof(double) * 36);
    double *ky = malloc(sizeof(double) * 36);
    irrep_lattice_k_grid(kg, kx, ky);
    /* k=0 lives at index 0 */
    IRREP_ASSERT_NEAR(kx[0], 0.0, 1e-15);
    IRREP_ASSERT_NEAR(ky[0], 0.0, 1e-15);
    /* All 36 k-points distinct (in the first BZ, equivalent points
     * are separated by reciprocal primitive vectors) */
    int unique = 1;
    for (int a = 0; a < 36 && unique; ++a) {
        for (int b = a + 1; b < 36 && unique; ++b) {
            double dkx = kx[a]-kx[b], dky = ky[a]-ky[b];
            if (dkx*dkx + dky*dky < 1e-16) unique = 0;
        }
    }
    IRREP_ASSERT(unique);
    free(kx); free(ky);

    /* ------------------------------------------------------------------ */
    /* Small-cluster sanity: 2×2 square — PBC collapses the system to a
     * square with 2 distinct neighbours per site.  This exercises the
     * dedup pass.                                                        */
    /* ------------------------------------------------------------------ */
    irrep_lattice_t *sq22 = irrep_lattice_build(IRREP_LATTICE_SQUARE, 2, 2);
    IRREP_ASSERT(sq22 != NULL);
    /* Each site has only 2 distinct NN (the left = right wrap, etc.) */
    IRREP_ASSERT(irrep_lattice_num_bonds_nn(sq22) == 4);
    irrep_lattice_free(sq22);

    /* Cleanup */
    free(ii); free(jj);
    irrep_lattice_free(sq);
    irrep_lattice_free(tri);
    irrep_lattice_free(hc);
    irrep_lattice_free(kg);

    return IRREP_TEST_END();
}
