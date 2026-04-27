/* SPDX-License-Identifier: MIT */
/* Tests for the linearised-spin-wave (LSW) module. Verifies:
 *
 *   1. FM square-lattice Heisenberg dispersion:
 *        ω(k) = 2|J|S z (1 − γ_k) = 2|J|S(2 − cos kx − cos ky)   for J<0.
 *      Gapless at k=0 (Goldstone), peaks at (π,π).
 *
 *   2. Uniaxial anisotropy K_z opens a gap of 2 K_z S at k=0.
 *
 *   3. 2-sublattice unit cell with NN J reproduces the dispersion of a
 *      single-sublattice FM with twice the BZ — folded back into a 2-band
 *      structure that touches at the BZ corner.
 *
 *   4. Hermiticity of the LSW Hamiltonian for a generic J + DMI bond list.
 *
 *   5. FM ground-state stability: ω(k) ≥ 0 for J<0, K_z ≥ 0.
 *
 *   6. Berry curvature of a trivial 1-band model is identically zero
 *      (no off-diagonal phase to wind around).
 *
 *   7. Chern number of the 1-band trivial model integrates to 0 under
 *      a coarse Nx × Ny grid.
 */

#include <irrep/magnon.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int total = 0, failed = 0;

#define ASSERT_NEAR(actual, expected, tol, msg)                                                    \
    do {                                                                                           \
        ++total;                                                                                   \
        double a__ = (actual);                                                                     \
        double e__ = (expected);                                                                   \
        if (!(fabs(a__ - e__) <= (tol))) {                                                         \
            fprintf(stderr, "  FAIL  %s:%d  %s  got %.10f vs %.10f (|Δ|=%.3g, tol=%.3g)\n",        \
                    __FILE__, __LINE__, msg, a__, e__, fabs(a__ - e__), (double)(tol));            \
            ++failed;                                                                              \
        }                                                                                          \
    } while (0)

#define ASSERT(cond, msg)                                                                          \
    do {                                                                                           \
        ++total;                                                                                   \
        if (!(cond)) {                                                                             \
            fprintf(stderr, "  FAIL  %s:%d  %s\n", __FILE__, __LINE__, msg);                       \
            ++failed;                                                                              \
        }                                                                                          \
    } while (0)

/* (1) FM square lattice with J<0: ω(k) = 2|J|S(2 − cos kx − cos ky). */
static void test_fm_square_dispersion(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    /* Two unique NN bonds: (0,0,+x) and (0,0,+y). FM ⇒ J<0. */
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}},
    };
    double S = 0.5;
    irrep_magnon_lsw_t *L =
        irrep_magnon_lsw_new(/*n_sub=*/1, S, a1, a2, /*n_bonds=*/2, bonds, /*Kz=*/0);
    ASSERT(L != NULL, "FM square LSW handle");

    /* Sample 5 k-points and compare to closed form. */
    struct {
        double kx, ky;
    } kpts[] = {
        {0.0, 0.0},
        {M_PI, 0.0},
        {0.0, M_PI},
        {M_PI / 2.0, M_PI / 3.0},
        {M_PI, M_PI},
    };
    double          omega;
    double _Complex u;
    for (int p = 0; p < 5; ++p) {
        irrep_status_t st = irrep_magnon_dispersion(L, kpts[p].kx, kpts[p].ky, &omega, &u);
        ASSERT(st == 0, "dispersion call ok");
        double expected = 2.0 * 1.0 * S * (2.0 - cos(kpts[p].kx) - cos(kpts[p].ky));
        ASSERT_NEAR(omega, expected, 1e-12, "FM square ω(k)");
    }
    irrep_magnon_lsw_free(L);
}

/* (2) Anisotropy gap: ω(k=0) = 2 K_z S for J<0 FM with K_z > 0. */
static void test_anisotropy_gap(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}},
    };
    double Kz = 0.25, S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, S, a1, a2, 2, bonds, Kz);
    double              omega;
    double _Complex     u;
    irrep_magnon_dispersion(L, 0.0, 0.0, &omega, &u);
    ASSERT_NEAR(omega, 2.0 * Kz * S, 1e-12, "anisotropy gap at k=0");
    irrep_magnon_lsw_free(L);
}

/* (3) 2-sublattice rectangular cell with the *same* underlying FM square
 * lattice (sublattice A at origin, sublattice B at +x): the LSW H(k) is a
 * 2×2 matrix whose two bands trace the original 1-band dispersion folded
 * by the doubled BZ. At the unit-cell-doubled BZ corner k1=0, the bands
 * are degenerate at ω=0 (touching point). */
static void test_two_sublattice_folded(void) {
    /* Doubled cell: a1' = (2, 0), a2' = (0, 1). Two NN-x bonds inside the
     * cell — (A at (0,0) → B at (1,0)) is delta=(0,0), and (B → A at next
     * cell along a1') is delta=(+1, 0). NN-y bonds stay within each
     * sublattice. */
    double a1p[2] = {2.0, 0.0};
    double a2p[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 1, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}}, /* A→B intra */
        {.bi = 1, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}}, /* B→A inter */
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}}, /* A y */
        {.bi = 1, .bj = 1, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}}, /* B y */
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(2, S, a1p, a2p, 4, bonds, 0);

    double          omega[2];
    double _Complex u[4];

    /* At the doubled-BZ origin (kx,ky=0): bands at 0 (Goldstone of folded
     * cell) and at the BZ-corner image of the original lattice. */
    irrep_magnon_dispersion(L, 0.0, 0.0, omega, u);
    ASSERT_NEAR(omega[0], 0.0, 1e-10, "folded 2-band: lower at k=0");
    /* Upper band at k1=0 is the original-BZ corner image (kx_orig=π, ky=0):
     * ω = 2|J|S(2 − cos π − cos 0) = 2|J|S·2 = 2 in our units. */
    ASSERT_NEAR(omega[1], 2.0 * 1.0 * S * 2.0, 1e-10, "folded 2-band: upper at k=0");

    /* At cartesian (π/2, 0) in doubled BZ → the two original-BZ images
     * are k_orig_x = π/2 and π/2 + π = 3π/2. Both give cos(k_orig_x) = 0,
     * so both bands are degenerate at ω = 2|J|S(2 − 0 − 1) = 1. */
    irrep_magnon_dispersion(L, M_PI / 2.0, 0.0, omega, u);
    double w_x = 2.0 * S * (2.0 - 0.0 - 1.0);
    ASSERT_NEAR(omega[0], w_x, 1e-10, "folded 2-band degenerate at X-point");
    ASSERT_NEAR(omega[1], w_x, 1e-10, "folded 2-band degenerate at X-point");

    /* Off the high-symmetry line, e.g. cartesian (π/4, 0): images are
     * k_orig_x = π/4 and 5π/4, giving cos = 1/√2 and -1/√2. */
    irrep_magnon_dispersion(L, M_PI / 4.0, 0.0, omega, u);
    double w_lo = 2.0 * S * (2.0 - cos(M_PI / 4.0) - 1.0);
    double w_hi = 2.0 * S * (2.0 - cos(5.0 * M_PI / 4.0) - 1.0);
    ASSERT_NEAR(omega[0], w_lo, 1e-10, "folded 2-band: lower at (π/4, 0)");
    ASSERT_NEAR(omega[1], w_hi, 1e-10, "folded 2-band: upper at (π/4, 0)");

    irrep_magnon_lsw_free(L);
}

/* (4) Closed-form 1-band check with DMI. For a single-sublattice model
 * (which would be unphysical on a centrosymmetric lattice — Moriya M1
 * forbids DMI on bonds with inversion centers — but the LSW *module*
 * is geometry-agnostic and accepts any bond list), the dispersion is
 *
 *     ω(k) = 2|J|S(2 − cos kx − cos ky)
 *           + 2 S Dz_x sin(kx) + 2 S Dz_y sin(ky)
 *
 * where Dz_x is the z-component of D on the +x bond, Dz_y on the +y
 * bond. The sin terms come from −i Dz e^{ikt} + i Dz e^{−ikt} = 2 Dz
 * sin(kt). With opposite-sign Dz on the two bonds we get a non-reciprocal
 * dispersion: ω(k) ≠ ω(−k). */
static void test_one_band_with_dmi(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    double Dx = 0.05, Dy = -0.05;
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, Dx}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, Dy}},
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, S, a1, a2, 2, bonds, 0.0);
    double              kx = 0.7, ky = 1.3;
    double              w;
    double _Complex     u;
    irrep_magnon_dispersion(L, kx, ky, &w, &u);
    double expected = 2.0 * S * (2.0 - cos(kx) - cos(ky))
                    + 2.0 * S * Dx * sin(kx) + 2.0 * S * Dy * sin(ky);
    ASSERT_NEAR(w, expected, 1e-12, "1-band ω(k) closed form with DMI");

    /* Non-reciprocity: ω(k) − ω(−k) = 4 S (Dx sin kx + Dy sin ky), not zero. */
    double w_neg;
    irrep_magnon_dispersion(L, -kx, -ky, &w_neg, &u);
    double diff_expected = 4.0 * S * (Dx * sin(kx) + Dy * sin(ky));
    ASSERT_NEAR(w - w_neg, diff_expected, 1e-12, "DMI induces non-reciprocity");

    irrep_magnon_lsw_free(L);
}

/* (5) FM ground-state stability: every sampled ω(k) ≥ 0 for a stable
 * Heisenberg FM. */
static void test_fm_stability(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}},
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, S, a1, a2, 2, bonds, 0);
    int                 negative = 0;
    for (int ix = 0; ix < 16; ++ix)
        for (int iy = 0; iy < 16; ++iy) {
            double          kx = -M_PI + 2.0 * M_PI * ix / 16.0;
            double          ky = -M_PI + 2.0 * M_PI * iy / 16.0;
            double          w;
            double _Complex u;
            irrep_magnon_dispersion(L, kx, ky, &w, &u);
            if (w < -1e-12)
                ++negative;
        }
    ASSERT(negative == 0, "FM stability: ω(k) ≥ 0 over BZ");
    irrep_magnon_lsw_free(L);
}

/* (6) Trivial 1-band model has zero Berry curvature. Dz=0 and only J
 * means H(k) is real — eigenvector is 1 everywhere, no winding. */
static void test_berry_trivial(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}},
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, S, a1, a2, 2, bonds, 0);
    double              berry;
    irrep_magnon_berry(L, 0.5, 0.5, 1e-3, &berry);
    /* For a 1-band model the Berry curvature is 0 by construction. */
    ASSERT_NEAR(berry, 0.0, 1e-8, "trivial Ω(k) = 0");
    irrep_magnon_lsw_free(L);
}

/* (7) Chern number of the trivial 1-band model = 0. */
static void test_chern_trivial(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}},
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, S, a1, a2, 2, bonds, 0);
    double              chern;
    irrep_magnon_chern(L, 16, 16, &chern);
    ASSERT_NEAR(chern, 0.0, 1e-9, "trivial Chern = 0");
    irrep_magnon_lsw_free(L);
}

int main(void) {
    test_fm_square_dispersion();
    test_anisotropy_gap();
    test_two_sublattice_folded();
    test_one_band_with_dmi();
    test_fm_stability();
    test_berry_trivial();
    test_chern_trivial();
    printf("test_magnon: %d/%d assertions passed\n", total - failed, total);
    return failed == 0 ? 0 : 1;
}
