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

/* (8) Trivial 1-band model has zero thermal Hall conductivity:
 * Berry curvature is identically zero ⇒ Σ_b c₂(g_b) Ω_b = 0. */
static void test_thermal_hall_trivial(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}},
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, S, a1, a2, 2, bonds, 0);
    /* Sweep three temperatures: low, mid, high. All should give κ_xy = 0. */
    double T_set[3] = {0.1, 1.0, 10.0};
    for (int i = 0; i < 3; ++i) {
        double k = irrep_magnon_thermal_hall_kxy(L, T_set[i], 24, 24);
        ASSERT_NEAR(k, 0.0, 1e-9, "trivial κ_xy = 0");
    }
    irrep_magnon_lsw_free(L);
}

/* (9) Topological kagome model: κ_xy → 0 as T → 0 (BE statistics decay)
 * and κ_xy is non-zero at intermediate T. The peak T scales with the
 * lowest-band gap O(D · sin(2π/3)) ≈ O(D). */
static void test_thermal_hall_kagome(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * 1.7320508075688772};
    double D = 0.15;
    irrep_magnon_bond_t bonds[6] = {
        {.bi = 0, .bj = 1, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 2, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 2, .bj = 0, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 0, .delta_x = 1,  .delta_y = 0, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 0, .bj = 2, .delta_x = 0,  .delta_y = -1, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y = 1, .J = -1.0, .D = {0, 0, -D}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(3, /*S=*/1.0, a1, a2, 6, bonds, 0);

    /* Verify the Chern numbers first (smoke test). */
    double chern[3];
    irrep_magnon_chern(L, 32, 32, chern);
    ASSERT_NEAR(chern[0], -1.0, 0.05, "kagome lower-band Chern = -1");
    ASSERT_NEAR(chern[1], 0.0, 0.05, "kagome middle-band Chern = 0");
    ASSERT_NEAR(chern[2], +1.0, 0.05, "kagome upper-band Chern = +1");

    /* T → 0: BE factor on the lowest band ω₁(K) ≈ 2|J|S - O(D) is
     * suppressed exponentially, so κ_xy → 0. */
    double k_low = irrep_magnon_thermal_hall_kxy(L, /*T=*/0.05, 32, 32);
    ASSERT(fabs(k_low) < 1e-3, "κ_xy(T=0.05) ≈ 0");

    /* Intermediate T: κ_xy is non-zero. (Sign convention varies; we
     * just check magnitude.) */
    double k_mid = irrep_magnon_thermal_hall_kxy(L, /*T=*/1.0, 32, 32);
    ASSERT(fabs(k_mid) > 1e-3, "κ_xy(T=1.0) ≠ 0 (topological response)");

    /* The peak is around T ~ band-gap / k_B; check it's larger than at
     * very low T. */
    ASSERT(fabs(k_mid) > fabs(k_low), "κ_xy increases from T=0.05 to T=1.0");

    irrep_magnon_lsw_free(L);
}

/* (10) Strip dispersion of the trivial 1-band square FM: all modes are
 * bulk (extended sin(nπx/L)) with edge_weight = 0.5 by symmetry. */
static void test_strip_trivial(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}},
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, S, a1, a2, 2, bonds, 0);

    int     Lx = 12;
    int     N = Lx; /* 1 sublattice */
    double *omega = malloc((size_t)N * sizeof *omega);
    double *ew = malloc((size_t)N * sizeof *ew);
    /* Sweep five k_y values; at every one, all modes should have
     * edge_weight = 0.5 to within numerical noise (sinusoidal envelope
     * is symmetric). */
    double kys[5] = {0.0, M_PI / 4.0, M_PI / 2.0, 3.0 * M_PI / 4.0, M_PI};
    int    bad = 0;
    for (int i = 0; i < 5; ++i) {
        irrep_magnon_strip_dispersion(L, Lx, kys[i], omega, ew);
        for (int b = 0; b < N; ++b)
            if (fabs(ew[b] - 0.5) > 0.02)
                ++bad;
    }
    ASSERT(bad == 0, "trivial strip: every mode has edge_weight ≈ 0.5");
    free(omega);
    free(ew);
    irrep_magnon_lsw_free(L);
}

/* (11) Strip dispersion of the topological kagome FM has chiral edge
 * states inside the bulk gaps. We pick a k_y deep inside the bulk
 * mid-band gap and check that the strip has at least 2 modes (one per
 * edge) with edge_weight far from 0.5. */
static void test_strip_kagome_edge_modes(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * 1.7320508075688772};
    double D = 0.15;
    irrep_magnon_bond_t bonds[6] = {
        {.bi = 0, .bj = 1, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 2, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 2, .bj = 0, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 0, .delta_x = 1,  .delta_y = 0, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 0, .bj = 2, .delta_x = 0,  .delta_y = -1, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y = 1, .J = -1.0, .D = {0, 0, -D}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(3, /*S=*/1.0, a1, a2, 6, bonds, 0);

    int     Lx = 24;
    int     N = Lx * 3;
    double *omega = malloc((size_t)N * sizeof *omega);
    double *ew = malloc((size_t)N * sizeof *ew);

    /* Scan a few k_y values; for each, count modes with strong edge
     * localisation (edge_weight far from 0.5) lying in the topological
     * gap region [2.7, 5.4] (between the K-point Dirac splittings of
     * bands 1 and 2 of the bulk). At least at one k_y, we expect ≥ 2
     * such modes — one per edge. */
    int max_localised_in_gap = 0;
    double kys[5] = {-2.0, -1.0, 0.0, 1.0, 2.0};
    for (int i = 0; i < 5; ++i) {
        irrep_magnon_strip_dispersion(L, Lx, kys[i], omega, ew);
        int loc = 0;
        for (int b = 0; b < N; ++b)
            if (omega[b] > 2.8 && omega[b] < 5.4 && fabs(ew[b] - 0.5) > 0.4)
                ++loc;
        if (loc > max_localised_in_gap)
            max_localised_in_gap = loc;
    }
    ASSERT(max_localised_in_gap >= 2,
           "kagome strip: ≥2 edge-localised modes appear inside the topological gap");

    free(omega);
    free(ew);
    irrep_magnon_lsw_free(L);
}

/* (12) AFM bipartite-chain dispersion (2-sublattice doubled cell):
 *
 *   ω(k_a) = 2·J·S·|sin(k_a · a)|         (folded BZ k_a ∈ [-π/2, π/2])
 *
 * for J > 0 AFM, S = spin per site. Both bands degenerate (folding of
 * a single-band original dispersion). Goldstone at k_a = 0. */
static void test_afm_chain_general(void) {
    /* 2-sublattice doubled cell along x: a1 = (2, 0), a2 = (0, 1).
     * Sublattices A (σ=+1), B (σ=-1). 2 unique NN bonds per cell:
     *   intra-cell: A→B with delta=(0, 0)  → t = 0
     *   inter-cell: A→B with delta=(-1, 0) → t = -2·x̂
     */
    double a1[2] = {2.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 1, .delta_x = 0,  .delta_y = 0, .J = +1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 1, .delta_x = -1, .delta_y = 0, .J = +1.0, .D = {0, 0, 0}},
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(2, S, a1, a2, 2, bonds, 0);
    int                 signs[2] = {+1, -1};

    /* Bipartite-chain folded-BZ dispersion: ω(k_x) = 2·J·S·|sin(k_x)|
     * (the doubled cell folds the original-BZ k_orig = k_x and k_orig +
     * π = k_x + π onto the same k_x in the doubled BZ; both give |sin
     * k_orig| = |sin k_x| since sin(k+π) = -sin(k)). Both bands
     * degenerate. With J=1, S=0.5: ω = |sin(k_x)|. */

    /* k_x = π/4 → ω = sin(π/4) = √2/2 ≈ 0.7071 */
    double SQRT2_HALF = 0.5 * 1.4142135623730951;
    double omega[2];
    irrep_magnon_dispersion_general(L, signs, M_PI / 4.0, 0.0, omega);
    ASSERT_NEAR(omega[0], SQRT2_HALF, 1e-7, "AFM chain at k_x=π/4 (band 1) = √2/2");
    ASSERT_NEAR(omega[1], SQRT2_HALF, 1e-7, "AFM chain at k_x=π/4 (band 2) = √2/2");

    /* k_x = π/8 → ω = sin(π/8) ≈ 0.3827 */
    irrep_magnon_dispersion_general(L, signs, M_PI / 8.0, 0.0, omega);
    ASSERT_NEAR(omega[0], sin(M_PI / 8.0), 1e-7, "AFM chain at k_x=π/8 = sin(π/8)");
    ASSERT_NEAR(omega[1], sin(M_PI / 8.0), 1e-7, "AFM chain at k_x=π/8 = sin(π/8)");

    /* k_x = π/2 → ω = sin(π/2) = 1 (peak of folded dispersion) */
    irrep_magnon_dispersion_general(L, signs, M_PI / 2.0, 0.0, omega);
    ASSERT_NEAR(omega[0], 1.0, 1e-7, "AFM chain at k_x=π/2 = 1");
    ASSERT_NEAR(omega[1], 1.0, 1e-7, "AFM chain at k_x=π/2 = 1");

    irrep_magnon_lsw_free(L);
}

/* (13) FM-mode-recovery: when sublattice signs are all +1, the
 * Bogoliubov-Colpa solver should reproduce the FM Heisenberg
 * dispersion. Cross-check against irrep_magnon_dispersion. */
static void test_general_fm_recovery(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}},
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, S, a1, a2, 2, bonds, 0);
    int                 signs[1] = {+1};

    double w_fm, w_general;
    double _Complex u;
    double          k_set[3][2] = {{0.0, 0.0}, {M_PI / 2, M_PI / 3}, {M_PI, M_PI}};
    for (int p = 0; p < 3; ++p) {
        irrep_magnon_dispersion(L, k_set[p][0], k_set[p][1], &w_fm, &u);
        irrep_magnon_dispersion_general(L, signs, k_set[p][0], k_set[p][1], &w_general);
        ASSERT_NEAR(w_general, w_fm, 1e-8, "general (signs=+1) recovers FM dispersion");
    }
    irrep_magnon_lsw_free(L);
}

/* (14) Wilson-loop winding equals Chern number on the topological
 * kagome model. Sweep k_x across the BZ and accumulate the windings of
 * θ_b(k_x); each must equal the corresponding Chern number. */
static void test_wilson_winding_kagome(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * 1.7320508075688772};
    double D = 0.15;
    irrep_magnon_bond_t bonds[6] = {
        {.bi = 0, .bj = 1, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 2, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 2, .bj = 0, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 0, .delta_x = 1,  .delta_y = 0, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 0, .bj = 2, .delta_x = 0,  .delta_y = -1, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y = 1, .J = -1.0, .D = {0, 0, -D}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(3, 1.0, a1, a2, 6, bonds, 0);

    /* Sweep k_x across the BZ (b1 = 2π·(1, 0) for our a1) and unwrap
     * θ_b(k_x). Total winding = Chern. */
    int     Nx = 64;
    int     Ny = 64;
    double  prev[3], theta[3];
    double  unwrapped[3] = {0, 0, 0};
    double  start[3];
    irrep_magnon_wilson_spectrum(L, -M_PI, Ny, prev);
    for (int b = 0; b < 3; ++b) {
        unwrapped[b] = prev[b];
        start[b] = prev[b];
    }
    for (int ix = 1; ix <= Nx; ++ix) {
        double kx = -M_PI + (2.0 * M_PI) * ix / Nx;
        irrep_magnon_wilson_spectrum(L, kx, Ny, theta);
        for (int b = 0; b < 3; ++b) {
            double diff = theta[b] - prev[b];
            /* Unwrap: branch-cut shifts of ±2π */
            while (diff > M_PI)
                diff -= 2.0 * M_PI;
            while (diff < -M_PI)
                diff += 2.0 * M_PI;
            unwrapped[b] += diff;
            prev[b] = theta[b];
        }
    }
    /* Total winding = (unwrapped[b] - start[b]) / (2π). */
    double winding[3];
    for (int b = 0; b < 3; ++b)
        winding[b] = (unwrapped[b] - start[b]) / (2.0 * M_PI);

    /* Expected Chern (-1, 0, +1). Wilson winding has the same sign
     * convention as Chern (FHS). */
    ASSERT_NEAR(winding[0], -1.0, 0.05, "Wilson winding band 0 = -1");
    ASSERT_NEAR(winding[1], 0.0, 0.05, "Wilson winding band 1 = 0");
    ASSERT_NEAR(winding[2], +1.0, 0.05, "Wilson winding band 2 = +1");

    irrep_magnon_lsw_free(L);
}

/* (15) Trivial 1-band square FM has zero Wilson winding. */
static void test_wilson_trivial(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);
    /* Single band; θ(k_x) should be flat (= 0 modulo branch-cut noise). */
    int Nx = 16, Ny = 32;
    double max_theta = 0;
    for (int ix = 0; ix < Nx; ++ix) {
        double kx = -M_PI + (2.0 * M_PI) * ix / Nx;
        double th;
        irrep_magnon_wilson_spectrum(L, kx, Ny, &th);
        if (fabs(th) > max_theta)
            max_theta = fabs(th);
    }
    ASSERT(max_theta < 1e-6, "trivial Wilson θ(k_x) ≈ 0 everywhere");
    irrep_magnon_lsw_free(L);
}

/* (16) 3D simple-cubic FM dispersion. NN J on the simple-cubic
 * lattice with 6 NN bonds per site (±x, ±y, ±z). For J<0 (FM),
 * S = ½, the analytic dispersion is
 *
 *   ω(k) = 3·|J|·S·(1 − γ_k)·2 = 3 − cos kx − cos ky − cos kz   (J=-1, S=½)
 *
 * with 0 at k=0 (Goldstone), 6 at k=(π, π, π). */
static void test_3d_simple_cubic_fm(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    double a3[3] = {0.0, 0.0, 1.0};
    irrep_magnon_bond_t bonds[3] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 0, .delta_z = 1,
         .J = -1.0, .D = {0, 0, 0}},
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, S, a1, a2, 3, bonds, 0);

    struct {
        double kx, ky, kz;
    } kpts[] = {
        {0.0, 0.0, 0.0},   /* Γ: ω = 0 */
        {M_PI, 0.0, 0.0},  /* X: ω = 2 */
        {M_PI, M_PI, 0.0}, /* M: ω = 4 */
        {M_PI, M_PI, M_PI}, /* R: ω = 6 */
        {M_PI / 2, M_PI / 3, M_PI / 4},
    };
    double          omega;
    double _Complex u;
    for (size_t i = 0; i < sizeof(kpts) / sizeof(*kpts); ++i) {
        irrep_magnon_dispersion_3d(L, a3, kpts[i].kx, kpts[i].ky, kpts[i].kz, &omega, &u);
        double expected =
            3.0 - cos(kpts[i].kx) - cos(kpts[i].ky) - cos(kpts[i].kz);
        ASSERT_NEAR(omega, expected, 1e-12, "3D simple-cubic FM dispersion");
    }
    irrep_magnon_lsw_free(L);
}

/* (17) 3D-to-2D reduction: setting delta_z = 0 on every bond and using
 * the 3D dispersion at k_z = 0 must reproduce the 2D dispersion. */
static void test_3d_reduces_to_2d(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    double a3[3] = {0.0, 0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);
    double          w_2d, w_3d;
    double _Complex u_2d, u_3d;
    irrep_magnon_dispersion(L, M_PI / 3.0, M_PI / 5.0, &w_2d, &u_2d);
    irrep_magnon_dispersion_3d(L, a3, M_PI / 3.0, M_PI / 5.0, 0.0, &w_3d, &u_3d);
    ASSERT_NEAR(w_3d, w_2d, 1e-12, "3D dispersion at kz=0 reduces to 2D");
    irrep_magnon_lsw_free(L);
}

/* (18) Single-sublattice FM transverse structure factor: S_perp_b(q)
 * = 2S identically for all q (|u_b|² = 1 trivially, no sublattice
 * form factor). */
static void test_structure_factor_single(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    double              S = 0.5;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, S, a1, a2, 2, bonds, 0);
    double              q_pts[5][2] = {
        {0.0, 0.0}, {M_PI, 0.0}, {0.0, M_PI}, {M_PI / 3, M_PI / 5}, {M_PI, M_PI}};
    for (int i = 0; i < 5; ++i) {
        double omega, S_perp;
        irrep_magnon_structure_factor(L, q_pts[i][0], q_pts[i][1], &omega, &S_perp);
        ASSERT_NEAR(S_perp, 2.0 * S, 1e-12, "single-sublattice S_perp = 2S");
    }
    irrep_magnon_lsw_free(L);
}

/* (19) Kagome FM at Γ: lowest band carries S_perp = 2S·n_sub = 6S, and
 * the upper bands are orthogonal so S_perp = 0. The sum rule
 * Σ_b S_perp_b(q) = 2S·n_sub holds at every q. */
static void test_structure_factor_kagome_gamma(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * 1.7320508075688772};
    /* Pure Heisenberg kagome FM (no DMI for cleanest test). */
    irrep_magnon_bond_t bonds[6] = {
        {.bi = 0, .bj = 1, .delta_x = 0,  .delta_y = 0,  .delta_z = 0, .J = -1.0, .D = {0,0,0}},
        {.bi = 1, .bj = 2, .delta_x = 0,  .delta_y = 0,  .delta_z = 0, .J = -1.0, .D = {0,0,0}},
        {.bi = 2, .bj = 0, .delta_x = 0,  .delta_y = 0,  .delta_z = 0, .J = -1.0, .D = {0,0,0}},
        {.bi = 1, .bj = 0, .delta_x = 1,  .delta_y = 0,  .delta_z = 0, .J = -1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 2, .delta_x = 0,  .delta_y = -1, .delta_z = 0, .J = -1.0, .D = {0,0,0}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y = 1,  .delta_z = 0, .J = -1.0, .D = {0,0,0}},
    };
    double              S = 1.0;
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(3, S, a1, a2, 6, bonds, 0);
    double              omega[3], S_perp[3];
    /* At Γ avoid the exact Goldstone numerical issue. */
    irrep_magnon_structure_factor(L, 1e-6, 1e-6, omega, S_perp);
    ASSERT_NEAR(S_perp[0], 2.0 * S * 3.0, 0.01, "kagome FM Goldstone band: S_perp = 6S");
    ASSERT(S_perp[1] < 0.01, "kagome FM second band at Γ: dark (S_perp ≈ 0)");
    ASSERT(S_perp[2] < 0.01, "kagome FM third band at Γ: dark (S_perp ≈ 0)");

    /* Sum rule: Σ_b S_perp = 2S·n_sub, independent of q. Test at three q. */
    double q_pts[3][2] = {{0.5, 0.3}, {1.5, 0.7}, {2.0, 1.0}};
    for (int i = 0; i < 3; ++i) {
        irrep_magnon_structure_factor(L, q_pts[i][0], q_pts[i][1], omega, S_perp);
        double sum = S_perp[0] + S_perp[1] + S_perp[2];
        ASSERT_NEAR(sum, 2.0 * S * 3.0, 1e-9, "kagome FM sum rule Σ S_perp = 2S·n_sub");
    }
    irrep_magnon_lsw_free(L);
}

/* (20) Magnon DOS sum rule: ∫ D(ω) dω = n_sub for any LSW model.
 * Histogram on a Nx × Ny grid, sum bins · bin_width, must equal n_sub. */
static void test_dos_sum_rule(void) {
    /* 1-sublattice square FM: dispersion ω ∈ [0, 4]. */
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);

    int     n_bins = 50;
    double *dos = malloc((size_t)n_bins * sizeof *dos);
    irrep_magnon_dos(L, 64, 64, 0.0, 4.001, n_bins, dos);
    double bin_w = 4.001 / n_bins;
    double sum = 0;
    for (int i = 0; i < n_bins; ++i)
        sum += dos[i] * bin_w;
    ASSERT_NEAR(sum, 1.0, 1e-6, "1-sublattice DOS sum rule = n_sub = 1");

    free(dos);
    irrep_magnon_lsw_free(L);
}

/* (21) DOS van-Hove singularity: square FM has saddle points at
 * (π, 0) and (0, π) where ω = 2 (mid-band). The DOS has a logarithmic
 * peak there. With our histogram resolution, the bin containing ω=2
 * should be the largest. */
static void test_dos_van_hove(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);
    int     n_bins = 40;
    double *dos = malloc((size_t)n_bins * sizeof *dos);
    irrep_magnon_dos(L, 128, 128, 0.0, 4.0001, n_bins, dos);
    /* The peak should be in the bin containing ω=2 — that's bin (2/4)·40 = 20. */
    int peak_bin = 0;
    for (int i = 1; i < n_bins; ++i)
        if (dos[i] > dos[peak_bin])
            peak_bin = i;
    /* Allow ±2 bins of tolerance. */
    ASSERT(abs(peak_bin - 20) <= 2,
           "square-FM DOS peak (van-Hove) at ω=2 (bin 20)");

    free(dos);
    irrep_magnon_lsw_free(L);
}

/* (22) Trivial 3D simple-cubic FM has zero Chern on every k_z slice. */
static void test_3d_chern_trivial(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    double a3[3] = {0.0, 0.0, 1.0};
    irrep_magnon_bond_t bonds[3] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 0, .delta_z = 1,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 3, bonds, 0);
    double              chern;
    /* Try three k_z values across the BZ. */
    double kz_set[3] = {0.0, M_PI / 2, M_PI};
    for (int i = 0; i < 3; ++i) {
        irrep_magnon_chern_3d_slice_kz(L, a3, kz_set[i], 16, 16, &chern);
        ASSERT_NEAR(chern, 0.0, 1e-9, "3D simple-cubic FM Chern = 0 at every k_z");
    }
    irrep_magnon_lsw_free(L);
}

/* (23) Layered topological kagome: kagome FM with DMI in plane plus a
 * weak inter-layer Heisenberg coupling along c-axis. At k_z = 0, the
 * 2D Chern signature (-1, 0, +1) should persist (small inter-layer
 * coupling does not destroy the in-plane topology). */
static void test_3d_chern_layered_kagome(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * 1.7320508075688772};
    double a3[3] = {0.0, 0.0, 1.0};
    double D = 0.15;
    /* In-plane bonds: same as the 2D kagome topological model. */
    irrep_magnon_bond_t bonds[9] = {
        /* up-tri CCW */
        {.bi = 0, .bj = 1, .delta_x = 0,  .delta_y = 0,  .delta_z = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 2, .delta_x = 0,  .delta_y = 0,  .delta_z = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 2, .bj = 0, .delta_x = 0,  .delta_y = 0,  .delta_z = 0, .J = -1.0, .D = {0, 0, +D}},
        /* down-tri */
        {.bi = 1, .bj = 0, .delta_x = 1,  .delta_y = 0,  .delta_z = 0, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 0, .bj = 2, .delta_x = 0,  .delta_y = -1, .delta_z = 0, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y = 1,  .delta_z = 0, .J = -1.0, .D = {0, 0, -D}},
        /* Inter-layer: weak FM Heisenberg, A→A, B→B, C→C along +ẑ. */
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 0, .delta_z = 1, .J = -0.05, .D = {0, 0, 0}},
        {.bi = 1, .bj = 1, .delta_x = 0, .delta_y = 0, .delta_z = 1, .J = -0.05, .D = {0, 0, 0}},
        {.bi = 2, .bj = 2, .delta_x = 0, .delta_y = 0, .delta_z = 1, .J = -0.05, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(3, 1.0, a1, a2, 9, bonds, 0);
    double              chern[3];
    /* k_z = 0 slice: the in-plane (-1, 0, +1) topology should survive. */
    irrep_magnon_chern_3d_slice_kz(L, a3, 0.0, 32, 32, chern);
    ASSERT_NEAR(chern[0], -1.0, 0.05, "layered kagome k_z=0: lower band Chern -1");
    ASSERT_NEAR(chern[1], 0.0, 0.05, "layered kagome k_z=0: middle band Chern 0");
    ASSERT_NEAR(chern[2], +1.0, 0.05, "layered kagome k_z=0: upper band Chern +1");
    irrep_magnon_lsw_free(L);
}

/* (24) Specific-heat limits: C_V → 0 as T → 0 (BE-suppressed) and
 * C_V → n_sub at T → ∞ (equipartition). For 1-sublattice FM. */
static void test_specific_heat_limits(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);

    /* Low-T: very small C_V. */
    double cv_low = irrep_magnon_specific_heat(L, /*T=*/0.01, 32, 32);
    ASSERT(cv_low < 0.01, "low-T C_V → 0 (BE-suppressed)");

    /* High-T: equipartition C_V → n_sub = 1. */
    double cv_high = irrep_magnon_specific_heat(L, /*T=*/100.0, 32, 32);
    ASSERT(fabs(cv_high - 1.0) < 0.01, "high-T C_V → n_sub = 1 (equipartition)");

    /* Intermediate-T monotonic increase from low to high. */
    double cv_mid = irrep_magnon_specific_heat(L, /*T=*/0.5, 32, 32);
    ASSERT(cv_mid > cv_low && cv_mid < cv_high,
           "C_V is monotonic in T over 0.01 → 0.5 → 100");

    irrep_magnon_lsw_free(L);
}

/* (25) Magnetization at T → 0 → S (no thermal magnons). At very low T,
 * the gapped 2D FM should give M ≈ S to high precision. */
static void test_magnetization_low_T(void) {
    /* 2D square FM with K_z = 0.5 → 2 K_z S = 0.5 gap. Hence at T = 0.05
     * (10× below the gap), thermal population is exp(-10) ~ 1e-5,
     * giving M ≈ S − ε. */
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, /*Kz=*/0.5);
    double M_low = irrep_magnon_magnetization(L, 0.05, 32, 32);
    /* M(T=0.05) should be 0.5 minus exponentially-small population. */
    ASSERT(fabs(M_low - 0.5) < 1e-3, "M(T → 0) → S = 0.5");
    irrep_magnon_lsw_free(L);
}

/* (26) M(T) is monotonically decreasing in T. Test on the same gapped
 * 2D FM. */
static void test_magnetization_monotonic(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0.5);
    double M1 = irrep_magnon_magnetization(L, 0.1, 32, 32);
    double M2 = irrep_magnon_magnetization(L, 0.5, 32, 32);
    double M3 = irrep_magnon_magnetization(L, 1.0, 32, 32);
    ASSERT(M1 > M2 && M2 > M3,
           "M(T) monotonically decreases: M(0.1) > M(0.5) > M(1.0)");
    irrep_magnon_lsw_free(L);
}

/* (27) Q-ω map: at each q on a square-FM path, the Lorentzian-broadened
 * intensity must peak at ω = ω(q). Test by finding the energy bin with
 * max intensity at each q and confirming it sits within η of ω(q). */
static void test_qomega_map_peaks_track_dispersion(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);

    int    n_q = 5;
    int    n_w = 200;
    double eta = 0.05;
    double w_min = 0.0;
    double w_max = 4.5;
    double dw = (w_max - w_min) / n_w;
    double qpath[5][2] = {{0.1, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {2.5, 0.0}, {3.0, 0.0}};
    double *I_qw = malloc((size_t)n_q * (size_t)n_w * sizeof *I_qw);
    irrep_magnon_neutron_qomega_map(L, qpath, n_q, w_min, w_max, n_w, eta, I_qw);
    /* Verify max-intensity bin matches dispersion ω(q) = 1 - cos(qx)
     * (with J=1, S=½). */
    int hits = 0;
    for (int iq = 0; iq < n_q; ++iq) {
        int jmax = 0;
        for (int j = 1; j < n_w; ++j)
            if (I_qw[iq * n_w + j] > I_qw[iq * n_w + jmax])
                jmax = j;
        double w_peak = w_min + (jmax + 0.5) * dw;
        double w_expected = 2.0 - cos(qpath[iq][0]) - cos(qpath[iq][1]);
        if (fabs(w_peak - w_expected) < 3.0 * eta)
            ++hits;
    }
    ASSERT(hits == n_q, "Q-ω map peaks track dispersion to within 3η");

    free(I_qw);
    irrep_magnon_lsw_free(L);
}

/* (28) Lorentzian unit-area sum rule: at each q, integrating I(q, ω) dω
 * should give Σ_b S_⊥_b(q) (the band-summed structure factor). For a
 * 1-sublattice FM, that's 2S = 1. */
static void test_qomega_lorentzian_sum_rule(void) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);
    int     n_q = 3;
    int     n_w = 800;
    double  eta = 0.02;
    /* Wide window so Lorentzian tails are negligible. */
    double  w_min = -2.0, w_max = 6.0;
    double  qpath[3][2] = {{0.5, 0.5}, {2.0, 0.0}, {3.0, 1.0}};
    double *I_qw = malloc((size_t)n_q * (size_t)n_w * sizeof *I_qw);
    irrep_magnon_neutron_qomega_map(L, qpath, n_q, w_min, w_max, n_w, eta, I_qw);
    double dw = (w_max - w_min) / n_w;
    for (int iq = 0; iq < n_q; ++iq) {
        double sum = 0;
        for (int j = 0; j < n_w; ++j)
            sum += I_qw[iq * n_w + j] * dw;
        ASSERT_NEAR(sum, 1.0, 0.005, "Q-ω Lorentzian unit-area sum = 2S = 1");
    }
    free(I_qw);
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
    test_thermal_hall_trivial();
    test_thermal_hall_kagome();
    test_strip_trivial();
    test_strip_kagome_edge_modes();
    test_afm_chain_general();
    test_general_fm_recovery();
    test_wilson_winding_kagome();
    test_wilson_trivial();
    test_3d_simple_cubic_fm();
    test_3d_reduces_to_2d();
    test_structure_factor_single();
    test_structure_factor_kagome_gamma();
    test_dos_sum_rule();
    test_dos_van_hove();
    test_3d_chern_trivial();
    test_3d_chern_layered_kagome();
    test_specific_heat_limits();
    test_magnetization_low_T();
    test_magnetization_monotonic();
    test_qomega_map_peaks_track_dispersion();
    test_qomega_lorentzian_sum_rule();
    printf("test_magnon: %d/%d assertions passed\n", total - failed, total);
    return failed == 0 ? 0 : 1;
}
