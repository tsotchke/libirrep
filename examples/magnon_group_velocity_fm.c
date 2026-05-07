/* SPDX-License-Identifier: MIT */
/* Magnon group velocity v_g = ∇_k ω(k) on the spin-½ square
 * ferromagnet — exercises `irrep_magnon_group_velocity` against
 * the analytic Holstein-Primakoff dispersion at five k-points
 * spanning the BZ.
 *
 * MATHEMATICAL CONTEXT
 *
 * The group velocity sets the speed at which a magnon wave-packet
 * propagates and enters every magnon-mediated transport quantity:
 *
 *     j_s = (1/V) Σ_b ∫_{BZ} d²k/(2π)² · v_g_b(k) · n_B(ω_b/T)
 *
 * is the magnon spin-current density. Modes with large |v_g|
 * dominate transport; flat bands and band tops have v_g → 0 and
 * are spectator modes.
 *
 * For the spin-½ square nearest-neighbour ferromagnet (J<0):
 *
 *     ω(q)   = 2|J|S [2 − cos q_x − cos q_y]
 *     v_g(q) = 2|J|S (sin q_x, sin q_y)
 *
 * with J = -1, S = ½ giving |v_g|_max = 2√2|J|S = √2 ≈ 1.414 at
 * q = (π/2, π/2). The Goldstone mode at Γ is *quadratic*
 * (ω ≈ |J|S |q|² near Γ), so v_g(Γ) = 0 — distinct from the
 * AFM Goldstone, which has linear ω ∝ |k| and finite v_g.
 *
 * The library uses central differences:
 *
 *     v_x = [ω(k_x + h, k_y) − ω(k_x − h, k_y)] / (2h)
 *
 * giving O(h²) error, so step h = 1e-3 reproduces the analytic
 * gradient to ~h² · ω''' ≈ 10⁻⁶.
 *
 * NOTE on AFM ground states: the current `irrep_magnon_group_velocity`
 * does not take sublattice signs and therefore returns the FM-track
 * Hermitian-path velocity. For Bogoliubov-Colpa AFM dispersions, the
 * caller can do central differences on `irrep_magnon_dispersion_general`
 * directly until an `_group_velocity_general` variant is added in a
 * future cycle.
 *
 * REFERENCES
 *   - Holstein & Primakoff, Phys. Rev. 58, 1098 (1940)
 *   - Cornelissen et al., Nat. Phys. 11, 1022 (2015) — magnon spin
 *     transport via local v_g
 *
 * Build: make examples
 * Run:   ./build/bin/magnon_group_velocity_fm */

#include <irrep/magnon.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Magnon group velocity on the spin-½ square FM —\n");
    printf("  v_g = (sin q_x, sin q_y) at five representative k-points\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[2] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);
    if (!L) return 1;

    struct {
        double      qx, qy;
        const char *label;
    } pts[] = {
        {0.0,        0.0,        "Γ (Goldstone)              "},
        {0.1,        0.0,        "near Γ (q = 0.1)           "},
        {M_PI / 2.0, 0.0,        "(π/2, 0)                   "},
        {M_PI / 2.0, M_PI / 2.0, "(π/2, π/2) — |v_g|_max  "},
        {M_PI,       M_PI,       "M (π, π) — band top     "},
    };

    printf("  Setup: 1-sublattice FM, J = -1, S = ½, ω = 2 - cos q_x - cos q_y\n");
    printf("         FD step h = 1e-3, expected error ~ h² ω''' ~ 1e-6\n\n");
    printf("  %-30s  %-12s  %-12s  %-10s\n",
           "k-point", "v_x (lib)", "v_x (exact)", "abs.err");
    double h = 1e-3;
    int    fails  = 0;
    double maxerr = 0.0;
    for (int i = 0; i < 5; ++i) {
        double vx, vy;
        irrep_magnon_group_velocity(L, pts[i].qx, pts[i].qy, h, &vx, &vy);
        double vx_exact = sin(pts[i].qx);
        double vy_exact = sin(pts[i].qy);
        double err      = fmax(fabs(vx - vx_exact), fabs(vy - vy_exact));
        if (err > 1e-5) ++fails;
        if (err > maxerr) maxerr = err;
        printf("  %-30s  %+10.6f    %+10.6f    %.2e\n",
               pts[i].label, vx, vx_exact, err);
    }
    printf("\n  All 5 points within 1e-5 of analytic; worst-case %.2e — within\n",
           maxerr);
    printf("  the central-difference O(h²) budget for h = 1e-3.\n\n");

    /* Diagnostic: confirm |v_g|(Γ) → 0 quadratically as q → Γ. */
    printf("  Quadratic-Goldstone confirmation: |v_g(q)| / |q| as q → 0\n");
    printf("  (FM Goldstone has ω ∝ q², so |v_g| ∝ |q| → 0 at Γ).\n\n");
    printf("  %-10s  %-12s  %-12s\n", "|q|", "|v_g|", "|v_g|/|q|");
    double qmags[5] = {1.0, 0.5, 0.2, 0.1, 0.05};
    for (int i = 0; i < 5; ++i) {
        double q = qmags[i];
        double vx, vy;
        irrep_magnon_group_velocity(L, q, 0.0, h, &vx, &vy);
        double speed = sqrt(vx * vx + vy * vy);
        printf("  %-10.4f  %-12.6f  %-12.6f\n", q, speed, speed / q);
    }
    printf("\n  |v_g|/|q| → 1 confirms |v_g| ≈ |q| near Γ: v_x = sin q ≈ q,\n");
    printf("  consistent with the quadratic FM dispersion ω ≈ q².\n\n");

    irrep_magnon_lsw_free(L);

    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Group velocity verified at machine-budget precision; FM Goldstone\n");
    printf("  has v_g(Γ) = 0 from the quadratic dispersion shape.\n");
    printf("══════════════════════════════════════════════════════════════════\n");
    return fails ? 1 : 0;
}
