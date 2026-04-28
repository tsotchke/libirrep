/* SPDX-License-Identifier: MIT */
/* Polarised 1-magnon structure factor demo on a square FM. Shows the
 * channel split S^xx, S^yy, S^zz with three different ground-state
 * orientations:
 *   1. Collinear FM along ẑ — S^xx = S^yy, S^zz = 0
 *   2. Collinear FM along x̂ — S^yy = S^zz, S^xx = 0
 *   3. Non-collinear (single sublattice tilted to (1,1,1)/√3) — all
 *      three channels carry weight
 *
 * For polarised-INS experiments the channel split distinguishes
 * spin-flip scattering (transverse) from non-spin-flip
 * (longitudinal). LSW + the new _polarized_structure_factor gives
 * this directly.
 *
 * Build: make examples
 * Run:   ./build/bin/polarized_sf_demo */

#include <irrep/magnon.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Polarised 1-magnon SF demo — square FM under three orientations\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[2] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0,0,0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);

    /* Three orientations:
     *   A: collinear ẑ (default)
     *   B: collinear x̂
     *   C: tilted (1,1,1)/√3
     */
    double n_z[1][3]   = {{0, 0, 1}};
    double n_x[1][3]   = {{1, 0, 0}};
    double inv_root3   = 1.0 / sqrt(3.0);
    double n_xyz[1][3] = {{inv_root3, inv_root3, inv_root3}};

    struct {
        const char  *name;
        double     (*n_hat)[3];
    } cases[] = {
        {"FM ẑ        (collinear z)",     NULL},   /* NULL → assume ẑ */
        {"FM ẑ explicit",                 n_z},
        {"FM x̂        (collinear x)",     n_x},
        {"non-collinear (1,1,1)/√3",      n_xyz},
    };
    int n_cases = 4;

    /* Probe at a generic q = (π/2, π/3) — interior of BZ, away from Γ. */
    double qx = M_PI / 2;
    double qy = M_PI / 3;

    printf("  At q = (π/2, π/3):\n\n");
    printf("  %-32s   %-9s   %-9s   %-9s   %-9s\n", "ground state", "ω_b", "S^xx", "S^yy",
           "S^zz");
    printf("  ─────────────────────────────────┼───────────┼───────────┼───────────┼──────────\n");
    for (int c = 0; c < n_cases; ++c) {
        double omega[1], Sxx[1], Syy[1], Szz[1];
        irrep_magnon_polarized_structure_factor(L, cases[c].n_hat, qx, qy, omega, Sxx, Syy,
                                                  Szz);
        printf("  %-32s   %-9.4f   %-9.4f   %-9.4f   %-9.4f\n", cases[c].name, omega[0],
               Sxx[0], Syy[0], Szz[0]);
    }
    printf("\n");

    printf("  PHYSICS:\n");
    printf("  - FM along ẑ: spin-flip scattering carries weight in S^xx and\n");
    printf("    S^yy by axial symmetry around the ordering axis. S^zz = 0\n");
    printf("    because longitudinal spin fluctuations vanish at the\n");
    printf("    1-magnon level for collinear FMs.\n");
    printf("  - FM along x̂: same physics, but the lab axes that carry\n");
    printf("    1-magnon weight are now (y, z) — perpendicular to x̂.\n");
    printf("  - Tilted (1,1,1)/√3: all three channels are non-zero. The\n");
    printf("    longitudinal (S^zz_lab) channel acquires 1-magnon weight\n");
    printf("    when the ordering axis tilts away from ẑ. The total\n");
    printf("    transverse-to-ordering-axis weight is conserved\n");
    printf("    (frame-invariant); the channel SPLIT depends on the\n");
    printf("    rotation matrix.\n");
    printf("  - Useful for: polarised INS where spin-flip / non-spin-flip\n");
    printf("    discrimination distinguishes transverse from longitudinal\n");
    printf("    fluctuations directly, and for ML pipelines that need\n");
    printf("    polarised observables as training labels.\n\n");

    irrep_magnon_lsw_free(L);
    return 0;
}
