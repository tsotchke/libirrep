/* SPDX-License-Identifier: MIT */
/* Chiral edge modes on a kagome FM strip — bulk-boundary correspondence
 * for topological magnons.
 *
 * Same model as kagome_topological_magnons.c (FM Heisenberg J = -1
 * with NN DMI Dz alternating by triangle parity, S = 1, Mook-Henk-
 * Mertig (-1, 0, +1) Chern signature). Open boundary conditions
 * along a₁ (Lx unit cells stacked), periodic along a₂. Dispersion
 * ω(k_y) on the strip exhibits gapless chiral edge modes that
 * traverse every bulk gap.
 *
 * The Hatsugai bulk-boundary correspondence theorem (1993) says: at
 * each edge of the strip, the *net* number of chiral edge modes in
 * a given bulk gap equals the sum of bulk Chern numbers of all
 * bands below that gap. For our (-1, 0, +1) signature:
 *
 *   - lower gap (between bands 1 and 2):  C_below = -1
 *     → 1 chiral mode per edge, opposite chirality on the two edges
 *
 *   - upper gap (between bands 2 and 3):  C_below = -1 + 0 = -1
 *     → 1 chiral mode per edge
 *
 * Each chiral mode carries an unidirectional magnon current —
 * left-edge modes flow one way, right-edge the other. The mode is
 * topologically protected: small symmetry-preserving perturbations
 * cannot back-scatter or gap it out, only push the dispersion.
 *
 * This is the bosonic analog of quantum-Hall edge channels and the
 * direct mechanism behind the chiral magnon transport observed at
 * the edge of Cu(1,3-bdc) flakes (Akazawa 2020).
 *
 * Build: `make examples`
 * Run:   `./build/bin/kagome_chiral_edge_modes` */

#include <irrep/magnon.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== libirrep — chiral edge modes on a kagome-FM strip ===\n\n");
    printf("    Setup: 3-sublattice kagome, NN J=-1 (FM), Dz=±0.15 by triangle\n");
    printf("           parity, S=1. Strip: Lx=24 cells (open BC on a₁), PBC on\n");
    printf("           a₂. Bulk Chern signature: (-1, 0, +1).\n\n");

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
    if (!L) {
        fprintf(stderr, "magnon LSW handle allocation failed\n");
        return 1;
    }

    int     Lx = 24;
    int     N = Lx * 3;
    double *omega = malloc((size_t)N * sizeof *omega);
    double *ew = malloc((size_t)N * sizeof *ew);

    /* k_y scan: for each k_y, count and locate edge-localised modes
     * inside the lower topological gap (ω ∈ [2.8, 3.4] — between the
     * upper edge of the lower bulk band and the lower edge of the
     * middle bulk band, both K-point-split by O(D)). */
    printf("    Edge-mode count inside lower gap (ω ∈ [2.8, 3.4]) vs k_y:\n");
    printf("    %-9s  %-7s  %-7s  %-s\n", "k_y", "n_left", "n_right", "ω of edge modes");
    printf("    ──────────────────────────────────────────────────────────\n");
    for (int i = 0; i <= 12; ++i) {
        double ky = -M_PI + (2.0 * M_PI) * i / 12.0;
        irrep_magnon_strip_dispersion(L, Lx, ky, omega, ew);
        int    n_left = 0, n_right = 0;
        double l_ws[8] = {0}, r_ws[8] = {0};
        for (int b = 0; b < N; ++b) {
            if (omega[b] < 2.7 || omega[b] > 3.5)
                continue;
            if (ew[b] > 0.85 && n_left < 8)
                l_ws[n_left++] = omega[b];
            else if (ew[b] < 0.15 && n_right < 8)
                r_ws[n_right++] = omega[b];
        }
        printf("    %+9.3f  %-7d  %-7d  ", ky, n_left, n_right);
        for (int j = 0; j < n_left; ++j)
            printf("L=%.3f ", l_ws[j]);
        for (int j = 0; j < n_right; ++j)
            printf("R=%.3f ", r_ws[j]);
        printf("\n");
    }

    /* Snapshot at k_y where the test_magnon assertion fires —
     * confirm at least one k_y carries 2 edge-localised modes. */
    printf("\n    Snapshot scan: at which k_y does the gap host both an L and R mode?\n");
    int n_kpoints_with_both = 0;
    for (int i = 0; i <= 60; ++i) {
        double ky = -M_PI + (2.0 * M_PI) * i / 60.0;
        irrep_magnon_strip_dispersion(L, Lx, ky, omega, ew);
        int    has_L = 0, has_R = 0;
        for (int b = 0; b < N; ++b) {
            if (omega[b] < 2.7 || omega[b] > 3.5)
                continue;
            if (ew[b] > 0.85)
                has_L = 1;
            if (ew[b] < 0.15)
                has_R = 1;
        }
        if (has_L && has_R)
            ++n_kpoints_with_both;
    }
    printf("    %d / 61 sampled k_y points host BOTH a left- and right-edge mode\n",
           n_kpoints_with_both);
    printf("    inside the lower topological gap. (At other k_y, the edge modes\n");
    printf("    have dispersed into the bulk bands above/below the gap window.)\n");

    printf("\n  ━ Interpretation ━\n");
    printf("    The right-edge mode (R) disperses ω(k_y) ↗ as k_y ↘ — it\n");
    printf("    propagates *one direction* along the right edge (right-moving\n");
    printf("    convention). The left-edge mode (L) has the opposite slope:\n");
    printf("    ω(k_y) ↘ as k_y ↘, so it propagates the other direction along\n");
    printf("    the left edge. The two edges carry counter-propagating chiral\n");
    printf("    magnon currents — the bosonic analog of quantum-Hall edge\n");
    printf("    channels.\n\n");
    printf("    Counted per BZ winding: the edge-mode dispersion crosses any\n");
    printf("    horizontal line ω=ε∈[2.8, 3.4] *exactly once* on each edge.\n");
    printf("    The signed crossing count on a given edge is the Chern number\n");
    printf("    of the band below the gap (here C₁ = -1) — the Hatsugai 1993\n");
    printf("    bulk-boundary correspondence theorem.\n\n");
    printf("    This is the bulk-boundary correspondence theorem (Hatsugai\n");
    printf("    1993) operating on magnons. The same mechanism that gives\n");
    printf("    quantum-Hall systems their robust edge currents gives this\n");
    printf("    kagome FM unidirectional magnon transport at its physical\n");
    printf("    boundary — observed in Cu(1,3-bdc) (Akazawa 2020) and in\n");
    printf("    Lu₂V₂O₇ pyrochlore (Hirschberger 2015).\n\n");
    printf("    Closes another loop: from a libirrep DMI bond list to the\n");
    printf("    *direct geometric signature* of the topological invariant —\n");
    printf("    edge modes that the bulk Chern numbers predict by the\n");
    printf("    Hatsugai correspondence.\n");

    free(omega);
    free(ew);
    irrep_magnon_lsw_free(L);
    return 0;
}
