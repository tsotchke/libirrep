/* SPDX-License-Identifier: MIT */
/* Kagome FM topological transition under DM-coefficient sign flip —
 * demonstrating the libirrep public Chern-sweep + boundary-detection
 * workflow.
 *
 * MATHEMATICAL CONTEXT
 *
 * The kagome FM with out-of-plane DMI Dz has band Chern numbers
 *
 *   sign(D) = +1: (C_1, C_2, C_3) = (-1, 0, +1)
 *   sign(D) = -1: (C_1, C_2, C_3) = (+1, 0, -1)
 *
 * The two regimes are separated by D = 0, where the bands touch at the
 * K point of the magnetic Brillouin zone (Dirac cones, gap = 0). This
 * is the simplest topological phase transition on the magnon side: the
 * Dirac mass changes sign through the gap-closing point, and each band
 * exchanges a unit of Berry-curvature flux with the band it touches.
 *
 * THIS EXAMPLE
 *
 * Walks the workflow for finding such transitions automatically:
 *
 *   1. Define a factory callback that constructs the LSW handle from a
 *      single scalar parameter (here D).
 *   2. Call irrep_magnon_chern_sweep with a coarse parameter grid
 *      across the suspected transition region.
 *   3. Call irrep_magnon_chern_detect_boundaries on the resulting
 *      Chern table to flag the indices where the integer signature
 *      changes between adjacent samples.
 *   4. Refine: rerun the sweep on a tighter grid around each flagged
 *      boundary to bracket the transition more precisely.
 *
 * Both the coarse and refined sweeps avoid D = 0 itself (the bands
 * touch there, so the FHS-plaquette method on a finite k-grid returns
 * a non-integer "Chern" value due to the band-degeneracy ambiguity).
 *
 * Build: make examples
 * Run:   ./build/bin/kagome_dmi_topological_transition */

#include <irrep/magnon.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static irrep_magnon_lsw_t *kagome_fm_dmi_factory(double D, void *ud) {
    (void)ud;
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * sqrt(3.0)};
    irrep_magnon_bond_t bonds[6] = {
        {.bi = 0, .bj = 1, .delta_x =  0, .delta_y =  0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 2, .delta_x =  0, .delta_y =  0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 2, .bj = 0, .delta_x =  0, .delta_y =  0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 0, .delta_x =  1, .delta_y =  0, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 0, .bj = 2, .delta_x =  0, .delta_y = -1, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y =  1, .J = -1.0, .D = {0, 0, -D}},
    };
    return irrep_magnon_lsw_new(/*n_sub=*/3, /*S=*/1.0, a1, a2, /*n_bonds=*/6, bonds, /*Kz=*/0);
}

static void print_chern_row(double D, const double *chern) {
    printf("    D = %+7.4f   (% .4f, % .4f, % .4f)   Σ = %+8.1e\n",
           D, chern[0], chern[1], chern[2], chern[0] + chern[1] + chern[2]);
}

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Kagome FM topological transition under D → -D\n");
    printf("  Demonstrating chern_sweep + chern_detect_boundaries workflow\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    /* === Stage 1: coarse sweep across the transition region === */
    printf("  STAGE 1 — Coarse parameter sweep, |D| ∈ {0.10, 0.20, 0.30}\n\n");
    double D_coarse[6] = {-0.30, -0.20, -0.10, +0.10, +0.20, +0.30};
    int    n_coarse    = 6;
    double chern_coarse[18];

    irrep_status_t st = irrep_magnon_chern_sweep(kagome_fm_dmi_factory, NULL, D_coarse, n_coarse,
                                                  /*n_sub=*/3, /*Nx=*/64, /*Ny=*/64, chern_coarse);
    if (st != IRREP_OK) {
        fprintf(stderr, "chern_sweep failed: %s\n", irrep_strerror(st));
        return 1;
    }
    for (int i = 0; i < n_coarse; ++i) print_chern_row(D_coarse[i], chern_coarse + i * 3);

    int boundaries_coarse[6];
    int n_boundaries_coarse = -1;
    st = irrep_magnon_chern_detect_boundaries(chern_coarse, n_coarse, /*n_sub=*/3,
                                                boundaries_coarse, /*max=*/6,
                                                &n_boundaries_coarse);
    if (st != IRREP_OK) {
        fprintf(stderr, "detect_boundaries failed\n");
        return 1;
    }
    printf("\n  Detected %d boundary(ies):\n", n_boundaries_coarse);
    for (int k = 0; k < n_boundaries_coarse; ++k) {
        int    i = boundaries_coarse[k];
        double Dl = D_coarse[i], Dr = D_coarse[i + 1];
        printf("    [%2d]  bracketed by D ∈ [%+.4f, %+.4f]  (gap = %.4f)\n",
               k, Dl, Dr, Dr - Dl);
    }

    /* === Stage 2: refine each detected boundary === */
    printf("\n  STAGE 2 — Refinement: tighter sweep on the bracketed interval\n\n");
    for (int k = 0; k < n_boundaries_coarse; ++k) {
        int    i  = boundaries_coarse[k];
        double Dl = D_coarse[i], Dr = D_coarse[i + 1];
        printf("  Refining boundary [%d] in [%+.4f, %+.4f]:\n", k, Dl, Dr);

        /* Six points spanning the bracket; skip exact zero (bands touch). */
        double D_fine[6];
        for (int j = 0; j < 6; ++j) {
            double t = (j + 0.5) / 6.0;
            D_fine[j] = Dl + t * (Dr - Dl);
            if (fabs(D_fine[j]) < 5e-3) D_fine[j] = (D_fine[j] > 0 ? 1 : -1) * 0.005;
        }

        double chern_fine[18];
        st = irrep_magnon_chern_sweep(kagome_fm_dmi_factory, NULL, D_fine, 6, 3, 96, 96,
                                       chern_fine);
        if (st != IRREP_OK) {
            fprintf(stderr, "  refinement chern_sweep failed: %s\n", irrep_strerror(st));
            continue;
        }
        for (int j = 0; j < 6; ++j) print_chern_row(D_fine[j], chern_fine + j * 3);

        int boundaries_fine[6];
        int n_boundaries_fine = -1;
        st = irrep_magnon_chern_detect_boundaries(chern_fine, 6, 3, boundaries_fine, 6,
                                                    &n_boundaries_fine);
        if (st == IRREP_OK && n_boundaries_fine > 0) {
            int    j  = boundaries_fine[0];
            double DL = D_fine[j], DR = D_fine[j + 1];
            printf("    refined bracket: D ∈ [%+.4f, %+.4f]  (gap = %.4f)\n\n",
                   DL, DR, DR - DL);
        } else {
            printf("    (no refined boundary detected — try a denser sub-grid)\n\n");
        }
    }

    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  PHYSICAL INTERPRETATION\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");
    printf("  The detected transition at D = 0 is the canonical magnon\n");
    printf("  Dirac-cone gap closing: at D > 0 the lower band has Chern\n");
    printf("  number C₁ = -1 (Berry-curvature flux of -2π through the BZ);\n");
    printf("  at D < 0 the same band has C₁ = +1. The gap-closing point is\n");
    printf("  the magnon analogue of the Haldane-model Dirac-cone gap-\n");
    printf("  closing under haldane-mass sign flip.\n\n");
    printf("  Workflow generality: any factory callback that builds an LSW\n");
    printf("  handle from a scalar parameter — vary J, J', K_z, applied\n");
    printf("  field, or any combination — drives the same chern_sweep +\n");
    printf("  detect_boundaries pipeline. The library returns the boundary\n");
    printf("  brackets; the user refines.\n\n");

    return 0;
}
