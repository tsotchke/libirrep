/* SPDX-License-Identifier: MIT */
/* Magnon decay rate Γ_b(k) on the kagome 120° Néel — canonical
 * usage of `irrep_magnon_heisenberg_decay_rate` with explicit
 * `gap_floor` for the IR cutoff.
 *
 * The Born self-energy is logarithmically IR-divergent on
 * Goldstone-having models like the kagome 120° Néel. The function
 * requires the caller to supply `gap_floor` representing the
 * physical model gap (anisotropy, applied field, dissipation
 * width, finite-T self-consistent gap, etc.). At any FIXED
 * gap_floor the BZ integral converges as Nx, Ny → ∞.
 *
 * This example sweeps several gap_floor values to show the
 * regularization-dependence explicitly so the user understands
 * what's being computed.
 *
 * Build: make examples
 * Run:   ./build/bin/kagome_120_neel_decay */

#include <irrep/magnon.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Magnon decay Γ_b(k) on kagome 120° Néel\n");
    printf("  (Heisenberg cubic-vertex Born self-energy with gap_floor IR cutoff)\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    /* Standard kagome bond list, AFM J=1, S=1/2. */
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * sqrt(3.0)};
    irrep_magnon_bond_t bonds[6] = {
        {.bi = 0, .bj = 1, .delta_x = 0,  .delta_y = 0,  .J = +1.0, .D = {0,0,0}},
        {.bi = 1, .bj = 2, .delta_x = 0,  .delta_y = 0,  .J = +1.0, .D = {0,0,0}},
        {.bi = 2, .bj = 0, .delta_x = 0,  .delta_y = 0,  .J = +1.0, .D = {0,0,0}},
        {.bi = 1, .bj = 0, .delta_x = 1,  .delta_y = 0,  .J = +1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 2, .delta_x = 0,  .delta_y = -1, .J = +1.0, .D = {0,0,0}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y = 1,  .J = +1.0, .D = {0,0,0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(3, 0.5, a1, a2, 6, bonds, 0);
    /* 120° Néel ordering vectors in the lattice plane (Chubukov-Senthil). */
    double n_vecs[9] = {
        +1.0,           0.0,           0.0,
        -0.5,           0.5*sqrt(3.0), 0.0,
        -0.5,          -0.5*sqrt(3.0), 0.0,
    };

    /* Probe at three high-symmetry k-points. */
    double k_pts[3][2] = {
        {1.5, 0.7},                   /* zone-interior */
        {M_PI, M_PI/sqrt(3.0)},       /* M-point */
        {4*M_PI/3, 0},                /* K-point */
    };
    const char *labels[3] = {"interior", "M-point", "K-point"};
    int n_k = 3;

    /* Sweep gap_floor at fixed (N, η) to expose IR-divergence behavior. */
    double gap_floors[5] = {0.5, 0.2, 0.1, 0.05, 0.025};
    int n_gaps = 5;

    printf("  Sweep over gap_floor at N=24, η=0.1 (Born linewidth in J units):\n\n");
    printf("  %-9s   %-12s   %-12s   %-12s\n", "gap_floor", "Γ_band 0", "Γ_band 1", "Γ_band 2");
    for (int q = 0; q < n_k; ++q) {
        printf("\n  q = %s :\n", labels[q]);
        printf("  ─────────┼──────────────┼──────────────┼──────────────\n");
        for (int g = 0; g < n_gaps; ++g) {
            double k_one[1][2] = {{k_pts[q][0], k_pts[q][1]}};
            double Gamma[3];
            irrep_magnon_heisenberg_decay_rate(L, n_vecs, k_one, 1, 24, 24, 0.1,
                                                 gap_floors[g], Gamma);
            printf("  %-9.3f   %-12.4f   %-12.4f   %-12.4f\n", gap_floors[g], Gamma[0],
                   Gamma[1], Gamma[2]);
        }
    }

    printf("\n══════════════════════════════════════════════════════════════════\n");
    printf("  How to use this output\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");
    printf("  At any FIXED gap_floor, Γ is well-defined and converges as the\n");
    printf("  BZ grid is refined. The gap_floor → 0 limit is genuinely\n");
    printf("  divergent because the Born self-energy on a truly-gapless\n");
    printf("  Goldstone model is IR-divergent.\n\n");
    printf("  Pick gap_floor based on physical context:\n");
    printf("    - For a model with single-ion anisotropy K, use gap_floor\n");
    printf("      ≪ K (e.g., 1e-3 if K = 0.1·J).\n");
    printf("    - For a model in an applied field h, use gap_floor ≪ h·g·μ_B.\n");
    printf("    - For comparing to experimental linewidths, use the\n");
    printf("      experimental gap (instrument resolution + intrinsic T-\n");
    printf("      dependent gap from disorder + magnon-magnon scattering).\n");
    printf("    - For matching published Mourigal-Chernyshev-style maps,\n");
    printf("      use gap_floor consistent with their self-consistent or\n");
    printf("      experimental IR treatment (typically ~0.05-0.1·J).\n\n");
    printf("  The function is symmetry-correct under all standard checks:\n");
    printf("  Γ(k) = Γ(-k), bond-orientation invariant, U(1)-protected\n");
    printf("  zero on collinear-and-symmetric models, unitarity (Γ ≥ 0).\n\n");

    irrep_magnon_lsw_free(L);
    return 0;
}
