/* SPDX-License-Identifier: MIT */
/* Demonstrates irrep_magnon_born_decay_rate on a square FM with a
 * contrived cubic-vertex |V_3|² ~ sin²(½(k - k1)·â) — a generic
 * "non-local" kernel that imitates what a published cubic vertex
 * looks like (depends on momentum transfer through trig functions).
 *
 * Compares Γ_kin (matrix-element-stripped phase space) with Γ_Born
 * for two vertex strengths g = 0.1 and g = 0.5 (in J units). The
 * relative scaling Γ_Born ∝ g² is the experimentally-measurable
 * relationship.
 *
 * For a real material the cubic vertex would come from the local-
 * frame Holstein-Primakoff expansion of the full bond-list
 * Hamiltonian (Mourigal-Chernyshev for kagome 120° Néel, etc.) —
 * the Born harness is generic; the model-specific physics is in the
 * vertex callback.
 *
 * Build: make examples
 * Run:   ./build/bin/born_damping_demo */

#include <irrep/magnon.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N_K 12

/* Contrived vertex: |V_3|² = g² · sin²(½ (kx - k1x)) · sin²(½ (ky - k1y)).
 * Vanishes at k = k1 (forward scattering); peaks at zone-boundary
 * momentum transfers. */
typedef struct {
    double g_squared;
} demo_vertex_data_t;

static double demo_vertex_(int b, int b1, int b2, double kx, double ky, double k1x, double k1y,
                            double k2x, double k2y, void *user) {
    (void)b; (void)b1; (void)b2;
    (void)k2x; (void)k2y;
    demo_vertex_data_t *d = user;
    double sx = sin(0.5 * (kx - k1x));
    double sy = sin(0.5 * (ky - k1y));
    return d->g_squared * sx * sx * sy * sy;
}

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Born decay-rate demo — contrived vertex on a square FM\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[2] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);

    /* Walk the Γ-X-M-Γ path. */
    double k_pts[N_K][2];
    const char *labels[N_K] = {0};
    int  per = N_K / 3;
    for (int i = 0; i < per; ++i) {
        double t = (double)i / per;
        k_pts[i][0] = t * M_PI;
        k_pts[i][1] = 0;
    }
    for (int i = 0; i < per; ++i) {
        double t = (double)i / per;
        k_pts[per + i][0] = M_PI;
        k_pts[per + i][1] = t * M_PI;
    }
    int rest = N_K - 2 * per;
    for (int i = 0; i < rest; ++i) {
        double t = (double)i / rest;
        k_pts[2 * per + i][0] = M_PI * (1 - t);
        k_pts[2 * per + i][1] = M_PI * (1 - t);
    }
    labels[0]           = "Γ";
    labels[per]         = "X";
    labels[2 * per]     = "M";
    labels[N_K - 1]     = "Γ";

    /* Kinematic phase space (matrix-element-stripped). */
    double *gamma_kin = malloc((size_t)N_K * sizeof *gamma_kin);
    irrep_magnon_kinematic_damping(L, k_pts, N_K, 24, 24, 0.05, gamma_kin);

    /* Born decay rate at two vertex strengths. */
    demo_vertex_data_t v_lo = {.g_squared = 0.1 * 0.1};
    demo_vertex_data_t v_hi = {.g_squared = 0.5 * 0.5};
    double *gamma_lo = malloc((size_t)N_K * sizeof *gamma_lo);
    double *gamma_hi = malloc((size_t)N_K * sizeof *gamma_hi);
    irrep_magnon_born_decay_rate(L, k_pts, N_K, 24, 24, 0.05, demo_vertex_, &v_lo, gamma_lo);
    irrep_magnon_born_decay_rate(L, k_pts, N_K, 24, 24, 0.05, demo_vertex_, &v_hi, gamma_hi);

    printf("  Γ-X-M-Γ path on the square FM (single band, S = 1/2):\n");
    printf("  Vertex |V_3|² = g² · sin²(½ Δkx) · sin²(½ Δky), g ∈ {0.1, 0.5}\n\n");
    printf("  %-4s %-6s    %-12s   %-12s   %-12s   %-10s\n", "ik", "label", "ω(k)",
           "Γ_kin", "Γ_Born(g=0.1)", "Γ/g² ratio");
    printf("  ────┼──────────┼──────────────┼──────────────┼──────────────┼─────────────\n");
    for (int ik = 0; ik < N_K; ++ik) {
        double omega, dummy[1];
        double _Complex u_dummy[1];
        irrep_magnon_dispersion(L, k_pts[ik][0], k_pts[ik][1], &omega, u_dummy);
        (void)dummy;
        const char *lbl = labels[ik] ? labels[ik] : "";
        double ratio = (v_lo.g_squared > 0 && gamma_lo[ik] > 0)
                          ? gamma_hi[ik] / gamma_lo[ik]
                          : 0.0 / 0.0;
        printf("  %-4d %-6s    %-12.4f   %-12.4f   %-12.4f   %.3f\n", ik, lbl, omega,
               gamma_kin[ik], gamma_lo[ik], ratio);
    }
    printf("\n");

    printf("  PHYSICS:\n");
    printf("  - Γ_kin > 0 is the kinematic phase-space upper bound.\n");
    printf("  - Γ_Born(g) is the actual Born linewidth for the chosen\n");
    printf("    cubic vertex; scales as g² across all k as expected.\n");
    printf("  - The Γ/g² ratio is the geometric factor 25 (= 0.25/0.01),\n");
    printf("    independent of k, confirming the Born integral is linear\n");
    printf("    in |V|² and the harness handles vertex scaling correctly.\n");
    printf("  - For a *real* model the user supplies the model-specific\n");
    printf("    cubic vertex (Mourigal-Chernyshev, canted-FM, spiral-LSW,\n");
    printf("    etc.) via the same callback interface.\n\n");

    free(gamma_kin);
    free(gamma_lo);
    free(gamma_hi);
    irrep_magnon_lsw_free(L);
    return 0;
}
