/* SPDX-License-Identifier: MIT */
/* Full 1+2+3 magnon spectroscopy on a NN-only bipartite square AFM —
 * exercises all three orders of the new dynamical-structure-factor
 * stack on the same model.
 *
 * Method: at three q-points (zone-interior, zone-boundary, Néel
 * wavevector) compute:
 *   - 1-magnon S⁽¹⁾(q, ω) via _one_magnon_qomega_general
 *   - 2-magnon S⁽²⁾(q, ω) via _two_magnon_qomega_general
 *   - 3-magnon S⁽³⁾(q, ω) via _three_magnon_qomega_general
 * Reports peak energy, peak intensity, and ω-integrated weight for
 * each order.
 *
 * Physics interpretation: the orders dominate in different ω regions
 * - ω ∈ [0, ω_max ≈ 4]: 1-magnon dominates (sharp peaks at the AFM
 *   spinwave dispersion, captured exactly by LSW Bogoliubov)
 * - ω ∈ [ω_max, 2·ω_max]: 2-magnon continuum extends here
 *   (typically a shoulder visible in real INS data)
 * - ω ∈ [2·ω_max, 3·ω_max]: 3-magnon continuum is the only LSW-level
 *   contribution; the integrated weight scales as 1/(2S)² so it's
 *   small for S = 1/2 but non-zero
 *
 * Build: make examples
 * Run:   ./build/bin/square_afm_three_magnon */

#include <irrep/magnon.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const struct {
    const char *name;
    double      qx, qy;
} q_pts[] = {
    {"(π/4, 0)",  M_PI / 4, 0},
    {"(π/2, 0)",  M_PI / 2, 0},
    {"(π, 0)",    M_PI,     0},
};
static const int N_Q = 3;

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Square AFM — full 1+2+3 magnon spectroscopy on a single model\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    /* Bipartite supercell with NN J = +1, S = 1/2. Bandwidth ≈ 4. */
    double a1[2] = {1.0, 1.0};
    double a2[2] = {1.0, -1.0};
    irrep_magnon_bond_t bonds[4] = {
        {.bi = 0, .bj = 1, .delta_x = 0,  .delta_y = 0,  .J = +1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 1, .delta_x = -1, .delta_y = -1, .J = +1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 1, .delta_x = 0,  .delta_y = -1, .J = +1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 1, .delta_x = -1, .delta_y = 0,  .J = +1.0, .D = {0,0,0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(2, 0.5, a1, a2, 4, bonds, 0);
    int signs[2] = {+1, -1};

    int    n_omega = 100;
    double w_max   = 12.0;     /* spans 1, 2, and 3 magnon ranges */
    double dw      = w_max / n_omega;
    double eta     = 0.10;

    double q_arr[3][2];
    for (int i = 0; i < N_Q; ++i) {
        q_arr[i][0] = q_pts[i].qx;
        q_arr[i][1] = q_pts[i].qy;
    }

    double *I1 = malloc((size_t)N_Q * n_omega * sizeof *I1);
    double *I2 = malloc((size_t)N_Q * n_omega * sizeof *I2);
    double *I3 = malloc((size_t)N_Q * n_omega * sizeof *I3);

    irrep_magnon_one_magnon_qomega_general(L, signs, q_arr, N_Q, 0.0, w_max, n_omega, eta, I1);
    irrep_magnon_two_magnon_qomega_general(L, signs, q_arr, N_Q, 24, 24, 0.0, w_max, n_omega,
                                             eta, I2);
    /* 3-magnon costs O(N^4); use small grid. */
    irrep_magnon_three_magnon_qomega_general(L, signs, q_arr, N_Q, 8, 8, 0.0, w_max, n_omega,
                                               eta, I3);

    printf("  Energy axis: ω ∈ [0, %.1f] in J units, %d bins, η = %.2f.\n", w_max, n_omega,
           eta);
    printf("  Bandwidths:  1-magnon ≈ 4 = 2·z·J·S\n");
    printf("               2-magnon ≈ 8 = 2 × 1-magnon band-top\n");
    printf("               3-magnon ≈ 12 = 3 × 1-magnon band-top\n\n");

    for (int iq = 0; iq < N_Q; ++iq) {
        printf("  ─── q = %s ───\n", q_pts[iq].name);
        double total[3] = {0, 0, 0};
        int    peak[3]  = {0, 0, 0};
        for (int j = 0; j < n_omega; ++j) {
            total[0] += I1[iq * n_omega + j] * dw;
            total[1] += I2[iq * n_omega + j] * dw;
            total[2] += I3[iq * n_omega + j] * dw;
            if (I1[iq * n_omega + j] > I1[iq * n_omega + peak[0]]) peak[0] = j;
            if (I2[iq * n_omega + j] > I2[iq * n_omega + peak[1]]) peak[1] = j;
            if (I3[iq * n_omega + j] > I3[iq * n_omega + peak[2]]) peak[2] = j;
        }
        double w_peak[3];
        for (int k = 0; k < 3; ++k) w_peak[k] = (peak[k] + 0.5) * dw;
        printf("    %-12s  %-15s  %-15s\n", "order", "peak (ω, S)", "∫ S(q, ω) dω");
        printf("    ──────────────┼─────────────────┼─────────────────\n");
        printf("    1-magnon      (%5.2f, %6.3f)  %7.4f\n", w_peak[0],
               I1[iq * n_omega + peak[0]], total[0]);
        printf("    2-magnon      (%5.2f, %6.3f)  %7.4f\n", w_peak[1],
               I2[iq * n_omega + peak[1]], total[1]);
        printf("    3-magnon      (%5.2f, %6.3f)  %7.4f\n", w_peak[2],
               I3[iq * n_omega + peak[2]], total[2]);
        printf("    total weight ratio 1:2:3 = 1.000 : %.3f : %.4f\n\n",
               total[1] / total[0], total[2] / total[0]);
    }

    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  PHYSICS\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");
    printf("  The 1-, 2-, and 3-magnon orders dominate different parts of\n");
    printf("  the (q, ω) plane on the square AFM:\n\n");
    printf("  - 1-magnon: sharp peaks at ω_b(q) ≤ 4 — the conventional\n");
    printf("    spinwave dispersion captured by LSW Bogoliubov. Total\n");
    printf("    spectral weight 2S·n_sub = 1 per q (sum rule).\n");
    printf("  - 2-magnon: continuum extending to 2·ω_max ≈ 8. In real INS\n");
    printf("    on cuprates this appears as a shoulder above the sharp\n");
    printf("    1-magnon peaks. Weight scales as 1/(2S) relative to\n");
    printf("    1-magnon for moderate S — visible but smaller.\n");
    printf("  - 3-magnon: continuum to 3·ω_max ≈ 12. For S = 1/2 the\n");
    printf("    weight scales as 1/(2S)² so it's smaller still but\n");
    printf("    explains residual high-energy intensity in INS data.\n\n");
    printf("  Together these give the full LSW prediction for inelastic\n");
    printf("  neutron scattering, comparable directly to experiment after\n");
    printf("  intra-cell form-factor corrections via\n");
    printf("  `_structure_factor_general_with_form_factor`.\n\n");

    free(I1);
    free(I2);
    free(I3);
    irrep_magnon_lsw_free(L);
    return 0;
}
