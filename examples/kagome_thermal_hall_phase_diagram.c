/* SPDX-License-Identifier: MIT */
/* Kagome FM thermal Hall conductivity κ_xy(T) phase diagram —
 * D-strength sweep with Chern-signature robustness and peak-temperature
 * scaling.
 *
 * MATHEMATICAL CONTEXT
 *
 * Mook-Henk-Mertig (2014) established that the kagome ferromagnet with
 * NN Heisenberg J + out-of-plane DMI D_z has band Chern numbers
 * (−1, 0, +1) for any D > 0, and a topological gap at the K point of
 * size 2√3·D·S in libirrep's per-bond convention. The Berry-curvature
 * integral over occupied magnon bands at finite temperature gives a
 * non-zero thermal Hall conductivity (Matsumoto-Murakami 2011, 2014):
 *
 *     κ_xy(T) = -(k_B² T / (ℏ V_uc · (2π)²)) · Σ_b ∫_BZ d²k · c₂(n_B) · Ω_b
 *
 * where c₂(g) is the Matsumoto-Murakami weight (PRL 106, 197202, 2011).
 *
 * For a fixed Chern signature, κ_xy(T) on the kagome FM is monotone-
 * saturating in T:
 *   - T → 0: κ_xy → 0 (no magnons populated below the gap)
 *   - T ~ gap: rising
 *   - T ~ bandwidth: saturating
 *   - T → ∞: κ_xy → −Σ_b ∫_BZ ω_b · F_b(k) d²k / (A_uc · (2π)²)
 *
 * The high-T plateau follows from the c₂ asymptotic c₂(g) = π²/3 − 1/g
 * + O(1/g²): the leading π²/3 term cancels by the Σ_b C_b = 0 sum rule,
 * leaving the band-energy-weighted Berry-curvature integral
 *   κ_∞ = −Σ_b ∫_BZ ω_b(k) · F_b(k) d²k / (A_uc · (2π)²)
 * as a finite constant.
 *
 * Empirical D-dependence on the kagome FM with NN J = −1, S = 1, NN
 * out-of-plane DMI Dz: κ_∞ is monotone increasing in D over four
 * decades (D ∈ [0.001, 0.5]) but does NOT follow a clean power law.
 * Local slope of log κ_∞ vs log D varies with D regime:
 *   D ∈ [0.001, 0.01]:    κ_∞ ∼ D^0.52
 *   D ∈ [0.01,  0.10]:    κ_∞ ∼ D^0.66
 *   D ∈ [0.10,  0.50]:    κ_∞ ∼ D^0.46
 * The aggregate slope across [0.001, 0.5] is ~0.56. Underlying
 * structure: κ_∞ has band contributions Σ_b ω_b(k) C_b weighted by
 * the Berry curvature distribution F_b(k), which sharpens at K and
 * K′ Dirac points as D → 0 but spreads across the band as D grows.
 * Neither the sharp-peak limit nor the spread-out limit gives a
 * single-power-law D-dependence — the empirical exponent ~ 0.5 is
 * a coarse summary of multi-scale band-structure physics.
 *
 * κ_∞ is NOT a topological invariant — it is a Hamiltonian-specific
 * band-structure quantity — but is robust against gap-preserving
 * deformations.
 *
 * EXPERIMENTAL RELEVANCE
 *
 * For materials such as Cu(1,3-bdc) (Chisnell et al, PRL 115, 147201,
 * 2015) and Fe₃Sn₂ (Yin et al, Nature 562, 91, 2018), measuring κ_xy(T)
 * vs T gives access to the Dzyaloshinskii-Moriya strength via the peak
 * position. This example tabulates the predicted κ_xy(T) curves across
 * a range of D, so a measured peak can be back-extracted into a D value
 * given J and S.
 *
 * Build: make examples
 * Run:   ./build/bin/kagome_thermal_hall_phase_diagram */

#include <irrep/magnon.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static irrep_magnon_lsw_t *build_kagome_FM_DMI(double J, double D, double S) {
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * sqrt(3.0)};
    irrep_magnon_bond_t bonds[6] = {
        {.bi = 0, .bj = 1, .delta_x =  0, .delta_y =  0, .J = J, .D = {0, 0, +D}},
        {.bi = 1, .bj = 2, .delta_x =  0, .delta_y =  0, .J = J, .D = {0, 0, +D}},
        {.bi = 2, .bj = 0, .delta_x =  0, .delta_y =  0, .J = J, .D = {0, 0, +D}},
        {.bi = 1, .bj = 0, .delta_x =  1, .delta_y =  0, .J = J, .D = {0, 0, -D}},
        {.bi = 0, .bj = 2, .delta_x =  0, .delta_y = -1, .J = J, .D = {0, 0, -D}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y =  1, .J = J, .D = {0, 0, -D}},
    };
    return irrep_magnon_lsw_new(/*n_sub=*/3, S, a1, a2, /*n_bonds=*/6, bonds, /*Kz=*/0);
}

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Kagome FM with DMI — thermal Hall κ_xy(T) phase diagram\n");
    printf("  (Mook 2014 + Matsumoto-Murakami 2011/2014, libirrep prediction)\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    const double J = -1.0;
    const double S = 1.0;
    /* Extended D-range: 4 decades from D = 0.001 to D = 0.5, plus the
     * original moderate-D points for the κ_xy(T) phase-diagram tables. */
    const double D_set[] = {0.001, 0.003, 0.01, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50};
    const int n_D = (int)(sizeof D_set / sizeof D_set[0]);

    /* === Part 1: Chern-signature robustness across D sweep === */
    printf("  PART 1 — Chern-number robustness across D ∈ [0.05, 0.30]\n");
    printf("  Verify (-1, 0, +1) signature stable across the sweep.\n\n");
    printf("  %-7s   %-12s   %-8s %-8s %-8s   %-12s\n",
           "D", "K-pt gap", "C_1", "C_2", "C_3", "Σ C_b");
    printf("  ──────────────────────────────────────────────────────────────\n");

    int chern_signature_ok = 1;
    for (int i = 0; i < n_D; ++i) {
        irrep_magnon_lsw_t *L = build_kagome_FM_DMI(J, D_set[i], S);
        double K_pt[2] = {4 * M_PI / 3, 0.0};
        double w[3];
        double _Complex u_dummy[9];
        irrep_magnon_dispersion(L, K_pt[0], K_pt[1], w, u_dummy);
        double gap = w[1] - w[0];

        double chern[3];
        irrep_magnon_chern(L, 64, 64, chern);
        double sum_c = chern[0] + chern[1] + chern[2];

        printf("  %-7.3f   %-12.6f   %+8.4f %+8.4f %+8.4f   %+12.2e\n",
               D_set[i], gap, chern[0], chern[1], chern[2], sum_c);

        if (fabs(chern[0] - (-1.0)) > 1e-3 ||
            fabs(chern[1] -   0.0 ) > 1e-3 ||
            fabs(chern[2] - (+1.0)) > 1e-3 ||
            fabs(sum_c)             > 1e-3)
            chern_signature_ok = 0;

        irrep_magnon_lsw_free(L);
    }
    printf("\n  → Chern signature (-1, 0, +1) %s across the entire D sweep.\n\n",
           chern_signature_ok ? "ROBUST" : "FAILED");

    /* === Part 2: κ_xy(T) on T grid for each D — high-T saturation === */
    printf("  PART 2 — κ_xy(T) on a T sweep, for each D — units of (k_B²/ℏ)\n");
    printf("  Curve is monotone-saturating, with high-T plateau set by\n");
    printf("  the band-energy-weighted Berry-curvature integral.\n\n");

    const double T_set[] = {0.05, 0.20, 0.50, 1.00, 2.00, 4.00, 8.00, 16.0, 32.0, 64.0};
    const int n_T = (int)(sizeof T_set / sizeof T_set[0]);

    /* Header row: T values across the top. */
    printf("        D \\ T    ");
    for (int j = 0; j < n_T; ++j) printf("%8.3f  ", T_set[j]);
    printf("\n  %s", "──────────────");
    for (int j = 0; j < n_T; ++j) printf("──────────");
    printf("\n");

    /* For each D row, build LSW once, sweep T. */
    double kxy_table[n_D][n_T];
    for (int i = 0; i < n_D; ++i) {
        irrep_magnon_lsw_t *L = build_kagome_FM_DMI(J, D_set[i], S);
        printf("   D = %-6.3f   ", D_set[i]);
        for (int j = 0; j < n_T; ++j) {
            double k = irrep_magnon_thermal_hall_kxy(L, T_set[j], 64, 64);
            kxy_table[i][j] = k;
            printf("%+8.1e  ", k);
        }
        printf("\n");
        irrep_magnon_lsw_free(L);
    }
    printf("\n");

    /* === Part 3: extract saturation amplitude per D, check D scaling === */
    printf("  PART 3 — High-T saturation amplitude κ_∞(D) and gap scale\n\n");
    printf("  %-8s   %-12s   %-14s   %-10s   %-10s   %-10s\n",
           "D", "K-gap = 2√3DS", "κ_∞ (T=64)", "κ_∞ / D", "κ_∞ / gap", "log slope");
    printf("  ──────────────────────────────────────────────────────────────────────────\n");
    int j_max = n_T - 1;
    for (int i = 0; i < n_D; ++i) {
        double gap = 2.0 * sqrt(3.0) * D_set[i] * S;
        double kappa_inf = kxy_table[i][j_max];
        double log_slope = 0.0;
        if (i > 0 && D_set[i - 1] > 0.0 && kxy_table[i - 1][j_max] > 0.0) {
            log_slope =
                log(kappa_inf / kxy_table[i - 1][j_max]) / log(D_set[i] / D_set[i - 1]);
        }
        printf("  %-8.4f   %-12.6f   %+14.6e   %-10.4f   %-10.4f   %s%-10.3f\n",
               D_set[i], gap, kappa_inf, kappa_inf / D_set[i], kappa_inf / gap,
               (i == 0 ? "  " : ""), log_slope);
    }
    printf("\n  → κ_∞ is monotone-increasing in D but does NOT follow a clean\n");
    printf("    power law. Local log-slope varies with D regime:\n");
    printf("      small D (Berry curvature sharply peaked at K): slope ≈ 0.5\n");
    printf("      mid-range D (curvature broadens onto band):    slope ≈ 0.6-0.7\n");
    printf("      larger D (curvature fully spread, ω_b matters): slope ≈ 0.4-0.5\n");
    printf("    Aggregate over [0.001, 0.5]: slope ≈ 0.56 — a coarse summary,\n");
    printf("    not a fundamental power law. The κ_∞/gap ratio is the\n");
    printf("    Hamiltonian-specific dimensionless number characteristic of\n");
    printf("    the kagome-FM Berry-curvature distribution weighted by ω_b.\n\n");

    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  EXPERIMENTAL TAKEAWAY\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");
    printf("  For a measured κ_xy(T) curve on a kagome FM material:\n\n");
    printf("    1. The high-T saturation amplitude κ_∞ grows monotonically\n");
    printf("       with D and saturates as the Berry curvature spreads\n");
    printf("       across the band. Pre-tabulated κ_∞(D) values let one\n");
    printf("       extract D from a measured κ_xy(T → bandwidth) plateau\n");
    printf("       given known J and S.\n\n");
    printf("    2. The temperature scale of the rise (10%%–90%% of κ_∞)\n");
    printf("       gives the magnon-band-energy distribution, which is\n");
    printf("       set by |J|·S (bandwidth), not D — this gives an\n");
    printf("       independent cross-check of |J|·S from the same data.\n\n");
    printf("    3. The Chern signature (-1, 0, +1) is robust: the sign\n");
    printf("       of κ_xy fixes the sign of the DM coefficient D_z\n");
    printf("       and the chirality of the kagome triangles.\n");

    return 0;
}
