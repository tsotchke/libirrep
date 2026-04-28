/* SPDX-License-Identifier: MIT */
/* Kinematic magnon-decay phase space on the kagome topological FM —
 * exercises irrep_magnon_kinematic_damping along Γ-K-M-Γ.
 *
 * Same model as kagome_topological_magnons.c (Cu(1,3-bdc)-like FM
 * kagome with NN J=-1 and out-of-plane DMI Dz = ±0.15 on alternating
 * triangles). Three magnon bands: the lowest is the gapless Goldstone,
 * the upper two get gapped by DMI.
 *
 * What's printed:
 *   - the 1-magnon dispersion ω_b(k) along Γ-K-M-Γ
 *   - the kinematic damping Γ_kin_b(k) at each k
 *   - regions where Γ_kin > 0 are *kinematically allowed* decay
 *     channels; regions where Γ_kin = 0 are kinematically forbidden
 *
 * Physics interpretation: even where Γ_kin > 0, the actual linewidth
 * is Γ = Γ_kin · |V_3|² where V_3 is the cubic vertex. For a
 * collinear FM with U(1) symmetry, V_3 = 0 by symmetry → no real
 * damping despite finite phase space. DMI breaks U(1) → V_3 ~ D, so
 * real linewidth scales as (D/J)² · Γ_kin.
 *
 * Build: make examples
 * Run:   ./build/bin/kagome_kinematic_damping */

#include <irrep/magnon.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N_K 30  /* Γ → K → M → Γ in 30 points */

static void build_path(double (*kpath)[2], const char **label, double *t_axis) {
    double K[2] = {4.0 * M_PI / 3.0, 0.0};
    double M[2] = {M_PI, M_PI / sqrt(3.0)};
    int    n_per = N_K / 3;
    for (int i = 0; i < n_per; ++i) {
        double t = (double)i / n_per;
        kpath[i][0] = t * K[0];
        kpath[i][1] = t * K[1];
    }
    for (int i = 0; i < n_per; ++i) {
        double t = (double)i / n_per;
        kpath[n_per + i][0] = K[0] + t * (M[0] - K[0]);
        kpath[n_per + i][1] = K[1] + t * (M[1] - K[1]);
    }
    int rest = N_K - 2 * n_per;
    for (int i = 0; i < rest; ++i) {
        double t = (double)i / rest;
        kpath[2 * n_per + i][0] = M[0] * (1 - t);
        kpath[2 * n_per + i][1] = M[1] * (1 - t);
    }
    for (int i = 0; i < N_K; ++i) t_axis[i] = (double)i / (N_K - 1);
    label[0]              = "Γ";
    label[n_per]          = "K";
    label[2 * n_per]      = "M";
    label[N_K - 1]        = "Γ";
}

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Kinematic magnon damping — kagome topological FM with DMI\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * sqrt(3.0)};
    double D = 0.15;
    irrep_magnon_bond_t bonds[6] = {
        {.bi = 0, .bj = 1, .delta_x = 0, .delta_y = 0,  .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 2, .delta_x = 0, .delta_y = 0,  .J = -1.0, .D = {0, 0, +D}},
        {.bi = 2, .bj = 0, .delta_x = 0, .delta_y = 0,  .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 0, .delta_x = 1,  .delta_y = 0,  .J = -1.0, .D = {0, 0, -D}},
        {.bi = 0, .bj = 2, .delta_x = 0,  .delta_y = -1, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y = 1,  .J = -1.0, .D = {0, 0, -D}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(3, 1.0, a1, a2, 6, bonds, 0);

    double (*kpath)[2] = malloc((size_t)N_K * sizeof *kpath);
    double *t_axis     = malloc((size_t)N_K * sizeof *t_axis);
    const char **label = calloc(N_K, sizeof *label);
    build_path(kpath, label, t_axis);

    /* Dispersion. */
    double *omega_all = malloc((size_t)N_K * 3 * sizeof *omega_all);
    double _Complex *u_all = malloc((size_t)N_K * 9 * sizeof *u_all);
    for (int ik = 0; ik < N_K; ++ik)
        irrep_magnon_dispersion(L, kpath[ik][0], kpath[ik][1], omega_all + ik * 3,
                                 u_all + ik * 9);

    /* Kinematic damping (slow — quadratic in BZ grid). 24×24 grid. */
    double *gamma_all = malloc((size_t)N_K * 3 * sizeof *gamma_all);
    irrep_magnon_kinematic_damping(L, kpath, N_K, 24, 24, 0.05, gamma_all);

    printf("  Path: Γ → K → M → Γ on the kagome BZ.\n");
    printf("  3 bands per k. Γ_kin in units of J = 1.\n\n");
    printf("  %-4s %-8s    %-10s   %-10s   %-10s   %-22s\n", "ik", "label", "ω₁(k)", "ω₂(k)",
           "ω₃(k)", "Γ_kin (band 1, 2, 3)");
    printf("  ─────┼──────────┼────────────┼────────────┼────────────┼──────────────────────\n");
    for (int ik = 0; ik < N_K; ++ik) {
        const char *lbl = label[ik] ? label[ik] : "";
        printf("  %-4d %-8s    %-10.4f   %-10.4f   %-10.4f   %.4f, %.4f, %.4f\n", ik, lbl,
               omega_all[ik * 3 + 0], omega_all[ik * 3 + 1], omega_all[ik * 3 + 2],
               gamma_all[ik * 3 + 0], gamma_all[ik * 3 + 1], gamma_all[ik * 3 + 2]);
    }
    printf("\n");

    /* Summary: where is Γ_kin maximised? */
    int    ik_max = 0, b_max = 0;
    double g_max = 0;
    for (int ik = 0; ik < N_K; ++ik)
        for (int b = 0; b < 3; ++b)
            if (gamma_all[ik * 3 + b] > g_max) {
                g_max = gamma_all[ik * 3 + b];
                ik_max = ik;
                b_max  = b;
            }
    printf("  Maximal Γ_kin = %.4f at ik=%d (label %s), band %d, ω = %.3f.\n", g_max,
           ik_max, label[ik_max] ? label[ik_max] : "—", b_max,
           omega_all[ik_max * 3 + b_max]);

    int    ik_min = 0, b_min = 0;
    double g_min = 1e9;
    for (int ik = 0; ik < N_K; ++ik)
        for (int b = 0; b < 3; ++b)
            if (gamma_all[ik * 3 + b] < g_min) {
                g_min  = gamma_all[ik * 3 + b];
                ik_min = ik;
                b_min  = b;
            }
    printf("  Minimal Γ_kin = %.4f at ik=%d (label %s), band %d, ω = %.3f.\n\n", g_min,
           ik_min, label[ik_min] ? label[ik_min] : "—", b_min,
           omega_all[ik_min * 3 + b_min]);

    printf("  PHYSICS:\n");
    printf("  - Γ_kin > 0 marks regions with kinematic decay phase space.\n");
    printf("    These bands have *available* 1→2 channels, but the actual\n");
    printf("    linewidth needs the cubic vertex too: Γ = Γ_kin · |V_3|².\n");
    printf("  - This particular model (collinear FM along ẑ with out-of-\n");
    printf("    plane DMI ẑ) RETAINS U(1)_ẑ symmetry — DMI ẑ commutes\n");
    printf("    with spin rotations around ẑ. So V_3 = 0 by symmetry and\n");
    printf("    the real linewidth is identically zero, despite the finite\n");
    printf("    kinematic phase space shown above. The kinematic estimate\n");
    printf("    is an UPPER BOUND, not the actual linewidth.\n");
    printf("  - To get non-zero V_3 on the kagome FM one needs U(1)\n");
    printf("    breaking: in-plane DMI components, single-ion anisotropy\n");
    printf("    K_x or K_y, an external in-plane field, or a non-collinear\n");
    printf("    ground state (canted, spiral, 120° kagome AFM).\n");
    printf("  - For materials with in-plane DMI Dxy ≪ Dz (e.g., Cu(1,3-bdc)\n");
    printf("    where in-plane components are weak), Γ_real ~ (Dxy/J)² ·\n");
    printf("    Γ_kin. Upper bound for Dxy = 0.15: (0.15)² · Γ_kin_max\n");
    printf("    = %.4f · J for the most-damped band/k.\n", 0.15 * 0.15 * g_max);
    printf("  - This example demonstrates the kinematic *infrastructure*;\n");
    printf("    real damping prediction needs the matrix-element side as\n");
    printf("    well, which is model-specific and not yet built into the\n");
    printf("    library.\n\n");

    free(kpath);
    free(t_axis);
    free(label);
    free(omega_all);
    free(u_all);
    free(gamma_all);
    irrep_magnon_lsw_free(L);
    return 0;
}
