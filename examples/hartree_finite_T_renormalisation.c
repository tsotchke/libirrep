/* SPDX-License-Identifier: MIT */
/* Bloch-Dyson Hartree-Fock renormalisation factor Z(T) of the
 * spin-½ square ferromagnet — exercises
 * `irrep_magnon_hartree_renormalisation` and demonstrates the
 * standard finite-T LSW-validity diagnostic.
 *
 * MATHEMATICAL CONTEXT
 *
 * Linear spin-wave theory truncates the Holstein-Primakoff boson
 * expansion at quadratic order. Including the leading quartic
 * (Hartree-Fock) self-energy gives the Bloch-Dyson renormalised
 * dispersion
 *
 *     ω_HF(k, T) = ω_LSW(k) · Z(T)
 *
 * where the scalar renormalisation factor is
 *
 *     Z(T) = 1 − ⟨n⟩(T) / S,
 *     ⟨n⟩(T) = (1 / N_BZ) Σ_b ∫_{BZ} n_BE(ω_b / T) d²k.
 *
 * Limits:
 *
 *   - T → 0:          Z → 1 (LSW exact)
 *   - T near T_c:     Z → 0 (LSW breakdown — the magnon bands
 *                      collapse onto themselves and quartic
 *                      interactions are no longer perturbative)
 *
 * Z is monotone decreasing in T and provides a *diagnostic*: if
 * Z > 0.7 at the temperature of interest, LSW is quantitatively
 * trustworthy; if Z < 0.5, beyond-LSW machinery (Schwinger-boson
 * MF, classical MC) is needed.
 *
 * The example sweeps Z(T) for a square FM with S = ½, J = -1 over
 * T ∈ [0.01, 2.0] and verifies:
 *
 *   1. Z(T → 0) → 1 monotonically.
 *   2. Z is strictly decreasing in T.
 *   3. Z drops below 0.7 by T ≈ 1 (T ≈ |J|S = 0.5 is the rough
 *      LSW-validity scale).
 *
 * REFERENCES
 *   - Bloch, Z. Phys. 61, 206 (1930) — magnetisation T^{3/2} law
 *   - Dyson, Phys. Rev. 102, 1217 (1956) — leading 1/S beyond-LSW
 *   - Holstein & Primakoff, Phys. Rev. 58, 1098 (1940)
 *
 * Build: make examples
 * Run:   ./build/bin/hartree_finite_T_renormalisation */

#include <irrep/magnon.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Bloch-Dyson Z(T) on the spin-½ square FM — finite-T LSW\n");
    printf("  renormalisation factor and validity diagnostic.\n");
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

    /* Temperature grid: span low-T (LSW exact) to T ~ |J|S where
     * Z is expected to drop into the warning range. */
    double T_vals[8] = {0.01, 0.05, 0.10, 0.25, 0.50, 1.00, 1.50, 2.00};
    int    Nx = 64, Ny = 64;

    printf("  Setup: square FM, S = ½, J = -1, BZ %d×%d.\n\n", Nx, Ny);
    printf("  %-8s   %-12s   %-30s\n", "T", "Z(T)", "regime");

    double Z_prev = 1.0 + 1e-12;
    int    monotone_fails = 0;
    double Z_at_T_unity = NAN;
    for (int i = 0; i < 8; ++i) {
        double Z = irrep_magnon_hartree_renormalisation(L, T_vals[i], Nx, Ny);
        const char *regime;
        if      (Z > 0.95) regime = "LSW exact (Z ≈ 1)";
        else if (Z > 0.70) regime = "LSW reliable (Z > 0.7)";
        else if (Z > 0.50) regime = "LSW marginal (0.5 < Z < 0.7)";
        else if (Z > 0.0)  regime = "LSW BREAKDOWN (Z < 0.5)";
        else               regime = "post-LSW (Z ≤ 0)";
        printf("  %-8.4f   %-12.6f   %-30s\n", T_vals[i], Z, regime);

        if (Z > Z_prev + 1e-9) ++monotone_fails;
        Z_prev = Z;
        if (fabs(T_vals[i] - 1.0) < 1e-9) Z_at_T_unity = Z;
    }
    printf("\n");

    /* Verifications. */
    double Z_lowT = irrep_magnon_hartree_renormalisation(L, 0.01, Nx, Ny);
    int    pass_lowT     = (Z_lowT > 0.99);
    int    pass_monotone = (monotone_fails == 0);
    int    pass_t_unity  = (!isnan(Z_at_T_unity) && Z_at_T_unity < 0.7);

    printf("  Verifications:\n");
    printf("    Z(0.01)  = %.6f  → expected > 0.99   %s\n",
           Z_lowT, pass_lowT ? "PASS" : "FAIL");
    printf("    monotone-decrease in T               %s (%d violations)\n",
           pass_monotone ? "PASS" : "FAIL", monotone_fails);
    printf("    Z(T = 1.0) = %.6f → expected < 0.70  %s\n",
           Z_at_T_unity, pass_t_unity ? "PASS" : "FAIL");
    printf("\n");

    /* Convergence sanity: Z should be near-independent of BZ grid
     * size at any fixed T (modulo low-k Goldstone-mode integration
     * floor — for 2D FM ⟨n⟩(T) has a soft IR contribution that
     * converges slowly with BZ grid in the very low-T regime). */
    printf("  BZ-grid convergence at T = 0.5:\n");
    for (int N = 32; N <= 128; N *= 2) {
        double Z = irrep_magnon_hartree_renormalisation(L, 0.5, N, N);
        printf("    Nx = Ny = %3d   Z = %.6f\n", N, Z);
    }
    printf("\n");

    irrep_magnon_lsw_free(L);

    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Z(T) reproduces Bloch-Dyson finite-T renormalisation: → 1 at\n");
    printf("  low T, monotone decreasing, drops into LSW-breakdown range\n");
    printf("  near T ≈ |J|S, providing a quantitative LSW-validity check.\n");
    printf("══════════════════════════════════════════════════════════════════\n");
    return (pass_lowT && pass_monotone && pass_t_unity) ? 0 : 1;
}
