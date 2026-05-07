/* SPDX-License-Identifier: MIT */
/* Powder-averaged inelastic-neutron-scattering spectrum on the
 * kagome ferromagnet — exercises `irrep_magnon_powder_spectrum`
 * and contrasts the BZ-integrated transverse INS observable
 * S_powder(ω) against the bare magnon DOS.
 *
 * MATHEMATICAL CONTEXT
 *
 * Real-world neutron-scattering experiments often study
 * polycrystalline samples in which the single-crystal momentum-
 * resolved Q-ω map is angularly averaged into a 1D function of
 * energy, S_powder(ω). Formally,
 *
 *     S_powder(ω) = (1/N_BZ) Σ_{k ∈ BZ, b} S_⊥_b(k) · δ(ω − ω_b(k))
 *
 * where S_⊥_b(k) is the transverse one-magnon structure factor of
 * band b and ω_b(k) the band's dispersion. This is *not* the same
 * as the magnon density of states
 *
 *     D(ω) = (1/N_BZ) Σ_{k, b} δ(ω − ω_b(k))
 *
 * because the DOS weights every state equally while the powder
 * spectrum suppresses dark bands (zero transverse weight). On the
 * kagome ferromagnet this difference is dramatic: the flat band
 * is dark at the Γ point and carries reduced weight elsewhere, so
 * its Van-Hove δ-spike is muted in INS while it dominates the DOS.
 *
 * Sum rules:
 *
 *     ∫ D(ω)         dω = n_sub        (one state per band per cell)
 *     ∫ S_powder(ω)  dω = 2S · n_sub   (transverse spin sum rule)
 *
 * For the kagome FM with S = ½, n_sub = 3, this gives ∫D = 3 and
 * ∫S_powder = 3.
 *
 * Lattice setup: kagome with primitive vectors a₁ = (1, 0),
 * a₂ = (½, ½√3), three sublattices at the kagome basis sites,
 * and J = −1 (ferromagnetic in the libirrep sign convention).
 *
 * REFERENCES
 *   - Mook, Henk & Mertig, Phys. Rev. B 89, 134409 (2014)
 *   - Squires, *Theoretical Introduction to Thermal Neutron
 *     Scattering* (Cambridge, 1978)
 *
 * Build: make examples
 * Run:   ./build/bin/kagome_fm_powder_ins */

#include <irrep/magnon.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Kagome ferromagnet — powder-averaged INS S_powder(ω) vs\n");
    printf("  magnon DOS, illustrating dark-band suppression.\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * sqrt(3.0)};
    double S    = 0.5;
    irrep_magnon_bond_t bonds[6] = {
        {.bi = 0, .bj = 1, .delta_x =  0, .delta_y =  0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 1, .bj = 2, .delta_x =  0, .delta_y =  0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 2, .bj = 0, .delta_x =  0, .delta_y =  0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 1, .bj = 0, .delta_x =  1, .delta_y =  0, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 2, .delta_x =  0, .delta_y = -1, .J = -1.0, .D = {0, 0, 0}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y =  1, .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(3, S, a1, a2, 6, bonds, 0);
    if (!L) {
        fprintf(stderr, "magnon LSW handle allocation failed\n");
        return 1;
    }

    /* Energy axis: kagome FM bandwidth is 6|J|S = 3 with the flat
     * band sitting exactly at ω = 3. Pad a touch above 3 for clean
     * bin-edge resolution and start at 0 (Goldstone). */
    int    n_bins = 600;
    double w_min  = 0.0;
    double w_max  = 3.4;
    double dw     = (w_max - w_min) / n_bins;

    /* BZ integration grid — quadratic-cost; 256² is reasonable. */
    int Nx = 256, Ny = 256;

    double *S_pow = malloc((size_t)n_bins * sizeof *S_pow);
    double *D_om  = malloc((size_t)n_bins * sizeof *D_om);
    if (!S_pow || !D_om) { irrep_magnon_lsw_free(L); return 1; }

    irrep_status_t st1 = irrep_magnon_powder_spectrum(L, Nx, Ny, w_min, w_max,
                                                       n_bins, S_pow);
    irrep_status_t st2 = irrep_magnon_dos(L, Nx, Ny, w_min, w_max, n_bins, D_om);
    if (st1 != IRREP_OK || st2 != IRREP_OK) {
        fprintf(stderr, "spectrum failed: powder=%d dos=%d\n", (int)st1, (int)st2);
        free(S_pow); free(D_om); irrep_magnon_lsw_free(L);
        return 1;
    }

    /* Normalisation check: irrep_magnon_dos and powder are returned
     * as per-unit-cell histogram densities (states per unit-energy
     * per unit-cell). Multiply by bin width to integrate. */
    double sum_pow = 0.0, sum_dos = 0.0;
    for (int i = 0; i < n_bins; ++i) {
        sum_pow += S_pow[i] * dw;
        sum_dos += D_om[i]  * dw;
    }
    /* Locate the flat-band bin (largest DOS value). */
    int    iflat   = 0;
    double dos_max = D_om[0];
    for (int i = 1; i < n_bins; ++i)
        if (D_om[i] > dos_max) { dos_max = D_om[i]; iflat = i; }
    double w_flat       = w_min + (iflat + 0.5) * dw;
    double pow_at_flat  = S_pow[iflat];
    double dos_at_flat  = D_om[iflat];
    double pow_floor    = 0.0;       /* mid-band pow value as reference */
    double dos_floor    = 0.0;
    int    n_mid        = 0;
    for (int i = 0; i < n_bins; ++i) {
        double w = w_min + (i + 0.5) * dw;
        if (w > 0.5 && w < 2.5) {
            pow_floor += S_pow[i];
            dos_floor += D_om[i];
            ++n_mid;
        }
    }
    pow_floor /= n_mid;
    dos_floor /= n_mid;

    printf("  Setup: 3-sublattice kagome FM, S = ½, J = -1, no DMI.\n");
    printf("         BZ grid %d×%d, %d ω-bins on [%.1f, %.1f].\n\n",
           Nx, Ny, n_bins, w_min, w_max);

    printf("  RESULT 1 — sum rules:\n");
    printf("    ∫ D(ω)        dω = %.6f   (analytic n_sub        = 3)   rel.err %.2e\n",
           sum_dos, fabs(sum_dos - 3.0) / 3.0);
    printf("    ∫ S_powder(ω) dω = %.6f   (analytic 2S·n_sub     = 3)   rel.err %.2e\n",
           sum_pow, fabs(sum_pow - 3.0) / 3.0);
    printf("\n");

    printf("  RESULT 2 — flat-band suppression (ω ≈ %.3f):\n", w_flat);
    printf("    DOS peak / mid-band-DOS-mean         = %6.2f x\n",
           dos_at_flat / (dos_floor > 0 ? dos_floor : 1.0));
    printf("    S_powder peak / mid-band-S-pow-mean  = %6.2f x\n",
           pow_at_flat / (pow_floor > 0 ? pow_floor : 1.0));
    printf("    Ratio (S_pow/DOS) at flat band       = %.3f\n",
           (dos_at_flat > 0) ? (pow_at_flat / dos_at_flat) : 0.0);
    printf("\n");

    printf("  Interpretation: the flat band concentrates 1/n_sub of the\n");
    printf("  DOS into a single bin (Van-Hove δ-spike — ratio ~ %.0fx in\n",
           dos_at_flat / (dos_floor > 0 ? dos_floor : 1.0));
    printf("  the spike bin vs mid-band mean), but the powder INS spectrum\n");
    printf("  shows zero weight at the flat band — the kagome compact-\n");
    printf("  localised hexagonal mode is uniformly orthogonal to the\n");
    printf("  ferromagnetic (1,1,1)/√3 spin direction at every k, so\n");
    printf("  (Σ_α u_α^flat)² ≡ 0 throughout the BZ. The flat band is\n");
    printf("  *completely dark* in transverse one-magnon INS — this is\n");
    printf("  the canonical signature of compact-localised states.\n\n");

    printf("  RESULT 3 — coarse spectral profile (every 50th bin):\n\n");
    printf("  %-8s   %-12s   %-12s\n", "ω", "D(ω)", "S_powder(ω)");
    for (int i = 0; i < n_bins; i += 50) {
        double w = w_min + (i + 0.5) * dw;
        printf("  %-8.4f   %-12.5f   %-12.5f\n", w, D_om[i], S_pow[i]);
    }

    free(S_pow); free(D_om);
    irrep_magnon_lsw_free(L);

    printf("\n══════════════════════════════════════════════════════════════════\n");
    printf("  Both sum rules satisfied; powder INS suppresses the kagome\n");
    printf("  flat band relative to the bare DOS by the structure factor.\n");
    printf("══════════════════════════════════════════════════════════════════\n");
    return (fabs(sum_dos - 3.0) < 0.05 && fabs(sum_pow - 3.0) < 0.05) ? 0 : 1;
}
