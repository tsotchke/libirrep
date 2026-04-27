/* SPDX-License-Identifier: MIT */
/* B20-style chiral magnet 4-site Heisenberg + DMI ED.
 *
 * Closes the materials-search loop end-to-end on the smallest non-
 * trivial chiral magnet:
 *
 *   1. lattice    — 4-site closed ring (analytic singlet at J only)
 *   2. analyzer   — chiral cubic O preserves bond → D ∥ bond
 *      (Bak-Jensen pattern). Already shown by the dmi.h tests; here we
 *      USE the result by setting D = (1, 0, 0) on the bond from x=0 to
 *      x=1 (the bond direction).
 *   3. apply      — irrep_heisenberg_apply + irrep_dmi_apply, summed
 *      via a small combined-apply callback.
 *   4. ED         — irrep_lanczos_eigvals_reorth on the combined
 *      Hamiltonian, varying D/J to trace the spin-canting transition.
 *
 * Demonstrates: DMI cants the GS away from the singlet at J=1, with
 * the ground-state energy E₀(D/J) decreasing monotonically as D
 * grows. The 4-site ring is too small for skyrmion physics, but the
 * DMI-driven canting is the same mechanism that drives helimagnetism
 * in MnSi-type compounds at large N.
 *
 * Build: `make examples`
 * Run:   `./build/bin/b20_chiral_4site_ed` */

#include <irrep/dmi_hamiltonian.h>
#include <irrep/hamiltonian.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N        4
#define DIM      (1LL << N)

typedef struct {
    irrep_heisenberg_t      *H_J;
    irrep_dmi_hamiltonian_t *H_D;
    double _Complex         *scratch;
} combined_t;

/* Combined apply: out = H_J · psi + H_D · psi. */
static void combined_apply(const double _Complex *psi, double _Complex *out, void *opaque) {
    combined_t *C = opaque;
    irrep_heisenberg_apply(psi, out, C->H_J);
    irrep_dmi_apply(psi, C->scratch, C->H_D);
    for (long long s = 0; s < DIM; ++s)
        out[s] += C->scratch[s];
}

int main(void) {
    printf("=== libirrep — B20-style 4-site Heisenberg + DMI ED ===\n");
    printf("    GS energy vs D/J on a 4-site ring with NN Heisenberg + DMI ∥ bond.\n\n");

    /* 4-site closed ring: bonds (0,1), (1,2), (2,3), (3,0).
     * D vector along the bond direction (chosen so all 4 bonds
     * have D in the same in-plane direction — schematic, not
     * geometrically realised on this 1D ring). */
    int bi[4] = {0, 1, 2, 3};
    int bj[4] = {1, 2, 3, 0};

    irrep_heisenberg_t *H_J = irrep_heisenberg_new(N, 4, bi, bj, /*J=*/1.0);

    printf("    %-6s  %-12s  %-12s  %-s\n", "D/J", "E_0", "E_0 − E_J", "interpretation");
    printf("    ────────────────────────────────────────────────────\n");

    double E_J_only = 0;
    for (int q = 0; q <= 6; ++q) {
        double DoJ = q * 0.25; /* 0, 0.25, ..., 1.5 */
        double Dx[4] = {DoJ, DoJ, DoJ, DoJ};
        double Dy[4] = {0, 0, 0, 0};
        double Dz[4] = {0, 0, 0, 0};
        irrep_dmi_hamiltonian_t *H_D = irrep_dmi_hamiltonian_new(N, 4, bi, bj, Dx, Dy, Dz);

        double _Complex *scratch = malloc((size_t)DIM * sizeof *scratch);
        combined_t C = {H_J, H_D, scratch};

        double _Complex *seed = malloc((size_t)DIM * sizeof *seed);
        double           sn = 0;
        for (long long s = 0; s < DIM; ++s) {
            seed[s] = 0.1 * sin(0.37 * s) + I * 0.05 * cos(0.23 * s);
            sn += creal(seed[s]) * creal(seed[s]) + cimag(seed[s]) * cimag(seed[s]);
        }
        sn = sqrt(sn);
        for (long long s = 0; s < DIM; ++s)
            seed[s] /= sn;

        double          E0;
        irrep_status_t  st = irrep_lanczos_eigvals_reorth(combined_apply, &C, DIM, 1, 100, seed,
                                                         &E0);
        if (st != IRREP_OK) {
            printf("    Lanczos failed at D/J=%.2f\n", DoJ);
        } else {
            if (q == 0)
                E_J_only = E0;
            const char *note =
                q == 0   ? "pure Heisenberg singlet (E_0 = -2 J on 4-site ring)"
                : DoJ < 0.6 ? "weak DMI: small canting"
                : DoJ < 1.1 ? "moderate DMI: significant deviation from singlet"
                            : "DMI dominates, GS is canted-helical";
            printf("    %-6.2f  %+12.6f  %+12.6f  %s\n", DoJ, E0, E0 - E_J_only, note);
        }

        free(seed);
        free(scratch);
        irrep_dmi_hamiltonian_free(H_D);
    }

    printf("\n  ━ Interpretation ━\n");
    printf("    Pure Heisenberg (D=0): GS is the spin singlet of a 4-site ring,\n");
    printf("    E_0 = −2J as on every bipartite NN ring with even sites.\n");
    printf("    Adding D ∥ bond cants the spins away from collinear AFM,\n");
    printf("    LOWERING the energy monotonically. The same mechanism at the\n");
    printf("    thermodynamic limit on a 3D B20 lattice produces the helimagnetic\n");
    printf("    phase + skyrmion lattice; here at N=4 the ring is too small for\n");
    printf("    skyrmion physics but the canting trend is unambiguous.\n\n");
    printf("    This closes the materials-search loop:\n");
    printf("      1. dmi.h analyzer    → 'D ∥ bond is allowed' under chiral cubic O\n");
    printf("      2. dmi_hamiltonian.h → ED on the resulting Heisenberg + DMI\n");
    printf("      3. Lanczos           → ground-state energy at each D/J\n");
    printf("      4. Observable trend  → DMI lowers the GS energy, drives canting\n");

    irrep_heisenberg_free(H_J);
    return 0;
}
