/* SPDX-License-Identifier: MIT */
/* Inelastic-neutron-scattering Q-ω map for the spin-½ square
 * ferromagnet — analytic-dispersion verification of
 * `irrep_magnon_neutron_qomega_map` along the standard high-
 * symmetry k-path, plus the LSW Lorentzian sum-rule.
 *
 * MATHEMATICAL CONTEXT
 *
 * The single-sublattice nearest-neighbour Heisenberg ferromagnet
 * on the square lattice with J < 0 (in our convention; ferromagnetic
 * means negative bond, since the LSW handle minimises ⟨S_i · S_j⟩
 * along bonds with positive J) has the textbook Holstein-Primakoff
 * dispersion
 *
 *     ω(q) = 2 |J| S [ 2 − cos(q_x) − cos(q_y) ]
 *
 * which is quadratic near the Goldstone point Γ (ω ≈ |J|·S·|q|² in
 * 2D), in contrast to the linear AFM dispersion. Spin S = ½ and
 * J = 1 give bandwidth 8|J|S = 4 with a maximum at the M point
 * (q = (π, π)).
 *
 * The 1-magnon dynamical structure factor for this state, after
 * Lorentzian broadening with half-width η, is
 *
 *     S(q, ω) = (1/π) · [ S_⊥(q) · η / ((ω − ω(q))² + η²) ]
 *
 * with S_⊥(q) = 2S = 1 at every q (single-band saturated FM,
 * polarised along ẑ). The neutron form factor for an unpolarised
 * sample is the geometric average along the polarisation-orthogonal
 * direction, which collapses to the perpendicular structure factor
 * up to a kinematic prefactor.
 *
 * Two analytic checks the example performs:
 *
 *   1. **Peak tracking**: at each q, the energy bin with maximum
 *      intensity must coincide with ω(q) to within ~η.
 *   2. **Lorentzian unit-area sum rule**: integrating I(q, ω) dω
 *      over a window wide compared to η must give 2S = 1 to high
 *      precision.
 *
 * REFERENCES
 *   - Holstein & Primakoff, Phys. Rev. 58, 1098 (1940)
 *   - Anderson, Phys. Rev. 86, 694 (1952) — AFM analog
 *   - Squires, *Theoretical Introduction to Thermal Neutron
 *     Scattering* (Cambridge, 1978) — INS sum rules
 *
 * Build: make examples
 * Run:   ./build/bin/square_fm_neutron_ins */

#include <irrep/magnon.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Square-lattice S=½ FM — INS Q-ω map vs analytic dispersion\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    /* Single-sublattice ferromagnet — primitive square cell. */
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    irrep_magnon_bond_t bonds[2] = {
        {.bi = 0, .bj = 0, .delta_x = 1, .delta_y = 0, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
        {.bi = 0, .bj = 0, .delta_x = 0, .delta_y = 1, .delta_z = 0,
         .J = -1.0, .D = {0, 0, 0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 2, bonds, 0);
    if (!L) {
        fprintf(stderr, "magnon LSW handle allocation failed\n");
        return 1;
    }

    /* Standard high-symmetry k-path Γ → X → M → Γ on the 2D BZ.
     * Γ = (0,0), X = (π,0), M = (π,π) — 7 sample points. */
    enum { N_Q = 7 };
    double qpath[N_Q][2] = {
        {0.10,         0.0       },  /* near Γ */
        {0.5 * M_PI,   0.0       },
        {        M_PI, 0.0       },  /* X */
        {        M_PI, 0.5 * M_PI},
        {        M_PI,        M_PI},  /* M */
        {0.5 * M_PI,   0.5 * M_PI},
        {0.10,         0.10      },  /* near Γ along (1,1) */
    };
    const char *labels[N_Q] = {
        "≈ Γ",
        "Γ→X (½)",
        "X",
        "X→M (½)",
        "M",
        "M→Γ (½)",
        "≈ Γ (1,1)",
    };

    /* Energy axis: bandwidth is 8|J|S = 4. We pad the window
     * generously on both sides so Lorentzian tails are negligible:
     * with η = 0.02, the truncation error of a unit Lorentzian
     * outside [-2, 6] is (1/π)·atan(η/(w_pad)) summed over both
     * tails — at the band-edge ω(M) = 4 the lower tail is 6/η = 300
     * and upper tail 2/η = 100 from peak, both well into asymptotic
     * arctan ≈ π/2 region; net truncation ≈ 3e-3. */
    int    n_w   = 1600;
    double w_min = -2.0;
    double w_max = 6.0;
    double dw    = (w_max - w_min) / n_w;
    double eta   = 0.02;

    double *I_qw = malloc((size_t)N_Q * (size_t)n_w * sizeof *I_qw);
    if (!I_qw) {
        irrep_magnon_lsw_free(L);
        return 1;
    }

    irrep_status_t st = irrep_magnon_neutron_qomega_map(
        L, qpath, N_Q, w_min, w_max, n_w, eta, I_qw);
    if (st != IRREP_OK) {
        fprintf(stderr, "neutron_qomega_map failed: status %d\n", (int)st);
        free(I_qw);
        irrep_magnon_lsw_free(L);
        return 1;
    }

    printf("  Setup: 1-sublattice ferromagnet, J = -1 (FM convention),\n");
    printf("         S = 1/2, bandwidth 8|J|S = 4.\n");
    printf("         Q-ω grid: %d ω-bins on [%.1f, %.1f], η = %.3f.\n\n",
           n_w, w_min, w_max, eta);

    printf("  RESULT 1 — peak tracking vs analytic ω(q) = 2 - cos qx - cos qy:\n\n");
    printf("  %-12s  %-18s  %-10s  %-10s  %-12s\n",
           "k-point", "(qx, qy)", "ω_peak", "ω_exact", "|Δω| / η");
    int peak_hits = 0;
    for (int iq = 0; iq < N_Q; ++iq) {
        int jmax = 0;
        for (int j = 1; j < n_w; ++j)
            if (I_qw[iq * n_w + j] > I_qw[iq * n_w + jmax])
                jmax = j;
        double w_peak     = w_min + (jmax + 0.5) * dw;
        double w_expected = 2.0 - cos(qpath[iq][0]) - cos(qpath[iq][1]);
        double err_in_eta = fabs(w_peak - w_expected) / eta;
        if (err_in_eta < 3.0)
            ++peak_hits;
        printf("  %-12s  (%5.3f, %5.3f)     %-10.4f  %-10.4f  %-12.3f\n",
               labels[iq], qpath[iq][0], qpath[iq][1],
               w_peak, w_expected, err_in_eta);
    }
    printf("\n  All %d/%d peaks within 3η of analytic ω(q): %s\n\n",
           peak_hits, (int)N_Q,
           peak_hits == (int)N_Q ? "PASS" : "FAIL");

    printf("  RESULT 2 — Lorentzian unit-area sum rule:\n");
    printf("             ∫ I(q, ω) dω = Σ_b S_⊥_b(q) = 2S = 1 (single-band FM)\n\n");
    printf("  %-12s  %-18s  %-12s  %-12s\n",
           "k-point", "(qx, qy)", "∫I dω", "rel. err");
    int    sum_hits = 0;
    double max_err  = 0.0;
    for (int iq = 0; iq < N_Q; ++iq) {
        double sum = 0.0;
        for (int j = 0; j < n_w; ++j)
            sum += I_qw[iq * n_w + j] * dw;
        double rel = fabs(sum - 1.0);
        if (rel > max_err) max_err = rel;
        if (rel < 5e-3) ++sum_hits;
        printf("  %-12s  (%5.3f, %5.3f)     %-12.6f  %-12.2e\n",
               labels[iq], qpath[iq][0], qpath[iq][1], sum, rel);
    }
    printf("\n  All %d/%d q-points satisfy ∫ I dω = 1 to 5e-3 (max %.2e): %s\n\n",
           sum_hits, (int)N_Q, max_err,
           sum_hits == (int)N_Q ? "PASS" : "FAIL");

    free(I_qw);
    irrep_magnon_lsw_free(L);

    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  INS prediction matches analytic Holstein-Primakoff dispersion;\n");
    printf("  Lorentzian sum rule recovers 2S = 1 across the full k-path.\n");
    printf("══════════════════════════════════════════════════════════════════\n");
    return (peak_hits == (int)N_Q && sum_hits == (int)N_Q) ? 0 : 1;
}
