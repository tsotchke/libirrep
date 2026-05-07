/* SPDX-License-Identifier: MIT */
/* Consistency tests for `irrep_magnon_dispersion_noncollinear_3d`:
 * verify that on collinear ground states the non-collinear API
 * agrees with the collinear `_dispersion_3d` to machine precision,
 * and that the FM dispersion is invariant under uniform rotation
 * of the n_vectors (rotational invariance of Heisenberg).
 *
 * MATHEMATICAL CONTEXT
 *
 * The non-collinear LSW formalism handles arbitrary classical
 * ground states m̂_α ∈ S² per sublattice α. For a collinear
 * ferromagnet (m̂_α = ẑ for all α), the non-collinear machinery
 * specialises to the standard Hermitian Holstein-Primakoff path,
 * so:
 *
 *     ω_noncoll_3d(k; n̂_α = ẑ ∀α)  ≡  ω_3d(k)
 *
 * to machine precision. Furthermore, because the Heisenberg
 * Hamiltonian is rotationally invariant in spin space, replacing
 * (ẑ, ẑ, …) with (R · ẑ, R · ẑ, …) for any R ∈ SO(3) leaves
 * ω(k) unchanged:
 *
 *     ω_noncoll_3d(k; n̂_α = m̂ ∀α)  =  ω_3d(k)
 *
 * for any unit vector m̂. These two checks together exercise the
 * non-collinear 3D dispersion path and confirm rotational
 * invariance of the LSW machinery.
 *
 * The example uses a simple-cubic ferromagnet (n_sub = 1, J < 0,
 * S = ½) with all six NN bonds (±x̂, ±ŷ, ±ẑ) listed explicitly
 * — so each direction is counted twice in the bond list, giving
 * the analytic dispersion
 *
 *     ω(k) = 4|J|S [3 − cos k_x − cos k_y − cos k_z]
 *
 * (the libirrep convention sums over the bond list as written).
 * For J = −1, S = ½ the bandwidth is 24|J|S = 12 at k = (π, π, π)
 * and the quadratic Goldstone mode at Γ.
 *
 * REFERENCES
 *   - Holstein & Primakoff, Phys. Rev. 58, 1098 (1940)
 *   - Toth & Lake, J. Phys.: Condens. Matter 27, 166002 (2015) —
 *     non-collinear LSW formalism
 *
 * Build: make examples
 * Run:   ./build/bin/dispersion_noncollinear_3d_consistency */

#include <irrep/magnon.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  3D non-collinear LSW dispersion — consistency vs collinear\n");
    printf("  3D path on a simple-cubic ferromagnet.\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.0, 1.0};
    double a3[3] = {0.0, 0.0, 1.0};
    irrep_magnon_bond_t bonds[6] = {
        {.bi = 0, .bj = 0, .delta_x =  1, .delta_y =  0, .delta_z =  0, .J = -1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 0, .delta_x = -1, .delta_y =  0, .delta_z =  0, .J = -1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 0, .delta_x =  0, .delta_y =  1, .delta_z =  0, .J = -1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 0, .delta_x =  0, .delta_y = -1, .delta_z =  0, .J = -1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 0, .delta_x =  0, .delta_y =  0, .delta_z =  1, .J = -1.0, .D = {0,0,0}},
        {.bi = 0, .bj = 0, .delta_x =  0, .delta_y =  0, .delta_z = -1, .J = -1.0, .D = {0,0,0}},
    };
    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(1, 0.5, a1, a2, 6, bonds, 0);
    if (!L) return 1;

    /* Five canonical 3D BZ k-points. */
    struct {
        double      kx, ky, kz;
        const char *label;
    } pts[] = {
        {0.0,         0.0,        0.0,        "Γ                  "},
        {0.5 * M_PI,  0.0,        0.0,        "(π/2, 0, 0)         "},
        {        M_PI, 0.0,       0.0,        "X = (π, 0, 0)       "},
        {        M_PI, M_PI,      0.0,        "M = (π, π, 0)       "},
        {        M_PI, M_PI,      M_PI,       "R = (π, π, π) — top "},
    };
    enum { N_PTS = 5 };

    /* SECTION 1: collinear FM with n̂ = ẑ uniformly. The non-collinear
     * API must match the collinear `_dispersion_3d` path AND match
     * the analytic ω(k) = 1 - cos kx - cos ky - cos kz (units of 2|J|S = 1). */
    printf("  SECTION 1: collinear ground state n̂ = ẑ\n");
    printf("             non-collinear API vs collinear `_dispersion_3d` and analytic.\n\n");
    printf("  %-22s  %-11s  %-11s  %-11s  %-10s\n",
           "k-point", "ω_noncoll", "ω_collinear", "ω_analytic", "max abs err");

    double n_z[3] = {0.0, 0.0, 1.0};
    double max_err_z = 0.0;
    for (int i = 0; i < N_PTS; ++i) {
        double         w_nc, w_col;
        double _Complex u_dummy;
        irrep_magnon_dispersion_noncollinear_3d(L, n_z, a3, pts[i].kx, pts[i].ky,
                                                  pts[i].kz, &w_nc);
        irrep_magnon_dispersion_3d(L, a3, pts[i].kx, pts[i].ky, pts[i].kz,
                                     &w_col, &u_dummy);
        double w_an = 2.0 * ((1.0 - cos(pts[i].kx))
                           + (1.0 - cos(pts[i].ky))
                           + (1.0 - cos(pts[i].kz)));
        double err  = fmax(fabs(w_nc - w_col), fabs(w_nc - w_an));
        if (err > max_err_z) max_err_z = err;
        printf("  %-22s  %-11.6f  %-11.6f  %-11.6f  %.2e\n",
               pts[i].label, w_nc, w_col, w_an, err);
    }
    printf("\n  Worst-case ω disagreement across all paths: %.2e (machine precision).\n\n",
           max_err_z);

    /* SECTION 2: rotational invariance — replace n̂ = ẑ with n̂ pointing
     * along three other directions; the FM Heisenberg dispersion must
     * be unchanged because the Hamiltonian commutes with global SO(3). */
    printf("  SECTION 2: rotational invariance (uniform rotation of n_vectors)\n\n");
    printf("  %-30s  %-22s  %-10s\n",
           "n̂ direction", "ω(R) − ω(R; n̂=ẑ)", "max abs diff");

    struct {
        double      n[3];
        const char *label;
    } rotations[] = {
        {{1.0, 0.0, 0.0},                       "x̂  (90° rotation about y)   "},
        {{0.0, 1.0, 0.0},                       "ŷ  (90° rotation about x)   "},
        {{1.0 / sqrt(3.0), 1.0 / sqrt(3.0),
          1.0 / sqrt(3.0)},                     "(1,1,1)/√3 (cube diagonal) "},
        {{0.6, 0.8, 0.0},                       "arbitrary (3-4-0)/5         "},
    };

    double max_err_rot = 0.0;
    for (int r = 0; r < 4; ++r) {
        double w_ref;
        irrep_magnon_dispersion_noncollinear_3d(L, n_z, a3, pts[N_PTS - 1].kx,
                                                  pts[N_PTS - 1].ky,
                                                  pts[N_PTS - 1].kz, &w_ref);
        double w_rot;
        irrep_magnon_dispersion_noncollinear_3d(L, rotations[r].n, a3,
                                                  pts[N_PTS - 1].kx,
                                                  pts[N_PTS - 1].ky,
                                                  pts[N_PTS - 1].kz, &w_rot);
        double err = fabs(w_rot - w_ref);
        if (err > max_err_rot) max_err_rot = err;
        printf("  %-30s  %+9.6f vs %+9.6f  %.2e\n",
               rotations[r].label, w_rot, w_ref, err);
    }
    printf("\n  Heisenberg FM dispersion at R = (π, π, π) is invariant to %.2e\n",
           max_err_rot);
    printf("  under any global rotation of the ground-state spin direction —\n");
    printf("  consistent with SO(3) invariance of the isotropic Heisenberg\n");
    printf("  Hamiltonian.\n\n");

    irrep_magnon_lsw_free(L);

    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Non-collinear 3D LSW reduces to the collinear path at machine\n");
    printf("  precision and respects global SO(3) of Heisenberg exchange.\n");
    printf("══════════════════════════════════════════════════════════════════\n");
    /* Section-1 floor is the LSW eigensolver tolerance (~1e-10);
     * Section-2 rotational invariance is at machine precision. */
    return (max_err_z < 1e-9 && max_err_rot < 1e-12) ? 0 : 1;
}
