/* SPDX-License-Identifier: MIT */
/* Topological magnons on a ferromagnetic kagome lattice with
 * Dzyaloshinskii-Moriya interaction.
 *
 * The kagome FM with NN Heisenberg J (FM, J < 0) and an out-of-plane
 * DMI Dz on every NN bond is the simplest model that hosts:
 *
 *   1. Three magnon bands (one per sublattice).
 *   2. Dirac touching points at K (and K') in the unfrustrated Heisenberg
 *      limit (Dz = 0).
 *   3. Gap opening at K with non-zero Berry curvature when Dz ≠ 0.
 *   4. Non-trivial band Chern numbers (typically (-1, 0, +1) for the
 *      three bands), giving a thermal-Hall effect κ_xy(T).
 *
 * This is the magnon analog of the Haldane model. Materials known to
 * realise it include Cu(1,3-bdc) (Chisnell et al. 2015) and Fe₃Sn₂
 * (Yin et al. 2018; both kagome FM with measurable thermal Hall).
 *
 * Convention for Dz signs on NN bonds: traverse each triangle
 * counterclockwise. Up-triangle (in-cell A→B→C→A): D = +Dz on each
 * bond. Down-triangle (CCW around the triangle that sits between unit
 * cells): D = −Dz when expressed in the sublattice ordering of the
 * source-in-cell-(0,0) convention. This sign assignment is consistent
 * with the C_3 rotation around each triangle centre and produces the
 * canonical (-1, 0, +1) Chern signature.
 *
 * Build: `make examples`
 * Run:   `./build/bin/kagome_topological_magnons` */

#include <irrep/magnon.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== libirrep — topological magnons on the kagome ferromagnet ===\n\n");
    printf("    Setup: 3-sublattice kagome (A, B, C), NN Heisenberg J=-1 (FM),\n");
    printf("           DMI Dz = ±D on alternating triangles (Mook 2014 convention).\n");
    printf("           Spin S = 1.\n\n");

    /* Kagome primitive cell: a1 = (1, 0), a2 = (1/2, √3/2). Sublattices:
     *   A at (0, 0), B at (1/2, 0), C at (1/4, √3/4). */
    double a1[2] = {1.0, 0.0};
    double a2[2] = {0.5, 0.5 * 1.7320508075688772};

    /* All 6 NN bonds (length 1/2):
     *   in-cell up-triangle, source-in-cell-(0,0):
     *     1. A → B    (delta=(0,0))
     *     2. B → C    (delta=(0,0))
     *     3. C → A    (delta=(0,0))    [ closes the up-triangle CCW ]
     *   inter-cell, going to the down-triangle:
     *     4. B → A    (delta=(+1, 0))   [ B(0,0) → A(1,0) ]
     *     5. A → C    (delta=(0,-1))    [ A(0,0) → C(0,-1) ]
     *     6. C → B    (delta=(-1, 1))   [ C(0,0) → B(-1,1) ]
     *
     * Up-tri CCW gets +Dz, down-tri CCW gets -Dz. */
    double D = 0.15;
    irrep_magnon_bond_t bonds[6] = {
        /* up-tri CCW: A→B→C→A */
        {.bi = 0, .bj = 1, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 1, .bj = 2, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        {.bi = 2, .bj = 0, .delta_x = 0, .delta_y = 0, .J = -1.0, .D = {0, 0, +D}},
        /* down-tri CCW (at the boundary, expressed source-in-cell-(0,0)) */
        {.bi = 1, .bj = 0, .delta_x = 1,  .delta_y = 0, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 0, .bj = 2, .delta_x = 0,  .delta_y = -1, .J = -1.0, .D = {0, 0, -D}},
        {.bi = 2, .bj = 1, .delta_x = -1, .delta_y = 1, .J = -1.0, .D = {0, 0, -D}},
    };

    irrep_magnon_lsw_t *L = irrep_magnon_lsw_new(/*n_sub=*/3, /*S=*/1.0, a1, a2,
                                                  /*n_bonds=*/6, bonds, /*Kz=*/0);
    if (!L) {
        fprintf(stderr, "magnon LSW handle allocation failed\n");
        return 1;
    }

    /* (1) Dispersion along Γ → K → M → Γ.
     *   K  = (4π/3, 0)
     *   M  = (π, π/√3)
     *   Γ  = (0, 0) */
    printf("    Dispersion along Γ → K → M → Γ:\n");
    printf("    %-8s  %-9s  %-9s  %-9s  %-9s\n", "label", "kx", "ky", "ω₁", "ω₂ ω₃");
    printf("    ────────────────────────────────────────────────────────────\n");
    struct {
        const char *name;
        double      kx, ky;
    } kpath[] = {
        {"Γ", 0.0, 0.0},
        {"K", 4.0 * M_PI / 3.0, 0.0},
        {"M", M_PI, M_PI / 1.7320508075688772},
        {"Γ", 0.0, 0.0},
    };
    double          omega[3];
    double _Complex u[9];
    for (size_t i = 0; i < sizeof(kpath) / sizeof(*kpath); ++i) {
        irrep_magnon_dispersion(L, kpath[i].kx, kpath[i].ky, omega, u);
        printf("    %-8s  %+9.4f  %+9.4f  %9.4f  %9.4f  %9.4f\n", kpath[i].name, kpath[i].kx,
               kpath[i].ky, omega[0], omega[1], omega[2]);
    }

    /* (2) Berry curvature concentration near K. With D ≠ 0 the Dirac
     * touching at K is gapped, so Berry curvature accumulates around K. */
    printf("\n    Berry curvature Ω(k) near K = (4π/3, 0):\n");
    printf("    %-12s  %-9s  %-9s  %-9s  %-9s\n", "Δk from K", "Ω₁", "Ω₂", "Ω₃", "|Ω|/π²");
    printf("    ────────────────────────────────────────────────────────\n");
    double K_kx = 4.0 * M_PI / 3.0;
    double K_ky = 0.0;
    for (int q = 0; q <= 4; ++q) {
        double dq = q * 0.1;
        double berry[3];
        irrep_magnon_berry(L, K_kx + dq, K_ky, /*delta_k=*/1e-3, berry);
        printf("    %+12.4f  %+9.4f  %+9.4f  %+9.4f  %+9.4f\n", dq, berry[0], berry[1], berry[2],
               (berry[0] + berry[1] + berry[2]) / (M_PI * M_PI));
    }

    /* (3) Chern numbers — integrate Ω over the BZ. */
    printf("\n    Chern numbers (Nx × Ny = 64 × 64):\n");
    double chern[3];
    irrep_magnon_chern(L, 64, 64, chern);
    printf("    Band 1 (lower): C₁ = %+8.4f\n", chern[0]);
    printf("    Band 2 (mid):   C₂ = %+8.4f\n", chern[1]);
    printf("    Band 3 (upper): C₃ = %+8.4f\n", chern[2]);
    double sum = chern[0] + chern[1] + chern[2];
    printf("    Sum (should be 0 within rounding): %+8.4f\n", sum);

    /* (4) Thermal Hall conductivity κ_xy(T) — Matsumoto-Murakami. The
     * peak T scales with the topological gap O(D). */
    printf("\n    Thermal Hall conductivity κ_xy(T) (natural units, Nx×Ny = 32×32):\n");
    printf("    %-9s  %-14s\n", "T/J", "κ_xy");
    printf("    ─────────────────────────\n");
    double Ts[] = {0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0};
    for (size_t i = 0; i < sizeof(Ts) / sizeof(*Ts); ++i) {
        double k = irrep_magnon_thermal_hall_kxy(L, Ts[i], 32, 32);
        printf("    %-9.3f  %+14.6e\n", Ts[i], k);
    }

    printf("\n  ━ Interpretation ━\n");
    printf("    With Dz ≠ 0, the Dirac touchings at K, K' are gapped and each\n");
    printf("    band carries a non-zero integer Chern number. The sum is zero\n");
    printf("    (any closed Hilbert space has total Chern = 0); the (-1, 0, +1)\n");
    printf("    pattern is the canonical Mook-Henk-Mertig signature observed in\n");
    printf("    Cu(1,3-bdc) (Chisnell 2015) and predicted in Fe₃Sn₂ (Yin 2018).\n\n");
    printf("    The thermal Hall κ_xy(T) decays exponentially as T → 0 (no\n");
    printf("    magnons populated below the topological gap of O(D·sin(2π/3))\n");
    printf("    ≈ 0.13 in J units), grows steeply through the gap-crossing\n");
    printf("    regime, and asymptotes to a finite plateau at T >> bandwidth.\n");
    printf("    The plateau value is the magnon-Berry-energy integral\n");
    printf("    (1/V_uc) Σ_b ∫ Ω_b·ω_b — the leading O(T) term in κ_xy\n");
    printf("    cancels because Σ_b C_b = 0, but the next-order term ~ ω/T\n");
    printf("    leaves a finite residual.\n\n");
    printf("    The crossover temperature is set by the topological gap, so\n");
    printf("    a measurement of κ_xy(T) directly fixes the gap and thereby\n");
    printf("    the DMI strength D. Inelastic-neutron + thermal-transport\n");
    printf("    are the two halves of the experimental fingerprint.\n\n");
    printf("    Closes another loop: from a libirrep DMI bond list (J + Dz on\n");
    printf("    each NN bond, signs alternating by triangle parity) to a\n");
    printf("    *measurable* observable — κ_xy(T) — whose entire shape is\n");
    printf("    fixed by the magnon Chern numbers + dispersion + DMI gap.\n");

    irrep_magnon_lsw_free(L);
    return 0;
}
