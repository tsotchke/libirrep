/* SPDX-License-Identifier: MIT */
/* MTC consistency-proof demo: prints the runtime residuals for the
 * two MTCs shipped by libirrep — Ising (non-abelian) and Z₂×Z₂
 * (abelian) — against every axiom they should satisfy.
 *
 * # What this demonstrates
 *
 * Each MTC has a closed web of algebraic / topological identities
 * that must hold to define a consistent anyon theory. The library
 * computes runtime "residuals" for each identity — the max
 * absolute deviation from the expected value. A residual < 1e-12
 * is a numerical proof that the identity holds at machine
 * precision.
 *
 * The demo prints the residuals side-by-side, showing the same
 * audit structure across both MTCs.
 */
#include <irrep/crane_yetter_ising.h>
#include <irrep/toric_mtc.h>

#include <stdio.h>

int main(void) {
    printf("MTC consistency-proof audit (residuals; ≤ 1e-12 indicates passing proof)\n");
    printf("========================================================================\n");

    printf("\nIsing MTC (non-abelian, 3 simples {1, σ, ψ}, D = 2, c = 1/2):\n");
    printf("  F-symbol unitarity on σσσσ block      : %12.3e\n",
           irrep_ising_F_unitarity_residual());
    printf("  twist θ_a = R^{aa}_? trace formula     : %12.3e\n",
           irrep_ising_twist_from_R_residual());
    printf("  S-matrix from twists                   : %12.3e\n",
           irrep_ising_S_from_twist_residual());
    printf("  Verlinde formula N = Σ S·S·S*/S₀       : %12.3e\n",
           irrep_ising_verlinde_residual());
    printf("  modular S² = I (self-dual ⇒ C = I)     : %12.3e\n",
           irrep_ising_S_squared_residual());
    printf("  CY connected-sum multiplicativity      : %12.3e\n",
           irrep_crane_yetter_connected_sum_residual());

    printf("\nZ₂×Z₂ toric MTC (abelian, 4 simples {1, e, m, ψ}, D = 2, c = 0):\n");
    printf("  twist θ_a = R^{aa}_1 (abelian formula) : %12.3e\n",
           irrep_toric_mtc_twist_from_R_residual());
    printf("  Verlinde formula                       : %12.3e\n",
           irrep_toric_mtc_verlinde_residual());
    printf("  modular S² = I                         : %12.3e\n",
           irrep_toric_mtc_S_squared_residual());
    printf("  CY connected-sum multiplicativity      : %12.3e\n",
           irrep_toric_mtc_connected_sum_residual());

    printf("\nWalker-Wang ground-state dim (state-sum count on 3D polytopes):\n");
    printf("                       Ising    Z₂×Z₂\n");
    printf("  3-simplex            %5lld    %5lld   (= D-driven for both)\n",
           irrep_ising_walker_wang_simplex3_full_count(),
           irrep_toric_mtc_walker_wang_simplex3_full_count());
    printf("  cube                 %5lld    %5lld   (MTC algebra visible)\n",
           irrep_ising_walker_wang_cube_full_count(),
           irrep_toric_mtc_walker_wang_cube_full_count());
    printf("  octahedron (dual)    %5lld    %5lld   (= cube count: WW duality)\n",
           irrep_ising_walker_wang_octahedron_full_count(),
           irrep_toric_mtc_walker_wang_octahedron_full_count());

    printf("\nAll residuals are runtime-computed; a value > 1e-12 would indicate\n");
    printf("an MTC axiom violation. The polyhedral-duality equality is exact\n");
    printf("(integer counts; same value on dual polyhedra is the runtime witness).\n");
    printf("========================================================================\n");
    return 0;
}
