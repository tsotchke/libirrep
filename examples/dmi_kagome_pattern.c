/* SPDX-License-Identifier: MIT */
/* Bond-exchange-tensor decomposition for the 2D kagome lattice under
 * the in-plane D_6 site symmetry — directly applicable to the
 * room-temperature kagome magnet family:
 *
 *   Fe₃Sn₂          T_C ≈ 660 K, RT skyrmion-like bubbles, abundant Fe+Sn.
 *   Mn₃Sn / Mn₃Ge   T_N ≈ 420 / 380 K, AFM kagome with topological Hall effect.
 *   Co₃Sn₂S₂        T_C ≈ 177 K, magnetic Weyl semimetal, kagome layers.
 *
 * The 2D in-plane symmetry of these layered crystals is D_6 (or D_6h in
 * 3D, but the σ_h adds only the same constraint on J^s as the in-plane
 * elements for the symmetric tensor; for DMI it permits the out-of-plane
 * D_z component or forbids it depending on layer stacking). This demo
 * uses the 2D D_6 group from `point_group.h`, applied to a kagome NN
 * bond with positions padded to z = 0.
 *
 * Result derived from group theory alone — no input from material data:
 *
 *   - D_6 (chiral hexagonal, σ_v's containing the bond): see output below.
 *   - For a centrosymmetric kagome (D_6h, includes inversion), the in-
 *     plane DMI is forbidden by the inversion through the bond midpoint
 *     when the kagome plane has a true σ_h. In the layered metallic
 *     kagome magnets (Fe₃Sn₂ etc.), σ_h is broken by the stacking, so
 *     in-plane DMI survives — the toolkit predicts it from the broken-
 *     symmetry subgroup. */

#include <irrep/dmi.h>
#include <irrep/lattice.h>
#include <irrep/point_group.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void print_J_matrix(const char *label, int b, const double *J) {
    printf("      %s_%d = [[%+0.4f %+0.4f %+0.4f]\n", label, b + 1, J[0], J[1], J[2]);
    printf("                 [%+0.4f %+0.4f %+0.4f]\n", J[3], J[4], J[5]);
    printf("                 [%+0.4f %+0.4f %+0.4f]]\n", J[6], J[7], J[8]);
}

int main(void) {
    printf("=== libirrep — kagome NN bond exchange-tensor decomposition ===\n");
    printf("    Site symmetry under D_6 (in-plane).  Relevant to:\n");
    printf("      Fe₃Sn₂          (T_C ≈ 660 K — RT)\n");
    printf("      Mn₃Sn / Mn₃Ge   (T_N ≈ 420 / 380 K — RT)\n");
    printf("      Co₃Sn₂S₂        (T_C = 177 K, magnetic Weyl)\n\n");

    /* Kagome 2D lattice, pick NN bond between sublattice 0 (=(0,0)) and
     * sublattice 1 (=(1,0)) of the same cell. Bond direction along x. */
    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
    double           xy0[2], xy1[2];
    irrep_lattice_site_position(L, 0, xy0); /* sub 0 in cell (0,0) */
    irrep_lattice_site_position(L, 1, xy1); /* sub 1 in cell (0,0) */

    /* Pad positions to 3D (z = 0); centre the bond at the origin so the
     * libirrep point groups (operating around (0,0,0)) realise the bond's
     * site stabiliser in D_6. */
    double mid[3] = {0.5 * (xy0[0] + xy1[0]), 0.5 * (xy0[1] + xy1[1]), 0.0};
    double r_a[3] = {xy0[0] - mid[0], xy0[1] - mid[1], 0.0};
    double r_b[3] = {xy1[0] - mid[0], xy1[1] - mid[1], 0.0};
    double bondvec[3] = {xy1[0] - xy0[0], xy1[1] - xy0[1], 0.0};
    double bondlen = sqrt(bondvec[0] * bondvec[0] + bondvec[1] * bondvec[1]);
    double bondhat[3] = {bondvec[0] / bondlen, bondvec[1] / bondlen, 0.0};

    printf("  Kagome NN bond (sub 0 — sub 1):\n");
    printf("    r_a (cell origin)  = (%+0.4f, %+0.4f, %+0.4f)\n", xy0[0], xy0[1], 0.0);
    printf("    r_b                = (%+0.4f, %+0.4f, %+0.4f)\n", xy1[0], xy1[1], 0.0);
    printf("    bond direction     = (%+0.4f, %+0.4f, %+0.4f)   |bond| = %.4f\n", bondhat[0],
           bondhat[1], bondhat[2], bondlen);
    printf("\n");

    irrep_pg_table_t *d6 = irrep_pg_table_build(IRREP_PG_D6);
    irrep_pg_table_t *c3v = irrep_pg_table_build(IRREP_PG_C3V);
    irrep_pg_table_t *d3 = irrep_pg_table_build(IRREP_PG_D3);

    /* ---- Antisymmetric (DMI) part ---- */
    printf("  ━ DMI vector D under each candidate site symmetry ━\n\n");

    double basis_d6[9] = {0}, basis_c3v[9] = {0}, basis_d3[9] = {0};
    int    n_d6 = irrep_dmi_allowed_basis_from_pg(r_a, r_b, d6, 1e-6, basis_d6);
    int    n_c3v = irrep_dmi_allowed_basis_from_pg(r_a, r_b, c3v, 1e-6, basis_c3v);
    int    n_d3 = irrep_dmi_allowed_basis_from_pg(r_a, r_b, d3, 1e-6, basis_d3);

    printf("    D_6  (chiral hex, vertical mirrors): %d-dim", n_d6);
    if (n_d6 > 0) {
        printf("    D_1 = (%+0.4f, %+0.4f, %+0.4f)", basis_d6[0], basis_d6[1], basis_d6[2]);
    }
    printf("\n");

    printf("    C_3v (3-fold + σ_v's):             %d-dim", n_c3v);
    if (n_c3v > 0)
        printf("    D_1 = (%+0.4f, %+0.4f, %+0.4f)", basis_c3v[0], basis_c3v[1], basis_c3v[2]);
    printf("\n");

    printf("    D_3  (3-fold proper-only):         %d-dim", n_d3);
    if (n_d3 > 0)
        printf("    D_1 = (%+0.4f, %+0.4f, %+0.4f)", basis_d3[0], basis_d3[1], basis_d3[2]);
    printf("\n\n");

    /* ---- Symmetric (Heisenberg + Γ-type) part ---- */
    printf("  ━ Symmetric exchange tensor J^s under each candidate site symmetry ━\n\n");

    double J_basis[54] = {0};

    int n_J_d6 = irrep_exchange_symmetric_basis_from_pg(r_a, r_b, d6, 1e-6, J_basis);
    printf("    D_6 J^s: %d-dim\n", n_J_d6);
    for (int b = 0; b < n_J_d6; ++b)
        print_J_matrix("J", b, J_basis + b * 9);

    int n_J_d3 = irrep_exchange_symmetric_basis_from_pg(r_a, r_b, d3, 1e-6, J_basis);
    printf("\n    D_3 J^s: %d-dim\n", n_J_d3);
    for (int b = 0; b < n_J_d3; ++b)
        print_J_matrix("J", b, J_basis + b * 9);

    /* ---- Summary ---- */
    printf("\n  ━ Summary table ━\n");
    printf("  %-6s  %-8s  %-8s  %s\n", "sym", "DMI dim", "J^s dim", "physical content");
    printf("  %-6s  %-8d  %-8d  %s\n", "D_6", n_d6, n_J_d6,
           n_d6 == 1 ? "D ⊥ kagome plane (D_z component)" : "");
    printf("  %-6s  %-8d  %-8d  %s\n", "C_3v", n_c3v, 0 /* not computed but predictable */,
           "polar 3-fold; allows in-plane D");
    printf("  %-6s  %-8d  %-8d  %s\n", "D_3", n_d3, n_J_d3,
           n_d3 == 1 ? "no improper to flip parity (chiral)" : "");
    printf("\n");

    /* Adversarial interpretation: for layered kagome metallic ferromagnets
     * (Fe₃Sn₂, Co₃Sn₂S₂), the relevant 3D sym is D_6h (with σ_h). The σ_h
     * forbids the in-plane D vector (Moriya rule c — mirror containing
     * bond → D ⊥ that plane → D ∥ z). So 3D-D_6h NN bond DMI is
     * Z-only — and that's exactly what's measured in those compounds.
     *
     * For STACKING-BROKEN kagome (e.g., Fe₃Sn₂ has bilayer stacking that
     * breaks σ_h between non-equivalent layers), the inversion is also
     * broken. The site symmetry drops to D_6 (or lower); in-plane DMI
     * components become allowed, and the system enters the helimagnet /
     * skyrmion-bubble regime — the actual mechanism behind RT skyrmion-
     * like textures in Fe₃Sn₂. The toolkit predicts when this occurs
     * directly from the broken-symmetry point group. */
    printf("  Adversarial interpretation:\n");
    printf("    - D_6 alone (in-plane symmetry, no σ_h): DMI has 1 component along z.\n");
    printf("      This is the in-plane bond's out-of-plane DMI.\n");
    printf("    - For 3D D_6h (with σ_h between layers, applies to a non-broken-stacking\n");
    printf("      crystal): the σ_h would force D ∥ z too (mirror containing bond).\n");
    printf("    - When stacking breaks σ_h (Fe₃Sn₂ bilayers): DMI gains in-plane components,\n");
    printf("      driving the helimagnetic / RT skyrmion-bubble texture observed there.\n");
    printf("    - The materials hunt for RT chiral kagome magnets distils to: find compounds\n");
    printf("      that break σ_h and inversion through the kagome plane while preserving the\n");
    printf("      in-plane D_6 (or break it too, via Jahn-Teller distortion).\n");

    irrep_pg_table_free(d6);
    irrep_pg_table_free(c3v);
    irrep_pg_table_free(d3);
    irrep_lattice_free(L);
    return 0;
}
