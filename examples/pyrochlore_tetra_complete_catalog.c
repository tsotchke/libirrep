/* SPDX-License-Identifier: MIT */
/* Complete exchange-tensor catalog for the pyrochlore "up tetrahedron".
 *
 * The pyrochlore conventional cell holds 4 corner-sharing tetrahedra; the
 * "up tetrahedron" centred on FCC sublattice 0 has vertices at
 *   (0,0,0),  (¼,¼,0),  (¼,0,¼),  (0,¼,¼)
 * which are libirrep's pyrochlore sublattices 0, 1, 2, 3. The tetrahedron
 * has 6 NN bonds (all 6 edges) and 4 triangular faces. This example
 * tabulates the COMPLETE symmetry-allowed bilinear + trilinear coupling
 * structure under each candidate site symmetry:
 *
 *    For each bond:    DMI dim    + J^s dim    (D = 0, 1, 2, 3 components;
 *                                                J^s 1..6 components)
 *    For each face:    χ allowed (0/1)
 *
 * This is the parameter scaffolding a downstream DFT or micromagnetic
 * code needs to model the pyrochlore family — from group theory alone,
 * with no input from the materials literature. The Curnoe-Ross-Kao
 * parametrisation of Yb₂Ti₂O₇ / Er₂Ti₂O₇ / Tb₂Ti₂O₇ phase diagrams
 * (Ross 2011, PRX 1, 021002) is the form the J^s output reproduces.
 *
 * Build / run:
 *   make examples
 *   ./build/bin/pyrochlore_tetra_complete_catalog */

#include <irrep/dmi.h>
#include <irrep/lattice3d.h>
#include <irrep/point_group.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const int TET_VERTICES[4] = {0, 1, 2, 3};
static const int TET_BONDS[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
static const int TET_FACES[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};

static void centre_bond(const double r0[3], const double r1[3], double r_a_out[3],
                        double r_b_out[3]) {
    double mid[3] = {0.5 * (r0[0] + r1[0]), 0.5 * (r0[1] + r1[1]), 0.5 * (r0[2] + r1[2])};
    for (int k = 0; k < 3; ++k) {
        r_a_out[k] = r0[k] - mid[k];
        r_b_out[k] = r1[k] - mid[k];
    }
}

static void centre_triangle(const double r0[3], const double r1[3], const double r2[3],
                            double r_a_out[3], double r_b_out[3], double r_c_out[3]) {
    double cent[3] = {(r0[0] + r1[0] + r2[0]) / 3.0, (r0[1] + r1[1] + r2[1]) / 3.0,
                      (r0[2] + r1[2] + r2[2]) / 3.0};
    for (int k = 0; k < 3; ++k) {
        r_a_out[k] = r0[k] - cent[k];
        r_b_out[k] = r1[k] - cent[k];
        r_c_out[k] = r2[k] - cent[k];
    }
}

static void analyze_under_pg(const char *pg_name, irrep_pg_table_t *pg,
                             const double pos[4][3]) {
    printf("\n  ━ Site symmetry: %s ━\n", pg_name);

    /* DMI + J^s per bond. */
    printf("    Bonds (6 NN edges):\n");
    printf("    %-9s  %-9s  %-9s  %-s\n", "endpoints", "DMI dim", "J^s dim", "interpretation");
    for (int b = 0; b < 6; ++b) {
        double r_a[3], r_b[3];
        centre_bond(pos[TET_BONDS[b][0]], pos[TET_BONDS[b][1]], r_a, r_b);

        double dmi_basis[9] = {0};
        int    n_dmi = irrep_dmi_allowed_basis_from_pg(r_a, r_b, pg, 1e-6, dmi_basis);

        double jsym_basis[54] = {0};
        int    n_jsym = irrep_exchange_symmetric_basis_from_pg(r_a, r_b, pg, 1e-6, jsym_basis);

        const char *interp = "";
        if (n_dmi == 0 && n_jsym == 1)
            interp = "Heisenberg only (full inversion-stabilised)";
        else if (n_dmi == 0)
            interp = "no DMI; nontrivial J^s anisotropy";
        else if (n_dmi == 1)
            interp = "DMI 1-dim + J^s anisotropy";
        else
            interp = "high DMI freedom (low-symmetry bond)";

        printf("    %d → %d      %d-dim      %d-dim      %s\n", TET_BONDS[b][0], TET_BONDS[b][1],
               n_dmi, n_jsym, interp);
    }

    /* Chirality per triangle face. */
    printf("    Triangle faces (4):\n");
    printf("    %-13s  %-s\n", "vertices", "scalar chirality χ");
    for (int f = 0; f < 4; ++f) {
        double r_a[3], r_b[3], r_c[3];
        centre_triangle(pos[TET_FACES[f][0]], pos[TET_FACES[f][1]], pos[TET_FACES[f][2]], r_a,
                        r_b, r_c);
        int rc = irrep_chirality_allowed_from_pg(r_a, r_b, r_c, pg, 1e-6);
        printf("    {%d, %d, %d}      %s\n", TET_FACES[f][0], TET_FACES[f][1], TET_FACES[f][2],
               rc == 1 ? "ALLOWED"
                       : rc == 0 ? "FORBIDDEN" : "(no preserving op)");
    }
}

int main(void) {
    printf("=== libirrep — pyrochlore up-tetrahedron complete exchange catalog ===\n");
    printf("    All bond DMI + J^s + face chirality, derived from group theory alone.\n");

    /* Pyrochlore 1×1×1 — pick the "up tetrahedron" at FCC sublattice 0. */
    irrep_lattice3d_t *L = irrep_lattice3d_build(IRREP_LATTICE3D_PYROCHLORE, 1, 1, 1);

    double pos[4][3];
    for (int v = 0; v < 4; ++v)
        irrep_lattice3d_site_position(L, TET_VERTICES[v], pos[v]);

    /* Centre the tetrahedron at the origin so the libirrep point groups
     * (operating around (0,0,0)) realise the bond / triangle stabilisers. */
    double cent[3] = {0, 0, 0};
    for (int v = 0; v < 4; ++v)
        for (int k = 0; k < 3; ++k)
            cent[k] += pos[v][k] / 4.0;
    for (int v = 0; v < 4; ++v)
        for (int k = 0; k < 3; ++k)
            pos[v][k] -= cent[k];

    printf("\n  Tetrahedron vertices (centred at origin):\n");
    for (int v = 0; v < 4; ++v)
        printf("    sub %d:  (%+0.4f, %+0.4f, %+0.4f)\n", TET_VERTICES[v], pos[v][0],
               pos[v][1], pos[v][2]);

    irrep_pg_table_t *Oh = irrep_pg_table_build(IRREP_PG_OH);
    irrep_pg_table_t *O = irrep_pg_table_build(IRREP_PG_O);
    irrep_pg_table_t *Td = irrep_pg_table_build(IRREP_PG_TD);

    analyze_under_pg("O_h (full cubic, with inversion)", Oh, pos);
    analyze_under_pg("O   (chiral cubic, no inversion)", O, pos);
    analyze_under_pg("T_d (tetrahedral with σ_d)", Td, pos);

    printf("\n  ━ Interpretation ━\n");
    printf("\n  Pyrochlore site symmetry varies by sublattice and is famously\n");
    printf("  D_3d at the magnetic 16d Wyckoff position; the ABOVE analysis\n");
    printf("  uses cubic point groups to assess what the parent O_h imposes\n");
    printf("  globally rather than the per-site D_3d. The spinel structure\n");
    printf("  Fd-3m has full O_h with inversion centres at the 4-coordinate\n");
    printf("  (1/2, 1/2, 1/2)-equivalent positions; bonds whose midpoints sit\n");
    printf("  on inversion centres carry zero DMI by Moriya rule (a). Other\n");
    printf("  bonds (where the midpoint sits between two NN positions but\n");
    printf("  NOT on an inversion centre) acquire non-zero DMI direction —\n");
    printf("  the Elhajal-Canals-Lacroix DMI pattern (PRB 2002).\n");
    printf("\n  For materials-design: this catalog is the parameter scaffold\n");
    printf("  that DFT magnitudes plug into. Every DMI / J^s / χ entry\n");
    printf("  marked ALLOWED with `n > 0` is a free parameter that DFT must\n");
    printf("  determine from first principles; entries marked 0 / FORBIDDEN\n");
    printf("  are zeroed by symmetry and should NOT be Wannier-fitted (any\n");
    printf("  nonzero DFT result on them is numerical noise).\n");

    irrep_pg_table_free(Oh);
    irrep_pg_table_free(O);
    irrep_pg_table_free(Td);
    irrep_lattice3d_free(L);
    return 0;
}
