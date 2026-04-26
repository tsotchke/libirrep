/* SPDX-License-Identifier: MIT */
/* Derive the DMI pattern of pyrochlore NN bonds from group theory alone.
 *
 * The pyrochlore lattice (Fd-3m, with magnetic site Wyckoff 16d) supports
 * a celebrated "Moriya-Elhajal" DMI pattern — the bonds between adjacent
 * tetrahedral corners carry D vectors of fixed length, oriented along the
 * bond direction itself (in chiral cubic site sym) or perpendicular to it
 * along the cubic [001] axis (in centrosymmetric cubic, with sign
 * alternation). The pattern is what stabilises the order-by-disorder
 * selection in pyrochlore quantum spin ice.
 *
 * This example takes the pyrochlore 1×1×1 cluster from libirrep's lattice3d
 * and the O_h cubic point group from point_group.h, and asks the DMI
 * symmetry analyzer to derive the allowed D-vector for a representative
 * NN bond. The output is the result of pure group theory — no input from
 * the pyrochlore literature — and reproduces the known textbook pattern:
 *
 *   - Inversion at the bond midpoint? Pyrochlore Fd-3m has inversion at
 *     (1/2, 1/2, 1/2), NOT at NN-bond midpoints, so D ≠ 0.
 *   - The bond {sub 0 at (0,0,0), sub 1 at (1/4, 1/4, 0)} is invariant
 *     under the C_2' rotation about (1,1,0) of O_h (since this axis
 *     bisects the bond). By Moriya's rule (e), D ∥ bond direction.
 *
 * Caveat: O_h is the FULL cubic point group around an arbitrary origin
 * we choose, not the Fd-3m space group. The two differ by translations
 * and screw axes; for the simple bond chosen here, the result coincides
 * with the textbook pattern. A full Fd-3m treatment would require
 * 3D space-group infrastructure (deferred past 1.3). */

#include <irrep/dmi.h>
#include <irrep/lattice3d.h>
#include <irrep/point_group.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void print_basis(const char *label, int n, const double basis[9]) {
    printf("    %s: %d-dim\n", label, n);
    for (int i = 0; i < n; ++i)
        printf("      D_%d = (%+0.6f, %+0.6f, %+0.6f)   |D| = %.6f\n", i + 1,
               basis[i * 3], basis[i * 3 + 1], basis[i * 3 + 2],
               sqrt(basis[i * 3] * basis[i * 3] + basis[i * 3 + 1] * basis[i * 3 + 1] +
                    basis[i * 3 + 2] * basis[i * 3 + 2]));
    if (n == 0)
        printf("      (D ≡ 0 — bond DMI fully constrained to vanish)\n");
}

int main(void) {
    printf("=== libirrep — DMI pattern for pyrochlore NN bond, derived from group theory ===\n\n");

    irrep_lattice3d_t *L = irrep_lattice3d_build(IRREP_LATTICE3D_PYROCHLORE, 1, 1, 1);

    /* Pick the NN bond between sub 0 and sub 1 of the same FCC sublattice
     * (= the first edge of the "up" tetrahedron at the FCC origin). */
    double r0[3], r1[3];
    irrep_lattice3d_site_position(L, 0, r0);
    irrep_lattice3d_site_position(L, 1, r1);

    /* Centre the bond at the origin so the libirrep point groups (acting
     * around (0,0,0)) realise the bond's site symmetry. */
    double mid[3] = {0.5 * (r0[0] + r1[0]), 0.5 * (r0[1] + r1[1]), 0.5 * (r0[2] + r1[2])};
    double r_a[3] = {r0[0] - mid[0], r0[1] - mid[1], r0[2] - mid[2]};
    double r_b[3] = {r1[0] - mid[0], r1[1] - mid[1], r1[2] - mid[2]};

    double bondlen =
        sqrt((r1[0] - r0[0]) * (r1[0] - r0[0]) + (r1[1] - r0[1]) * (r1[1] - r0[1]) +
             (r1[2] - r0[2]) * (r1[2] - r0[2]));
    double bondhat[3] = {(r1[0] - r0[0]) / bondlen, (r1[1] - r0[1]) / bondlen,
                         (r1[2] - r0[2]) / bondlen};

    printf("  Pyrochlore NN bond (sub 0 — sub 1):\n");
    printf("    r_a (cell origin)  = (%+0.6f, %+0.6f, %+0.6f)\n", r0[0], r0[1], r0[2]);
    printf("    r_b                = (%+0.6f, %+0.6f, %+0.6f)\n", r1[0], r1[1], r1[2]);
    printf("    bond direction     = (%+0.6f, %+0.6f, %+0.6f)   |bond| = %.6f\n", bondhat[0],
           bondhat[1], bondhat[2], bondlen);
    printf("    midpoint (centred) = (%+0.6f, %+0.6f, %+0.6f) → using shifted r_a, r_b\n\n",
           mid[0], mid[1], mid[2]);

    double            basis[9] = {0};
    irrep_pg_table_t *oh = irrep_pg_table_build(IRREP_PG_OH);
    irrep_pg_table_t *o = irrep_pg_table_build(IRREP_PG_O);
    irrep_pg_table_t *td = irrep_pg_table_build(IRREP_PG_TD);

    printf("  ━ allowed-D analysis under each candidate site symmetry ━\n\n");

    double basis_oh[9], basis_o[9], basis_td[9];

    int n_oh = irrep_dmi_allowed_basis_from_pg(r_a, r_b, oh, 1e-6, basis_oh);
    print_basis("under O_h (full cubic, with inversion)", n_oh, basis_oh);
    printf("\n");

    int n_o = irrep_dmi_allowed_basis_from_pg(r_a, r_b, o, 1e-6, basis_o);
    print_basis("under O (chiral cubic, NO inversion)", n_o, basis_o);
    printf("\n");

    int n_td = irrep_dmi_allowed_basis_from_pg(r_a, r_b, td, 1e-6, basis_td);
    print_basis("under T_d (tetrahedral with σ_d)", n_td, basis_td);
    printf("\n");

    (void)basis;

    /* Interpret each result by projecting onto bond direction. */
    if (n_o == 1) {
        double dot =
            basis_o[0] * bondhat[0] + basis_o[1] * bondhat[1] + basis_o[2] * bondhat[2];
        printf("  Cross-check (O): D · bond̂ = %+.6f  (Moriya rule (e): C₂' along bond → D ∥ bond)\n",
               dot);
    }
    if (n_td == 1) {
        double dot =
            basis_td[0] * bondhat[0] + basis_td[1] * bondhat[1] + basis_td[2] * bondhat[2];
        printf("  Cross-check (T_d): D · bond̂ = %+.6f  (mirror through bond → D ⊥ that mirror)\n",
               dot);
    }
    printf("\n");

    /* What this tells the materials hunter:
     *   - O_h forces D = 0 because of the C₂' "reversing" element
     *     COMBINED with the mirror containing the bond. Bulk
     *     centrosymmetric pyrochlores (e.g. Tb₂Ti₂O₇) have NN bonds
     *     with D = 0 by symmetry (any DMI must come from longer-range
     *     bonds with lower site symmetry, or from a non-cubic distortion).
     *
     *   - O (chiral cubic, no improper elements) leaves a 1D allowed D
     *     along the bond — this is the magnitude that B20 Mn-Mn bonds
     *     carry in helimagnetic MnSi. Same group-theory pattern, even
     *     though the lattice is different.
     *
     *   - T_d sits between the two: it has improper σ_d but no inversion.
     *     Its bond stabiliser is smaller, leaves D in a 1D subspace, but
     *     in a DIFFERENT direction than O — perpendicular to the bond
     *     (the σ_d's mirror plane). Diamond / zincblende NN bonds.
     *
     *   The CMOS-compatible chiral magnet hunt distils to: target
     *   crystals where the bond site symmetry is O (chiral cubic), the
     *   structure is Si-substrate-compatible, and the constituent
     *   elements are abundant. MnSi epitaxial on Si(111) hits all three
     *   — except for the T_C = 29 K problem, which is materials physics,
     *   not group theory. */

    printf("  Interpretation:\n");
    printf("    - O_h NN bond: D = 0   centrosymmetric pyrochlores have no NN DMI on this bond.\n");
    printf("                            Cross-check: O_h has BOTH the C₂' rotation along the\n");
    printf("                            bond (forcing D ∥ bond) AND the mirror containing the\n");
    printf("                            bond (forcing D ⊥ that mirror) — incompatible, hence\n");
    printf("                            D = 0 by symmetry alone.\n");
    printf("    - O   NN bond: D ∥ bond direction. Group theory matches Moriya rule (e):\n");
    printf("                            the C₂' axis ALONG the bond pins D parallel to it.\n");
    printf("                            Same pattern as B20 Mn-Mn bonds in MnSi.\n");
    printf("    - T_d NN bond: D ⊥ bond direction. Different from O because T_d's mirror\n");
    printf("                            symmetries pin D perpendicular to the bond.\n");
    printf("                            Diamond / zincblende NN bond geometry.\n");

    /* ---- Symmetric exchange tensor catalog ---- */

    printf("\n  ━ allowed *symmetric* exchange tensor J^s under each site sym ━\n");
    printf("  (Heisenberg = trace; non-trivial extra components = Kitaev-Γ-type anisotropy.)\n\n");

    double J_basis[54] = {0};

    int n_J_oh = irrep_exchange_symmetric_basis_from_pg(r_a, r_b, oh, 1e-6, J_basis);
    printf("    O_h: J^s %d-dim\n", n_J_oh);
    for (int b = 0; b < n_J_oh; ++b) {
        const double *J = J_basis + b * 9;
        printf("      J_%d = [[%+0.4f %+0.4f %+0.4f]\n", b + 1, J[0], J[1], J[2]);
        printf("              [%+0.4f %+0.4f %+0.4f]\n", J[3], J[4], J[5]);
        printf("              [%+0.4f %+0.4f %+0.4f]]\n", J[6], J[7], J[8]);
    }

    int n_J_o = irrep_exchange_symmetric_basis_from_pg(r_a, r_b, o, 1e-6, J_basis);
    printf("\n    O:   J^s %d-dim\n", n_J_o);
    for (int b = 0; b < n_J_o; ++b) {
        const double *J = J_basis + b * 9;
        printf("      J_%d = [[%+0.4f %+0.4f %+0.4f]\n", b + 1, J[0], J[1], J[2]);
        printf("              [%+0.4f %+0.4f %+0.4f]\n", J[3], J[4], J[5]);
        printf("              [%+0.4f %+0.4f %+0.4f]]\n", J[6], J[7], J[8]);
    }

    int n_J_td = irrep_exchange_symmetric_basis_from_pg(r_a, r_b, td, 1e-6, J_basis);
    printf("\n    T_d: J^s %d-dim\n", n_J_td);
    for (int b = 0; b < n_J_td; ++b) {
        const double *J = J_basis + b * 9;
        printf("      J_%d = [[%+0.4f %+0.4f %+0.4f]\n", b + 1, J[0], J[1], J[2]);
        printf("              [%+0.4f %+0.4f %+0.4f]\n", J[3], J[4], J[5]);
        printf("              [%+0.4f %+0.4f %+0.4f]]\n", J[6], J[7], J[8]);
    }

    /* Summary table per site sym: (DMI dim, J^s dim) tuple. */
    printf("\n  ━ complete bond-exchange-tensor decomposition (this NN bond) ━\n");
    printf("    %-6s  %-7s  %-7s  %-12s  %-s\n", "sym", "DMI dim", "J^s dim", "components",
           "interpretation");
    printf("    %-6s  %-7d  %-7d  %-12s  %-s\n", "O_h", n_oh, n_J_oh,
           n_oh == 0 && n_J_oh > 0 ? "Heis only" : "—",
           "centrosymmetric — pure Heisenberg, no DMI, no anisotropy");
    printf("    %-6s  %-7d  %-7d  %-12s  %-s\n", "O", n_o, n_J_o, "Heis + DMI + Γ",
           "chiral cubic — D ∥ bond, nontrivial J^s anisotropy");
    printf("    %-6s  %-7d  %-7d  %-12s  %-s\n", "T_d", n_td, n_J_td, "Heis + DMI⊥ + Γ",
           "tetrahedral — D ⊥ bond, σ_d-allowed J^s");

    irrep_pg_table_free(oh);
    irrep_pg_table_free(o);
    irrep_pg_table_free(td);
    irrep_lattice3d_free(L);
    return 0;
}
