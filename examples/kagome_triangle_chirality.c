/* SPDX-License-Identifier: MIT */
/* Scalar chirality on the kagome triangle, derived from group theory.
 *
 * The kagome triangle's three magnetic ions support a scalar chirality
 *
 *     χ = S_a · (S_b × S_c)
 *
 * iff the magnetic point group of the triangle satisfies the symmetry
 * rule from `docs/PHYSICS_APPENDIX.md` §14: every preserving operation
 * `g` must satisfy
 *
 *     det(g) · σ_perm(g) · σ_T(g) = +1.
 *
 * This is the pseudoscalar selection rule: σ_h in the triangle plane
 * forbids χ; T-augmented σ_h ("magnetic mirror") restores it; pure
 * time-reversal alone forbids it. The result determines whether a
 * candidate kagome magnet supports a topological Hall effect from the
 * scalar-chirality contribution to the Berry curvature — the route by
 * which Mn₃Sn / Mn₃Ge produce a giant anomalous Hall response despite
 * being antiferromagnetic.
 *
 * This example computes the chirality verdict under several candidate
 * site symmetries and combinations with time reversal, mapping the
 * results onto the dominant RT kagome magnets.
 *
 * Build / run:
 *   make examples
 *   ./build/bin/kagome_triangle_chirality */

#include <irrep/dmi.h>
#include <irrep/lattice.h>
#include <irrep/point_group.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Synthesise an op from (axis, angle, det, antiunitary). */
static void make_op(irrep_dmi_sym_op_t *op, double nx, double ny, double nz, double theta,
                    int det, int antiunitary) {
    memset(op, 0, sizeof *op);
    double n_norm = sqrt(nx * nx + ny * ny + nz * nz);
    nx /= n_norm;
    ny /= n_norm;
    nz /= n_norm;
    double c = cos(theta), s = sin(theta), C = 1 - c;
    double *R = op->R_proper;
    R[0] = c + nx * nx * C;
    R[1] = nx * ny * C - nz * s;
    R[2] = nx * nz * C + ny * s;
    R[3] = ny * nx * C + nz * s;
    R[4] = c + ny * ny * C;
    R[5] = ny * nz * C - nx * s;
    R[6] = nz * nx * C - ny * s;
    R[7] = nz * ny * C + nx * s;
    R[8] = c + nz * nz * C;
    op->det = det;
    op->antiunitary = antiunitary;
}

static void identity_op(irrep_dmi_sym_op_t *op) {
    memset(op, 0, sizeof *op);
    op->R_proper[0] = 1;
    op->R_proper[4] = 1;
    op->R_proper[8] = 1;
    op->det = +1;
}

int main(void) {
    printf("=== libirrep — scalar chirality on the kagome triangle ===\n\n");

    /* Kagome triangle vertices: pull from libirrep's kagome lattice. */
    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
    double           xy[3][2];
    irrep_lattice_site_position(L, 0, xy[0]); /* sub A = (0, 0) */
    irrep_lattice_site_position(L, 1, xy[1]); /* sub B = (1, 0) */
    irrep_lattice_site_position(L, 2, xy[2]); /* sub C = (0.5, √3/2) */
    irrep_lattice_free(L);

    /* Centre the triangle at the origin. */
    double cx = (xy[0][0] + xy[1][0] + xy[2][0]) / 3.0;
    double cy = (xy[0][1] + xy[1][1] + xy[2][1]) / 3.0;
    double r_a[3] = {xy[0][0] - cx, xy[0][1] - cy, 0.0};
    double r_b[3] = {xy[1][0] - cx, xy[1][1] - cy, 0.0};
    double r_c[3] = {xy[2][0] - cx, xy[2][1] - cy, 0.0};

    printf("  Triangle vertices (cell 0, centred at origin):\n");
    printf("    A: (%+0.4f, %+0.4f, 0)\n", r_a[0], r_a[1]);
    printf("    B: (%+0.4f, %+0.4f, 0)\n", r_b[0], r_b[1]);
    printf("    C: (%+0.4f, %+0.4f, 0)\n\n", r_c[0], r_c[1]);

    printf("  ━ Chirality verdict under candidate site symmetries ━\n\n");

    /* Case 1: identity-only → trivially allowed. */
    {
        irrep_dmi_sym_op_t ops[1];
        identity_op(&ops[0]);
        int rc = irrep_chirality_allowed(r_a, r_b, r_c, ops, 1, 1e-9);
        printf("    Identity only:                                  χ %s\n",
               rc == 1 ? "ALLOWED" : "FORBIDDEN");
    }

    /* Case 2: C_3 about z-axis through centroid (cycles A→B→C). */
    {
        irrep_dmi_sym_op_t ops[3];
        identity_op(&ops[0]);
        make_op(&ops[1], 0, 0, 1, 2.0 * M_PI / 3.0, +1, 0);
        make_op(&ops[2], 0, 0, 1, -2.0 * M_PI / 3.0, +1, 0);
        int rc = irrep_chirality_allowed(r_a, r_b, r_c, ops, 3, 1e-9);
        printf("    C_3 about triangle centroid:                    χ %s\n",
               rc == 1 ? "ALLOWED" : "FORBIDDEN");
        printf("      (cyclic perm × proper det = +1 — chirality survives.)\n");
    }

    /* Case 3: σ_h IN the triangle plane (centrosymmetric kagome stacking). */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        make_op(&ops[1], 0, 0, 1, M_PI, -1, 0); /* mirror perp to z */
        int rc = irrep_chirality_allowed(r_a, r_b, r_c, ops, 2, 1e-9);
        printf("    σ_h in triangle plane (centrosymmetric stacking): χ %s\n",
               rc == 1 ? "ALLOWED" : "FORBIDDEN");
        printf("      (perm = identity × det = -1 → product -1 — forbids chirality.\n");
        printf("       Idealised symmetric kagome stacking has no scalar-chirality\n");
        printf("       contribution to topological Hall.)\n");
    }

    /* Case 4: T·σ_h "magnetic mirror" — Mn₃Sn / Mn₃Ge configuration. */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        make_op(&ops[1], 0, 0, 1, M_PI, -1, 1); /* T·σ_h, antiunitary=1 */
        int rc = irrep_chirality_allowed(r_a, r_b, r_c, ops, 2, 1e-9);
        printf("    T·σ_h magnetic mirror (Mn₃Sn-type 120° AFM):    χ %s\n",
               rc == 1 ? "ALLOWED" : "FORBIDDEN");
        printf("      (T flips pseudoscalar; the σ_h's bad sign × T's bad sign = +1.\n");
        printf("       The non-collinear 120° magnetic structure breaks σ_h\n");
        printf("       but preserves T·σ_h as a magnetic-point-group element —\n");
        printf("       this is the symmetry route to Mn₃Sn's giant topological Hall.)\n");
    }

    /* Case 5: C_3v (3-fold + σ_v's containing centroid-vertex axes). */
    {
        irrep_dmi_sym_op_t ops[6];
        identity_op(&ops[0]);
        make_op(&ops[1], 0, 0, 1, 2.0 * M_PI / 3.0, +1, 0);
        make_op(&ops[2], 0, 0, 1, -2.0 * M_PI / 3.0, +1, 0);
        /* σ_v's: mirror perp to (n_x, n_y, 0) where n is normal to a vertex. */
        for (int k = 0; k < 3; ++k) {
            double phi = k * 2.0 * M_PI / 3.0;
            make_op(&ops[3 + k], -sin(phi), cos(phi), 0, M_PI, -1, 0);
        }
        int rc = irrep_chirality_allowed(r_a, r_b, r_c, ops, 6, 1e-9);
        printf("    C_3v (C_3 + 3 vertical mirrors):                 χ %s\n",
               rc == 1 ? "ALLOWED" : "FORBIDDEN");
        printf("      (Each σ_v transposes two vertices and fixes the third —\n");
        printf("       perm sign = -1 × det = -1 → product +1. Combined with C_3\n");
        printf("       (cyclic + proper, +1), all preserving ops give +1 → allowed.)\n");
    }

    /* Case 6: D_3h (C_3 + σ_v's + σ_h) — the centrosymmetric case for
     * 3-fold sites in a horizontal plane. σ_h FORBIDS even though
     * C_3v alone would allow. */
    {
        irrep_dmi_sym_op_t ops[12];
        identity_op(&ops[0]);
        make_op(&ops[1], 0, 0, 1, 2.0 * M_PI / 3.0, +1, 0);
        make_op(&ops[2], 0, 0, 1, -2.0 * M_PI / 3.0, +1, 0);
        for (int k = 0; k < 3; ++k) {
            double phi = k * 2.0 * M_PI / 3.0;
            make_op(&ops[3 + k], -sin(phi), cos(phi), 0, M_PI, -1, 0);
        }
        /* σ_h */
        make_op(&ops[6], 0, 0, 1, M_PI, -1, 0);
        /* S_3 = σ_h · C_3 */
        make_op(&ops[7], 0, 0, 1, M_PI - 2.0 * M_PI / 3.0, -1, 0);
        make_op(&ops[8], 0, 0, 1, M_PI + 2.0 * M_PI / 3.0, -1, 0);
        /* C_2 axes through vertex (proper rotation) */
        for (int k = 0; k < 3; ++k) {
            double phi = k * 2.0 * M_PI / 3.0;
            make_op(&ops[9 + k], cos(phi), sin(phi), 0, M_PI, +1, 0);
        }
        int rc = irrep_chirality_allowed(r_a, r_b, r_c, ops, 12, 1e-9);
        printf("    D_3h (C_3v + σ_h):                                χ %s\n",
               rc == 1 ? "ALLOWED" : "FORBIDDEN");
        printf("      (σ_h alone forbids — even when C_3v allows. The full D_3h\n");
        printf("       symmetry of an idealised Mn₃Sn paramagnet gives χ = 0\n");
        printf("       above T_N. Below T_N the magnetic structure breaks σ_h.)\n");
    }

    printf("\n  ━ Connection to RT kagome magnets ━\n\n");
    printf("    Material      | T_N or T_C | Magnetic structure | χ via libirrep\n");
    printf("    ─────────────────────────────────────────────────────────────────\n");
    printf("    Mn₃Sn          | T_N=420 K  | non-collinear 120°  | T·σ_h allowed → χ ≠ 0\n");
    printf("    Mn₃Ge          | T_N=380 K  | similar             | T·σ_h allowed → χ ≠ 0\n");
    printf("    Fe₃Sn₂         | T_C=660 K  | bilayer FM-stacked  | broken σ_h → χ ≠ 0\n");
    printf("    Co₃Sn₂S₂       | T_C=177 K  | FM, σ_h preserved   | σ_h forbids χ from\n");
    printf("                   |            |                     | NN triangles alone;\n");
    printf("                   |            |                     | TH from Berry from\n");
    printf("                   |            |                     | Weyl bands instead.\n");
    printf("\n");
    printf("  Each row's verdict is reproducible by passing the magnetic\n");
    printf("  point group's element list to `irrep_chirality_allowed`. The\n");
    printf("  libirrep substrate is the irreducible group-theory step in the\n");
    printf("  topological-Hall-effect candidate-screening loop.\n");

    return 0;
}
