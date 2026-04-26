/* SPDX-License-Identifier: MIT */
/* Crystal-field splitting of d-orbitals under cubic point groups.
 *
 * Demonstrates the O_h and T_d projectors on the l=2 (parity-even)
 * spherical-harmonic basis — the canonical demonstration of crystal-field
 * theory on transition-metal d-orbitals. The expected splittings:
 *
 *   O_h:  d  →  e_g (2)  +  t_2g (3)
 *   T_d:  d  →  e   (2)  +  t_2  (3)
 *   O:    d  →  e   (2)  +  t_2  (3)   (same numerical splitting as T_d
 *                                        — O is the proper subgroup)
 *
 * The implied energy split between e_g and t_2g (Δ_oct = 10Dq) is the
 * organising parameter of transition-metal complex spectra. libirrep
 * provides the *symmetry decomposition*; the energy gap itself is set
 * by the radial part of the crystal field, not by group theory.
 *
 * Build / run:
 *   make examples
 *   ./build/bin/cubic_crystal_field */

#include <irrep/multiset.h>
#include <irrep/point_group.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double frob_norm_sq(const double *v, int n) {
    double s = 0;
    for (int i = 0; i < n; ++i)
        s += v[i] * v[i];
    return s;
}

static void print_one(const irrep_pg_table_t *t, const irrep_multiset_t *spec, const double *in,
                      int dim) {
    double *out = calloc((size_t)dim, sizeof(double));
    int     n_irreps = irrep_pg_num_irreps(t);
    double  total = 0;
    printf("    %-6s ", "irrep");
    for (int mu = 0; mu < n_irreps; ++mu)
        printf(" %-6s", irrep_pg_irrep_label(t, mu));
    printf("\n    %-6s ", "‖P_μ‖²");
    for (int mu = 0; mu < n_irreps; ++mu) {
        memset(out, 0, (size_t)dim * sizeof(double));
        irrep_pg_project(t, mu, spec, in, out);
        double w = frob_norm_sq(out, dim);
        total += w;
        printf(" %-6.3f", w);
    }
    printf("\n    Σ ‖P_μ‖² = %.6f  (input ‖·‖² = %.6f → projector completeness)\n", total,
           frob_norm_sq(in, dim));
    free(out);
}

static void demo_d_orbital(irrep_point_group_t group, const char *group_label) {
    printf("\n━ %s — d-orbital (l=2 even) splitting ━\n", group_label);

    irrep_pg_table_t *t = irrep_pg_table_build(group);
    irrep_multiset_t *spec = irrep_multiset_parse("1x2e");
    int               dim = 5; /* 2l+1 = 5 d-orbital components */

    /* Reduction multiplicities: which irreps appear in 1x2e ? */
    int *mult = calloc((size_t)irrep_pg_num_irreps(t), sizeof(int));
    irrep_pg_reduce(t, spec, mult);
    printf("    decomposition (multiplicities):");
    for (int mu = 0; mu < irrep_pg_num_irreps(t); ++mu)
        if (mult[mu] > 0)
            printf("  %dx%s", mult[mu], irrep_pg_irrep_label(t, mu));
    printf("\n");
    free(mult);

    /* Probe basis: each of the 5 components in turn. The l=2 real-SH order
     * within libirrep is the standard m = -2, -1, 0, +1, +2 mapping which,
     * in cartesian form, corresponds to (xy, yz, 3z²-r², xz, x²-y²). The
     * crystal-field literature labels these as (d_xy, d_yz, d_z², d_xz, d_x²-y²).
     * Under O_h: d_z² and d_x²-y² → e_g; d_xy, d_yz, d_xz → t_2g. */
    const char *basis_names[5] = {"d_xy", "d_yz", "d_z²", "d_xz", "d_x²-y²"};
    for (int b = 0; b < dim; ++b) {
        double in[5] = {0};
        in[b] = 1.0;
        printf("  basis component %d (%s):\n", b, basis_names[b]);
        print_one(t, spec, in, dim);
    }

    irrep_multiset_free(spec);
    irrep_pg_table_free(t);
}

static void demo_p_orbital(irrep_point_group_t group, const char *group_label) {
    printf("\n━ %s — p-orbital (l=1 odd) splitting ━\n", group_label);

    irrep_pg_table_t *t = irrep_pg_table_build(group);
    irrep_multiset_t *spec = irrep_multiset_parse("1x1o");
    int               dim = 3;

    int *mult = calloc((size_t)irrep_pg_num_irreps(t), sizeof(int));
    irrep_pg_reduce(t, spec, mult);
    printf("    decomposition (multiplicities):");
    for (int mu = 0; mu < irrep_pg_num_irreps(t); ++mu)
        if (mult[mu] > 0)
            printf("  %dx%s", mult[mu], irrep_pg_irrep_label(t, mu));
    printf("\n");
    free(mult);

    /* In O_h the polar vector p = (x, y, z) is a single 3-dim irrep T_1u —
     * the three p-orbitals do NOT split. In T_d (no inversion), the same
     * polar vector lives in T_2 — likewise no splitting. */
    const char *basis_names[3] = {"p_y", "p_z", "p_x"};
    for (int b = 0; b < dim; ++b) {
        double in[3] = {0};
        in[b] = 1.0;
        printf("  basis component %d (%s):\n", b, basis_names[b]);
        print_one(t, spec, in, dim);
    }

    irrep_multiset_free(spec);
    irrep_pg_table_free(t);
}

int main(void) {
    printf("=== libirrep — cubic crystal-field splittings ===\n");
    printf("    Each |P_μ ψ|² is the weight of |ψ⟩ in irrep μ.\n");
    printf("    Σ_μ |P_μ ψ|² = |ψ|² is projector completeness.\n");

    demo_p_orbital(IRREP_PG_OH, "O_h");
    demo_d_orbital(IRREP_PG_OH, "O_h");
    demo_d_orbital(IRREP_PG_TD, "T_d");

    return 0;
}
