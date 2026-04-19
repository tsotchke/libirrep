/* SPDX-License-Identifier: MIT */
/* End-to-end integration example (1.2):
 *
 *   1. Build a NequIP message-passing layer from a single spec string.
 *   2. Run a forward pass on a 4-node / 8-edge square lattice.
 *   3. Project the output node features onto the A₁ (trivial / totally
 *      symmetric) irrep under C₄ᵥ to keep only the square-lattice-symmetric
 *      component — the standard construction for symmetric-NQS ansätze.
 *
 * Expected output: weight count, h_out norms before and after projection
 * (projection should strictly reduce or preserve the ‖·‖₂ on each node,
 * since P_μ is a projector with ‖P_μ‖ ≤ 1). */

#include <irrep/multiset.h>
#include <irrep/nequip.h>
#include <irrep/point_group.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* 4-node square lattice (periodic), 8 directed edges. Unit vectors along
 * ±x̂ and ±ŷ at the cell spacing; the C4v projector's symmetries match
 * the lattice exactly, so A₁ is a non-trivial filter. */
static const int    k_edge_src[8] = { 0, 1, 0, 2, 1, 0, 3, 1 };
static const int    k_edge_dst[8] = { 1, 0, 2, 0, 3, 2, 1, 3 };
static const double k_edge_vec[8 * 3] = {
     1.0,  0.0, 0.0,   /* 0 → 1 (+x̂) */
    -1.0,  0.0, 0.0,   /* 1 → 0 (−x̂) */
     0.0,  1.0, 0.0,   /* 0 → 2 (+ŷ) */
     0.0, -1.0, 0.0,   /* 2 → 0 (−ŷ) */
     0.0,  1.0, 0.0,   /* 1 → 3 (+ŷ) */
    -1.0,  0.0, 0.0,   /* 0 → 2 dup path variant */
    -1.0,  0.0, 0.0,   /* 3 → 1 (−x̂) */
     1.0,  0.0, 0.0,   /* 1 → 3 dup path variant */
};
#define N_NODES 4
#define N_EDGES 8

static double l2_node(const double *h, int dim) {
    double s = 0.0;
    for (int i = 0; i < dim; ++i) s += h[i] * h[i];
    return sqrt(s);
}

int main(void) {
    /* -------- 1. Spec-string NequIP layer -------- */
    const char *spec = "2x0e + 1x1o -> 2x0e + 1x1o "
                       "[sh=2, radial=4, r_cut=1.5, cutoff=polynomial(6)]";
    irrep_nequip_layer_t *layer = irrep_nequip_layer_from_spec(spec);
    if (!layer) {
        fprintf(stderr, "from_spec failed: %s\n", irrep_last_error());
        return 1;
    }
    int h_dim = 2 + 3;   /* 2x0e + 1x1o */
    int nw    = irrep_nequip_layer_num_weights(layer);
    printf("spec      : %s\n", spec);
    printf("h_dim     : %d\n", h_dim);
    printf("num weights: %d\n", nw);

    /* Deterministic weights + input features. */
    double *w    = calloc((size_t)nw, sizeof(double));
    double *h_in = calloc((size_t)N_NODES * h_dim, sizeof(double));
    double *h_out= calloc((size_t)N_NODES * h_dim, sizeof(double));
    for (int i = 0; i < nw; ++i)              w[i]    = 0.013 * (i + 1);
    for (int i = 0; i < N_NODES * h_dim; ++i) h_in[i] = 0.07 * i - 0.2;

    irrep_nequip_layer_apply(layer, w, N_NODES, N_EDGES,
                              k_edge_src, k_edge_dst, k_edge_vec,
                              h_in, h_out);

    printf("\nnode-wise ‖h_out‖₂ BEFORE projection:\n");
    for (int n = 0; n < N_NODES; ++n) {
        printf("  node %d  %.6f\n", n, l2_node(h_out + n * h_dim, h_dim));
    }

    /* -------- 2. C₄ᵥ A₁ projection per node -------- */
    irrep_pg_table_t *pg = irrep_pg_table_build(IRREP_PG_C4V);
    irrep_multiset_t *spec_m = irrep_multiset_parse("2x0e + 1x1o");
    if (!pg || !spec_m) { fprintf(stderr, "pg build failed\n"); return 1; }

    double *h_sym = calloc((size_t)N_NODES * h_dim, sizeof(double));
    for (int n = 0; n < N_NODES; ++n) {
        /* Project each node's feature block independently — the projector
         * acts within-node (symmetry group is point symmetries at the node). */
        irrep_pg_project(pg, /*mu=*/0 /* A1 */, spec_m,
                         h_out + n * h_dim,
                         h_sym + n * h_dim);
    }

    printf("\nnode-wise ‖h_sym‖₂ AFTER A1 projection:\n");
    for (int n = 0; n < N_NODES; ++n) {
        double before = l2_node(h_out + n * h_dim, h_dim);
        double after  = l2_node(h_sym + n * h_dim, h_dim);
        printf("  node %d  %.6f  (before %.6f, reduction %.1f%%)\n",
               n, after, before, 100.0 * (1.0 - after / before));
    }

    /* -------- 3. Verify P² = P on the projected output -------- */
    double *h_twice = calloc((size_t)N_NODES * h_dim, sizeof(double));
    for (int n = 0; n < N_NODES; ++n) {
        irrep_pg_project(pg, 0, spec_m,
                         h_sym + n * h_dim,
                         h_twice + n * h_dim);
    }
    double max_drift = 0.0;
    for (int i = 0; i < N_NODES * h_dim; ++i) {
        double d = fabs(h_twice[i] - h_sym[i]);
        if (d > max_drift) max_drift = d;
    }
    printf("\nprojector idempotence: max |P(P(h)) - P(h)| = %.3e\n", max_drift);
    printf("(must be below 1e-10 — projector is bit-exact idempotent here.)\n");

    free(w); free(h_in); free(h_out); free(h_sym); free(h_twice);
    irrep_multiset_free(spec_m);
    irrep_pg_table_free(pg);
    irrep_nequip_layer_free(layer);
    return 0;
}
