/* SPDX-License-Identifier: MIT */
/* NequIP descriptor benchmarks on the shared downstream lattice shapes —
 * see `tests/test_downstream_compat/lattice_connectivity.json`.
 *
 * Sweeps edge counts {64, 256, 1024, 4096, 16384} matching the downstream
 * torque-net's operational envelope from unit-test through
 * Paper-4-device-simulation scale. Ships apply / backward / forces timings
 * per shape.
 *
 * Deliberately does NOT run the 262144-edge micromagnetic shape by default
 * — per the coordination doc that's a µMAG follow-up regime and forcing it
 * into every `make bench` run makes CI slow without buying much signal.
 * Flip `IRREP_BENCH_INCLUDE_MUMAG=1` in the environment to include it.
 */
#include "harness.h"
#include <irrep/multiset.h>
#include <irrep/nequip.h>
#include "../tests/test_downstream_compat/lattice_loader.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char *size_name;        /* "small" / "medium" / "large" */
    const char *h_in_spec;
    const char *h_out_spec;
    int         l_sh_max;
    int         n_radial;
    int         cutoff_p;
} shape_cfg_t;

typedef struct {
    const char *shape_key;         /* key in lattice_connectivity.json */
    long        iterations;        /* scale iteration count with edge count */
} lattice_cfg_t;

static const char *k_json_path =
    "tests/test_downstream_compat/lattice_connectivity.json";

static const shape_cfg_t k_shapes[] = {
    { "small",
      "2x0e + 1x1o",
      "1x1o",
      /*sh=*/2, /*radial=*/4, /*cutoff_p=*/6 },
    { "medium",
      "4x0e + 2x1o + 1x2e",
      "2x0e + 1x1o",
      /*sh=*/3, /*radial=*/6, /*cutoff_p=*/6 },
    { "large",
      "8x0e + 4x1o + 2x2e + 1x3o",
      "4x0e + 2x1o + 1x2e",
      /*sh=*/4, /*radial=*/8, /*cutoff_p=*/6 },
};
#define N_SHAPES ((int)(sizeof(k_shapes) / sizeof(k_shapes[0])))

/* (iterations chosen so each call completes in ~1–3 s on Apple Silicon) */
static const lattice_cfg_t k_base_lattices[] = {
    { "4x4_periodic_64_edges",      2000 },
    { "8x8_periodic_256_edges",      500 },
    { "16x16_periodic_1024_edges",   100 },
    { "32x32_periodic_4096_edges",    20 },
    { "64x64_periodic_16384_edges",    5 },
};
#define N_BASE_LATTICES ((int)(sizeof(k_base_lattices) / sizeof(k_base_lattices[0])))

static const lattice_cfg_t k_mumag_lattice =
    { "256x256_periodic_262144_edges", 1 };

/* Fill `buf` with deterministic per-node features. Range-limited to avoid
 * FP over/underflow and to keep the TP short-circuit `if (wt == 0.0)` path
 * from triggering on sparse input signal. */
static void fill_h_in(double *buf, int n_nodes, int h_in_dim) {
    for (int i = 0; i < n_nodes * h_in_dim; ++i) {
        buf[i] = 0.021 * (double)i - 0.4;
    }
}

static void fill_weights(double *w, int nw) {
    /* Placeholder init, dense — matches the 1.1 baseline so regressions on
     * the apply/backward/forces paths stay comparable. A realistic-init
     * suite using scripts/compute_nequip_init.c can be added later. */
    for (int i = 0; i < nw; ++i) w[i] = 0.013 * (double)(i + 1);
}

static void run_one(const shape_cfg_t *shape, const lattice_cfg_t *lat) {
    lattice_shape_t lattice;
    if (lattice_shape_load(k_json_path, lat->shape_key, &lattice) != 0) {
        fprintf(stderr, "failed to load %s\n", lat->shape_key);
        return;
    }

    irrep_multiset_t *h_in  = irrep_multiset_parse(shape->h_in_spec);
    irrep_multiset_t *h_out = irrep_multiset_parse(shape->h_out_spec);
    if (!h_in || !h_out) { lattice_shape_free(&lattice); return; }

    irrep_nequip_layer_t *layer = irrep_nequip_layer_build(
        h_in, shape->l_sh_max, shape->n_radial, /*r_cut=*/3.0,
        IRREP_NEQUIP_CUTOFF_POLYNOMIAL, shape->cutoff_p, h_out);
    if (!layer) {
        irrep_multiset_free(h_in); irrep_multiset_free(h_out);
        lattice_shape_free(&lattice); return;
    }

    int nw = irrep_nequip_layer_num_weights(layer);
    double *w = calloc((size_t)nw, sizeof(double));
    fill_weights(w, nw);

    double *h_in_buf  = calloc((size_t)lattice.num_nodes * h_in->total_dim,
                               sizeof(double));
    double *h_out_buf = calloc((size_t)lattice.num_nodes * h_out->total_dim,
                               sizeof(double));
    fill_h_in(h_in_buf, lattice.num_nodes, h_in->total_dim);

    const long iters = lat->iterations;

    /* Forward. */
    double t0 = irrep_bench_now_ns();
    for (long i = 0; i < iters; ++i) {
        irrep_nequip_layer_apply(layer, w,
                                  lattice.num_nodes, lattice.num_edges,
                                  lattice.edge_src, lattice.edge_dst,
                                  lattice.edge_vec, h_in_buf, h_out_buf);
    }
    double t1 = irrep_bench_now_ns();
    char name[192];
    snprintf(name, sizeof(name), "nequip_apply_%s_%dE",
             shape->size_name, lattice.num_edges);
    irrep_bench_report(name, iters, t1 - t0, lattice.num_edges);

    /* Backward. */
    double *gho = calloc((size_t)lattice.num_nodes * h_out->total_dim,
                         sizeof(double));
    double *ghi = calloc((size_t)lattice.num_nodes * h_in->total_dim,
                         sizeof(double));
    double *gw  = calloc((size_t)nw, sizeof(double));
    for (int i = 0; i < lattice.num_nodes * h_out->total_dim; ++i) {
        gho[i] = 0.011 * (double)i + 0.2;
    }
    t0 = irrep_bench_now_ns();
    for (long i = 0; i < iters; ++i) {
        memset(ghi, 0, (size_t)lattice.num_nodes * h_in->total_dim * sizeof(double));
        memset(gw,  0, (size_t)nw * sizeof(double));
        irrep_nequip_layer_apply_backward(layer, w,
                                           lattice.num_nodes, lattice.num_edges,
                                           lattice.edge_src, lattice.edge_dst,
                                           lattice.edge_vec,
                                           h_in_buf, gho, ghi, gw);
    }
    t1 = irrep_bench_now_ns();
    snprintf(name, sizeof(name), "nequip_backward_%s_%dE",
             shape->size_name, lattice.num_edges);
    irrep_bench_report(name, iters, t1 - t0, lattice.num_edges);

    /* Forces — edge-geometry gradient. */
    double *grad_edge = calloc((size_t)lattice.num_edges * 3, sizeof(double));
    t0 = irrep_bench_now_ns();
    for (long i = 0; i < iters; ++i) {
        memset(grad_edge, 0, (size_t)lattice.num_edges * 3 * sizeof(double));
        irrep_nequip_layer_apply_forces(layer, w,
                                         lattice.num_nodes, lattice.num_edges,
                                         lattice.edge_src, lattice.edge_dst,
                                         lattice.edge_vec,
                                         h_in_buf, gho, grad_edge);
    }
    t1 = irrep_bench_now_ns();
    snprintf(name, sizeof(name), "nequip_forces_%s_%dE",
             shape->size_name, lattice.num_edges);
    irrep_bench_report(name, iters, t1 - t0, lattice.num_edges);

    free(w); free(h_in_buf); free(h_out_buf);
    free(gho); free(ghi); free(gw); free(grad_edge);
    irrep_nequip_layer_free(layer);
    irrep_multiset_free(h_in);
    irrep_multiset_free(h_out);
    lattice_shape_free(&lattice);
}

int main(void) {
    int include_mumag = 0;
    const char *env = getenv("IRREP_BENCH_INCLUDE_MUMAG");
    if (env && env[0] && env[0] != '0') include_mumag = 1;

    for (int s = 0; s < N_SHAPES; ++s) {
        for (int l = 0; l < N_BASE_LATTICES; ++l) {
            run_one(&k_shapes[s], &k_base_lattices[l]);
        }
        if (include_mumag) {
            run_one(&k_shapes[s], &k_mumag_lattice);
        }
    }
    return 0;
}
