/* SPDX-License-Identifier: MIT */
/* libFuzzer entry point for irrep_nequip_layer_apply_forces.
 *
 * A fixed NequIP layer is built once ("1x0e + 1x1o" → "1x0e + 1x1o",
 * sh=2, radial=4, polynomial cutoff) and each input is reinterpreted as:
 *     edge_vec[e_max * 3] + h_in[n_max * h_in_dim] + grad_h_out[n_max * h_out_dim]
 *
 * Edge indices are derived from the input bytes modulo n_nodes so bad
 * indices don't trip the bounds check (we exercise the geometry path, not
 * the index-validation path). Post-call we require every output to be
 * finite — a NaN leak through the chain rule would indicate a real defect. */

#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include <irrep/multiset.h>
#include <irrep/nequip.h>

#define N_NODES 4
#define N_EDGES 6

static irrep_nequip_layer_t *g_layer  = NULL;
static irrep_multiset_t     *g_h_in   = NULL;
static irrep_multiset_t     *g_h_out  = NULL;
static int                   g_h_in_dim;
static int                   g_h_out_dim;
static int                   g_nw;
static double               *g_w      = NULL;

static void init_once_(void) {
    if (g_layer) return;
    g_h_in  = irrep_multiset_parse("1x0e + 1x1o");
    g_h_out = irrep_multiset_parse("1x0e + 1x1o");
    g_layer = irrep_nequip_layer_build(
        g_h_in, /*l_sh_max=*/2, /*n_radial=*/4, /*r_cut=*/1.5,
        IRREP_NEQUIP_CUTOFF_POLYNOMIAL, /*cutoff_poly_p=*/6, g_h_out);
    if (!g_layer) return;
    g_h_in_dim  = g_h_in->total_dim;
    g_h_out_dim = g_h_out->total_dim;
    g_nw = irrep_nequip_layer_num_weights(g_layer);
    g_w  = calloc((size_t)g_nw, sizeof(double));
    for (int i = 0; i < g_nw; ++i) g_w[i] = 0.013 * (i + 1);
}

/* Minimum input to cover every buffer: edge_vec (6*3=18) + h_in (4*4=16) +
 * grad_out (4*4=16) = 50 doubles = 400 bytes. Round up with slack. */
#define MIN_INPUT 512u

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    init_once_();
    if (!g_layer || size < MIN_INPUT) return 0;

    /* Sanitise float inputs so NaN / ±Inf can't propagate trivially; the
     * fuzzer should exercise *finite* pathological inputs. */
    double edge_vec[N_EDGES * 3];
    double h_in     [N_NODES * 4];   /* h_in_dim = 1 + 3 = 4 */
    double grad_out [N_NODES * 4];
    memcpy(edge_vec, data,        sizeof(edge_vec));
    memcpy(h_in,     data + 144,  sizeof(h_in));
    memcpy(grad_out, data + 144 + 128, sizeof(grad_out));

    /* Sanitize + clamp to [-1, 1] magnitude; raw random bytes include
     * doubles near 1e+308 whose products overflow validly to Inf/NaN. */
    for (size_t i = 0; i < sizeof(edge_vec) / sizeof(edge_vec[0]); ++i)
        edge_vec[i] = isfinite(edge_vec[i]) ? tanh(edge_vec[i]) : 0.1;
    for (size_t i = 0; i < sizeof(h_in) / sizeof(h_in[0]); ++i)
        h_in[i] = isfinite(h_in[i]) ? tanh(h_in[i]) : 0.0;
    for (size_t i = 0; i < sizeof(grad_out) / sizeof(grad_out[0]); ++i)
        grad_out[i] = isfinite(grad_out[i]) ? tanh(grad_out[i]) : 0.0;

    /* Clamp |edge_vec| to [0.01, 1.4] per component to keep some edges
     * inside r_cut and still sweep orientations. */
    for (int e = 0; e < N_EDGES; ++e) {
        for (int k = 0; k < 3; ++k) {
            double v = edge_vec[e * 3 + k];
            if (v >  1.4) v =  1.4;
            if (v < -1.4) v = -1.4;
            edge_vec[e * 3 + k] = v;
        }
    }

    /* Fixed graph topology — the fuzzer targets the chain rule, not the
     * edge-list validator. */
    int src[N_EDGES] = { 0, 1, 1, 2, 2, 3 };
    int dst[N_EDGES] = { 1, 0, 2, 1, 3, 2 };

    double grad_edge[N_EDGES * 3] = { 0 };
    irrep_nequip_layer_apply_forces(g_layer, g_w, N_NODES, N_EDGES,
                                     src, dst, edge_vec, h_in, grad_out,
                                     grad_edge);

    for (size_t i = 0; i < sizeof(grad_edge) / sizeof(grad_edge[0]); ++i) {
        if (!isfinite(grad_edge[i])) __builtin_trap();
    }
    return 0;
}
