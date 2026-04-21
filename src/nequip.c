/* SPDX-License-Identifier: MIT */
/* NequIP message-passing layer: the first-class composition of edge SH,
 * radial basis, cutoff, UVW tensor product, and neighbour aggregation. */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/multiset.h>
#include <irrep/nequip.h>
#include <irrep/radial.h>
#include <irrep/spherical_harmonics.h>
#include <irrep/tensor_product.h>

extern void irrep_set_error_(const char *fmt, ...);

struct irrep_nequip_layer {
    int                   l_sh_max;
    int                   n_radial;
    double                r_cut;
    irrep_nequip_cutoff_t cutoff_kind;
    int                   cutoff_poly_p;

    irrep_multiset_t     *sh_mset; /* "1x0e + 1x1o + ... + 1xLpar" */
    irrep_multiset_t     *hidden_in_copy;
    irrep_multiset_t     *hidden_out_copy;

    tp_descriptor_t      *tp;
    int                   sh_dim; /* (l_sh_max+1)² */
    int                   h_in_dim;
    int                   h_out_dim;
};

/* Build "1x0e + 1x1o + 1x2e + ..." up to l_sh_max — parity follows (-1)^l
 * so that the physical cartesian SH transforms as a natural O(3) irrep. */
static irrep_multiset_t *build_sh_multiset_(int l_sh_max) {
    irrep_multiset_t *m = irrep_multiset_new(l_sh_max + 1);
    if (!m)
        return NULL;
    for (int l = 0; l <= l_sh_max; ++l) {
        irrep_label_t lbl = {.l = l, .parity = (l & 1) ? IRREP_ODD : IRREP_EVEN};
        if (irrep_multiset_append(m, lbl, 1) != IRREP_OK) {
            irrep_multiset_free(m);
            return NULL;
        }
    }
    return m;
}

static irrep_multiset_t *clone_multiset_(const irrep_multiset_t *src) {
    if (!src)
        return NULL;
    irrep_multiset_t *dup = irrep_multiset_new(src->num_terms);
    if (!dup)
        return NULL;
    for (int i = 0; i < src->num_terms; ++i) {
        if (irrep_multiset_append(dup, src->labels[i], src->multiplicities[i]) != IRREP_OK) {
            irrep_multiset_free(dup);
            return NULL;
        }
    }
    return dup;
}

irrep_nequip_layer_t *irrep_nequip_layer_build(const irrep_multiset_t *hidden_in, int l_sh_max,
                                               int n_radial, double r_cut,
                                               irrep_nequip_cutoff_t cutoff_kind, int cutoff_poly_p,
                                               const irrep_multiset_t *hidden_out) {

    if (!hidden_in || !hidden_out) {
        irrep_set_error_("irrep_nequip_layer_build: NULL multiset input");
        return NULL;
    }
    if (l_sh_max < 0 || r_cut <= 0.0 || n_radial <= 0) {
        irrep_set_error_("irrep_nequip_layer_build: invalid numeric arg");
        return NULL;
    }

    irrep_nequip_layer_t *layer = calloc(1, sizeof(*layer));
    if (!layer)
        return NULL;
    layer->l_sh_max = l_sh_max;
    layer->n_radial = n_radial;
    layer->r_cut = r_cut;
    layer->cutoff_kind = cutoff_kind;
    layer->cutoff_poly_p = cutoff_poly_p > 0 ? cutoff_poly_p : 6;

    layer->sh_mset = build_sh_multiset_(l_sh_max);
    layer->hidden_in_copy = clone_multiset_(hidden_in);
    layer->hidden_out_copy = clone_multiset_(hidden_out);
    if (!layer->sh_mset || !layer->hidden_in_copy || !layer->hidden_out_copy) {
        irrep_nequip_layer_free(layer);
        return NULL;
    }
    layer->sh_dim = (l_sh_max + 1) * (l_sh_max + 1);
    layer->h_in_dim = hidden_in->total_dim;
    layer->h_out_dim = hidden_out->total_dim;

    /* Enumerate and build the UVW tensor product: hidden_in ⊗ sh → hidden_out. */
    int  n_paths = irrep_tp_enumerate_paths(layer->hidden_in_copy, layer->sh_mset,
                                            layer->hidden_out_copy, NULL, 0);
    int *paths = malloc((size_t)(n_paths > 0 ? n_paths : 1) * 3 * sizeof(int));
    if (!paths) {
        irrep_nequip_layer_free(layer);
        return NULL;
    }
    irrep_tp_enumerate_paths(layer->hidden_in_copy, layer->sh_mset, layer->hidden_out_copy, paths,
                             n_paths);
    layer->tp = irrep_tp_build_uvw(layer->hidden_in_copy, layer->sh_mset, layer->hidden_out_copy,
                                   paths, n_paths);
    free(paths);
    if (!layer->tp) {
        irrep_nequip_layer_free(layer);
        return NULL;
    }
    return layer;
}

void irrep_nequip_layer_free(irrep_nequip_layer_t *layer) {
    if (!layer)
        return;
    irrep_tp_free(layer->tp);
    irrep_multiset_free(layer->sh_mset);
    irrep_multiset_free(layer->hidden_in_copy);
    irrep_multiset_free(layer->hidden_out_copy);
    free(layer);
}

int irrep_nequip_layer_num_weights(const irrep_nequip_layer_t *layer) {
    return layer ? irrep_tp_num_weights_uvw(layer->tp) : 0;
}

static double eval_cutoff_(const irrep_nequip_layer_t *layer, double r) {
    if (layer->cutoff_kind == IRREP_NEQUIP_CUTOFF_COSINE) {
        return irrep_cutoff_cosine(r, layer->r_cut);
    }
    return irrep_cutoff_polynomial(r, layer->r_cut, layer->cutoff_poly_p);
}

static double eval_cutoff_d_(const irrep_nequip_layer_t *layer, double r) {
    if (layer->cutoff_kind == IRREP_NEQUIP_CUTOFF_COSINE) {
        return irrep_cutoff_cosine_d(r, layer->r_cut);
    }
    return irrep_cutoff_polynomial_d(r, layer->r_cut, layer->cutoff_poly_p);
}

void irrep_nequip_layer_apply(const irrep_nequip_layer_t *layer, const double *tp_weights,
                              int n_nodes, int n_edges, const int *edge_src, const int *edge_dst,
                              const double *edge_vec, const double *h_in, double *h_out) {

    if (!layer)
        return;
    memset(h_out, 0, (size_t)n_nodes * (size_t)layer->h_out_dim * sizeof(double));
    if (n_edges <= 0)
        return;

    /* Per-edge scratch — rebuilt each iteration for simplicity. Batched
     * NEON/AVX variants land in the SIMD sprint. */
    double *sh_buf = malloc((size_t)layer->sh_dim * sizeof(double));
    double *msg = malloc((size_t)layer->h_out_dim * sizeof(double));
    if (!sh_buf || !msg) {
        free(sh_buf);
        free(msg);
        return;
    }

    for (int e = 0; e < n_edges; ++e) {
        int src = edge_src[e];
        int dst = edge_dst[e];
        if (src < 0 || src >= n_nodes)
            continue;
        if (dst < 0 || dst >= n_nodes)
            continue;

        double rx = edge_vec[e * 3 + 0];
        double ry = edge_vec[e * 3 + 1];
        double rz = edge_vec[e * 3 + 2];
        double r = sqrt(rx * rx + ry * ry + rz * rz);
        if (r <= 0.0 || r >= layer->r_cut)
            continue;
        double inv_r = 1.0 / r;
        double r_hat[3] = {rx * inv_r, ry * inv_r, rz * inv_r};

        irrep_sph_harm_cart_all(layer->l_sh_max, sh_buf, r_hat);

        double cutoff = eval_cutoff_(layer, r);
        double rbf_sum = 0.0;
        for (int n = 1; n <= layer->n_radial; ++n) {
            rbf_sum += irrep_rbf_bessel(n, r, layer->r_cut);
        }
        double scale = cutoff * rbf_sum;
        if (scale == 0.0)
            continue;

        /* msg = tp(h_in[src], Y(r̂_ij)) */
        irrep_tp_apply_uvw(layer->tp, tp_weights, h_in + (size_t)src * layer->h_in_dim, sh_buf,
                           msg);
        /* accumulate into h_out[dst]: h_out[dst] += scale · msg */
        double *target = h_out + (size_t)dst * layer->h_out_dim;
        for (int k = 0; k < layer->h_out_dim; ++k)
            target[k] += scale * msg[k];
    }

    free(sh_buf);
    free(msg);
}

void irrep_nequip_layer_apply_backward(const irrep_nequip_layer_t *layer, const double *tp_weights,
                                       int n_nodes, int n_edges, const int *edge_src,
                                       const int *edge_dst, const double *edge_vec,
                                       const double *h_in, const double *grad_h_out,
                                       double *grad_h_in, double *grad_tp_weights) {

    if (!layer)
        return;
    (void)n_nodes;

    /* Scratch per edge. */
    double *sh_buf = calloc((size_t)layer->sh_dim, sizeof(double));
    double *scaled_g_out = calloc((size_t)layer->h_out_dim, sizeof(double));
    double *ga_local = calloc((size_t)layer->h_in_dim, sizeof(double));
    if (!sh_buf || !scaled_g_out || !ga_local) {
        free(sh_buf);
        free(scaled_g_out);
        free(ga_local);
        return;
    }

    for (int e = 0; e < n_edges; ++e) {
        int src = edge_src[e];
        int dst = edge_dst[e];
        if (src < 0 || src >= n_nodes)
            continue;
        if (dst < 0 || dst >= n_nodes)
            continue;

        double rx = edge_vec[e * 3 + 0];
        double ry = edge_vec[e * 3 + 1];
        double rz = edge_vec[e * 3 + 2];
        double r = sqrt(rx * rx + ry * ry + rz * rz);
        if (r <= 0.0 || r >= layer->r_cut)
            continue;
        double inv_r = 1.0 / r;
        double r_hat[3] = {rx * inv_r, ry * inv_r, rz * inv_r};

        irrep_sph_harm_cart_all(layer->l_sh_max, sh_buf, r_hat);

        double cutoff = eval_cutoff_(layer, r);
        double rbf_sum = 0.0;
        for (int n = 1; n <= layer->n_radial; ++n) {
            rbf_sum += irrep_rbf_bessel(n, r, layer->r_cut);
        }
        double scale = cutoff * rbf_sum;
        if (scale == 0.0)
            continue;

        /* h_out[dst] += scale · msg(src, e); so
         * ∂L/∂msg = scale · ∂L/∂h_out[dst] and we pass that into TP backward. */
        const double *go_dst = grad_h_out + (size_t)dst * layer->h_out_dim;
        for (int k = 0; k < layer->h_out_dim; ++k) {
            scaled_g_out[k] = scale * go_dst[k];
        }

        /* Zero the per-edge a-gradient scratch; accumulate into grad_h_in[src]. */
        memset(ga_local, 0, (size_t)layer->h_in_dim * sizeof(double));
        irrep_tp_apply_uvw_backward(layer->tp, tp_weights, h_in + (size_t)src * layer->h_in_dim,
                                    sh_buf, scaled_g_out, grad_h_in ? ga_local : NULL,
                                    NULL, /* skip grad of SH */
                                    grad_tp_weights);
        if (grad_h_in) {
            double *gh = grad_h_in + (size_t)src * layer->h_in_dim;
            for (int k = 0; k < layer->h_in_dim; ++k)
                gh[k] += ga_local[k];
        }
    }

    free(sh_buf);
    free(scaled_g_out);
    free(ga_local);
}

/* -------------------------------------------------------------------------- *
 * Edge-geometry gradient                                                     *
 *                                                                            *
 *   h_out[dst, k] = Σ_e  scale(r_e) · tp(h_src(e), Y(r̂_e))[k]                *
 *                                                                            *
 *   ∂h_out[dst]/∂r_e = (∂scale/∂r_e) · tp_out(e) + scale(r_e) · (∂tp_out/∂r_e) *
 *                                                                            *
 *   ∂r/∂edge_vec[e,axis]  = r̂[axis]                                          *
 *   ∂r̂[axis']/∂edge_vec[e,axis] = (δ_{axis,axis'} − r̂[axis]·r̂[axis']) / r.  *
 *                                                                            *
 * `Y(r̂)`'s gradient is tangent to S² (test_spherical_harmonics asserts this  *
 * to 1e-8 for l ≤ 3), so the radial projection term `r̂[axis]·(r̂·∇_r̂ Y_m)`  *
 * vanishes and `∂Y_m/∂edge_vec[e,axis] = (1/r) · ∂Y_m/∂r̂[axis]`.             *
 * -------------------------------------------------------------------------- */

void irrep_nequip_layer_apply_forces(const irrep_nequip_layer_t *layer, const double *tp_weights,
                                     int n_nodes, int n_edges, const int *edge_src,
                                     const int *edge_dst, const double *edge_vec,
                                     const double *h_in, const double *grad_h_out,
                                     double *grad_edge_vec) {

    if (!layer || !grad_edge_vec)
        return;
    (void)n_nodes;

    const int sh_dim = layer->sh_dim;
    const int h_out_dim = layer->h_out_dim;
    const int h_in_dim = layer->h_in_dim;
    const int l_sh_max = layer->l_sh_max;

    double   *sh_buf = calloc((size_t)sh_dim, sizeof(double));
    double   *grad_sh = calloc((size_t)sh_dim, sizeof(double));
    double   *msg = calloc((size_t)h_out_dim, sizeof(double));
    double   *scaled_g_out = calloc((size_t)h_out_dim, sizeof(double));
    /* Gradient of each SH component w.r.t. r̂ axes. Layout:
     *   sh_grad[axis * sh_dim + (l² + m + l)] = ∂Y_{l,m}/∂r̂[axis]. */
    double *sh_grad = calloc(3u * (size_t)sh_dim, sizeof(double));
    if (!sh_buf || !grad_sh || !msg || !scaled_g_out || !sh_grad) {
        free(sh_buf);
        free(grad_sh);
        free(msg);
        free(scaled_g_out);
        free(sh_grad);
        return;
    }

    for (int e = 0; e < n_edges; ++e) {
        int src = edge_src[e];
        int dst = edge_dst[e];
        if (src < 0 || src >= n_nodes)
            continue;
        if (dst < 0 || dst >= n_nodes)
            continue;

        double rx = edge_vec[e * 3 + 0];
        double ry = edge_vec[e * 3 + 1];
        double rz = edge_vec[e * 3 + 2];
        double r = sqrt(rx * rx + ry * ry + rz * rz);
        if (r <= 0.0 || r >= layer->r_cut)
            continue;
        double inv_r = 1.0 / r;
        double r_hat[3] = {rx * inv_r, ry * inv_r, rz * inv_r};

        /* Recompute per-edge forward state. */
        irrep_sph_harm_cart_all(l_sh_max, sh_buf, r_hat);
        irrep_sph_harm_cart_all_grad_batch(l_sh_max, 1, r_hat, sh_grad);

        double cutoff = eval_cutoff_(layer, r);
        double cutoff_d = eval_cutoff_d_(layer, r);
        double rbf_sum = 0.0;
        double rbf_sum_d = 0.0;
        for (int n = 1; n <= layer->n_radial; ++n) {
            rbf_sum += irrep_rbf_bessel(n, r, layer->r_cut);
            rbf_sum_d += irrep_rbf_bessel_d(n, r, layer->r_cut);
        }
        double scale = cutoff * rbf_sum;
        double scale_d = cutoff_d * rbf_sum + cutoff * rbf_sum_d;
        if (scale == 0.0 && scale_d == 0.0)
            continue;

        /* Forward message so we can contract it with grad_h_out for the r̂·scale'·msg term. */
        irrep_tp_apply_uvw(layer->tp, tp_weights, h_in + (size_t)src * h_in_dim, sh_buf, msg);

        const double *go_dst = grad_h_out + (size_t)dst * h_out_dim;

        /* First term: (grad_h_out · msg) · scale' · r̂. */
        double dot_gm = 0.0;
        for (int k = 0; k < h_out_dim; ++k)
            dot_gm += go_dst[k] * msg[k];
        double coeff_r = dot_gm * scale_d;

        /* Second term: SH-side pullback. scaled_g_out = scale · grad_h_out[dst]. */
        for (int k = 0; k < h_out_dim; ++k)
            scaled_g_out[k] = scale * go_dst[k];

        /* Zero grad_sh then accumulate ∂(msg · scaled_g_out)/∂Y via TP backward. */
        memset(grad_sh, 0, (size_t)sh_dim * sizeof(double));
        irrep_tp_apply_uvw_backward(layer->tp, tp_weights, h_in + (size_t)src * h_in_dim, sh_buf,
                                    scaled_g_out,
                                    /*grad_a=*/NULL,
                                    /*grad_b=*/grad_sh,
                                    /*grad_w=*/NULL);

        /* Accumulate into grad_edge_vec[e]:
         *   coeff_r · r̂[axis] + (1/r) · Σ_m grad_sh[m] · sh_grad[axis, m]. */
        for (int axis = 0; axis < 3; ++axis) {
            double        s = 0.0;
            const double *gb_axis = sh_grad + (size_t)axis * sh_dim;
            for (int m = 0; m < sh_dim; ++m)
                s += grad_sh[m] * gb_axis[m];
            grad_edge_vec[e * 3 + axis] += coeff_r * r_hat[axis] + inv_r * s;
        }
    }

    free(sh_buf);
    free(grad_sh);
    free(msg);
    free(scaled_g_out);
    free(sh_grad);
}
