/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/nequip.h>
#include <irrep/multiset.h>
#include <irrep/so3.h>
#include <irrep/wigner_d.h>
#include <irrep/spherical_harmonics.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

static void build_real_rot_l_(int l, double *R_out,
                              double alpha, double beta, double gamma) {
    int d = 2 * l + 1;
    double _Complex U[(2 * 8 + 1) * (2 * 8 + 1)];
    double _Complex D[(2 * 8 + 1) * (2 * 8 + 1)];
    double _Complex tmp[(2 * 8 + 1) * (2 * 8 + 1)];
    irrep_sph_harm_complex_to_real(l, U);
    irrep_wigner_D_matrix(l, D, alpha, beta, gamma);
    for (int i = 0; i < d; ++i)
        for (int j = 0; j < d; ++j) {
            double _Complex s = 0;
            for (int k = 0; k < d; ++k) s += U[i * d + k] * D[k * d + j];
            tmp[i * d + j] = s;
        }
    for (int i = 0; i < d; ++i)
        for (int j = 0; j < d; ++j) {
            double _Complex s = 0;
            for (int k = 0; k < d; ++k) s += tmp[i * d + k] * conj(U[j * d + k]);
            R_out[i * d + j] = creal(s);
        }
}

/* Apply the real-basis D(R) on a multiset feature vector in place. */
static void rotate_multiset_(const irrep_multiset_t *m,
                             double alpha, double beta, double gamma,
                             const double *in, double *out) {
    int offset = 0;
    for (int t = 0; t < m->num_terms; ++t) {
        int l    = m->labels[t].l;
        int d    = 2 * l + 1;
        int mult = m->multiplicities[t];
        double R[(2 * 8 + 1) * (2 * 8 + 1)];
        build_real_rot_l_(l, R, alpha, beta, gamma);
        for (int u = 0; u < mult; ++u) {
            const double *src = in  + offset + u * d;
            double       *dst = out + offset + u * d;
            for (int i = 0; i < d; ++i) {
                double s = 0;
                for (int j = 0; j < d; ++j) s += R[i * d + j] * src[j];
                dst[i] = s;
            }
        }
        offset += mult * d;
    }
}

int main(void) {
    IRREP_TEST_START("nequip");

    irrep_multiset_t *h_in  = irrep_multiset_parse("2x0e + 1x1o");
    irrep_multiset_t *h_out = irrep_multiset_parse("2x0e + 1x1o");

    irrep_nequip_layer_t *layer = irrep_nequip_layer_build(
        h_in, /*l_sh_max=*/2, /*n_radial=*/4, /*r_cut=*/3.0,
        IRREP_NEQUIP_CUTOFF_POLYNOMIAL, /*cutoff_poly_p=*/6, h_out);
    IRREP_ASSERT(layer != NULL);
    int nw = irrep_nequip_layer_num_weights(layer);
    IRREP_ASSERT(nw > 0);

    /* Simple 3-node line graph: 0 → 1, 1 → 2, plus back edges. */
    const int n_nodes = 3;
    const int n_edges = 4;
    int edge_src[4] = { 0, 1, 1, 2 };
    int edge_dst[4] = { 1, 0, 2, 1 };
    double edge_vec[12] = {
         0.5, -0.3,  0.8,     /* 0 → 1 */
        -0.5,  0.3, -0.8,     /* 1 → 0 */
         0.2,  0.7, -0.4,     /* 1 → 2 */
        -0.2, -0.7,  0.4,     /* 2 → 1 */
    };
    int h_in_dim = h_in->total_dim;
    int h_out_dim = h_out->total_dim;

    double *h_in_buf  = calloc((size_t)n_nodes * h_in_dim,  sizeof(double));
    double *h_out_buf = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
    double *w         = calloc((size_t)nw, sizeof(double));
    for (int i = 0; i < n_nodes * h_in_dim; ++i) h_in_buf[i] = 0.17 * i - 0.3;
    for (int i = 0; i < nw; ++i)                 w[i]        = 0.23 * i - 0.5;

    irrep_nequip_layer_apply(layer, w, n_nodes, n_edges,
                             edge_src, edge_dst, edge_vec, h_in_buf, h_out_buf);

    /* Output should be finite and at least the first dst node gets something. */
    for (int i = 0; i < n_nodes * h_out_dim; ++i) IRREP_ASSERT(isfinite(h_out_buf[i]));

    /* -------- Equivariance: rotate input + edge vectors, apply, compare ------- */
    double alpha = 0.4, beta = 0.9, gamma = 1.1;
    /* Rotate h_in features on each node. */
    double *h_in_rot  = calloc((size_t)n_nodes * h_in_dim,  sizeof(double));
    for (int i = 0; i < n_nodes; ++i) {
        rotate_multiset_(h_in, alpha, beta, gamma,
                         h_in_buf + (size_t)i * h_in_dim,
                         h_in_rot + (size_t)i * h_in_dim);
    }
    /* Rotate edge vectors: v' = R · v in cartesian. */
    double R3[9];
    build_real_rot_l_(1, R3, alpha, beta, gamma);
    /* R3 is the real-basis l=1 rotation in (y, z, x) layout. Convert to the
     * cartesian (x, y, z) basis for edge_vec: permute axes. */
    /* Real-basis index 0=y, 1=z, 2=x → cartesian permutation matrix P. */
    double edge_vec_rot[12];
    for (int e = 0; e < n_edges; ++e) {
        double vx = edge_vec[e * 3 + 0];
        double vy = edge_vec[e * 3 + 1];
        double vz = edge_vec[e * 3 + 2];
        /* Express in real-basis layout (y, z, x): */
        double v_real[3] = { vy, vz, vx };
        /* Rotate: */
        double w_real[3] = { 0, 0, 0 };
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) w_real[i] += R3[i * 3 + j] * v_real[j];
        /* Back to (x, y, z) layout: */
        edge_vec_rot[e * 3 + 0] = w_real[2];   /* x */
        edge_vec_rot[e * 3 + 1] = w_real[0];   /* y */
        edge_vec_rot[e * 3 + 2] = w_real[1];   /* z */
    }

    double *h_out_rot_in = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
    irrep_nequip_layer_apply(layer, w, n_nodes, n_edges,
                             edge_src, edge_dst, edge_vec_rot, h_in_rot,
                             h_out_rot_in);

    /* Rotate the direct output. */
    double *h_out_rotated = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
    for (int i = 0; i < n_nodes; ++i) {
        rotate_multiset_(h_out, alpha, beta, gamma,
                         h_out_buf     + (size_t)i * h_out_dim,
                         h_out_rotated + (size_t)i * h_out_dim);
    }

    double max_err = 0.0;
    for (int i = 0; i < n_nodes * h_out_dim; ++i) {
        double d = fabs(h_out_rot_in[i] - h_out_rotated[i]);
        if (d > max_err) max_err = d;
    }
    IRREP_ASSERT(max_err < 1e-10);

    /* -------- Backward pass: finite-difference cross-check on one node ------- */
    {
        double *grad_h_out = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
        double *grad_h_in  = calloc((size_t)n_nodes * h_in_dim,  sizeof(double));
        double *grad_w     = calloc((size_t)nw, sizeof(double));
        for (int i = 0; i < n_nodes * h_out_dim; ++i) grad_h_out[i] = 0.11 * i + 0.2;

        irrep_nequip_layer_apply_backward(layer, w, n_nodes, n_edges,
                                           edge_src, edge_dst, edge_vec,
                                           h_in_buf, grad_h_out,
                                           grad_h_in, grad_w);

        double h = 1e-7;
        /* FD on h_in[0][0]. */
        {
            double save = h_in_buf[0];
            double *cp = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
            double *cm = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
            h_in_buf[0] = save + h;
            irrep_nequip_layer_apply(layer, w, n_nodes, n_edges,
                                      edge_src, edge_dst, edge_vec, h_in_buf, cp);
            h_in_buf[0] = save - h;
            irrep_nequip_layer_apply(layer, w, n_nodes, n_edges,
                                      edge_src, edge_dst, edge_vec, h_in_buf, cm);
            h_in_buf[0] = save;
            double fd = 0;
            for (int j = 0; j < n_nodes * h_out_dim; ++j)
                fd += grad_h_out[j] * (cp[j] - cm[j]) / (2 * h);
            IRREP_ASSERT(fabs(grad_h_in[0] - fd) < 1e-6);
            free(cp); free(cm);
        }
        /* FD on weights[0]. */
        {
            double save = w[0];
            double *cp = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
            double *cm = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
            w[0] = save + h;
            irrep_nequip_layer_apply(layer, w, n_nodes, n_edges,
                                      edge_src, edge_dst, edge_vec, h_in_buf, cp);
            w[0] = save - h;
            irrep_nequip_layer_apply(layer, w, n_nodes, n_edges,
                                      edge_src, edge_dst, edge_vec, h_in_buf, cm);
            w[0] = save;
            double fd = 0;
            for (int j = 0; j < n_nodes * h_out_dim; ++j)
                fd += grad_h_out[j] * (cp[j] - cm[j]) / (2 * h);
            IRREP_ASSERT(fabs(grad_w[0] - fd) < 1e-6);
            free(cp); free(cm);
        }
        free(grad_h_out); free(grad_h_in); free(grad_w);
    }

    /* -------- Edge-geometry gradient via _apply_forces vs finite difference. --------
     *
     * Exercises the full chain rule (SH gradient, RBF derivative, cutoff
     * derivative) across every edge and axis, so a sign or factor mistake on
     * any of the three legs surfaces here rather than hiding behind a
     * single-edge spot check. */
    {
        double *grad_h_out  = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
        double *grad_edge   = calloc((size_t)n_edges * 3,         sizeof(double));
        for (int i = 0; i < n_nodes * h_out_dim; ++i) grad_h_out[i] = 0.07 * i + 0.1;

        irrep_nequip_layer_apply_forces(layer, w, n_nodes, n_edges,
                                         edge_src, edge_dst, edge_vec,
                                         h_in_buf, grad_h_out, grad_edge);

        double h = 1e-6;
        double edge_fd[12];
        /* Every edge × every axis — 12 FD checks, not just the first edge. */
        for (int e = 0; e < n_edges; ++e) {
            for (int axis = 0; axis < 3; ++axis) {
                memcpy(edge_fd, edge_vec, sizeof(edge_fd));
                double save = edge_fd[e * 3 + axis];

                double *hp = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
                double *hm = calloc((size_t)n_nodes * h_out_dim, sizeof(double));
                edge_fd[e * 3 + axis] = save + h;
                irrep_nequip_layer_apply(layer, w, n_nodes, n_edges,
                                          edge_src, edge_dst, edge_fd, h_in_buf, hp);
                edge_fd[e * 3 + axis] = save - h;
                irrep_nequip_layer_apply(layer, w, n_nodes, n_edges,
                                          edge_src, edge_dst, edge_fd, h_in_buf, hm);

                double fd = 0.0;
                for (int k = 0; k < n_nodes * h_out_dim; ++k) {
                    fd += grad_h_out[k] * (hp[k] - hm[k]) / (2 * h);
                }
                IRREP_ASSERT(fabs(grad_edge[e * 3 + axis] - fd) < 1e-6);
                free(hp); free(hm);
            }
        }

        free(grad_h_out); free(grad_edge);
    }

    /* -------- Forces path across both cutoff variants and l_sh_max ∈ {0, 3}.
     *
     * The polynomial-cutoff test above covers `l_sh_max=2`. Rebuild tiny
     * layers with the cosine cutoff and with l_sh_max at the boundaries
     * (l=0 has no tangential SH; l=3 recomputes gradients with longer chain)
     * so regressions in those code paths don't slip through.
     */
    {
        struct { irrep_nequip_cutoff_t kind; int l_sh; int poly_p; } cases[] = {
            { IRREP_NEQUIP_CUTOFF_COSINE,     1, 0 },
            { IRREP_NEQUIP_CUTOFF_COSINE,     3, 0 },
            { IRREP_NEQUIP_CUTOFF_POLYNOMIAL, 0, 6 },
            { IRREP_NEQUIP_CUTOFF_POLYNOMIAL, 3, 6 },
        };
        for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
            irrep_multiset_t *hi = irrep_multiset_parse("1x0e + 1x1o");
            irrep_multiset_t *ho = irrep_multiset_parse("1x0e + 1x1o");
            irrep_nequip_layer_t *lc = irrep_nequip_layer_build(
                hi, cases[c].l_sh, /*n_radial=*/3, /*r_cut=*/3.0,
                cases[c].kind, cases[c].poly_p, ho);
            IRREP_ASSERT(lc != NULL);

            int lc_nw = irrep_nequip_layer_num_weights(lc);
            double *lw = calloc((size_t)lc_nw, sizeof(double));
            for (int i = 0; i < lc_nw; ++i) lw[i] = 0.12 * i - 0.3;

            const int nn = 2, ne = 2;
            int es[2] = { 0, 1 }, ed[2] = { 1, 0 };
            double ev[6] = { 0.4, -0.2, 0.5,  -0.4, 0.2, -0.5 };
            double *hin  = calloc((size_t)nn * hi->total_dim, sizeof(double));
            for (int i = 0; i < nn * hi->total_dim; ++i) hin[i] = 0.11 * i + 0.2;

            double *gho = calloc((size_t)nn * ho->total_dim, sizeof(double));
            for (int i = 0; i < nn * ho->total_dim; ++i) gho[i] = 0.09 * i - 0.1;
            double *ge = calloc((size_t)ne * 3, sizeof(double));

            irrep_nequip_layer_apply_forces(lc, lw, nn, ne, es, ed, ev, hin, gho, ge);

            double h = 1e-6;
            for (int e = 0; e < ne; ++e) {
                for (int axis = 0; axis < 3; ++axis) {
                    double ev_fd[6]; memcpy(ev_fd, ev, sizeof(ev));
                    double save = ev_fd[e * 3 + axis];
                    double *hp = calloc((size_t)nn * ho->total_dim, sizeof(double));
                    double *hm = calloc((size_t)nn * ho->total_dim, sizeof(double));
                    ev_fd[e * 3 + axis] = save + h;
                    irrep_nequip_layer_apply(lc, lw, nn, ne, es, ed, ev_fd, hin, hp);
                    ev_fd[e * 3 + axis] = save - h;
                    irrep_nequip_layer_apply(lc, lw, nn, ne, es, ed, ev_fd, hin, hm);
                    double fd = 0.0;
                    for (int k = 0; k < nn * ho->total_dim; ++k) {
                        fd += gho[k] * (hp[k] - hm[k]) / (2 * h);
                    }
                    IRREP_ASSERT(fabs(ge[e * 3 + axis] - fd) < 1e-6);
                    free(hp); free(hm);
                }
            }

            free(lw); free(hin); free(gho); free(ge);
            irrep_nequip_layer_free(lc);
            irrep_multiset_free(hi);
            irrep_multiset_free(ho);
        }
    }

    free(h_in_buf); free(h_out_buf); free(w);
    free(h_in_rot); free(h_out_rot_in); free(h_out_rotated);
    irrep_nequip_layer_free(layer);
    irrep_multiset_free(h_in);
    irrep_multiset_free(h_out);

    return IRREP_TEST_END();
}
