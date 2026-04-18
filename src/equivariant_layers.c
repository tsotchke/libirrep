/* SPDX-License-Identifier: MIT */
/* M9: minimal equivariant-NN primitives.
 *
 * Linear-on-irreps: per matched irrep label (l, parity) there is an
 * (mult_out × mult_in) weight block that mixes channels. Blocks for
 * different (l, parity) are independent — preserves equivariance.
 *
 * RMS norm per (term, copy): normalizes the (2l+1)-dimensional block by
 * its own RMS and scales by a learnable scalar. Rotation-equivariant
 * because the RMS is a scalar (l=0) invariant of the block.
 *
 * Gate: each (term, copy) block is multiplied by a caller-supplied scalar
 * gate (typically produced by passing a scalar feature through σ or tanh).
 *
 * `in_channels` / `out_channels` in `irrep_linear_build` are accepted for
 * future uvw-mode use; M9 ignores them (set to 1 by convention).
 */

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/equivariant_layers.h>

#define UNUSED(x) ((void)(x))

/* -------------------------------------------------------------------------- *
 * Linear-on-irreps                                                           *
 * -------------------------------------------------------------------------- */

struct linear_recipe {
    int t_in, t_out;
    int l;
    int d;                 /* 2l + 1 */
    int mult_in, mult_out;
    int in_offset, out_offset;
    int weight_offset;     /* into the flat weights[] */
};

struct irrep_linear {
    int in_dim, out_dim;
    int in_channels, out_channels;
    int num_recipes;
    int num_weights;
    struct linear_recipe *recipes;
};

static int multiset_offset_(const irrep_multiset_t *m, int t) {
    int off = 0;
    for (int i = 0; i < t; ++i) off += m->multiplicities[i] * (2 * m->labels[i].l + 1);
    return off;
}

irrep_linear_t *irrep_linear_build(const irrep_multiset_t *in,
                                   const irrep_multiset_t *out,
                                   int in_channels, int out_channels) {
    if (!in || !out) return NULL;
    irrep_linear_t *lin = calloc(1, sizeof(*lin));
    if (!lin) return NULL;
    lin->in_dim        = in->total_dim;
    lin->out_dim       = out->total_dim;
    lin->in_channels   = in_channels;
    lin->out_channels  = out_channels;

    /* Count matched (l, parity) label pairs. */
    int max_pairs = in->num_terms * out->num_terms;
    lin->recipes = calloc((size_t)max_pairs, sizeof(*lin->recipes));
    if (!lin->recipes) { free(lin); return NULL; }

    int count = 0;
    int wsum  = 0;
    for (int ti = 0; ti < in->num_terms; ++ti) {
        for (int to = 0; to < out->num_terms; ++to) {
            if (in->labels[ti].l      != out->labels[to].l)      continue;
            if (in->labels[ti].parity != out->labels[to].parity) continue;
            struct linear_recipe *r = &lin->recipes[count++];
            r->t_in          = ti;
            r->t_out         = to;
            r->l             = in->labels[ti].l;
            r->d             = 2 * r->l + 1;
            r->mult_in       = in->multiplicities[ti];
            r->mult_out      = out->multiplicities[to];
            r->in_offset     = multiset_offset_(in,  ti);
            r->out_offset    = multiset_offset_(out, to);
            r->weight_offset = wsum;
            wsum += r->mult_in * r->mult_out;
        }
    }
    lin->num_recipes = count;
    lin->num_weights = wsum;
    return lin;
}

void irrep_linear_free(irrep_linear_t *lin) {
    if (!lin) return;
    free(lin->recipes);
    free(lin);
}

int irrep_linear_num_weights(const irrep_linear_t *lin) {
    return lin ? lin->num_weights : 0;
}

void irrep_linear_apply(const irrep_linear_t *lin,
                        const double *weights, const double *in, double *out) {
    if (!lin) return;
    for (int i = 0; i < lin->out_dim; ++i) out[i] = 0.0;

    for (int k = 0; k < lin->num_recipes; ++k) {
        const struct linear_recipe *r = &lin->recipes[k];
        const double *W = weights + r->weight_offset;
        for (int u_out = 0; u_out < r->mult_out; ++u_out) {
            double *o = out + r->out_offset + u_out * r->d;
            for (int u_in = 0; u_in < r->mult_in; ++u_in) {
                double w_uv = W[u_out * r->mult_in + u_in];
                if (w_uv == 0.0) continue;
                const double *v = in + r->in_offset + u_in * r->d;
                for (int m = 0; m < r->d; ++m) o[m] += w_uv * v[m];
            }
        }
    }
}

void irrep_linear_backward(const irrep_linear_t *lin,
                           const double *weights, const double *in,
                           const double *grad_out,
                           double *grad_weights, double *grad_in) {
    if (!lin) return;
    for (int k = 0; k < lin->num_recipes; ++k) {
        const struct linear_recipe *r = &lin->recipes[k];
        const double *W  = weights;
        double       *gW = grad_weights ? grad_weights + r->weight_offset : NULL;

        for (int u_out = 0; u_out < r->mult_out; ++u_out) {
            const double *go = grad_out + r->out_offset + u_out * r->d;
            for (int u_in = 0; u_in < r->mult_in; ++u_in) {
                const double *v = in + r->in_offset + u_in * r->d;
                double dw = 0.0;
                double w_uv = W[r->weight_offset + u_out * r->mult_in + u_in];
                double *gi = grad_in ? grad_in + r->in_offset + u_in * r->d : NULL;
                for (int m = 0; m < r->d; ++m) {
                    dw += v[m] * go[m];
                    if (gi) gi[m] += w_uv * go[m];
                }
                if (gW) gW[u_out * r->mult_in + u_in] += dw;
            }
        }
    }
    UNUSED(weights);   /* quiet unused warning in the branch where gW is NULL */
}

/* -------------------------------------------------------------------------- *
 * RMS norm per (term, copy)                                                  *
 *                                                                            *
 * `channels` is accepted for future uvw extension (M9 treats it as 1).       *
 * Scales layout: one scalar per (term, copy), flat in multiset order.        *
 * -------------------------------------------------------------------------- */

void irrep_norm_rms(const irrep_multiset_t *m, int channels,
                    const double *scales, const double *in, double *out) {
    UNUSED(channels);
    if (!m) return;
    int scale_idx = 0;
    int off       = 0;
    for (int t = 0; t < m->num_terms; ++t) {
        int d    = 2 * m->labels[t].l + 1;
        int mult = m->multiplicities[t];
        for (int u = 0; u < mult; ++u) {
            const double *x = in + off;
            double        *y = out + off;
            double sum_sq = 0.0;
            for (int k = 0; k < d; ++k) sum_sq += x[k] * x[k];
            double rms = sqrt(sum_sq / (double)d);
            double inv = rms > 1e-14 ? 1.0 / rms : 0.0;
            double s   = scales ? scales[scale_idx] : 1.0;
            double g   = s * inv;
            for (int k = 0; k < d; ++k) y[k] = g * x[k];
            off += d;
            scale_idx++;
        }
    }
}

void irrep_norm_rms_backward(const irrep_multiset_t *m, int channels,
                             const double *scales, const double *in,
                             const double *grad_out,
                             double *grad_scales, double *grad_in) {
    UNUSED(channels);
    if (!m) return;
    int scale_idx = 0;
    int off       = 0;
    for (int t = 0; t < m->num_terms; ++t) {
        int d    = 2 * m->labels[t].l + 1;
        int mult = m->multiplicities[t];
        for (int u = 0; u < mult; ++u) {
            const double *x  = in       + off;
            const double *go = grad_out + off;
            double       *gi = grad_in ? grad_in + off : NULL;

            double sum_sq = 0.0;
            for (int k = 0; k < d; ++k) sum_sq += x[k] * x[k];
            double rms = sqrt(sum_sq / (double)d);
            if (rms < 1e-14) {
                /* gradient undefined; leave grad_in as zero update */
                off += d; scale_idx++; continue;
            }
            double s    = scales ? scales[scale_idx] : 1.0;
            double inv  = 1.0 / rms;
            double inv3 = inv * inv * inv;

            /* d rms / d x_k = x_k / (d · rms)
             * d y_k / d x_j = s · (δ_{kj} / rms − x_k · x_j / (d · rms³))
             * grad_x_j = Σ_k go_k · s · (δ_{kj}/rms − x_k x_j / (d rms³))
             *         = s · (go_j / rms − (Σ_k go_k · x_k) · x_j / (d · rms³)) */
            double go_dot_x = 0.0;
            for (int k = 0; k < d; ++k) go_dot_x += go[k] * x[k];

            if (gi) {
                double a = s * inv;
                double b = s * go_dot_x * inv3 / (double)d;
                for (int j = 0; j < d; ++j) gi[j] += a * go[j] - b * x[j];
            }
            if (grad_scales) {
                double gs = 0.0;
                for (int k = 0; k < d; ++k) gs += go[k] * (x[k] * inv);
                grad_scales[scale_idx] += gs;
            }
            off += d;
            scale_idx++;
        }
    }
}

/* -------------------------------------------------------------------------- *
 * Gate — scalar · irrep block                                                *
 * -------------------------------------------------------------------------- */

void irrep_gate_apply(const irrep_multiset_t *m, int channels,
                      const double *scalar_gates,
                      const double *in, double *out) {
    UNUSED(channels);
    if (!m) return;
    int gate_idx = 0;
    int off      = 0;
    for (int t = 0; t < m->num_terms; ++t) {
        int d    = 2 * m->labels[t].l + 1;
        int mult = m->multiplicities[t];
        for (int u = 0; u < mult; ++u) {
            double g = scalar_gates ? scalar_gates[gate_idx] : 1.0;
            for (int k = 0; k < d; ++k) out[off + k] = g * in[off + k];
            off += d;
            gate_idx++;
        }
    }
}
