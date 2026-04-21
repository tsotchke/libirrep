/* SPDX-License-Identifier: MIT */
/* e3nn-style path-indexed tensor products.
 *
 * Two channel modes supported (see IRREP_TP_MODE_UUU / IRREP_TP_MODE_UVW):
 *
 *   UUU mode (`irrep_tp_build` / `_apply` / `_apply_weighted`):
 *     For each selected path (i_a, i_b, i_c) we require
 *     mult_a = mult_b = mult_c, and pair channels one-to-one. One scalar
 *     weight per path. Copy index u ∈ [0, mult_p):
 *
 *         c[i_c, u, m_c] += w_p · Σ_{m_a, m_b}
 *             W_p[m_a, m_b, m_c] · a[i_a, u, m_a] · b[i_b, u, m_b]
 *
 *   UVW mode (`irrep_tp_build_uvw` / `_apply_uvw` / `_apply_uvw_backward`):
 *     Independent multiplicities. Full [W, V, U] weight tensor per path.
 *
 *         c[i_c, w, m_c] += Σ_{u, v} weight_p[w, v, u] · Σ_{m_a, m_b}
 *             W_p[m_a, m_b, m_c] · a[i_a, u, m_a] · b[i_b, v, m_b]
 *
 * The path coupling tensor W_p is built from the real-basis change-of-
 * basis of the complex Clebsch-Gordan coefficients, with an
 * i^{l_a + l_b − l_c} phase correction so odd-l-sum paths (e.g. the
 * cross-product (1, 1, 1) path) produce guaranteed-real output.
 * Parity filter `p_a · p_b = p_c` is enforced by irrep_tp_enumerate_paths.
 *
 * Backward accumulates gradients into grad_a / grad_b / grad_w via the
 * transposed recipe; see irrep_tp_apply_backward, irrep_tp_apply_uvw_backward.
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/clebsch_gordan.h>
#include <irrep/parity.h>
#include <irrep/quadrature.h>
#include <irrep/spherical_harmonics.h>
#include <irrep/tensor_product.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern void irrep_set_error_(const char *fmt, ...);

#define UNUSED(x) ((void)(x))

/* -------------------------------------------------------------------------- *
 * Descriptor                                                                 *
 * -------------------------------------------------------------------------- */

/* Sparse representation of a single non-zero CG tensor entry after the
 * real-basis change. Packs three small indices + a weight in 16 bytes so
 * four entries fit a cache line. Real-basis CG tensors are ~20 % dense
 * (a tight analogue of the complex-basis m_a + m_b = m_c selection rule),
 * so iterating the sparse list is 4–5× less work than sweeping the dense
 * d_a · d_b · d_c block. */
struct tp_nz_entry {
    short  ima;
    short  imb;
    short  imc;
    short  pad_;
    double cg;
};

struct tp_path {
    int i_a, i_b, i_c;
    int l_a, l_b, l_c;
    int d_a, d_b, d_c;          /* 2l + 1 */
    int mult_u, mult_v, mult_w; /* u = a-channels, v = b-channels,
                                 * w = c-channels; for UUU all equal. */
    int     offset_a, offset_b, offset_c;
    int     weight_offset;  /* into the flat weights buffer */
    double *cg;             /* d_a · d_b · d_c doubles (dense, kept
                             * for the backward-by-index paths) */
    struct tp_nz_entry *nz; /* sparse non-zero entries */
    int                 n_nz;
};

struct tp_descriptor {
    irrep_tp_mode_t mode;
    int             num_paths;
    int             num_weights;
    int             a_dim, b_dim, c_dim;
    struct tp_path *paths;
};

/* -------------------------------------------------------------------------- *
 * Helpers                                                                    *
 * -------------------------------------------------------------------------- */

static int multiset_offset_(const irrep_multiset_t *m, int term_idx) {
    int off = 0;
    for (int i = 0; i < term_idx; ++i) {
        off += m->multiplicities[i] * (2 * m->labels[i].l + 1);
    }
    return off;
}

/* -------------------------------------------------------------------------- *
 * Path enumeration                                                           *
 * -------------------------------------------------------------------------- */

int irrep_tp_enumerate_paths(const irrep_multiset_t *a, const irrep_multiset_t *b,
                             const irrep_multiset_t *c, int *out_paths, int max_paths) {
    if (!a || !b || !c)
        return 0;
    int count = 0;
    for (int ia = 0; ia < a->num_terms; ++ia) {
        int la = a->labels[ia].l;
        int pa = a->labels[ia].parity;
        for (int ib = 0; ib < b->num_terms; ++ib) {
            int lb = b->labels[ib].l;
            int pb = b->labels[ib].parity;
            int pab = pa * pb;
            int lo = la > lb ? la - lb : lb - la;
            int hi = la + lb;
            for (int ic = 0; ic < c->num_terms; ++ic) {
                int lc = c->labels[ic].l;
                int pc = c->labels[ic].parity;
                if (lc < lo || lc > hi)
                    continue;
                if (pab != pc)
                    continue;
                /* Both even and odd l-sum paths are supported — the coupling
                 * tensor is computed as a real-basis Gaunt integral, which is
                 * real by construction for every triangle-valid triple. */
                if (out_paths && count < max_paths) {
                    out_paths[count * 3 + 0] = ia;
                    out_paths[count * 3 + 1] = ib;
                    out_paths[count * 3 + 2] = ic;
                }
                count++;
            }
        }
    }
    return count;
}

/* -------------------------------------------------------------------------- *
 * Build                                                                      *
 * -------------------------------------------------------------------------- */

tp_descriptor_t *irrep_tp_build(const irrep_multiset_t *a, const irrep_multiset_t *b,
                                const irrep_multiset_t *c, const int *selected_paths,
                                int num_selected_paths) {
    if (!a || !b || !c) {
        irrep_set_error_("irrep_tp_build: NULL multiset input");
        return NULL;
    }
    if (num_selected_paths < 0) {
        irrep_set_error_("irrep_tp_build: negative path count");
        return NULL;
    }

    tp_descriptor_t *desc = calloc(1, sizeof(*desc));
    if (!desc)
        return NULL;
    desc->mode = IRREP_TP_MODE_UUU;
    desc->a_dim = a->total_dim;
    desc->b_dim = b->total_dim;
    desc->c_dim = c->total_dim;
    desc->num_paths = num_selected_paths;
    desc->num_weights = num_selected_paths;

    if (num_selected_paths == 0)
        return desc;

    desc->paths = calloc((size_t)num_selected_paths, sizeof(*desc->paths));
    if (!desc->paths) {
        free(desc);
        return NULL;
    }

    for (int k = 0; k < num_selected_paths; ++k) {
        int ia = selected_paths[k * 3 + 0];
        int ib = selected_paths[k * 3 + 1];
        int ic = selected_paths[k * 3 + 2];
        if (ia < 0 || ia >= a->num_terms || ib < 0 || ib >= b->num_terms || ic < 0 ||
            ic >= c->num_terms) {
            irrep_set_error_("irrep_tp_build: path %d has out-of-range index", k);
            goto fail;
        }
        int la = a->labels[ia].l, pa = a->labels[ia].parity;
        int lb = b->labels[ib].l, pb = b->labels[ib].parity;
        int lc = c->labels[ic].l, pc = c->labels[ic].parity;
        int lo = la > lb ? la - lb : lb - la;
        int hi = la + lb;
        if (lc < lo || lc > hi) {
            irrep_set_error_("irrep_tp_build: path %d violates triangle (|%d-%d|<=%d<=%d+%d)", k,
                             la, lb, lc, la, lb);
            goto fail;
        }
        if (pa * pb != pc) {
            irrep_set_error_("irrep_tp_build: path %d violates parity (%d * %d = %d, got %d)", k,
                             pa, pb, pa * pb, pc);
            goto fail;
        }
        int mult_a = a->multiplicities[ia];
        int mult_b = b->multiplicities[ib];
        int mult_c = c->multiplicities[ic];
        if (mult_a != mult_b || mult_a != mult_c) {
            irrep_set_error_("irrep_tp_build: path %d needs matching multiplicities "
                             "(got %d, %d, %d)",
                             k, mult_a, mult_b, mult_c);
            goto fail;
        }

        struct tp_path *p = &desc->paths[k];
        p->i_a = ia;
        p->i_b = ib;
        p->i_c = ic;
        p->l_a = la;
        p->l_b = lb;
        p->l_c = lc;
        p->d_a = 2 * la + 1;
        p->d_b = 2 * lb + 1;
        p->d_c = 2 * lc + 1;
        p->mult_u = p->mult_v = p->mult_w = mult_a;
        p->weight_offset = k; /* 1 weight per path */
        p->offset_a = multiset_offset_(a, ia);
        p->offset_b = multiset_offset_(b, ib);
        p->offset_c = multiset_offset_(c, ic);

        /* Real-basis coupling tensor W via complex CG + change-of-basis,
         * with the e3nn "i^{l_a + l_b − l_c}" phase that keeps W real for
         * every triangle-valid (l_a, l_b, l_c) — including the antisymmetric
         * odd-l-sum couplings (e.g. 1 ⊗ 1 → 1, the cross product):
         *
         *   Q_l[n, m] = i^l · U_l[n, m]     (complex-to-real with extra phase)
         *   W_C[n_a, n_b, n_c]
         *     = Σ_{m_a, m_b, m_c}
         *         conj(Q_c[n_c, m_c]) · CG[m_a, m_b, m_c]
         *           · Q_a[n_a, m_a] · Q_b[n_b, m_b]
         *     = i^{l_a + l_b − l_c} · [old U-only formula]
         *
         * The result is purely real; we store its real part directly.
         *
         * References:
         *   • Limpanuparb & Milthorpe (2014), arXiv:1410.1748, §3.
         *   • e3nn `o3.wigner_3j`: the same convention, repackaged as a 3j. */
        size_t           nW = (size_t)p->d_a * (size_t)p->d_b * (size_t)p->d_c;
        double          *cg_c = calloc(nW, sizeof(double));
        double _Complex *U_a = calloc((size_t)(p->d_a * p->d_a), sizeof(double _Complex));
        double _Complex *U_b = calloc((size_t)(p->d_b * p->d_b), sizeof(double _Complex));
        double _Complex *U_c = calloc((size_t)(p->d_c * p->d_c), sizeof(double _Complex));
        p->cg = calloc(nW, sizeof(double));
        if (!cg_c || !U_a || !U_b || !U_c || !p->cg) {
            free(cg_c);
            free(U_a);
            free(U_b);
            free(U_c);
            goto fail;
        }
        irrep_cg_block(la, lb, lc, cg_c);
        irrep_sph_harm_complex_to_real(la, U_a);
        irrep_sph_harm_complex_to_real(lb, U_b);
        irrep_sph_harm_complex_to_real(lc, U_c);

        /* Phase exponent i^{l_a + l_b − l_c} ∈ {1, i, −1, −i}. */
        int phase_mod = (((la + lb - lc) % 4) + 4) % 4;
        double _Complex phase;
        switch (phase_mod) {
        case 0:
            phase = 1.0;
            break;
        case 1:
            phase = 1.0 * I;
            break;
        case 2:
            phase = -1.0;
            break;
        default:
            phase = -1.0 * I;
            break;
        }

        for (int na = 0; na < p->d_a; ++na) {
            for (int nb = 0; nb < p->d_b; ++nb) {
                for (int nc = 0; nc < p->d_c; ++nc) {
                    double _Complex s = 0.0;
                    for (int ma = 0; ma < p->d_a; ++ma) {
                        double _Complex ua = U_a[na * p->d_a + ma];
                        if (ua == 0.0)
                            continue;
                        for (int mb = 0; mb < p->d_b; ++mb) {
                            double _Complex ub = U_b[nb * p->d_b + mb];
                            if (ub == 0.0)
                                continue;
                            for (int mc = 0; mc < p->d_c; ++mc) {
                                double cg = cg_c[(ma * p->d_b + mb) * p->d_c + mc];
                                if (cg == 0.0)
                                    continue;
                                s += conj(U_c[nc * p->d_c + mc]) * cg * ua * ub;
                            }
                        }
                    }
                    p->cg[(na * p->d_b + nb) * p->d_c + nc] = creal(phase * s);
                }
            }
        }
        free(cg_c);
        free(U_a);
        free(U_b);
        free(U_c);

        /* Build the sparse non-zero index for the apply/backward hot
         * loops. Visit (ima, imb, imc) in the same order apply_core_
         * does so iteration is memory-sequential. */
        int n_nz = 0;
        for (size_t i = 0; i < nW; ++i)
            if (p->cg[i] != 0.0)
                ++n_nz;
        p->nz = (n_nz > 0) ? malloc((size_t)n_nz * sizeof *p->nz) : NULL;
        p->n_nz = n_nz;
        if (n_nz > 0 && !p->nz)
            goto fail;
        int e = 0;
        for (int ima = 0; ima < p->d_a; ++ima) {
            for (int imb = 0; imb < p->d_b; ++imb) {
                for (int imc = 0; imc < p->d_c; ++imc) {
                    double v = p->cg[(ima * p->d_b + imb) * p->d_c + imc];
                    if (v == 0.0)
                        continue;
                    p->nz[e].ima = (short)ima;
                    p->nz[e].imb = (short)imb;
                    p->nz[e].imc = (short)imc;
                    p->nz[e].pad_ = 0;
                    p->nz[e].cg = v;
                    ++e;
                }
            }
        }
    }
    return desc;

fail:
    irrep_tp_free(desc);
    return NULL;
}

void irrep_tp_free(tp_descriptor_t *desc) {
    if (!desc)
        return;
    if (desc->paths) {
        for (int k = 0; k < desc->num_paths; ++k) {
            free(desc->paths[k].cg);
            free(desc->paths[k].nz);
        }
        free(desc->paths);
    }
    free(desc);
}

int irrep_tp_output_dim(const tp_descriptor_t *d) {
    return d ? d->c_dim : 0;
}
int irrep_tp_num_paths(const tp_descriptor_t *d) {
    return d ? d->num_paths : 0;
}

irrep_tp_mode_t irrep_tp_mode(const tp_descriptor_t *d) {
    return d ? d->mode : IRREP_TP_MODE_UUU;
}

/* -------------------------------------------------------------------------- *
 * Forward (unweighted / weighted)                                            *
 *                                                                            *
 * NOTE: apply zeroes c_out before accumulating. Backward accumulates into    *
 * caller-provided grad buffers (caller pre-zeros).                           *
 * -------------------------------------------------------------------------- */

static void apply_core_(const tp_descriptor_t *desc, const double *weights, /* may be NULL */
                        const double *a_in, const double *b_in, double *c_out) {
    for (int i = 0; i < desc->c_dim; ++i)
        c_out[i] = 0.0;

    for (int k = 0; k < desc->num_paths; ++k) {
        const struct tp_path     *p = &desc->paths[k];
        double                    w = weights ? weights[k] : 1.0;
        int                       d_a = p->d_a, d_b = p->d_b, d_c = p->d_c;
        const struct tp_nz_entry *nz = p->nz;
        int                       n_nz = p->n_nz;

        for (int u = 0; u < p->mult_u; ++u) {
            const double *a_block = a_in + p->offset_a + u * d_a;
            const double *b_block = b_in + p->offset_b + u * d_b;
            /* */ double *c_block = c_out + p->offset_c + u * d_c;

            /* Sparse iteration: visit only the ~20 % of (ima, imb, imc)
             * triples with non-zero CG. Each entry is a straight-line
             * FMA with no inner-loop branch. */
            for (int e = 0; e < n_nz; ++e) {
                c_block[nz[e].imc] += w * nz[e].cg * a_block[nz[e].ima] * b_block[nz[e].imb];
            }
        }
    }
}

void irrep_tp_apply(const tp_descriptor_t *desc, const double *a_in, const double *b_in,
                    double *c_out) {
    if (!desc)
        return;
    apply_core_(desc, NULL, a_in, b_in, c_out);
}

void irrep_tp_apply_weighted(const tp_descriptor_t *desc, const double *weights, const double *a_in,
                             const double *b_in, double *c_out) {
    if (!desc)
        return;
    apply_core_(desc, weights, a_in, b_in, c_out);
}

/* -------------------------------------------------------------------------- *
 * Backward                                                                   *
 * -------------------------------------------------------------------------- */

static void backward_core_(const tp_descriptor_t *desc, const double *weights, /* NULL → all ones */
                           const double *a_in, const double *b_in, const double *grad_c,
                           double *grad_a, double *grad_b,
                           double *grad_w) /* NULL if not wanted */ {
    for (int k = 0; k < desc->num_paths; ++k) {
        const struct tp_path     *p = &desc->paths[k];
        double                    w = weights ? weights[k] : 1.0;
        int                       d_a = p->d_a, d_b = p->d_b, d_c = p->d_c;
        const struct tp_nz_entry *nz = p->nz;
        int                       n_nz = p->n_nz;
        double                    dw = 0.0;

        for (int u = 0; u < p->mult_u; ++u) {
            const double *a_block = a_in + p->offset_a + u * d_a;
            const double *b_block = b_in + p->offset_b + u * d_b;
            const double *gc_block = grad_c + p->offset_c + u * d_c;
            double       *ga_block = grad_a ? grad_a + p->offset_a + u * d_a : NULL;
            double       *gb_block = grad_b ? grad_b + p->offset_b + u * d_b : NULL;

            /* Sparse backward: iterate only the non-zero (ima, imb, imc)
             * entries. Each entry contributes to three gradients. */
            for (int e = 0; e < n_nz; ++e) {
                double cg = nz[e].cg;
                int    ima = nz[e].ima, imb = nz[e].imb, imc = nz[e].imc;
                double gc = gc_block[imc];
                if (gc == 0.0)
                    continue;
                double a_v = a_block[ima];
                double b_v = b_block[imb];
                if (ga_block)
                    ga_block[ima] += w * cg * b_v * gc;
                if (gb_block)
                    gb_block[imb] += w * cg * a_v * gc;
                if (grad_w)
                    dw += cg * a_v * b_v * gc;
            }
        }
        if (grad_w)
            grad_w[k] += dw;
    }
}

void irrep_tp_apply_backward(const tp_descriptor_t *desc, const double *a_in, const double *b_in,
                             const double *grad_c_out, double *grad_a, double *grad_b) {
    if (!desc)
        return;
    backward_core_(desc, NULL, a_in, b_in, grad_c_out, grad_a, grad_b, NULL);
}

void irrep_tp_apply_backward_weighted(const tp_descriptor_t *desc, const double *weights,
                                      const double *a_in, const double *b_in,
                                      const double *grad_c_out, double *grad_a, double *grad_b,
                                      double *grad_w) {
    if (!desc)
        return;
    backward_core_(desc, weights, a_in, b_in, grad_c_out, grad_a, grad_b, grad_w);
}

/* -------------------------------------------------------------------------- *
 * Batched variants                                                           *
 * -------------------------------------------------------------------------- */

void irrep_tp_apply_batch(const tp_descriptor_t *desc, size_t batch, const double *a_in,
                          const double *b_in, double *c_out) {
    if (!desc)
        return;
    for (size_t bi = 0; bi < batch; ++bi) {
        apply_core_(desc, NULL, a_in + bi * (size_t)desc->a_dim, b_in + bi * (size_t)desc->b_dim,
                    c_out + bi * (size_t)desc->c_dim);
    }
}

void irrep_tp_apply_weighted_batch(const tp_descriptor_t *desc, size_t batch, const double *weights,
                                   const double *a_in, const double *b_in, double *c_out) {
    if (!desc)
        return;
    for (size_t bi = 0; bi < batch; ++bi) {
        apply_core_(desc, weights, a_in + bi * (size_t)desc->a_dim, b_in + bi * (size_t)desc->b_dim,
                    c_out + bi * (size_t)desc->c_dim);
    }
}

void irrep_tp_apply_backward_batch(const tp_descriptor_t *desc, size_t batch, const double *a_in,
                                   const double *b_in, const double *grad_c_out, double *grad_a,
                                   double *grad_b) {
    if (!desc)
        return;
    for (size_t bi = 0; bi < batch; ++bi) {
        backward_core_(desc, NULL, a_in + bi * (size_t)desc->a_dim, b_in + bi * (size_t)desc->b_dim,
                       grad_c_out + bi * (size_t)desc->c_dim,
                       grad_a ? grad_a + bi * (size_t)desc->a_dim : NULL,
                       grad_b ? grad_b + bi * (size_t)desc->b_dim : NULL, NULL);
    }
}

void irrep_tp_apply_backward_weighted_batch(const tp_descriptor_t *desc, size_t batch,
                                            const double *weights, const double *a_in,
                                            const double *b_in, const double *grad_c_out,
                                            double *grad_a, double *grad_b, double *grad_w) {
    if (!desc)
        return;
    for (size_t bi = 0; bi < batch; ++bi) {
        backward_core_(desc, weights, a_in + bi * (size_t)desc->a_dim,
                       b_in + bi * (size_t)desc->b_dim, grad_c_out + bi * (size_t)desc->c_dim,
                       grad_a ? grad_a + bi * (size_t)desc->a_dim : NULL,
                       grad_b ? grad_b + bi * (size_t)desc->b_dim : NULL,
                       grad_w); /* grad_w accumulates across batch */
    }
}

/* ========================================================================== *
 * UVW mode: independent multiplicities, per-(path, w, v, u) weights.         *
 * ========================================================================== */

tp_descriptor_t *irrep_tp_build_uvw(const irrep_multiset_t *a, const irrep_multiset_t *b,
                                    const irrep_multiset_t *c, const int *selected_paths,
                                    int num_selected_paths) {
    if (!a || !b || !c) {
        irrep_set_error_("irrep_tp_build_uvw: NULL multiset input");
        return NULL;
    }
    if (num_selected_paths < 0) {
        irrep_set_error_("irrep_tp_build_uvw: negative path count");
        return NULL;
    }

    tp_descriptor_t *desc = calloc(1, sizeof(*desc));
    if (!desc)
        return NULL;
    desc->mode = IRREP_TP_MODE_UVW;
    desc->a_dim = a->total_dim;
    desc->b_dim = b->total_dim;
    desc->c_dim = c->total_dim;
    desc->num_paths = num_selected_paths;
    desc->num_weights = 0;

    if (num_selected_paths == 0)
        return desc;
    desc->paths = calloc((size_t)num_selected_paths, sizeof(*desc->paths));
    if (!desc->paths) {
        free(desc);
        return NULL;
    }

    int wsum = 0;
    for (int k = 0; k < num_selected_paths; ++k) {
        int ia = selected_paths[k * 3 + 0];
        int ib = selected_paths[k * 3 + 1];
        int ic = selected_paths[k * 3 + 2];
        if (ia < 0 || ia >= a->num_terms || ib < 0 || ib >= b->num_terms || ic < 0 ||
            ic >= c->num_terms) {
            irrep_set_error_("irrep_tp_build_uvw: path %d out-of-range index", k);
            goto fail;
        }
        int la = a->labels[ia].l, pa = a->labels[ia].parity;
        int lb = b->labels[ib].l, pb = b->labels[ib].parity;
        int lc = c->labels[ic].l, pc = c->labels[ic].parity;
        int lo = la > lb ? la - lb : lb - la;
        int hi = la + lb;
        if (lc < lo || lc > hi) {
            irrep_set_error_("irrep_tp_build_uvw: path %d violates triangle", k);
            goto fail;
        }
        if (pa * pb != pc) {
            irrep_set_error_("irrep_tp_build_uvw: path %d violates parity", k);
            goto fail;
        }

        struct tp_path *p = &desc->paths[k];
        p->i_a = ia;
        p->i_b = ib;
        p->i_c = ic;
        p->l_a = la;
        p->l_b = lb;
        p->l_c = lc;
        p->d_a = 2 * la + 1;
        p->d_b = 2 * lb + 1;
        p->d_c = 2 * lc + 1;
        p->mult_u = a->multiplicities[ia];
        p->mult_v = b->multiplicities[ib];
        p->mult_w = c->multiplicities[ic];
        p->offset_a = multiset_offset_(a, ia);
        p->offset_b = multiset_offset_(b, ib);
        p->offset_c = multiset_offset_(c, ic);
        p->weight_offset = wsum;
        wsum += p->mult_w * p->mult_v * p->mult_u;

        /* Reuse the same Gaunt / CG+phase coupling tensor as UUU build. */
        size_t           nW = (size_t)p->d_a * (size_t)p->d_b * (size_t)p->d_c;
        double          *cg_c = calloc(nW, sizeof(double));
        double _Complex *U_a = calloc((size_t)(p->d_a * p->d_a), sizeof(double _Complex));
        double _Complex *U_b = calloc((size_t)(p->d_b * p->d_b), sizeof(double _Complex));
        double _Complex *U_c = calloc((size_t)(p->d_c * p->d_c), sizeof(double _Complex));
        p->cg = calloc(nW, sizeof(double));
        if (!cg_c || !U_a || !U_b || !U_c || !p->cg) {
            free(cg_c);
            free(U_a);
            free(U_b);
            free(U_c);
            goto fail;
        }
        irrep_cg_block(la, lb, lc, cg_c);
        irrep_sph_harm_complex_to_real(la, U_a);
        irrep_sph_harm_complex_to_real(lb, U_b);
        irrep_sph_harm_complex_to_real(lc, U_c);

        int phase_mod = (((la + lb - lc) % 4) + 4) % 4;
        double _Complex phase;
        switch (phase_mod) {
        case 0:
            phase = 1.0;
            break;
        case 1:
            phase = 1.0 * I;
            break;
        case 2:
            phase = -1.0;
            break;
        default:
            phase = -1.0 * I;
            break;
        }

        for (int na = 0; na < p->d_a; ++na)
            for (int nb = 0; nb < p->d_b; ++nb)
                for (int nc = 0; nc < p->d_c; ++nc) {
                    double _Complex s = 0.0;
                    for (int ma = 0; ma < p->d_a; ++ma) {
                        double _Complex ua = U_a[na * p->d_a + ma];
                        if (ua == 0.0)
                            continue;
                        for (int mb = 0; mb < p->d_b; ++mb) {
                            double _Complex ub = U_b[nb * p->d_b + mb];
                            if (ub == 0.0)
                                continue;
                            for (int mc = 0; mc < p->d_c; ++mc) {
                                double cg = cg_c[(ma * p->d_b + mb) * p->d_c + mc];
                                if (cg == 0.0)
                                    continue;
                                s += conj(U_c[nc * p->d_c + mc]) * cg * ua * ub;
                            }
                        }
                    }
                    p->cg[(na * p->d_b + nb) * p->d_c + nc] = creal(phase * s);
                }
        free(cg_c);
        free(U_a);
        free(U_b);
        free(U_c);

        /* Sparse non-zero list (same structure as UUU path). */
        int n_nz_u = 0;
        for (size_t i = 0; i < nW; ++i)
            if (p->cg[i] != 0.0)
                ++n_nz_u;
        p->nz = (n_nz_u > 0) ? malloc((size_t)n_nz_u * sizeof *p->nz) : NULL;
        p->n_nz = n_nz_u;
        if (n_nz_u > 0 && !p->nz)
            goto fail;
        int e2 = 0;
        for (int ima = 0; ima < p->d_a; ++ima) {
            for (int imb = 0; imb < p->d_b; ++imb) {
                for (int imc = 0; imc < p->d_c; ++imc) {
                    double v = p->cg[(ima * p->d_b + imb) * p->d_c + imc];
                    if (v == 0.0)
                        continue;
                    p->nz[e2].ima = (short)ima;
                    p->nz[e2].imb = (short)imb;
                    p->nz[e2].imc = (short)imc;
                    p->nz[e2].pad_ = 0;
                    p->nz[e2].cg = v;
                    ++e2;
                }
            }
        }
    }
    desc->num_weights = wsum;
    return desc;

fail:
    irrep_tp_free(desc);
    return NULL;
}

int irrep_tp_num_weights_uvw(const tp_descriptor_t *d) {
    return (d && d->mode == IRREP_TP_MODE_UVW) ? d->num_weights : 0;
}

void irrep_tp_apply_uvw(const tp_descriptor_t *desc, const double *weights, const double *a_in,
                        const double *b_in, double *c_out) {
    if (!desc || desc->mode != IRREP_TP_MODE_UVW)
        return;
    for (int i = 0; i < desc->c_dim; ++i)
        c_out[i] = 0.0;

    for (int k = 0; k < desc->num_paths; ++k) {
        const struct tp_path *p = &desc->paths[k];
        int                   d_a = p->d_a, d_b = p->d_b, d_c = p->d_c;
        int                   U = p->mult_u, V = p->mult_v, W = p->mult_w;
        const double         *wp = weights + p->weight_offset; /* [W, V, U] */

        for (int w = 0; w < W; ++w) {
            double *c_block = c_out + p->offset_c + w * d_c;
            for (int v = 0; v < V; ++v) {
                const double *b_block = b_in + p->offset_b + v * d_b;
                for (int u = 0; u < U; ++u) {
                    double wt = wp[(w * V + v) * U + u];
                    if (wt == 0.0)
                        continue;
                    const double *a_block = a_in + p->offset_a + u * d_a;
                    for (int imc = 0; imc < d_c; ++imc) {
                        double s = 0.0;
                        for (int ima = 0; ima < d_a; ++ima) {
                            double a_v = a_block[ima];
                            if (a_v == 0.0)
                                continue;
                            for (int imb = 0; imb < d_b; ++imb) {
                                double cg = p->cg[(ima * d_b + imb) * d_c + imc];
                                if (cg == 0.0)
                                    continue;
                                s += cg * a_v * b_block[imb];
                            }
                        }
                        c_block[imc] += wt * s;
                    }
                }
            }
        }
    }
}

void irrep_tp_apply_uvw_backward(const tp_descriptor_t *desc, const double *weights,
                                 const double *a_in, const double *b_in, const double *grad_c_out,
                                 double *grad_a, double *grad_b, double *grad_w) {
    if (!desc || desc->mode != IRREP_TP_MODE_UVW)
        return;
    for (int k = 0; k < desc->num_paths; ++k) {
        const struct tp_path *p = &desc->paths[k];
        int                   d_a = p->d_a, d_b = p->d_b, d_c = p->d_c;
        int                   U = p->mult_u, V = p->mult_v, W = p->mult_w;
        const double         *wp = weights + p->weight_offset;
        double               *gwp = grad_w ? grad_w + p->weight_offset : NULL;

        for (int w = 0; w < W; ++w) {
            const double *gc = grad_c_out + p->offset_c + w * d_c;
            for (int v = 0; v < V; ++v) {
                const double *b_block = b_in + p->offset_b + v * d_b;
                double       *gb_block = grad_b ? grad_b + p->offset_b + v * d_b : NULL;
                for (int u = 0; u < U; ++u) {
                    const double *a_block = a_in + p->offset_a + u * d_a;
                    double       *ga_block = grad_a ? grad_a + p->offset_a + u * d_a : NULL;
                    double        wt = wp[(w * V + v) * U + u];
                    double        dw = 0.0;

                    for (int imc = 0; imc < d_c; ++imc) {
                        double gv = gc[imc];
                        if (gv == 0.0)
                            continue;
                        for (int ima = 0; ima < d_a; ++ima) {
                            double a_v = a_block[ima];
                            for (int imb = 0; imb < d_b; ++imb) {
                                double cg = p->cg[(ima * d_b + imb) * d_c + imc];
                                if (cg == 0.0)
                                    continue;
                                double b_v = b_block[imb];
                                if (ga_block)
                                    ga_block[ima] += wt * cg * b_v * gv;
                                if (gb_block)
                                    gb_block[imb] += wt * cg * a_v * gv;
                                if (gwp)
                                    dw += cg * a_v * b_v * gv;
                            }
                        }
                    }
                    if (gwp)
                        gwp[(w * V + v) * U + u] += dw;
                }
            }
        }
    }
}

/* -------------------------------------------------------------------------- *
 * Per-path weight L2 — useful for SR-style training where different paths    *
 * on the variational manifold have different curvatures and need per-path    *
 * regularisation strengths.                                                  *
 * -------------------------------------------------------------------------- */

void irrep_tp_weight_l2_per_path_uvw(const tp_descriptor_t *desc, const double *weights,
                                     double *out) {
    if (!desc || desc->mode != IRREP_TP_MODE_UVW || !weights || !out)
        return;
    for (int k = 0; k < desc->num_paths; ++k) {
        const struct tp_path *p = &desc->paths[k];
        int                   U = p->mult_u, V = p->mult_v, W = p->mult_w;
        const double         *wp = weights + p->weight_offset;
        double                s = 0.0;
        int                   n = W * V * U;
        for (int i = 0; i < n; ++i)
            s += wp[i] * wp[i];
        out[k] = s;
    }
}

void irrep_tp_weight_l2_per_path_uvw_backward(const tp_descriptor_t *desc, const double *weights,
                                              const double *grad_per_path, double *grad_weights) {
    if (!desc || desc->mode != IRREP_TP_MODE_UVW || !weights || !grad_per_path || !grad_weights)
        return;
    for (int k = 0; k < desc->num_paths; ++k) {
        const struct tp_path *p = &desc->paths[k];
        int                   U = p->mult_u, V = p->mult_v, W = p->mult_w;
        const double         *wp = weights + p->weight_offset;
        double               *gwp = grad_weights + p->weight_offset;
        double                coeff = 2.0 * grad_per_path[k];
        int                   n = W * V * U;
        for (int i = 0; i < n; ++i)
            gwp[i] += coeff * wp[i];
    }
}
