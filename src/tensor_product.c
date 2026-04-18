/* SPDX-License-Identifier: MIT */
/* M7: e3nn-style path-indexed tensor products.
 *
 * Channel mode: "uuu". For each selected path (i_a, i_b, i_c) we require
 * mult_a = mult_b = mult_c, and pair channels one-to-one across the three
 * multisets. One scalar weight per path.
 *
 * For each path p with copy index u ∈ [0, mult_p):
 *
 *     c[i_c, u, m_c] += w_p · Σ_{m_a, m_b}
 *         ⟨l_a m_a; l_b m_b | l_c m_c⟩ · a[i_a, u, m_a] · b[i_b, u, m_b]
 *
 * Build precomputes the CG block once; apply is a loop over (m_c, m_a, m_b).
 *
 * Backward:
 *     ∂L/∂a[i_a, u, m_a] += w_p · Σ_{m_b, m_c} CG · b[i_b, u, m_b] · ∂L/∂c[i_c, u, m_c]
 *     ∂L/∂b[i_b, u, m_b] += w_p · Σ_{m_a, m_c} CG · a[i_a, u, m_a] · ∂L/∂c[i_c, u, m_c]
 *     ∂L/∂w_p            += Σ_u Σ_{m_a,m_b,m_c} CG · a · b · ∂L/∂c
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/clebsch_gordan.h>
#include <irrep/parity.h>
#include <irrep/spherical_harmonics.h>
#include <irrep/tensor_product.h>

extern void irrep_set_error_(const char *fmt, ...);

#define UNUSED(x) ((void)(x))

/* -------------------------------------------------------------------------- *
 * Descriptor                                                                 *
 * -------------------------------------------------------------------------- */

struct tp_path {
    int     i_a, i_b, i_c;
    int     l_a, l_b, l_c;
    int     d_a, d_b, d_c;        /* 2l + 1 */
    int     mult;                  /* mult_a = mult_b = mult_c */
    int     offset_a, offset_b, offset_c;
    double *cg;                    /* d_a · d_b · d_c doubles */
};

struct tp_descriptor {
    int             num_paths;
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

int irrep_tp_enumerate_paths(const irrep_multiset_t *a,
                             const irrep_multiset_t *b,
                             const irrep_multiset_t *c,
                             int *out_paths, int max_paths) {
    if (!a || !b || !c) return 0;
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
                if (lc < lo || lc > hi) continue;
                if (pab != pc)          continue;
                /* v1.0 supports "parity-natural" couplings only: the real-basis
                 * coupling tensor is real iff l_a + l_b + l_c is even. Odd-sum
                 * paths are antisymmetric in the complex basis and would need
                 * an explicit i-factor convention (planned for v1.1). */
                if (((la + lb + lc) & 1) != 0) continue;
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

tp_descriptor_t *irrep_tp_build(const irrep_multiset_t *a,
                                const irrep_multiset_t *b,
                                const irrep_multiset_t *c,
                                const int *selected_paths,
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
    if (!desc) return NULL;
    desc->a_dim     = a->total_dim;
    desc->b_dim     = b->total_dim;
    desc->c_dim     = c->total_dim;
    desc->num_paths = num_selected_paths;

    if (num_selected_paths == 0) return desc;

    desc->paths = calloc((size_t)num_selected_paths, sizeof(*desc->paths));
    if (!desc->paths) { free(desc); return NULL; }

    for (int k = 0; k < num_selected_paths; ++k) {
        int ia = selected_paths[k * 3 + 0];
        int ib = selected_paths[k * 3 + 1];
        int ic = selected_paths[k * 3 + 2];
        if (ia < 0 || ia >= a->num_terms ||
            ib < 0 || ib >= b->num_terms ||
            ic < 0 || ic >= c->num_terms) {
            irrep_set_error_("irrep_tp_build: path %d has out-of-range index", k);
            goto fail;
        }
        int la = a->labels[ia].l, pa = a->labels[ia].parity;
        int lb = b->labels[ib].l, pb = b->labels[ib].parity;
        int lc = c->labels[ic].l, pc = c->labels[ic].parity;
        int lo = la > lb ? la - lb : lb - la;
        int hi = la + lb;
        if (lc < lo || lc > hi) {
            irrep_set_error_("irrep_tp_build: path %d violates triangle (|%d-%d|<=%d<=%d+%d)",
                             k, la, lb, lc, la, lb);
            goto fail;
        }
        if (pa * pb != pc) {
            irrep_set_error_("irrep_tp_build: path %d violates parity (%d * %d = %d, got %d)",
                             k, pa, pb, pa * pb, pc);
            goto fail;
        }
        if (((la + lb + lc) & 1) != 0) {
            irrep_set_error_("irrep_tp_build: path %d has odd l-sum (%d + %d + %d); "
                             "the real-basis coupling is imaginary, not yet supported",
                             k, la, lb, lc);
            goto fail;
        }
        int mult_a = a->multiplicities[ia];
        int mult_b = b->multiplicities[ib];
        int mult_c = c->multiplicities[ic];
        if (mult_a != mult_b || mult_a != mult_c) {
            irrep_set_error_("irrep_tp_build: path %d needs matching multiplicities "
                             "(got %d, %d, %d)", k, mult_a, mult_b, mult_c);
            goto fail;
        }

        struct tp_path *p = &desc->paths[k];
        p->i_a = ia; p->i_b = ib; p->i_c = ic;
        p->l_a = la; p->l_b = lb; p->l_c = lc;
        p->d_a = 2 * la + 1; p->d_b = 2 * lb + 1; p->d_c = 2 * lc + 1;
        p->mult = mult_a;
        p->offset_a = multiset_offset_(a, ia);
        p->offset_b = multiset_offset_(b, ib);
        p->offset_c = multiset_offset_(c, ic);

        /* Compute the real-basis coupling tensor
         *   W[n_a, n_b, n_c] = Σ_{m_a, m_b, m_c}
         *       U_c[n_c, m_c] · CG[m_a, m_b, m_c]
         *         · conj(U_a[n_a, m_a]) · conj(U_b[n_b, m_b])
         * where U_l is the complex-to-real SH basis change. The sum is
         * identically real when l_a + l_b + l_c is even (the only case
         * allowed by the parity selection plus triangle + parity above),
         * so we drop the tiny imaginary residual that remains from
         * floating-point round-off. */
        size_t nW = (size_t)p->d_a * (size_t)p->d_b * (size_t)p->d_c;
        double          *cg_c = calloc(nW, sizeof(double));
        double _Complex *U_a  = calloc((size_t)(p->d_a * p->d_a), sizeof(double _Complex));
        double _Complex *U_b  = calloc((size_t)(p->d_b * p->d_b), sizeof(double _Complex));
        double _Complex *U_c  = calloc((size_t)(p->d_c * p->d_c), sizeof(double _Complex));
        p->cg = calloc(nW, sizeof(double));
        if (!cg_c || !U_a || !U_b || !U_c || !p->cg) {
            free(cg_c); free(U_a); free(U_b); free(U_c);
            goto fail;
        }
        irrep_cg_block(la, lb, lc, cg_c);
        irrep_sph_harm_complex_to_real(la, U_a);
        irrep_sph_harm_complex_to_real(lb, U_b);
        irrep_sph_harm_complex_to_real(lc, U_c);

        for (int na = 0; na < p->d_a; ++na) {
            for (int nb = 0; nb < p->d_b; ++nb) {
                for (int nc = 0; nc < p->d_c; ++nc) {
                    double _Complex s = 0.0;
                    for (int ma = 0; ma < p->d_a; ++ma) {
                        double _Complex ua_star = conj(U_a[na * p->d_a + ma]);
                        if (ua_star == 0.0) continue;
                        for (int mb = 0; mb < p->d_b; ++mb) {
                            double _Complex ub_star = conj(U_b[nb * p->d_b + mb]);
                            if (ub_star == 0.0) continue;
                            for (int mc = 0; mc < p->d_c; ++mc) {
                                double cg = cg_c[(ma * p->d_b + mb) * p->d_c + mc];
                                if (cg == 0.0) continue;
                                s += U_c[nc * p->d_c + mc] * cg * ua_star * ub_star;
                            }
                        }
                    }
                    p->cg[(na * p->d_b + nb) * p->d_c + nc] = creal(s);
                }
            }
        }
        free(cg_c); free(U_a); free(U_b); free(U_c);
    }
    return desc;

fail:
    irrep_tp_free(desc);
    return NULL;
}

void irrep_tp_free(tp_descriptor_t *desc) {
    if (!desc) return;
    if (desc->paths) {
        for (int k = 0; k < desc->num_paths; ++k) free(desc->paths[k].cg);
        free(desc->paths);
    }
    free(desc);
}

int irrep_tp_output_dim(const tp_descriptor_t *d) { return d ? d->c_dim    : 0; }
int irrep_tp_num_paths (const tp_descriptor_t *d) { return d ? d->num_paths : 0; }

/* -------------------------------------------------------------------------- *
 * Forward (unweighted / weighted)                                            *
 *                                                                            *
 * NOTE: apply zeroes c_out before accumulating. Backward accumulates into    *
 * caller-provided grad buffers (caller pre-zeros).                           *
 * -------------------------------------------------------------------------- */

static void apply_core_(const tp_descriptor_t *desc,
                        const double *weights, /* may be NULL */
                        const double *a_in,
                        const double *b_in,
                        double       *c_out) {
    for (int i = 0; i < desc->c_dim; ++i) c_out[i] = 0.0;

    for (int k = 0; k < desc->num_paths; ++k) {
        const struct tp_path *p = &desc->paths[k];
        double w = weights ? weights[k] : 1.0;
        int d_a = p->d_a, d_b = p->d_b, d_c = p->d_c;

        for (int u = 0; u < p->mult; ++u) {
            const double *a_block = a_in + p->offset_a + u * d_a;
            const double *b_block = b_in + p->offset_b + u * d_b;
            /* */  double *c_block = c_out + p->offset_c + u * d_c;

            for (int imc = 0; imc < d_c; ++imc) {
                double s = 0.0;
                for (int ima = 0; ima < d_a; ++ima) {
                    double a_v = a_block[ima];
                    if (a_v == 0.0) continue;
                    for (int imb = 0; imb < d_b; ++imb) {
                        double cg = p->cg[(ima * d_b + imb) * d_c + imc];
                        if (cg == 0.0) continue;
                        s += cg * a_v * b_block[imb];
                    }
                }
                c_block[imc] += w * s;
            }
        }
    }
}

void irrep_tp_apply(const tp_descriptor_t *desc,
                    const double *a_in, const double *b_in, double *c_out) {
    if (!desc) return;
    apply_core_(desc, NULL, a_in, b_in, c_out);
}

void irrep_tp_apply_weighted(const tp_descriptor_t *desc,
                             const double *weights,
                             const double *a_in, const double *b_in,
                             double       *c_out) {
    if (!desc) return;
    apply_core_(desc, weights, a_in, b_in, c_out);
}

/* -------------------------------------------------------------------------- *
 * Backward                                                                   *
 * -------------------------------------------------------------------------- */

static void backward_core_(const tp_descriptor_t *desc,
                           const double *weights,    /* NULL → all ones */
                           const double *a_in,
                           const double *b_in,
                           const double *grad_c,
                           double       *grad_a,
                           double       *grad_b,
                           double       *grad_w)     /* NULL if not wanted */ {
    for (int k = 0; k < desc->num_paths; ++k) {
        const struct tp_path *p = &desc->paths[k];
        double w = weights ? weights[k] : 1.0;
        int d_a = p->d_a, d_b = p->d_b, d_c = p->d_c;
        double dw = 0.0;

        for (int u = 0; u < p->mult; ++u) {
            const double *a_block  = a_in  + p->offset_a + u * d_a;
            const double *b_block  = b_in  + p->offset_b + u * d_b;
            const double *gc_block = grad_c + p->offset_c + u * d_c;
            double       *ga_block = grad_a ? grad_a + p->offset_a + u * d_a : NULL;
            double       *gb_block = grad_b ? grad_b + p->offset_b + u * d_b : NULL;

            for (int imc = 0; imc < d_c; ++imc) {
                double gc = gc_block[imc];
                if (gc == 0.0) continue;
                for (int ima = 0; ima < d_a; ++ima) {
                    double a_v = a_block[ima];
                    for (int imb = 0; imb < d_b; ++imb) {
                        double cg = p->cg[(ima * d_b + imb) * d_c + imc];
                        if (cg == 0.0) continue;
                        double b_v = b_block[imb];
                        if (ga_block) ga_block[ima] += w * cg * b_v * gc;
                        if (gb_block) gb_block[imb] += w * cg * a_v * gc;
                        if (grad_w)   dw            +=     cg * a_v * b_v * gc;
                    }
                }
            }
        }
        if (grad_w) grad_w[k] += dw;
    }
}

void irrep_tp_apply_backward(const tp_descriptor_t *desc,
                             const double *a_in, const double *b_in,
                             const double *grad_c_out,
                             double *grad_a, double *grad_b) {
    if (!desc) return;
    backward_core_(desc, NULL, a_in, b_in, grad_c_out, grad_a, grad_b, NULL);
}

void irrep_tp_apply_backward_weighted(const tp_descriptor_t *desc,
                                      const double *weights,
                                      const double *a_in, const double *b_in,
                                      const double *grad_c_out,
                                      double *grad_a, double *grad_b, double *grad_w) {
    if (!desc) return;
    backward_core_(desc, weights, a_in, b_in, grad_c_out, grad_a, grad_b, grad_w);
}

/* -------------------------------------------------------------------------- *
 * Batched variants                                                           *
 * -------------------------------------------------------------------------- */

void irrep_tp_apply_batch(const tp_descriptor_t *desc, size_t batch,
                          const double *a_in, const double *b_in, double *c_out) {
    if (!desc) return;
    for (size_t bi = 0; bi < batch; ++bi) {
        apply_core_(desc, NULL,
                    a_in  + bi * (size_t)desc->a_dim,
                    b_in  + bi * (size_t)desc->b_dim,
                    c_out + bi * (size_t)desc->c_dim);
    }
}

void irrep_tp_apply_weighted_batch(const tp_descriptor_t *desc, size_t batch,
                                   const double *weights,
                                   const double *a_in, const double *b_in,
                                   double *c_out) {
    if (!desc) return;
    for (size_t bi = 0; bi < batch; ++bi) {
        apply_core_(desc, weights,
                    a_in  + bi * (size_t)desc->a_dim,
                    b_in  + bi * (size_t)desc->b_dim,
                    c_out + bi * (size_t)desc->c_dim);
    }
}

void irrep_tp_apply_backward_batch(const tp_descriptor_t *desc, size_t batch,
                                   const double *a_in, const double *b_in,
                                   const double *grad_c_out,
                                   double *grad_a, double *grad_b) {
    if (!desc) return;
    for (size_t bi = 0; bi < batch; ++bi) {
        backward_core_(desc, NULL,
                       a_in      + bi * (size_t)desc->a_dim,
                       b_in      + bi * (size_t)desc->b_dim,
                       grad_c_out + bi * (size_t)desc->c_dim,
                       grad_a ? grad_a + bi * (size_t)desc->a_dim : NULL,
                       grad_b ? grad_b + bi * (size_t)desc->b_dim : NULL,
                       NULL);
    }
}

void irrep_tp_apply_backward_weighted_batch(const tp_descriptor_t *desc, size_t batch,
                                            const double *weights,
                                            const double *a_in, const double *b_in,
                                            const double *grad_c_out,
                                            double *grad_a, double *grad_b, double *grad_w) {
    if (!desc) return;
    for (size_t bi = 0; bi < batch; ++bi) {
        backward_core_(desc, weights,
                       a_in      + bi * (size_t)desc->a_dim,
                       b_in      + bi * (size_t)desc->b_dim,
                       grad_c_out + bi * (size_t)desc->c_dim,
                       grad_a ? grad_a + bi * (size_t)desc->a_dim : NULL,
                       grad_b ? grad_b + bi * (size_t)desc->b_dim : NULL,
                       grad_w);   /* grad_w accumulates across batch */
    }
}
