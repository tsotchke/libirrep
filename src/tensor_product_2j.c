/* SPDX-License-Identifier: MIT */
/* Half-integer (spinor) tensor-product backend 
 *
 * Operates on `irrep_multiset_2j_t` with **complex** amplitudes, because
 * half-integer irreps carry no natural real basis. Selection rules:
 *   - Triangle on doubled-j: |two_j_a − two_j_b| ≤ two_j_c ≤ two_j_a + two_j_b.
 *   - Parity multiplication: p_a · p_b = p_c.
 *   - two_j_a + two_j_b + two_j_c even (automatic under the triangle).
 *
 * Clebsch-Gordan coefficients come from irrep_cg_2j() (Racah single-sum,
 * log-gamma-stabilised). The apply kernel is a straightforward 3-D sum:
 *
 *   c[i_c, w, m_c] += Σ_{u, v} weight_p[w, v, u] · Σ_{m_a, m_b}
 *                      CG_p[m_a, m_b, m_c] · a[i_a, u, m_a] · b[i_b, v, m_b]
 *
 * with m_a, m_b, m_c running in doubled-integer steps of 2 from −two_j to
 * +two_j (so the indexing is the same as in the integer-l case: storage index
 * = (m + two_j) / 2, dimension = two_j + 1).
 *
 * Backward is the transposed recipe and is straightforward.
 */

#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/clebsch_gordan.h>
#include <irrep/tensor_product.h>
#include <irrep/multiset_2j.h>
#include <irrep/parity.h>

extern void irrep_set_error_(const char *fmt, ...);

struct tp_2j_path {
    int     i_a, i_b, i_c;
    int     two_j_a, two_j_b, two_j_c;
    int     d_a, d_b, d_c;                /* two_j + 1 */
    int     mult_u, mult_v, mult_w;
    int     offset_a, offset_b, offset_c;
    int     weight_offset;                 /* into the per-descriptor weights */
    double *cg;                            /* d_a · d_b · d_c doubles */
};

struct tp_2j_descriptor {
    int                 num_paths;
    int                 num_weights;       /* Σ w·v·u per path */
    int                 a_dim, b_dim, c_dim;
    struct tp_2j_path  *paths;
};

/* -------------------------------------------------------------------------- *
 * Helpers                                                                    *
 * -------------------------------------------------------------------------- */

static int multiset_2j_offset_(const irrep_multiset_2j_t *m, int term_idx) {
    int off = 0;
    for (int i = 0; i < term_idx; ++i) {
        off += m->multiplicities[i] * (m->labels[i].two_j + 1);
    }
    return off;
}

static int triangle_ok_2j_(int two_j_a, int two_j_b, int two_j_c) {
    int lo = (two_j_a > two_j_b) ? (two_j_a - two_j_b) : (two_j_b - two_j_a);
    int hi = two_j_a + two_j_b;
    if (two_j_c < lo || two_j_c > hi) return 0;
    if (((two_j_a + two_j_b + two_j_c) & 1) != 0) return 0;
    return 1;
}

/* -------------------------------------------------------------------------- *
 * Path enumeration                                                           *
 * -------------------------------------------------------------------------- */

int irrep_tp_2j_enumerate_paths(const irrep_multiset_2j_t *a,
                                const irrep_multiset_2j_t *b,
                                const irrep_multiset_2j_t *c,
                                int *out_paths, int max_paths) {
    if (!a || !b || !c) return 0;
    int count = 0;
    for (int ia = 0; ia < a->num_terms; ++ia) {
        int ja = a->labels[ia].two_j, pa = a->labels[ia].parity;
        for (int ib = 0; ib < b->num_terms; ++ib) {
            int jb = b->labels[ib].two_j, pb = b->labels[ib].parity;
            int pc = pa * pb;
            for (int ic = 0; ic < c->num_terms; ++ic) {
                int jc = c->labels[ic].two_j;
                if (c->labels[ic].parity != pc) continue;
                if (!triangle_ok_2j_(ja, jb, jc)) continue;
                if (out_paths && count < max_paths) {
                    out_paths[3*count + 0] = ia;
                    out_paths[3*count + 1] = ib;
                    out_paths[3*count + 2] = ic;
                }
                ++count;
            }
        }
    }
    return count;
}

/* -------------------------------------------------------------------------- *
 * Build                                                                      *
 * -------------------------------------------------------------------------- */

tp_2j_descriptor_t *irrep_tp_2j_build(const irrep_multiset_2j_t *a,
                                      const irrep_multiset_2j_t *b,
                                      const irrep_multiset_2j_t *c,
                                      const int *selected_paths,
                                      int num_selected_paths) {
    if (!a || !b || !c) {
        irrep_set_error_("irrep_tp_2j_build: NULL multiset");
        return NULL;
    }

    /* Auto-enumerate paths if none supplied. */
    int *owned_paths = NULL;
    int  num_paths   = num_selected_paths;
    const int *paths_src = selected_paths;
    if (!selected_paths) {
        int total = irrep_tp_2j_enumerate_paths(a, b, c, NULL, 0);
        if (total == 0) {
            irrep_set_error_("irrep_tp_2j_build: no valid paths");
            return NULL;
        }
        owned_paths = malloc((size_t)total * 3 * sizeof(int));
        if (!owned_paths) {
            irrep_set_error_("irrep_tp_2j_build: out of memory");
            return NULL;
        }
        irrep_tp_2j_enumerate_paths(a, b, c, owned_paths, total);
        num_paths = total;
        paths_src = owned_paths;
    }

    tp_2j_descriptor_t *desc = calloc(1, sizeof(*desc));
    if (!desc) {
        free(owned_paths);
        irrep_set_error_("irrep_tp_2j_build: out of memory");
        return NULL;
    }
    desc->paths = calloc((size_t)num_paths, sizeof(*desc->paths));
    if (!desc->paths) {
        free(desc); free(owned_paths);
        irrep_set_error_("irrep_tp_2j_build: out of memory");
        return NULL;
    }
    desc->num_paths = num_paths;
    desc->a_dim = irrep_multiset_2j_dim(a);
    desc->b_dim = irrep_multiset_2j_dim(b);
    desc->c_dim = irrep_multiset_2j_dim(c);

    int wbase = 0;
    for (int p = 0; p < num_paths; ++p) {
        int ia = paths_src[3*p + 0];
        int ib = paths_src[3*p + 1];
        int ic = paths_src[3*p + 2];
        if (ia < 0 || ia >= a->num_terms ||
            ib < 0 || ib >= b->num_terms ||
            ic < 0 || ic >= c->num_terms) {
            irrep_set_error_("irrep_tp_2j_build: path %d out of range", p);
            goto fail;
        }
        struct tp_2j_path *pp = &desc->paths[p];
        pp->i_a = ia; pp->i_b = ib; pp->i_c = ic;
        pp->two_j_a = a->labels[ia].two_j;
        pp->two_j_b = b->labels[ib].two_j;
        pp->two_j_c = c->labels[ic].two_j;
        pp->d_a = pp->two_j_a + 1;
        pp->d_b = pp->two_j_b + 1;
        pp->d_c = pp->two_j_c + 1;

        if (!triangle_ok_2j_(pp->two_j_a, pp->two_j_b, pp->two_j_c)) {
            irrep_set_error_("irrep_tp_2j_build: path %d violates triangle", p);
            goto fail;
        }
        if (a->labels[ia].parity * b->labels[ib].parity != c->labels[ic].parity) {
            irrep_set_error_("irrep_tp_2j_build: path %d violates parity", p);
            goto fail;
        }

        pp->mult_u = a->multiplicities[ia];
        pp->mult_v = b->multiplicities[ib];
        pp->mult_w = c->multiplicities[ic];

        pp->offset_a = multiset_2j_offset_(a, ia);
        pp->offset_b = multiset_2j_offset_(b, ib);
        pp->offset_c = multiset_2j_offset_(c, ic);

        pp->weight_offset = wbase;
        wbase += pp->mult_u * pp->mult_v * pp->mult_w;

        /* CG block: cg[m_a * d_b * d_c + m_b * d_c + m_c] */
        pp->cg = malloc((size_t)pp->d_a * pp->d_b * pp->d_c * sizeof(double));
        if (!pp->cg) {
            irrep_set_error_("irrep_tp_2j_build: out of memory (CG block)");
            goto fail;
        }
        for (int ma_i = 0; ma_i < pp->d_a; ++ma_i) {
            int two_ma = 2 * ma_i - pp->two_j_a;
            for (int mb_i = 0; mb_i < pp->d_b; ++mb_i) {
                int two_mb = 2 * mb_i - pp->two_j_b;
                for (int mc_i = 0; mc_i < pp->d_c; ++mc_i) {
                    int two_mc = 2 * mc_i - pp->two_j_c;
                    double cg = (two_ma + two_mb == two_mc)
                              ? irrep_cg_2j(pp->two_j_a, two_ma,
                                            pp->two_j_b, two_mb,
                                            pp->two_j_c, two_mc)
                              : 0.0;
                    pp->cg[(size_t)ma_i * pp->d_b * pp->d_c
                         + (size_t)mb_i * pp->d_c + mc_i] = cg;
                }
            }
        }
    }
    desc->num_weights = wbase;

    free(owned_paths);
    return desc;

fail:
    if (desc) {
        if (desc->paths) {
            for (int p = 0; p < desc->num_paths; ++p) free(desc->paths[p].cg);
            free(desc->paths);
        }
        free(desc);
    }
    free(owned_paths);
    return NULL;
}

void irrep_tp_2j_free(tp_2j_descriptor_t *desc) {
    if (!desc) return;
    if (desc->paths) {
        for (int p = 0; p < desc->num_paths; ++p) free(desc->paths[p].cg);
        free(desc->paths);
    }
    free(desc);
}

int irrep_tp_2j_output_dim(const tp_2j_descriptor_t *desc) {
    return desc ? desc->c_dim : 0;
}
int irrep_tp_2j_num_paths(const tp_2j_descriptor_t *desc) {
    return desc ? desc->num_paths : 0;
}

/* -------------------------------------------------------------------------- *
 * Apply: c = a ⊗ b with unit weights (each path contributes w=v=u=1 copy).   *
 * -------------------------------------------------------------------------- */

void irrep_tp_2j_apply(const tp_2j_descriptor_t *desc,
                       const double _Complex *a_in,
                       const double _Complex *b_in,
                       double _Complex *c_out) {
    if (!desc || !a_in || !b_in || !c_out) return;
    memset(c_out, 0, (size_t)desc->c_dim * sizeof(double _Complex));

    for (int p = 0; p < desc->num_paths; ++p) {
        const struct tp_2j_path *pp = &desc->paths[p];
        int d_a = pp->d_a, d_b = pp->d_b, d_c = pp->d_c;

        /* If mults disagree, treat as uvw with all-ones weight over the
         * overlapping copy-channels u==v==w. For true UVW, callers should
         * use apply_weighted with an explicit weight tensor. */
        int mult = pp->mult_u;
        if (pp->mult_v < mult) mult = pp->mult_v;
        if (pp->mult_w < mult) mult = pp->mult_w;

        for (int u = 0; u < mult; ++u) {
            const double _Complex *a_block =
                a_in + pp->offset_a + (size_t)u * d_a;
            const double _Complex *b_block =
                b_in + pp->offset_b + (size_t)u * d_b;
            double _Complex *c_block =
                c_out + pp->offset_c + (size_t)u * d_c;

            for (int mc = 0; mc < d_c; ++mc) {
                double _Complex acc = 0.0;
                for (int ma = 0; ma < d_a; ++ma) {
                    double _Complex a_v = a_block[ma];
                    for (int mb = 0; mb < d_b; ++mb) {
                        double cg = pp->cg[(size_t)ma * d_b * d_c
                                         + (size_t)mb * d_c + mc];
                        if (cg == 0.0) continue;
                        acc += cg * a_v * b_block[mb];
                    }
                }
                c_block[mc] += acc;
            }
        }
    }
}

/* -------------------------------------------------------------------------- *
 * Weighted apply:                                                            *
 *   c[i_c, w_ch, m_c] += Σ_{u, v} weight[i_p, w_ch, v, u] · Σ_{m_a, m_b}     *
 *                         CG · a[i_a, u, m_a] · b[i_b, v, m_b]                *
 * weights are complex; shape path-major, (w, v, u) innermost.                *
 * -------------------------------------------------------------------------- */

void irrep_tp_2j_apply_weighted(const tp_2j_descriptor_t *desc,
                                const double _Complex *weights,
                                const double _Complex *a_in,
                                const double _Complex *b_in,
                                double _Complex *c_out) {
    if (!desc || !a_in || !b_in || !c_out || !weights) return;
    memset(c_out, 0, (size_t)desc->c_dim * sizeof(double _Complex));

    for (int p = 0; p < desc->num_paths; ++p) {
        const struct tp_2j_path *pp = &desc->paths[p];
        int d_a = pp->d_a, d_b = pp->d_b, d_c = pp->d_c;
        int U = pp->mult_u, V = pp->mult_v, W = pp->mult_w;
        const double _Complex *wblk = weights + pp->weight_offset;

        for (int wch = 0; wch < W; ++wch) {
            double _Complex *c_block =
                c_out + pp->offset_c + (size_t)wch * d_c;
            for (int v = 0; v < V; ++v) {
                const double _Complex *b_block =
                    b_in + pp->offset_b + (size_t)v * d_b;
                for (int u = 0; u < U; ++u) {
                    const double _Complex *a_block =
                        a_in + pp->offset_a + (size_t)u * d_a;
                    double _Complex w = wblk[(size_t)wch * V * U
                                           + (size_t)v  * U + u];
                    if (w == 0.0) continue;

                    for (int mc = 0; mc < d_c; ++mc) {
                        double _Complex acc = 0.0;
                        for (int ma = 0; ma < d_a; ++ma) {
                            double _Complex a_v = a_block[ma];
                            for (int mb = 0; mb < d_b; ++mb) {
                                double cg = pp->cg[(size_t)ma * d_b * d_c
                                                 + (size_t)mb * d_c + mc];
                                if (cg == 0.0) continue;
                                acc += cg * a_v * b_block[mb];
                            }
                        }
                        c_block[mc] += w * acc;
                    }
                }
            }
        }
    }
}

/* -------------------------------------------------------------------------- *
 * Backward: transpose the recipe against grad_c_out.                         *
 *   grad_w[i_p, w, v, u] += Σ_{m_c} conj(a ⊗ b reduced) · grad_c_out[m_c]    *
 *   grad_a[i_a, u, m_a]   += Σ_{v, w, m_b, m_c} weight · CG ·                *
 *                            conj(b[i_b, v, m_b]) · grad_c_out[m_c]          *
 *   (symmetric for grad_b)                                                   *
 *                                                                            *
 * Complex conjugation: since TP is bilinear in (a, b) and complex-linear in  *
 * weights, the backward uses conj on the *other* factors when pulling back   *
 * against grad_c (Wirtinger convention; caller's autograd expects            *
 * ∂/∂conj(x) = 0). For all-real amplitudes the conj is a no-op and this     *
 * reduces to the integer-l real TP backward.                                 *
 * -------------------------------------------------------------------------- */

void irrep_tp_2j_apply_backward(const tp_2j_descriptor_t *desc,
                                const double _Complex *weights,
                                const double _Complex *a_in,
                                const double _Complex *b_in,
                                const double _Complex *grad_c_out,
                                double _Complex *grad_a,
                                double _Complex *grad_b,
                                double _Complex *grad_w) {
    if (!desc || !weights || !a_in || !b_in || !grad_c_out) return;

    for (int p = 0; p < desc->num_paths; ++p) {
        const struct tp_2j_path *pp = &desc->paths[p];
        int d_a = pp->d_a, d_b = pp->d_b, d_c = pp->d_c;
        int U = pp->mult_u, V = pp->mult_v, W = pp->mult_w;
        const double _Complex *wblk       = weights    + pp->weight_offset;
        double _Complex       *gwblk      = grad_w ? grad_w + pp->weight_offset : NULL;

        for (int wch = 0; wch < W; ++wch) {
            const double _Complex *gc_block =
                grad_c_out + pp->offset_c + (size_t)wch * d_c;

            for (int v = 0; v < V; ++v) {
                const double _Complex *b_block =
                    b_in + pp->offset_b + (size_t)v * d_b;
                double _Complex *gb_block = grad_b
                    ? grad_b + pp->offset_b + (size_t)v * d_b : NULL;

                for (int u = 0; u < U; ++u) {
                    const double _Complex *a_block =
                        a_in + pp->offset_a + (size_t)u * d_a;
                    double _Complex *ga_block = grad_a
                        ? grad_a + pp->offset_a + (size_t)u * d_a : NULL;

                    double _Complex w_val = wblk[(size_t)wch * V * U
                                               + (size_t)v  * U + u];

                    /* Accumulate scalar for grad_w[p, wch, v, u] plus
                     * contribute to grad_a / grad_b. */
                    double _Complex gw_acc = 0.0;
                    for (int mc = 0; mc < d_c; ++mc) {
                        double _Complex gc = gc_block[mc];
                        for (int ma = 0; ma < d_a; ++ma) {
                            for (int mb = 0; mb < d_b; ++mb) {
                                double cg = pp->cg[(size_t)ma * d_b * d_c
                                                 + (size_t)mb * d_c + mc];
                                if (cg == 0.0) continue;
                                gw_acc += cg * a_block[ma] * b_block[mb] * gc;

                                if (ga_block) {
                                    ga_block[ma] += w_val * cg * b_block[mb] * gc;
                                }
                                if (gb_block) {
                                    gb_block[mb] += w_val * cg * a_block[ma] * gc;
                                }
                            }
                        }
                    }
                    if (gwblk) {
                        gwblk[(size_t)wch * V * U + (size_t)v * U + u] += gw_acc;
                    }
                }
            }
        }
    }
}
