/* SPDX-License-Identifier: MIT */
/* Point-group character tables and projector.
 *
 * Elements are stored as `(proper rotation matrix, det)` pairs — an improper
 * element has `det = -1` and its full action on an `(l, parity)` irrep is
 * `parity · D^l(R_proper)`. Character tables follow Bradley-Cracknell Vol. 1;
 * they are laid out element-by-element (not class-by-class) so the
 * projector loop is index-indexed and doesn't need a class mapping table.
 *
 * The real-basis Wigner-D matrix `D^l(R)` is built lazily per call via
 *   D_real = U · D_complex · U†,  U = complex-to-real-SH basis change
 * using the existing @ref irrep_wigner_D_matrix + @ref irrep_sph_harm_complex_to_real
 * primitives. At v1.2 scope (C₄ᵥ, D₆ only; ℓ ≤ 8) the cost is negligible
 * next to the matrix-vector application. */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/point_group.h>
#include <irrep/so3.h>
#include <irrep/spherical_harmonics.h>
#include <irrep/types.h>
#include <irrep/wigner_d.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern void irrep_set_error_(const char *fmt, ...);

/* -------------------------------------------------------------------------- *
 * Table layout                                                               *
 * -------------------------------------------------------------------------- */

typedef struct {
    irrep_rot_matrix_t R;   /* proper rotation part */
    int                det; /* +1 (proper) or -1 (improper / reflection) */
} pg_element_t;

/* Real-basis D^l matrices are cached at table-build time for l ≤ this bound.
 * Covers the NequIP / MACE / Allegro `l_sh_max` sweet spot without bloating
 * the table footprint. Larger l falls back to on-the-fly recomputation in
 * irrep_pg_project — no regression, just no cache hit. */
#define IRREP_PG_CACHED_L_MAX 4

/* Per-l block size (2l+1)² summed across l = 0..L, giving the cache stride
 * per group element. At L = 4: 1 + 9 + 25 + 49 + 81 = 165 doubles. */
static int cache_block_size_(int l_max) {
    int s = 0;
    for (int l = 0; l <= l_max; ++l)
        s += (2 * l + 1) * (2 * l + 1);
    return s;
}

static int cache_block_offset_(int l) {
    int s = 0;
    for (int ll = 0; ll < l; ++ll)
        s += (2 * ll + 1) * (2 * ll + 1);
    return s;
}

struct irrep_pg_table {
    irrep_point_group_t group;
    int                 num_irreps;
    int                 order;
    int                *irrep_dims;   /* [num_irreps] */
    const char        **irrep_labels; /* [num_irreps] — static pointers */
    double             *characters;   /* [num_irreps * order] row-major */
    pg_element_t       *elements;     /* [order] */
    /* Flat cache: D^l(g) for l = 0..IRREP_PG_CACHED_L_MAX, every g ∈ G.
     * Indexing: `D_cache + g * stride + cache_block_offset_(l) + i`. */
    double *D_cache;
    int     D_cache_stride; /* = cache_block_size_(IRREP_PG_CACHED_L_MAX) */
};

/* Forward decl: real-basis D^l(R) matrix, defined below. */
static void real_d_matrix_l_(int l, irrep_rot_matrix_t R, double *out);

/* -------------------------------------------------------------------------- *
 * Rotation helpers                                                           *
 * -------------------------------------------------------------------------- */

static irrep_rot_matrix_t rot_z_(double theta) {
    irrep_axis_angle_t aa = {.axis = {0.0, 0.0, 1.0}, .angle = theta};
    return irrep_rot_from_axis_angle(aa);
}

static irrep_rot_matrix_t rot_axis_(double ax, double ay, double az, double theta) {
    double             n = sqrt(ax * ax + ay * ay + az * az);
    irrep_axis_angle_t aa = {.axis = {ax / n, ay / n, az / n}, .angle = theta};
    return irrep_rot_from_axis_angle(aa);
}

/* -------------------------------------------------------------------------- *
 * C₄ᵥ  —  order 8, 5 irreps                                                  *
 * Elements: {E, C₄, C₂, C₄³, σ_v(xz), σ_v(yz), σ_d(x=y), σ_d(x=-y)}.         *
 * Character table (Bradley-Cracknell).                                       *
 * -------------------------------------------------------------------------- */

/*                          E   C4   C2   C4³   σvx  σvy  σd1  σd2  */
static const double C4V_CHAR[5 * 8] = {
    /* A1 */ 1, 1,  1,  1,  1,  1,  1,  1,
    /* A2 */ 1, 1,  1,  1,  -1, -1, -1, -1,
    /* B1 */ 1, -1, 1,  -1, 1,  1,  -1, -1,
    /* B2 */ 1, -1, 1,  -1, -1, -1, 1,  1,
    /* E  */ 2, 0,  -2, 0,  0,  0,  0,  0,
};
static const int         C4V_DIMS[5] = {1, 1, 1, 1, 2};
static const char *const C4V_LABELS[5] = {"A1", "A2", "B1", "B2", "E"};

static void              fill_c4v_elements_(pg_element_t *e) {
    e[0].R = irrep_rot_identity();
    e[0].det = +1; /* E   */
    e[1].R = rot_z_(M_PI / 2.0);
    e[1].det = +1; /* C4  */
    e[2].R = rot_z_(M_PI);
    e[2].det = +1; /* C2  */
    e[3].R = rot_z_(-M_PI / 2.0);
    e[3].det = +1; /* C4³ */
    /* σ_v through xz-plane ≡ improper · Ry(π).
     * σ_v through yz-plane ≡ improper · Rx(π).
     * σ_d diagonals likewise rotate (x, y) → (y, x) etc. composed with σ_h.
     * Using proper rotation = (reflection) · (inversion); det = -1 marks
     * improper. The "proper part" for σ_v(xz) is Ry(π): (x,y,z) → (-x,y,-z),
     * then times inversion gives (x,-y,z) — the xz-plane reflection. */
    e[4].R = rot_axis_(0.0, 1.0, 0.0, M_PI);
    e[4].det = -1; /* σ_v(xz) */
    e[5].R = rot_axis_(1.0, 0.0, 0.0, M_PI);
    e[5].det = -1; /* σ_v(yz) */
    /* σ_d(x=y): maps (x,y,z) → (y,x,-z). That equals inversion · Rn(π) where
     * n is the axis along (1, -1, 0)/√2. */
    e[6].R = rot_axis_(1.0, -1.0, 0.0, M_PI);
    e[6].det = -1; /* σ_d(x=y) */
    e[7].R = rot_axis_(1.0, 1.0, 0.0, M_PI);
    e[7].det = -1; /* σ_d(x=-y) */
}

/* -------------------------------------------------------------------------- *
 * D₆  —  order 12, 6 irreps                                                  *
 * Elements: {E, C6, C3, C2, C3², C6⁵} ∪ {3 C₂′ (edge midpoints) + 3 C₂″ (vertices)} *
 * Bradley-Cracknell table.                                                   *
 * -------------------------------------------------------------------------- */

/*                          E    C6   C3   C2   C3²  C6⁵  C2′1 C2′2 C2′3 C2″1 C2″2 C2″3 */
static const double D6_CHAR[6 * 12] = {
    /* A1 */ 1, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    /* A2 */ 1, 1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, -1,
    /* B1 */ 1, -1, 1,  -1, 1,  -1, 1,  1,  1,  -1, -1, -1,
    /* B2 */ 1, -1, 1,  -1, 1,  -1, -1, -1, -1, 1,  1,  1,
    /* E1 */ 2, 1,  -1, -2, -1, 1,  0,  0,  0,  0,  0,  0,
    /* E2 */ 2, -1, -1, 2,  -1, -1, 0,  0,  0,  0,  0,  0,
};
static const int         D6_DIMS[6] = {1, 1, 1, 1, 2, 2};
static const char *const D6_LABELS[6] = {"A1", "A2", "B1", "B2", "E1", "E2"};

static void              fill_d6_elements_(pg_element_t *e) {
    /* 6 proper rotations about z. */
    for (int k = 0; k < 6; ++k) {
        e[k].R = rot_z_((double)k * M_PI / 3.0); /* 0°, 60°, 120°, 180°, 240°, 300° */
        e[k].det = +1;
    }
    /* 3 C₂′ axes in xy-plane through edge midpoints (angle 30°, 90°, 150° from x). */
    const double primed_angles[3] = {M_PI / 6.0, M_PI / 2.0, 5.0 * M_PI / 6.0};
    for (int k = 0; k < 3; ++k) {
        double a = primed_angles[k];
        e[6 + k].R = rot_axis_(cos(a), sin(a), 0.0, M_PI);
        e[6 + k].det = +1; /* C₂′ in D₆ is a PROPER rotation by π about a horizontal axis. */
    }
    /* 3 C₂″ axes in xy-plane through vertices (angle 0°, 60°, 120° from x). */
    const double doubled_angles[3] = {0.0, M_PI / 3.0, 2.0 * M_PI / 3.0};
    for (int k = 0; k < 3; ++k) {
        double a = doubled_angles[k];
        e[9 + k].R = rot_axis_(cos(a), sin(a), 0.0, M_PI);
        e[9 + k].det = +1;
    }
}

/* -------------------------------------------------------------------------- *
 * C₃ᵥ  —  order 6, 3 irreps                                                  *
 * Elements: {E, C₃, C₃², σ_v(1), σ_v(2), σ_v(3)}.                             *
 * Mirror planes at angles φ ∈ {0, 2π/3, 4π/3} containing the z-axis.         *
 * Bradley-Cracknell table.                                                   *
 * -------------------------------------------------------------------------- */

/*                          E    C3   C3²  σv1  σv2  σv3 */
static const double C3V_CHAR[3 * 6] = {
    /* A1 */ 1, 1,  1,  1,  1,  1,
    /* A2 */ 1, 1,  1,  -1, -1, -1,
    /* E  */ 2, -1, -1, 0,  0,  0,
};
static const int         C3V_DIMS[3] = {1, 1, 2};
static const char *const C3V_LABELS[3] = {"A1", "A2", "E"};

static void              fill_c3v_elements_(pg_element_t *e) {
    e[0].R = irrep_rot_identity();
    e[0].det = +1; /* E   */
    e[1].R = rot_z_(2.0 * M_PI / 3.0);
    e[1].det = +1; /* C3  */
    e[2].R = rot_z_(-2.0 * M_PI / 3.0);
    e[2].det = +1; /* C3² */
    /* σ_v reflections through vertical planes at φ = 0, 2π/3, 4π/3.
     * Following the C4v precedent: improper part = inversion ·
     * (proper rotation by π about the normal to the mirror plane).
     * Mirror plane at angle φ has normal (−sin φ, cos φ, 0) (rotated 90°
     * from the plane's in-plane direction). */
    const double sv_angles[3] = {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0};
    for (int k = 0; k < 3; ++k) {
        double a = sv_angles[k];
        e[3 + k].R = rot_axis_(-sin(a), cos(a), 0.0, M_PI);
        e[3 + k].det = -1;
    }
}

/* -------------------------------------------------------------------------- *
 * D₃  —  order 6, 3 irreps                                                   *
 * Elements: {E, C₃, C₃², C₂(1), C₂(2), C₂(3)}.                                *
 * Three C₂ axes in the xy-plane at φ ∈ {0, 2π/3, 4π/3}.                       *
 * Bradley-Cracknell table.                                                   *
 * Note: D₃ is isomorphic to C₃ᵥ as an abstract group; the tables are the    *
 * same, but the geometric realisation (purely proper vs. with improper      *
 * reflections) differs, so `_reduce` behaves differently on parity-odd inputs. *
 * -------------------------------------------------------------------------- */

/*                          E    C3   C3²  C2(1) C2(2) C2(3) */
static const double D3_CHAR[3 * 6] = {
    /* A1 */ 1, 1,  1,  1,  1,  1,
    /* A2 */ 1, 1,  1,  -1, -1, -1,
    /* E  */ 2, -1, -1, 0,  0,  0,
};
static const int         D3_DIMS[3] = {1, 1, 2};
static const char *const D3_LABELS[3] = {"A1", "A2", "E"};

static void              fill_d3_elements_(pg_element_t *e) {
    e[0].R = irrep_rot_identity();
    e[0].det = +1; /* E   */
    e[1].R = rot_z_(2.0 * M_PI / 3.0);
    e[1].det = +1; /* C3  */
    e[2].R = rot_z_(-2.0 * M_PI / 3.0);
    e[2].det = +1; /* C3² */
    /* Three C₂ axes in xy-plane at angles 0, 2π/3, 4π/3. All proper. */
    const double c2_angles[3] = {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0};
    for (int k = 0; k < 3; ++k) {
        double a = c2_angles[k];
        e[3 + k].R = rot_axis_(cos(a), sin(a), 0.0, M_PI);
        e[3 + k].det = +1;
    }
}

/* -------------------------------------------------------------------------- *
 * Table build / destroy                                                      *
 * -------------------------------------------------------------------------- */

irrep_pg_table_t *irrep_pg_table_build(irrep_point_group_t g) {
    irrep_pg_table_t *t = calloc(1, sizeof(*t));
    if (!t)
        return NULL;
    t->group = g;

    const double *chars_src = NULL;
    switch (g) {
    case IRREP_PG_C4V:
        t->num_irreps = 5;
        t->order = 8;
        t->irrep_dims = (int *)C4V_DIMS;
        t->irrep_labels = (const char **)C4V_LABELS;
        chars_src = C4V_CHAR;
        t->elements = calloc(8, sizeof(*t->elements));
        if (!t->elements) {
            free(t);
            return NULL;
        }
        fill_c4v_elements_(t->elements);
        break;
    case IRREP_PG_D6:
        t->num_irreps = 6;
        t->order = 12;
        t->irrep_dims = (int *)D6_DIMS;
        t->irrep_labels = (const char **)D6_LABELS;
        chars_src = D6_CHAR;
        t->elements = calloc(12, sizeof(*t->elements));
        if (!t->elements) {
            free(t);
            return NULL;
        }
        fill_d6_elements_(t->elements);
        break;
    case IRREP_PG_C3V:
        t->num_irreps = 3;
        t->order = 6;
        t->irrep_dims = (int *)C3V_DIMS;
        t->irrep_labels = (const char **)C3V_LABELS;
        chars_src = C3V_CHAR;
        t->elements = calloc(6, sizeof(*t->elements));
        if (!t->elements) {
            free(t);
            return NULL;
        }
        fill_c3v_elements_(t->elements);
        break;
    case IRREP_PG_D3:
        t->num_irreps = 3;
        t->order = 6;
        t->irrep_dims = (int *)D3_DIMS;
        t->irrep_labels = (const char **)D3_LABELS;
        chars_src = D3_CHAR;
        t->elements = calloc(6, sizeof(*t->elements));
        if (!t->elements) {
            free(t);
            return NULL;
        }
        fill_d3_elements_(t->elements);
        break;
    default:
        irrep_set_error_("irrep_pg_table_build: unknown point group %d", (int)g);
        free(t);
        return NULL;
    }

    /* Character table: caller may mutate in principle, so duplicate. */
    size_t ch_n = (size_t)t->num_irreps * (size_t)t->order;
    t->characters = malloc(ch_n * sizeof(double));
    if (!t->characters) {
        free(t->elements);
        free(t);
        return NULL;
    }
    memcpy(t->characters, chars_src, ch_n * sizeof(double));

    /* Pre-compute real-basis D^l matrices per element for l ≤ cached max.
     * Cost: `order × block_size × 8 B`. At order=12 (D6) and cached_L=4,
     * that's 12 · 165 · 8 = ~16 KB. On failure, fall back gracefully — the
     * projector's recompute path still works. */
    t->D_cache_stride = cache_block_size_(IRREP_PG_CACHED_L_MAX);
    size_t cache_n = (size_t)t->order * (size_t)t->D_cache_stride;
    t->D_cache = malloc(cache_n * sizeof(double));
    if (t->D_cache) {
        for (int e = 0; e < t->order; ++e) {
            for (int l = 0; l <= IRREP_PG_CACHED_L_MAX; ++l) {
                double *slot =
                    t->D_cache + (size_t)e * (size_t)t->D_cache_stride + cache_block_offset_(l);
                real_d_matrix_l_(l, t->elements[e].R, slot);
            }
        }
    }

    return t;
}

void irrep_pg_table_free(irrep_pg_table_t *t) {
    if (!t)
        return;
    free(t->D_cache);
    free(t->elements);
    free(t->characters);
    /* irrep_dims and irrep_labels point into static storage — do not free. */
    free(t);
}

int irrep_pg_num_irreps(const irrep_pg_table_t *t) {
    return t ? t->num_irreps : 0;
}
int irrep_pg_order(const irrep_pg_table_t *t) {
    return t ? t->order : 0;
}
const char *irrep_pg_irrep_label(const irrep_pg_table_t *t, int mu) {
    if (!t || mu < 0 || mu >= t->num_irreps)
        return NULL;
    return t->irrep_labels[mu];
}

/* -------------------------------------------------------------------------- *
 * Real-basis Wigner-D^l(R) — wraps the complex full-D + U · D · U† bridge.   *
 *                                                                            *
 * Output is row-major `(2l+1) × (2l+1)` real. For l = 0 this is trivially 1. *
 * -------------------------------------------------------------------------- */

static void real_d_matrix_l_(int l, irrep_rot_matrix_t R, double *out) {
    if (l == 0) {
        out[0] = 1.0;
        return;
    }
    int               d = 2 * l + 1;

    irrep_euler_zyz_t e = irrep_euler_zyz_from_rot(R);
    /* Stack space — safe at our l ≤ IRREP_L_MAX. */
    double _Complex D[(2 * IRREP_L_MAX + 1) * (2 * IRREP_L_MAX + 1)];
    double _Complex U[(2 * IRREP_L_MAX + 1) * (2 * IRREP_L_MAX + 1)];
    double _Complex tmp[(2 * IRREP_L_MAX + 1) * (2 * IRREP_L_MAX + 1)];

    irrep_wigner_D_matrix(l, D, e.alpha, e.beta, e.gamma);
    irrep_sph_harm_complex_to_real(l, U);

    /* tmp = U · D */
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            double _Complex s = 0.0;
            for (int k = 0; k < d; ++k)
                s += U[i * d + k] * D[k * d + j];
            tmp[i * d + j] = s;
        }
    }
    /* out = tmp · U† — guaranteed real by construction. */
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            double _Complex s = 0.0;
            for (int k = 0; k < d; ++k)
                s += tmp[i * d + k] * conj(U[j * d + k]);
            out[i * d + j] = creal(s);
        }
    }
}

/* -------------------------------------------------------------------------- *
 * Projection                                                                 *
 * -------------------------------------------------------------------------- */

void irrep_pg_project(const irrep_pg_table_t *t, int mu, const irrep_multiset_t *spec,
                      const double *in, double *out) {
    if (!t || !spec || !in || !out)
        return;
    if (mu < 0 || mu >= t->num_irreps)
        return;
    const int total = spec->total_dim;
    for (int i = 0; i < total; ++i)
        out[i] = 0.0;

    const double  prefactor = (double)t->irrep_dims[mu] / (double)t->order;
    const double *chi = t->characters + (size_t)mu * (size_t)t->order;

    /* Accumulator sized to the largest block, allocated once outside the loop. */
    double block_tmp[2 * IRREP_L_MAX + 1];

    for (int e_idx = 0; e_idx < t->order; ++e_idx) {
        const pg_element_t *g = &t->elements[e_idx];
        double              w = prefactor * chi[e_idx];
        if (w == 0.0)
            continue;

        int offset = 0;
        for (int ti = 0; ti < spec->num_terms; ++ti) {
            int l = spec->labels[ti].l;
            int par = spec->labels[ti].parity;
            int mult = spec->multiplicities[ti];
            int d = 2 * l + 1;

            /* Use the cached D^l when possible; fall back to per-call build
             * for l > cached bound. Cache hit skips Euler extraction +
             * complex-D build + U·D·U† — the dominant `_project` cost. */
            const double *D_block;
            double        D_fallback[(2 * IRREP_L_MAX + 1) * (2 * IRREP_L_MAX + 1)];
            if (l <= IRREP_PG_CACHED_L_MAX && t->D_cache) {
                D_block =
                    t->D_cache + (size_t)e_idx * (size_t)t->D_cache_stride + cache_block_offset_(l);
            } else {
                real_d_matrix_l_(l, g->R, D_fallback);
                D_block = D_fallback;
            }

            double flip = (g->det < 0) ? (double)par : 1.0;
            double coeff = w * flip;
            if (coeff == 0.0) {
                offset += mult * d;
                continue;
            }

            for (int u = 0; u < mult; ++u) {
                const double *src = in + offset + u * d;
                double       *dst = out + offset + u * d;
                /* block_tmp = D_block · src */
                for (int i = 0; i < d; ++i) {
                    double s = 0.0;
                    for (int j = 0; j < d; ++j)
                        s += D_block[i * d + j] * src[j];
                    block_tmp[i] = s;
                }
                for (int i = 0; i < d; ++i)
                    dst[i] += coeff * block_tmp[i];
            }
            offset += mult * d;
        }
    }
}

/* -------------------------------------------------------------------------- *
 * Reduction via character formula                                            *
 *                                                                            *
 *   m_μ = (1/|G|) Σ_g χ*_μ(g) · χ_V(g)                                       *
 *                                                                            *
 * where χ_V(g) = Σ_term mult_term · χ_block(g),                              *
 *       χ_block(g) = trace of real-D_l(R_proper(g)) · (g.det<0 ? parity : 1) *
 *                  = χ_l(θ_g) · parity_flip                                  *
 *       χ_l(θ) = sin((l + 1/2)·θ) / sin(θ/2)                                 *
 * -------------------------------------------------------------------------- */

/* χ_l(θ) — the character of an SO(3) irrep of order l at a rotation by θ.
 * At θ = 0 the limit is 2l + 1. */
static double chi_l_(int l, double theta) {
    if (fabs(theta) < 1e-14)
        return (double)(2 * l + 1);
    double num = sin((double)(2 * l + 1) * theta * 0.5);
    double den = sin(theta * 0.5);
    if (fabs(den) < 1e-14) {
        /* θ = 2π — round trip; also 2l+1. */
        return (double)(2 * l + 1);
    }
    return num / den;
}

/* Extract rotation angle from the proper rotation matrix (abs value). */
static double rot_angle_(irrep_rot_matrix_t R) {
    irrep_axis_angle_t aa = irrep_axis_angle_from_rot(R);
    return fabs(aa.angle);
}

void irrep_pg_reduce(const irrep_pg_table_t *t, const irrep_multiset_t *spec, int *out_mult) {
    if (!t || !spec || !out_mult)
        return;

    /* Precompute χ_V(g) over all group elements. */
    double *chi_V = calloc((size_t)t->order, sizeof(double));
    if (!chi_V) {
        for (int m = 0; m < t->num_irreps; ++m)
            out_mult[m] = 0;
        return;
    }
    for (int e_idx = 0; e_idx < t->order; ++e_idx) {
        const pg_element_t *g = &t->elements[e_idx];
        double              theta = rot_angle_(g->R);
        double              s = 0.0;
        for (int ti = 0; ti < spec->num_terms; ++ti) {
            int    l = spec->labels[ti].l;
            int    par = spec->labels[ti].parity;
            int    mult = spec->multiplicities[ti];
            double chi_block = chi_l_(l, theta);
            if (g->det < 0)
                chi_block *= (double)par;
            s += (double)mult * chi_block;
        }
        chi_V[e_idx] = s;
    }

    for (int mu = 0; mu < t->num_irreps; ++mu) {
        const double *chi = t->characters + (size_t)mu * (size_t)t->order;
        double        acc = 0.0;
        for (int e_idx = 0; e_idx < t->order; ++e_idx)
            acc += chi[e_idx] * chi_V[e_idx];
        double raw = acc / (double)t->order;
        /* Multiplicity must be a non-negative integer; round to the nearest. */
        out_mult[mu] = (int)lround(raw);
    }
    free(chi_V);
}
