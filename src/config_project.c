/* SPDX-License-Identifier: MIT */
/* Configuration-space projection onto a space-group irrep.
 *
 * The heavy lifting — computing ψ(g·σ) for each g — is the caller's job, as
 * it is wavefunction-specific. This module is the group-theoretic reducer:
 * given the amplitudes ψ_g for each group element and a character row for
 * the target irrep, return P_μ ψ(σ).
 */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/config_project.h>
#include <irrep/lattice.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern void irrep_set_error_(const char *fmt, ...);

struct irrep_sg_irrep {
    const irrep_space_group_t *G; /* borrowed */
    int                        order;
    int                        irrep_dim;
    double _Complex           *characters;
};

irrep_sg_irrep_t *irrep_sg_irrep_new(const irrep_space_group_t *G,
                                     const double _Complex *characters, int irrep_dim) {
    if (!G || !characters || irrep_dim <= 0) {
        irrep_set_error_("irrep_sg_irrep_new: invalid arguments");
        return NULL;
    }
    int               order = irrep_space_group_order(G);

    irrep_sg_irrep_t *mu = calloc(1, sizeof(*mu));
    if (!mu) {
        irrep_set_error_("irrep_sg_irrep_new: out of memory");
        return NULL;
    }
    mu->characters = malloc((size_t)order * sizeof(double _Complex));
    if (!mu->characters) {
        free(mu);
        irrep_set_error_("irrep_sg_irrep_new: out of memory");
        return NULL;
    }
    memcpy(mu->characters, characters, (size_t)order * sizeof(double _Complex));
    mu->G = G;
    mu->order = order;
    mu->irrep_dim = irrep_dim;
    return mu;
}

void irrep_sg_irrep_free(irrep_sg_irrep_t *mu) {
    if (!mu)
        return;
    free(mu->characters);
    free(mu);
}

irrep_sg_irrep_t *irrep_sg_trivial(const irrep_space_group_t *G) {
    if (!G)
        return NULL;
    int              order = irrep_space_group_order(G);
    double _Complex *chi = malloc((size_t)order * sizeof(double _Complex));
    if (!chi) {
        irrep_set_error_("irrep_sg_trivial: out of memory");
        return NULL;
    }
    for (int g = 0; g < order; ++g)
        chi[g] = 1.0 + 0.0 * I;
    irrep_sg_irrep_t *mu = irrep_sg_irrep_new(G, chi, 1);
    free(chi);
    return mu;
}

irrep_sg_irrep_t *irrep_sg_sign_rep(const irrep_space_group_t *G) {
    if (!G)
        return NULL;
    int order = irrep_space_group_order(G);
    int point_order = irrep_space_group_point_order(G);
    /* p4mm / p6mm: the builder enumerates rotations first, then mirrors —
     * so index p ∈ [0, point_order/2) is a rotation and p ∈ [point_order/2,
     * point_order) is a reflection. For p1 (point_order = 1) every element
     * is a translation, sign +1 throughout. Full element index is
     * (cell_index · point_order + p); characters are constant across cells
     * at the Γ-point.                                                      */
    int              half = point_order / 2;
    double _Complex *chi = malloc((size_t)order * sizeof(double _Complex));
    if (!chi) {
        irrep_set_error_("irrep_sg_sign_rep: out of memory");
        return NULL;
    }
    for (int g = 0; g < order; ++g) {
        int p = g % point_order;
        int sign = (point_order > 1 && p >= half) ? -1 : +1;
        chi[g] = (double)sign + 0.0 * I;
    }
    irrep_sg_irrep_t *mu = irrep_sg_irrep_new(G, chi, 1);
    free(chi);
    return mu;
}

double _Complex irrep_sg_project_amplitude(const irrep_sg_irrep_t *mu,
                                           const double _Complex  *psi_of_g) {
    if (!mu || !psi_of_g)
        return NAN + NAN * I;
    double _Complex acc = 0.0 + 0.0 * I;
    for (int g = 0; g < mu->order; ++g) {
        acc += conj(mu->characters[g]) * psi_of_g[g];
    }
    return acc * ((double)mu->irrep_dim / (double)mu->order);
}

double _Complex irrep_sg_project_A1(const irrep_space_group_t *G, const double _Complex *psi_of_g) {
    if (!G || !psi_of_g)
        return NAN + NAN * I;
    int order = irrep_space_group_order(G);
    double _Complex acc = 0.0 + 0.0 * I;
    for (int g = 0; g < order; ++g)
        acc += psi_of_g[g];
    return acc / (double)order;
}

void irrep_sg_enumerate_orbit(const irrep_space_group_t *G, const double *sigma,
                              double *out_orbit) {
    if (!G || !sigma || !out_orbit)
        return;
    int  n = irrep_space_group_num_sites(G);
    int  order = irrep_space_group_order(G);
    int *inv = malloc((size_t)n * sizeof(int));
    if (!inv)
        return;
    for (int g = 0; g < order; ++g) {
        irrep_space_group_permutation_inverse(G, g, inv);
        double *dst = out_orbit + (size_t)g * n;
        for (int s = 0; s < n; ++s)
            dst[s] = sigma[inv[s]];
    }
    free(inv);
}

/* -------------------------------------------------------------------------- *
 * Symmetry-adapted basis builder                                             *
 * -------------------------------------------------------------------------- */

static long long ipow_ll_(int base, int exp) {
    long long r = 1;
    for (int k = 0; k < exp; ++k)
        r *= (long long)base;
    return r;
}

/* Apply a site permutation to the digit-encoded index `s`: bit/digit j of
 * the output = bit/digit perm[j] of the input. */
static long long permute_digits_(long long s, int Nsites, int d, const int *perm) {
    if (d == 2) {
        long long out = 0;
        for (int j = 0; j < Nsites; ++j) {
            long long bit = (s >> perm[j]) & 1LL;
            out |= bit << j;
        }
        return out;
    }
    /* Generic base-d case. Extract digits once, then reassemble. */
    int       digits[64];
    long long tmp = s;
    for (int j = 0; j < Nsites; ++j) {
        digits[j] = (int)(tmp % d);
        tmp /= d;
    }
    long long out = 0, w = 1;
    for (int j = 0; j < Nsites; ++j) {
        out += (long long)digits[perm[j]] * w;
        w *= d;
    }
    return out;
}

int irrep_sg_adapted_basis(const irrep_space_group_t *G, const irrep_sg_irrep_t *mu, int num_sites,
                           int local_dim, double _Complex *basis_out, int n_max) {
    if (!G || !mu || !basis_out || n_max <= 0 || num_sites < 1 || local_dim < 2)
        return -1;
    if (num_sites != irrep_space_group_num_sites(G))
        return -1;

    int order = irrep_space_group_order(G);
    if (mu->order != order)
        return -1;

    long long D = ipow_ll_(local_dim, num_sites);
    if (D <= 0)
        return -1;

    /* Cache every group element's inverse permutation once, upfront. Turns
     * the inner loop's permutation lookup into a pointer dereference rather
     * than a memcpy via irrep_space_group_permutation_inverse. */
    int *perm_inv_table = malloc((size_t)order * (size_t)num_sites * sizeof(int));
    if (!perm_inv_table)
        return -1;
    for (int g = 0; g < order; ++g) {
        irrep_space_group_permutation_inverse(G, g, perm_inv_table + (size_t)g * num_sites);
    }

    double _Complex *v = malloc((size_t)D * sizeof(double _Complex));
    if (!v) {
        free(perm_inv_table);
        return -1;
    }

    const double tol = 1e-9;
    int          n_basis = 0;
    double _Complex scale = (double)mu->irrep_dim / (double)order;

    for (long long s = 0; s < D; ++s) {
        /* Orbit-representative filter: skip |s⟩ if some element of its orbit
         * has a smaller index. Within each G-orbit we process exactly one
         * seed — the minimum. This cuts the apply-phase cost by a factor of
         * |G| on typical orbits. Correctness: every orbit rep gets processed,
         * so every basis vector of the μ-isotypic subspace is still reached
         * (Gram-Schmidt rejects the redundant ones; orbits that project to
         * zero contribute nothing, as with the all-seeds loop). */
        int is_min = 1;
        for (int g = 1; g < order && is_min; ++g) {
            const int *perm = perm_inv_table + (size_t)g * num_sites;
            long long  s_g = permute_digits_(s, num_sites, local_dim, perm);
            if (s_g < s)
                is_min = 0;
        }
        if (!is_min)
            continue;

        /* v = P_μ |s⟩ = (d_μ / |G|) Σ_g conj(χ(g)) · g|s⟩ */
        memset(v, 0, (size_t)D * sizeof(double _Complex));
        for (int g = 0; g < order; ++g) {
            const int *perm = perm_inv_table + (size_t)g * num_sites;
            long long  s_g = permute_digits_(s, num_sites, local_dim, perm);
            v[s_g] += scale * conj(mu->characters[g]);
        }

        /* Orthogonalise against accepted basis. */
        for (int k = 0; k < n_basis; ++k) {
            const double _Complex *b_k = basis_out + (size_t)k * D;
            double _Complex overlap = 0.0;
            for (long long t = 0; t < D; ++t)
                overlap += conj(b_k[t]) * v[t];
            for (long long t = 0; t < D; ++t)
                v[t] -= overlap * b_k[t];
        }

        double norm2 = 0.0;
        for (long long t = 0; t < D; ++t)
            norm2 += creal(v[t]) * creal(v[t]) + cimag(v[t]) * cimag(v[t]);
        double norm = sqrt(norm2);

        if (norm > tol && n_basis < n_max) {
            double _Complex *dst = basis_out + (size_t)n_basis * D;
            for (long long t = 0; t < D; ++t)
                dst[t] = v[t] / norm;
            ++n_basis;
        }
        if (n_basis >= n_max)
            break;
    }

    free(v);
    free(perm_inv_table);
    return n_basis;
}

/* -------------------------------------------------------------------------- *
 * Bloch-momentum projection (translation subgroup only)                      *
 *                                                                            *
 * Convention: k = (kx/Lx) b1 + (ky/Ly) b2 with b_i · a_j = 2π δ_ij, so the   *
 * phase on translation t = (tx, ty) is                                       *
 *     e^{-i k·t} = exp(-2π i (kx tx / Lx + ky ty / Ly)).                     *
 * -------------------------------------------------------------------------- */

double _Complex irrep_sg_bloch_amplitude(const irrep_space_group_t *G, int kx, int ky,
                                         const double _Complex *psi_of_g) {
    if (!G || !psi_of_g)
        return NAN + NAN * I;
    const irrep_lattice_t *L = irrep_space_group_lattice(G);
    if (!L)
        return NAN + NAN * I;
    int Lx = irrep_lattice_Lx(L);
    int Ly = irrep_lattice_Ly(L);
    /* Bloch indices live mod (Lx, Ly). Canonicalise so callers can pass any
     * integer (including negatives for signed BZ conventions) without hitting
     * the orbit-rep filter's [0, Lx) assumption in bloch_basis. */
    kx = ((kx % Lx) + Lx) % Lx;
    ky = ((ky % Ly) + Ly) % Ly;
    int point_order = irrep_space_group_point_order(G);
    double _Complex acc = 0.0 + 0.0 * I;
    double inv_Lx = 1.0 / (double)Lx;
    double inv_Ly = 1.0 / (double)Ly;
    for (int ty = 0; ty < Ly; ++ty) {
        for (int tx = 0; tx < Lx; ++tx) {
            int    tidx = ty * Lx + tx;
            int    g = tidx * point_order + 0; /* identity point element */
            double phase = -2.0 * M_PI * ((double)(kx * tx) * inv_Lx + (double)(ky * ty) * inv_Ly);
            double _Complex w = cos(phase) + I * sin(phase);
            acc += w * psi_of_g[g];
        }
    }
    return acc / (double)(Lx * Ly);
}

int irrep_sg_bloch_basis(const irrep_space_group_t *G, int kx, int ky, int num_sites, int local_dim,
                         double _Complex *basis_out, int n_max) {
    if (!G || !basis_out || n_max <= 0 || num_sites < 1 || local_dim < 2)
        return -1;
    if (num_sites != irrep_space_group_num_sites(G))
        return -1;
    const irrep_lattice_t *L = irrep_space_group_lattice(G);
    if (!L)
        return -1;
    int Lx = irrep_lattice_Lx(L);
    int Ly = irrep_lattice_Ly(L);
    kx = ((kx % Lx) + Lx) % Lx;
    ky = ((ky % Ly) + Ly) % Ly;
    int       point_order = irrep_space_group_point_order(G);
    int       num_trans = Lx * Ly;

    long long D = ipow_ll_(local_dim, num_sites);
    if (D <= 0)
        return -1;

    /* Cache translation-only inverse permutations. Group element
     * g = tidx · point_order + 0 has permutation at row (tidx · point_order). */
    int *trans_inv = malloc((size_t)num_trans * (size_t)num_sites * sizeof(int));
    if (!trans_inv)
        return -1;
    int *scratch = malloc((size_t)num_sites * sizeof(int));
    if (!scratch) {
        free(trans_inv);
        return -1;
    }
    for (int tidx = 0; tidx < num_trans; ++tidx) {
        int g = tidx * point_order + 0;
        irrep_space_group_permutation_inverse(G, g, scratch);
        memcpy(trans_inv + (size_t)tidx * num_sites, scratch, (size_t)num_sites * sizeof(int));
    }
    free(scratch);

    /* Precompute Bloch phases for each translation index. */
    double _Complex *phases = malloc((size_t)num_trans * sizeof(double _Complex));
    if (!phases) {
        free(trans_inv);
        return -1;
    }
    double inv_Lx = 1.0 / (double)Lx;
    double inv_Ly = 1.0 / (double)Ly;
    for (int ty = 0; ty < Ly; ++ty) {
        for (int tx = 0; tx < Lx; ++tx) {
            int    tidx = ty * Lx + tx;
            double ph = -2.0 * M_PI * ((double)(kx * tx) * inv_Lx + (double)(ky * ty) * inv_Ly);
            phases[tidx] = cos(ph) + I * sin(ph);
        }
    }

    double _Complex *v = malloc((size_t)D * sizeof(double _Complex));
    if (!v) {
        free(phases);
        free(trans_inv);
        return -1;
    }

    const double tol = 1e-9;
    int          n_basis = 0;
    double _Complex scale = 1.0 / (double)num_trans;

    for (long long s = 0; s < D; ++s) {
        /* Orbit-representative filter over translations. For non-Γ k the
         * orbit phase structure still makes each orbit produce one vector;
         * keep the minimum translate as seed to avoid duplicated work. */
        int is_min = 1;
        for (int tidx = 1; tidx < num_trans && is_min; ++tidx) {
            const int *perm = trans_inv + (size_t)tidx * num_sites;
            long long  s_g = permute_digits_(s, num_sites, local_dim, perm);
            if (s_g < s)
                is_min = 0;
        }
        if (!is_min)
            continue;

        memset(v, 0, (size_t)D * sizeof(double _Complex));
        for (int tidx = 0; tidx < num_trans; ++tidx) {
            const int *perm = trans_inv + (size_t)tidx * num_sites;
            long long  s_g = permute_digits_(s, num_sites, local_dim, perm);
            v[s_g] += scale * phases[tidx];
        }

        for (int k = 0; k < n_basis; ++k) {
            const double _Complex *b_k = basis_out + (size_t)k * D;
            double _Complex overlap = 0.0;
            for (long long t = 0; t < D; ++t)
                overlap += conj(b_k[t]) * v[t];
            for (long long t = 0; t < D; ++t)
                v[t] -= overlap * b_k[t];
        }

        double norm2 = 0.0;
        for (long long t = 0; t < D; ++t)
            norm2 += creal(v[t]) * creal(v[t]) + cimag(v[t]) * cimag(v[t]);
        double norm = sqrt(norm2);

        if (norm > tol && n_basis < n_max) {
            double _Complex *dst = basis_out + (size_t)n_basis * D;
            for (long long t = 0; t < D; ++t)
                dst[t] = v[t] / norm;
            ++n_basis;
        }
        if (n_basis >= n_max)
            break;
    }

    free(v);
    free(phases);
    free(trans_inv);
    return n_basis;
}

/* -------------------------------------------------------------------------- *
 * Little group at Bloch momentum k                                           *
 *                                                                            *
 * For each pure-point-group element p ∈ [0, point_order) we need to decide   *
 * whether R_p fixes k (mod reciprocal lattice). We extract R_p as a 2x2      *
 * integer matrix by applying p to sublattice-0 sites at cells (1, 0) and    *
 * (0, 1) and reading the cell decomposition of the image. This relies on    *
 * the wallpaper-group builder placing sublattice 0 at the rotation centre   *
 * (true for the p4mm / p6mm constructors currently shipped), so that a      *
 * point operation around the origin never induces an extra cell shift on    *
 * top of the Bravais rotation.                                              *
 *                                                                            *
 * Once R_p is in hand, the fix-point test is                                *
 *     R_p · (kx, ky) ≡ (kx, ky)  (mod (Lx, Ly)).                             *
 * -------------------------------------------------------------------------- */

struct irrep_sg_little_group {
    const irrep_space_group_t *G; /* borrowed */
    int                        kx, ky;
    int                        point_order;   /* cardinality of little point group */
    int                       *point_ops;     /* length = point_order */
    int                        num_translations;
    /* Cached 2×2 real-space rotation matrix per little-point element (in
     * signed-representative form, same as little_group_extract_rotation_).
     * 4 × point_order ints. */
    int                       *matrices;
};

struct irrep_sg_little_group_irrep {
    const irrep_sg_little_group_t *lg; /* borrowed */
    int                            dim;
    int                            point_order;
    double _Complex               *characters; /* length = point_order */
    /* Optional full D_μ(g) matrix per element. NULL for 1D irreps (characters
     * suffice) or for hand-supplied _irrep_new calls where only characters are
     * known. For 2D named irreps, filled by _irrep_named from the element
     * geometric action (rotation angle / mirror axis) in cartesian. Stored
     * row-major as [point_order][dim * dim] contiguous. */
    double _Complex               *matrices;
};

/* Extract the 2x2 integer matrix R_p (in lattice-cell basis) by measuring
 * image-cell displacements. For a pure point element `p` centred at the
 * wallpaper-group origin O (possibly not a lattice site, e.g. hexagon
 * centre on kagome), images of same-sublattice sites transform as
 *
 *     cell(p · site_at(cell=v, sub=0)) = R · v + c_p
 *
 * where `c_p` depends on how p maps sublattice 0 to some sublattice k_p
 * (and is independent of v). Taking the difference cancels c_p:
 *
 *     R · (1, 0) = cell(p · (1, 0)) - cell(p · (0, 0))
 *     R · (0, 1) = cell(p · (0, 1)) - cell(p · (0, 0))
 *
 * Differences are reduced modulo (Lx, Ly) to fold wraparound back onto the
 * torus, which is all we need for the subsequent fix-point test. */
static int little_group_extract_rotation_(const irrep_space_group_t *G, int p, int R[2][2]) {
    const irrep_lattice_t *L = irrep_space_group_lattice(G);
    if (!L)
        return -1;
    int Lx = irrep_lattice_Lx(L);
    int Ly = irrep_lattice_Ly(L);
    if (Lx < 2 || Ly < 2)
        return -1;

    int g = p;
    int s00 = irrep_lattice_site_index(L, 0, 0, 0);
    int s10 = irrep_lattice_site_index(L, 1, 0, 0);
    int s01 = irrep_lattice_site_index(L, 0, 1, 0);
    if (s00 < 0 || s10 < 0 || s01 < 0)
        return -1;

    int s00_img = irrep_space_group_apply(G, g, s00);
    int s10_img = irrep_space_group_apply(G, g, s10);
    int s01_img = irrep_space_group_apply(G, g, s01);
    if (s00_img < 0 || s10_img < 0 || s01_img < 0)
        return -1;

    int ix00, iy00, ix10, iy10, ix01, iy01;
    if (irrep_lattice_cell_of(L, s00_img, &ix00, &iy00) != IRREP_OK ||
        irrep_lattice_cell_of(L, s10_img, &ix10, &iy10) != IRREP_OK ||
        irrep_lattice_cell_of(L, s01_img, &ix01, &iy01) != IRREP_OK)
        return -1;

    R[0][0] = ((ix10 - ix00) % Lx + Lx) % Lx;
    R[1][0] = ((iy10 - iy00) % Ly + Ly) % Ly;
    R[0][1] = ((ix01 - ix00) % Lx + Lx) % Lx;
    R[1][1] = ((iy01 - iy00) % Ly + Ly) % Ly;

    /* Lift to signed representatives in (-L/2, L/2] so the det computation
     * and M^{-T} that follow are exact integer arithmetic rather than
     * mod-(Lx·Ly) algebra. */
    if (R[0][0] > Lx / 2)
        R[0][0] -= Lx;
    if (R[0][1] > Lx / 2)
        R[0][1] -= Lx;
    if (R[1][0] > Ly / 2)
        R[1][0] -= Ly;
    if (R[1][1] > Ly / 2)
        R[1][1] -= Ly;
    return 0;
}

/* Signed modulo into [0, m). `m` must be positive. */
static int smod_(int v, int m) { return ((v % m) + m) % m; }

/* Does point op `p` fix (kx, ky) modulo (Lx, Ly) on this lattice?
 *
 * `little_group_extract_rotation_` returns the real-space action M on lattice
 * vectors: (R_p · a_j) in a-basis = column j of M. Under R_p, the reciprocal
 * basis (b1, b2) transforms by the inverse transpose, so on a k-vector stored
 * as (kx, ky) ∈ (b1, b2)-coords the action is M^{-T}:
 *
 *     k_new = M^{-T} · k.
 *
 * For |det M| = 1 (every wallpaper-group lattice symmetry), M^{-T} has
 * integer entries, so the check stays in exact integer arithmetic. For square
 * lattices M happens to be orthogonal, so M^{-T} = M and the distinction
 * doesn't matter — but on triangular / hexagonal lattices the two differ
 * and only M^{-T} gives the correct K-point stabiliser. */
static int point_op_fixes_k_(const irrep_space_group_t *G, int p, int kx, int ky, int Lx, int Ly) {
    int M[2][2];
    if (little_group_extract_rotation_(G, p, M) != 0)
        return 0;
    int det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
    if (det != 1 && det != -1)
        return 0; /* not a lattice symmetry — shouldn't happen for p4mm/p6mm */
    /* M^{-T} = (1/det) · [[M[1][1], -M[1][0]], [-M[0][1], M[0][0]]] */
    int MinvT[2][2] = {{det * M[1][1], -det * M[1][0]}, {-det * M[0][1], det * M[0][0]}};
    int nkx = smod_(MinvT[0][0] * kx + MinvT[0][1] * ky, Lx);
    int nky = smod_(MinvT[1][0] * kx + MinvT[1][1] * ky, Ly);
    return (nkx == kx) && (nky == ky);
}

irrep_sg_little_group_t *irrep_sg_little_group_build(const irrep_space_group_t *G, int kx, int ky) {
    if (!G) {
        irrep_set_error_("irrep_sg_little_group_build: NULL space group");
        return NULL;
    }
    const irrep_lattice_t *L = irrep_space_group_lattice(G);
    if (!L) {
        irrep_set_error_("irrep_sg_little_group_build: space group has no lattice handle");
        return NULL;
    }
    int Lx = irrep_lattice_Lx(L);
    int Ly = irrep_lattice_Ly(L);
    kx = smod_(kx, Lx);
    ky = smod_(ky, Ly);

    int point_order = irrep_space_group_point_order(G);

    /* p1: only identity in the point group; trivially fixes every k. */
    /* p4mm / p6mm: run the fix-point test on each of the 8 / 12 elements. */
    int *ops = malloc((size_t)point_order * sizeof(int));
    if (!ops) {
        irrep_set_error_("irrep_sg_little_group_build: out of memory");
        return NULL;
    }
    int n_ops = 0;
    for (int p = 0; p < point_order; ++p) {
        if (point_op_fixes_k_(G, p, kx, ky, Lx, Ly)) {
            ops[n_ops++] = p;
        }
    }

    irrep_sg_little_group_t *lg = calloc(1, sizeof(*lg));
    if (!lg) {
        free(ops);
        irrep_set_error_("irrep_sg_little_group_build: out of memory");
        return NULL;
    }
    /* Shrink to fit. */
    int *shrunk = realloc(ops, (size_t)(n_ops > 0 ? n_ops : 1) * sizeof(int));
    lg->G = G;
    lg->kx = kx;
    lg->ky = ky;
    lg->point_order = n_ops;
    lg->point_ops = shrunk ? shrunk : ops;
    lg->num_translations = Lx * Ly;

    /* Cache the rotation matrix of each selected point op. A fresh call to
     * _element_matrix would re-do the extraction on every query; a few
     * tens of ints here keep that O(1). */
    lg->matrices = malloc((size_t)n_ops * 4 * sizeof(int));
    if (lg->matrices) {
        for (int i = 0; i < n_ops; ++i) {
            int M[2][2];
            if (little_group_extract_rotation_(G, lg->point_ops[i], M) == 0) {
                int *base = lg->matrices + (size_t)i * 4;
                base[0] = M[0][0];
                base[1] = M[0][1];
                base[2] = M[1][0];
                base[3] = M[1][1];
            } else {
                /* Should not happen given build-time filtering. */
                int *base = lg->matrices + (size_t)i * 4;
                base[0] = base[3] = 1;
                base[1] = base[2] = 0;
            }
        }
    }
    return lg;
}

void irrep_sg_little_group_free(irrep_sg_little_group_t *lg) {
    if (!lg)
        return;
    free(lg->matrices);
    free(lg->point_ops);
    free(lg);
}

void irrep_sg_little_group_element_matrix(const irrep_sg_little_group_t *lg, int i,
                                          int out_M[2][2]) {
    if (!lg || i < 0 || i >= lg->point_order || !out_M || !lg->matrices)
        return;
    const int *base = lg->matrices + (size_t)i * 4;
    out_M[0][0] = base[0];
    out_M[0][1] = base[1];
    out_M[1][0] = base[2];
    out_M[1][1] = base[3];
}

/* ----- Little-group-irrep handle + composite projector ----- */

irrep_sg_little_group_irrep_t *irrep_sg_little_group_irrep_new(const irrep_sg_little_group_t *lg,
                                                               const double _Complex *characters,
                                                               int                    dim) {
    if (!lg || !characters || dim <= 0) {
        irrep_set_error_("irrep_sg_little_group_irrep_new: invalid arguments");
        return NULL;
    }
    irrep_sg_little_group_irrep_t *mu = calloc(1, sizeof(*mu));
    if (!mu) {
        irrep_set_error_("irrep_sg_little_group_irrep_new: out of memory");
        return NULL;
    }
    mu->characters = malloc((size_t)lg->point_order * sizeof(double _Complex));
    if (!mu->characters) {
        free(mu);
        irrep_set_error_("irrep_sg_little_group_irrep_new: out of memory");
        return NULL;
    }
    memcpy(mu->characters, characters, (size_t)lg->point_order * sizeof(double _Complex));
    mu->lg = lg;
    mu->dim = dim;
    mu->point_order = lg->point_order;
    return mu;
}

int irrep_sg_little_group_irrep_matrix(const irrep_sg_little_group_irrep_t *mu_k,
                                       int i, double _Complex *D_out) {
    if (!mu_k || !D_out || i < 0 || i >= mu_k->point_order)
        return -1;
    int d = mu_k->dim;
    if (d == 1) {
        D_out[0] = mu_k->characters[i];
        return 0;
    }
    if (!mu_k->matrices)
        return -1; /* 2D irrep built via _irrep_new (no matrix data) */
    memcpy(D_out, mu_k->matrices + (size_t)i * d * d,
           (size_t)(d * d) * sizeof(double _Complex));
    return 0;
}

void irrep_sg_little_group_irrep_free(irrep_sg_little_group_irrep_t *mu_k) {
    if (!mu_k)
        return;
    free(mu_k->matrices);
    free(mu_k->characters);
    free(mu_k);
}

int irrep_sg_little_group_irrep_dim(const irrep_sg_little_group_irrep_t *mu_k) {
    return mu_k ? mu_k->dim : 0;
}

/* ---- p6mm parent-op → conjugacy-class table ------------------------
 * Enumerated once on 6×6 kagome (where the lattice-rotation-matrix
 * extraction is unambiguous) and cached here. Valid for any cluster
 * size: the space-group builder indexes point ops the same way
 * regardless of `(Lx, Ly)`.
 *
 * Classes (in the order the C_6v character table expects):
 *   0 = E, 1 = C_6, 2 = C_3, 3 = C_2, 4 = σ_v, 5 = σ_d.
 *
 * Mirror split {6, 8, 10} σ_v vs {7, 9, 11} σ_d verified by explicit
 * conjugation `C_3 · M_6 · C_3⁻¹ = M_8` — geometric axes at
 * 0°/60°/120° (σ_v) vs 30°/90°/150° (σ_d) in cartesian. */
static const int p6mm_point_class_[12] = {
    0,                   /* parent 0: E         */
    1, 2, 3, 2, 1,       /* parent 1..5: C_6, C_3, C_2, C_3², C_6⁵ */
    4, 5, 4, 5, 4, 5,    /* parent 6..11: σ_v, σ_d, σ_v, σ_d, σ_v, σ_d */
};

/* C_6v character table, rows indexed by named irrep subset, columns by
 * conjugacy class (E, C_6, C_3, C_2, σ_v, σ_d). Unused named irreps are
 * flagged with INT_MIN in the dim-lookup so caller gets a clean error. */
static const double c6v_chi_[6][6] = {
    /*            E    C_6   C_3   C_2   σ_v   σ_d  */
    /* A_1 */  { +1,  +1,  +1,  +1,  +1,  +1 },
    /* A_2 */  { +1,  +1,  +1,  +1,  -1,  -1 },
    /* B_1 */  { +1,  -1,  +1,  -1,  +1,  -1 },
    /* B_2 */  { +1,  -1,  +1,  -1,  -1,  +1 },
    /* E_1 */  { +2,  +1,  -1,  -2,   0,   0 },
    /* E_2 */  { +2,  -1,  -1,  +2,   0,   0 },
};
static const int c6v_dim_[6] = {1, 1, 1, 1, 2, 2};

/* C_3v at the K-point on p6mm kagome. Little group {E, C_3, C_3², 3σ}
 * with parent indices {0, 2, 4, σ_1, σ_2, σ_3}. Classes (E, C_3, σ_v).
 *
 * Which three of the six C_6v mirrors appear at a given K depends on the
 * K-point, but all three are in the SAME C_6v class (σ_v or σ_d
 * exclusively). The C_3v classification below uses that invariant to
 * pick the class-2 slot. */
static const double c3v_chi_[3][3] = {
    /*          E    2C_3   3σ   */
    /* A_1 */  { 1,    1,    1 },
    /* A_2 */  { 1,    1,   -1 },
    /* E   */  { 2,   -1,    0 },
};
static const int c3v_dim_[3] = {1, 1, 2};

/* C_2v (order 4). Appears at M-points on p6mm and at X-points on p4mm.
 * The four elements always enumerate in ascending parent-index order as
 * (E, C_2, σ, σ'). σ vs σ' distinguishes the two mirror-axis classes —
 * we take σ = class 2 and σ' = class 3, consistent throughout the
 * library's test expectations. */
static const double c2v_chi_[4][4] = {
    /*          E    C_2   σ    σ'  */
    /* A_1 */  { 1,    1,   1,   1 },
    /* A_2 */  { 1,    1,  -1,  -1 },
    /* B_1 */  { 1,   -1,   1,  -1 },
    /* B_2 */  { 1,   -1,  -1,   1 },
};
static const int c2v_dim_[4] = {1, 1, 1, 1};

/* p4mm parent-op → conjugacy-class table (enumerated on a 4×4 square):
 *   0 = E, 1 = C_4, 2 = C_2, 3 = σ_v, 4 = σ_d. */
static const int p4mm_point_class_[8] = {
    0,       /* 0: E                  */
    1, 2, 1, /* 1: C_4, 2: C_2, 3: C_4³ */
    3, 3,    /* 4, 5: σ_v (axes along x / y) */
    4, 4,    /* 6, 7: σ_d (diagonal axes)    */
};

/* C_4v character table, rows indexed by A_1 / A_2 / B_1 / B_2 / E. */
static const double c4v_chi_[5][5] = {
    /*            E   2C_4   C_2   2σ_v   2σ_d */
    /* A_1 */  { +1,   +1,   +1,   +1,   +1 },
    /* A_2 */  { +1,   +1,   +1,   -1,   -1 },
    /* B_1 */  { +1,   -1,   +1,   +1,   -1 },
    /* B_2 */  { +1,   -1,   +1,   -1,   +1 },
    /* E   */  { +2,    0,   -2,    0,    0 },
};
static const int c4v_dim_[5] = {1, 1, 1, 1, 2};

/* =======================================================================
 * 2D-irrep matrix representations D_μ(g) per parent op.
 *
 * For 2D named irreps (C_3v E, C_6v E_1/E_2, C_4v E), store per-parent-op
 * 2×2 real orthogonal matrices. Indexed by the parent-group op index,
 * not the little-group slot index — the irrep builder then selects via
 * lg->point_ops[i] to get the matrix for the i-th little-group element.
 *
 * Conventions:
 *   Rotation by angle θ:  R(θ) = [[cos θ, -sin θ], [sin θ, cos θ]]
 *   Mirror at axis α:     M(α) = [[cos 2α, sin 2α], [sin 2α, -cos 2α]]
 *
 * p6mm parent ops (12): E, C_6, C_3, C_2, C_3², C_6⁵, σ_v(0°), σ_d(30°),
 *                       σ_v(60°), σ_d(90°), σ_v(120°), σ_d(150°).
 * p4mm parent ops (8):  E, C_4, C_2, C_4³, σ_v(x), σ_v(y), σ_d(45°), σ_d(135°).
 * ======================================================================= */

/* sqrt(3)/2 as a compile-time-ish constant (runtime since sqrt() isn't
 * constexpr in C). Used throughout the 2D-irrep tables. */
static double S3H(void) { return 0.5 * 1.7320508075688772; /* √3 / 2 */ }

/* Fill D[2][2] for the C_6v E_1 irrep on parent op `p ∈ [0, 12)`.
 * Rotation angle doubles: E_1 uses D(C_6) = R(π/3), D(C_3) = R(2π/3) etc.
 * Mirror angles single: D(σ at α) = [[c 2α, s 2α], [s 2α, -c 2α]]. */
static void fill_D_c6v_E1_(int p, double _Complex D[4]) {
    double c = 0.0, s = 0.0;
    double s3h = S3H();
    switch (p) {
        case 0:  /* E */             c = 1.0;  s = 0.0;  break;
        case 1:  /* C_6, θ=π/3 */    c = 0.5;  s = s3h;  break;
        case 2:  /* C_3, θ=2π/3 */   c = -0.5; s = s3h;  break;
        case 3:  /* C_2, θ=π */      c = -1.0; s = 0.0;  break;
        case 4:  /* C_3², θ=4π/3 */  c = -0.5; s = -s3h; break;
        case 5:  /* C_6⁵, θ=5π/3 */  c = 0.5;  s = -s3h; break;
        default: break; /* mirrors handled below */
    }
    if (p < 6) {
        D[0] = c; D[1] = -s; D[2] = s; D[3] = c;
        return;
    }
    /* Mirrors: σ_v axes at 0°, 60°, 120°; σ_d axes at 30°, 90°, 150°.
     * Slot 6,7,8,9,10,11 → α = 0, π/6, π/3, π/2, 2π/3, 5π/6. */
    double alpha;
    switch (p) {
        case 6:  alpha = 0.0;                          break;
        case 7:  alpha = 0.5235987755982988;            break; /* π/6 */
        case 8:  alpha = 1.0471975511965976;            break; /* π/3 */
        case 9:  alpha = 1.5707963267948966;            break; /* π/2 */
        case 10: alpha = 2.0943951023931953;            break; /* 2π/3 */
        case 11: alpha = 2.6179938779914944;            break; /* 5π/6 */
        default: alpha = 0.0; break;
    }
    double c2 = cos(2.0 * alpha), s2 = sin(2.0 * alpha);
    D[0] = c2;  D[1] = s2;
    D[2] = s2;  D[3] = -c2;
}

/* C_6v E_2: rotations doubled (D(C_6) = R(2π/3)), mirror axes doubled too
 * (D(σ at α) = [[c 4α, s 4α], [s 4α, -c 4α]]).
 *
 * Characters check: tr(R(2π/3)) = -1 (C_6 in E_2); tr(R(4π/3)) = -1 (C_3);
 * tr(R(2π)) = 2 (C_2); tr(M(·)) = 0. All match the C_6v E_2 character
 * row {2, -1, -1, 2, 0, 0} exactly. */
static void fill_D_c6v_E2_(int p, double _Complex D[4]) {
    double c = 0.0, s = 0.0;
    double s3h = S3H();
    switch (p) {
        case 0:  c = 1.0;  s = 0.0;  break; /* E: θ=0 */
        case 1:  c = -0.5; s = s3h;  break; /* C_6: θ=2π/3 */
        case 2:  c = -0.5; s = -s3h; break; /* C_3: θ=4π/3 */
        case 3:  c = 1.0;  s = 0.0;  break; /* C_2: θ=2π ≡ 0 */
        case 4:  c = -0.5; s = s3h;  break; /* C_3²: θ=8π/3 ≡ 2π/3 */
        case 5:  c = -0.5; s = -s3h; break; /* C_6⁵: θ=10π/3 ≡ 4π/3 */
        default: break;
    }
    if (p < 6) {
        D[0] = c; D[1] = -s; D[2] = s; D[3] = c;
        return;
    }
    /* Mirror axes at doubled angle: α → 2α, so reflection matrix uses 4α. */
    double alpha;
    switch (p) {
        case 6:  alpha = 0.0;                          break;
        case 7:  alpha = 0.5235987755982988;            break;
        case 8:  alpha = 1.0471975511965976;            break;
        case 9:  alpha = 1.5707963267948966;            break;
        case 10: alpha = 2.0943951023931953;            break;
        case 11: alpha = 2.6179938779914944;            break;
        default: alpha = 0.0; break;
    }
    double c4 = cos(4.0 * alpha), s4 = sin(4.0 * alpha);
    D[0] = c4;  D[1] = s4;
    D[2] = s4;  D[3] = -c4;
}

/* C_3v E: subgroup of C_6v restricted to {E, C_3, C_3², 3σ}. The three
 * mirrors at a K-point on kagome are those at 0°, 60°, 120° (σ_v-type).
 * Takes a LITTLE-group slot index i and fills D accordingly. Parent ops
 * at kagome K-point are {0, 2, 4, 6, 8, 10}. */
static void fill_D_c3v_E_(int p, double _Complex D[4]) {
    double s3h = S3H();
    switch (p) {
        case 0:  /* E */
            D[0] = 1.0; D[1] = 0.0; D[2] = 0.0; D[3] = 1.0;
            return;
        case 2:  /* C_3 */
            D[0] = -0.5; D[1] = -s3h; D[2] = s3h; D[3] = -0.5;
            return;
        case 4:  /* C_3² */
            D[0] = -0.5; D[1] = s3h; D[2] = -s3h; D[3] = -0.5;
            return;
        /* σ_v mirrors at 0°, 60°, 120° */
        case 6:  /* α=0 */
            D[0] = 1.0; D[1] = 0.0; D[2] = 0.0; D[3] = -1.0; return;
        case 8:  /* α=π/3 → M(2π/3) = [[-1/2, √3/2], [√3/2, 1/2]] */
            D[0] = -0.5; D[1] = s3h; D[2] = s3h; D[3] = 0.5; return;
        case 10: /* α=2π/3 → M(4π/3) = [[-1/2, -√3/2], [-√3/2, 1/2]] */
            D[0] = -0.5; D[1] = -s3h; D[2] = -s3h; D[3] = 0.5; return;
        /* σ_d subset (if the K-point picks those instead) — handled via
         * same formula applied at odd parent indices. */
        case 7:  D[0] = 0.5;  D[1] = s3h; D[2] = s3h;  D[3] = -0.5; return;
        case 9:  D[0] = -1.0; D[1] = 0.0; D[2] = 0.0;  D[3] = 1.0;  return;
        case 11: D[0] = 0.5;  D[1] = -s3h; D[2] = -s3h; D[3] = -0.5; return;
        default:
            D[0] = 1.0; D[1] = 0.0; D[2] = 0.0; D[3] = 1.0;
            return;
    }
}

/* C_4v E on p4mm: parent ops {0=E, 1=C_4, 2=C_2, 3=C_4³, 4=σ_v x-axis,
 * 5=σ_v y-axis, 6=σ_d at 45°, 7=σ_d at 135°}. */
static void fill_D_c4v_E_(int p, double _Complex D[4]) {
    switch (p) {
        case 0: D[0] = 1.0;  D[1] = 0.0;  D[2] = 0.0;  D[3] = 1.0;  return; /* E */
        case 1: D[0] = 0.0;  D[1] = -1.0; D[2] = 1.0;  D[3] = 0.0;  return; /* C_4: R(π/2) */
        case 2: D[0] = -1.0; D[1] = 0.0;  D[2] = 0.0;  D[3] = -1.0; return; /* C_2 */
        case 3: D[0] = 0.0;  D[1] = 1.0;  D[2] = -1.0; D[3] = 0.0;  return; /* C_4³: R(3π/2) */
        case 4: D[0] = 1.0;  D[1] = 0.0;  D[2] = 0.0;  D[3] = -1.0; return; /* σ_v(x) */
        case 5: D[0] = -1.0; D[1] = 0.0;  D[2] = 0.0;  D[3] = 1.0;  return; /* σ_v(y) */
        case 6: D[0] = 0.0;  D[1] = 1.0;  D[2] = 1.0;  D[3] = 0.0;  return; /* σ_d(π/4) */
        case 7: D[0] = 0.0;  D[1] = -1.0; D[2] = -1.0; D[3] = 0.0;  return; /* σ_d(3π/4) */
        default:
            D[0] = 1.0; D[1] = 0.0; D[2] = 0.0; D[3] = 1.0;
            return;
    }
}

/* Dispatch helper: map (abstract little-group type, parent-class, slot
 * index, name) → the character-table column (`lg_cls`). Returns -1 on
 * failure and writes an error. The little-group type is inferred from the
 * order + parent space-group kind. */
static int classify_lg_slot_(int parent_group_order, int lg_order, int parent_cls, int slot,
                             int parent_op, int *lg_cls) {
    if (parent_group_order == 12) {
        /* p6mm. Sub-little-groups: C_6v (12), C_3v (6), C_2v (4), C_s/C_1 (small). */
        if (lg_order == 12) {
            *lg_cls = parent_cls; /* C_6v uses 6 classes keyed 0..5 */
            return 0;
        }
        if (lg_order == 6) {
            /* C_3v: three rotations (E, 2 C_3) + three mirrors (single class). */
            if (parent_cls == 0)        { *lg_cls = 0; return 0; }
            if (parent_cls == 2)        { *lg_cls = 1; return 0; }
            if (parent_cls == 4 || parent_cls == 5) { *lg_cls = 2; return 0; }
            irrep_set_error_("irrep_sg_little_group_irrep_named: unexpected element at C_3v "
                             "slot %d (parent %d, p6mm class %d)",
                             slot, parent_op, parent_cls);
            return -1;
        }
        if (lg_order == 4) {
            /* C_2v at M on p6mm: {E, C_2, σ, σ'}. Parent classes for the
             * four slots are (E=0, C_2=3, σ_v=4, σ_d=5) in some order; the
             * natural slot-ordering from `_little_group_point_ops` is
             * ascending parent index which gives exactly {E, C_2, σ_v, σ_d}
             * for M_a=(1,0)->[0,3,7,10], M_b=(0,1)->[0,3,6,9],
             * M_c=(1,1)->[0,3,8,11] on the 2×2 kagome probe. We map the
             * two mirror-class parents {σ_v, σ_d} into C_2v classes {σ=2,
             * σ'=3} by the ORDER the slot appears (first mirror → class 2,
             * second mirror → class 3). */
            if (parent_cls == 0)                       { *lg_cls = 0; return 0; }
            if (parent_cls == 3)                       { *lg_cls = 1; return 0; }
            if (parent_cls == 4 || parent_cls == 5) {
                /* Determine σ vs σ' by mirror index within this little group. */
                *lg_cls = (slot == 2) ? 2 : 3;
                return 0;
            }
            irrep_set_error_("irrep_sg_little_group_irrep_named: unexpected element at C_2v "
                             "slot %d (parent %d, p6mm class %d)",
                             slot, parent_op, parent_cls);
            return -1;
        }
    }
    if (parent_group_order == 8) {
        /* p4mm. Sub-little-groups: C_4v (8), C_2v (4). */
        if (lg_order == 8) {
            *lg_cls = parent_cls;
            return 0;
        }
        if (lg_order == 4) {
            /* C_2v at X on p4mm. Parent classes (E=0, C_4=1, C_2=2, σ_v=3,
             * σ_d=4). The four slots always enumerate as {E, C_2, σ, σ'};
             * rotations are E and C_2 (classes 0 and 2); mirrors are σ or σ_d. */
            if (parent_cls == 0) { *lg_cls = 0; return 0; }
            if (parent_cls == 2) { *lg_cls = 1; return 0; }
            if (parent_cls == 3 || parent_cls == 4) {
                *lg_cls = (slot == 2) ? 2 : 3;
                return 0;
            }
            irrep_set_error_("irrep_sg_little_group_irrep_named: unexpected element at C_2v "
                             "slot %d (parent %d, p4mm class %d)",
                             slot, parent_op, parent_cls);
            return -1;
        }
    }
    irrep_set_error_("irrep_sg_little_group_irrep_named: unsupported (parent order=%d, "
                     "little-group order=%d) combination",
                     parent_group_order, lg_order);
    return -1;
}

irrep_sg_little_group_irrep_t *
irrep_sg_little_group_irrep_named(const irrep_sg_little_group_t *lg, irrep_lg_named_irrep_t name) {
    if (!lg) {
        irrep_set_error_("irrep_sg_little_group_irrep_named: NULL little group");
        return NULL;
    }
    int n = lg->point_order;
    int parent_order = irrep_space_group_point_order(lg->G);

    /* Pick the correct character table and validate name vs shape. */
    int             idx_in_table;
    int             dim;
    const double   *row;
    int             classes;

    if (n == 12 && parent_order == 12) {
        /* C_6v */
        switch (name) {
            case IRREP_LG_IRREP_A1: idx_in_table = 0; break;
            case IRREP_LG_IRREP_A2: idx_in_table = 1; break;
            case IRREP_LG_IRREP_B1: idx_in_table = 2; break;
            case IRREP_LG_IRREP_B2: idx_in_table = 3; break;
            case IRREP_LG_IRREP_E1: idx_in_table = 4; break;
            case IRREP_LG_IRREP_E2: idx_in_table = 5; break;
            default:
                irrep_set_error_("irrep_sg_little_group_irrep_named: named irrep not valid on C_6v");
                return NULL;
        }
        dim = c6v_dim_[idx_in_table];
        row = c6v_chi_[idx_in_table];
        classes = 6;
    } else if (n == 6 && parent_order == 12) {
        /* C_3v at K */
        switch (name) {
            case IRREP_LG_IRREP_A1: idx_in_table = 0; break;
            case IRREP_LG_IRREP_A2: idx_in_table = 1; break;
            case IRREP_LG_IRREP_E:  idx_in_table = 2; break;
            default:
                irrep_set_error_("irrep_sg_little_group_irrep_named: named irrep not valid on C_3v");
                return NULL;
        }
        dim = c3v_dim_[idx_in_table];
        row = c3v_chi_[idx_in_table];
        classes = 3;
    } else if (n == 4) {
        /* C_2v — appears at M on p6mm AND at X on p4mm. Classes: E, C_2, σ, σ'. */
        switch (name) {
            case IRREP_LG_IRREP_A1: idx_in_table = 0; break;
            case IRREP_LG_IRREP_A2: idx_in_table = 1; break;
            case IRREP_LG_IRREP_B1: idx_in_table = 2; break;
            case IRREP_LG_IRREP_B2: idx_in_table = 3; break;
            default:
                irrep_set_error_("irrep_sg_little_group_irrep_named: named irrep not valid on C_2v");
                return NULL;
        }
        dim = c2v_dim_[idx_in_table];
        row = c2v_chi_[idx_in_table];
        classes = 4;
    } else if (n == 8 && parent_order == 8) {
        /* C_4v */
        switch (name) {
            case IRREP_LG_IRREP_A1:    idx_in_table = 0; break;
            case IRREP_LG_IRREP_A2:    idx_in_table = 1; break;
            case IRREP_LG_IRREP_B1:    idx_in_table = 2; break;
            case IRREP_LG_IRREP_B2:    idx_in_table = 3; break;
            case IRREP_LG_IRREP_E_C4V: idx_in_table = 4; break;
            default:
                irrep_set_error_("irrep_sg_little_group_irrep_named: named irrep not valid on C_4v");
                return NULL;
        }
        dim = c4v_dim_[idx_in_table];
        row = c4v_chi_[idx_in_table];
        classes = 5;
    } else {
        irrep_set_error_("irrep_sg_little_group_irrep_named: no builtin for (parent order=%d, "
                         "little-group order=%d); use irrep_sg_little_group_irrep_new with a "
                         "hand-supplied character row",
                         parent_order, n);
        return NULL;
    }

    /* Build the per-element character row. */
    double _Complex *chi = malloc((size_t)n * sizeof(double _Complex));
    if (!chi) {
        irrep_set_error_("irrep_sg_little_group_irrep_named: out of memory");
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        int parent = lg->point_ops[i];
        if (parent < 0 || parent >= parent_order) {
            free(chi);
            irrep_set_error_("irrep_sg_little_group_irrep_named: parent op out of range");
            return NULL;
        }
        int parent_cls =
            (parent_order == 12) ? p6mm_point_class_[parent] : p4mm_point_class_[parent];
        int lg_cls;
        if (classify_lg_slot_(parent_order, n, parent_cls, i, parent, &lg_cls) != 0) {
            free(chi);
            return NULL;
        }
        if (lg_cls < 0 || lg_cls >= classes) {
            free(chi);
            irrep_set_error_("irrep_sg_little_group_irrep_named: class index out of range");
            return NULL;
        }
        chi[i] = row[lg_cls] + 0.0 * I;
    }

    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi, dim);
    free(chi);
    if (!mu)
        return NULL;

    /* Fill D-matrix table for 2D named irreps. */
    if (dim == 2) {
        mu->matrices = malloc((size_t)n * 4 * sizeof(double _Complex));
        if (!mu->matrices) {
            irrep_sg_little_group_irrep_free(mu);
            irrep_set_error_("irrep_sg_little_group_irrep_named: OOM matrix alloc");
            return NULL;
        }
        for (int i = 0; i < n; ++i) {
            int p = lg->point_ops[i];
            double _Complex D[4];
            if (parent_order == 12 && n == 12) {
                /* C_6v — E_1 or E_2. */
                if (name == IRREP_LG_IRREP_E1) fill_D_c6v_E1_(p, D);
                else /* E_2 */                 fill_D_c6v_E2_(p, D);
            } else if (parent_order == 12 && n == 6) {
                /* C_3v E */
                fill_D_c3v_E_(p, D);
            } else if (parent_order == 8 && n == 8) {
                /* C_4v E */
                fill_D_c4v_E_(p, D);
            } else {
                /* C_2v 2D-irrep shouldn't reach here (no E on C_2v). */
                D[0] = 1; D[1] = 0; D[2] = 0; D[3] = 1;
            }
            memcpy(mu->matrices + (size_t)i * 4, D, 4 * sizeof(double _Complex));
        }
    }
    return mu;
}

double _Complex irrep_sg_project_at_k(const irrep_sg_little_group_t       *lg,
                                      const irrep_sg_little_group_irrep_t *mu_k,
                                      const double _Complex               *psi_of_g) {
    if (!lg || !mu_k || !psi_of_g || mu_k->lg != lg)
        return NAN + NAN * I;
    const irrep_lattice_t *L = irrep_space_group_lattice(lg->G);
    if (!L)
        return NAN + NAN * I;
    int Lx = irrep_lattice_Lx(L);
    int Ly = irrep_lattice_Ly(L);
    int point_order = irrep_space_group_point_order(lg->G);
    int num_trans = lg->num_translations;

    double inv_Lx = 1.0 / (double)Lx;
    double inv_Ly = 1.0 / (double)Ly;

    double _Complex acc = 0.0 + 0.0 * I;
    for (int ty = 0; ty < Ly; ++ty) {
        for (int tx = 0; tx < Lx; ++tx) {
            int    tidx = ty * Lx + tx;
            double ph =
                -2.0 * M_PI * ((double)(lg->kx * tx) * inv_Lx + (double)(lg->ky * ty) * inv_Ly);
            double _Complex bloch = cos(ph) + I * sin(ph);
            for (int i = 0; i < lg->point_order; ++i) {
                int             p = lg->point_ops[i];
                int             g = tidx * point_order + p;
                double _Complex chi_conj = conj(mu_k->characters[i]);
                acc += bloch * chi_conj * psi_of_g[g];
            }
        }
    }
    int group_order = num_trans * lg->point_order;
    return acc * ((double)mu_k->dim / (double)group_order);
}

void irrep_sg_canonicalise(const irrep_space_group_t *G, uint64_t config_in,
                           uint64_t *rep_out, int *g_idx_out) {
    if (!G) {
        if (rep_out)   *rep_out   = config_in;
        if (g_idx_out) *g_idx_out = 0;
        return;
    }
    int order = irrep_space_group_order(G);
    uint64_t best  = config_in;
    int      best_g = 0;
    for (int g = 0; g < order; ++g) {
        uint64_t cand = irrep_space_group_apply_bits(G, g, config_in);
        if (cand < best) {
            best   = cand;
            best_g = g;
        }
    }
    if (rep_out)   *rep_out   = best;
    if (g_idx_out) *g_idx_out = best_g;
}

int irrep_sg_stabiliser(const irrep_space_group_t *G, uint64_t config, int *out_indices) {
    if (!G || !out_indices)
        return 0;
    int order = irrep_space_group_order(G);
    int count = 0;
    for (int g = 0; g < order; ++g) {
        if (irrep_space_group_apply_bits(G, g, config) == config)
            out_indices[count++] = g;
    }
    return count;
}

int irrep_sg_orbit_size(const irrep_space_group_t *G, uint64_t config) {
    if (!G)
        return 0;
    int order = irrep_space_group_order(G);
    int stabiliser = 0;
    for (int g = 0; g < order; ++g) {
        if (irrep_space_group_apply_bits(G, g, config) == config)
            ++stabiliser;
    }
    if (stabiliser == 0)
        return 0;
    return order / stabiliser;
}

int irrep_sg_projector_weights(const irrep_sg_little_group_t       *lg,
                               const irrep_sg_little_group_irrep_t *mu_k,
                               double _Complex                     *weights_out) {
    if (!lg || !mu_k || !weights_out || mu_k->lg != lg)
        return -1;
    const irrep_lattice_t *L = irrep_space_group_lattice(lg->G);
    if (!L)
        return -1;
    int Lx = irrep_lattice_Lx(L);
    int Ly = irrep_lattice_Ly(L);
    int point_order = irrep_space_group_point_order(lg->G);
    int num_trans = lg->num_translations;
    int full_point_denom = Lx * Ly * point_order;

    for (int g = 0; g < full_point_denom; ++g)
        weights_out[g] = 0.0 + 0.0 * I;

    double inv_Lx = 1.0 / (double)Lx;
    double inv_Ly = 1.0 / (double)Ly;
    int    group_order = num_trans * lg->point_order;
    double scale = (double)mu_k->dim / (double)group_order;

    for (int ty = 0; ty < Ly; ++ty) {
        for (int tx = 0; tx < Lx; ++tx) {
            int    tidx = ty * Lx + tx;
            double ph =
                -2.0 * M_PI * ((double)(lg->kx * tx) * inv_Lx + (double)(lg->ky * ty) * inv_Ly);
            double _Complex bloch = cos(ph) + I * sin(ph);
            for (int i = 0; i < lg->point_order; ++i) {
                int p = lg->point_ops[i];
                int g = tidx * point_order + p;
                weights_out[g] = scale * bloch * conj(mu_k->characters[i]);
            }
        }
    }
    return 0;
}

int irrep_sg_adapted_basis_at_k(const irrep_sg_little_group_t       *lg,
                                const irrep_sg_little_group_irrep_t *mu_k, int num_sites,
                                int local_dim, double _Complex *basis_out, int n_max) {
    if (!lg || !mu_k || !basis_out || n_max <= 0 || num_sites < 1 || local_dim < 2 ||
        mu_k->lg != lg)
        return -1;
    if (num_sites != irrep_space_group_num_sites(lg->G))
        return -1;
    const irrep_lattice_t *L = irrep_space_group_lattice(lg->G);
    if (!L)
        return -1;
    int Lx = irrep_lattice_Lx(L);
    int Ly = irrep_lattice_Ly(L);
    int point_order = irrep_space_group_point_order(lg->G);
    int num_trans = Lx * Ly;
    int lg_point = lg->point_order;
    int group_order = num_trans * lg_point;

    long long D = ipow_ll_(local_dim, num_sites);
    if (D <= 0)
        return -1;

    /* Cache the inverse permutation of every space-group element that lies
     * inside the little group — i.e. (tidx, p) with p ∈ lg->point_ops. Stored
     * as a flat [num_trans × lg_point × num_sites] int tensor. */
    int *perm_inv = malloc((size_t)num_trans * (size_t)lg_point * (size_t)num_sites * sizeof(int));
    if (!perm_inv)
        return -1;
    int *scratch = malloc((size_t)num_sites * sizeof(int));
    if (!scratch) {
        free(perm_inv);
        return -1;
    }
    for (int tidx = 0; tidx < num_trans; ++tidx) {
        for (int i = 0; i < lg_point; ++i) {
            int p = lg->point_ops[i];
            int g = tidx * point_order + p;
            irrep_space_group_permutation_inverse(lg->G, g, scratch);
            size_t off = ((size_t)tidx * lg_point + i) * num_sites;
            memcpy(perm_inv + off, scratch, (size_t)num_sites * sizeof(int));
        }
    }
    free(scratch);

    /* Precompute per-(tidx, i) phase weights: e^{-i k·t} · conj(χ_μ(p_i)). */
    double _Complex *weights =
        malloc((size_t)num_trans * (size_t)lg_point * sizeof(double _Complex));
    if (!weights) {
        free(perm_inv);
        return -1;
    }
    double inv_Lx = 1.0 / (double)Lx;
    double inv_Ly = 1.0 / (double)Ly;
    for (int ty = 0; ty < Ly; ++ty) {
        for (int tx = 0; tx < Lx; ++tx) {
            int    tidx = ty * Lx + tx;
            double ph =
                -2.0 * M_PI * ((double)(lg->kx * tx) * inv_Lx + (double)(lg->ky * ty) * inv_Ly);
            double _Complex bloch = cos(ph) + I * sin(ph);
            for (int i = 0; i < lg_point; ++i) {
                weights[(size_t)tidx * lg_point + i] = bloch * conj(mu_k->characters[i]);
            }
        }
    }

    double _Complex *v = malloc((size_t)D * sizeof(double _Complex));
    if (!v) {
        free(weights);
        free(perm_inv);
        return -1;
    }

    const double tol = 1e-9;
    int          n_basis = 0;
    double _Complex scale = (double)mu_k->dim / (double)group_order;

    for (long long s = 0; s < D; ++s) {
        /* Orbit-representative filter over the FULL little group: skip |s⟩
         * if any (tidx, p ∈ P_k) maps it to a smaller index — that smaller
         * index is (or will be) the seed for this orbit. */
        int is_min = 1;
        for (int tidx = 0; tidx < num_trans && is_min; ++tidx) {
            for (int i = 0; i < lg_point && is_min; ++i) {
                if (tidx == 0 && i == 0)
                    continue; /* identity */
                const int *perm = perm_inv + ((size_t)tidx * lg_point + i) * num_sites;
                long long  s_g = permute_digits_(s, num_sites, local_dim, perm);
                if (s_g < s)
                    is_min = 0;
            }
        }
        if (!is_min)
            continue;

        memset(v, 0, (size_t)D * sizeof(double _Complex));
        for (int tidx = 0; tidx < num_trans; ++tidx) {
            for (int i = 0; i < lg_point; ++i) {
                const int      *perm = perm_inv + ((size_t)tidx * lg_point + i) * num_sites;
                long long       s_g = permute_digits_(s, num_sites, local_dim, perm);
                double _Complex w = weights[(size_t)tidx * lg_point + i];
                v[s_g] += scale * w;
            }
        }

        for (int k = 0; k < n_basis; ++k) {
            const double _Complex *b_k = basis_out + (size_t)k * D;
            double _Complex        overlap = 0.0;
            for (long long t = 0; t < D; ++t)
                overlap += conj(b_k[t]) * v[t];
            for (long long t = 0; t < D; ++t)
                v[t] -= overlap * b_k[t];
        }

        double norm2 = 0.0;
        for (long long t = 0; t < D; ++t)
            norm2 += creal(v[t]) * creal(v[t]) + cimag(v[t]) * cimag(v[t]);
        double norm = sqrt(norm2);

        if (norm > tol && n_basis < n_max) {
            double _Complex *dst = basis_out + (size_t)n_basis * D;
            for (long long t = 0; t < D; ++t)
                dst[t] = v[t] / norm;
            ++n_basis;
        }
        if (n_basis >= n_max)
            break;
    }

    free(v);
    free(weights);
    free(perm_inv);
    return n_basis;
}

int irrep_sg_little_group_point_order(const irrep_sg_little_group_t *lg) {
    return lg ? lg->point_order : 0;
}

int irrep_sg_little_group_order(const irrep_sg_little_group_t *lg) {
    return lg ? lg->point_order * lg->num_translations : 0;
}

void irrep_sg_little_group_point_ops(const irrep_sg_little_group_t *lg, int *out_indices) {
    if (!lg || !out_indices)
        return;
    memcpy(out_indices, lg->point_ops, (size_t)lg->point_order * sizeof(int));
}

void irrep_sg_little_group_k(const irrep_sg_little_group_t *lg, int *out_kx, int *out_ky) {
    if (!lg)
        return;
    if (out_kx)
        *out_kx = lg->kx;
    if (out_ky)
        *out_ky = lg->ky;
}

const irrep_space_group_t *irrep_sg_little_group_parent(const irrep_sg_little_group_t *lg) {
    return lg ? lg->G : NULL;
}
