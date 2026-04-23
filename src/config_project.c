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

void irrep_sg_little_group_irrep_free(irrep_sg_little_group_irrep_t *mu_k) {
    if (!mu_k)
        return;
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

irrep_sg_little_group_irrep_t *
irrep_sg_little_group_irrep_named(const irrep_sg_little_group_t *lg, irrep_lg_named_irrep_t name) {
    if (!lg) {
        irrep_set_error_("irrep_sg_little_group_irrep_named: NULL little group");
        return NULL;
    }
    int n = lg->point_order;
    if (n != 12 && n != 6) {
        irrep_set_error_("irrep_sg_little_group_irrep_named: only C_6v (order 12) and C_3v "
                         "(order 6) supported; got order %d",
                         n);
        return NULL;
    }

    /* Validate name vs little-group shape and pick the character table. */
    int             idx_in_table;
    int             dim;
    const double   *table_row;
    int             classes; /* for the table's column count */
    const double (*full_table)[6];
    const double (*c3v_table)[3];

    if (n == 12) {
        /* C_6v. */
        switch (name) {
            case IRREP_LG_IRREP_A1: idx_in_table = 0; break;
            case IRREP_LG_IRREP_A2: idx_in_table = 1; break;
            case IRREP_LG_IRREP_B1: idx_in_table = 2; break;
            case IRREP_LG_IRREP_B2: idx_in_table = 3; break;
            case IRREP_LG_IRREP_E1: idx_in_table = 4; break;
            case IRREP_LG_IRREP_E2: idx_in_table = 5; break;
            default:
                irrep_set_error_(
                    "irrep_sg_little_group_irrep_named: named irrep not valid on C_6v");
                return NULL;
        }
        dim        = c6v_dim_[idx_in_table];
        full_table = c6v_chi_;
        table_row  = c6v_chi_[idx_in_table];
        classes    = 6;
        c3v_table  = NULL;
        (void)table_row;
        (void)c3v_table;
    } else {
        /* n == 6; C_3v. */
        switch (name) {
            case IRREP_LG_IRREP_A1: idx_in_table = 0; break;
            case IRREP_LG_IRREP_A2: idx_in_table = 1; break;
            case IRREP_LG_IRREP_E:  idx_in_table = 2; break;
            default:
                irrep_set_error_(
                    "irrep_sg_little_group_irrep_named: named irrep not valid on C_3v");
                return NULL;
        }
        dim       = c3v_dim_[idx_in_table];
        full_table = NULL;
        c3v_table  = c3v_chi_;
        classes    = 3;
        (void)full_table;
    }

    /* Build the per-element character row. */
    double _Complex *chi = malloc((size_t)n * sizeof(double _Complex));
    if (!chi) {
        irrep_set_error_("irrep_sg_little_group_irrep_named: out of memory");
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        int parent = lg->point_ops[i];
        if (parent < 0 || parent >= 12) {
            free(chi);
            irrep_set_error_("irrep_sg_little_group_irrep_named: parent op out of range");
            return NULL;
        }
        int p6mm_cls = p6mm_point_class_[parent];
        int lg_cls;
        if (n == 12) {
            /* C_6v: parent class maps 1-to-1 to table column. */
            lg_cls = p6mm_cls;
        } else {
            /* C_3v at K. Classes {E=0, C_3=1, σ=2}. The 6-element little
             * group contains parent classes {E, C_3, σ_v OR σ_d}. The
             * three rotations (including E) map to C_3v classes 0 / 1;
             * the three mirrors map to class 2. */
            if (p6mm_cls == 0)
                lg_cls = 0; /* E */
            else if (p6mm_cls == 2)
                lg_cls = 1; /* C_3 or C_3² */
            else if (p6mm_cls == 4 || p6mm_cls == 5)
                lg_cls = 2; /* mirror */
            else {
                free(chi);
                irrep_set_error_(
                    "irrep_sg_little_group_irrep_named: unexpected element at C_3v slot %d (parent %d, "
                    "p6mm class %d)",
                    i, parent, p6mm_cls);
                return NULL;
            }
        }
        if (lg_cls < 0 || lg_cls >= classes) {
            free(chi);
            irrep_set_error_(
                "irrep_sg_little_group_irrep_named: internal class index out of range");
            return NULL;
        }
        double v = (n == 12) ? full_table[idx_in_table][lg_cls] : c3v_table[idx_in_table][lg_cls];
        chi[i]   = v + 0.0 * I;
    }

    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi, dim);
    free(chi);
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
