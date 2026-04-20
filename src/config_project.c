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
#  define M_PI 3.14159265358979323846
#endif

extern void irrep_set_error_(const char *fmt, ...);

struct irrep_sg_irrep {
    const irrep_space_group_t *G;        /* borrowed */
    int                        order;
    int                        irrep_dim;
    double _Complex           *characters;
};

irrep_sg_irrep_t *irrep_sg_irrep_new(const irrep_space_group_t *G,
                                     const double _Complex *characters,
                                     int irrep_dim) {
    if (!G || !characters || irrep_dim <= 0) {
        irrep_set_error_("irrep_sg_irrep_new: invalid arguments");
        return NULL;
    }
    int order = irrep_space_group_order(G);

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
    mu->G         = G;
    mu->order     = order;
    mu->irrep_dim = irrep_dim;
    return mu;
}

void irrep_sg_irrep_free(irrep_sg_irrep_t *mu) {
    if (!mu) return;
    free(mu->characters);
    free(mu);
}

irrep_sg_irrep_t *irrep_sg_trivial(const irrep_space_group_t *G) {
    if (!G) return NULL;
    int order = irrep_space_group_order(G);
    double _Complex *chi = malloc((size_t)order * sizeof(double _Complex));
    if (!chi) {
        irrep_set_error_("irrep_sg_trivial: out of memory");
        return NULL;
    }
    for (int g = 0; g < order; ++g) chi[g] = 1.0 + 0.0*I;
    irrep_sg_irrep_t *mu = irrep_sg_irrep_new(G, chi, 1);
    free(chi);
    return mu;
}

irrep_sg_irrep_t *irrep_sg_sign_rep(const irrep_space_group_t *G) {
    if (!G) return NULL;
    int order       = irrep_space_group_order(G);
    int point_order = irrep_space_group_point_order(G);
    /* p4mm / p6mm: the builder enumerates rotations first, then mirrors —
     * so index p ∈ [0, point_order/2) is a rotation and p ∈ [point_order/2,
     * point_order) is a reflection. For p1 (point_order = 1) every element
     * is a translation, sign +1 throughout. Full element index is
     * (cell_index · point_order + p); characters are constant across cells
     * at the Γ-point.                                                      */
    int half = point_order / 2;
    double _Complex *chi = malloc((size_t)order * sizeof(double _Complex));
    if (!chi) {
        irrep_set_error_("irrep_sg_sign_rep: out of memory");
        return NULL;
    }
    for (int g = 0; g < order; ++g) {
        int p    = g % point_order;
        int sign = (point_order > 1 && p >= half) ? -1 : +1;
        chi[g] = (double)sign + 0.0*I;
    }
    irrep_sg_irrep_t *mu = irrep_sg_irrep_new(G, chi, 1);
    free(chi);
    return mu;
}

double _Complex
irrep_sg_project_amplitude(const irrep_sg_irrep_t *mu,
                           const double _Complex *psi_of_g) {
    if (!mu || !psi_of_g) return 0.0 + 0.0*I;
    double _Complex acc = 0.0 + 0.0*I;
    for (int g = 0; g < mu->order; ++g) {
        acc += conj(mu->characters[g]) * psi_of_g[g];
    }
    return acc * ((double)mu->irrep_dim / (double)mu->order);
}

double _Complex
irrep_sg_project_A1(const irrep_space_group_t *G,
                    const double _Complex *psi_of_g) {
    if (!G || !psi_of_g) return 0.0 + 0.0*I;
    int order = irrep_space_group_order(G);
    double _Complex acc = 0.0 + 0.0*I;
    for (int g = 0; g < order; ++g) acc += psi_of_g[g];
    return acc / (double)order;
}

void irrep_sg_enumerate_orbit(const irrep_space_group_t *G,
                              const double *sigma,
                              double *out_orbit) {
    if (!G || !sigma || !out_orbit) return;
    int n     = irrep_space_group_num_sites(G);
    int order = irrep_space_group_order(G);
    int *inv  = malloc((size_t)n * sizeof(int));
    if (!inv) return;
    for (int g = 0; g < order; ++g) {
        irrep_space_group_permutation_inverse(G, g, inv);
        double *dst = out_orbit + (size_t)g * n;
        for (int s = 0; s < n; ++s) dst[s] = sigma[inv[s]];
    }
    free(inv);
}

/* -------------------------------------------------------------------------- *
 * Symmetry-adapted basis builder                                             *
 * -------------------------------------------------------------------------- */

static long long ipow_ll_(int base, int exp) {
    long long r = 1;
    for (int k = 0; k < exp; ++k) r *= (long long)base;
    return r;
}

/* Apply a site permutation to the digit-encoded index `s`: bit/digit j of
 * the output = bit/digit perm[j] of the input. */
static long long permute_digits_(long long s, int Nsites, int d,
                                 const int *perm) {
    if (d == 2) {
        long long out = 0;
        for (int j = 0; j < Nsites; ++j) {
            long long bit = (s >> perm[j]) & 1LL;
            out |= bit << j;
        }
        return out;
    }
    /* Generic base-d case. Extract digits once, then reassemble. */
    int digits[64];
    long long tmp = s;
    for (int j = 0; j < Nsites; ++j) { digits[j] = (int)(tmp % d); tmp /= d; }
    long long out = 0, w = 1;
    for (int j = 0; j < Nsites; ++j) { out += (long long)digits[perm[j]] * w; w *= d; }
    return out;
}

int irrep_sg_adapted_basis(const irrep_space_group_t *G,
                           const irrep_sg_irrep_t *mu,
                           int num_sites, int local_dim,
                           double _Complex *basis_out,
                           int n_max) {
    if (!G || !mu || !basis_out || n_max <= 0 || num_sites < 1 || local_dim < 2)
        return -1;
    if (num_sites != irrep_space_group_num_sites(G)) return -1;

    int order = irrep_space_group_order(G);
    if (mu->order != order) return -1;

    long long D = ipow_ll_(local_dim, num_sites);
    if (D <= 0) return -1;

    /* Cache every group element's inverse permutation once, upfront. Turns
     * the inner loop's permutation lookup into a pointer dereference rather
     * than a memcpy via irrep_space_group_permutation_inverse. */
    int *perm_inv_table = malloc((size_t)order * (size_t)num_sites * sizeof(int));
    if (!perm_inv_table) return -1;
    for (int g = 0; g < order; ++g) {
        irrep_space_group_permutation_inverse(G, g,
                                              perm_inv_table + (size_t)g * num_sites);
    }

    double _Complex *v = malloc((size_t)D * sizeof(double _Complex));
    if (!v) { free(perm_inv_table); return -1; }

    const double tol = 1e-9;
    int n_basis = 0;
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
            long long s_g = permute_digits_(s, num_sites, local_dim, perm);
            if (s_g < s) is_min = 0;
        }
        if (!is_min) continue;

        /* v = P_μ |s⟩ = (d_μ / |G|) Σ_g conj(χ(g)) · g|s⟩ */
        memset(v, 0, (size_t)D * sizeof(double _Complex));
        for (int g = 0; g < order; ++g) {
            const int *perm = perm_inv_table + (size_t)g * num_sites;
            long long s_g = permute_digits_(s, num_sites, local_dim, perm);
            v[s_g] += scale * conj(mu->characters[g]);
        }

        /* Orthogonalise against accepted basis. */
        for (int k = 0; k < n_basis; ++k) {
            const double _Complex *b_k = basis_out + (size_t)k * D;
            double _Complex overlap = 0.0;
            for (long long t = 0; t < D; ++t) overlap += conj(b_k[t]) * v[t];
            for (long long t = 0; t < D; ++t) v[t] -= overlap * b_k[t];
        }

        double norm2 = 0.0;
        for (long long t = 0; t < D; ++t)
            norm2 += creal(v[t]) * creal(v[t]) + cimag(v[t]) * cimag(v[t]);
        double norm = sqrt(norm2);

        if (norm > tol && n_basis < n_max) {
            double _Complex *dst = basis_out + (size_t)n_basis * D;
            for (long long t = 0; t < D; ++t) dst[t] = v[t] / norm;
            ++n_basis;
        }
        if (n_basis >= n_max) break;
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

double _Complex
irrep_sg_bloch_amplitude(const irrep_space_group_t *G,
                         int kx, int ky,
                         const double _Complex *psi_of_g) {
    if (!G || !psi_of_g) return 0.0 + 0.0*I;
    const irrep_lattice_t *L = irrep_space_group_lattice(G);
    if (!L) return 0.0 + 0.0*I;
    int Lx = irrep_lattice_Lx(L);
    int Ly = irrep_lattice_Ly(L);
    /* Bloch indices live mod (Lx, Ly). Canonicalise so callers can pass any
     * integer (including negatives for signed BZ conventions) without hitting
     * the orbit-rep filter's [0, Lx) assumption in bloch_basis. */
    kx = ((kx % Lx) + Lx) % Lx;
    ky = ((ky % Ly) + Ly) % Ly;
    int point_order = irrep_space_group_point_order(G);
    double _Complex acc = 0.0 + 0.0*I;
    double inv_Lx = 1.0 / (double)Lx;
    double inv_Ly = 1.0 / (double)Ly;
    for (int ty = 0; ty < Ly; ++ty) {
        for (int tx = 0; tx < Lx; ++tx) {
            int tidx = ty * Lx + tx;
            int g = tidx * point_order + 0;   /* identity point element */
            double phase = -2.0 * M_PI *
                ((double)(kx * tx) * inv_Lx + (double)(ky * ty) * inv_Ly);
            double _Complex w = cos(phase) + I * sin(phase);
            acc += w * psi_of_g[g];
        }
    }
    return acc / (double)(Lx * Ly);
}

int irrep_sg_bloch_basis(const irrep_space_group_t *G,
                         int kx, int ky,
                         int num_sites, int local_dim,
                         double _Complex *basis_out,
                         int n_max) {
    if (!G || !basis_out || n_max <= 0 || num_sites < 1 || local_dim < 2)
        return -1;
    if (num_sites != irrep_space_group_num_sites(G)) return -1;
    const irrep_lattice_t *L = irrep_space_group_lattice(G);
    if (!L) return -1;
    int Lx = irrep_lattice_Lx(L);
    int Ly = irrep_lattice_Ly(L);
    kx = ((kx % Lx) + Lx) % Lx;
    ky = ((ky % Ly) + Ly) % Ly;
    int point_order = irrep_space_group_point_order(G);
    int num_trans = Lx * Ly;

    long long D = ipow_ll_(local_dim, num_sites);
    if (D <= 0) return -1;

    /* Cache translation-only inverse permutations. Group element
     * g = tidx · point_order + 0 has permutation at row (tidx · point_order). */
    int *trans_inv = malloc((size_t)num_trans * (size_t)num_sites * sizeof(int));
    if (!trans_inv) return -1;
    int *scratch = malloc((size_t)num_sites * sizeof(int));
    if (!scratch) { free(trans_inv); return -1; }
    for (int tidx = 0; tidx < num_trans; ++tidx) {
        int g = tidx * point_order + 0;
        irrep_space_group_permutation_inverse(G, g, scratch);
        memcpy(trans_inv + (size_t)tidx * num_sites, scratch,
               (size_t)num_sites * sizeof(int));
    }
    free(scratch);

    /* Precompute Bloch phases for each translation index. */
    double _Complex *phases = malloc((size_t)num_trans * sizeof(double _Complex));
    if (!phases) { free(trans_inv); return -1; }
    double inv_Lx = 1.0 / (double)Lx;
    double inv_Ly = 1.0 / (double)Ly;
    for (int ty = 0; ty < Ly; ++ty) {
        for (int tx = 0; tx < Lx; ++tx) {
            int tidx = ty * Lx + tx;
            double ph = -2.0 * M_PI *
                ((double)(kx * tx) * inv_Lx + (double)(ky * ty) * inv_Ly);
            phases[tidx] = cos(ph) + I * sin(ph);
        }
    }

    double _Complex *v = malloc((size_t)D * sizeof(double _Complex));
    if (!v) { free(phases); free(trans_inv); return -1; }

    const double tol   = 1e-9;
    int   n_basis      = 0;
    double _Complex scale = 1.0 / (double)num_trans;

    for (long long s = 0; s < D; ++s) {
        /* Orbit-representative filter over translations. For non-Γ k the
         * orbit phase structure still makes each orbit produce one vector;
         * keep the minimum translate as seed to avoid duplicated work. */
        int is_min = 1;
        for (int tidx = 1; tidx < num_trans && is_min; ++tidx) {
            const int *perm = trans_inv + (size_t)tidx * num_sites;
            long long s_g = permute_digits_(s, num_sites, local_dim, perm);
            if (s_g < s) is_min = 0;
        }
        if (!is_min) continue;

        memset(v, 0, (size_t)D * sizeof(double _Complex));
        for (int tidx = 0; tidx < num_trans; ++tidx) {
            const int *perm = trans_inv + (size_t)tidx * num_sites;
            long long s_g = permute_digits_(s, num_sites, local_dim, perm);
            v[s_g] += scale * phases[tidx];
        }

        for (int k = 0; k < n_basis; ++k) {
            const double _Complex *b_k = basis_out + (size_t)k * D;
            double _Complex overlap = 0.0;
            for (long long t = 0; t < D; ++t) overlap += conj(b_k[t]) * v[t];
            for (long long t = 0; t < D; ++t) v[t] -= overlap * b_k[t];
        }

        double norm2 = 0.0;
        for (long long t = 0; t < D; ++t)
            norm2 += creal(v[t]) * creal(v[t]) + cimag(v[t]) * cimag(v[t]);
        double norm = sqrt(norm2);

        if (norm > tol && n_basis < n_max) {
            double _Complex *dst = basis_out + (size_t)n_basis * D;
            for (long long t = 0; t < D; ++t) dst[t] = v[t] / norm;
            ++n_basis;
        }
        if (n_basis >= n_max) break;
    }

    free(v);
    free(phases);
    free(trans_inv);
    return n_basis;
}
