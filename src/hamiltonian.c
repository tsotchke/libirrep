/* SPDX-License-Identifier: MIT */
/* Spin-½ Hamiltonian apply operators.
 *
 * Internal representation: a flat bond list with *per-bond* diagonal
 * (S^z S^z) and off-diagonal (½·(S^+ S^- + h.c.)) couplings. Every
 * exchange-type model in this module is a specialisation:
 *
 *   Heisenberg:  coeff_zz[b] = J/4,   coeff_pm[b] = J/2
 *   XY:          coeff_zz[b] = 0,     coeff_pm[b] = J/2
 *   J₁–J₂:       coeff_zz / coeff_pm vary between NN and NNN segments
 *
 * Same apply loop for all three — data-driven, no branching on model
 * type inside the hot path. This is the canonical payoff of promoting
 * the ad-hoc apply_H_heisenberg to a library primitive: extending the
 * model surface does not cost a new apply, just a new constructor.
 *
 * Hot loop (irrep_heisenberg_apply):
 *   for each bond b, endpoints (i, j):
 *     for each basis state s:
 *       zz_sign = (bit_i ^ bit_j) ? −1 : +1
 *       out[s]             += coeff_zz[b] · zz_sign · psi[s]          (diag)
 *       if bit_i != bit_j and coeff_pm[b] != 0:
 *         out[s ⊕ (1<<i) ⊕ (1<<j)] += coeff_pm[b] · psi[s]             (flip)
 */

#include <complex.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/hamiltonian.h>

extern void irrep_set_error_(const char *fmt, ...);

struct irrep_heisenberg {
    int     num_sites;
    int     num_bonds;
    int    *bi;
    int    *bj;
    double *coeff_zz;        /* per-bond S^z_i S^z_j coefficient  */
    double *coeff_pm;        /* per-bond ½(S^+_i S^−_j + h.c.)    */
};

/* Shared allocator. Copies the bond arrays + per-bond coefficient
 * arrays. The caller owns the input arrays. */
static irrep_heisenberg_t *
build_(int num_sites, int num_bonds,
       const int *bi, const int *bj,
       const double *coeff_zz, const double *coeff_pm) {
    if (num_sites < 1 || num_sites > 62 || num_bonds < 0 || !bi || !bj) {
        irrep_set_error_("irrep_heisenberg build: invalid arguments");
        return NULL;
    }
    for (int b = 0; b < num_bonds; ++b) {
        if (bi[b] < 0 || bi[b] >= num_sites ||
            bj[b] < 0 || bj[b] >= num_sites ||
            bi[b] == bj[b]) {
            irrep_set_error_("irrep_heisenberg build: bond %d out of range", b);
            return NULL;
        }
    }
    irrep_heisenberg_t *H = calloc(1, sizeof *H);
    if (!H) { irrep_set_error_("irrep_heisenberg build: OOM"); return NULL; }
    H->num_sites = num_sites;
    H->num_bonds = num_bonds;
    if (num_bonds > 0) {
        H->bi       = malloc((size_t)num_bonds * sizeof(int));
        H->bj       = malloc((size_t)num_bonds * sizeof(int));
        H->coeff_zz = malloc((size_t)num_bonds * sizeof(double));
        H->coeff_pm = malloc((size_t)num_bonds * sizeof(double));
        if (!H->bi || !H->bj || !H->coeff_zz || !H->coeff_pm) {
            free(H->bi); free(H->bj); free(H->coeff_zz); free(H->coeff_pm);
            free(H);
            irrep_set_error_("irrep_heisenberg build: OOM (bond arrays)");
            return NULL;
        }
        memcpy(H->bi, bi, (size_t)num_bonds * sizeof(int));
        memcpy(H->bj, bj, (size_t)num_bonds * sizeof(int));
        memcpy(H->coeff_zz, coeff_zz, (size_t)num_bonds * sizeof(double));
        memcpy(H->coeff_pm, coeff_pm, (size_t)num_bonds * sizeof(double));
    }
    return H;
}

irrep_heisenberg_t *
irrep_heisenberg_new(int num_sites, int num_bonds,
                     const int *bi, const int *bj,
                     double J) {
    if (num_bonds < 0) num_bonds = 0;
    double *czz = malloc((size_t)(num_bonds ? num_bonds : 1) * sizeof(double));
    double *cpm = malloc((size_t)(num_bonds ? num_bonds : 1) * sizeof(double));
    if (!czz || !cpm) { free(czz); free(cpm); return NULL; }
    for (int b = 0; b < num_bonds; ++b) {
        czz[b] = 0.25 * J;
        cpm[b] = 0.5  * J;
    }
    irrep_heisenberg_t *H = build_(num_sites, num_bonds, bi, bj, czz, cpm);
    free(czz); free(cpm);
    return H;
}

irrep_heisenberg_t *
irrep_heisenberg_j1j2_new(int num_sites,
                          int num_bonds_nn,  const int *nn_i,  const int *nn_j,  double J1,
                          int num_bonds_nnn, const int *nnn_i, const int *nnn_j, double J2) {
    if (num_bonds_nn < 0 || num_bonds_nnn < 0) {
        irrep_set_error_("irrep_heisenberg_j1j2_new: negative bond count");
        return NULL;
    }
    int total = num_bonds_nn + num_bonds_nnn;
    int    *bi  = malloc((size_t)(total ? total : 1) * sizeof(int));
    int    *bj  = malloc((size_t)(total ? total : 1) * sizeof(int));
    double *czz = malloc((size_t)(total ? total : 1) * sizeof(double));
    double *cpm = malloc((size_t)(total ? total : 1) * sizeof(double));
    if (!bi || !bj || !czz || !cpm) {
        free(bi); free(bj); free(czz); free(cpm);
        irrep_set_error_("irrep_heisenberg_j1j2_new: OOM");
        return NULL;
    }
    for (int b = 0; b < num_bonds_nn; ++b) {
        bi[b]  = nn_i[b]; bj[b] = nn_j[b];
        czz[b] = 0.25 * J1;
        cpm[b] = 0.5  * J1;
    }
    for (int b = 0; b < num_bonds_nnn; ++b) {
        bi[num_bonds_nn + b]  = nnn_i[b];  bj[num_bonds_nn + b] = nnn_j[b];
        czz[num_bonds_nn + b] = 0.25 * J2;
        cpm[num_bonds_nn + b] = 0.5  * J2;
    }
    irrep_heisenberg_t *H = build_(num_sites, total, bi, bj, czz, cpm);
    free(bi); free(bj); free(czz); free(cpm);
    return H;
}

irrep_heisenberg_t *
irrep_xy_new(int num_sites, int num_bonds,
             const int *bi, const int *bj,
             double J) {
    if (num_bonds < 0) num_bonds = 0;
    double *czz = malloc((size_t)(num_bonds ? num_bonds : 1) * sizeof(double));
    double *cpm = malloc((size_t)(num_bonds ? num_bonds : 1) * sizeof(double));
    if (!czz || !cpm) { free(czz); free(cpm); return NULL; }
    for (int b = 0; b < num_bonds; ++b) {
        czz[b] = 0.0;                   /* XY has no S^z S^z term */
        cpm[b] = 0.5 * J;
    }
    irrep_heisenberg_t *H = build_(num_sites, num_bonds, bi, bj, czz, cpm);
    free(czz); free(cpm);
    return H;
}

void irrep_heisenberg_free(irrep_heisenberg_t *H) {
    if (!H) return;
    free(H->bi);
    free(H->bj);
    free(H->coeff_zz);
    free(H->coeff_pm);
    free(H);
}

int irrep_heisenberg_num_sites(const irrep_heisenberg_t *H) {
    return H ? H->num_sites : 0;
}

long long irrep_heisenberg_dim(const irrep_heisenberg_t *H) {
    return H ? (1LL << H->num_sites) : 0;
}

void irrep_heisenberg_apply(const double _Complex *psi,
                            double _Complex       *out,
                            void                  *opaque) {
    if (!psi || !out || !opaque) return;
    const irrep_heisenberg_t *H = (const irrep_heisenberg_t *)opaque;
    const long long dim = 1LL << H->num_sites;

    memset(out, 0, (size_t)dim * sizeof(double _Complex));
    for (int b = 0; b < H->num_bonds; ++b) {
        const int    i   = H->bi[b], j = H->bj[b];
        const long long mi = 1LL << i, mj = 1LL << j;
        const double czz = H->coeff_zz[b];
        const double cpm = H->coeff_pm[b];
        for (long long s = 0; s < dim; ++s) {
            const int zi = (int)((s >> i) & 1);
            const int zj = (int)((s >> j) & 1);
            const double zz_sign = (zi ^ zj) ? -1.0 : +1.0;
            out[s] += (czz * zz_sign) * psi[s];
            if (zi != zj && cpm != 0.0) {
                const long long s_flip = s ^ mi ^ mj;
                out[s_flip] += cpm * psi[s];
            }
        }
    }
}
