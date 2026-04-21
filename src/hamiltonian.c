/* SPDX-License-Identifier: MIT */
/* Spin-½ Heisenberg Hamiltonian apply-operator.
 *
 * Hot loop:
 *   for each bond (i, j):
 *     for each basis state s:
 *       out[s]      += J/4 · (+1 or −1) · psi[s]       (diagonal S^z S^z)
 *       if zi != zj: out[s ⊕ (1<<i) ⊕ (1<<j)] += J/2 · psi[s]
 *
 * Bit-twiddling used: each spin lives in a single bit of s, so flipping
 * two bits swaps an up-down pair to down-up (the off-diagonal
 * ½(S_i^+ S_j^- + h.c.) matrix element).
 *
 * This is a direct port of the ad-hoc apply_H_heisenberg previously
 * inlined in examples/kagome24_ed.c. Making it a library primitive
 * closes the audit item "every ED example re-implements Heisenberg".
 */

#include <complex.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/hamiltonian.h>

extern void irrep_set_error_(const char *fmt, ...);

struct irrep_heisenberg {
    int    num_sites;
    int    num_bonds;
    int   *bi;
    int   *bj;
    double J;
};

irrep_heisenberg_t *
irrep_heisenberg_new(int num_sites, int num_bonds,
                     const int *bi, const int *bj,
                     double J) {
    if (num_sites < 1 || num_sites > 62 || num_bonds < 0 || !bi || !bj) {
        irrep_set_error_("irrep_heisenberg_new: invalid arguments");
        return NULL;
    }
    for (int b = 0; b < num_bonds; ++b) {
        if (bi[b] < 0 || bi[b] >= num_sites ||
            bj[b] < 0 || bj[b] >= num_sites ||
            bi[b] == bj[b]) {
            irrep_set_error_("irrep_heisenberg_new: bond %d out of range", b);
            return NULL;
        }
    }
    irrep_heisenberg_t *H = calloc(1, sizeof *H);
    if (!H) { irrep_set_error_("irrep_heisenberg_new: OOM"); return NULL; }
    H->num_sites = num_sites;
    H->num_bonds = num_bonds;
    H->J         = J;
    if (num_bonds > 0) {
        H->bi = malloc((size_t)num_bonds * sizeof(int));
        H->bj = malloc((size_t)num_bonds * sizeof(int));
        if (!H->bi || !H->bj) {
            free(H->bi); free(H->bj); free(H);
            irrep_set_error_("irrep_heisenberg_new: OOM (bonds)");
            return NULL;
        }
        memcpy(H->bi, bi, (size_t)num_bonds * sizeof(int));
        memcpy(H->bj, bj, (size_t)num_bonds * sizeof(int));
    }
    return H;
}

void irrep_heisenberg_free(irrep_heisenberg_t *H) {
    if (!H) return;
    free(H->bi);
    free(H->bj);
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
    const double    J   = H->J;

    memset(out, 0, (size_t)dim * sizeof(double _Complex));
    for (int b = 0; b < H->num_bonds; ++b) {
        const int i = H->bi[b], j = H->bj[b];
        const long long mi = 1LL << i, mj = 1LL << j;
        for (long long s = 0; s < dim; ++s) {
            const int zi = (int)((s >> i) & 1);
            const int zj = (int)((s >> j) & 1);
            const double sign = (zi ^ zj) ? -1.0 : +1.0;
            out[s] += (J * 0.25 * sign) * psi[s];
            if (zi != zj) {
                const long long s_flip = s ^ mi ^ mj;
                out[s_flip] += (J * 0.5) * psi[s];
            }
        }
    }
}
