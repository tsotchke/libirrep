/* SPDX-License-Identifier: MIT */
/* Spin-½ Dzyaloshinskii-Moriya Hamiltonian apply operator.
 *
 * Algorithm: for each bond, the DMI term `D · (S_a × S_b)` decomposes
 * into three pair-flip operators (see dmi_hamiltonian.h for the full
 * derivation). Each operator flips at most two bits with a known
 * complex coefficient depending on (D_x, D_y, D_z) and the bit
 * pattern of the source state. */

#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/dmi_hamiltonian.h>

extern void irrep_set_error_(const char *fmt, ...);

struct irrep_dmi_hamiltonian {
    int      num_sites;
    int      num_bonds;
    int     *bi;
    int     *bj;
    double  *D_x;
    double  *D_y;
    double  *D_z;
    long long dim;
};

irrep_dmi_hamiltonian_t *irrep_dmi_hamiltonian_new(int num_sites, int num_bonds, const int *bi,
                                                   const int *bj, const double *D_x,
                                                   const double *D_y, const double *D_z) {
    if (num_sites <= 0 || num_sites > 30 || num_bonds < 0 || (num_bonds > 0 && (!bi || !bj))) {
        irrep_set_error_("irrep_dmi_hamiltonian_new: invalid arguments");
        return NULL;
    }
    irrep_dmi_hamiltonian_t *H = calloc(1, sizeof(*H));
    if (!H) {
        irrep_set_error_("irrep_dmi_hamiltonian_new: out of memory");
        return NULL;
    }
    H->num_sites = num_sites;
    H->num_bonds = num_bonds;
    H->dim = 1LL << num_sites;
    if (num_bonds > 0) {
        H->bi = malloc((size_t)num_bonds * sizeof(int));
        H->bj = malloc((size_t)num_bonds * sizeof(int));
        H->D_x = malloc((size_t)num_bonds * sizeof(double));
        H->D_y = malloc((size_t)num_bonds * sizeof(double));
        H->D_z = malloc((size_t)num_bonds * sizeof(double));
        if (!H->bi || !H->bj || !H->D_x || !H->D_y || !H->D_z) {
            irrep_dmi_hamiltonian_free(H);
            irrep_set_error_("irrep_dmi_hamiltonian_new: out of memory");
            return NULL;
        }
        memcpy(H->bi, bi, (size_t)num_bonds * sizeof(int));
        memcpy(H->bj, bj, (size_t)num_bonds * sizeof(int));
        if (D_x)
            memcpy(H->D_x, D_x, (size_t)num_bonds * sizeof(double));
        else
            for (int b = 0; b < num_bonds; ++b)
                H->D_x[b] = 0;
        if (D_y)
            memcpy(H->D_y, D_y, (size_t)num_bonds * sizeof(double));
        else
            for (int b = 0; b < num_bonds; ++b)
                H->D_y[b] = 0;
        if (D_z)
            memcpy(H->D_z, D_z, (size_t)num_bonds * sizeof(double));
        else
            for (int b = 0; b < num_bonds; ++b)
                H->D_z[b] = 0;
    }
    return H;
}

void irrep_dmi_hamiltonian_free(irrep_dmi_hamiltonian_t *H) {
    if (!H)
        return;
    free(H->bi);
    free(H->bj);
    free(H->D_x);
    free(H->D_y);
    free(H->D_z);
    free(H);
}

int irrep_dmi_hamiltonian_num_sites(const irrep_dmi_hamiltonian_t *H) {
    return H ? H->num_sites : 0;
}

long long irrep_dmi_hamiltonian_dim(const irrep_dmi_hamiltonian_t *H) {
    return H ? H->dim : 0;
}

void irrep_dmi_apply(const double _Complex *psi, double _Complex *out, void *opaque) {
    const irrep_dmi_hamiltonian_t *H = opaque;
    long long dim = H->dim;
    memset(out, 0, (size_t)dim * sizeof *out);

    /* Spin-½ matrix elements of D · (S_a × S_b):
     *
     *   coeff_pp = ½ (D_x − i D_y)   on   |↑_a ↓_b⟩ → contributes
     *                                       z-component flips on b-only
     *                                     (S_a^z S_b^+ − S_a^+ S_b^z)
     *
     * Working it out term-by-term in the |bit_a bit_b⟩ basis (other bits
     * untouched). Let s_a = bit a of s, s_b = bit b of s. Define
     * z_a = ½ if s_a = 1 else −½; same for z_b.
     *
     * Term 1: D_x (S_a^y S_b^z − S_a^z S_b^y)
     *   S_a^y |s⟩ flips bit a, multiplies by ±i/2 (sign by old s_a)
     *   Combined: D_x · (i/2) · z_b · (1 − 2 s_a) [bit-a flip]
     *           − D_x · z_a · (i/2) · (1 − 2 s_b) [bit-b flip]
     *
     * Term 2: D_y (S_a^z S_b^x − S_a^x S_b^z)
     *   D_y · z_a · (½) · [bit-b flip]
     *   − D_y · (½) · z_b · [bit-a flip]
     *
     * Term 3: D_z (S_a^x S_b^y − S_a^y S_b^x)
     *   Both x and y are bit-flips. Together:
     *   S^x_a S^y_b − S^y_a S^x_b on basis |s_a s_b⟩:
     *      = (1/4)[(σ^a_+ + σ^a_−)(σ^b_+ − σ^b_−)/i − (σ^a_+ − σ^a_−)/i (σ^b_+ + σ^b_−)]
     *   This flips both bits a and b. Coefficient depends on (s_a, s_b):
     *     |↓↓⟩ → only S_+S_+ term — 0 (S_+ S_+ raises both, both already
     *           up; (1/4)[1·(−1)/i − 1/i·1] = (1/4)(−2/i) = i/2)
     *           Actually more carefully:
     *           S^a_+ S^b_+ on |↓↓⟩ = |↑↑⟩
     *           Coefficient: (1/4)[(σ^a_+)(σ^b_+ − σ^b_-)/i + (σ^a_- − σ^a_+)/i (σ^b_+) ]
     *                      = (1/4)[ (σ^a_+ σ^b_+)/i − (σ^a_+ σ^b_-)/i − (σ^a_+ σ^b_+)/i + ... ]
     *   This algebra gets messy. The cleaner way: enumerate the 4 cases
     *   explicitly (↓↓, ↓↑, ↑↓, ↑↑) and tabulate the matrix elements. */

    for (int b = 0; b < H->num_bonds; ++b) {
        int      i = H->bi[b], j = H->bj[b];
        long long mi = 1LL << i, mj = 1LL << j;
        double Dx = H->D_x[b], Dy = H->D_y[b], Dz = H->D_z[b];

        for (long long s = 0; s < dim; ++s) {
            int sa = (int)((s >> i) & 1);  /* bit of site a */
            int sb = (int)((s >> j) & 1);
            double za = sa ? +0.5 : -0.5;
            double zb = sb ? +0.5 : -0.5;
            double sgn_a = sa ? -1.0 : +1.0;  /* sign for S^y on bit-a flip */
            double sgn_b = sb ? -1.0 : +1.0;

            /* Term 1 (D_x): two bit-a-only flips and bit-b-only flips */
            /* Bit-a flip: coefficient on |s⟩ → |s ⊕ mi⟩
             *   D_x · (i/2)·z_b·sgn_a    contribution to out[s ⊕ mi]
             */
            out[s ^ mi] += (Dx * 0.5 * I * sgn_a * zb) * psi[s];
            out[s ^ mj] += (-Dx * 0.5 * I * sgn_b * za) * psi[s];

            /* Term 2 (D_y): bit-b flip and bit-a flip */
            out[s ^ mj] += (Dy * 0.5 * za) * psi[s];
            out[s ^ mi] += (-Dy * 0.5 * zb) * psi[s];

            /* Term 3 (D_z): bit-a AND bit-b both flip; only nonzero
             * when bits are anti-aligned (S_+ S_- or S_- S_+).
             *   D_z · (i/2) · [S_+_a S_-_b − S_-_a S_+_b]
             *   On |↑↓⟩: S_- gives ½ on a, S_+ gives ½ on b. The
             *     S_+_a S_-_b term acts as |↓↑⟩ × ½ × ½. Hmm wait,
             *     S_+_a |↑⟩ = 0, so on |↑↓⟩ only S_-_a S_+_b is nonzero.
             *   Let me just enumerate: */
            if (sa != sb) {
                long long sf = s ^ mi ^ mj;
                /* On |↑_a ↓_b⟩ = (sa=1, sb=0): only S_-_a S_+_b nonzero.
                 * Coeff: (i/2) · [S_+_a S_-_b − S_-_a S_+_b] · |↑↓⟩
                 *      = (i/2) · [0 − ½ · ½] · |↓↑⟩ = (-i/8) |↓↑⟩
                 * On |↓_a ↑_b⟩ = (sa=0, sb=1): only S_+_a S_-_b nonzero.
                 *      = (i/2) · [½ · ½ − 0] · |↑↓⟩ = (+i/8) |↑↓⟩
                 *
                 * S_+ |↓⟩ = |↑⟩ (factor 1, since spin-½ S_+_z gives ℏ√(s(s+1) - m(m+1)) = √(¾ - ¼) = √(½) = ½√2 ... wait that's not 1)
                 *
                 * Convention check: the Pauli operator σ_+ |↓⟩ = |↑⟩ with
                 * coefficient 1, and S_+ = (1/2) σ_+ for spin-½, so
                 * S_+ |↓⟩ = (1/2) |↑⟩. Hmm. So S_+ S_- has coefficient
                 * 1/4 between bit-flipped states. Let me recompute.
                 *
                 *   S_+_a S_-_b |↓_a ↑_b⟩ = S_+_a (½ |↓_a ↓_b⟩) = ½ · ½ |↑_a ↓_b⟩
                 *                          = (1/4) |↑↓⟩
                 *   S_-_a S_+_b |↑_a ↓_b⟩ = (1/4) |↓↑⟩
                 *
                 * Then D_z (i/2) [S_+_a S_-_b − S_-_a S_+_b]:
                 *   On |↓_a ↑_b⟩:  D_z (i/2) (1/4) |↑↓⟩ = (i D_z / 8) |↑↓⟩
                 *   On |↑_a ↓_b⟩:  −D_z (i/2) (1/4) |↓↑⟩ = (−i D_z / 8) |↓↑⟩
                 */
                double sign = (sa == 0 && sb == 1) ? +1.0 : -1.0;
                out[sf] += (Dz * 0.25 * I * sign) * psi[s];
            }
        }
    }
}
