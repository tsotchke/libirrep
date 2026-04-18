/* SPDX-License-Identifier: MIT */
/* M6: time-reversal operator.
 *
 * T is antiunitary: T (α |ψ⟩) = α* T |ψ⟩.  The library exposes the *linear*
 * factor U so that T ψ = U · ψ*, where ψ* is the component-wise complex
 * conjugate of the input vector. Using the standard convention
 *
 *     T |j, m⟩ = (−1)^{j − m} |j, −m⟩,
 *
 * the matrix elements work out to
 *
 *     U_{m', m} = (−1)^{j + m'} δ_{m, −m'}.
 *
 * Under T² we get (U U*)_{m'', m} = (−1)^{2j} δ_{m, m''}, so T² = +1 for
 * integer j (O(3) irreps) and T² = −1 for half-integer j (Kramers pairs).
 */

#include <complex.h>

#include <irrep/time_reversal.h>
#include <irrep/types.h>

#define UNUSED(x) ((void)(x))

/* helper: compute (−1)^n safely for any int (two's-complement AND works) */
static inline double sign_from_(int n) {
    return (n & 1) ? -1.0 : 1.0;
}

void irrep_time_reversal_integer(int l, double _Complex *out) {
    if (l < 0 || !out) return;
    int d = 2 * l + 1;
    for (int i = 0; i < d * d; ++i) out[i] = 0.0;
    for (int i_mp = 0; i_mp < d; ++i_mp) {
        int mp   = i_mp - l;
        int i_m  = 2 * l - i_mp;   /* m = −m' → i_m = l + m = l − mp = 2l − i_mp */
        out[i_mp * d + i_m] = sign_from_(l + mp);
    }
}

void irrep_time_reversal_half_integer(int two_j, double _Complex *out) {
    if (two_j < 0 || !out) return;
    int d = two_j + 1;
    for (int i = 0; i < d * d; ++i) out[i] = 0.0;
    for (int i_mp = 0; i_mp < d; ++i_mp) {
        int two_mp = 2 * i_mp - two_j;
        int two_m  = -two_mp;
        int i_m    = (two_m + two_j) / 2;
        int phase  = (two_j + two_mp) / 2;   /* always integer by parity */
        out[i_mp * d + i_m] = sign_from_(phase);
    }
}

void irrep_time_reversal_multiset(const irrep_multiset_t *m, double _Complex *out) {
    if (!m || !out) return;
    int total = m->total_dim;
    for (int i = 0; i < total * total; ++i) out[i] = 0.0;

    double _Complex block[(2 * IRREP_L_MAX + 1) * (2 * IRREP_L_MAX + 1)];
    int offset = 0;
    for (int t = 0; t < m->num_terms; ++t) {
        int l    = m->labels[t].l;
        int mult = m->multiplicities[t];
        int d    = 2 * l + 1;
        irrep_time_reversal_integer(l, block);
        for (int c = 0; c < mult; ++c) {
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    out[(offset + i) * total + (offset + j)] = block[i * d + j];
                }
            }
            offset += d;
        }
    }
}

int irrep_time_reversal_square_sign(const irrep_multiset_t *m) {
    if (!m || m->num_terms == 0) return 0;
    /* irrep_multiset_t carries only integer-l labels, so T² = +1 always.
     * Half-integer multisets require a forthcoming irrep_multiset_2j_t type. */
    return +1;
}
