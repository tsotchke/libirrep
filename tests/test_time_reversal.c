/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/time_reversal.h>
#include <irrep/multiset.h>
#include <complex.h>
#include <math.h>

static int approx_eq_c(double _Complex a, double _Complex b, double tol) {
    return cabs(a - b) <= tol;
}

int main(void) {
    IRREP_TEST_START("time_reversal");

    const double tol = 1e-12;

    /* -------- integer l = 1: U permutes |m⟩ ↔ |−m⟩ with phase (−1)^{l+m'} -------- */
    {
        int l = 1, d = 3;
        double _Complex U[9];
        irrep_time_reversal_integer(l, U);
        /* Expected:
         *  row m'=-1 (i=0): [0, 0, 1]        (phase = (-1)^{1+(-1)} = +1)
         *  row m'= 0 (i=1): [0, -1, 0]       (phase = (-1)^{1+0}    = -1)
         *  row m'=+1 (i=2): [1, 0, 0]        (phase = (-1)^{1+1}    = +1) */
        IRREP_ASSERT(approx_eq_c(U[0*d + 2],  1.0, tol));
        IRREP_ASSERT(approx_eq_c(U[1*d + 1], -1.0, tol));
        IRREP_ASSERT(approx_eq_c(U[2*d + 0],  1.0, tol));
        /* All other entries zero */
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                int is_antidiag = (i + j == d - 1);
                if (!is_antidiag) {
                    IRREP_ASSERT(cabs(U[i*d + j]) < tol);
                }
            }
        }
    }

    /* -------- T² = +I for integer l (various l) -------- */
    {
        for (int l = 0; l <= 4; ++l) {
            int d = 2*l + 1;
            double _Complex U[9*9];
            irrep_time_reversal_integer(l, U);
            /* T² = U · U* (complex conjugate of U) */
            for (int i = 0; i < d; ++i) {
                for (int k = 0; k < d; ++k) {
                    double _Complex s = 0;
                    for (int n = 0; n < d; ++n) s += U[i*d + n] * conj(U[n*d + k]);
                    double target = (i == k) ? 1.0 : 0.0;
                    IRREP_ASSERT(cabs(s - target) < tol);
                }
            }
        }
    }

    /* -------- half-integer: T² = −I (Kramers) -------- */
    {
        for (int two_j = 1; two_j <= 5; two_j += 2) {
            int d = two_j + 1;
            double _Complex U[6*6];
            irrep_time_reversal_half_integer(two_j, U);
            for (int i = 0; i < d; ++i) {
                for (int k = 0; k < d; ++k) {
                    double _Complex s = 0;
                    for (int n = 0; n < d; ++n) s += U[i*d + n] * conj(U[n*d + k]);
                    double target = (i == k) ? -1.0 : 0.0;
                    IRREP_ASSERT(cabs(s - target) < tol);
                }
            }
        }
    }

    /* -------- spin-1/2: matrix matches i σ_y convention -------- */
    {
        /* T |m=−1/2⟩ = -|+1/2⟩,   T |m=+1/2⟩ = |-1/2⟩. */
        double _Complex U[4];
        irrep_time_reversal_half_integer(1, U);   /* two_j = 1 */
        /* Indexing: i=0 ↔ m=-1/2, i=1 ↔ m=+1/2. */
        IRREP_ASSERT(approx_eq_c(U[0*2 + 1],  1.0, tol));    /* m=+1/2 → |-1/2⟩ (row 0) */
        IRREP_ASSERT(approx_eq_c(U[1*2 + 0], -1.0, tol));    /* m=-1/2 → -|+1/2⟩ (row 1) */
        IRREP_ASSERT(cabs(U[0]) < tol);
        IRREP_ASSERT(cabs(U[3]) < tol);
    }

    /* -------- multiset block-diagonal -------- */
    {
        irrep_multiset_t *m = irrep_multiset_parse("1x1e + 1x2o");
        IRREP_ASSERT(m != NULL);
        int n = m->total_dim;          /* 3 + 5 = 8 */
        IRREP_ASSERT(n == 8);
        double _Complex T[8 * 8];
        irrep_time_reversal_multiset(m, T);

        /* Verify the (l=1) block is at rows/cols 0..2 and matches
         * irrep_time_reversal_integer(1, ...) */
        double _Complex T1[9];
        irrep_time_reversal_integer(1, T1);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                IRREP_ASSERT(cabs(T[i*8 + j] - T1[i*3 + j]) < tol);
            }
        }

        /* Off-block (rows 0..2, cols 3..7) must be zero */
        for (int i = 0; i < 3; ++i) {
            for (int j = 3; j < 8; ++j) {
                IRREP_ASSERT(cabs(T[i*8 + j]) < tol);
            }
        }

        /* T² = +I for this purely integer-l multiset */
        int sign = irrep_time_reversal_square_sign(m);
        IRREP_ASSERT(sign == +1);

        irrep_multiset_free(m);
    }

    return IRREP_TEST_END();
}
