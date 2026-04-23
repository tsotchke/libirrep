/* SPDX-License-Identifier: MIT */
/* Tests for `irrep_wigner_d_matrix_batch`, the vectorised per-β batch.
 *
 * Coverage:
 *   - Bit-exact agreement between the SIMD-dispatched `_batch` and the
 *     scalar per-β loop across j ∈ {0, 1, 2, 3, 4, 6, 10, 20} and a
 *     diverse set of β values (including boundaries β = 0, π).
 *   - Odd n_betas: the final un-paired β falls through to the scalar
 *     path; result must still agree.
 *   - n_betas = 1 (edge case): single-β batch equals a plain
 *     `_matrix` call.
 *   - Equivariance pass-through: the batched output at β = 0 is exactly
 *     the identity matrix at every j.
 */

#include "harness.h"
#include <irrep/wigner_d.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    IRREP_TEST_START("wigner_d_batch");

    const int js[] = {0, 1, 2, 3, 4, 6, 10, 20};
    const double betas_bank[] = {
        0.0, 0.123, 0.5, M_PI / 6.0, M_PI / 4.0, M_PI / 3.0,
        1.1, 1.5707963, 2.0, 2.5, 3.0, M_PI,
    };
    const size_t n_bank = sizeof(betas_bank) / sizeof(betas_bank[0]);

    /* ---- Bit-exact SIMD batch vs scalar per-β loop on N β values ------- */
    for (size_t ij = 0; ij < sizeof(js) / sizeof(js[0]); ++ij) {
        int j = js[ij];
        int d = 2 * j + 1;
        size_t stride = (size_t)d * (size_t)d;
        size_t n_betas = n_bank;
        double *batched = malloc(n_betas * stride * sizeof(double));
        double *scalar  = malloc(n_betas * stride * sizeof(double));
        IRREP_ASSERT(batched != NULL && scalar != NULL);

        irrep_wigner_d_matrix_batch(j, n_betas, betas_bank, batched);
        for (size_t b = 0; b < n_betas; ++b)
            irrep_wigner_d_matrix(j, scalar + b * stride, betas_bank[b]);

        /* Bit-exact agreement. The SIMD kernel runs the same arithmetic
         * in the same order lane-wise; pragma FP_CONTRACT OFF on both
         * sides keeps fma-fusion out of the picture. Tolerance 1e-14
         * to absorb any residual lane-scheduling reassociation without
         * hiding real regressions. */
        for (size_t k = 0; k < n_betas * stride; ++k) {
            double diff = fabs(batched[k] - scalar[k]);
            IRREP_ASSERT(diff < 1e-14);
        }

        free(batched);
        free(scalar);
    }

    /* ---- Odd n_betas: tail element through the scalar path ------------- */
    {
        int    j = 5;
        int    d = 11;
        size_t stride = 121;
        double betas[5] = {0.3, 0.7, 1.1, 1.5, 1.9};
        double batched[5 * 121];
        double scalar[5 * 121];
        irrep_wigner_d_matrix_batch(j, 5, betas, batched);
        for (int b = 0; b < 5; ++b)
            irrep_wigner_d_matrix(j, scalar + b * stride, betas[b]);
        for (size_t k = 0; k < 5 * stride; ++k)
            IRREP_ASSERT(fabs(batched[k] - scalar[k]) < 1e-14);
    }

    /* ---- n_betas = 1 (edge) equals scalar `_matrix` ------------------- */
    {
        int    j = 4;
        int    d = 9;
        size_t stride = 81;
        double batched[81];
        double scalar[81];
        double beta = 0.777;
        irrep_wigner_d_matrix_batch(j, 1, &beta, batched);
        irrep_wigner_d_matrix(j, scalar, beta);
        for (size_t k = 0; k < stride; ++k)
            IRREP_ASSERT(fabs(batched[k] - scalar[k]) < 1e-14);
    }

    /* ---- β = 0 gives identity for every batched lane ------------------- */
    {
        int    j = 3;
        int    d = 7;
        double betas[4] = {0.0, 0.0, 0.0, 0.0};
        double out[4 * 49];
        irrep_wigner_d_matrix_batch(j, 4, betas, out);
        for (int b = 0; b < 4; ++b) {
            double *m = out + b * 49;
            for (int i = 0; i < d; ++i)
                for (int k = 0; k < d; ++k) {
                    double expected = (i == k) ? 1.0 : 0.0;
                    IRREP_ASSERT_NEAR(m[i * d + k], expected, 1e-14);
                }
        }
    }

    /* ---- Invalid input: NULL buffers / n_betas = 0 are no-ops --------- */
    {
        double dummy = 0.0;
        irrep_wigner_d_matrix_batch(-1, 1, &dummy, &dummy); /* j < 0 */
        irrep_wigner_d_matrix_batch(2, 0, &dummy, &dummy);  /* n_betas = 0 */
        irrep_wigner_d_matrix_batch(2, 1, NULL, &dummy);
        irrep_wigner_d_matrix_batch(2, 1, &dummy, NULL);
        IRREP_ASSERT(1); /* no crash */
    }

    return IRREP_TEST_END();
}
