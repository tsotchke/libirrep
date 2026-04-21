/* SPDX-License-Identifier: MIT */
/* Tests for the primitive equivariant layers used by NequIP-style networks:
 * `irrep_linear`, block RMS-norm, and gated nonlinearity.
 *
 * Coverage:
 *   - `linear`: matched multisets route weights as Σ mult_out · mult_in.
 *   - `linear`: irrep labels with no corresponding out channel contribute
 *     nothing.
 *   - `linear`: SO(3) equivariance on a matched multiset (apply `linear`
 *     before and after a random rotation; compare).
 *   - `linear` backward pass matches finite difference.
 *   - RMS norm: output has unit block-RMS when scale = 1.
 *   - RMS norm backward matches finite difference.
 *   - Gated nonlinearity: multiplies each irrep block by its gate scalar,
 *     leaving non-gated blocks untouched.
 */
#include "harness.h"
#include <irrep/equivariant_layers.h>
#include <irrep/multiset.h>
#include <irrep/so3.h>
#include <irrep/wigner_d.h>
#include <irrep/spherical_harmonics.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    IRREP_TEST_START("equivariant_layers");

    const double tol = 1e-12;

    /* -------- linear: matched multisets count weights as Σ mult_out · mult_in ------- */
    {
        irrep_multiset_t *in = irrep_multiset_parse("2x0e + 3x1o");
        irrep_multiset_t *out = irrep_multiset_parse("4x0e + 5x1o");
        irrep_linear_t   *lin = irrep_linear_build(in, out, 1, 1);
        IRREP_ASSERT(lin != NULL);
        IRREP_ASSERT(irrep_linear_num_weights(lin) == 2 * 4 + 3 * 5); /* 8 + 15 = 23 */
        irrep_linear_free(lin);
        irrep_multiset_free(in);
        irrep_multiset_free(out);
    }

    /* -------- linear: unmatched labels contribute nothing ---------- */
    {
        irrep_multiset_t *in = irrep_multiset_parse("1x0e + 1x1o");
        irrep_multiset_t *out = irrep_multiset_parse("1x0e + 1x2e"); /* 1o ≠ 2e */
        irrep_linear_t   *lin = irrep_linear_build(in, out, 1, 1);
        IRREP_ASSERT(irrep_linear_num_weights(lin) == 1); /* only 0e matches */

        double W[1] = {2.0};
        double x[4] = {7.0, 1.0, 2.0, 3.0}; /* 0e + 1o */
        double y[6];                        /* 0e + 2e */
        irrep_linear_apply(lin, W, x, y);
        IRREP_ASSERT_NEAR(y[0], 14.0, tol); /* scalar mixed */
        /* l=2 block untouched (zero) */
        for (int i = 1; i < 6; ++i)
            IRREP_ASSERT_NEAR(y[i], 0.0, tol);

        irrep_linear_free(lin);
        irrep_multiset_free(in);
        irrep_multiset_free(out);
    }

    /* -------- linear: SO(3) equivariance on a matched multiset ------- */
    {
        irrep_multiset_t *M = irrep_multiset_parse("1x0e + 2x1o");
        irrep_linear_t   *lin = irrep_linear_build(M, M, 1, 1);
        int               nw = irrep_linear_num_weights(lin);
        IRREP_ASSERT(nw == 1 + 2 * 2);

        double weights[1 + 4] = {1.3, 0.7, -0.5, 0.2, 0.9};
        double x[7] = {1.0, 0.3, 0.1, -0.2, 0.5, -0.7, 0.4};
        double y[7];
        irrep_linear_apply(lin, weights, x, y);

        /* Rotate x with real-basis D, apply linear, should equal D · y. */
        double _Complex U1[9];
        irrep_sph_harm_complex_to_real(1, U1);
        double _Complex D1[9];
        irrep_wigner_D_matrix(1, D1, 0.3, 0.6, 1.1);
        double _Complex tmp[9], Dr1[9];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double _Complex s = 0;
                for (int k = 0; k < 3; ++k)
                    s += U1[i * 3 + k] * D1[k * 3 + j];
                tmp[i * 3 + j] = s;
            }
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double _Complex s = 0;
                for (int k = 0; k < 3; ++k)
                    s += tmp[i * 3 + k] * conj(U1[j * 3 + k]);
                Dr1[i * 3 + j] = s;
            }

        double x_rot[7];
        x_rot[0] = x[0]; /* scalar block unchanged */
        for (int u = 0; u < 2; ++u) {
            double out3[3] = {0, 0, 0};
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    out3[i] += creal(Dr1[i * 3 + j]) * x[1 + u * 3 + j];
            for (int i = 0; i < 3; ++i)
                x_rot[1 + u * 3 + i] = out3[i];
        }

        double y_rot_in[7];
        irrep_linear_apply(lin, weights, x_rot, y_rot_in);

        double y_out_rot[7];
        y_out_rot[0] = y[0];
        for (int u = 0; u < 2; ++u) {
            double out3[3] = {0, 0, 0};
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    out3[i] += creal(Dr1[i * 3 + j]) * y[1 + u * 3 + j];
            for (int i = 0; i < 3; ++i)
                y_out_rot[1 + u * 3 + i] = out3[i];
        }

        for (int i = 0; i < 7; ++i)
            IRREP_ASSERT(fabs(y_rot_in[i] - y_out_rot[i]) < 1e-10);

        irrep_linear_free(lin);
        irrep_multiset_free(M);
    }

    /* -------- linear backward vs finite difference ---------- */
    {
        irrep_multiset_t *in = irrep_multiset_parse("2x0e + 1x1o");
        irrep_multiset_t *out = irrep_multiset_parse("1x0e + 1x1o");
        irrep_linear_t   *lin = irrep_linear_build(in, out, 1, 1);
        int               nw = irrep_linear_num_weights(lin);
        double            W[10];
        for (int i = 0; i < nw; ++i)
            W[i] = 0.3 * (i + 1);
        double x[5] = {1.2, -0.4, 0.7, 0.9, -0.3};
        double grad_out[4] = {1.0, 0.5, -0.2, 0.8};
        double grad_W[10] = {0}, grad_in[5] = {0};
        irrep_linear_backward(lin, W, x, grad_out, grad_W, grad_in);

        double h = 1e-7;
        for (int i = 0; i < nw; ++i) {
            double save = W[i];
            double yp[4], yn[4];
            W[i] = save + h;
            irrep_linear_apply(lin, W, x, yp);
            W[i] = save - h;
            irrep_linear_apply(lin, W, x, yn);
            W[i] = save;
            double fd = 0;
            for (int k = 0; k < 4; ++k)
                fd += grad_out[k] * (yp[k] - yn[k]) / (2 * h);
            IRREP_ASSERT(fabs(grad_W[i] - fd) < 1e-6);
        }
        for (int i = 0; i < 5; ++i) {
            double save = x[i];
            double yp[4], yn[4];
            x[i] = save + h;
            irrep_linear_apply(lin, W, x, yp);
            x[i] = save - h;
            irrep_linear_apply(lin, W, x, yn);
            x[i] = save;
            double fd = 0;
            for (int k = 0; k < 4; ++k)
                fd += grad_out[k] * (yp[k] - yn[k]) / (2 * h);
            IRREP_ASSERT(fabs(grad_in[i] - fd) < 1e-6);
        }

        irrep_linear_free(lin);
        irrep_multiset_free(in);
        irrep_multiset_free(out);
    }

    /* -------- RMS norm: output has unit block-RMS when scale=1 ---------- */
    {
        irrep_multiset_t *M = irrep_multiset_parse("1x1o");
        double            in[3] = {3.0, 4.0, 0.0}; /* |x| = 5, d=3 → rms = √(25/3) */
        double            scales[1] = {1.0};
        double            out[3];
        irrep_norm_rms(M, 1, scales, in, out);
        /* out has rms = scale = 1; and direction preserved */
        double out_rms = sqrt((out[0] * out[0] + out[1] * out[1] + out[2] * out[2]) / 3.0);
        IRREP_ASSERT_NEAR(out_rms, 1.0, 1e-12);
        /* direction: out / |out| == in / |in| */
        double in_norm = sqrt(in[0] * in[0] + in[1] * in[1] + in[2] * in[2]);
        for (int i = 0; i < 3; ++i) {
            IRREP_ASSERT_NEAR(out[i] / sqrt(3.0), in[i] / in_norm, 1e-12);
        }
        irrep_multiset_free(M);
    }

    /* -------- RMS norm backward vs FD ---------- */
    {
        irrep_multiset_t *M = irrep_multiset_parse("1x1o");
        double            in[3] = {0.6, -0.3, 0.9};
        double            scales[1] = {1.4};
        double            go[3] = {0.7, 0.2, -0.5};
        double            grad_s[1] = {0}, grad_in[3] = {0};
        irrep_norm_rms_backward(M, 1, scales, in, go, grad_s, grad_in);

        double h = 1e-7;
        for (int i = 0; i < 3; ++i) {
            double save = in[i];
            double yp[3], yn[3];
            in[i] = save + h;
            irrep_norm_rms(M, 1, scales, in, yp);
            in[i] = save - h;
            irrep_norm_rms(M, 1, scales, in, yn);
            in[i] = save;
            double fd = 0;
            for (int k = 0; k < 3; ++k)
                fd += go[k] * (yp[k] - yn[k]) / (2 * h);
            IRREP_ASSERT(fabs(grad_in[i] - fd) < 1e-6);
        }
        double save = scales[0];
        double yp[3], yn[3];
        scales[0] = save + h;
        irrep_norm_rms(M, 1, scales, in, yp);
        scales[0] = save - h;
        irrep_norm_rms(M, 1, scales, in, yn);
        scales[0] = save;
        double fd_s = 0;
        for (int k = 0; k < 3; ++k)
            fd_s += go[k] * (yp[k] - yn[k]) / (2 * h);
        IRREP_ASSERT(fabs(grad_s[0] - fd_s) < 1e-6);
        irrep_multiset_free(M);
    }

    /* -------- gate: multiplies each block by its gate ---------- */
    {
        irrep_multiset_t *M = irrep_multiset_parse("2x1o"); /* 2 copies, d=3, 2 gates */
        double            in[6] = {1, 2, 3, 4, 5, 6};
        double            gates[2] = {-1.0, 2.5};
        double            out[6];
        irrep_gate_apply(M, 1, gates, in, out);
        for (int i = 0; i < 3; ++i)
            IRREP_ASSERT_NEAR(out[i], -in[i], tol);
        for (int i = 0; i < 3; ++i)
            IRREP_ASSERT_NEAR(out[3 + i], 2.5 * in[3 + i], tol);
        irrep_multiset_free(M);
    }

    return IRREP_TEST_END();
}
