/* SPDX-License-Identifier: MIT */
/* A two-layer E(3)-equivariant MLP on irrep feature vectors.
 *
 *   in  → Linear₁ → RMS-norm → Gate → Linear₂ → out
 *
 * Uses the M9 primitives. The gate is SiLU-like (sigmoid of a l=0 scalar).
 *
 * The point is that the whole thing is equivariant by construction —
 * rotating the input by R must rotate the output by the same R. We verify
 * that numerically at the end. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/multiset.h>
#include <irrep/equivariant_layers.h>
#include <irrep/so3.h>
#include <irrep/wigner_d.h>
#include <irrep/spherical_harmonics.h>
#include <complex.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

static double sigmoid(double x) { return 1.0 / (1.0 + exp(-x)); }

static void build_real_rot_1(double R[9],
                             double alpha, double beta, double gamma) {
    double _Complex U[9];
    double _Complex D[9];
    irrep_sph_harm_complex_to_real(1, U);
    irrep_wigner_D_matrix(1, D, alpha, beta, gamma);
    double _Complex tmp[9];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double _Complex s = 0;
            for (int k = 0; k < 3; ++k) s += U[i*3+k] * D[k*3+j];
            tmp[i*3+j] = s;
        }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double _Complex s = 0;
            for (int k = 0; k < 3; ++k) s += tmp[i*3+k] * conj(U[j*3+k]);
            R[i*3+j] = creal(s);
        }
}

static void apply_layer(const irrep_linear_t *lin, const double *W,
                        const irrep_multiset_t *middle,
                        const double *norm_scales,
                        const double *gates,
                        const double *in_vec, double *middle_vec) {
    /* lin: in → middle; then RMS-norm on middle; then gate on middle. */
    double tmp[32], tmp2[32];
    irrep_linear_apply(lin, W, in_vec, tmp);
    irrep_norm_rms(middle, 1, norm_scales, tmp, tmp2);
    irrep_gate_apply(middle, 1, gates, tmp2, middle_vec);
}

int main(void) {
    irrep_multiset_t *in_m     = irrep_multiset_parse("1x0e + 1x1o");
    irrep_multiset_t *middle_m = irrep_multiset_parse("2x0e + 2x1o");
    irrep_multiset_t *out_m    = irrep_multiset_parse("1x0e + 1x1o");

    irrep_linear_t *L1 = irrep_linear_build(in_m, middle_m, 1, 1);
    irrep_linear_t *L2 = irrep_linear_build(middle_m, out_m, 1, 1);

    /* Seeded deterministic weights. */
    int nw1 = irrep_linear_num_weights(L1);
    int nw2 = irrep_linear_num_weights(L2);
    double W1[32], W2[32];
    for (int i = 0; i < nw1; ++i) W1[i] = 0.31 * (i + 1) - 0.8;
    for (int i = 0; i < nw2; ++i) W2[i] = 0.27 * (i + 1) - 0.5;

    /* One scale per (term, copy) in the middle = 2 + 2 = 4 scales.
     * Gate layout same — 4 scalars; computed from a pre-tanh logit. */
    double norm_scales[4] = { 1.0, 1.1, 0.9, 1.2 };
    double gate_logits[4] = { 0.5, -0.3, 1.2, -0.1 };
    double gates[4];
    for (int i = 0; i < 4; ++i) gates[i] = sigmoid(gate_logits[i]);

    /* Input feature: scalar + vector. */
    double x[4]     = { 0.7,   0.3, -0.5, 0.2 };
    double middle[8];
    double y[4];
    apply_layer(L1, W1, middle_m, norm_scales, gates, x, middle);
    irrep_linear_apply(L2, W2, middle, y);

    /* Equivariance check: rotate input, re-apply, compare to rotated output. */
    double R[9];
    build_real_rot_1(R, 0.4, 0.9, 1.3);

    double x_rot[4];
    x_rot[0] = x[0];
    for (int i = 0; i < 3; ++i) {
        double s = 0;
        for (int j = 0; j < 3; ++j) s += R[i*3+j] * x[1+j];
        x_rot[1+i] = s;
    }
    double middle_rot[8], y_rot[4];
    apply_layer(L1, W1, middle_m, norm_scales, gates, x_rot, middle_rot);
    irrep_linear_apply(L2, W2, middle_rot, y_rot);

    double y_direct_rot[4];
    y_direct_rot[0] = y[0];
    for (int i = 0; i < 3; ++i) {
        double s = 0;
        for (int j = 0; j < 3; ++j) s += R[i*3+j] * y[1+j];
        y_direct_rot[1+i] = s;
    }

    printf("equivariant_mlp\n");
    printf("  input  x = (% .4f   % .4f % .4f % .4f)\n", x[0], x[1], x[2], x[3]);
    printf("  output y = (% .4f   % .4f % .4f % .4f)\n", y[0], y[1], y[2], y[3]);
    printf("  equivariance check: f(R·x) vs R·f(x)\n");
    double max_err = 0.0;
    for (int i = 0; i < 4; ++i) {
        double e = fabs(y_rot[i] - y_direct_rot[i]);
        if (e > max_err) max_err = e;
    }
    printf("  max |f(R·x) − R·f(x)| = %.3e\n", max_err);
    int ok = max_err < 1e-10;
    printf("  %s\n", ok ? "OK (equivariant)" : "FAIL");

    irrep_linear_free(L1);
    irrep_linear_free(L2);
    irrep_multiset_free(in_m);
    irrep_multiset_free(middle_m);
    irrep_multiset_free(out_m);
    return ok ? 0 : 1;
}
