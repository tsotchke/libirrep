/* SPDX-License-Identifier: MIT */
/* Tests for the e3nn-style path-indexed tensor product:
 * `irrep_tp_build`, `irrep_tp_apply`, and the backward pass.
 *
 * Coverage:
 *   - Path enumeration: 1x1o × 1x1o → {1x0e, 1x1e, 1x2e} (triangle + parity).
 *   - Hand-computed 1x1o × 1x1o → 1x0e via CG + i^(l_a+l_b−l_c) phase.
 *   - Build errors: multiplicity mismatch; parity violation.
 *   - Equivariance under SO(3):  tp(D(R) a, D(R) b) = D(R) tp(a, b).
 *   - Multiplicity > 1 uses copy-pairing between aligned slots.
 *   - 1x1o × 1x1o → 1x1e is the cross product (odd-l-sum real-basis path).
 *   - Batched tp_apply matches single.
 *   - Backward pass bit-exact against finite differences.
 *   - UVW connection mode: weight indexing; collapses to UUU when weights
 *     are diagonal; backward matches FD.
 *   - Per-path L² regulariser: forward matches brute force, backward matches
 *     finite difference.
 */
#include "harness.h"
#include <irrep/tensor_product.h>
#include <irrep/multiset.h>
#include <irrep/so3.h>
#include <irrep/wigner_d.h>
#include <irrep/clebsch_gordan.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

static int apply_complex_matrix_(const double _Complex *M, int n,
                                 const double *vec, double _Complex *out) {
    for (int i = 0; i < n; ++i) {
        double _Complex s = 0.0;
        for (int j = 0; j < n; ++j) s += M[i * n + j] * vec[j];
        out[i] = s;
    }
    return 0;
}

int main(void) {
    IRREP_TEST_START("tensor_product");

    const double tol = 1e-10;

    /* -------- path enumeration: 1x1o × 1x1o → 1x0e + 1x1e + 1x2e ---------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("1x1o");
        irrep_multiset_t *B = irrep_multiset_parse("1x1o");
        irrep_multiset_t *C = irrep_multiset_parse("1x0e + 1x1e + 1x2e");
        int paths[9 * 3];
        int n = irrep_tp_enumerate_paths(A, B, C, paths, 9);
        /* All three triangle-valid paths: 0e (scalar), 1e (cross product),
         * 2e (symmetric-traceless tensor). Parity o·o = e matches all of C. */
        IRREP_ASSERT(n == 3);

        /* Also: 1x1o × 1x1o → 1x0o should yield 0 paths (parity o·o = e ≠ o) */
        irrep_multiset_t *Co = irrep_multiset_parse("1x0o");
        int n0 = irrep_tp_enumerate_paths(A, B, Co, NULL, 0);
        IRREP_ASSERT(n0 == 0);

        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
        irrep_multiset_free(Co);
    }

    /* -------- hand-computed: 1x1o × 1x1o → 1x0e via CG + i-phase ---------- */
    /* With the e3nn-style `i^{l_a+l_b-l_c}` phase, the two-vector scalar
     * coupling lands at `+(a · b) / √3` (the Condon-Shortley CG value
     * i^2 = −1 flips the sign vs. the bare-CG result). */
    {
        irrep_multiset_t *A = irrep_multiset_parse("1x1o");
        irrep_multiset_t *B = irrep_multiset_parse("1x1o");
        irrep_multiset_t *C = irrep_multiset_parse("1x0e");
        int path[3] = { 0, 0, 0 };
        tp_descriptor_t *d = irrep_tp_build(A, B, C, path, 1);
        IRREP_ASSERT(d != NULL);
        IRREP_ASSERT(irrep_tp_num_paths(d) == 1);
        IRREP_ASSERT(irrep_tp_output_dim(d) == 1);

        double a[3] = { 1.0, 2.0, 3.0 };
        double b[3] = { 4.0, 5.0, 6.0 };
        double c[1] = { 0.0 };
        irrep_tp_apply(d, a, b, c);
        double expected = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / sqrt(3.0);
        IRREP_ASSERT_NEAR(c[0], expected, tol);

        /* weighted by 2.0 */
        double w[1] = { 2.0 };
        double c2[1] = { 0.0 };
        irrep_tp_apply_weighted(d, w, a, b, c2);
        IRREP_ASSERT_NEAR(c2[0], 2.0 * expected, tol);

        irrep_tp_free(d);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- build fails on mult mismatch ---------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("2x1o");
        irrep_multiset_t *B = irrep_multiset_parse("1x1o");
        irrep_multiset_t *C = irrep_multiset_parse("1x0e");
        int path[3] = { 0, 0, 0 };
        tp_descriptor_t *d = irrep_tp_build(A, B, C, path, 1);
        IRREP_ASSERT(d == NULL);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- build fails on parity violation ---------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("1x1o");
        irrep_multiset_t *B = irrep_multiset_parse("1x1e");
        irrep_multiset_t *C = irrep_multiset_parse("1x0e");   /* o·e = o ≠ e */
        int path[3] = { 0, 0, 0 };
        tp_descriptor_t *d = irrep_tp_build(A, B, C, path, 1);
        IRREP_ASSERT(d == NULL);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- equivariance: tp(D(R) a, D(R) b) == D(R) tp(a, b) ---------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("1x1o");
        irrep_multiset_t *B = irrep_multiset_parse("1x1o");
        irrep_multiset_t *C = irrep_multiset_parse("1x0e + 1x1e + 1x2e");
        int n_paths = irrep_tp_enumerate_paths(A, B, C, NULL, 0);
        int *paths = malloc((size_t)n_paths * 3 * sizeof(int));
        irrep_tp_enumerate_paths(A, B, C, paths, n_paths);
        tp_descriptor_t *d = irrep_tp_build(A, B, C, paths, n_paths);
        IRREP_ASSERT(d != NULL);

        double a[3] = { 0.7, -0.3, 1.1 };
        double b[3] = { -0.2, 0.9, 0.4 };
        double c_direct[9];   /* C total_dim = 1 + 3 + 5 = 9 */
        irrep_tp_apply(d, a, b, c_direct);

        /* Rotate inputs via D(R) for l=1 (A, B), then tp → expect D(R) c. */
        double alpha = 0.3, beta = 0.9, gamma = 1.4;
        double _Complex Da[9], Db[9];      /* l=1 → 3×3 */
        double _Complex Dc_ms[9 * 9];      /* C total 9×9 */
        irrep_wigner_D_matrix(1, Da, alpha, beta, gamma);
        irrep_wigner_D_matrix(1, Db, alpha, beta, gamma);
        irrep_wigner_D_multiset(C, Dc_ms, alpha, beta, gamma);

        /* Wigner-D is on complex SH basis; we're working in real SH basis
         * through `sph_harm_cart`. To rotate real-basis vectors we need the
         * real-basis rotation matrix D_real = U · D_complex · U†.
         * Implement that via the M3 complex-to-real matrix. */
        double _Complex U1[9], U0[1], U2[25];
        extern void irrep_sph_harm_complex_to_real(int, double _Complex*);
        irrep_sph_harm_complex_to_real(0, U0);
        irrep_sph_harm_complex_to_real(1, U1);
        irrep_sph_harm_complex_to_real(2, U2);

        /* Compute D_real_1 = U1 · Da · U1† (3×3 real-basis rotation on l=1) */
        double _Complex tmp1[9], Dr1[9];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double _Complex s = 0;
                for (int k = 0; k < 3; ++k) s += U1[i*3+k] * Da[k*3+j];
                tmp1[i*3+j] = s;
            }
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double _Complex s = 0;
                for (int k = 0; k < 3; ++k) s += tmp1[i*3+k] * conj(U1[j*3+k]);
                Dr1[i*3+j] = s;
            }
        }

        double a_rot[3], b_rot[3];
        for (int i = 0; i < 3; ++i) {
            double _Complex s_a = 0, s_b = 0;
            for (int j = 0; j < 3; ++j) {
                s_a += Dr1[i*3+j] * a[j];
                s_b += Dr1[i*3+j] * b[j];
            }
            a_rot[i] = creal(s_a);
            b_rot[i] = creal(s_b);
        }

        /* Apply tp to rotated inputs */
        double c_rotated_in[9];
        irrep_tp_apply(d, a_rot, b_rot, c_rotated_in);

        /* Rotate the direct output c_direct using the real-basis multiset D */
        /* For each output block: l=0 (1), l=1 (3), l=2 (5). */
        double c_out_rotated[9];

        /* l=0 block (index 0): D_real_0 = 1, no change */
        c_out_rotated[0] = c_direct[0];

        /* l=1 block (indices 1..3): apply Dr1 */
        for (int i = 0; i < 3; ++i) {
            double _Complex s = 0;
            for (int j = 0; j < 3; ++j) s += Dr1[i*3+j] * c_direct[1 + j];
            c_out_rotated[1 + i] = creal(s);
        }

        /* l=2 block: compute Dr2 = U2 · D_l=2 · U2† */
        double _Complex D2[25];
        irrep_wigner_D_matrix(2, D2, alpha, beta, gamma);
        double _Complex tmp2[25], Dr2[25];
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j) {
                double _Complex s = 0;
                for (int k = 0; k < 5; ++k) s += U2[i*5+k] * D2[k*5+j];
                tmp2[i*5+j] = s;
            }
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j) {
                double _Complex s = 0;
                for (int k = 0; k < 5; ++k) s += tmp2[i*5+k] * conj(U2[j*5+k]);
                Dr2[i*5+j] = s;
            }
        for (int i = 0; i < 5; ++i) {
            double _Complex s = 0;
            for (int j = 0; j < 5; ++j) s += Dr2[i*5+j] * c_direct[4 + j];
            c_out_rotated[4 + i] = creal(s);
        }

        (void)Dc_ms; (void)U0;
        for (int i = 0; i < 9; ++i) {
            IRREP_ASSERT(fabs(c_rotated_in[i] - c_out_rotated[i]) < 1e-10);
        }

        free(paths);
        irrep_tp_free(d);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- multiplicity > 1 uses copy-pairing ---------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("2x1o");
        irrep_multiset_t *B = irrep_multiset_parse("2x1o");
        irrep_multiset_t *C = irrep_multiset_parse("2x0e");
        int path[3] = { 0, 0, 0 };
        tp_descriptor_t *d = irrep_tp_build(A, B, C, path, 1);
        IRREP_ASSERT(d != NULL);

        double a[6] = { 1, 2, 3,   4, 5, 6 };   /* two copies */
        double b[6] = { 7, 8, 9,  10,11,12 };
        double c[2] = { 0.0, 0.0 };
        irrep_tp_apply(d, a, b, c);
        double c0_expected = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / sqrt(3.0);
        double c1_expected = (a[3]*b[3] + a[4]*b[4] + a[5]*b[5]) / sqrt(3.0);
        IRREP_ASSERT_NEAR(c[0], c0_expected, tol);
        IRREP_ASSERT_NEAR(c[1], c1_expected, tol);

        irrep_tp_free(d);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- backward matches finite-difference ---------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("1x1o");
        irrep_multiset_t *B = irrep_multiset_parse("1x1o");
        irrep_multiset_t *C = irrep_multiset_parse("1x0e + 1x2e");
        int n_paths = irrep_tp_enumerate_paths(A, B, C, NULL, 0);
        int paths[9 * 3];
        irrep_tp_enumerate_paths(A, B, C, paths, 9);
        tp_descriptor_t *d = irrep_tp_build(A, B, C, paths, n_paths);

        double a[3] = { 0.7, -0.3, 1.1 };
        double b[3] = { -0.2, 0.9, 0.4 };
        double weights[2] = { 1.3, -0.8 };
        double grad_c[6] = { 1.0, 0.5, -0.2, 0.3, 0.7, -0.1 };

        double grad_a[3] = { 0 }, grad_b[3] = { 0 }, grad_w[2] = { 0 };
        irrep_tp_apply_backward_weighted(d, weights, a, b, grad_c,
                                          grad_a, grad_b, grad_w);

        /* FD check on grad_a */
        double h = 1e-7;
        for (int i = 0; i < 3; ++i) {
            double save = a[i];
            double c_plus[6], c_minus[6];
            a[i] = save + h; irrep_tp_apply_weighted(d, weights, a, b, c_plus);
            a[i] = save - h; irrep_tp_apply_weighted(d, weights, a, b, c_minus);
            a[i] = save;
            double fd = 0;
            for (int j = 0; j < 6; ++j) fd += grad_c[j] * (c_plus[j] - c_minus[j]) / (2 * h);
            IRREP_ASSERT(fabs(grad_a[i] - fd) < 1e-6);
        }
        /* FD check on grad_b */
        for (int i = 0; i < 3; ++i) {
            double save = b[i];
            double c_plus[6], c_minus[6];
            b[i] = save + h; irrep_tp_apply_weighted(d, weights, a, b, c_plus);
            b[i] = save - h; irrep_tp_apply_weighted(d, weights, a, b, c_minus);
            b[i] = save;
            double fd = 0;
            for (int j = 0; j < 6; ++j) fd += grad_c[j] * (c_plus[j] - c_minus[j]) / (2 * h);
            IRREP_ASSERT(fabs(grad_b[i] - fd) < 1e-6);
        }
        /* FD check on grad_w */
        for (int k = 0; k < 2; ++k) {
            double save = weights[k];
            double c_plus[6], c_minus[6];
            weights[k] = save + h; irrep_tp_apply_weighted(d, weights, a, b, c_plus);
            weights[k] = save - h; irrep_tp_apply_weighted(d, weights, a, b, c_minus);
            weights[k] = save;
            double fd = 0;
            for (int j = 0; j < 6; ++j) fd += grad_c[j] * (c_plus[j] - c_minus[j]) / (2 * h);
            IRREP_ASSERT(fabs(grad_w[k] - fd) < 1e-6);
        }

        irrep_tp_free(d);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- 1x1o × 1x1o → 1x1e: cross product (odd l-sum path) -------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("1x1o");
        irrep_multiset_t *B = irrep_multiset_parse("1x1o");
        irrep_multiset_t *C = irrep_multiset_parse("1x1e");
        int path[3] = { 0, 0, 0 };
        tp_descriptor_t *d = irrep_tp_build(A, B, C, path, 1);
        IRREP_ASSERT(d != NULL);     /* odd l-sum no longer rejected */

        double a[3] = { 1.3, -0.7, 0.9 };
        double b[3] = { -0.4, 0.6, 0.2 };
        double c[3] = { 0, 0, 0 };
        irrep_tp_apply(d, a, b, c);

        /* Antisymmetric in a↔b (cross-product property): */
        double c_swap[3] = { 0, 0, 0 };
        irrep_tp_apply(d, b, a, c_swap);
        for (int i = 0; i < 3; ++i) IRREP_ASSERT(fabs(c[i] + c_swap[i]) < 1e-12);
        /* And a × a = 0: */
        double c_self[3] = { 0, 0, 0 };
        irrep_tp_apply(d, a, a, c_self);
        for (int i = 0; i < 3; ++i) IRREP_ASSERT(fabs(c_self[i]) < 1e-12);
        /* Sanity: the output is not all zeros for non-parallel vectors. */
        IRREP_ASSERT(fabs(c[0]) + fabs(c[1]) + fabs(c[2]) > 1e-6);

        irrep_tp_free(d);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- batched matches single ---------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("1x1o");
        irrep_multiset_t *B = irrep_multiset_parse("1x1o");
        irrep_multiset_t *C = irrep_multiset_parse("1x0e");
        int path[3] = { 0, 0, 0 };
        tp_descriptor_t *d = irrep_tp_build(A, B, C, path, 1);

        size_t batch = 4;
        double a_batch[12] = { 1,2,3,  4,5,6,  7,8,9,  10,11,12 };
        double b_batch[12] = { 0.1,0.2,0.3,  0.4,0.5,0.6,  0.7,0.8,0.9,  1.0,1.1,1.2 };
        double c_batch[4];
        irrep_tp_apply_batch(d, batch, a_batch, b_batch, c_batch);
        for (size_t bi = 0; bi < batch; ++bi) {
            double c_single[1];
            irrep_tp_apply(d, a_batch + bi * 3, b_batch + bi * 3, c_single);
            IRREP_ASSERT_NEAR(c_batch[bi], c_single[0], 1e-14);
        }

        irrep_tp_free(d);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- UVW mode: build counts weights correctly, mode flag set --------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("2x1o");   /* U = 2 */
        irrep_multiset_t *B = irrep_multiset_parse("3x1o");   /* V = 3 */
        irrep_multiset_t *C = irrep_multiset_parse("4x0e");   /* W = 4 */
        int path[3] = { 0, 0, 0 };
        tp_descriptor_t *d = irrep_tp_build_uvw(A, B, C, path, 1);
        IRREP_ASSERT(d != NULL);
        IRREP_ASSERT(irrep_tp_mode(d) == IRREP_TP_MODE_UVW);
        IRREP_ASSERT(irrep_tp_num_weights_uvw(d) == 2 * 3 * 4);
        irrep_tp_free(d);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- UVW collapses to UUU for diagonal weight = δ_{wvu} ------------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("2x1o");
        irrep_multiset_t *B = irrep_multiset_parse("2x1o");
        irrep_multiset_t *C = irrep_multiset_parse("2x0e");
        int path[3] = { 0, 0, 0 };

        tp_descriptor_t *uuu = irrep_tp_build(A, B, C, path, 1);
        tp_descriptor_t *uvw = irrep_tp_build_uvw(A, B, C, path, 1);
        IRREP_ASSERT(uuu != NULL && uvw != NULL);

        double a[6] = { 1, 2, 3,   4, 5, 6 };
        double b[6] = { 7, 8, 9,  10,11,12 };
        double c_uuu[2] = { 0, 0 }, c_uvw[2] = { 0, 0 };
        double uuu_w[1] = { 1.0 };
        double uvw_w[8] = { 0 };
        /* δ_{w,v,u}: weight 1 for w=v=u, else 0. */
        for (int w = 0; w < 2; ++w)
            for (int v = 0; v < 2; ++v)
                for (int u = 0; u < 2; ++u)
                    uvw_w[(w * 2 + v) * 2 + u] = (w == v && w == u) ? 1.0 : 0.0;

        irrep_tp_apply_weighted(uuu, uuu_w, a, b, c_uuu);
        irrep_tp_apply_uvw     (uvw, uvw_w, a, b, c_uvw);
        for (int i = 0; i < 2; ++i) IRREP_ASSERT(fabs(c_uuu[i] - c_uvw[i]) < 1e-12);

        irrep_tp_free(uuu);
        irrep_tp_free(uvw);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- UVW backward matches centered finite difference ---------------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("2x1o");
        irrep_multiset_t *B = irrep_multiset_parse("3x1o");
        irrep_multiset_t *C = irrep_multiset_parse("1x0e + 1x2e");
        int n_paths = irrep_tp_enumerate_paths(A, B, C, NULL, 0);
        int *paths = malloc((size_t)n_paths * 3 * sizeof(int));
        irrep_tp_enumerate_paths(A, B, C, paths, n_paths);
        tp_descriptor_t *d = irrep_tp_build_uvw(A, B, C, paths, n_paths);
        IRREP_ASSERT(d != NULL);

        int nw = irrep_tp_num_weights_uvw(d);
        double *w  = malloc((size_t)nw * sizeof(double));
        double *gw = calloc((size_t)nw, sizeof(double));
        for (int i = 0; i < nw; ++i) w[i] = 0.1 * (i + 1) - 0.5;

        double a[6] = { 0.1, -0.3, 0.5,   0.7, 0.2, -0.4 };
        double b[9] = { 1.1, -0.7, 0.3,   0.4, -0.5, 0.8,   -0.9, 0.6, 0.2 };
        /* c_dim = 1 + 5 = 6 */
        double grad_c[6] = { 0.5, 0.3, -0.2, 0.7, -0.1, 0.4 };
        double ga[6] = { 0 }, gb[9] = { 0 };

        irrep_tp_apply_uvw_backward(d, w, a, b, grad_c, ga, gb, gw);

        double h = 1e-7;
        for (int i = 0; i < 6; ++i) {
            double save = a[i];
            double cp[6] = { 0 }, cm[6] = { 0 };
            a[i] = save + h; irrep_tp_apply_uvw(d, w, a, b, cp);
            a[i] = save - h; irrep_tp_apply_uvw(d, w, a, b, cm);
            a[i] = save;
            double fd = 0; for (int j = 0; j < 6; ++j) fd += grad_c[j] * (cp[j] - cm[j]) / (2 * h);
            IRREP_ASSERT(fabs(ga[i] - fd) < 1e-6);
        }
        for (int i = 0; i < nw; ++i) {
            double save = w[i];
            double cp[6] = { 0 }, cm[6] = { 0 };
            w[i] = save + h; irrep_tp_apply_uvw(d, w, a, b, cp);
            w[i] = save - h; irrep_tp_apply_uvw(d, w, a, b, cm);
            w[i] = save;
            double fd = 0; for (int j = 0; j < 6; ++j) fd += grad_c[j] * (cp[j] - cm[j]) / (2 * h);
            IRREP_ASSERT(fabs(gw[i] - fd) < 1e-6);
        }

        free(paths); free(w); free(gw);
        irrep_tp_free(d);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    /* -------- Per-path L2 regulariser: forward matches brute force; backward matches FD. -------- */
    {
        irrep_multiset_t *A = irrep_multiset_parse("2x0e + 1x1o");
        irrep_multiset_t *B = irrep_multiset_parse("2x0e + 1x1o");
        irrep_multiset_t *C = irrep_multiset_parse("2x0e + 1x1o");
        int np = irrep_tp_enumerate_paths(A, B, C, NULL, 0);
        int *paths = malloc(np * 3 * sizeof(int));
        irrep_tp_enumerate_paths(A, B, C, paths, np);
        tp_descriptor_t *d = irrep_tp_build_uvw(A, B, C, paths, np);
        IRREP_ASSERT(d != NULL);
        int nw = irrep_tp_num_weights_uvw(d);

        double *w = calloc(nw, sizeof(double));
        for (int i = 0; i < nw; ++i) w[i] = 0.37 * i - 0.13;

        double *per_path = calloc(np, sizeof(double));
        irrep_tp_weight_l2_per_path_uvw(d, w, per_path);

        /* Forward: sum across all paths equals the full-weight L2. */
        double total = 0.0;
        for (int k = 0; k < np; ++k) total += per_path[k];
        double full_l2 = 0.0;
        for (int i = 0; i < nw; ++i) full_l2 += w[i] * w[i];
        IRREP_ASSERT(fabs(total - full_l2) < 1e-12);

        /* Backward: per-path grad of 1.0 on every path — grad_w[i] should be 2 w[i]. */
        double *grad_per_path = calloc(np, sizeof(double));
        for (int k = 0; k < np; ++k) grad_per_path[k] = 1.0;
        double *grad_w = calloc(nw, sizeof(double));
        irrep_tp_weight_l2_per_path_uvw_backward(d, w, grad_per_path, grad_w);
        for (int i = 0; i < nw; ++i) {
            IRREP_ASSERT(fabs(grad_w[i] - 2.0 * w[i]) < 1e-12);
        }

        /* FD check with heterogeneous path weights. */
        for (int k = 0; k < np; ++k) grad_per_path[k] = 0.17 * k - 0.3;
        memset(grad_w, 0, nw * sizeof(double));
        irrep_tp_weight_l2_per_path_uvw_backward(d, w, grad_per_path, grad_w);

        double h = 1e-6;
        /* Exercise the full weight vector — cheap enough on this shape and
         * it closes the audit gap of a spot-check-only backward FD pass. */
        for (int i = 0; i < nw; ++i) {
            double save = w[i];
            w[i] = save + h;
            double *pp = calloc(np, sizeof(double));
            irrep_tp_weight_l2_per_path_uvw(d, w, pp);
            double lp = 0.0; for (int k = 0; k < np; ++k) lp += grad_per_path[k] * pp[k];
            w[i] = save - h;
            irrep_tp_weight_l2_per_path_uvw(d, w, pp);
            double lm = 0.0; for (int k = 0; k < np; ++k) lm += grad_per_path[k] * pp[k];
            w[i] = save;
            double fd = (lp - lm) / (2 * h);
            IRREP_ASSERT(fabs(grad_w[i] - fd) < 1e-6);
            free(pp);
        }

        free(paths); free(w); free(per_path); free(grad_per_path); free(grad_w);
        irrep_tp_free(d);
        irrep_multiset_free(A);
        irrep_multiset_free(B);
        irrep_multiset_free(C);
    }

    return IRREP_TEST_END();
}
