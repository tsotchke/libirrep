/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/spherical_harmonics.h>
#include <irrep/quadrature.h>
#include <complex.h>
#include <math.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

int main(void) {
    IRREP_TEST_START("spherical_harmonics");

    const double tol = 1e-12;

    /* ---- Y_0^0(θ, φ) = 1 / √(4π) ---- */
    {
        double _Complex y = irrep_sph_harm(0, 0, 0.5, 1.2);
        IRREP_ASSERT_NEAR(creal(y), 1.0 / sqrt(4.0 * M_PI), tol);
        IRREP_ASSERT_NEAR(cimag(y), 0.0, tol);
    }

    /* ---- Y_1^0 at north pole = √(3/(4π)) · cos 0 = √(3/(4π)) ---- */
    {
        double _Complex y = irrep_sph_harm(1, 0, 0.0, 0.7);
        IRREP_ASSERT_NEAR(creal(y), sqrt(3.0 / (4.0 * M_PI)), tol);
    }

    /* ---- Y_1^{real} in cartesian: Y_{1,0} ∝ z, Y_{1,±1} ∝ (x, y) ---- */
    {
        double theta = 0.5, phi = 0.3;
        double xv = sin(theta) * cos(phi);
        double yv = sin(theta) * sin(phi);
        double zv = cos(theta);
        double k = sqrt(3.0 / (4.0 * M_PI));
        IRREP_ASSERT_NEAR(irrep_sph_harm_real(1,  0, theta, phi), k * zv, tol);
        IRREP_ASSERT_NEAR(irrep_sph_harm_real(1, +1, theta, phi), k * xv, tol);
        IRREP_ASSERT_NEAR(irrep_sph_harm_real(1, -1, theta, phi), k * yv, tol);
    }

    /* ---- cartesian ≡ polar for a fully populated direction, every l,m ---- */
    {
        double theta = 0.73, phi = 1.44;
        double r_hat[3] = { sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };
        for (int l = 0; l <= 6; ++l) {
            double cart[2 * 6 + 1];
            irrep_sph_harm_cart(l, cart, r_hat);
            for (int m = -l; m <= l; ++m) {
                double polar = irrep_sph_harm_real(l, m, theta, phi);
                IRREP_ASSERT_NEAR(cart[m + l], polar, tol);
            }
        }
    }

    /* ---- Condon-Shortley: Y_l^{-m} = (-1)^m conj(Y_l^m) ---- */
    {
        double theta = 0.4, phi = 0.6;
        for (int l = 1; l <= 5; ++l) {
            for (int m = 1; m <= l; ++m) {
                double _Complex yp = irrep_sph_harm(l,  m, theta, phi);
                double _Complex yn = irrep_sph_harm(l, -m, theta, phi);
                double          s  = (m & 1) ? -1.0 : 1.0;
                IRREP_ASSERT(cabs(yn - s * conj(yp)) < tol);
            }
        }
    }

    /* ---- sum rule: Σ_m |Y_l^m|² = (2l+1)/(4π) ---- */
    {
        double theta = 0.4, phi = 0.6;
        for (int l = 0; l <= 8; ++l) {
            double sum = 0.0;
            for (int m = -l; m <= l; ++m) {
                double _Complex y = irrep_sph_harm(l, m, theta, phi);
                sum += cabs(y) * cabs(y);
            }
            IRREP_ASSERT_NEAR(sum, (2 * l + 1) / (4.0 * M_PI), 1e-10);
        }
    }

    /* ---- sum rule (real basis) via cartesian ---- */
    {
        double r_hat[3] = { 0.3, 0.4, sqrt(1.0 - 0.09 - 0.16) };
        for (int l = 0; l <= 6; ++l) {
            double cart[2 * 6 + 1];
            irrep_sph_harm_cart(l, cart, r_hat);
            double sum = 0.0;
            for (int i = 0; i < 2 * l + 1; ++i) sum += cart[i] * cart[i];
            IRREP_ASSERT_NEAR(sum, (2 * l + 1) / (4.0 * M_PI), 1e-10);
        }
    }

    /* ---- addition theorem:  (4π/(2l+1)) Σ_m Y^real(r1) Y^real(r2) = P_l(r1·r2) ---- */
    {
        double r1[3] = { 0.5, 0.3, 0.8 };
        double r2[3] = { -0.2, 0.7, 0.68 };
        double n1 = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
        double n2 = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
        for (int i = 0; i < 3; ++i) { r1[i] /= n1; r2[i] /= n2; }
        double dot = r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2];
        /* Tolerance here is *machine-precision* — a degradation to ~1e-10
         * would mean the real-SH normalisation constants have drifted (seen
         * when a downstream observed 2.5e-10 on Y_2^0(π/4) against the stale
         * release/1.0.0 pkg). */
        for (int l = 0; l <= 8; ++l) {
            double ha[2 * 8 + 1], hb[2 * 8 + 1];
            irrep_sph_harm_cart(l, ha, r1);
            irrep_sph_harm_cart(l, hb, r2);
            double sum = 0.0;
            for (int i = 0; i < 2 * l + 1; ++i) sum += ha[i] * hb[i];
            sum *= 4.0 * M_PI / (2 * l + 1);
            double pl = irrep_legendre_assoc(l, 0, dot);
            IRREP_ASSERT_NEAR(sum, pl, 5e-14);
        }
    }

    /* ---- Closed-form Y_2^0(π/4) sanity check ----
     * Downstream cross-check against hand-expanded √(5/(16π)) · (3 cos²θ − 1)
     * flagged a 2.5e-10 drift on an older binary. This assertion locks the
     * current bit-exact behaviour in. */
    {
        double theta = M_PI / 4.0;
        double phi   = 0.0;
        double c     = cos(theta);
        double expected = sqrt(5.0 / (16.0 * M_PI)) * (3.0 * c * c - 1.0);

        /* Polar form. */
        IRREP_ASSERT(fabs(irrep_sph_harm_real(2, 0, theta, phi) - expected) < 1e-14);

        /* Cartesian form — should agree bit-for-bit with polar. */
        double r_hat[3] = { sin(theta) * cos(phi),
                            sin(theta) * sin(phi),
                            cos(theta) };
        double buf[5];
        irrep_sph_harm_cart(2, buf, r_hat);
        IRREP_ASSERT(fabs(buf[2] - expected) < 1e-14);
    }

    /* ---- radial derivative of Y(r̂) is zero ---- */
    {
        double r_hat[3] = { 0.4, 0.6, sqrt(1.0 - 0.16 - 0.36) };
        for (int l = 1; l <= 3; ++l) {
            double grad[3 * (2 * 3 + 1)];
            irrep_sph_harm_cart_grad(l, grad, r_hat);
            int d = 2 * l + 1;
            for (int m = 0; m < d; ++m) {
                double radial = r_hat[0] * grad[0 * d + m]
                              + r_hat[1] * grad[1 * d + m]
                              + r_hat[2] * grad[2 * d + m];
                IRREP_ASSERT(fabs(radial) < 1e-8);
            }
        }
    }

    /* ---- complex-to-real matrix is unitary ---- */
    {
        for (int l = 0; l <= 4; ++l) {
            int d = 2 * l + 1;
            double _Complex U[(2 * 4 + 1) * (2 * 4 + 1)];
            irrep_sph_harm_complex_to_real(l, U);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    double _Complex s = 0.0;
                    for (int k = 0; k < d; ++k) s += U[i * d + k] * conj(U[j * d + k]);
                    double target = (i == j) ? 1.0 : 0.0;
                    IRREP_ASSERT(cabs(s - target) < 1e-12);
                }
            }
        }
    }

    /* ---- sph_harm_cart_all layout matches concatenated sph_harm_cart ---- */
    {
        double r_hat[3] = { 0.3, 0.4, sqrt(1.0 - 0.25) };
        double flat[(6 + 1) * (6 + 1)];
        irrep_sph_harm_cart_all(6, flat, r_hat);
        int offset = 0;
        for (int l = 0; l <= 6; ++l) {
            double single[2 * 6 + 1];
            irrep_sph_harm_cart(l, single, r_hat);
            for (int m = 0; m < 2 * l + 1; ++m) {
                IRREP_ASSERT_NEAR(flat[offset + m], single[m], 1e-14);
            }
            offset += 2 * l + 1;
        }
    }

    /* ---- poles: only Y_{l,0} is nonzero; Y_{l,±m}(±ẑ) = 0 for m > 0 ---- */
    {
        double north[3] = { 0, 0,  1 };
        double south[3] = { 0, 0, -1 };
        for (int l = 1; l <= 5; ++l) {
            double cart[2 * 5 + 1];
            irrep_sph_harm_cart(l, cart, north);
            for (int m = 1; m <= l; ++m) {
                IRREP_ASSERT(fabs(cart[l + m]) < 1e-12);
                IRREP_ASSERT(fabs(cart[l - m]) < 1e-12);
            }
            /* Y_{l,0}(north) = N P_l(1) = √((2l+1)/(4π)) */
            IRREP_ASSERT_NEAR(cart[l], sqrt((2 * l + 1) / (4.0 * M_PI)), tol);
            /* Y_{l,0}(south) = √((2l+1)/(4π)) · (-1)^l */
            irrep_sph_harm_cart(l, cart, south);
            double sgn = (l & 1) ? -1.0 : 1.0;
            IRREP_ASSERT_NEAR(cart[l], sgn * sqrt((2 * l + 1) / (4.0 * M_PI)), tol);
        }
    }

    /* ---- orthonormality via Gauss-Legendre (θ) × uniform (φ) tensor product
     *      ∫ Y_l^m Y_{l'}^{m'}* dΩ = δ_{ll'} δ_{mm'} ---- */
    {
        int n_theta = 32;
        double nodes[32], weights[32];
        irrep_gauss_legendre(n_theta, nodes, weights);
        int n_phi = 40;
        double dphi = 2.0 * M_PI / (double)n_phi;
        for (int l1 = 0; l1 <= 3; ++l1) {
            for (int l2 = 0; l2 <= 3; ++l2) {
                for (int m1 = -l1; m1 <= l1; ++m1) {
                    for (int m2 = -l2; m2 <= l2; ++m2) {
                        double _Complex s = 0.0;
                        for (int i = 0; i < n_theta; ++i) {
                            double cos_theta = nodes[i];
                            double theta = acos(cos_theta);
                            double w_t = weights[i];
                            for (int j = 0; j < n_phi; ++j) {
                                double phi = (j + 0.5) * dphi;
                                double _Complex y1 = irrep_sph_harm(l1, m1, theta, phi);
                                double _Complex y2 = irrep_sph_harm(l2, m2, theta, phi);
                                s += w_t * dphi * y1 * conj(y2);
                            }
                        }
                        double expected = (l1 == l2 && m1 == m2) ? 1.0 : 0.0;
                        IRREP_ASSERT(cabs(s - expected) < 1e-10);
                    }
                }
            }
        }
    }

    /* ---- f32 wrapper agrees with double within float precision ---- */
    {
        float r_hat[3] = { 0.3f, 0.4f, (float)sqrt(1.0 - 0.25) };
        float cart_f[2 * 4 + 1];
        double r_hat_d[3] = { r_hat[0], r_hat[1], r_hat[2] };
        double cart_d[2 * 4 + 1];
        for (int l = 0; l <= 4; ++l) {
            irrep_sph_harm_cart_f32(l, cart_f, r_hat);
            irrep_sph_harm_cart(l, cart_d, r_hat_d);
            for (int i = 0; i < 2 * l + 1; ++i) {
                IRREP_ASSERT(fabs((double)cart_f[i] - cart_d[i]) < 1e-5);
            }
        }
    }

    /* -------- NEON batch kernel is bit-exact against scalar --------
     *
     * Sweeps l_max ∈ {0..6} over N = 23 unit vectors, deliberately seeded
     * with pole cases and near-pole inputs to cover the r_xy → 0 branch.
     * Odd N exercises the tail path. Dispatched kernel (NEON on aarch64,
     * scalar elsewhere) must match the per-edge scalar call byte-for-byte. */
    {
        const size_t N = 23;
        double rhats[23 * 3];
        /* Deterministic generic inputs for the first 17 entries. */
        for (size_t i = 0; i < 17; ++i) {
            double t = 0.37 * (double)i + 0.11;
            double x = sin(t), y = cos(0.7 * t), z = sin(1.3 * t);
            double n = sqrt(x*x + y*y + z*z);
            rhats[3*i+0] = x / n;
            rhats[3*i+1] = y / n;
            rhats[3*i+2] = z / n;
        }
        /* Edge-case stressors: exact poles and r_xy ≪ 1. */
        rhats[17*3+0] = 0.0;    rhats[17*3+1] = 0.0;    rhats[17*3+2] =  1.0;  /* +ẑ pole */
        rhats[18*3+0] = 0.0;    rhats[18*3+1] = 0.0;    rhats[18*3+2] = -1.0;  /* −ẑ pole */
        {   /* near-pole: r_xy just above the scalar threshold (1e-14). */
            double rxy = 1e-12;
            double z   = sqrt(1.0 - rxy * rxy);
            rhats[19*3+0] = rxy; rhats[19*3+1] = 0.0; rhats[19*3+2] = z;
        }
        {   /* near-pole: r_xy just below threshold — scalar falls into (cp, sp) = (1, 0). */
            double rxy = 1e-16;
            double z   = sqrt(1.0 - rxy * rxy);
            rhats[20*3+0] = rxy; rhats[20*3+1] = 0.0; rhats[20*3+2] = z;
        }
        rhats[21*3+0] = 1.0;    rhats[21*3+1] = 0.0;    rhats[21*3+2] = 0.0;  /* +x̂ equator */
        rhats[22*3+0] = 0.0;    rhats[22*3+1] = 1.0;    rhats[22*3+2] = 0.0;  /* +ŷ equator */

        for (int l_max = 0; l_max <= 6; ++l_max) {
            int block = (l_max + 1) * (l_max + 1);
            double *batch_out = calloc(N * (size_t)block, sizeof(double));
            double *ref_out   = calloc(N * (size_t)block, sizeof(double));
            IRREP_ASSERT(batch_out != NULL && ref_out != NULL);

            irrep_sph_harm_cart_all_batch(l_max, N, rhats, batch_out);
            for (size_t i = 0; i < N; ++i) {
                irrep_sph_harm_cart_all(l_max, ref_out + i * (size_t)block,
                                        rhats + i * 3);
            }
            for (size_t i = 0; i < N * (size_t)block; ++i) {
                IRREP_ASSERT(batch_out[i] == ref_out[i]);
            }
            free(batch_out); free(ref_out);
        }
    }

    /* -------- Batched gradient matches the per-l gradient --------
     * For each edge and axis and l, irrep_sph_harm_cart_all_grad_batch should
     * produce the same values as a direct per-l call to _grad. Bit-exact
     * (single call path on both sides). */
    {
        const int l_max = 4;
        const int block = (l_max + 1) * (l_max + 1);
        const size_t N = 3;
        double r_hats[3 * 3] = {
            0.6, -0.3,  0.7,
           -0.4,  0.8,  0.2,
            0.0,  0.0,  1.0,
        };
        /* Normalise. */
        for (size_t e = 0; e < N; ++e) {
            double *v = r_hats + e * 3;
            double n = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            v[0] /= n; v[1] /= n; v[2] /= n;
        }
        double batch_grad[3 * 3 * (4 + 1) * (4 + 1)];
        irrep_sph_harm_cart_all_grad_batch(l_max, N, r_hats, batch_grad);

        double ref_grad[3 * (2 * 4 + 1)];
        for (size_t e = 0; e < N; ++e) {
            int off = 0;
            for (int l = 0; l <= l_max; ++l) {
                int d = 2 * l + 1;
                irrep_sph_harm_cart_grad(l, ref_grad, r_hats + e * 3);
                for (int axis = 0; axis < 3; ++axis) {
                    for (int i = 0; i < d; ++i) {
                        double got = batch_grad[(e * 3 + axis) * block + off + i];
                        double ref = ref_grad[axis * d + i];
                        IRREP_ASSERT(got == ref);
                    }
                }
                off += d;
            }
        }
    }

    return IRREP_TEST_END();
}
