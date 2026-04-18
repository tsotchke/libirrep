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
        for (int l = 0; l <= 5; ++l) {
            double ha[2 * 5 + 1], hb[2 * 5 + 1];
            irrep_sph_harm_cart(l, ha, r1);
            irrep_sph_harm_cart(l, hb, r2);
            double sum = 0.0;
            for (int i = 0; i < 2 * l + 1; ++i) sum += ha[i] * hb[i];
            sum *= 4.0 * M_PI / (2 * l + 1);
            double pl = irrep_legendre_assoc(l, 0, dot);
            IRREP_ASSERT_NEAR(sum, pl, 1e-10);
        }
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

    return IRREP_TEST_END();
}
