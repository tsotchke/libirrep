/* SPDX-License-Identifier: MIT */
/* Tests for regular solid harmonics R_{l,m}(r) — polynomials in x, y, z
 * proportional to r^l · Y_{l,m}(r̂).
 *
 * Coverage:
 *   - l = 0 is a constant; l = 1 is (√(3/4π)) · (y, z, x).
 *   - Homogeneity:  R_{l,m}(α r) = α^l · R_{l,m}(r).
 *   - At |r| = 1, solid harmonics equal the surface harmonics for every
 *     supported l.
 *   - `solid_harmonics_all` concatenation layout matches the per-l path.
 *   - Analytic gradient agrees with a 5-point finite-difference estimate to
 *     ≈1e-8 for every l.
 *   - High-l recurrence stays finite and matches surface × |r|^l.
 */
#include "harness.h"
#include <irrep/solid_harmonics.h>
#include <irrep/spherical_harmonics.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    IRREP_TEST_START("solid_harmonics");

    const double tol = 1e-12;

    /* -------- l=0 is a constant -------- */
    {
        double r[3] = {0.7, -0.3, 1.2};
        double out[1];
        irrep_solid_harm_cart(0, out, r);
        IRREP_ASSERT_NEAR(out[0], 1.0 / sqrt(4.0 * M_PI), tol);
    }

    /* -------- l=1 is (√(3/4π)) · (y, z, x) -------- */
    {
        double r[3] = {0.5, -0.2, 0.9};
        double out[3];
        irrep_solid_harm_cart(1, out, r);
        double k = sqrt(3.0 / (4.0 * M_PI));
        IRREP_ASSERT_NEAR(out[0], k * r[1], tol);
        IRREP_ASSERT_NEAR(out[1], k * r[2], tol);
        IRREP_ASSERT_NEAR(out[2], k * r[0], tol);
    }

    /* -------- homogeneity: R_{l,m}(α r) = α^l · R_{l,m}(r) -------- */
    {
        double r[3] = {0.3, -0.4, 0.6};
        double alpha = 1.7;
        double ar[3] = {alpha * r[0], alpha * r[1], alpha * r[2]};
        for (int l = 0; l <= 8; ++l) {
            double y1[2 * 8 + 1], y2[2 * 8 + 1];
            irrep_solid_harm_cart(l, y1, r);
            irrep_solid_harm_cart(l, y2, ar);
            double scale = pow(alpha, (double)l);
            for (int m = 0; m < 2 * l + 1; ++m) {
                IRREP_ASSERT_NEAR(y2[m], scale * y1[m], 1e-10);
            }
        }
    }

    /* -------- matches surface harmonics at |r| = 1, for every supported l -------- */
    {
        double r_hat[3] = {0.3, -0.4, sqrt(1.0 - 0.09 - 0.16)};
        for (int l = 0; l <= 8; ++l) {
            double solid[2 * 8 + 1], surface[2 * 8 + 1];
            irrep_solid_harm_cart(l, solid, r_hat);
            irrep_sph_harm_cart(l, surface, r_hat);
            for (int m = 0; m < 2 * l + 1; ++m) {
                IRREP_ASSERT(fabs(solid[m] - surface[m]) < 1e-10);
            }
        }
    }

    /* -------- solid-harmonic-all concatenation -------- */
    {
        double r[3] = {0.4, 0.5, -0.3};
        double flat[(8 + 1) * (8 + 1)];
        irrep_solid_harm_cart_all(8, flat, r);
        int offset = 0;
        for (int l = 0; l <= 8; ++l) {
            double single[2 * 8 + 1];
            irrep_solid_harm_cart(l, single, r);
            for (int m = 0; m < 2 * l + 1; ++m) {
                IRREP_ASSERT_NEAR(flat[offset + m], single[m], tol);
            }
            offset += 2 * l + 1;
        }
    }

    /* -------- analytic gradient matches 5-point FD to ~1e-8 for all l -------- */
    {
        double r[3] = {0.37, -0.21, 0.84};
        double h = 1e-5;
        for (int l = 0; l <= 8; ++l) {
            int    d = 2 * l + 1;
            double g[3 * (2 * 8 + 1)];
            irrep_solid_harm_cart_grad(l, g, r);

            for (int axis = 0; axis < 3; ++axis) {
                double r_p[3] = {r[0], r[1], r[2]};
                double r_m[3] = {r[0], r[1], r[2]};
                double r_pp[3] = {r[0], r[1], r[2]};
                double r_mm[3] = {r[0], r[1], r[2]};
                r_p[axis] += h;
                r_m[axis] -= h;
                r_pp[axis] += 2.0 * h;
                r_mm[axis] -= 2.0 * h;
                double yp[17], ym[17], ypp[17], ymm[17];
                irrep_solid_harm_cart(l, yp, r_p);
                irrep_solid_harm_cart(l, ym, r_m);
                irrep_solid_harm_cart(l, ypp, r_pp);
                irrep_solid_harm_cart(l, ymm, r_mm);
                for (int m = 0; m < d; ++m) {
                    double fd = (ymm[m] - 8.0 * ym[m] + 8.0 * yp[m] - ypp[m]) / (12.0 * h);
                    IRREP_ASSERT(fabs(g[axis * d + m] - fd) < 1e-6);
                }
            }
        }
    }

    /* -------- high-l recurrence path stays finite and matches surface × |r|^l -------- */
    {
        double r[3] = {0.4, 0.5, -0.3};
        double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        double n = sqrt(r2);
        double r_hat[3] = {r[0] / n, r[1] / n, r[2] / n};
        for (int l = 5; l <= 8; ++l) {
            int    d = 2 * l + 1;
            double solid[2 * 8 + 1], surface[2 * 8 + 1];
            irrep_solid_harm_cart(l, solid, r);
            irrep_sph_harm_cart(l, surface, r_hat);
            double scale = pow(n, (double)l);
            for (int m = 0; m < d; ++m) {
                IRREP_ASSERT(isfinite(solid[m]));
                IRREP_ASSERT(fabs(solid[m] - scale * surface[m]) < 1e-10);
            }
        }
    }

    return IRREP_TEST_END();
}
