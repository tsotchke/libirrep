/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/quadrature.h>
#include <math.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

int main(void) {
    IRREP_TEST_START("quadrature");

    /* -------- Gauss-Legendre weights sum to 2 -------- */
    for (int n = 1; n <= 16; ++n) {
        double nodes[32], weights[32];
        IRREP_ASSERT(irrep_gauss_legendre(n, nodes, weights));
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += weights[i];
        IRREP_ASSERT(fabs(s - 2.0) < 1e-12);
    }

    /* -------- Gauss-Legendre integrates polynomial degree ≤ 2n − 1 exactly -------- */
    for (int n = 1; n <= 8; ++n) {
        double nodes[16], weights[16];
        irrep_gauss_legendre(n, nodes, weights);
        for (int d = 0; d <= 2 * n - 1; ++d) {
            double s = 0.0;
            for (int i = 0; i < n; ++i) s += weights[i] * pow(nodes[i], d);
            double expected = (d & 1) ? 0.0 : 2.0 / (d + 1);    /* ∫_{-1}^{1} x^d dx */
            IRREP_ASSERT(fabs(s - expected) < 1e-10);
        }
    }

    /* -------- Lebedev: sizes match the Lebedev 1999 counts -------- */
    IRREP_ASSERT(irrep_lebedev_size(3) == 6);
    IRREP_ASSERT(irrep_lebedev_size(5) == 14);
    IRREP_ASSERT(irrep_lebedev_size(7) == 26);
    IRREP_ASSERT(irrep_lebedev_size(9) == 0);    /* not yet tabulated */

    /* -------- Lebedev: weights sum to 1 (unit-sphere-normalized) -------- */
    for (int order = 3; order <= 7; order += 2) {
        int sz = irrep_lebedev_size(order);
        if (sz == 0) continue;
        double *buf = malloc((size_t)sz * 4 * sizeof(double));
        IRREP_ASSERT(irrep_lebedev_fill(order, buf));
        double s = 0.0;
        for (int i = 0; i < sz; ++i) s += buf[i * 4 + 3];
        IRREP_ASSERT(fabs(s - 1.0) < 1e-12);
        /* Each point must be on the unit sphere */
        for (int i = 0; i < sz; ++i) {
            double nx = buf[i*4+0], ny = buf[i*4+1], nz = buf[i*4+2];
            double n2 = nx*nx + ny*ny + nz*nz;
            IRREP_ASSERT(fabs(n2 - 1.0) < 1e-12);
        }
        free(buf);
    }

    /* -------- Lebedev rule exactness: ∫ x² dΩ = 4π/3, ∫ x⁴ dΩ = 4π/5, etc.
     *          Our weights are normalized so Σ w_i f(x_i) gives ∫/(4π).         */
    for (int order = 3; order <= 7; order += 2) {
        int sz = irrep_lebedev_size(order);
        double *buf = malloc((size_t)sz * 4 * sizeof(double));
        irrep_lebedev_fill(order, buf);

        /* ⟨x²⟩_{S²} = 1/3 */
        double s = 0.0;
        for (int i = 0; i < sz; ++i) {
            double x = buf[i * 4 + 0], w = buf[i * 4 + 3];
            s += w * x * x;
        }
        IRREP_ASSERT(fabs(s - 1.0 / 3.0) < 1e-12);

        /* ⟨x² + y² + z²⟩ = 1 */
        double s_full = 0.0;
        for (int i = 0; i < sz; ++i) {
            double x = buf[i*4+0], y = buf[i*4+1], z = buf[i*4+2], w = buf[i*4+3];
            s_full += w * (x*x + y*y + z*z);
        }
        IRREP_ASSERT(fabs(s_full - 1.0) < 1e-12);
        free(buf);
    }

    /* -------- Tensor-product sphere quadrature: exactness across all orders -------- */
    {
        for (int D = 2; D <= 20; ++D) {
            int sz = irrep_quadrature_sphere_size(D);
            double *buf = malloc((size_t)sz * 4 * sizeof(double));
            IRREP_ASSERT(irrep_quadrature_sphere_fill(D, buf));

            /* Weights sum to 1 and each point is unit. */
            double s_w = 0.0;
            for (int i = 0; i < sz; ++i) {
                double x = buf[i*4+0], y = buf[i*4+1], z = buf[i*4+2], w = buf[i*4+3];
                s_w += w;
                double n2 = x*x + y*y + z*z;
                IRREP_ASSERT(fabs(n2 - 1.0) < 1e-12);
            }
            IRREP_ASSERT(fabs(s_w - 1.0) < 1e-12);

            /* Integrate x^2 — expect 1/3 exactly for D ≥ 2. */
            double s_x2 = 0.0;
            for (int i = 0; i < sz; ++i) {
                double x = buf[i*4+0], w = buf[i*4+3];
                s_x2 += w * x * x;
            }
            IRREP_ASSERT(fabs(s_x2 - 1.0 / 3.0) < 1e-12);

            /* Integrate x^4 — expect 1/5 for D ≥ 4. */
            if (D >= 4) {
                double s_x4 = 0.0;
                for (int i = 0; i < sz; ++i) {
                    double x = buf[i*4+0], w = buf[i*4+3];
                    s_x4 += w * x * x * x * x;
                }
                IRREP_ASSERT(fabs(s_x4 - 1.0 / 5.0) < 1e-12);
            }

            free(buf);
        }
    }

    /* -------- Lebedev order 5: integrate degree-4 monomials exactly -------- */
    /* ⟨x⁴⟩_{S²} = 1/5;  ⟨x²y²⟩ = 1/15 */
    {
        int sz = irrep_lebedev_size(5);
        double buf[14 * 4];
        irrep_lebedev_fill(5, buf);
        double s4 = 0.0, s22 = 0.0;
        for (int i = 0; i < sz; ++i) {
            double x = buf[i*4+0], y = buf[i*4+1], w = buf[i*4+3];
            s4  += w * x*x*x*x;
            s22 += w * x*x*y*y;
        }
        IRREP_ASSERT(fabs(s4  - 1.0 / 5.0 ) < 1e-12);
        IRREP_ASSERT(fabs(s22 - 1.0 / 15.0) < 1e-12);
    }

    return IRREP_TEST_END();
}
