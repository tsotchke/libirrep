/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/radial.h>
#include <irrep/quadrature.h>
#include <math.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

int main(void) {
    IRREP_TEST_START("radial");

    const double tol = 1e-12;

    /* -------- Bessel RBF orthogonality on [0, r_cut] with r² weight -------- */
    /* Use Gauss-Legendre on [0, r_cut]: change of variables x = 2r/r_cut − 1. */
    {
        double r_cut = 2.5;
        int n = 64;
        double nodes[64], weights[64];
        irrep_gauss_legendre(n, nodes, weights);

        for (int n_a = 1; n_a <= 4; ++n_a) {
            for (int n_b = 1; n_b <= 4; ++n_b) {
                double s = 0.0;
                for (int i = 0; i < n; ++i) {
                    double r = 0.5 * r_cut * (nodes[i] + 1.0);
                    double jac = 0.5 * r_cut;
                    double pa = irrep_rbf_bessel(n_a, r, r_cut);
                    double pb = irrep_rbf_bessel(n_b, r, r_cut);
                    s += weights[i] * jac * pa * pb * r * r;
                }
                double expected = (n_a == n_b) ? 1.0 : 0.0;
                IRREP_ASSERT(fabs(s - expected) < 1e-8);
            }
        }
    }

    /* -------- Bessel RBF batched matches single -------- */
    {
        double buf[4];
        irrep_rbf_bessel_all(4, buf, 0.5, 2.0);
        for (int n = 1; n <= 4; ++n) {
            IRREP_ASSERT_NEAR(buf[n - 1], irrep_rbf_bessel(n, 0.5, 2.0), tol);
        }
    }

    /* -------- Gaussian RBF peak at μ, grid samples at right positions -------- */
    {
        IRREP_ASSERT_NEAR(irrep_rbf_gaussian(0.5, 0.5, 0.1), 1.0, tol);
        IRREP_ASSERT(irrep_rbf_gaussian(0.7, 0.5, 0.1) < 1.0);

        double grid[5];
        irrep_rbf_gaussian_grid(5, grid, 0.0, 0.0, 1.0, 0.1);
        IRREP_ASSERT_NEAR(grid[0], 1.0, tol);
        IRREP_ASSERT(grid[1] < grid[0]);
        IRREP_ASSERT(grid[4] < grid[1]);
    }

    /* -------- cosine cutoff: c(0) = 1, c(r_cut) = 0, c'(r_cut^−) = 0 -------- */
    {
        IRREP_ASSERT_NEAR(irrep_cutoff_cosine(0.0, 2.0), 1.0, tol);
        IRREP_ASSERT_NEAR(irrep_cutoff_cosine(2.0, 2.0), 0.0, tol);
        /* c(r_cut/2) = 1/2 */
        IRREP_ASSERT_NEAR(irrep_cutoff_cosine(1.0, 2.0), 0.5, tol);
        /* c'(0) = 0 */
        IRREP_ASSERT_NEAR(irrep_cutoff_cosine_d(0.0, 2.0), 0.0, tol);
        /* c'(r_cut) → 0 just below */
        IRREP_ASSERT(fabs(irrep_cutoff_cosine_d(1.999, 2.0)) < 2e-3);
    }

    /* -------- cosine cutoff derivative matches finite difference -------- */
    {
        double r_cut = 2.0, r = 0.7, h = 1e-6;
        double fd = (irrep_cutoff_cosine(r + h, r_cut)
                   - irrep_cutoff_cosine(r - h, r_cut)) / (2 * h);
        IRREP_ASSERT(fabs(fd - irrep_cutoff_cosine_d(r, r_cut)) < 1e-8);
    }

    /* -------- polynomial cutoff: f_p(0) = 1, f_p(r_cut) = 0, smooth at r_cut -------- */
    {
        for (int p = 1; p <= 6; ++p) {
            IRREP_ASSERT_NEAR(irrep_cutoff_polynomial(0.0, 2.0, p), 1.0, tol);
            /* f_p(r→r_cut^−) → 0 */
            IRREP_ASSERT(fabs(irrep_cutoff_polynomial(1.999, 2.0, p)) < 1e-4);
            /* f_p(r ≥ r_cut) = 0 */
            IRREP_ASSERT_NEAR(irrep_cutoff_polynomial(2.0, 2.0, p), 0.0, tol);
        }
    }

    /* -------- polynomial cutoff derivative matches FD -------- */
    {
        double r_cut = 2.0;
        double h = 1e-6;
        for (int p = 2; p <= 6; ++p) {
            double r = 1.3;
            double fd = (irrep_cutoff_polynomial(r + h, r_cut, p)
                       - irrep_cutoff_polynomial(r - h, r_cut, p)) / (2 * h);
            double an = irrep_cutoff_polynomial_d(r, r_cut, p);
            IRREP_ASSERT(fabs(fd - an) < 1e-7);
        }
    }

    return IRREP_TEST_END();
}
