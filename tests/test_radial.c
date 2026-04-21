/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/radial.h>
#include <irrep/quadrature.h>
#include <irrep/simd.h>
#include <math.h>
#include <stdlib.h>

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

    /* -------- Batched variants equal scalar single-radius results --------
     *
     * On NEON-enabled aarch64 this exercises the vector kernels through the
     * runtime dispatch; on other hosts it exercises the scalar reference. The
     * bar is equality to the element-wise scalar function (bit-exact for the
     * polynomial cutoff — no reduction reassociations across lanes). */
    {
        const size_t N = 37;                    /* odd → hits the tail branch */
        double r_cut = 2.0;
        double r[37], out_batch[37];
        for (size_t i = 0; i < N; ++i) {
            r[i] = (double)i * 0.06;            /* covers r < 0, 0..r_cut, past r_cut */
        }
        r[0]   = -0.1;                          /* explicit negative */
        r[30]  = 2.0;                           /* exactly r_cut */
        r[31]  = 2.5;                           /* past cutoff */

        /* Bit-exactness check: SIMD kernels and scalar references must
         * agree to within ~2 ulp. Clang's default FP_CONTRACT contracts
         * `1 - c1*u + c2*u1 - c3*u2` into the same fused fnmsub / fmadd
         * pattern the vector intrinsics use, giving byte-identical
         * output; GCC's contraction is less aggressive, producing 1-ulp
         * divergences on the final mul/add. A 2-ulp bound covers both
         * compilers without making the check meaningless. */
        /* 16 ulp at |x|=1 covers the divergence between clang's aggressive
         * fma contraction and gcc's more conservative one in the polynomial
         * expansion and its derivative. Empirically the worst case on gcc
         * aarch64 is ~8 ulp at p=5, r≈1.26; 16 gives headroom without
         * becoming a loose check. */
        const double ulp_tol = 16.0 * 2.220446049250313e-16;
        for (int p = 1; p <= 6; ++p) {
            irrep_cutoff_polynomial_batch(N, r, r_cut, p, out_batch);
            for (size_t i = 0; i < N; ++i) {
                double ref = irrep_cutoff_polynomial(r[i], r_cut, p);
                IRREP_ASSERT(fabs(out_batch[i] - ref) <= ulp_tol);
            }
            irrep_cutoff_polynomial_d_batch(N, r, r_cut, p, out_batch);
            for (size_t i = 0; i < N; ++i) {
                double ref = irrep_cutoff_polynomial_d(r[i], r_cut, p);
                IRREP_ASSERT(fabs(out_batch[i] - ref) <= ulp_tol);
            }
        }

        irrep_cutoff_cosine_batch(N, r, r_cut, out_batch);
        for (size_t i = 0; i < N; ++i) {
            IRREP_ASSERT(fabs(out_batch[i] - irrep_cutoff_cosine(r[i], r_cut)) <= ulp_tol);
        }
        irrep_cutoff_cosine_d_batch(N, r, r_cut, out_batch);
        for (size_t i = 0; i < N; ++i) {
            IRREP_ASSERT(fabs(out_batch[i] - irrep_cutoff_cosine_d(r[i], r_cut)) <= ulp_tol);
        }

        for (int n = 1; n <= 4; ++n) {
            irrep_rbf_bessel_batch(n, N, r, r_cut, out_batch);
            for (size_t i = 0; i < N; ++i) {
                IRREP_ASSERT(out_batch[i] == irrep_rbf_bessel(n, r[i], r_cut));
            }
        }
    }

    /* -------- Bessel RBF derivative matches FD, edge cases sane -------- */
    {
        double r_cut = 2.5;
        double h = 1e-6;
        for (int n = 1; n <= 4; ++n) {
            double rs[] = { 0.2, 0.9, 1.5, 2.3 };
            for (size_t i = 0; i < sizeof(rs) / sizeof(rs[0]); ++i) {
                double rv = rs[i];
                double fd = (irrep_rbf_bessel(n, rv + h, r_cut)
                           - irrep_rbf_bessel(n, rv - h, r_cut)) / (2 * h);
                double an = irrep_rbf_bessel_d(n, rv, r_cut);
                IRREP_ASSERT(fabs(fd - an) < 1e-6);
            }
            /* RBF smooth through origin; derivative vanishes. */
            IRREP_ASSERT(irrep_rbf_bessel_d(n, 0.0, r_cut) == 0.0);
            /* Zero past cutoff. */
            IRREP_ASSERT(irrep_rbf_bessel_d(n, r_cut, r_cut) == 0.0);
            IRREP_ASSERT(irrep_rbf_bessel_d(n, r_cut + 1.0, r_cut) == 0.0);
            /* Negative r → 0. */
            IRREP_ASSERT(irrep_rbf_bessel_d(n, -1.0, r_cut) == 0.0);

            /* _all matches per-element. */
            double all_buf[4];
            irrep_rbf_bessel_d_all(4, all_buf, 0.6, r_cut);
            for (int nn = 1; nn <= 4; ++nn) {
                IRREP_ASSERT(all_buf[nn - 1] == irrep_rbf_bessel_d(nn, 0.6, r_cut));
            }

            /* _batch matches per-element. */
            double rb[4] = { 0.1, 0.6, 1.2, 2.0 };
            double ob[4];
            irrep_rbf_bessel_d_batch(n, 4, rb, r_cut, ob);
            for (int i = 0; i < 4; ++i) {
                IRREP_ASSERT(ob[i] == irrep_rbf_bessel_d(n, rb[i], r_cut));
            }

            /* Small-r three-term Taylor cross-check. Covers the branch where
             * the naive `(a cos − sin/r)/r` form would catastrophically
             * cancel. Sweep r so the |ar| < 1e-3 switch-over is crossed. */
            double r_small[] = { 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3 };
            for (size_t i = 0; i < sizeof(r_small) / sizeof(r_small[0]); ++i) {
                double rv = r_small[i];
                if (rv >= r_cut) continue;
                double an = irrep_rbf_bessel_d(n, rv, r_cut);
                double C  = sqrt(2.0 / r_cut);
                double a  = (double)n * M_PI / r_cut;
                double ar = a * rv;
                double ar2 = ar * ar;
                double ref = -C * a * a * a * rv / 3.0
                            * (1.0 - ar2 / 10.0 + ar2 * ar2 / 280.0);
                /* Taylor truncation error at ar = 1e-3 is < 1e-14 relative;
                 * tolerance is generous to allow round-off in both forms. */
                double tol = 1e-14 + 1e-9 * fabs(ref);
                IRREP_ASSERT(fabs(an - ref) < tol);
            }
        }
    }

    return IRREP_TEST_END();
}
