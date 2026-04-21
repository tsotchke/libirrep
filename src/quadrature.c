/* SPDX-License-Identifier: MIT */
/* M8: quadrature — Gauss-Legendre (any n) + Lebedev (orders 3/5/7 only;
 * higher orders deferred to the Lebedev-Laikov 1999 data import). */

#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include <irrep/quadrature.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT2
#  define M_SQRT2 1.41421356237309504880
#endif

/* -------------------------------------------------------------------------- *
 * Gauss-Legendre on [−1, 1] via Newton iteration on P_n(x).                  *
 * Nodes returned in descending order (Chebyshev-guess convention).           *
 * -------------------------------------------------------------------------- */

bool irrep_gauss_legendre(int n, double *nodes, double *weights) {
    if (n <= 0 || !nodes || !weights) return false;

    if (n == 1) {
        nodes[0] = 0.0;
        weights[0] = 2.0;
        return true;
    }

    const int    max_iter = 100;
    const double tol      = 1e-15;

    for (int k = 0; k < n; ++k) {
        double x = cos(M_PI * ((double)k + 0.75) / ((double)n + 0.5));

        double p0 = 1.0, p1 = x, dp = 0.0;
        for (int iter = 0; iter < max_iter; ++iter) {
            p0 = 1.0; p1 = x;
            for (int j = 1; j < n; ++j) {
                double p2 = ((2.0 * j + 1.0) * x * p1 - (double)j * p0) / ((double)j + 1.0);
                p0 = p1;
                p1 = p2;
            }
            /* p1 = P_n(x), p0 = P_{n-1}(x); derivative via the standard identity. */
            dp = (double)n * (x * p1 - p0) / (x * x - 1.0);
            double dx = p1 / dp;
            x -= dx;
            if (fabs(dx) < tol) break;
        }

        nodes[k]   = x;
        weights[k] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
    return true;
}

/* -------------------------------------------------------------------------- *
 * Lebedev quadrature on S².                                                  *
 *                                                                            *
 * Rules 3 / 5 / 7 are hand-encoded here. Remaining rules                     *
 * (9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41) land with the      *
 * full Lebedev-Laikov 1999 data import in a future polish pass.              *
 * -------------------------------------------------------------------------- */

int irrep_lebedev_size(int order) {
    switch (order) {
        case 3:  return 6;
        case 5:  return 14;
        case 7:  return 26;
        default: return 0;
    }
}

static void add_a1_(double *buf, size_t *idx, double w) {
    /* 6 axis points. */
    const double axes[6][3] = {
        { +1, 0, 0 }, { -1, 0, 0 },
        { 0, +1, 0 }, { 0, -1, 0 },
        { 0, 0, +1 }, { 0, 0, -1 },
    };
    for (int i = 0; i < 6; ++i) {
        buf[*idx * 4 + 0] = axes[i][0];
        buf[*idx * 4 + 1] = axes[i][1];
        buf[*idx * 4 + 2] = axes[i][2];
        buf[*idx * 4 + 3] = w;
        (*idx)++;
    }
}

static void add_a2_(double *buf, size_t *idx, double w) {
    /* 12 edge points: (0, ±1/√2, ±1/√2) and the two cyclic positions for 0. */
    const double r = 1.0 / M_SQRT2;
    const double pts[12][3] = {
        { 0, +r, +r }, { 0, +r, -r }, { 0, -r, +r }, { 0, -r, -r },
        { +r, 0, +r }, { +r, 0, -r }, { -r, 0, +r }, { -r, 0, -r },
        { +r, +r, 0 }, { +r, -r, 0 }, { -r, +r, 0 }, { -r, -r, 0 },
    };
    for (int i = 0; i < 12; ++i) {
        buf[*idx * 4 + 0] = pts[i][0];
        buf[*idx * 4 + 1] = pts[i][1];
        buf[*idx * 4 + 2] = pts[i][2];
        buf[*idx * 4 + 3] = w;
        (*idx)++;
    }
}

static void add_a3_(double *buf, size_t *idx, double w) {
    /* 8 corner points: ±1/√3 in each component. */
    const double r = 1.0 / 1.7320508075688772935;    /* 1/√3 */
    for (int sx = 0; sx < 2; ++sx) {
        for (int sy = 0; sy < 2; ++sy) {
            for (int sz = 0; sz < 2; ++sz) {
                double x = sx ? -r : r;
                double y = sy ? -r : r;
                double z = sz ? -r : r;
                buf[*idx * 4 + 0] = x;
                buf[*idx * 4 + 1] = y;
                buf[*idx * 4 + 2] = z;
                buf[*idx * 4 + 3] = w;
                (*idx)++;
            }
        }
    }
}

bool irrep_lebedev_fill(int order, double *xyz_weights) {
    if (!xyz_weights) return false;
    size_t idx = 0;
    switch (order) {
        case 3:
            add_a1_(xyz_weights, &idx, 1.0 / 6.0);
            return true;
        case 5:
            add_a1_(xyz_weights, &idx, 1.0 / 15.0);
            add_a3_(xyz_weights, &idx, 3.0 / 40.0);
            return true;
        case 7:
            add_a1_(xyz_weights, &idx, 1.0 / 21.0);
            add_a2_(xyz_weights, &idx, 4.0 / 105.0);
            add_a3_(xyz_weights, &idx, 9.0 / 280.0);
            return true;
        default:
            return false;
    }
}

/* -------------------------------------------------------------------------- *
 * Tensor-product S² quadrature: Gauss-Legendre in cos θ × uniform in φ.      *
 * Covers arbitrary exactness at the cost of ≈ 2× more points than Lebedev.   *
 * -------------------------------------------------------------------------- */

#include <stdlib.h>

int irrep_quadrature_sphere_size(int exactness_deg) {
    if (exactness_deg < 0) return 0;
    int n_theta = exactness_deg / 2 + 1;
    int n_phi   = exactness_deg + 1;
    if (n_theta < 1) n_theta = 1;
    if (n_phi   < 1) n_phi   = 1;
    return n_theta * n_phi;
}

bool irrep_quadrature_sphere_fill(int exactness_deg, double *xyz_weights) {
    if (!xyz_weights || exactness_deg < 0) return false;
    int n_theta = exactness_deg / 2 + 1;
    int n_phi   = exactness_deg + 1;
    if (n_theta < 1) n_theta = 1;
    if (n_phi   < 1) n_phi   = 1;

    double *nodes   = malloc((size_t)n_theta * sizeof(double));
    double *weights = malloc((size_t)n_theta * sizeof(double));
    if (!nodes || !weights) { free(nodes); free(weights); return false; }
    if (!irrep_gauss_legendre(n_theta, nodes, weights)) {
        free(nodes); free(weights); return false;
    }

    double dphi    = 2.0 * M_PI / (double)n_phi;
    double inv_4pi = 1.0 / (4.0 * M_PI);
    int idx = 0;
    for (int i = 0; i < n_theta; ++i) {
        double ct  = nodes[i];
        double st2 = 1.0 - ct * ct;
        double st  = st2 > 0.0 ? sqrt(st2) : 0.0;
        double wt  = weights[i];
        for (int j = 0; j < n_phi; ++j) {
            double ph = ((double)j + 0.5) * dphi;
            xyz_weights[idx * 4 + 0] = st * cos(ph);
            xyz_weights[idx * 4 + 1] = st * sin(ph);
            xyz_weights[idx * 4 + 2] = ct;
            xyz_weights[idx * 4 + 3] = wt * dphi * inv_4pi;
            idx++;
        }
    }
    free(nodes); free(weights);
    return true;
}
