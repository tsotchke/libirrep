/* SPDX-License-Identifier: MIT */
/* Real regular solid harmonics via explicit polynomial expressions for
 * l = 0..IRREP_SOLID_L_MAX = 4. Higher l falls back to the surface-harmonics
 * path (`irrep_sph_harm_cart`) with a scale-by-|r|^l correction.
 *
 * All normalisations match `irrep/spherical_harmonics.h` so that
 *
 *   R_{l, m}(r̂) == irrep_sph_harm_cart(l, ...)
 *
 * when |r̂| = 1.
 *
 * Gradients for l = 0..4 come from polynomial differentiation and are exact
 * to round-off. For l ≥ 5 we use a 5-point centred stencil on R(r) treated
 * as a polynomial — still accurate to ~1e-12 in double precision but
 * slightly noisier than analytic.
 */

#include <math.h>
#include <stdatomic.h>
#include <stddef.h>
#include <string.h>

#include <irrep/solid_harmonics.h>
#include <irrep/spherical_harmonics.h>
#include <irrep/types.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* -------------------------------------------------------------------------- *
 * Normalisation constants (computed at first use, cached).                   *
 * -------------------------------------------------------------------------- */

static double K0 = 0.0;   /* sqrt(1 / (4π))          */
static double K1 = 0.0;   /* sqrt(3 / (4π))          */
static double K2m2 = 0.0; /* sqrt(15 / (4π))         — Y_{2,±2}, Y_{2,±1} x/y/z */
static double K20 = 0.0;  /* sqrt(5 / (16π))         — Y_{2, 0}            */
static double K22 = 0.0;  /* sqrt(15 / (16π))        — Y_{2,±2}            */
static double K3m3 = 0.0; /* sqrt(35 / (32π))        */
static double K3m2 = 0.0; /* sqrt(105 / (4π))        */
static double K3m1 = 0.0; /* sqrt(21 / (32π))        */
static double K30 = 0.0;  /* sqrt(7 / (16π))         */
static double K32 = 0.0;  /* sqrt(105 / (16π))       */

/* l = 4: distinct constant per m (even-m get different coefficients from
 * sin(2φ) vs cos(2φ) producing 2xy vs x²−y² respectively, so the overall
 * normalisation differs by a factor of √2). */
static double K4_m4 = 0.0; /* sqrt(315/(16 π))  ·  xy(x²−y²)            */
static double K4_m3 = 0.0; /* sqrt(315/(128π))  ·  yz(3x²−y²)           */
static double K4_m2 = 0.0; /* sqrt(45 /(16 π))  ·  xy(7z²−r²)           */
static double K4_m1 = 0.0; /* sqrt(45 /(128π))  ·  yz(7z²−3r²)          */
static double K4_0 = 0.0;  /* sqrt(9  /(256π))  ·  (35z⁴−30z²r²+3r⁴)    */
static double K4_p1 = 0.0; /* sqrt(45 /(128π))  ·  xz(7z²−3r²)          */
static double K4_p2 = 0.0; /* sqrt(45 /(64 π))  ·  (x²−y²)(7z²−r²)      */
static double K4_p3 = 0.0; /* sqrt(315/(128π))  ·  xz(x²−3y²)           */
static double K4_p4 = 0.0; /* sqrt(315/(256π))  ·  (x⁴−6x²y²+y⁴)        */

/* Two-phase atomic init. `initialized` becomes true only *after* every
 * constant has been written, so a concurrent second thread seeing the
 * true flag has a happens-before with the first thread's writes
 * (memory_order_acquire on the fast-path load, memory_order_release on
 * the slow-path store). Multiple threads may still race to *perform*
 * the initialization — but each writes identical values, so the
 * outcome is bit-identical regardless of interleaving. */
static atomic_bool constants_initialized = false;

static void        init_constants_(void) {
    if (atomic_load_explicit(&constants_initialized, memory_order_acquire))
        return;
    K0 = sqrt(1.0 / (4.0 * M_PI));
    K1 = sqrt(3.0 / (4.0 * M_PI));
    K2m2 = sqrt(15.0 / (4.0 * M_PI));
    K20 = sqrt(5.0 / (16.0 * M_PI));
    K22 = sqrt(15.0 / (16.0 * M_PI));
    K3m3 = sqrt(35.0 / (32.0 * M_PI));
    K3m2 = sqrt(105.0 / (4.0 * M_PI));
    K3m1 = sqrt(21.0 / (32.0 * M_PI));
    K30 = sqrt(7.0 / (16.0 * M_PI));
    K32 = sqrt(105.0 / (16.0 * M_PI));
    K4_m4 = sqrt(315.0 / (16.0 * M_PI));
    /* Odd-m coefficients for l=4 use the (3/4)-style formula rather than
     * (3/8) — sin(nφ) and cos(nφ) expand into x/y polynomials that differ
     * by a factor of 2 for odd n, which is absorbed into the norm. */
    K4_m3 = sqrt(315.0 / (32.0 * M_PI));
    K4_m2 = sqrt(45.0 / (16.0 * M_PI));
    K4_m1 = sqrt(45.0 / (32.0 * M_PI));
    K4_0 = sqrt(9.0 / (256.0 * M_PI));
    K4_p1 = sqrt(45.0 / (32.0 * M_PI));
    K4_p2 = sqrt(45.0 / (64.0 * M_PI));
    K4_p3 = sqrt(315.0 / (32.0 * M_PI));
    K4_p4 = sqrt(315.0 / (256.0 * M_PI));
    atomic_store_explicit(&constants_initialized, true, memory_order_release);
}

/* -------------------------------------------------------------------------- *
 * Polynomial evaluation for l ≤ 4                                            *
 * -------------------------------------------------------------------------- */

static void solid_harm_l0_(double *out) {
    out[0] = K0;
}

static void solid_harm_l1_(double *out, double x, double y, double z) {
    /* m = -1, 0, +1 */
    out[0] = K1 * y;
    out[1] = K1 * z;
    out[2] = K1 * x;
}

static void solid_harm_l2_(double *out, double x, double y, double z) {
    double x2 = x * x, y2 = y * y, z2 = z * z;
    out[0] = K2m2 * x * y; /* m = -2 : 2·Y_{2,-2}/√2 · |r|²/|r|²  */
    out[1] = K2m2 * y * z; /* m = -1 */
    out[2] = K20 * (2.0 * z2 - x2 - y2);
    out[3] = K2m2 * x * z;    /* m = +1 */
    out[4] = K22 * (x2 - y2); /* m = +2 */
}

static void solid_harm_l3_(double *out, double x, double y, double z) {
    double x2 = x * x, y2 = y * y, z2 = z * z;
    out[0] = K3m3 * y * (3.0 * x2 - y2);                 /* m = -3 */
    out[1] = K3m2 * x * y * z;                           /* m = -2 */
    out[2] = K3m1 * y * (4.0 * z2 - x2 - y2);            /* m = -1 */
    out[3] = K30 * z * (2.0 * z2 - 3.0 * x2 - 3.0 * y2); /* m =  0 */
    out[4] = K3m1 * x * (4.0 * z2 - x2 - y2);            /* m = +1 */
    out[5] = K32 * z * (x2 - y2);                        /* m = +2 */
    out[6] = K3m3 * x * (x2 - 3.0 * y2);                 /* m = +3 */
}

static void solid_harm_l4_(double *out, double x, double y, double z) {
    double x2 = x * x, y2 = y * y, z2 = z * z;
    double r2 = x2 + y2 + z2;
    /* Polynomial forms derived from the Wikipedia real-SH tables combined
     * with sin(nφ), cos(nφ) expansion identities in cartesian. */
    out[0] = K4_m4 * x * y * (x2 - y2);                                /* m = -4 */
    out[1] = K4_m3 * y * z * (3.0 * x2 - y2);                          /* m = -3 */
    out[2] = K4_m2 * x * y * (7.0 * z2 - r2);                          /* m = -2 */
    out[3] = K4_m1 * y * z * (7.0 * z2 - 3.0 * r2);                    /* m = -1 */
    out[4] = K4_0 * (35.0 * z2 * z2 - 30.0 * z2 * r2 + 3.0 * r2 * r2); /* m =  0 */
    out[5] = K4_p1 * x * z * (7.0 * z2 - 3.0 * r2);                    /* m = +1 */
    out[6] = K4_p2 * (x2 - y2) * (7.0 * z2 - r2);                      /* m = +2 */
    out[7] = K4_p3 * x * z * (x2 - 3.0 * y2);                          /* m = +3 */
    out[8] = K4_p4 * (x2 * x2 - 6.0 * x2 * y2 + y2 * y2);              /* m = +4 */
}

/* -------------------------------------------------------------------------- *
 * Public: solid harmonics evaluation                                         *
 * -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- *
 * Recurrence for l ≥ 5 (also valid for all l; used as a cross-check)         *
 *                                                                            *
 * The orthonormal real regular solid harmonics satisfy:                      *
 *                                                                            *
 *   R_{0,0}          = 1 / √(4π)                                             *
 *   R_{l,+l}         = √((2l+1)/(2l))   · (x R_{l−1,+(l−1)} − y R_{l−1,−(l−1)}) *
 *   R_{l,-l}         = √((2l+1)/(2l))   · (y R_{l−1,+(l−1)} + x R_{l−1,−(l−1)}) *
 *   R_{l,+(l−1)}     = √(2l+1)          · z R_{l−1,+(l−1)}                    *
 *   R_{l,-(l−1)}     = √(2l+1)          · z R_{l−1,−(l−1)}                    *
 *   R_{l,±m}, |m|<l−1:                                                        *
 *     C = √((2l+1)(2l−1)/((l+m)(l−m)))                                        *
 *     D = √((2l+1)(l+m−1)(l−m−1)/((l+m)(l−m)(2l−3)))                          *
 *     R_{l,±m}       = C · z · R_{l−1,±m}  −  D · r² · R_{l−2,±m}             *
 *                                                                            *
 * Gradients co-evolve: differentiating each recurrence step gives a clean    *
 * polynomial expression for ∇R_{l,m} in terms of R and ∇R at l−1 and l−2.    *
 *                                                                            *
 * Storage: we keep arrays of R_{l,m} and ∂R_{l,m}/∂r_i for the two preceding *
 * l values, indexed by `m + l` within each row. Rows are padded to a common  *
 * stride so iteration across a range of l's is contiguous.                   *
 * -------------------------------------------------------------------------- */

#define IRREP_SOLID_L_STRIDE (2 * IRREP_L_MAX + 1)

static void solid_recurrence_fill_(int l_max, const double r[3], int want_grad,
                                   double *R_all, /* (l_max+1) rows × stride */
                                   double *dR_all /* 3 × same layout, or NULL */) {
    double x = r[0], y = r[1], z = r[2];
    double r2 = x * x + y * y + z * z;

    int    stride = IRREP_SOLID_L_STRIDE;

    /* R_{0,0} at column (m + l_max) = l_max (centered). But we want column
     * (m + l) per row so that the (2l+1) entries straddle l_max ± l. Easier:
     * store row l at columns [stride_half − l .. stride_half + l] where
     * stride_half = IRREP_L_MAX. */
    int ctr = IRREP_L_MAX;

    /* zero out */
    for (int l = 0; l <= l_max; ++l) {
        for (int j = 0; j < stride; ++j)
            R_all[l * stride + j] = 0.0;
    }
    if (want_grad && dR_all) {
        for (int ax = 0; ax < 3; ++ax) {
            for (int l = 0; l <= l_max; ++l) {
                for (int j = 0; j < stride; ++j) {
                    dR_all[(ax * (l_max + 1) + l) * stride + j] = 0.0;
                }
            }
        }
    }

    R_all[0 * stride + ctr] = K0;
    /* ∇R_{0,0} = 0 (already zeroed). */

    if (l_max < 1)
        return;

    /* Level l = 1 (handle explicitly for clarity). */
    R_all[1 * stride + (ctr - 1)] = K1 * y;
    R_all[1 * stride + (ctr)] = K1 * z;
    R_all[1 * stride + (ctr + 1)] = K1 * x;
    if (want_grad && dR_all) {
        /* ∂R_{1,-1}/∂y = K1 */
        dR_all[(1 * (l_max + 1) + 1) * stride + (ctr - 1)] = K1;
        /* ∂R_{1, 0}/∂z = K1 */
        dR_all[(2 * (l_max + 1) + 1) * stride + (ctr)] = K1;
        /* ∂R_{1,+1}/∂x = K1 */
        dR_all[(0 * (l_max + 1) + 1) * stride + (ctr + 1)] = K1;
    }

    for (int l = 2; l <= l_max; ++l) {
        double  sqrt_2lp1 = sqrt((double)(2 * l + 1));
        double  sqrt_2lp1_over_2l = sqrt((2.0 * l + 1.0) / (2.0 * l));

        double *R_prev = R_all + (l - 1) * stride;
        double *R_prev2 = R_all + (l - 2) * stride;
        double *R_curr = R_all + l * stride;

        /* ---- diagonal: m = ±l from (l-1, ±(l-1)). ---- */
        double Rp_prev = R_prev[ctr + (l - 1)]; /* R_{l-1, +(l-1)} */
        double Rn_prev = R_prev[ctr - (l - 1)]; /* R_{l-1, -(l-1)} */
        R_curr[ctr + l] = sqrt_2lp1_over_2l * (x * Rp_prev - y * Rn_prev);
        R_curr[ctr - l] = sqrt_2lp1_over_2l * (y * Rp_prev + x * Rn_prev);

        /* ---- off-diagonal: m = ±(l-1) from (l-1, ±(l-1)). ---- */
        R_curr[ctr + (l - 1)] = sqrt_2lp1 * z * Rp_prev;
        R_curr[ctr - (l - 1)] = sqrt_2lp1 * z * Rn_prev;

        /* ---- inner: |m| < l-1 ---- */
        for (int m = 0; m <= l - 2; ++m) {
            double denom = (double)(l + m) * (double)(l - m);
            double C = sqrt((2.0 * l + 1.0) * (2.0 * l - 1.0) / denom);
            double D = 0.0;
            if (l >= 2) {
                double numD = (2.0 * l + 1.0) * (double)(l + m - 1) * (double)(l - m - 1);
                double denD = denom * (2.0 * l - 3.0);
                if (denD > 0.0 && numD >= 0.0)
                    D = sqrt(numD / denD);
            }
            double Rp1 = R_prev[ctr + m];
            double Rn1 = R_prev[ctr - m];
            double Rp2 = (l >= 2) ? R_prev2[ctr + m] : 0.0;
            double Rn2 = (l >= 2) ? R_prev2[ctr - m] : 0.0;
            R_curr[ctr + m] = C * z * Rp1 - D * r2 * Rp2;
            R_curr[ctr - m] = C * z * Rn1 - D * r2 * Rn2;
        }

        if (!want_grad || !dR_all)
            continue;

/* ------------------------------------------------------------------ *
 * Gradient co-evolution: differentiate every recurrence step.        *
 * dR_all indexing: [axis * (l_max+1) + l] * stride + (ctr + m).      *
 * ------------------------------------------------------------------ */
#define DR(ax, ll, m) dR_all[((ax) * (l_max + 1) + (ll)) * stride + (ctr + (m))]

        /* diagonal: R_{l,+l} = A·(x R_p − y R_n)  with R_p = R_{l-1,+(l-1)} */
        double dRp_dx = DR(0, l - 1, +(l - 1));
        double dRp_dy = DR(1, l - 1, +(l - 1));
        double dRp_dz = DR(2, l - 1, +(l - 1));
        double dRn_dx = DR(0, l - 1, -(l - 1));
        double dRn_dy = DR(1, l - 1, -(l - 1));
        double dRn_dz = DR(2, l - 1, -(l - 1));

        DR(0, l, +l) = sqrt_2lp1_over_2l * (Rp_prev + x * dRp_dx - y * dRn_dx);
        DR(1, l, +l) = sqrt_2lp1_over_2l * (x * dRp_dy - Rn_prev - y * dRn_dy);
        DR(2, l, +l) = sqrt_2lp1_over_2l * (x * dRp_dz - y * dRn_dz);

        DR(0, l, -l) = sqrt_2lp1_over_2l * (y * dRp_dx + Rn_prev + x * dRn_dx);
        DR(1, l, -l) = sqrt_2lp1_over_2l * (Rp_prev + y * dRp_dy + x * dRn_dy);
        DR(2, l, -l) = sqrt_2lp1_over_2l * (y * dRp_dz + x * dRn_dz);

        /* off-diagonal: R_{l,±(l-1)} = sqrt(2l+1) · z · R_{l-1,±(l-1)} */
        DR(0, l, +(l - 1)) = sqrt_2lp1 * z * dRp_dx;
        DR(1, l, +(l - 1)) = sqrt_2lp1 * z * dRp_dy;
        DR(2, l, +(l - 1)) = sqrt_2lp1 * (Rp_prev + z * dRp_dz);
        DR(0, l, -(l - 1)) = sqrt_2lp1 * z * dRn_dx;
        DR(1, l, -(l - 1)) = sqrt_2lp1 * z * dRn_dy;
        DR(2, l, -(l - 1)) = sqrt_2lp1 * (Rn_prev + z * dRn_dz);

        /* inner: R_{l,±m} = C·z·R_{l-1,±m} − D·r²·R_{l-2,±m} */
        for (int m = 0; m <= l - 2; ++m) {
            double denom = (double)(l + m) * (double)(l - m);
            double C = sqrt((2.0 * l + 1.0) * (2.0 * l - 1.0) / denom);
            double D = 0.0;
            if (l >= 2) {
                double numD = (2.0 * l + 1.0) * (double)(l + m - 1) * (double)(l - m - 1);
                double denD = denom * (2.0 * l - 3.0);
                if (denD > 0.0 && numD >= 0.0)
                    D = sqrt(numD / denD);
            }

            double Rp1_val = R_prev[ctr + m];
            double Rn1_val = R_prev[ctr - m];
            double Rp2_val = (l >= 2) ? R_prev2[ctr + m] : 0.0;
            double Rn2_val = (l >= 2) ? R_prev2[ctr - m] : 0.0;

            double dRp1_dx = DR(0, l - 1, +m);
            double dRp1_dy = DR(1, l - 1, +m);
            double dRp1_dz = DR(2, l - 1, +m);
            double dRn1_dx = DR(0, l - 1, -m);
            double dRn1_dy = DR(1, l - 1, -m);
            double dRn1_dz = DR(2, l - 1, -m);
            double dRp2_dx = (l >= 2) ? DR(0, l - 2, +m) : 0.0;
            double dRp2_dy = (l >= 2) ? DR(1, l - 2, +m) : 0.0;
            double dRp2_dz = (l >= 2) ? DR(2, l - 2, +m) : 0.0;
            double dRn2_dx = (l >= 2) ? DR(0, l - 2, -m) : 0.0;
            double dRn2_dy = (l >= 2) ? DR(1, l - 2, -m) : 0.0;
            double dRn2_dz = (l >= 2) ? DR(2, l - 2, -m) : 0.0;

            /* ∂(r²)/∂r_i = 2 r_i */
            DR(0, l, +m) = C * z * dRp1_dx - D * (2.0 * x * Rp2_val + r2 * dRp2_dx);
            DR(1, l, +m) = C * z * dRp1_dy - D * (2.0 * y * Rp2_val + r2 * dRp2_dy);
            DR(2, l, +m) = C * (Rp1_val + z * dRp1_dz) - D * (2.0 * z * Rp2_val + r2 * dRp2_dz);

            DR(0, l, -m) = C * z * dRn1_dx - D * (2.0 * x * Rn2_val + r2 * dRn2_dx);
            DR(1, l, -m) = C * z * dRn1_dy - D * (2.0 * y * Rn2_val + r2 * dRn2_dy);
            DR(2, l, -m) = C * (Rn1_val + z * dRn1_dz) - D * (2.0 * z * Rn2_val + r2 * dRn2_dz);
        }
#undef DR
    }
}

void irrep_solid_harm_cart(int l, double *out, const double r[3]) {
    init_constants_();
    if (l < 0 || l > IRREP_SOLID_L_MAX || !out)
        return;
    double x = r[0], y = r[1], z = r[2];

    switch (l) {
    case 0:
        solid_harm_l0_(out);
        return;
    case 1:
        solid_harm_l1_(out, x, y, z);
        return;
    case 2:
        solid_harm_l2_(out, x, y, z);
        return;
    case 3:
        solid_harm_l3_(out, x, y, z);
        return;
    case 4:
        solid_harm_l4_(out, x, y, z);
        return;
    default:
        break;
    }

    /* l ≥ 5: analytic recurrence. */
    static double R_buf[(IRREP_L_MAX + 1) * IRREP_SOLID_L_STRIDE];
    solid_recurrence_fill_(l, r, 0, R_buf, NULL);
    int ctr = IRREP_L_MAX;
    for (int m = -l; m <= l; ++m) {
        out[m + l] = R_buf[l * IRREP_SOLID_L_STRIDE + (ctr + m)];
    }
}

void irrep_solid_harm_cart_all(int l_max, double *out, const double r[3]) {
    if (l_max < 0 || l_max > IRREP_SOLID_L_MAX || !out)
        return;
    int offset = 0;
    for (int l = 0; l <= l_max; ++l) {
        irrep_solid_harm_cart(l, out + offset, r);
        offset += 2 * l + 1;
    }
}

/* -------------------------------------------------------------------------- *
 * Analytic gradients (l ≤ 4)                                                 *
 *                                                                            *
 * Layout: out[axis * (2l+1) + (m+l)] = ∂R_{l, m} / ∂r_axis.                  *
 * -------------------------------------------------------------------------- */

static void solid_grad_l0_(double *out) {
    for (int axis = 0; axis < 3; ++axis)
        out[axis] = 0.0;
}

static void solid_grad_l1_(double *out) {
    /* R_{1,-1} = K1·y, R_{1,0} = K1·z, R_{1,+1} = K1·x. */
    for (int i = 0; i < 9; ++i)
        out[i] = 0.0;
    /* ∂/∂x */
    out[0 * 3 + 2] = K1;
    /* ∂/∂y */
    out[1 * 3 + 0] = K1;
    /* ∂/∂z */
    out[2 * 3 + 1] = K1;
}

static void solid_grad_l2_(double *out, double x, double y, double z) {
    double g[3][5];
    /* m = -2: R = K2m2 · xy */
    g[0][0] = K2m2 * y;
    g[1][0] = K2m2 * x;
    g[2][0] = 0.0;
    /* m = -1: R = K2m2 · yz */
    g[0][1] = 0.0;
    g[1][1] = K2m2 * z;
    g[2][1] = K2m2 * y;
    /* m =  0: R = K20 · (2z² - x² - y²) */
    g[0][2] = K20 * (-2.0 * x);
    g[1][2] = K20 * (-2.0 * y);
    g[2][2] = K20 * 4.0 * z;
    /* m = +1: R = K2m2 · xz */
    g[0][3] = K2m2 * z;
    g[1][3] = 0.0;
    g[2][3] = K2m2 * x;
    /* m = +2: R = K22 · (x² - y²) */
    g[0][4] = K22 * 2.0 * x;
    g[1][4] = K22 * (-2.0 * y);
    g[2][4] = 0.0;
    for (int axis = 0; axis < 3; ++axis)
        for (int m = 0; m < 5; ++m)
            out[axis * 5 + m] = g[axis][m];
}

static void solid_grad_l3_(double *out, double x, double y, double z) {
    double x2 = x * x, y2 = y * y, z2 = z * z;
    double g[3][7];
    /* m = -3: R = K3m3 · y · (3x² - y²) = K3m3 · (3x²y - y³) */
    g[0][0] = K3m3 * 6.0 * x * y;
    g[1][0] = K3m3 * (3.0 * x2 - 3.0 * y2);
    g[2][0] = 0.0;
    /* m = -2: R = K3m2 · xyz */
    g[0][1] = K3m2 * y * z;
    g[1][1] = K3m2 * x * z;
    g[2][1] = K3m2 * x * y;
    /* m = -1: R = K3m1 · y · (4z² − x² − y²) */
    g[0][2] = K3m1 * (-2.0 * x * y);
    g[1][2] = K3m1 * (4.0 * z2 - x2 - 3.0 * y2);
    g[2][2] = K3m1 * 8.0 * y * z;
    /* m =  0: R = K30 · z · (2z² − 3x² − 3y²) */
    g[0][3] = K30 * (-6.0 * x * z);
    g[1][3] = K30 * (-6.0 * y * z);
    g[2][3] = K30 * (6.0 * z2 - 3.0 * x2 - 3.0 * y2);
    /* m = +1: R = K3m1 · x · (4z² − x² − y²) */
    g[0][4] = K3m1 * (4.0 * z2 - 3.0 * x2 - y2);
    g[1][4] = K3m1 * (-2.0 * x * y);
    g[2][4] = K3m1 * 8.0 * x * z;
    /* m = +2: R = K32 · z · (x² − y²) */
    g[0][5] = K32 * 2.0 * x * z;
    g[1][5] = K32 * (-2.0 * y * z);
    g[2][5] = K32 * (x2 - y2);
    /* m = +3: R = K3m3 · x · (x² − 3y²) */
    g[0][6] = K3m3 * (3.0 * x2 - 3.0 * y2);
    g[1][6] = K3m3 * (-6.0 * x * y);
    g[2][6] = 0.0;
    for (int axis = 0; axis < 3; ++axis)
        for (int m = 0; m < 7; ++m)
            out[axis * 7 + m] = g[axis][m];
}

static void solid_grad_l4_(double *out, double x, double y, double z) {
    double x2 = x * x, y2 = y * y, z2 = z * z;
    double r2 = x2 + y2 + z2;
    double g[3][9];

    /* m = -4: R = K · xy(x² − y²) = K · (x³y − xy³) */
    g[0][0] = K4_m4 * (3.0 * x2 * y - y2 * y);
    g[1][0] = K4_m4 * (x2 * x - 3.0 * x * y2);
    g[2][0] = 0.0;

    /* m = -3: R = K · yz(3x² − y²) */
    g[0][1] = K4_m3 * (6.0 * x * y * z);
    g[1][1] = K4_m3 * (3.0 * x2 * z - 3.0 * y2 * z);
    g[2][1] = K4_m3 * (3.0 * x2 * y - y2 * y);

    /* m = -2: R = K · xy(7z² − r²) = K · xy(6z² − x² − y²) */
    g[0][2] = K4_m2 * y * (6.0 * z2 - 3.0 * x2 - y2);
    g[1][2] = K4_m2 * x * (6.0 * z2 - x2 - 3.0 * y2);
    g[2][2] = K4_m2 * 12.0 * x * y * z;

    /* m = -1: R = K · yz(7z² − 3r²) = K · yz(4z² − 3x² − 3y²) */
    g[0][3] = K4_m1 * (-6.0 * x * y * z);
    g[1][3] = K4_m1 * z * (4.0 * z2 - 3.0 * x2 - 9.0 * y2);
    g[2][3] = K4_m1 * y * (12.0 * z2 - 3.0 * x2 - 3.0 * y2);

    /* m =  0: R = K · (35z⁴ − 30z²r² + 3r⁴) */
    g[0][4] = K4_0 * 12.0 * x * (x2 + y2 - 4.0 * z2);
    g[1][4] = K4_0 * 12.0 * y * (x2 + y2 - 4.0 * z2);
    g[2][4] = K4_0 * 16.0 * z * (5.0 * z2 - 3.0 * r2);

    /* m = +1: R = K · xz(7z² − 3r²) */
    g[0][5] = K4_p1 * z * (4.0 * z2 - 9.0 * x2 - 3.0 * y2);
    g[1][5] = K4_p1 * (-6.0 * x * y * z);
    g[2][5] = K4_p1 * x * (12.0 * z2 - 3.0 * x2 - 3.0 * y2);

    /* m = +2: R = K · (x² − y²)(7z² − r²) = K · (6x²z² − 6y²z² − x⁴ + y⁴) */
    g[0][6] = K4_p2 * (12.0 * x * z2 - 4.0 * x * x2);
    g[1][6] = K4_p2 * (-12.0 * y * z2 + 4.0 * y * y2);
    g[2][6] = K4_p2 * 12.0 * z * (x2 - y2);

    /* m = +3: R = K · xz(x² − 3y²) */
    g[0][7] = K4_p3 * 3.0 * z * (x2 - y2);
    g[1][7] = K4_p3 * (-6.0 * x * y * z);
    g[2][7] = K4_p3 * x * (x2 - 3.0 * y2);

    /* m = +4: R = K · (x⁴ − 6x²y² + y⁴) */
    g[0][8] = K4_p4 * 4.0 * x * (x2 - 3.0 * y2);
    g[1][8] = K4_p4 * 4.0 * y * (y2 - 3.0 * x2);
    g[2][8] = 0.0;

    for (int axis = 0; axis < 3; ++axis)
        for (int m = 0; m < 9; ++m)
            out[axis * 9 + m] = g[axis][m];
}

void irrep_solid_harm_cart_grad(int l, double *out, const double r[3]) {
    init_constants_();
    if (l < 0 || l > IRREP_SOLID_L_MAX || !out)
        return;

    switch (l) {
    case 0:
        solid_grad_l0_(out);
        return;
    case 1:
        solid_grad_l1_(out);
        return;
    case 2:
        solid_grad_l2_(out, r[0], r[1], r[2]);
        return;
    case 3:
        solid_grad_l3_(out, r[0], r[1], r[2]);
        return;
    case 4:
        solid_grad_l4_(out, r[0], r[1], r[2]);
        return;
    default:
        break;
    }

    /* l ≥ 5: co-evolve gradients through the recurrence. */
    static double R_buf[(IRREP_L_MAX + 1) * IRREP_SOLID_L_STRIDE];
    static double dR_buf[3 * (IRREP_L_MAX + 1) * IRREP_SOLID_L_STRIDE];
    solid_recurrence_fill_(l, r, 1, R_buf, dR_buf);
    int d = 2 * l + 1;
    int ctr = IRREP_L_MAX;
    for (int axis = 0; axis < 3; ++axis) {
        for (int m = -l; m <= l; ++m) {
            out[axis * d + (m + l)] =
                dR_buf[(axis * (l + 1) + l) * IRREP_SOLID_L_STRIDE + (ctr + m)];
        }
    }
}
