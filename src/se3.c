/* SPDX-License-Identifier: MIT */
/** @file se3.c
 *  @brief Implementation of SE(3) / E(3) group + se(3) Lie algebra. */
#include <irrep/se3.h>

#include <math.h>
#include <stdbool.h>
#include <string.h>

/* ---------- 3×3 matrix utilities (row-major) ---------- */

static inline void mat3_identity(double R[9]) {
    R[0]=1; R[1]=0; R[2]=0;
    R[3]=0; R[4]=1; R[5]=0;
    R[6]=0; R[7]=0; R[8]=1;
}

static inline void mat3_mul(const double A[9], const double B[9], double C[9]) {
    double t[9];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            t[i*3+j] = A[i*3+0]*B[0*3+j] + A[i*3+1]*B[1*3+j] + A[i*3+2]*B[2*3+j];
        }
    memcpy(C, t, sizeof t);
}

static inline void mat3_apply(const double R[9], const double v[3], double out[3]) {
    double t0 = R[0]*v[0] + R[1]*v[1] + R[2]*v[2];
    double t1 = R[3]*v[0] + R[4]*v[1] + R[5]*v[2];
    double t2 = R[6]*v[0] + R[7]*v[1] + R[8]*v[2];
    out[0] = t0; out[1] = t1; out[2] = t2;
}

static inline void mat3_transpose(const double A[9], double At[9]) {
    double t[9] = { A[0], A[3], A[6],
                    A[1], A[4], A[7],
                    A[2], A[5], A[8] };
    memcpy(At, t, sizeof t);
}

/* Skew-symmetric matrix `[v]×` from a 3-vector. */
static inline void skew3(const double v[3], double S[9]) {
    S[0] =  0;     S[1] = -v[2];  S[2] =  v[1];
    S[3] =  v[2];  S[4] =  0;     S[5] = -v[0];
    S[6] = -v[1];  S[7] =  v[0];  S[8] =  0;
}

/* ---------- SE(3) basic ops ---------- */

irrep_se3_t
irrep_se3_identity(void)
{
    irrep_se3_t g;
    mat3_identity(g.R);
    g.t[0] = g.t[1] = g.t[2] = 0.0;
    return g;
}

irrep_se3_t
irrep_se3_from_R_t(const double R[9], const double t[3])
{
    irrep_se3_t g;
    memcpy(g.R, R, sizeof g.R);
    g.t[0] = t[0]; g.t[1] = t[1]; g.t[2] = t[2];
    return g;
}

/* Rodrigues formula: R = I + sin(θ)·K + (1-cos θ)·K² where K = [axis]×
 * and axis is unit. */
static void
rodrigues_R(const double axis[3], double angle, double R[9])
{
    double c = cos(angle), s = sin(angle);
    double K[9];
    skew3(axis, K);
    double K2[9];
    mat3_mul(K, K, K2);
    for (int i = 0; i < 9; ++i) {
        R[i] = (i % 4 == 0 ? 1.0 : 0.0) + s * K[i] + (1.0 - c) * K2[i];
    }
}

irrep_se3_t
irrep_se3_from_axis_angle_t(const double axis[3], double angle, const double t[3])
{
    irrep_se3_t g;
    /* Normalise axis defensively (callers may pass slightly off-unit). */
    double nrm = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (nrm == 0.0) {
        mat3_identity(g.R);
    } else {
        double u[3] = { axis[0]/nrm, axis[1]/nrm, axis[2]/nrm };
        rodrigues_R(u, angle, g.R);
    }
    g.t[0] = t[0]; g.t[1] = t[1]; g.t[2] = t[2];
    return g;
}

irrep_se3_t
irrep_se3_from_quat_t(const double q[4], const double t[3])
{
    /* Unit quaternion (w, x, y, z) → rotation matrix. */
    double w = q[0], x = q[1], y = q[2], z = q[3];
    irrep_se3_t g;
    g.R[0] = 1 - 2*(y*y + z*z); g.R[1] = 2*(x*y - z*w);   g.R[2] = 2*(x*z + y*w);
    g.R[3] = 2*(x*y + z*w);     g.R[4] = 1 - 2*(x*x+z*z); g.R[5] = 2*(y*z - x*w);
    g.R[6] = 2*(x*z - y*w);     g.R[7] = 2*(y*z + x*w);   g.R[8] = 1 - 2*(x*x+y*y);
    g.t[0] = t[0]; g.t[1] = t[1]; g.t[2] = t[2];
    return g;
}

irrep_se3_t
irrep_se3_compose(const irrep_se3_t *a, const irrep_se3_t *b)
{
    irrep_se3_t c;
    mat3_mul(a->R, b->R, c.R);
    /* t_c = R_a · t_b + t_a. */
    double Rb_t[3];
    mat3_apply(a->R, b->t, Rb_t);
    c.t[0] = Rb_t[0] + a->t[0];
    c.t[1] = Rb_t[1] + a->t[1];
    c.t[2] = Rb_t[2] + a->t[2];
    return c;
}

irrep_se3_t
irrep_se3_inverse(const irrep_se3_t *g)
{
    irrep_se3_t inv;
    mat3_transpose(g->R, inv.R);
    /* t_inv = -R^T · t. */
    double minus_t[3] = { -g->t[0], -g->t[1], -g->t[2] };
    mat3_apply(inv.R, minus_t, inv.t);
    return inv;
}

void
irrep_se3_apply(const irrep_se3_t *g, const double p[3], double out[3])
{
    double Rp[3];
    mat3_apply(g->R, p, Rp);
    out[0] = Rp[0] + g->t[0];
    out[1] = Rp[1] + g->t[1];
    out[2] = Rp[2] + g->t[2];
}

void
irrep_se3_rotate(const irrep_se3_t *g, const double v[3], double out[3])
{
    mat3_apply(g->R, v, out);
}

/* Polar-decomposition iteration: R ← (R + R^{-T}) / 2 with R^{-T}
 * computed via the cofactor / det formula. Converges quadratically
 * to the nearest orthogonal matrix to the input. Five iterations
 * is enough at IEEE double precision for inputs within 1e-2 of
 * SO(3). */
void
irrep_se3_renormalise(irrep_se3_t *g)
{
    if (g == NULL) return;
    double R[9];
    memcpy(R, g->R, sizeof R);
    for (int iter = 0; iter < 8; ++iter) {
        /* R^{-T} via cofactors / det. */
        double cof[9] = {
            R[4]*R[8] - R[5]*R[7],
            R[5]*R[6] - R[3]*R[8],
            R[3]*R[7] - R[4]*R[6],
            R[2]*R[7] - R[1]*R[8],
            R[0]*R[8] - R[2]*R[6],
            R[1]*R[6] - R[0]*R[7],
            R[1]*R[5] - R[2]*R[4],
            R[2]*R[3] - R[0]*R[5],
            R[0]*R[4] - R[1]*R[3],
        };
        double det = R[0]*cof[0] + R[1]*cof[1] + R[2]*cof[2];
        if (det == 0.0) break;
        double inv_det = 1.0 / det;
        double R_inv_T[9];
        for (int i = 0; i < 9; ++i) R_inv_T[i] = cof[i] * inv_det;
        double delta = 0.0;
        for (int i = 0; i < 9; ++i) {
            double nx = 0.5 * (R[i] + R_inv_T[i]);
            double d = nx - R[i];
            R[i] = nx;
            delta += d * d;
        }
        if (delta < 1e-30) break;
    }
    memcpy(g->R, R, sizeof R);
}

/* ---------- Lie algebra: hat, vee, exp, log, adjoint ---------- */

void
irrep_se3_hat(const irrep_se3_twist_t *xi, double mat[16])
{
    /* 4×4 row-major: [ [ω]× | v ; 0 0 0 0 ]. */
    mat[0]  =  0;          mat[1]  = -xi->omega[2]; mat[2]  =  xi->omega[1]; mat[3]  = xi->v[0];
    mat[4]  =  xi->omega[2]; mat[5]  =  0;          mat[6]  = -xi->omega[0]; mat[7]  = xi->v[1];
    mat[8]  = -xi->omega[1]; mat[9]  =  xi->omega[0]; mat[10] =  0;          mat[11] = xi->v[2];
    mat[12] =  0;          mat[13] =  0;          mat[14] =  0;          mat[15] = 0;
}

void
irrep_se3_vee(const double mat[16], irrep_se3_twist_t *xi)
{
    /* ω = (mat[9], mat[2], mat[4]) by the skew convention above. */
    xi->omega[0] = mat[9];
    xi->omega[1] = mat[2];
    xi->omega[2] = mat[4];
    xi->v[0] = mat[3];
    xi->v[1] = mat[7];
    xi->v[2] = mat[11];
}

/* Helper coefficients of the exp / log series in `θ = |ω|`:
 *   A(θ) = sin θ / θ            (∼ 1 - θ²/6 + θ⁴/120)
 *   B(θ) = (1 - cos θ) / θ²      (∼ 1/2 - θ²/24 + θ⁴/720)
 *   C(θ) = (1 - A(θ)) / θ²       (∼ 1/6 - θ²/120 + θ⁴/5040)
 * For small θ we switch to Taylor to keep ~1e-15 accuracy. */
static inline double coef_A(double theta) {
    if (fabs(theta) < 1e-4) {
        double t2 = theta * theta;
        return 1.0 - t2/6.0 + t2*t2/120.0;
    }
    return sin(theta) / theta;
}
static inline double coef_B(double theta) {
    if (fabs(theta) < 1e-4) {
        double t2 = theta * theta;
        return 0.5 - t2/24.0 + t2*t2/720.0;
    }
    return (1.0 - cos(theta)) / (theta * theta);
}
static inline double coef_C(double theta) {
    if (fabs(theta) < 1e-4) {
        double t2 = theta * theta;
        return 1.0/6.0 - t2/120.0 + t2*t2/5040.0;
    }
    return (1.0 - coef_A(theta)) / (theta * theta);
}

irrep_se3_t
irrep_se3_exp(const irrep_se3_twist_t *xi)
{
    irrep_se3_t g;
    double theta = sqrt(xi->omega[0]*xi->omega[0]
                      + xi->omega[1]*xi->omega[1]
                      + xi->omega[2]*xi->omega[2]);

    /* Rotation: R = I + A(θ)·[ω]× + B(θ)·[ω]×². Note that [ω]× is the
     * unnormalised skew of ω, so the angle scaling is implicit in θ. */
    double K[9], K2[9];
    skew3(xi->omega, K);
    mat3_mul(K, K, K2);
    double A = coef_A(theta), B = coef_B(theta);
    for (int i = 0; i < 9; ++i) {
        g.R[i] = (i % 4 == 0 ? 1.0 : 0.0) + A * K[i] + B * K2[i];
    }

    /* Translation: t = V(ω) · v with V(ω) = I + B·[ω]× + C·[ω]×². */
    double C = coef_C(theta);
    double V[9];
    for (int i = 0; i < 9; ++i) {
        V[i] = (i % 4 == 0 ? 1.0 : 0.0) + B * K[i] + C * K2[i];
    }
    mat3_apply(V, xi->v, g.t);
    return g;
}

irrep_se3_twist_t
irrep_se3_log(const irrep_se3_t *g)
{
    irrep_se3_twist_t xi;

    /* Rotation log: θ from trace, axis from skew-part. */
    double tr = g->R[0] + g->R[4] + g->R[8];
    double cos_theta = 0.5 * (tr - 1.0);
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    double theta = acos(cos_theta);

    if (fabs(theta) < 1e-8) {
        /* Near identity: log R ≈ (R - I)/2 in skew. */
        xi.omega[0] = 0.5 * (g->R[7] - g->R[5]);
        xi.omega[1] = 0.5 * (g->R[2] - g->R[6]);
        xi.omega[2] = 0.5 * (g->R[3] - g->R[1]);
    } else {
        double s = theta / (2.0 * sin(theta));
        xi.omega[0] = s * (g->R[7] - g->R[5]);
        xi.omega[1] = s * (g->R[2] - g->R[6]);
        xi.omega[2] = s * (g->R[3] - g->R[1]);
    }

    /* Translation: v = V(ω)^{-1} · t where
     *   V^{-1} = I - (1/2)·[ω]× + ((1 - A/(2B)) / θ²) · [ω]×². */
    double K[9], K2[9];
    skew3(xi.omega, K);
    mat3_mul(K, K, K2);

    double Vinv[9];
    double half = 0.5;
    double coef2;
    if (fabs(theta) < 1e-4) {
        /* Series: ((1 - A/(2B))/θ²) = 1/12 + θ²/720 + ... */
        double t2 = theta * theta;
        coef2 = 1.0/12.0 + t2/720.0 + t2*t2/30240.0;
    } else {
        double A = sin(theta) / theta;
        double B = (1.0 - cos(theta)) / (theta * theta);
        coef2 = (1.0 - A / (2.0 * B)) / (theta * theta);
    }
    for (int i = 0; i < 9; ++i) {
        Vinv[i] = (i % 4 == 0 ? 1.0 : 0.0) - half * K[i] + coef2 * K2[i];
    }
    mat3_apply(Vinv, g->t, xi.v);
    return xi;
}

void
irrep_se3_adjoint(const irrep_se3_t *g, double Ad6x6[36])
{
    /* Ad_g = [[ R, 0 ], [ [t]× R, R ]]. Row-major 6×6. Order (ω₁ω₂ω₃ v₁v₂v₃). */
    double tx[9];
    skew3(g->t, tx);
    double txR[9];
    mat3_mul(tx, g->R, txR);

    /* Upper-left R, upper-right zero. */
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Ad6x6[i*6 + j]     = g->R[i*3 + j];   /* upper-left */
            Ad6x6[i*6 + (j+3)] = 0.0;             /* upper-right */
        }
    }
    /* Lower-left [t]×R, lower-right R. */
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Ad6x6[(i+3)*6 + j]     = txR[i*3 + j];  /* lower-left */
            Ad6x6[(i+3)*6 + (j+3)] = g->R[i*3 + j]; /* lower-right */
        }
    }
}

/* ---------- E(3) ---------- */

irrep_e3_t
irrep_e3_identity(void)
{
    irrep_e3_t g;
    g.se3 = irrep_se3_identity();
    g.parity = +1;
    return g;
}

irrep_e3_t
irrep_e3_inversion(void)
{
    irrep_e3_t g;
    g.se3 = irrep_se3_identity();
    /* Negate R so that det = -1. */
    for (int i = 0; i < 9; ++i) g.se3.R[i] = -g.se3.R[i];
    g.parity = -1;
    return g;
}

irrep_e3_t
irrep_e3_compose(const irrep_e3_t *a, const irrep_e3_t *b)
{
    irrep_e3_t c;
    c.se3 = irrep_se3_compose(&a->se3, &b->se3);
    c.parity = a->parity * b->parity;
    return c;
}

irrep_e3_t
irrep_e3_inverse(const irrep_e3_t *g)
{
    irrep_e3_t inv;
    inv.se3 = irrep_se3_inverse(&g->se3);
    inv.parity = g->parity;   /* parity is self-inverse in Z/2Z. */
    return inv;
}

void
irrep_e3_apply(const irrep_e3_t *g, const double p[3], double out[3])
{
    irrep_se3_apply(&g->se3, p, out);
}
