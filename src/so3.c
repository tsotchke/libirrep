/* SPDX-License-Identifier: MIT */
/* M2: SO(3) rotations — all conversion pairs, composition, exp/log, sampling,
 * SLERP, geodesic distance, Fréchet mean (chordal), shortest-arc, validity.
 *
 * Conventions: active rotations, right-handed, Euler ZYZ (physics),
 * quaternion layout {x, y, z, w} with w last. See docs/PHYSICS_APPENDIX.md.
 */

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/so3.h>

#define UNUSED(x) ((void)(x))

/* -------------------------------------------------------------------------- *
 * Local helpers                                                              *
 * -------------------------------------------------------------------------- */

static inline double vec3_norm(const double v[3]) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static inline double quat_norm_sq(irrep_quaternion_t q) {
    return q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
}

/* -------------------------------------------------------------------------- *
 * Quaternion basics                                                          *
 * -------------------------------------------------------------------------- */

irrep_quaternion_t irrep_quat_identity(void) {
    return (irrep_quaternion_t){ .x = 0, .y = 0, .z = 0, .w = 1 };
}

double irrep_quat_norm(irrep_quaternion_t q) {
    return sqrt(quat_norm_sq(q));
}

irrep_quaternion_t irrep_quat_normalize(irrep_quaternion_t q) {
    double n2 = quat_norm_sq(q);
    if (n2 == 0.0) return irrep_quat_identity();
    double inv = 1.0 / sqrt(n2);
    return (irrep_quaternion_t){ q.x * inv, q.y * inv, q.z * inv, q.w * inv };
}

irrep_quaternion_t irrep_quat_conjugate(irrep_quaternion_t q) {
    return (irrep_quaternion_t){ -q.x, -q.y, -q.z, q.w };
}

irrep_quaternion_t irrep_quat_inverse(irrep_quaternion_t q) {
    double n2 = quat_norm_sq(q);
    if (n2 == 0.0) return irrep_quat_identity();
    double inv = 1.0 / n2;
    return (irrep_quaternion_t){ -q.x * inv, -q.y * inv, -q.z * inv, q.w * inv };
}

irrep_quaternion_t irrep_quat_compose(irrep_quaternion_t a, irrep_quaternion_t b) {
    /* Hamilton product with xyzw layout; equivalent to rotation composition
     * "a then b" under the active convention: (ab) v (ab)* = b (a v a*) b*. */
    return (irrep_quaternion_t){
        .x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        .y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        .z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
        .w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
    };
}

/* -------------------------------------------------------------------------- *
 * Conversions                                                                *
 * -------------------------------------------------------------------------- */

irrep_quaternion_t irrep_quat_from_axis_angle(irrep_axis_angle_t aa) {
    double an = vec3_norm(aa.axis);
    if (an < 1e-15) return irrep_quat_identity();
    double half = aa.angle * 0.5;
    double s = sin(half) / an;
    return (irrep_quaternion_t){
        .x = aa.axis[0] * s,
        .y = aa.axis[1] * s,
        .z = aa.axis[2] * s,
        .w = cos(half),
    };
}

irrep_axis_angle_t irrep_axis_angle_from_quat(irrep_quaternion_t q) {
    q = irrep_quat_normalize(q);
    /* Canonicalize w ≥ 0 so angle ∈ [0, π]. */
    if (q.w < 0.0) { q.x = -q.x; q.y = -q.y; q.z = -q.z; q.w = -q.w; }
    double vn = sqrt(q.x * q.x + q.y * q.y + q.z * q.z);
    double angle = 2.0 * atan2(vn, q.w);
    if (vn < 1e-15) {
        return (irrep_axis_angle_t){ .axis = { 0, 0, 1 }, .angle = 0 };
    }
    double inv = 1.0 / vn;
    return (irrep_axis_angle_t){
        .axis  = { q.x * inv, q.y * inv, q.z * inv },
        .angle = angle,
    };
}

irrep_rot_matrix_t irrep_rot_from_quat(irrep_quaternion_t q) {
    q = irrep_quat_normalize(q);
    double xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
    double xy = q.x * q.y, xz = q.x * q.z, yz = q.y * q.z;
    double wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;
    irrep_rot_matrix_t R;
    R.m[0] = 1.0 - 2.0 * (yy + zz);
    R.m[1] = 2.0 * (xy - wz);
    R.m[2] = 2.0 * (xz + wy);
    R.m[3] = 2.0 * (xy + wz);
    R.m[4] = 1.0 - 2.0 * (xx + zz);
    R.m[5] = 2.0 * (yz - wx);
    R.m[6] = 2.0 * (xz - wy);
    R.m[7] = 2.0 * (yz + wx);
    R.m[8] = 1.0 - 2.0 * (xx + yy);
    return R;
}

irrep_quaternion_t irrep_quat_from_rot(irrep_rot_matrix_t R) {
    /* Shepperd 1978 branch switch on the largest diagonal contribution;
     * numerically stable across all rotations. */
    double m00 = R.m[0], m01 = R.m[1], m02 = R.m[2];
    double m10 = R.m[3], m11 = R.m[4], m12 = R.m[5];
    double m20 = R.m[6], m21 = R.m[7], m22 = R.m[8];
    double trace = m00 + m11 + m22;
    irrep_quaternion_t q;
    if (trace > 0.0) {
        double s = sqrt(trace + 1.0) * 2.0; /* s = 4 qw */
        q.w = 0.25 * s;
        q.x = (m21 - m12) / s;
        q.y = (m02 - m20) / s;
        q.z = (m10 - m01) / s;
    } else if (m00 > m11 && m00 > m22) {
        double s = sqrt(1.0 + m00 - m11 - m22) * 2.0; /* s = 4 qx */
        q.w = (m21 - m12) / s;
        q.x = 0.25 * s;
        q.y = (m01 + m10) / s;
        q.z = (m02 + m20) / s;
    } else if (m11 > m22) {
        double s = sqrt(1.0 + m11 - m00 - m22) * 2.0; /* s = 4 qy */
        q.w = (m02 - m20) / s;
        q.x = (m01 + m10) / s;
        q.y = 0.25 * s;
        q.z = (m12 + m21) / s;
    } else {
        double s = sqrt(1.0 + m22 - m00 - m11) * 2.0; /* s = 4 qz */
        q.w = (m10 - m01) / s;
        q.x = (m02 + m20) / s;
        q.y = (m12 + m21) / s;
        q.z = 0.25 * s;
    }
    /* Canonicalize sign for consistent round-trip. */
    if (q.w < 0.0) { q.x = -q.x; q.y = -q.y; q.z = -q.z; q.w = -q.w; }
    return q;
}

irrep_rot_matrix_t irrep_rot_from_euler_zyz(irrep_euler_zyz_t e) {
    double ca = cos(e.alpha), sa = sin(e.alpha);
    double cb = cos(e.beta),  sb = sin(e.beta);
    double cg = cos(e.gamma), sg = sin(e.gamma);
    /* R = Rz(α) Ry(β) Rz(γ). */
    irrep_rot_matrix_t R;
    R.m[0] = ca * cb * cg - sa * sg;
    R.m[1] = -ca * cb * sg - sa * cg;
    R.m[2] = ca * sb;
    R.m[3] = sa * cb * cg + ca * sg;
    R.m[4] = -sa * cb * sg + ca * cg;
    R.m[5] = sa * sb;
    R.m[6] = -sb * cg;
    R.m[7] = sb * sg;
    R.m[8] = cb;
    return R;
}

irrep_euler_zyz_t irrep_euler_zyz_from_rot(irrep_rot_matrix_t R) {
    /* Clamp R[2][2] to [-1, 1] to tolerate tiny float drift. */
    double r22 = R.m[8];
    if (r22 >  1.0) r22 =  1.0;
    if (r22 < -1.0) r22 = -1.0;
    double beta = acos(r22);
    double sb = sin(beta);
    double alpha, gamma;
    if (sb > 1e-12) {
        alpha = atan2(R.m[5],  R.m[2]);     /* (sin α sin β, cos α sin β) */
        gamma = atan2(R.m[7], -R.m[6]);     /* (sin β sin γ, sin β cos γ) */
    } else {
        /* Gimbal lock. */
        gamma = 0.0;
        if (r22 > 0.0) {
            /* β ≈ 0: R = Rz(α + γ); collapse to α. */
            alpha = atan2(R.m[3], R.m[0]);
        } else {
            /* β ≈ π: Ry(π) = diag(-1, 1, -1); R[0][0] = -cos(α-γ). */
            alpha = atan2(-R.m[3], -R.m[0]);
        }
    }
    return (irrep_euler_zyz_t){ .alpha = alpha, .beta = beta, .gamma = gamma };
}

irrep_rot_matrix_t irrep_rot_from_axis_angle(irrep_axis_angle_t aa) {
    return irrep_rot_from_quat(irrep_quat_from_axis_angle(aa));
}

irrep_axis_angle_t irrep_axis_angle_from_rot(irrep_rot_matrix_t R) {
    return irrep_axis_angle_from_quat(irrep_quat_from_rot(R));
}

irrep_quaternion_t irrep_quat_from_euler_zyz(irrep_euler_zyz_t e) {
    return irrep_quat_from_rot(irrep_rot_from_euler_zyz(e));
}

irrep_euler_zyz_t irrep_euler_zyz_from_quat(irrep_quaternion_t q) {
    return irrep_euler_zyz_from_rot(irrep_rot_from_quat(q));
}

/* -------------------------------------------------------------------------- *
 * Group operations                                                           *
 * -------------------------------------------------------------------------- */

irrep_rot_matrix_t irrep_rot_identity(void) {
    return (irrep_rot_matrix_t){ .m = { 1, 0, 0, 0, 1, 0, 0, 0, 1 } };
}

irrep_rot_matrix_t irrep_rot_compose(irrep_rot_matrix_t A, irrep_rot_matrix_t B) {
    irrep_rot_matrix_t C;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += A.m[i * 3 + k] * B.m[k * 3 + j];
            C.m[i * 3 + j] = s;
        }
    }
    return C;
}

irrep_rot_matrix_t irrep_rot_inverse(irrep_rot_matrix_t R) {
    /* SO(3) inverse is the transpose. */
    irrep_rot_matrix_t T;
    T.m[0] = R.m[0]; T.m[1] = R.m[3]; T.m[2] = R.m[6];
    T.m[3] = R.m[1]; T.m[4] = R.m[4]; T.m[5] = R.m[7];
    T.m[6] = R.m[2]; T.m[7] = R.m[5]; T.m[8] = R.m[8];
    return T;
}

/* -------------------------------------------------------------------------- *
 * Exp / log (Rodrigues + quaternion-routed log for π-safety)                 *
 * -------------------------------------------------------------------------- */

irrep_rot_matrix_t irrep_rot_exp(const double omega[3]) {
    double theta = vec3_norm(omega);
    irrep_rot_matrix_t R;
    if (theta < 1e-12) {
        /* R ≈ I + [ω]_× (+ O(θ²)); for our tolerance the linear term suffices. */
        R.m[0] = 1;          R.m[1] = -omega[2];  R.m[2] =  omega[1];
        R.m[3] =  omega[2];  R.m[4] = 1;          R.m[5] = -omega[0];
        R.m[6] = -omega[1];  R.m[7] =  omega[0];  R.m[8] = 1;
        return R;
    }
    double inv = 1.0 / theta;
    double nx = omega[0] * inv, ny = omega[1] * inv, nz = omega[2] * inv;
    double c = cos(theta), s = sin(theta), k = 1.0 - c;
    R.m[0] = c + nx * nx * k;
    R.m[1] = nx * ny * k - nz * s;
    R.m[2] = nx * nz * k + ny * s;
    R.m[3] = ny * nx * k + nz * s;
    R.m[4] = c + ny * ny * k;
    R.m[5] = ny * nz * k - nx * s;
    R.m[6] = nz * nx * k - ny * s;
    R.m[7] = nz * ny * k + nx * s;
    R.m[8] = c + nz * nz * k;
    return R;
}

void irrep_rot_log(irrep_rot_matrix_t R, double omega_out[3]) {
    /* Route through the quaternion: Shepperd is stable across the full
     * rotation range, including angles arbitrarily close to π. This replaces
     * the degenerate acos((tr - 1)/2) formula with a numerically robust path. */
    irrep_quaternion_t   q  = irrep_quat_from_rot(R);
    irrep_axis_angle_t   aa = irrep_axis_angle_from_quat(q);
    omega_out[0] = aa.axis[0] * aa.angle;
    omega_out[1] = aa.axis[1] * aa.angle;
    omega_out[2] = aa.axis[2] * aa.angle;
}

/* -------------------------------------------------------------------------- *
 * Apply rotation to a vector                                                 *
 * -------------------------------------------------------------------------- */

void irrep_rot_apply(irrep_rot_matrix_t R, const double v[3], double out[3]) {
    double x = R.m[0] * v[0] + R.m[1] * v[1] + R.m[2] * v[2];
    double y = R.m[3] * v[0] + R.m[4] * v[1] + R.m[5] * v[2];
    double z = R.m[6] * v[0] + R.m[7] * v[1] + R.m[8] * v[2];
    out[0] = x; out[1] = y; out[2] = z;
}

void irrep_quat_apply(irrep_quaternion_t q, const double v[3], double out[3]) {
    irrep_rot_matrix_t R = irrep_rot_from_quat(q);
    irrep_rot_apply(R, v, out);
}

/* -------------------------------------------------------------------------- *
 * Sampling, interpolation, distance                                          *
 * -------------------------------------------------------------------------- */

static inline uint64_t lcg64_next_(uint64_t *s) {
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return *s;
}

static inline double u01_(uint64_t x) {
    /* top 53 bits → uniform [0, 1) double */
    return (double)(x >> 11) * (1.0 / (double)(1ULL << 53));
}

irrep_quaternion_t irrep_quat_random(uint64_t *rng_state) {
    /* Shoemake 1992, Graphics Gems III, p. 124. */
    if (!rng_state) return irrep_quat_identity();
    double u1 = u01_(lcg64_next_(rng_state));
    double u2 = u01_(lcg64_next_(rng_state));
    double u3 = u01_(lcg64_next_(rng_state));
    const double two_pi = 6.283185307179586476925286766559;
    double s1 = sqrt(1.0 - u1), s2 = sqrt(u1);
    irrep_quaternion_t q = {
        .x = s1 * sin(two_pi * u2),
        .y = s1 * cos(two_pi * u2),
        .z = s2 * sin(two_pi * u3),
        .w = s2 * cos(two_pi * u3),
    };
    if (q.w < 0.0) { q.x = -q.x; q.y = -q.y; q.z = -q.z; q.w = -q.w; }
    return q;
}

irrep_quaternion_t irrep_quat_slerp(irrep_quaternion_t a, irrep_quaternion_t b, double t) {
    double dot = a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    if (dot < 0.0) {
        b.x = -b.x; b.y = -b.y; b.z = -b.z; b.w = -b.w;
        dot = -dot;
    }
    if (dot > 0.9995) {
        /* Linear + normalize for numerical stability near-parallel. */
        irrep_quaternion_t q = {
            .x = a.x * (1.0 - t) + b.x * t,
            .y = a.y * (1.0 - t) + b.y * t,
            .z = a.z * (1.0 - t) + b.z * t,
            .w = a.w * (1.0 - t) + b.w * t,
        };
        return irrep_quat_normalize(q);
    }
    if (dot > 1.0) dot = 1.0;
    double theta = acos(dot);
    double sin_theta = sin(theta);
    double wa = sin((1.0 - t) * theta) / sin_theta;
    double wb = sin(t * theta) / sin_theta;
    return (irrep_quaternion_t){
        .x = a.x * wa + b.x * wb,
        .y = a.y * wa + b.y * wb,
        .z = a.z * wa + b.z * wb,
        .w = a.w * wa + b.w * wb,
    };
}

double irrep_rot_geodesic_distance(irrep_rot_matrix_t A, irrep_rot_matrix_t B) {
    irrep_rot_matrix_t Rd = irrep_rot_compose(irrep_rot_inverse(A), B);
    irrep_axis_angle_t aa = irrep_axis_angle_from_rot(Rd);
    return aa.angle;
}

irrep_quaternion_t irrep_quat_frechet_mean(const irrep_quaternion_t *qs,
                                           const double *weights, size_t n) {
    /* Chordal (extrinsic) mean: sign-align to the first quaternion, then
     * weighted average + normalize. A good approximation of the Karcher mean
     * for quaternions close to each other; v1.1 adds the iterative variant. */
    if (n == 0 || !qs) return irrep_quat_identity();

    irrep_quaternion_t ref = qs[0];
    irrep_quaternion_t acc = { 0, 0, 0, 0 };
    double total = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double w = weights ? weights[i] : 1.0;
        irrep_quaternion_t q = qs[i];
        double d = ref.x * q.x + ref.y * q.y + ref.z * q.z + ref.w * q.w;
        if (d < 0.0) { q.x = -q.x; q.y = -q.y; q.z = -q.z; q.w = -q.w; }
        acc.x += w * q.x;
        acc.y += w * q.y;
        acc.z += w * q.z;
        acc.w += w * q.w;
        total += w;
    }
    if (total == 0.0) return irrep_quat_identity();
    double inv = 1.0 / total;
    acc.x *= inv; acc.y *= inv; acc.z *= inv; acc.w *= inv;
    return irrep_quat_normalize(acc);
}

irrep_quaternion_t irrep_quat_from_two_vectors(const double a[3], const double b[3]) {
    double na = vec3_norm(a), nb = vec3_norm(b);
    if (na < 1e-15 || nb < 1e-15) return irrep_quat_identity();
    double ax = a[0] / na, ay = a[1] / na, az = a[2] / na;
    double bx = b[0] / nb, by = b[1] / nb, bz = b[2] / nb;
    double dot = ax * bx + ay * by + az * bz;
    if (dot > 1.0 - 1e-15) return irrep_quat_identity();
    if (dot < -1.0 + 1e-9) {
        /* 180°: pick any axis perpendicular to a. */
        double axis[3];
        if (fabs(ax) > 0.9) { axis[0] = -ay; axis[1] =  ax; axis[2] = 0.0; }
        else                { axis[0] = 0.0; axis[1] = -az; axis[2] = ay; }
        double an = vec3_norm(axis);
        if (an < 1e-15) return irrep_quat_identity();
        return (irrep_quaternion_t){ axis[0] / an, axis[1] / an, axis[2] / an, 0.0 };
    }
    /* Shortest-arc: axis = a × b, scalar = 1 + a·b, then normalize. */
    double cx = ay * bz - az * by;
    double cy = az * bx - ax * bz;
    double cz = ax * by - ay * bx;
    return irrep_quat_normalize((irrep_quaternion_t){ cx, cy, cz, 1.0 + dot });
}

bool irrep_rot_is_valid(irrep_rot_matrix_t R, double tolerance) {
    /* R^T R = I (orthogonality) */
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += R.m[k * 3 + i] * R.m[k * 3 + j];
            double target = (i == j) ? 1.0 : 0.0;
            if (fabs(s - target) > tolerance) return false;
        }
    }
    /* det(R) = +1 */
    double det = R.m[0] * (R.m[4] * R.m[8] - R.m[5] * R.m[7])
               - R.m[1] * (R.m[3] * R.m[8] - R.m[5] * R.m[6])
               + R.m[2] * (R.m[3] * R.m[7] - R.m[4] * R.m[6]);
    return fabs(det - 1.0) <= tolerance;
}
