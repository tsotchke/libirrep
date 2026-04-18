/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/so3.h>
#include <math.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

static int approx_eq_rot(irrep_rot_matrix_t A, irrep_rot_matrix_t B, double tol) {
    for (int i = 0; i < 9; ++i) if (fabs(A.m[i] - B.m[i]) > tol) return 0;
    return 1;
}

static int approx_eq_vec3(const double a[3], const double b[3], double tol) {
    return fabs(a[0] - b[0]) <= tol && fabs(a[1] - b[1]) <= tol && fabs(a[2] - b[2]) <= tol;
}

static int approx_eq_quat(irrep_quaternion_t a, irrep_quaternion_t b, double tol) {
    /* Accept q ≡ -q. */
    double dpos = fabs(a.x - b.x) + fabs(a.y - b.y) + fabs(a.z - b.z) + fabs(a.w - b.w);
    double dneg = fabs(a.x + b.x) + fabs(a.y + b.y) + fabs(a.z + b.z) + fabs(a.w + b.w);
    return (dpos <= tol) || (dneg <= tol);
}

int main(void) {
    IRREP_TEST_START("so3");

    const double tol = 1e-12;

    /* ---- identity round trips ---- */
    irrep_rot_matrix_t I  = irrep_rot_identity();
    irrep_quaternion_t qi = irrep_quat_identity();
    IRREP_ASSERT(approx_eq_rot(irrep_rot_from_quat(qi), I, tol));
    IRREP_ASSERT(approx_eq_quat(irrep_quat_from_rot(I), qi, tol));

    /* ---- axis-angle ↔ quaternion ---- */
    {
        irrep_axis_angle_t aa = { .axis = { 0, 0, 1 }, .angle = M_PI / 3.0 };
        irrep_quaternion_t q  = irrep_quat_from_axis_angle(aa);
        irrep_axis_angle_t aa_back = irrep_axis_angle_from_quat(q);
        IRREP_ASSERT_NEAR(aa_back.angle, M_PI / 3.0, tol);
        IRREP_ASSERT_NEAR(aa_back.axis[2], 1.0, tol);
    }

    /* ---- rotation by θ ẑ sends (1,0,0) → (cos θ, sin θ, 0) ---- */
    {
        double theta = 0.7;
        irrep_axis_angle_t aa = { .axis = { 0, 0, 1 }, .angle = theta };
        irrep_rot_matrix_t R  = irrep_rot_from_axis_angle(aa);
        double v[3] = { 1, 0, 0 }, out[3];
        irrep_rot_apply(R, v, out);
        double expected[3] = { cos(theta), sin(theta), 0.0 };
        IRREP_ASSERT(approx_eq_vec3(out, expected, tol));
    }

    /* ---- rotation by θ x̂ sends (0,0,1) → (0, -sin θ, cos θ) ---- */
    {
        double theta = 0.5;
        irrep_axis_angle_t aa = { .axis = { 1, 0, 0 }, .angle = theta };
        irrep_rot_matrix_t R  = irrep_rot_from_axis_angle(aa);
        double v[3] = { 0, 0, 1 }, out[3];
        irrep_rot_apply(R, v, out);
        double expected[3] = { 0.0, -sin(theta), cos(theta) };
        IRREP_ASSERT(approx_eq_vec3(out, expected, tol));
    }

    /* ---- quaternion compose matches matrix compose ---- */
    {
        irrep_quaternion_t q1 = irrep_quat_from_axis_angle((irrep_axis_angle_t){ { 1, 0, 0 }, 0.3 });
        irrep_quaternion_t q2 = irrep_quat_from_axis_angle((irrep_axis_angle_t){ { 0, 1, 0 }, 0.7 });
        irrep_rot_matrix_t R1 = irrep_rot_from_quat(q1);
        irrep_rot_matrix_t R2 = irrep_rot_from_quat(q2);
        irrep_rot_matrix_t Rq = irrep_rot_from_quat(irrep_quat_compose(q1, q2));
        irrep_rot_matrix_t Rm = irrep_rot_compose(R1, R2);
        IRREP_ASSERT(approx_eq_rot(Rq, Rm, tol));
    }

    /* ---- compose(R, R^-1) = I ---- */
    {
        irrep_rot_matrix_t R = irrep_rot_from_axis_angle((irrep_axis_angle_t){ { 0.5, -0.2, 0.8 }, 1.1 });
        IRREP_ASSERT(approx_eq_rot(irrep_rot_compose(R, irrep_rot_inverse(R)), I, tol));
    }

    /* ---- Euler ZYZ round trip (avoiding gimbal lock) ---- */
    {
        irrep_euler_zyz_t e = { .alpha = 0.7, .beta = 1.1, .gamma = 2.3 };
        irrep_rot_matrix_t R = irrep_rot_from_euler_zyz(e);
        irrep_euler_zyz_t back = irrep_euler_zyz_from_rot(R);
        IRREP_ASSERT_NEAR(back.alpha, e.alpha, tol);
        IRREP_ASSERT_NEAR(back.beta,  e.beta,  tol);
        IRREP_ASSERT_NEAR(back.gamma, e.gamma, tol);
    }

    /* ---- rot_exp / rot_log round trip, including near π ---- */
    {
        double omega[3] = { 0.2, -0.4, 0.9 };
        irrep_rot_matrix_t R = irrep_rot_exp(omega);
        double back[3];
        irrep_rot_log(R, back);
        IRREP_ASSERT(approx_eq_vec3(omega, back, tol));
    }
    {
        /* Angle very close to π — quaternion path must be stable. */
        double angle = M_PI - 1e-6;
        double omega[3] = { angle * 0.6, angle * 0.8, 0.0 };
        irrep_rot_matrix_t R = irrep_rot_exp(omega);
        double back[3];
        irrep_rot_log(R, back);
        /* axis.angle is canonicalized to [0, π]; magnitudes should match */
        double mag_in  = sqrt(omega[0]*omega[0] + omega[1]*omega[1] + omega[2]*omega[2]);
        double mag_out = sqrt(back[0]*back[0] + back[1]*back[1] + back[2]*back[2]);
        IRREP_ASSERT_NEAR(mag_in, mag_out, 1e-10);
    }

    /* ---- Shoemake: unit-norm and w ≥ 0 ---- */
    {
        uint64_t st = 0x9E3779B97F4A7C15ULL;
        for (int i = 0; i < 128; ++i) {
            irrep_quaternion_t q = irrep_quat_random(&st);
            double n = irrep_quat_norm(q);
            IRREP_ASSERT(fabs(n - 1.0) <= 1e-14);
            IRREP_ASSERT(q.w >= 0.0);
        }
    }

    /* ---- SLERP endpoints ---- */
    {
        irrep_quaternion_t a = irrep_quat_identity();
        irrep_quaternion_t b = irrep_quat_from_axis_angle((irrep_axis_angle_t){ { 0, 0, 1 }, 0.9 });
        irrep_quaternion_t p0 = irrep_quat_slerp(a, b, 0.0);
        irrep_quaternion_t p1 = irrep_quat_slerp(a, b, 1.0);
        IRREP_ASSERT(approx_eq_quat(p0, a, tol));
        IRREP_ASSERT(approx_eq_quat(p1, b, tol));
    }

    /* ---- geodesic distance: from I to R(θ) equals θ ---- */
    {
        double theta = 1.2345;
        irrep_rot_matrix_t R = irrep_rot_from_axis_angle(
            (irrep_axis_angle_t){ { 0.0, 1.0, 0.0 }, theta });
        double d = irrep_rot_geodesic_distance(I, R);
        IRREP_ASSERT_NEAR(d, theta, tol);
    }

    /* ---- Fréchet mean of identical quaternions is that quaternion ---- */
    {
        irrep_quaternion_t q = irrep_quat_from_axis_angle((irrep_axis_angle_t){ { 1, 0, 0 }, 0.6 });
        irrep_quaternion_t qs[3] = { q, q, q };
        irrep_quaternion_t mean = irrep_quat_frechet_mean(qs, NULL, 3);
        IRREP_ASSERT(approx_eq_quat(mean, q, tol));
    }

    /* ---- shortest-arc from ẑ to x̂ is 90° about ŷ ---- */
    {
        double a[3] = { 0, 0, 1 }, b[3] = { 1, 0, 0 };
        irrep_quaternion_t q = irrep_quat_from_two_vectors(a, b);
        irrep_axis_angle_t aa = irrep_axis_angle_from_quat(q);
        IRREP_ASSERT_NEAR(aa.angle, M_PI / 2.0, tol);
        IRREP_ASSERT_NEAR(fabs(aa.axis[1]), 1.0, tol);
    }

    /* ---- validity check ---- */
    IRREP_ASSERT(irrep_rot_is_valid(I, 1e-14));
    {
        irrep_rot_matrix_t bad = I;
        bad.m[0] = 2.0;
        IRREP_ASSERT(!irrep_rot_is_valid(bad, 1e-10));
    }

    return IRREP_TEST_END();
}
