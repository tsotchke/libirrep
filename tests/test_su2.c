/* SPDX-License-Identifier: MIT */
/* Tests for SU(2): spin-½ matrices, Pauli algebra, the 2-to-1 cover over
 * SO(3), and the matrix exponential.
 *
 * Coverage:
 *   - Identity SU(2) = I.
 *   - SU(2) ↔ unit quaternion round-trip.
 *   - Double cover: ±q map to the same SO(3) matrix.
 *   - Pauli algebra: σ_x σ_x = I  and  σ_x σ_y = i σ_z  (anticommutation).
 *   - Berry phase: D^{1/2}(R(2π, ẑ)) |↑⟩ = −|↑⟩.
 *   - σ_x |↑⟩ = |↓⟩.
 *   - SU(2) compose agrees with quaternion compose.
 *   - SU(2) inverse equals conjugate transpose.
 *   - `su2_exp(−i π σ_z / 2) = −i σ_z`.
 */
#include "harness.h"
#include <irrep/su2.h>
#include <irrep/so3.h>
#include <math.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

static int approx_eq_c(double _Complex a, double _Complex b, double tol) {
    return cabs(a - b) <= tol;
}

static int approx_eq_rot_su2(irrep_rot_matrix_t A, irrep_rot_matrix_t B, double tol) {
    for (int i = 0; i < 9; ++i) if (fabs(A.m[i] - B.m[i]) > tol) return 0;
    return 1;
}

int main(void) {
    IRREP_TEST_START("su2");

    const double tol = 1e-12;

    /* ---- identity ---- */
    {
        irrep_su2_t Umat = irrep_su2_identity();
        IRREP_ASSERT(approx_eq_c(Umat.u00, 1.0, tol));
        IRREP_ASSERT(approx_eq_c(Umat.u11, 1.0, tol));
        IRREP_ASSERT(approx_eq_c(Umat.u01, 0.0, tol));
        IRREP_ASSERT(approx_eq_c(Umat.u10, 0.0, tol));
    }

    /* ---- SU(2) ↔ quaternion round trip ---- */
    {
        irrep_quaternion_t q = irrep_quat_from_axis_angle(
            (irrep_axis_angle_t){ { 0.3, -0.6, 0.2 }, 0.9 });
        irrep_su2_t U = irrep_su2_from_quat(q);
        irrep_quaternion_t q_back = irrep_quat_from_su2(U);
        IRREP_ASSERT_NEAR(q.x, q_back.x, tol);
        IRREP_ASSERT_NEAR(q.y, q_back.y, tol);
        IRREP_ASSERT_NEAR(q.z, q_back.z, tol);
        IRREP_ASSERT_NEAR(q.w, q_back.w, tol);
    }

    /* ---- double cover: q and -q give the same SO(3) matrix ---- */
    {
        irrep_quaternion_t q1 = irrep_quat_from_axis_angle(
            (irrep_axis_angle_t){ { 1, 0, 0 }, 0.7 });
        irrep_quaternion_t q2 = { -q1.x, -q1.y, -q1.z, -q1.w };
        irrep_su2_t U1 = irrep_su2_from_quat(q1);
        irrep_su2_t U2 = irrep_su2_from_quat(q2);
        irrep_rot_matrix_t R1 = irrep_rot_from_su2(U1);
        irrep_rot_matrix_t R2 = irrep_rot_from_su2(U2);
        IRREP_ASSERT(approx_eq_rot_su2(R1, R2, tol));
    }

    /* ---- Pauli algebra: σ_x σ_x = I ---- */
    {
        double _Complex sx[4];
        irrep_pauli_x(sx);
        /* σ_x² = [[0,1],[1,0]]² = [[1,0],[0,1]] */
        double _Complex p00 = sx[0]*sx[0] + sx[1]*sx[2];
        double _Complex p01 = sx[0]*sx[1] + sx[1]*sx[3];
        double _Complex p10 = sx[2]*sx[0] + sx[3]*sx[2];
        double _Complex p11 = sx[2]*sx[1] + sx[3]*sx[3];
        IRREP_ASSERT(approx_eq_c(p00, 1.0, tol));
        IRREP_ASSERT(approx_eq_c(p01, 0.0, tol));
        IRREP_ASSERT(approx_eq_c(p10, 0.0, tol));
        IRREP_ASSERT(approx_eq_c(p11, 1.0, tol));
    }

    /* ---- Pauli anticommutation: σ_x σ_y = i σ_z ---- */
    {
        double _Complex sx[4], sy[4], sz[4];
        irrep_pauli_x(sx);
        irrep_pauli_y(sy);
        irrep_pauli_z(sz);
        double _Complex p00 = sx[0]*sy[0] + sx[1]*sy[2];
        double _Complex p01 = sx[0]*sy[1] + sx[1]*sy[3];
        double _Complex p10 = sx[2]*sy[0] + sx[3]*sy[2];
        double _Complex p11 = sx[2]*sy[1] + sx[3]*sy[3];
        IRREP_ASSERT(approx_eq_c(p00, I * sz[0], tol));
        IRREP_ASSERT(approx_eq_c(p01, I * sz[1], tol));
        IRREP_ASSERT(approx_eq_c(p10, I * sz[2], tol));
        IRREP_ASSERT(approx_eq_c(p11, I * sz[3], tol));
    }

    /* ---- Berry phase: D^{1/2}(R(2π, ẑ)) |↑⟩ = -|↑⟩ ---- */
    {
        irrep_quaternion_t q = irrep_quat_from_axis_angle(
            (irrep_axis_angle_t){ { 0, 0, 1 }, 2.0 * M_PI });
        irrep_su2_t U = irrep_su2_from_quat(q);
        double _Complex up[2]  = { 1.0, 0.0 };
        double _Complex out[2];
        irrep_su2_apply(U, up, out);
        IRREP_ASSERT_NEAR(creal(out[0]), -1.0, 1e-10);
        IRREP_ASSERT_NEAR(cimag(out[0]),  0.0, 1e-10);
        IRREP_ASSERT(cabs(out[1]) < 1e-10);
    }

    /* ---- σ_x on |↑⟩ = |↓⟩ ---- */
    {
        double _Complex sx[4];
        irrep_pauli_x(sx);
        irrep_su2_t Sx = { sx[0], sx[1], sx[2], sx[3] };
        double _Complex up[2]   = { 1.0, 0.0 };
        double _Complex out[2];
        irrep_su2_apply(Sx, up, out);
        IRREP_ASSERT(approx_eq_c(out[0], 0.0, tol));
        IRREP_ASSERT(approx_eq_c(out[1], 1.0, tol));
    }

    /* ---- SU(2) compose vs. quaternion compose ---- */
    {
        irrep_quaternion_t q1 = irrep_quat_from_axis_angle(
            (irrep_axis_angle_t){ { 1, 0, 0 }, 0.3 });
        irrep_quaternion_t q2 = irrep_quat_from_axis_angle(
            (irrep_axis_angle_t){ { 0, 1, 0 }, 0.7 });
        irrep_su2_t U1 = irrep_su2_from_quat(q1);
        irrep_su2_t U2 = irrep_su2_from_quat(q2);
        irrep_su2_t Uc = irrep_su2_compose(U1, U2);
        irrep_quaternion_t qc = irrep_quat_compose(q1, q2);
        irrep_su2_t Uq = irrep_su2_from_quat(qc);
        IRREP_ASSERT(approx_eq_c(Uc.u00, Uq.u00, tol));
        IRREP_ASSERT(approx_eq_c(Uc.u01, Uq.u01, tol));
        IRREP_ASSERT(approx_eq_c(Uc.u10, Uq.u10, tol));
        IRREP_ASSERT(approx_eq_c(Uc.u11, Uq.u11, tol));
    }

    /* ---- SU(2) inverse is conjugate transpose ---- */
    {
        irrep_quaternion_t q = irrep_quat_from_axis_angle(
            (irrep_axis_angle_t){ { 0.1, 0.9, -0.3 }, 1.4 });
        irrep_su2_t U   = irrep_su2_from_quat(q);
        irrep_su2_t Uinv = irrep_su2_inverse(U);
        irrep_su2_t C = irrep_su2_compose(U, Uinv);
        IRREP_ASSERT(approx_eq_c(C.u00, 1.0, tol));
        IRREP_ASSERT(approx_eq_c(C.u11, 1.0, tol));
        IRREP_ASSERT(approx_eq_c(C.u01, 0.0, tol));
        IRREP_ASSERT(approx_eq_c(C.u10, 0.0, tol));
    }

    /* ---- SU(2) exp: exp(-i π σ_z / 2) = -i σ_z ---- */
    {
        double _Complex sz[4];
        irrep_pauli_z(sz);
        double _Complex gen[4];
        for (int i = 0; i < 4; ++i) gen[i] = -I * (M_PI / 2.0) * sz[i];
        irrep_su2_t U = irrep_su2_exp(gen);
        /* exp(-i π σ_z / 2) = cos(π/2) I - i sin(π/2) σ_z = -i σ_z */
        IRREP_ASSERT(approx_eq_c(U.u00, -I * sz[0], 1e-10));
        IRREP_ASSERT(approx_eq_c(U.u11, -I * sz[3], 1e-10));
    }

    return IRREP_TEST_END();
}
