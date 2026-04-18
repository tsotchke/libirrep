/* SPDX-License-Identifier: MIT */
/* M2: SU(2) — 2×2 complex unitary group, Pauli matrices, double-cover to SO(3).
 *
 * Convention: U(R(θ, n̂)) = cos(θ/2) I − i sin(θ/2) (n̂ · σ), so for a unit
 * quaternion q = (x, y, z, w),
 *
 *     U = w I − i (x σ_x + y σ_y + z σ_z)
 *       = [[w − i z,  −y − i x],
 *          [ y − i x,   w + i z]].
 */

#include <complex.h>
#include <math.h>

#include <irrep/so3.h>
#include <irrep/su2.h>

#define UNUSED(x) ((void)(x))

static inline irrep_su2_t su2_identity_(void) {
    return (irrep_su2_t){ .u00 = 1.0, .u01 = 0.0, .u10 = 0.0, .u11 = 1.0 };
}

irrep_su2_t irrep_su2_identity(void) { return su2_identity_(); }

irrep_su2_t irrep_su2_from_quat(irrep_quaternion_t q) {
    double n2 = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    if (n2 > 0.0) {
        double inv = 1.0 / sqrt(n2);
        q.x *= inv; q.y *= inv; q.z *= inv; q.w *= inv;
    }
    irrep_su2_t U;
    U.u00 = q.w - q.z * I;
    U.u01 = -q.y - q.x * I;
    U.u10 =  q.y - q.x * I;
    U.u11 = q.w + q.z * I;
    return U;
}

irrep_quaternion_t irrep_quat_from_su2(irrep_su2_t U) {
    double w =  0.5 * (creal(U.u00) + creal(U.u11));
    double z =  0.5 * (cimag(U.u11) - cimag(U.u00));
    double y =  0.5 * (creal(U.u10) - creal(U.u01));
    double x = -0.5 * (cimag(U.u01) + cimag(U.u10));
    return (irrep_quaternion_t){ x, y, z, w };
}

irrep_rot_matrix_t irrep_rot_from_su2(irrep_su2_t U) {
    return irrep_rot_from_quat(irrep_quat_from_su2(U));
}

irrep_su2_t irrep_su2_compose(irrep_su2_t A, irrep_su2_t B) {
    irrep_su2_t C;
    C.u00 = A.u00 * B.u00 + A.u01 * B.u10;
    C.u01 = A.u00 * B.u01 + A.u01 * B.u11;
    C.u10 = A.u10 * B.u00 + A.u11 * B.u10;
    C.u11 = A.u10 * B.u01 + A.u11 * B.u11;
    return C;
}

irrep_su2_t irrep_su2_inverse(irrep_su2_t U) {
    /* For SU(2), U^{-1} = U† (conjugate transpose). */
    return (irrep_su2_t){
        .u00 = conj(U.u00),
        .u01 = conj(U.u10),
        .u10 = conj(U.u01),
        .u11 = conj(U.u11),
    };
}

void irrep_pauli_x(double _Complex out[4]) {
    out[0] = 0.0; out[1] = 1.0;
    out[2] = 1.0; out[3] = 0.0;
}

void irrep_pauli_y(double _Complex out[4]) {
    out[0] = 0.0;     out[1] = -I;
    out[2] =  I;      out[3] = 0.0;
}

void irrep_pauli_z(double _Complex out[4]) {
    out[0] =  1.0; out[1] =  0.0;
    out[2] =  0.0; out[3] = -1.0;
}

void irrep_su2_apply(irrep_su2_t U, const double _Complex in[2], double _Complex out[2]) {
    double _Complex c0 = U.u00 * in[0] + U.u01 * in[1];
    double _Complex c1 = U.u10 * in[0] + U.u11 * in[1];
    out[0] = c0; out[1] = c1;
}

/* -------------------------------------------------------------------------- *
 * 2×2 complex matrix exponential                                             *
 * Scaling-and-squaring with Taylor series; sufficient for SU(2) generators.  *
 * -------------------------------------------------------------------------- */

irrep_su2_t irrep_su2_exp(const double _Complex generator[4]) {
    /* Find max-magnitude entry to decide scaling depth. */
    double norm = 0.0;
    for (int i = 0; i < 4; ++i) {
        double a = cabs(generator[i]);
        if (a > norm) norm = a;
    }
    int scale = 0;
    double scaled = norm;
    while (scaled > 0.5 && scale < 30) { scale++; scaled *= 0.5; }
    double factor = ldexp(1.0, -scale);

    double _Complex A[4] = {
        generator[0] * factor, generator[1] * factor,
        generator[2] * factor, generator[3] * factor,
    };

    /* Taylor: e^A = Σ_{k≥0} A^k / k! */
    double _Complex R[4] = { 1.0, 0.0, 0.0, 1.0 };
    double _Complex P[4] = { 1.0, 0.0, 0.0, 1.0 };
    for (int k = 1; k <= 30; ++k) {
        double _Complex Q0 = P[0] * A[0] + P[1] * A[2];
        double _Complex Q1 = P[0] * A[1] + P[1] * A[3];
        double _Complex Q2 = P[2] * A[0] + P[3] * A[2];
        double _Complex Q3 = P[2] * A[1] + P[3] * A[3];
        double inv_k = 1.0 / (double)k;
        P[0] = Q0 * inv_k; P[1] = Q1 * inv_k;
        P[2] = Q2 * inv_k; P[3] = Q3 * inv_k;
        R[0] += P[0]; R[1] += P[1];
        R[2] += P[2]; R[3] += P[3];
        double pn = cabs(P[0]) + cabs(P[1]) + cabs(P[2]) + cabs(P[3]);
        if (pn < 1e-18) break;
    }

    /* Square `scale` times to undo the rescaling. */
    for (int s = 0; s < scale; ++s) {
        double _Complex R0 = R[0] * R[0] + R[1] * R[2];
        double _Complex R1 = R[0] * R[1] + R[1] * R[3];
        double _Complex R2 = R[2] * R[0] + R[3] * R[2];
        double _Complex R3 = R[2] * R[1] + R[3] * R[3];
        R[0] = R0; R[1] = R1; R[2] = R2; R[3] = R3;
    }
    return (irrep_su2_t){ .u00 = R[0], .u01 = R[1], .u10 = R[2], .u11 = R[3] };
}
