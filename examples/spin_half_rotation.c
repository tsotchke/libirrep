/* SPDX-License-Identifier: MIT */
/* Demonstrates the SU(2) Berry phase of a 2π rotation on a spin-½ state.
 *
 * The state |↑⟩ sent through D^{1/2}(R(2π, ẑ)) picks up a minus sign —
 * a famously measurable consequence of SU(2) double-covering SO(3). */

#include <complex.h>
#include <math.h>
#include <stdio.h>

#include <irrep/so3.h>
#include <irrep/su2.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

int main(void) {
    /* Build U = U(q) where q = (0, 0, sin(θ/2), cos(θ/2)), θ = 2π. */
    irrep_quaternion_t q = irrep_quat_from_axis_angle(
        (irrep_axis_angle_t){ .axis = { 0, 0, 1 }, .angle = 2.0 * M_PI });
    irrep_su2_t U = irrep_su2_from_quat(q);

    double _Complex up[2]  = { 1.0, 0.0 };     /* |↑⟩ */
    double _Complex out[2];
    irrep_su2_apply(U, up, out);

    printf("spin_half_rotation\n");
    printf("  D^{1/2}(2π, ẑ) |↑⟩ = (% .15f %+.15fi,  % .15f %+.15fi)\n",
           (double)__real__ out[0], (double)__imag__ out[0],
           (double)__real__ out[1], (double)__imag__ out[1]);
    printf("  expected (Berry phase -1) ≈ (-1, 0).\n");

    int ok = fabs(creal(out[0]) + 1.0) < 1e-10
          && fabs(cimag(out[0]))       < 1e-10
          && cabs(out[1])              < 1e-10;
    printf("  %s\n", ok ? "OK" : "MISMATCH");
    return ok ? 0 : 1;
}
