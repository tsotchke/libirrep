/* SPDX-License-Identifier: MIT */
#ifndef IRREP_SU2_H
#define IRREP_SU2_H

#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 2x2 complex unitary with determinant 1; row-major {u00, u01, u10, u11}. */
typedef struct {
    double _Complex u00, u01, u10, u11;
} irrep_su2_t;

IRREP_API irrep_su2_t   irrep_su2_identity (void);
IRREP_API irrep_su2_t   irrep_su2_from_quat(irrep_quaternion_t);          /* SU(2) ≃ unit quat */
IRREP_API irrep_quaternion_t irrep_quat_from_su2(irrep_su2_t);
IRREP_API irrep_rot_matrix_t irrep_rot_from_su2 (irrep_su2_t);            /* double cover */

IRREP_API irrep_su2_t   irrep_su2_compose(irrep_su2_t a, irrep_su2_t b);
IRREP_API irrep_su2_t   irrep_su2_inverse(irrep_su2_t a);
IRREP_API irrep_su2_t   irrep_su2_exp    (const double _Complex generator[4]);

/* Pauli matrices (row-major 2x2, written into `out`). */
IRREP_API void irrep_pauli_x(double _Complex out[4]);
IRREP_API void irrep_pauli_y(double _Complex out[4]);
IRREP_API void irrep_pauli_z(double _Complex out[4]);

/* Apply an SU(2) element to a spin-1/2 state. */
IRREP_API void irrep_su2_apply(irrep_su2_t U,
                               const double _Complex in[2],
                               double _Complex out[2]);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SU2_H */
