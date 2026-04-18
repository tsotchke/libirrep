/* SPDX-License-Identifier: MIT */
#ifndef IRREP_SO3_H
#define IRREP_SO3_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- conversions (all pairs between quaternion / axis-angle / matrix / euler) ---- */

IRREP_API irrep_quaternion_t irrep_quat_from_axis_angle(irrep_axis_angle_t);
IRREP_API irrep_axis_angle_t irrep_axis_angle_from_quat(irrep_quaternion_t);
IRREP_API irrep_rot_matrix_t irrep_rot_from_quat(irrep_quaternion_t);
IRREP_API irrep_quaternion_t irrep_quat_from_rot(irrep_rot_matrix_t);
IRREP_API irrep_rot_matrix_t irrep_rot_from_euler_zyz(irrep_euler_zyz_t);
IRREP_API irrep_euler_zyz_t  irrep_euler_zyz_from_rot(irrep_rot_matrix_t);
IRREP_API irrep_rot_matrix_t irrep_rot_from_axis_angle(irrep_axis_angle_t);
IRREP_API irrep_axis_angle_t irrep_axis_angle_from_rot(irrep_rot_matrix_t);
IRREP_API irrep_quaternion_t irrep_quat_from_euler_zyz(irrep_euler_zyz_t);
IRREP_API irrep_euler_zyz_t  irrep_euler_zyz_from_quat(irrep_quaternion_t);

/* ---- group operations ---- */

IRREP_API irrep_rot_matrix_t irrep_rot_identity (void);
IRREP_API irrep_quaternion_t irrep_quat_identity(void);

IRREP_API irrep_rot_matrix_t irrep_rot_compose (irrep_rot_matrix_t a, irrep_rot_matrix_t b);
IRREP_API irrep_rot_matrix_t irrep_rot_inverse (irrep_rot_matrix_t a);

IRREP_API irrep_quaternion_t irrep_quat_compose  (irrep_quaternion_t a, irrep_quaternion_t b);
IRREP_API irrep_quaternion_t irrep_quat_inverse  (irrep_quaternion_t a);
IRREP_API irrep_quaternion_t irrep_quat_conjugate(irrep_quaternion_t a);
IRREP_API irrep_quaternion_t irrep_quat_normalize(irrep_quaternion_t a);
IRREP_API double             irrep_quat_norm     (irrep_quaternion_t a);

/* ---- exp / log on rotation matrices (Rodrigues; Markley branch-switching) ---- */

IRREP_API irrep_rot_matrix_t irrep_rot_exp(const double omega[3]);
IRREP_API void               irrep_rot_log(irrep_rot_matrix_t R, double omega_out[3]);

/* ---- vector rotation ---- */

IRREP_API void irrep_rot_apply   (irrep_rot_matrix_t R, const double v[3], double out[3]);
IRREP_API void irrep_quat_apply  (irrep_quaternion_t q, const double v[3], double out[3]);

/* ---- sampling, interpolation, distance ---- */

IRREP_API irrep_quaternion_t irrep_quat_random(uint64_t *rng_state);  /* Shoemake */
IRREP_API irrep_quaternion_t irrep_quat_slerp (irrep_quaternion_t a,
                                               irrep_quaternion_t b, double t);

IRREP_API double             irrep_rot_geodesic_distance(irrep_rot_matrix_t a,
                                                         irrep_rot_matrix_t b);

IRREP_API irrep_quaternion_t irrep_quat_frechet_mean(const irrep_quaternion_t *qs,
                                                     const double *weights,
                                                     size_t n);

/* Shortest-arc rotation sending unit vector `a` onto unit vector `b`. */
IRREP_API irrep_quaternion_t irrep_quat_from_two_vectors(const double a[3], const double b[3]);

/* Validity check: is this an orthogonal matrix with det = +1? */
IRREP_API bool irrep_rot_is_valid(irrep_rot_matrix_t R, double tolerance);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SO3_H */
