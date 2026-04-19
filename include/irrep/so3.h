/* SPDX-License-Identifier: MIT */
/** @file so3.h
 *  @brief Rotation-group primitives: conversions between quaternions, Euler
 *         ZYZ, axis-angle, and 3×3 rotation matrices; group operations, exp /
 *         log, sampling, SLERP, geodesic distance, and the Karcher mean.
 *
 *  All rotations are **active** (they rotate vectors, not frames) and use a
 *  **right-handed** basis. Euler angles are physics-convention **ZYZ**
 *  (α ∈ [0, 2π), β ∈ [0, π], γ ∈ [0, 2π)). Quaternions are stored as
 *  `{x, y, z, w}` with `w` the scalar component; unit-norm is assumed by the
 *  rotation ops (`irrep_quat_normalize` is provided if you need it).
 */
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

/** @name Conversions
 *  All pairs between quaternion, axis-angle, rotation matrix, and Euler ZYZ.
 *  Every pair round-trips to within `1e-12` on the unit sphere of inputs.
 *  @{ */

/** @brief Quaternion from axis-angle. */
IRREP_API irrep_quaternion_t irrep_quat_from_axis_angle(irrep_axis_angle_t aa);
/** @brief Axis-angle from quaternion (canonical with angle ∈ [0, π]). */
IRREP_API irrep_axis_angle_t irrep_axis_angle_from_quat(irrep_quaternion_t q);
/** @brief 3×3 row-major rotation matrix from quaternion. */
IRREP_API irrep_rot_matrix_t irrep_rot_from_quat(irrep_quaternion_t q);
/** @brief Quaternion from rotation matrix (Shepperd's method). */
IRREP_API irrep_quaternion_t irrep_quat_from_rot(irrep_rot_matrix_t R);
/** @brief Rotation matrix from Euler ZYZ triple. */
IRREP_API irrep_rot_matrix_t irrep_rot_from_euler_zyz(irrep_euler_zyz_t e);
/** @brief Euler ZYZ from rotation matrix (gimbal-lock guarded at |sin β| < 1e-12). */
IRREP_API irrep_euler_zyz_t  irrep_euler_zyz_from_rot(irrep_rot_matrix_t R);
/** @brief Rotation matrix from axis-angle (Rodrigues formula). */
IRREP_API irrep_rot_matrix_t irrep_rot_from_axis_angle(irrep_axis_angle_t aa);
/** @brief Axis-angle from rotation matrix (Markley branch-switching near π). */
IRREP_API irrep_axis_angle_t irrep_axis_angle_from_rot(irrep_rot_matrix_t R);
/** @brief Quaternion from Euler ZYZ triple. */
IRREP_API irrep_quaternion_t irrep_quat_from_euler_zyz(irrep_euler_zyz_t e);
/** @brief Euler ZYZ from quaternion. */
IRREP_API irrep_euler_zyz_t  irrep_euler_zyz_from_quat(irrep_quaternion_t q);
/** @} */

/** @name Group operations
 *  @{ */

/** @brief Identity rotation matrix. */
IRREP_API irrep_rot_matrix_t irrep_rot_identity (void);
/** @brief Identity quaternion {0,0,0,1}. */
IRREP_API irrep_quaternion_t irrep_quat_identity(void);

/** @brief Rotation composition: @p a then @p b (left-multiplication convention). */
IRREP_API irrep_rot_matrix_t irrep_rot_compose (irrep_rot_matrix_t a, irrep_rot_matrix_t b);
/** @brief Rotation inverse (transpose). */
IRREP_API irrep_rot_matrix_t irrep_rot_inverse (irrep_rot_matrix_t a);

/** @brief Quaternion Hamilton product. */
IRREP_API irrep_quaternion_t irrep_quat_compose  (irrep_quaternion_t a, irrep_quaternion_t b);
/** @brief Quaternion inverse (for unit input: conjugate / ‖q‖²). */
IRREP_API irrep_quaternion_t irrep_quat_inverse  (irrep_quaternion_t a);
/** @brief Quaternion conjugate {−x, −y, −z, w}. */
IRREP_API irrep_quaternion_t irrep_quat_conjugate(irrep_quaternion_t a);
/** @brief Project to the unit 3-sphere (idempotent on unit input). */
IRREP_API irrep_quaternion_t irrep_quat_normalize(irrep_quaternion_t a);
/** @brief Euclidean norm of the 4-vector. */
IRREP_API double             irrep_quat_norm     (irrep_quaternion_t a);
/** @} */

/** @name Exp / log on SO(3)
 *  @{ */

/** @brief Matrix exponential of the skew-symmetric generator ω̂ (Rodrigues).
 *  @param omega  3-vector rotation generator in ℝ³; `|ω|` is the angle in radians. */
IRREP_API irrep_rot_matrix_t irrep_rot_exp(const double omega[3]);

/** @brief Matrix logarithm on SO(3) with Markley branch-switching near π.
 *  @param R          valid rotation matrix (orthogonal, det = +1).
 *  @param omega_out  3-vector rotation generator; angle = ‖ω‖, axis = ω / ‖ω‖. */
IRREP_API void               irrep_rot_log(irrep_rot_matrix_t R, double omega_out[3]);
/** @} */

/** @name Vector rotation
 *  @{ */

/** @brief Rotate `v` by `R` into `out` (in-place permitted). */
IRREP_API void irrep_rot_apply   (irrep_rot_matrix_t R, const double v[3], double out[3]);
/** @brief Rotate `v` by unit quaternion `q` into `out` (in-place permitted). */
IRREP_API void irrep_quat_apply  (irrep_quaternion_t q, const double v[3], double out[3]);
/** @} */

/** @name Sampling, interpolation, distance
 *  @{ */

/** @brief Uniform-on-SO(3) quaternion sample (Shoemake 1992).
 *  @param rng_state  caller-owned 64-bit LCG state; updated in place. */
IRREP_API irrep_quaternion_t irrep_quat_random(uint64_t *rng_state);

/** @brief Shortest-arc SLERP from @p a to @p b at parameter `t ∈ [0, 1]`. */
IRREP_API irrep_quaternion_t irrep_quat_slerp (irrep_quaternion_t a,
                                               irrep_quaternion_t b, double t);

/** @brief Geodesic distance on SO(3): the rotation angle of @p a<sup>-1</sup>·@p b. */
IRREP_API double             irrep_rot_geodesic_distance(irrep_rot_matrix_t a,
                                                         irrep_rot_matrix_t b);

/** @brief Weighted Karcher / Fréchet mean on SO(3).
 *
 *  Iteratively solves `μ = argmin_μ Σ wᵢ · d(μ, qᵢ)²` via
 *  `μ_{k+1} = μ_k · exp((1/Σw) · Σ wᵢ · log(μ_k⁻¹ · qᵢ))` with a chordal
 *  initialisation. Converges quadratically near the minimum.
 *  @param qs       array of @p n unit quaternions.
 *  @param weights  array of @p n non-negative weights (may be NULL for uniform).
 *  @param n        number of samples. */
IRREP_API irrep_quaternion_t irrep_quat_frechet_mean(const irrep_quaternion_t *qs,
                                                     const double *weights,
                                                     size_t n);

/** @brief Shortest-arc rotation sending unit vector @p a onto unit vector @p b. */
IRREP_API irrep_quaternion_t irrep_quat_from_two_vectors(const double a[3], const double b[3]);

/** @brief Check that @p R is orthogonal with det = +1 within @p tolerance. */
IRREP_API bool irrep_rot_is_valid(irrep_rot_matrix_t R, double tolerance);
/** @} */

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SO3_H */
