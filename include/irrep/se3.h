/* SPDX-License-Identifier: MIT */
/** @file se3.h
 *  @brief SE(3) — orientation-preserving rigid motions of R³ — and the
 *         Euclidean extension E(3) with parity.
 *
 *  This is the natural completion of `<irrep/so3.h>` to the full Lie
 *  group of 3D pose, with translations. SE(3) underpins:
 *
 *    - 3D registration, ICP, factor graphs (robotics, SLAM).
 *    - SE(3)-equivariant message-passing networks (e3nn / NequIP /
 *      Equiformer family).
 *    - Continuous rigid-body simulation (rigid-body dynamics on
 *      manifolds, screw theory).
 *    - The Lie-algebra side of optical-flow / camera pose estimation.
 *
 *  ## Group structure
 *
 *  An SE(3) element is a pair `g = (R, t)` with `R ∈ SO(3)` and
 *  `t ∈ R³`, acting on a point `p ∈ R³` as `g·p = R p + t`.
 *  Composition is `(R_a, t_a) ∘ (R_b, t_b) = (R_a R_b, R_a t_b + t_a)`;
 *  inverse is `(R, t)^{-1} = (R^T, -R^T t)`.
 *
 *  Internally we store `R` row-major as 9 doubles plus `t` as 3
 *  doubles, total 12 doubles per element. The matrix is *not*
 *  re-projected to SO(3) on each operation — composition stays in
 *  the manifold to machine precision when inputs are; callers who
 *  need long-product stability should periodically re-orthonormalise
 *  via `irrep_se3_renormalise`.
 *
 *  ## Lie algebra se(3)
 *
 *  `se(3)` is 6-dimensional, parameterised by a twist
 *  `ξ = (ω, v) ∈ R³ × R³`. The hat operator gives the 4×4
 *  homogeneous-matrix form
 *
 *      ξ^hat = [ ω×  v ;
 *                0   0 ]
 *
 *  (a 3×3 skew-symmetric block plus a translation column). The
 *  exponential map `exp : se(3) → SE(3)` is the closed-form Rodrigues
 *  formula on the rotation part plus a `V(ω) v` correction on the
 *  translation:
 *
 *      R = I + sin(θ)/θ · [ω]× + (1-cos θ)/θ² · [ω]×²,
 *      t = V(ω) · v,
 *      V(ω) = I + (1-cos θ)/θ² · [ω]× + (θ-sin θ)/θ³ · [ω]×²,
 *
 *  with the standard Taylor expansion in the `θ → 0` limit (handled
 *  to ~1e-15 precision by switching to the series at `θ < 1e-4`).
 *  The logarithm `log : SE(3) → se(3)` is its inverse.
 *
 *  ## E(3) extension
 *
 *  The full Euclidean group E(3) = O(3) ⋉ R³ allows reflections. We
 *  expose a single parity bit `det = ±1`; group operations follow
 *  the semidirect-product rule with det multiplied. This matches
 *  the e3nn "parity" convention and supports E(3)-equivariant
 *  message passing.
 *
 *  ## Primary references
 *
 *  - Murray-Li-Sastry, *A Mathematical Introduction to Robotic
 *    Manipulation*, CRC Press (1994) — the canonical SE(3) reference.
 *  - Sola-Deray-Atchuthan, *A micro Lie theory for state estimation
 *    in robotics*, arXiv:1812.01537 (2018) — modern formulation.
 *  - Geiger et al., *e3nn: Euclidean Neural Networks*, arXiv:2207.09453
 *    (2022) — parity convention and ML usage.
 */
#ifndef IRREP_SE3_H
#define IRREP_SE3_H

#include <stdbool.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================
 * SE(3) — orientation-preserving rigid motions.
 * ==================================================================== */

/** @brief An element of SE(3): rotation matrix + translation vector. */
typedef struct {
    double R[9];   /**< Row-major 3×3 rotation, R · R^T = I, det R = +1. */
    double t[3];   /**< Translation. */
} irrep_se3_t;

/** @brief Identity element (R = I, t = 0). */
IRREP_API irrep_se3_t
irrep_se3_identity(void);

/** @brief Construct from an explicit 3×3 rotation (row-major) and 3-vector. */
IRREP_API irrep_se3_t
irrep_se3_from_R_t(const double R[9], const double t[3]);

/** @brief Construct from axis-angle (unit axis, scalar angle) and translation. */
IRREP_API irrep_se3_t
irrep_se3_from_axis_angle_t(const double axis[3], double angle, const double t[3]);

/** @brief Construct from a unit quaternion (w, x, y, z) and translation. */
IRREP_API irrep_se3_t
irrep_se3_from_quat_t(const double q[4], const double t[3]);

/** @brief Composition `c = a ∘ b`. */
IRREP_API irrep_se3_t
irrep_se3_compose(const irrep_se3_t *a, const irrep_se3_t *b);

/** @brief Group-theoretic inverse. */
IRREP_API irrep_se3_t
irrep_se3_inverse(const irrep_se3_t *g);

/** @brief Apply to a point: `out = R · p + t`. */
IRREP_API void
irrep_se3_apply(const irrep_se3_t *g, const double p[3], double out[3]);

/** @brief Apply only the rotation part to a vector: `out = R · v`. */
IRREP_API void
irrep_se3_rotate(const irrep_se3_t *g, const double v[3], double out[3]);

/** @brief Re-orthonormalise the rotation block via polar / SVD-like
 *  Newton iteration. Useful after many compositions accumulate drift
 *  in the orthonormality of R. Idempotent on already-orthonormal inputs. */
IRREP_API void
irrep_se3_renormalise(irrep_se3_t *g);

/* ====================================================================
 * se(3) Lie algebra and (exp, log) maps.
 * ==================================================================== */

/** @brief Twist in the Lie algebra se(3) ≃ R³ × R³. */
typedef struct {
    double omega[3];  /**< Angular velocity / rotation generator. */
    double v[3];      /**< Linear velocity / translation generator. */
} irrep_se3_twist_t;

/** @brief Exponential map se(3) → SE(3) in closed form.
 *
 *  Uses the Rodrigues formula on the rotation part and the `V(ω) v`
 *  correction on the translation, with Taylor-series fallback at
 *  small `|ω|` for numerical stability. */
IRREP_API irrep_se3_t
irrep_se3_exp(const irrep_se3_twist_t *xi);

/** @brief Logarithm SE(3) → se(3). Inverse of `_exp` for `|ω| < π`. */
IRREP_API irrep_se3_twist_t
irrep_se3_log(const irrep_se3_t *g);

/** @brief Hat operator: write a twist `(ω, v)` as a 4×4 matrix in the
 *  Lie algebra (row-major). Skew-symmetric 3×3 block plus translation
 *  column, with zero bottom row. */
IRREP_API void
irrep_se3_hat(const irrep_se3_twist_t *xi, double mat4x4[16]);

/** @brief Vee operator: inverse of `_hat`. */
IRREP_API void
irrep_se3_vee(const double mat4x4[16], irrep_se3_twist_t *xi);

/** @brief Adjoint representation `Ad_g : se(3) → se(3)` as a 6×6
 *  matrix (row-major). Order of the 6 components is `(ω, v)`.
 *
 *      Ad_g = [ R          0
 *               [t]× R     R ]
 *
 *  Property: `g · exp(ξ) · g^{-1} = exp(Ad_g · ξ)`. */
IRREP_API void
irrep_se3_adjoint(const irrep_se3_t *g, double Ad6x6[36]);

/* ====================================================================
 * E(3) — full Euclidean group with parity.
 * ==================================================================== */

/** @brief An element of E(3): SE(3) element + parity bit (det = ±1). */
typedef struct {
    irrep_se3_t se3;   /**< Underlying orientation-preserving part. */
    int parity;        /**< +1 (proper) or -1 (improper / reflection). */
} irrep_e3_t;

/** @brief Identity (parity +1). */
IRREP_API irrep_e3_t
irrep_e3_identity(void);

/** @brief Reflection through the origin: `(R = -I, t = 0)`, parity -1. */
IRREP_API irrep_e3_t
irrep_e3_inversion(void);

/** @brief E(3) composition: parities multiply, SE(3) parts compose. */
IRREP_API irrep_e3_t
irrep_e3_compose(const irrep_e3_t *a, const irrep_e3_t *b);

/** @brief Inverse. */
IRREP_API irrep_e3_t
irrep_e3_inverse(const irrep_e3_t *g);

/** @brief Apply to a point. For improper elements (parity = -1) the
 *  rotation part `R` already encodes the sign flip via `det R = -1`;
 *  the action is still `out = R · p + t`. */
IRREP_API void
irrep_e3_apply(const irrep_e3_t *g, const double p[3], double out[3]);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SE3_H */
