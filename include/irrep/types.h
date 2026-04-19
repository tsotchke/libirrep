/* SPDX-License-Identifier: MIT */
/** @file types.h
 *  @brief Scalar types, irrep labels, status codes, and library-wide limits
 *         shared across every other header.
 */
#ifndef IRREP_TYPES_H
#define IRREP_TYPES_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @name Rotation element types
 *  @{ */

/** @brief Quaternion with `xyzw` layout; `w` is the scalar component, unit
 *         norm is assumed by the rotation API (call
 *         @ref irrep_quat_normalize if needed). */
typedef struct {
    double x, y, z, w;
} irrep_quaternion_t;

/** @brief Physics-convention ZYZ Euler triple in radians. Ranges:
 *         `α ∈ [0, 2π)`, `β ∈ [0, π]`, `γ ∈ [0, 2π)`. */
typedef struct {
    double alpha, beta, gamma;
} irrep_euler_zyz_t;

/** @brief Axis-angle rotation: a unit-norm 3-axis and an angle in radians. */
typedef struct {
    double axis[3];
    double angle;
} irrep_axis_angle_t;

/** @brief Row-major 3×3 rotation matrix. */
typedef struct {
    double m[9];
} irrep_rot_matrix_t;
/** @} */

/** @name Irrep labels and multisets
 *  @{ */

#define IRREP_EVEN (+1)     /**< Parity-even O(3) irrep. */
#define IRREP_ODD  (-1)     /**< Parity-odd  O(3) irrep. */

/** @brief Irreducible-representation label for O(3): order `l ≥ 0` and parity. */
typedef struct {
    int l;
    int parity;
} irrep_label_t;

/** @brief Direct sum of irreps, e.g. `"1x0e + 2x1o + 1x2e"`.
 *
 *  Semantic invariant: `total_dim = Σ mult_i · (2 l_i + 1)`. Construct via
 *  @ref irrep_multiset_parse or @ref irrep_multiset_new; feature vectors are
 *  laid out term-major: block i spans `[offset_i, offset_i + mult_i · (2 l_i + 1))`. */
typedef struct {
    irrep_label_t *labels;
    int           *multiplicities;
    int            num_terms;
    int            capacity;
    int            total_dim;
} irrep_multiset_t;
/** @} */

/** @name Status / error
 *  @{ */

/** @brief Status code for validating / builder functions. Pure math functions
 *         return mathematical defaults (`0.0`) on bad input rather than status. */
typedef enum {
    IRREP_OK                  =  0,
    IRREP_ERR_INVALID_ARG     = -1,
    IRREP_ERR_OUT_OF_MEMORY   = -2,
    IRREP_ERR_SELECTION_RULE  = -3,
    IRREP_ERR_NOT_IMPLEMENTED = -4,
    IRREP_ERR_PRECONDITION    = -5,
    IRREP_ERR_PARSE           = -6
} irrep_status_t;

/** @brief Human-readable message for a status code. */
IRREP_API const char *irrep_strerror(irrep_status_t status);

/** @brief Thread-local last-error string, set by internal helpers when a
 *         NULL-returning builder fails. Read immediately after the failing
 *         call; the buffer is overwritten by the next failure on the same thread. */
IRREP_API const char *irrep_last_error(void);
/** @} */

/** @name Library limits
 *  @{ */

#define IRREP_L_MAX      16     /**< Maximum integer angular momentum l supported. */
#define IRREP_TWO_J_MAX  32     /**< Maximum doubled-integer spin supported (j = 16). */
/** @} */

#ifdef __cplusplus
}
#endif

#endif /* IRREP_TYPES_H */
