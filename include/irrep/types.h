/* SPDX-License-Identifier: MIT */
#ifndef IRREP_TYPES_H
#define IRREP_TYPES_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------- *
 * Rotation element types                                                     *
 * ------------------------------------------------------------------------- */

typedef struct {
    double x, y, z, w;                      /* xyzw layout; unit for SO(3) */
} irrep_quaternion_t;

typedef struct {
    double alpha, beta, gamma;              /* ZYZ, radians */
} irrep_euler_zyz_t;

typedef struct {
    double axis[3];                         /* unit-norm */
    double angle;                           /* radians */
} irrep_axis_angle_t;

typedef struct {
    double m[9];                            /* row-major 3x3 */
} irrep_rot_matrix_t;

/* ------------------------------------------------------------------------- *
 * Irrep labels and multisets                                                 *
 * ------------------------------------------------------------------------- */

#define IRREP_EVEN (+1)
#define IRREP_ODD  (-1)

typedef struct {
    int l;                                  /* non-negative integer */
    int parity;                             /* +1 (even) or -1 (odd) */
} irrep_label_t;

/* Direct sum of irreps, e.g. "1x0e + 2x1o + 1x2e".
 * `labels[i]` and `multiplicities[i]` run 0..num_terms-1.
 * `total_dim = sum_i multiplicities[i] * (2 * labels[i].l + 1)`. */
typedef struct {
    irrep_label_t *labels;
    int           *multiplicities;
    int            num_terms;
    int            capacity;
    int            total_dim;
} irrep_multiset_t;

/* ------------------------------------------------------------------------- *
 * Status / error                                                             *
 * ------------------------------------------------------------------------- */

typedef enum {
    IRREP_OK                  =  0,
    IRREP_ERR_INVALID_ARG     = -1,
    IRREP_ERR_OUT_OF_MEMORY   = -2,
    IRREP_ERR_SELECTION_RULE  = -3,
    IRREP_ERR_NOT_IMPLEMENTED = -4,
    IRREP_ERR_PRECONDITION    = -5,
    IRREP_ERR_PARSE           = -6
} irrep_status_t;

IRREP_API const char *irrep_strerror(irrep_status_t status);

/* Thread-local last-error string, set by internal helpers when a NULL-returning
 * builder fails. Callers may read this immediately after the failing call. */
IRREP_API const char *irrep_last_error(void);

/* ------------------------------------------------------------------------- *
 * Library limits                                                             *
 * ------------------------------------------------------------------------- */

#define IRREP_L_MAX      16
#define IRREP_TWO_J_MAX  32

#ifdef __cplusplus
}
#endif

#endif /* IRREP_TYPES_H */
