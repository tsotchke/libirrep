/* SPDX-License-Identifier: MIT */
/** @file su2.h
 *  @brief 2×2 complex unitary group SU(2), Pauli matrices, and the SU(2)
 *         ↔ SO(3) double cover.
 *
 *  SU(2) elements are the universal cover of SO(3): every rotation `R` has
 *  exactly two preimages `±U`. The library uses the isomorphism
 *  `SU(2) ≃ unit quaternions` throughout; @ref irrep_quat_from_su2 and
 *  @ref irrep_su2_from_quat are the canonical bridges.
 */
#ifndef IRREP_SU2_H
#define IRREP_SU2_H

#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief 2×2 complex unitary matrix with `det = 1`, row-major storage. */
typedef struct {
    double _Complex u00, u01, u10, u11;
} irrep_su2_t;

/** @brief Identity SU(2) element. */
IRREP_API irrep_su2_t irrep_su2_identity(void);
/** @brief SU(2) element matching a unit quaternion (canonical double-cover map). */
IRREP_API irrep_su2_t irrep_su2_from_quat(irrep_quaternion_t q);
/** @brief Inverse map — quaternion from an SU(2) element. */
IRREP_API irrep_quaternion_t irrep_quat_from_su2(irrep_su2_t U);
/** @brief Double-cover projection `SU(2) → SO(3)`. Two-to-one. */
IRREP_API irrep_rot_matrix_t irrep_rot_from_su2(irrep_su2_t U);

/** @brief Group composition (matrix product). */
IRREP_API irrep_su2_t irrep_su2_compose(irrep_su2_t a, irrep_su2_t b);
/** @brief Inverse (conjugate transpose, since U is unitary). */
IRREP_API irrep_su2_t irrep_su2_inverse(irrep_su2_t a);
/** @brief Matrix exponential of a complex 2×2 generator (row-major). */
IRREP_API irrep_su2_t irrep_su2_exp(const double _Complex generator[4]);

/** @brief Pauli-X matrix σ_x — `[[0,1],[1,0]]`. */
IRREP_API void irrep_pauli_x(double _Complex out[4]);
/** @brief Pauli-Y matrix σ_y — `[[0,−i],[i,0]]`. */
IRREP_API void irrep_pauli_y(double _Complex out[4]);
/** @brief Pauli-Z matrix σ_z — `[[1,0],[0,−1]]`. */
IRREP_API void irrep_pauli_z(double _Complex out[4]);

/** @brief Apply an SU(2) element to a spin-½ state. */
IRREP_API void irrep_su2_apply(irrep_su2_t U, const double _Complex in[2], double _Complex out[2]);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SU2_H */
