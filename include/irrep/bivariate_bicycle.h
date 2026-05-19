/* SPDX-License-Identifier: MIT */
/** @file bivariate_bicycle.h
 *  @brief Bivariate bicycle codes (Bravyi-Cross-Gambetta-Maslov-Rall-Yoder,
 *         *Nature* 627 (2024) 778).
 *
 *  Bivariate bicycle (BB) codes are CSS qLDPC codes parameterised by
 *  a pair of polynomials `A, B ∈ F₂[x, y] / (xˡ - 1, yᵐ - 1)`. The
 *  ring `R = F₂[x,y]/(xˡ-1, yᵐ-1)` has dimension `ℓ · m` over F₂; we
 *  identify `R ≅ F₂^{ℓm}` via the basis `{x^a y^b : 0 ≤ a < ℓ, 0 ≤ b < m}`,
 *  with linear ordering `index(a, b) = a + ℓ · b`.
 *
 *  Multiplication by `x` on this basis is a cyclic shift in the `a`
 *  coordinate; multiplication by `y` is a cyclic shift in the `b`
 *  coordinate. Each polynomial in `R` thus acts as an `(ℓm) × (ℓm)`
 *  binary matrix (a sparse circulant).
 *
 *  The BB code is then built from `n = 2 ℓ m` qubits:
 *
 *    H_X = [ A | B ]       (m_X = ℓm rows)
 *    H_Z = [ B^T | A^T ]   (m_Z = ℓm rows)
 *
 *  CSS orthogonality requires `A B^T + B A^T ≡ 0` in `R`; for `R`
 *  commutative (which it is, since `xy = yx` and `R` is a quotient
 *  of a commutative ring), this holds automatically because
 *  `B^T A = A B^T = (B^T A)^T = A^T B`, so `A B^T + B A^T = A B + B A = 0`
 *  in characteristic 2.
 *
 *  Concrete instances from the Nature 2024 paper:
 *
 *  | code           | ℓ  | m  | n   | k  | d  |
 *  |----------------|----|----|-----|----|----|
 *  | [[72, 12, 6]]  | 6  | 6  | 72  | 12 | 6  |
 *  | [[90, 8, 10]]  | 15 | 3  | 90  | 8  | 10 |
 *  | [[108, 8, 10]] | 9  | 6  | 108 | 8  | 10 |
 *  | [[144, 12, 12]]| 12 | 6  | 144 | 12 | 12 |
 *  | [[288, 12, 18]]| 12 | 12 | 288 | 12 | 18 |
 *
 *  The [[144, 12, 12]] is the headline result that beats the surface code
 *  on IBM hardware at threshold ~0.7%.
 *
 *  ## Polynomial encoding
 *
 *  A polynomial `p(x, y) ∈ R` is stored as the bit-vector of its
 *  coefficients in the canonical basis: bit `a + ℓ·b` is the coefficient
 *  of `x^a y^b`. Construction proceeds by listing the nonzero monomials.
 *
 *  ## Primary references
 *
 *  - Bravyi-Cross-Gambetta-Maslov-Rall-Yoder, *High-threshold and
 *    low-overhead fault-tolerant quantum memory*, Nature 627 (2024) 778
 *    [arXiv:2308.07915].
 *  - MacKay-Mitchison-Shokrollahi, *More sparse-graph codes for quantum
 *    error-correction*, IEEE Trans. Inf. Theory 50 (2004) — sparse
 *    cyclic-code precursor.
 */
#ifndef IRREP_BIVARIATE_BICYCLE_H
#define IRREP_BIVARIATE_BICYCLE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief A polynomial in `F₂[x, y] / (xˡ - 1, yᵐ - 1)`.
 *
 *  The bit at position `a + ℓ · b` is the coefficient of `x^a y^b`. */
typedef struct {
    int ell;        /**< Modulus in x: ring is F₂[x,y]/(xˡ-1, ...). */
    int m;          /**< Modulus in y. */
    int dim;        /**< ℓ · m. */
    int n_words;    /**< ceil(dim / 64). */
    uint64_t *coeffs; /**< Bit-packed coefficient vector. */
} irrep_bb_poly_t;

/** @brief Allocate a zero polynomial in F₂[x,y]/(xˡ-1, yᵐ-1). */
IRREP_API irrep_status_t
irrep_bb_poly_new(irrep_bb_poly_t *out, int ell, int m);

/** @brief Free internal storage. */
IRREP_API void
irrep_bb_poly_free(irrep_bb_poly_t *p);

/** @brief Add the monomial `x^a y^b` to the polynomial (XOR in F₂). */
IRREP_API irrep_status_t
irrep_bb_poly_add_monomial(irrep_bb_poly_t *p, int a, int b);

/** @brief Get the coefficient of `x^a y^b` (0 or 1, or -1 on error). */
IRREP_API int
irrep_bb_poly_get(const irrep_bb_poly_t *p, int a, int b);

/** @brief Build the BB code from polynomials `A` and `B` in the same ring.
 *
 *  Allocates `out` as a CSS code on `n = 2 ℓ m` qubits. The check
 *  matrices are
 *
 *    H_X[r, c] = A[r, c]            for c ∈ [0, ℓm)
 *    H_X[r, c] = B[r, c - ℓm]       for c ∈ [ℓm, 2ℓm)
 *    H_Z[r, c] = B^T[r, c]          for c ∈ [0, ℓm)
 *    H_Z[r, c] = A^T[r, c - ℓm]     for c ∈ [ℓm, 2ℓm)
 *
 *  where the `(ℓm) × (ℓm)` action of A (resp. B) on the basis monomials
 *  `x^a y^b` is computed via the cyclic-shift rule. CSS orthogonality
 *  (`H_X · H_Zᵀ = 0`) is automatic from commutativity of R but is not
 *  verified by this constructor — call `irrep_css_code_verify` to check. */
IRREP_API irrep_status_t
irrep_bb_code_build(const irrep_bb_poly_t *A,
                    const irrep_bb_poly_t *B,
                    irrep_css_code_t *out);

/* ====================================================================
 * Named IBM bivariate-bicycle instances
 *
 * These are the headline qLDPC codes from Bravyi et al. 2024 (Nature
 * 627, 778) that beat the surface code on IBM-Heron-class hardware
 * at physical-error rate ~0.7%. The polynomial-pair specifications
 * follow Table 3 of that paper:
 *
 *   [[72, 12, 6]]   : ℓ =  6, m = 6
 *                     A(x, y) = x³ + y + y²
 *                     B(x, y) = y³ + x + x²
 *
 *   [[144, 12, 12]] : ℓ = 12, m = 6
 *                     A(x, y) = x³ + y + y²
 *                     B(x, y) = y³ + x + x²
 *
 *   [[288, 12, 18]] : ℓ = 12, m = 12
 *                     A(x, y) = x³ + y² + y⁷
 *                     B(x, y) = y³ + x + x²
 *
 * Both have `n = 2 ℓ m` qubits and `k = 12` logical qubits, with
 * distance `d = 12` and `d = 18` respectively. The encoding rate
 * `k/n = 1/12` is the appeal versus the surface code's `1/d²`.
 * ==================================================================== */

/** @brief Build the [[72, 12, 6]] IBM bivariate-bicycle code.
 *
 *  Allocates `out` and fills it with the CSS code defined by
 *  `A = x³ + y + y²` and `B = y³ + x + x²` in `F₂[x, y]/(x⁶ - 1, y⁶ - 1)`. */
IRREP_API irrep_status_t
irrep_bb_code_ibm_72_12_6(irrep_css_code_t *out);

/** @brief Build the [[144, 12, 12]] IBM bivariate-bicycle code.
 *
 *  Allocates `out` and fills it with the CSS code defined by
 *  `A = x³ + y + y²` and `B = y³ + x + x²` in `F₂[x, y]/(x¹² - 1, y⁶ - 1)`.
 *  Caller must `irrep_css_code_free(out)` when done. */
IRREP_API irrep_status_t
irrep_bb_code_ibm_144_12_12(irrep_css_code_t *out);

/** @brief Build the [[288, 12, 18]] IBM bivariate-bicycle code.
 *
 *  Allocates `out` and fills it with the CSS code defined by
 *  `A = x³ + y² + y⁷` and `B = y³ + x + x²` in
 *  `F₂[x, y]/(x¹² - 1, y¹² - 1)`. */
IRREP_API irrep_status_t
irrep_bb_code_ibm_288_12_18(irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_BIVARIATE_BICYCLE_H */
