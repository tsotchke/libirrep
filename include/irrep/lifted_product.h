/* SPDX-License-Identifier: MIT */
/** @file lifted_product.h
 *  @brief Lifted-product CSS qLDPC codes (Panteleev-Kalachev 2022).
 *
 *  Lifted-product codes generalise the hypergraph product (Tillich-Zemor
 *  2014) by allowing the matrix entries to live in a *non-trivial group
 *  algebra* `F₂[G]` rather than `F₂` directly, with `G` a finite group.
 *  The "lift" operation replaces each entry by a `|G| × |G|` permutation
 *  matrix and gives a quantum code on a much larger alphabet.
 *
 *  This is the family that achieves the *asymptotically good* parameters
 *  `[[n, kn, dn]]` with `k, d → constants > 0` as `n → ∞` — the first
 *  asymptotically good quantum codes (Panteleev-Kalachev 2022).
 *
 *  ## What we implement
 *
 *  We implement the **abelian** lifted product where `G = ℤ/ℓℤ` is the
 *  cyclic group of order `ℓ`. The group algebra `F₂[ℤ/ℓ] ≅ F₂[x]/(xˡ-1)`
 *  consists of polynomials in one variable, mod `xˡ - 1`. Each
 *  polynomial acts as an `ℓ × ℓ` circulant matrix on `F₂^ℓ`.
 *
 *  Given:
 *    - A "base" matrix `A` of size `m_a × n_a` whose entries are
 *      polynomials in `R = F₂[x]/(xˡ-1)`.
 *    - A second base matrix `B` of size `m_b × n_b` over the same `R`.
 *
 *  the lifted product is the CSS code on
 *    `n = ℓ · (n_a · n_b + m_a · m_b)`  qubits,
 *  with check matrices
 *    `H_X = [ A ⊗ I_{n_b}     |  I_{m_a} ⊗ B^T   ]`   (over R)
 *    `H_Z = [ I_{n_a} ⊗ σ(B)  |  σ(A^T) ⊗ I_{m_b} ]`
 *  where `σ` is the antipode of the cyclic group, `σ(x^a) = x^{-a}`,
 *  acting entry-wise on the polynomial matrices. The σ-conjugate is
 *  required for CSS orthogonality at the *F₂-lifted* level: lifting
 *  `R → F₂^{ℓ × ℓ}` maps the F₂ matrix transpose to the σ-action
 *  (`(M_p)^T = M_{σ(p)}`), so a literal block-Kronecker formula must
 *  pre-apply σ in `H_Z` to make `H_X · H_Z^T = 0` at the F₂ level.
 *
 *  ## Specialisations
 *
 *  - `ℓ = 1`: `R = F₂`, recovers the **hypergraph product** of `A` and `B`.
 *  - `m_a = n_a = 1`, `m_b = n_b = 1`: a single polynomial entry on each
 *    side; `n = 2ℓ`, `m_X = m_Z = ℓ`. This is the **abelian
 *    bivariate-bicycle** family at `(ℓ, m=1)`.
 *  - General `m_a, n_a`: the lift gives genuine qLDPC codes with rate
 *    bounded away from 0 as `ℓ → ∞`, and (for carefully chosen `A, B`)
 *    distance scaling linearly in `n`.
 *
 *  ## Primary references
 *
 *  - Panteleev-Kalachev, *Asymptotically good quantum and locally
 *    testable classical LDPC codes*, IEEE Trans. Inf. Theory 68 (2022)
 *    [arXiv:2111.03654].
 *  - Leverrier-Zemor, *Quantum Tanner codes*, FOCS 2022 — alternative
 *    asymptotically-good construction.
 *  - Dinur-Hsieh-Lin-Vidick, *Good locally testable codes*, STOC 2022.
 */
#ifndef IRREP_LIFTED_PRODUCT_H
#define IRREP_LIFTED_PRODUCT_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief A matrix of polynomials in `R = F₂[x]/(xˡ-1)`.
 *
 *  Entries `a[i][j]` are bit-vectors of length `ℓ` (one polynomial per
 *  cell). Storage is row-major in `(rows × cols)` of `n_words = ceil(ℓ/64)`
 *  words each. */
typedef struct {
    int n_rows;     /**< Number of rows. */
    int n_cols;     /**< Number of columns. */
    int ell;        /**< Cyclic group order; modulus xˡ - 1. */
    int n_words;    /**< Words per polynomial = ceil(ℓ/64). */
    uint64_t *data; /**< Flat storage: row-major (n_rows · n_cols) cells. */
} irrep_poly_matrix_t;

/** @brief Allocate a polynomial matrix initialised to all-zero polynomials. */
IRREP_API irrep_status_t
irrep_poly_matrix_new(irrep_poly_matrix_t *out,
                      int n_rows, int n_cols, int ell);

/** @brief Free internal storage. */
IRREP_API void
irrep_poly_matrix_free(irrep_poly_matrix_t *m);

/** @brief Add the monomial `x^a` to entry `(row, col)` (XOR in F₂). */
IRREP_API irrep_status_t
irrep_poly_matrix_add_monomial(irrep_poly_matrix_t *m,
                               int row, int col, int a);

/** @brief Get the coefficient of `x^a` at entry `(row, col)`. */
IRREP_API int
irrep_poly_matrix_get(const irrep_poly_matrix_t *m,
                      int row, int col, int a);

/** @brief Build the lifted-product CSS code from polynomial-matrix
 *  base codes `A` and `B` (both over `F₂[x]/(xˡ-1)` for the same `ℓ`).
 *
 *  Caller must `irrep_css_code_free` when done. CSS orthogonality holds
 *  by construction for any `A, B` over a commutative group algebra. */
IRREP_API irrep_status_t
irrep_lifted_product_build(const irrep_poly_matrix_t *A,
                           const irrep_poly_matrix_t *B,
                           irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_LIFTED_PRODUCT_H */
