/* SPDX-License-Identifier: MIT */
/** @file css_code.h
 *  @brief Calderbank–Shor–Steane (CSS) quantum code from a pair of
 *         binary parity-check matrices.
 *
 *  A CSS code is a stabilizer code in which every stabilizer is *either*
 *  a tensor product of `X` and `I` operators, *or* a tensor product of
 *  `Z` and `I` operators. It is fully specified by two binary matrices:
 *
 *    - `H_X ∈ F₂^{m_X × n}` — each row is the support of an X-stabilizer.
 *    - `H_Z ∈ F₂^{m_Z × n}` — each row is the support of a Z-stabilizer.
 *
 *  These must satisfy the **CSS orthogonality condition**:
 *
 *    `H_X · H_Zᵀ ≡ 0 (mod 2)`
 *
 *  i.e. every row of `H_X` is orthogonal in the F₂ inner product to every
 *  row of `H_Z`. That is exactly the stabilizer commutativity condition
 *  specialised to CSS codes (the symplectic inner product reduces to
 *  `Σᵢ x_i z'_i + z_i x'_i = 0 + Σᵢ z_i x'_i` for the X-Z case).
 *
 *  ## Why this matters for libirrep
 *
 *  Almost every modern stabilizer code in the QEC literature is CSS:
 *
 *    - **Toric / surface code** — `H_X` from vertex stabilizers, `H_Z`
 *      from plaquette stabilizers.
 *    - **Color code** [Bombín-Martín-Delgado 2006] — `H_X = H_Z` from
 *      face stabilizers, with the special property that both X- and
 *      Z-checks on every face are simultaneously commuting.
 *    - **Bivariate bicycle codes** [Bravyi et al. Nature 2024] —
 *      `H_X = [A | B]`, `H_Z = [Bᵀ | Aᵀ]` for `A, B ∈ F₂[x,y]/(xˡ-1, yᵐ-1)`.
 *    - **Hypergraph product** [Tillich-Zemor 2014] — `H_X = (h_a ⊗ I |
 *      I ⊗ h_bᵀ)`, `H_Z = (I ⊗ h_b | h_aᵀ ⊗ I)` from classical codes.
 *    - **Lifted product** [Panteleev-Kalachev 2022] — generalises hypergraph
 *      product to non-abelian group lifts.
 *
 *  All of them collapse onto the same `(H_X, H_Z)` interface; the
 *  difference is just how the matrices are constructed.
 *
 *  ## Storage
 *
 *  Parity-check matrices are stored row-major, bit-packed. Each row uses
 *  `ceil(n / 64)` words. This matches the layout of `irrep_pauli_t.x`
 *  and `irrep_pauli_t.z`, so converting a CSS code to an
 *  `irrep_stabilizer_group_t` is just a memcpy per generator.
 *
 *  ## Primary references
 *
 *  - Calderbank-Shor, *Phys. Rev. A* 54 (1996) 1098.
 *  - Steane, *Phys. Rev. Lett.* 77 (1996) 793.
 *  - MacKay-Mitchison-McFadden, *IEEE Trans. Inf. Theory* 50 (2004) —
 *    sparse-quantum-code construction.
 */
#ifndef IRREP_CSS_CODE_H
#define IRREP_CSS_CODE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/stabilizer_group.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Bit-packed parity-check matrix.
 *
 *  Rows are F₂ vectors of length `n_cols`, packed into `n_words = ceil(n_cols/64)`
 *  words per row. `data` has length `n_rows * n_words`. */
typedef struct {
    int n_rows;       /**< Number of stabilizers / parity checks. */
    int n_cols;       /**< Number of qubits / classical bits. */
    int n_words;      /**< Words per row = ceil(n_cols / 64). */
    uint64_t *data;   /**< Row-major bit-packed storage. */
} irrep_parity_matrix_t;

/** @brief Allocate a parity-check matrix initialised to all-zero. */
IRREP_API irrep_status_t
irrep_parity_matrix_new(irrep_parity_matrix_t *out, int n_rows, int n_cols);

/** @brief Free internal storage. (Does not free the struct itself.) */
IRREP_API void
irrep_parity_matrix_free(irrep_parity_matrix_t *m);

/** @brief Set bit (row, col) to 1. */
IRREP_API irrep_status_t
irrep_parity_matrix_set(irrep_parity_matrix_t *m, int row, int col);

/** @brief Get bit (row, col). Returns 0 or 1, or -1 on error. */
IRREP_API int
irrep_parity_matrix_get(const irrep_parity_matrix_t *m, int row, int col);

/** @brief F₂ inner product of two rows from possibly-different matrices.
 *
 *  Computes `Σⱼ a[r_a, j] · b[r_b, j] (mod 2)`. Returns 0 or 1, or -1 on error. */
IRREP_API int
irrep_parity_matrix_row_inner(const irrep_parity_matrix_t *a, int row_a,
                              const irrep_parity_matrix_t *b, int row_b);

/* ====================================================================
 * CSS code
 * ==================================================================== */

/** @brief CSS code defined by an X-check matrix and a Z-check matrix.
 *
 *  Both matrices must have the same `n_cols`, which equals the number
 *  of physical qubits. The number of rows of `H_X` (resp. `H_Z`) is the
 *  number of X- (resp. Z-) stabilizer generators. The total number of
 *  stabilizer generators is `m_X + m_Z`; the number of logical qubits
 *  is `k = n - rank(H_X) - rank(H_Z)` (caller-computed; not stored). */
typedef struct {
    int n;                         /**< Number of physical qubits. */
    irrep_parity_matrix_t H_X;     /**< X-stabilizer support matrix. */
    irrep_parity_matrix_t H_Z;     /**< Z-stabilizer support matrix. */
} irrep_css_code_t;

/** @brief Allocate a CSS code with `n` qubits and the given check counts. */
IRREP_API irrep_status_t
irrep_css_code_new(irrep_css_code_t *out, int n, int m_X, int m_Z);

/** @brief Free internal storage. */
IRREP_API void
irrep_css_code_free(irrep_css_code_t *c);

/** @brief Verify the CSS orthogonality condition `H_X · H_Zᵀ ≡ 0 (mod 2)`.
 *
 *  Returns IRREP_OK if every X-stabilizer commutes with every Z-stabilizer
 *  (the only commutativity check needed in a CSS code, since pure-X
 *  stabilizers always commute with each other and pure-Z always commute
 *  with each other). Returns IRREP_ERR_PRECONDITION if the condition fails. */
IRREP_API irrep_status_t
irrep_css_code_verify(const irrep_css_code_t *c);

/** @brief Materialise the CSS code as an abstract stabilizer group.
 *
 *  Generators are ordered: X-stabilizers (rows of `H_X`) at indices
 *  `[0, m_X)`, then Z-stabilizers (rows of `H_Z`) at indices
 *  `[m_X, m_X + m_Z)`. Caller must `irrep_stabilizer_group_free` when done. */
IRREP_API irrep_status_t
irrep_css_code_to_stabilizer_group(const irrep_css_code_t *c,
                                   irrep_stabilizer_group_t *out);

/** @brief Compute the F₂-rank of a parity-check matrix via Gaussian
 *  elimination.
 *
 *  Runs in-place on a workspace copy of `H` so the input is unchanged.
 *  Complexity: `O(m · n · m / 64)` with bit-packed XOR. Use for
 *  rank-derived diagnostics (logical-qubit counts, code dimension,
 *  syndrome redundancies).
 *
 *  @param[in] H  Parity-check matrix; can have any dimensions, including
 *                degenerate (0-row or 0-col).
 *  @return  Rank of `H` over F₂, in `[0, min(m, n)]`. Returns -1 on
 *           allocation failure. */
IRREP_API int
irrep_parity_matrix_rank(const irrep_parity_matrix_t *H);

/** @brief Number of logical qubits of a CSS code.
 *
 *  Computed as `k = n - rank(H_X) - rank(H_Z)`. CSS-orthogonality
 *  guarantees `rank(H_X) + rank(H_Z) ≤ n`, so `k ≥ 0`. For
 *  well-formed codes this returns the analytical logical-qubit count
 *  without any code-family-specific arithmetic.
 *
 *  @param[in] c  CSS code (read-only).
 *  @return  Number of logical qubits, or -1 on error. */
IRREP_API int
irrep_css_code_logical_qubits(const irrep_css_code_t *c);

/** @brief Brute-force code distance via stabilizer-group enumeration.
 *
 *  Convenience wrapper that materialises the CSS code as a stabilizer
 *  group, calls `irrep_qec_distance_brute(g, max_weight)`, and frees
 *  the working storage. Useful for one-shot code-parameter audits.
 *
 *  Complexity: O(C(n, w) · 3^w · m · n) per weight class up to
 *  `max_weight`. Practical for `n ≤ 25` at distances `≤ 5`.
 *
 *  @param[in] c           CSS code (read-only).
 *  @param[in] max_weight  Largest Pauli weight to enumerate.
 *  @return  Smallest weight `w ≤ max_weight` for which a non-trivial
 *           logical operator exists; `max_weight + 1` if none found
 *           in range; -1 on error. */
IRREP_API int
irrep_css_code_distance(const irrep_css_code_t *c, int max_weight);

/** @brief Find a minimum-weight logical X̄ for any CSS code by brute-
 *  force search over X-only Paulis of weight up to `max_weight`.
 *
 *  Iterates over all weight-w supports for w = 1..max_weight and
 *  returns the FIRST X^v Pauli such that:
 *    1. `v` commutes with every Z-stabilizer (v · z = 0 for every z
 *       in the rowspace of H_Z), AND
 *    2. `v` is NOT in the rowspace of H_X (not a trivial stabilizer).
 *
 *  Such a `v` is by definition a non-trivial logical X̄ operator;
 *  the first one found has minimum weight = `d_X` (the X-side
 *  distance). The overall CSS code distance is `d = min(d_X, d_Z)`,
 *  so `d_X` may exceed `d` when the Z-side gives smaller logicals
 *  (e.g., the [[15, 1, 3]] Reed-Muller code has d_X = 7 but
 *  d_Z = 3, so its distance is 3 although this X-side search
 *  returns 7).
 *
 *  Complexity: O(C(n, w) · m · n) per weight class. Practical for
 *  `n ≤ 30` at small d_X.
 *
 *  @param[in]  c           CSS code (read-only).
 *  @param[in]  max_weight  Largest weight to search.
 *  @param[out] out_pauli   Caller-allocated `irrep_pauli_t`; on
 *                          success holds the discovered logical X̄.
 *  @return  Weight of the discovered logical (= d_X) on success;
 *           `max_weight + 1` if none found in range; -1 on error. */
IRREP_API int
irrep_css_code_compute_logical_X(const irrep_css_code_t *c, int max_weight,
                                  irrep_pauli_t *out_pauli);

/** @brief Find a minimum-weight logical Z̄ for any CSS code by brute-
 *  force search over Z-only Paulis (the Z-side companion to
 *  `irrep_css_code_compute_logical_X`).
 *
 *  Iterates over weight-w supports and returns the FIRST Z^v such
 *  that v commutes with all X-stabilizers (`v · x = 0` for every x
 *  in the rowspace of H_X) and v is not in the rowspace of H_Z.
 *
 *  Returns the Z-side distance `d_Z`. */
IRREP_API int
irrep_css_code_compute_logical_Z(const irrep_css_code_t *c, int max_weight,
                                  irrep_pauli_t *out_pauli);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_CSS_CODE_H */
