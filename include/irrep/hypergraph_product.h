/* SPDX-License-Identifier: MIT */
/** @file hypergraph_product.h
 *  @brief Hypergraph-product CSS qLDPC codes (Tillich-Zemor 2014).
 *
 *  Given two classical binary codes with parity-check matrices
 *  `H_a ∈ F₂^{m_1 × n_1}` and `H_b ∈ F₂^{m_2 × n_2}`, the **hypergraph
 *  product** is the CSS code with:
 *
 *    n_qubits = n_1·n_2 + m_1·m_2
 *    H_X      = [ H_a ⊗ I_{n_2}  |  I_{m_1} ⊗ H_b^T ]   (m_1·n_2 rows)
 *    H_Z      = [ I_{n_1} ⊗ H_b  |  H_a^T ⊗ I_{m_2} ]   (n_1·m_2 rows)
 *
 *  CSS orthogonality `H_X · H_Z^T = 0 (mod 2)` is automatic by
 *  Kronecker-product algebra:
 *
 *    H_X · H_Z^T = (H_a ⊗ I_{n_2})(I_{n_1} ⊗ H_b^T) + (I_{m_1} ⊗ H_b^T)(H_a ⊗ I_{m_2})
 *                = (H_a ⊗ H_b^T) + (H_a ⊗ H_b^T)
 *                = 0  (in characteristic 2).
 *
 *  ## Why this matters
 *
 *  The hypergraph product is *the* foundational construction for modern
 *  qLDPC codes. Known specialisations:
 *
 *    - If H_a = H_b = parity-check of a 1D classical repetition code on L bits:
 *      `H_a = (1 1 0 0 ... 0; 0 1 1 0 ... 0; ...; 0 0 ... 0 1 1)` (L-1 × L),
 *      then the hypergraph product is the **L × L toric code on a torus**
 *      (with periodic boundaries handled by replacing H_a with the cyclic
 *      circulant version).
 *
 *    - If H_a, H_b are LDPC parity-checks of classical codes with rate r
 *      and relative distance δ, the quantum code has rate ~ r² and
 *      relative distance Θ(δ √n). This was the first sub-linear distance
 *      qLDPC family; lifted product later improved to constant rate +
 *      constant relative distance.
 *
 *  Any classical LDPC parity-check matrix gives a hypergraph-product
 *  qLDPC code. This module accepts arbitrary `(H_a, H_b)` and emits the
 *  CSS code.
 *
 *  ## Primary references
 *
 *  - Tillich-Zemor, *Quantum LDPC codes with positive rate and minimum
 *    distance proportional to √n*, IEEE Trans. Inf. Theory 60 (2014) 1193
 *    [arXiv:0903.0566].
 *  - Bravyi-Hastings, *Homological product codes*, STOC 2014 — generalised
 *    homological framework.
 *  - Leverrier-Tillich-Zemor, *Quantum expander codes*, FOCS 2015 —
 *    decoding theorem.
 */
#ifndef IRREP_HYPERGRAPH_PRODUCT_H
#define IRREP_HYPERGRAPH_PRODUCT_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Build the hypergraph-product CSS code from two classical codes.
 *
 *  @param H_a  Classical parity-check matrix `m_1 × n_1`.
 *  @param H_b  Classical parity-check matrix `m_2 × n_2`.
 *  @param[out] out  CSS code with `n = n_1·n_2 + m_1·m_2` qubits.
 *
 *  Caller must `irrep_css_code_free` when done. */
IRREP_API irrep_status_t
irrep_hypergraph_product_build(const irrep_parity_matrix_t *H_a,
                               const irrep_parity_matrix_t *H_b,
                               irrep_css_code_t *out);

/* ====================================================================
 * Named hypergraph-product instances
 * ==================================================================== */

/** @brief Build the [[13, 1, 3]] hypergraph product of two [3, 1, 3]
 *  repetition codes.
 *
 *  Classical input: `H_a = H_b = [[1,1,0],[0,1,1]]` (2×3 parity-check
 *  matrix of the [3, 1, 3] repetition code). The HGP construction
 *  gives:
 *    - n_qubits = n_a · n_b + m_a · m_b = 9 + 4 = 13
 *    - m_X = m_a · n_b = 6
 *    - m_Z = n_a · m_b = 6
 *    - k_HGP = k_a · k_b + k_a^T · k_b^T = 1·1 + 0·0 = 1
 *    - d_HGP = min(d_a, d_b) = 3
 *
 *  This is the smallest non-trivial HGP instance and serves as a
 *  worked example of the construction. Caller must
 *  `irrep_css_code_free(out)` when done. */
IRREP_API irrep_status_t
irrep_hgp_repetition_3_13_1_3(irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_HYPERGRAPH_PRODUCT_H */
