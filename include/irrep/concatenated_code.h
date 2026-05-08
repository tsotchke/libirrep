/* SPDX-License-Identifier: MIT */
/** @file concatenated_code.h
 *  @brief Concatenated CSS codes (two-level code composition).
 *
 *  Concatenation is the standard recursive construction in fault-tolerant
 *  quantum computing: take an *inner* `[[n_in, 1, d_in]]` code, encode each
 *  of the *outer* code's `n_out` physical qubits into one inner block, and
 *  use the outer code's stabilizers expressed as tensor products of inner
 *  logical operators.
 *
 *  Resulting parameters (when both codes are CSS and the inner logical
 *  operators X̄ and Z̄ are pure-X and pure-Z respectively):
 *
 *      [[n_in · n_out, k_out, ≥ d_in · d_out]]
 *
 *  ## Stabilizer structure
 *
 *  - **Inner stabilizers, blockwise**: for each of the `n_out` physical
 *    qubits of the outer code, copy the inner code's `m_X_in` X-stabs and
 *    `m_Z_in` Z-stabs onto the corresponding `n_in`-qubit block.
 *    Contributes `n_out · (m_X_in + m_Z_in)` stabilizers.
 *
 *  - **Outer stabilizers, lifted**: for each outer X-stabilizer of weight
 *    `w` (acting on `w` outer-code qubits), the corresponding concatenated
 *    stabilizer is the tensor product of inner X̄ on each of those `w`
 *    blocks. Same for Z. Contributes `m_X_out + m_Z_out` stabilizers.
 *
 *  ## Why the dimensions work
 *
 *      n_X_total = n_out · m_X_in + m_X_out
 *      n_Z_total = n_out · m_Z_in + m_Z_out
 *      k_total  = n_in · n_out - (n_X_total + n_Z_total)
 *               = n_out · (n_in - m_X_in - m_Z_in) - (m_X_out + m_Z_out)
 *               = n_out · k_in - (m_X_out + m_Z_out)
 *               = n_out - (m_X_out + m_Z_out)   (when k_in = 1)
 *               = k_out.
 *
 *  ## Why the distance is at least `d_in · d_out`
 *
 *  Any non-trivial logical of the concatenated code, viewed at the
 *  *outer* level, is a non-trivial logical of the outer code, hence has
 *  outer-weight ≥ d_out. Each of those outer qubits is encoded as an
 *  inner block; flipping the inner logical of any one block requires at
 *  least d_in physical errors. So total weight ≥ d_in · d_out.
 *
 *  ## Restriction
 *
 *  This module assumes the inner logical operators are *pure-X* (X̄) and
 *  *pure-Z* (Z̄), so the concatenated code remains CSS. This holds for
 *  all standard CSS code families (Steane, surface, color, BB, etc.).
 *
 *  ## Primary references
 *
 *  - Knill-Laflamme, *Concatenated quantum codes*, arXiv:quant-ph/9608012
 *    (1996) — original definition.
 *  - Aharonov-Ben-Or, *Fault-tolerant quantum computation with constant
 *    error rate*, SIAM J. Comput. 38 (2008) 1207 — threshold theorem
 *    via concatenation.
 *  - Bombín, *Topological codes by concatenation*, Phys. Rev. A 79
 *    (2009) 042322.
 */
#ifndef IRREP_CONCATENATED_CODE_H
#define IRREP_CONCATENATED_CODE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/stabilizer_group.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Construct the concatenated CSS code `inner ⊗ outer` with
 *  caller-supplied inner logical operators.
 *
 *  Pre-conditions:
 *    - `inner_x_bar` is pure-X (no Z bits set).
 *    - `inner_z_bar` is pure-Z (no X bits set).
 *    - `inner_x_bar` and `inner_z_bar` are both in the centraliser of
 *      the inner stabilizer group and anti-commute with each other.
 *    - `inner` encodes exactly 1 logical qubit (i.e. the supplied X̄/Z̄
 *      are unique up to the stabilizer group).
 *
 *  Output: `[[n_in · n_out, k_out, ≥ d_in · d_out]]` CSS code. Caller
 *  must `irrep_css_code_free` when done. */
IRREP_API irrep_status_t
irrep_css_concatenate(const irrep_css_code_t *inner,
                      const irrep_pauli_t *inner_x_bar,
                      const irrep_pauli_t *inner_z_bar,
                      const irrep_css_code_t *outer,
                      irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_CONCATENATED_CODE_H */
