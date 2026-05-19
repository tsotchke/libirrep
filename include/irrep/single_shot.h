/* SPDX-License-Identifier: MIT */
/** @file single_shot.h
 *  @brief Single-shot qLDPC codes (Quintavalle-Vasmer-Roffe-Campbell 2021).
 *
 *  ## What "single-shot" means
 *
 *  In a generic stabilizer code, syndrome extraction is itself faulty:
 *  measuring the stabilizers gives noisy syndrome bits, so one ordinary
 *  *round* of measurements is not enough to reliably correct errors.
 *  Standard fault-tolerance protocols therefore repeat syndrome
 *  extraction `O(d)` times per logical clock, multiplying the
 *  experimental cost by a factor of `d`.
 *
 *  A **single-shot** decoder breaks this trade-off: a single noisy
 *  syndrome round suffices for correction, with the residual error
 *  bounded by a constant. The mechanism is **meta-checks** on the
 *  syndrome itself: linear combinations of stabilizer measurements that
 *  must equal zero in the absence of measurement errors.
 *
 *  Concretely, given a CSS code with check matrices `(H_X, H_Z)`, a
 *  single-shot lift adds **meta-check matrices** `(M_X, M_Z)` satisfying:
 *
 *      M_X · H_X = 0  (mod 2)
 *      M_Z · H_Z = 0  (mod 2)
 *
 *  Then any "syndrome vector" `s = H · e + ε` (where `e` is the data
 *  error and `ε` is the measurement error) satisfies `M · s = M · ε`,
 *  so `M · s` is a *witness* of the measurement error — independent of
 *  the data error.
 *
 *  ## What this module exposes
 *
 *  We provide:
 *    - `irrep_single_shot_code_t` — a CSS code plus its X- and Z-meta-
 *      check matrices.
 *    - `irrep_single_shot_verify_meta` — checks `M_X · H_X = 0` and
 *      `M_Z · H_Z = 0` over F₂.
 *    - `irrep_single_shot_meta_syndrome` — given a syndrome vector,
 *      compute its meta-syndrome (independent of any data error).
 *
 *  Constructions of `M_X, M_Z` for specific CSS codes (3D toric,
 *  hypergraph product, BB) are out-of-scope here — they require
 *  code-family-specific arguments. This module is the *framework* and
 *  is verified against the **3D toric code**, where `M_X` is the
 *  cube-redundancy matrix: each cube has 6 face stabilizers whose
 *  product equals identity, providing a natural meta-check.
 *
 *  ## Primary references
 *
 *  - Quintavalle-Vasmer-Roffe-Campbell, *Single-shot error correction
 *    of three-dimensional homological product codes*, PRX Quantum 2
 *    (2021) 020340.
 *  - Bombín, *Single-shot fault-tolerant quantum error correction*,
 *    Phys. Rev. X 5 (2015) 031043 — original single-shot framework.
 *  - Campbell, *A theory of single-shot error correction for adversarial
 *    noise*, Quantum Sci. Technol. 4 (2019) 025006.
 */
#ifndef IRREP_SINGLE_SHOT_H
#define IRREP_SINGLE_SHOT_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief A CSS code lifted to single-shot via X- and Z-meta-checks.
 *
 *  - `M_X` rows are linear combinations of `H_X` rows that sum to zero.
 *  - `M_Z` rows are linear combinations of `H_Z` rows that sum to zero. */
typedef struct {
    irrep_css_code_t css;        /**< Underlying CSS code. */
    irrep_parity_matrix_t M_X;   /**< X-meta-check: M_X · H_X = 0. */
    irrep_parity_matrix_t M_Z;   /**< Z-meta-check: M_Z · H_Z = 0. */
} irrep_single_shot_code_t;

/** @brief Allocate a single-shot code wrapping an existing CSS code.
 *
 *  Takes ownership of the CSS code's storage (consumer of `css` must
 *  not free it independently). M_X has `m_X_meta` rows × `H_X.n_rows`
 *  cols; M_Z has `m_Z_meta` rows × `H_Z.n_rows` cols. */
IRREP_API irrep_status_t
irrep_single_shot_code_new(irrep_single_shot_code_t *out,
                           irrep_css_code_t *css,
                           int m_X_meta, int m_Z_meta);

/** @brief Free internal storage. */
IRREP_API void
irrep_single_shot_code_free(irrep_single_shot_code_t *c);

/** @brief Verify the meta-check property: M_X · H_X = 0 and M_Z · H_Z = 0
 *  over F₂. Returns IRREP_OK iff both hold. */
IRREP_API irrep_status_t
irrep_single_shot_verify_meta(const irrep_single_shot_code_t *c);

/** @brief Apply M_X to an X-syndrome vector (length H_X.n_rows; bit-packed
 *  into `n_words = ceil(rows/64)` uint64_t words). Output `meta_syndrome`
 *  has length M_X.n_rows. */
IRREP_API irrep_status_t
irrep_single_shot_meta_syndrome_X(const irrep_single_shot_code_t *c,
                                  const uint64_t *syndrome,
                                  uint64_t *meta_syndrome);

/** @brief Same as above for Z-syndrome. */
IRREP_API irrep_status_t
irrep_single_shot_meta_syndrome_Z(const irrep_single_shot_code_t *c,
                                  const uint64_t *syndrome,
                                  uint64_t *meta_syndrome);

/* ====================================================================
 * Generic single-shot lifter
 *
 * For any CSS code, the meta-check matrices `(M_X, M_Z)` are the bases
 * of the left-nullspaces of `(H_X, H_Z)` respectively (i.e., the
 * F₂-linear redundancies between the stabilizer rows).
 *
 *   M_X · H_X = 0   ⟺   rows of M_X span { v ∈ F₂^{m_X} : v · H_X = 0 }
 *
 * For codes with no inter-stabilizer redundancy (Steane, [[5,1,3]],
 * minimal-distance surface code) the nullity is 0 and M_X / M_Z are
 * 0-row matrices — the code is not natively single-shot. For codes
 * with topological redundancy (3D toric, 4D toric, X-cube, 3D color
 * code) the nullity is positive and the lift gives a canonical
 * meta-check generator set.
 * ==================================================================== */

/** @brief Auto-compute meta-checks `(M_X, M_Z)` for a CSS code via
 *  F₂ left-nullspace.
 *
 *  Internally:
 *    - Deep-copies `css` into `out->css` (caller retains ownership of
 *      the original `css`, which can be independently freed).
 *    - Computes M_X as the basis of the left-nullspace of H_X using
 *      F₂ Gaussian elimination on the augmented matrix [H_X | I_{m_X}].
 *    - Computes M_Z analogously from H_Z.
 *    - `out->M_X.n_rows = m_X - rank(H_X)` (the nullity of H_X);
 *      `out->M_Z.n_rows = m_Z - rank(H_Z)`. Either may be 0.
 *
 *  Post-condition: `irrep_single_shot_verify_meta(out) == IRREP_OK`.
 *
 *  @param[in]  css  Source CSS code (read-only).
 *  @param[out] out  Caller-allocated; on success holds a deep copy of
 *                   `css` plus the computed meta-checks. Caller must
 *                   `irrep_single_shot_code_free(out)` when done.
 *  @return  `IRREP_OK` on success; `IRREP_ERR_OUT_OF_MEMORY` on allocation
 *           failure; `IRREP_ERR_INVALID_ARG` if either pointer is NULL. */
IRREP_API irrep_status_t
irrep_single_shot_lift(const irrep_css_code_t *css,
                       irrep_single_shot_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SINGLE_SHOT_H */
