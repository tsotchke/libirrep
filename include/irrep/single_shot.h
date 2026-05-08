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

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SINGLE_SHOT_H */
