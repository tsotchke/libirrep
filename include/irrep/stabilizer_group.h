/* SPDX-License-Identifier: MIT */
/** @file stabilizer_group.h
 *  @brief Abstract Pauli-string stabilizer group machinery.
 *
 *  This is the foundation layer for ALL stabilizer-code primitives in
 *  libirrep (toric, surface, color, bivariate bicycle, hypergraph
 *  product, Floquet, ...). Each individual code is just a specific
 *  choice of stabilizer support encoded in this representation.
 *
 *  ## Symplectic representation of Pauli operators
 *
 *  An n-qubit Pauli operator P is encoded as a pair of bit-vectors
 *  `(x, z) ∈ F_2^n × F_2^n`:
 *    - bit `x_i = 1` iff P has an X-component on qubit i
 *    - bit `z_i = 1` iff P has a Z-component on qubit i
 *
 *  In this *symplectic-only* convention `irrep_pauli_set(p, i, Y)` sets
 *  `x_i = z_i = 1` and **leaves the sign untouched** — so the resulting
 *  Pauli represents the operator `X_i · Z_i` (which equals `-i·Y_i` in
 *  physics convention). For signed Pauli arithmetic that recovers the
 *  textbook Y identities (e.g. Y² = +I), callers must apply the `+i`
 *  factor manually after setting Y bits, or work with X and Z directly.
 *  This trade-off keeps `irrep_pauli_set` reversible and side-effect
 *  free at the bit level.
 *
 *  Two Paulis P = (x_P, z_P) and Q = (x_Q, z_Q) commute if and only
 *  if the symplectic inner product is zero:
 *    `⟨P, Q⟩_symp = x_P · z_Q + z_P · x_Q (mod 2) = 0`
 *  where `·` is the F_2 inner product. This is the central fact this
 *  module exposes.
 *
 *  ## CSS codes
 *
 *  A *CSS code* is a stabilizer code where every stabilizer is either
 *  pure-X (i.e., `z = 0`) or pure-Z (`x = 0`). Equivalently, the
 *  parity-check matrices `H_X` (rows are X-stabilizer supports) and
 *  `H_Z` (rows are Z-stabilizer supports) satisfy
 *  `H_X · H_Zᵀ = 0 (mod 2)` — the row supports of `H_X` and `H_Z` are
 *  pairwise orthogonal in the F_2 sense.
 *
 *  Toric, surface, color, BB, hypergraph-product, lifted-product codes
 *  are all CSS. Honeycomb Floquet is *not* CSS in the static sense
 *  (its measurements are XX, YY, ZZ on different rounds), but each
 *  round's instantaneous stabilizer group is a CSS subgroup.
 *
 *  ## Primary references
 *
 *  - Calderbank-Shor, Phys. Rev. A 54 (1996) 1098 — CSS construction
 *  - Steane, Phys. Rev. Lett. 77 (1996) 793 — CSS construction
 *  - Gottesman PhD thesis (1997) — symplectic Pauli framework
 *  - Bravyi-Kitaev (1998) — toric code → CSS instance
 */
#ifndef IRREP_STABILIZER_GROUP_H
#define IRREP_STABILIZER_GROUP_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Sign of a Pauli operator in {+1, -1, +i, -i}.
 *
 *  Stored as an integer in {0, 1, 2, 3} corresponding to the power of i:
 *  0 → +1, 1 → +i, 2 → -1, 3 → -i. */
typedef enum {
    IRREP_PAULI_SIGN_POS = 0,  /**< +1 */
    IRREP_PAULI_SIGN_I   = 1,  /**< +i */
    IRREP_PAULI_SIGN_NEG = 2,  /**< -1 */
    IRREP_PAULI_SIGN_NEGI = 3, /**< -i */
} irrep_pauli_sign_t;

/** @brief A Pauli operator on n qubits in symplectic representation.
 *
 *  The bit-vectors `x` and `z` are stored as packed `uint64_t` arrays;
 *  bit i of word w corresponds to qubit `i + 64*w`. For n ≤ 64, this
 *  is a single word; for larger n, multiple words. -*/
typedef struct {
    int n;                /**< Number of qubits. */
    int n_words;          /**< Number of uint64_t words = ceil(n/64). */
    uint64_t *x;          /**< X-support (bit-packed). */
    uint64_t *z;          /**< Z-support (bit-packed). */
    irrep_pauli_sign_t sign; /**< Sign in {+1, +i, -1, -i}. */
} irrep_pauli_t;

/** @brief Allocate a new identity Pauli on n qubits (x = z = 0, sign = +1). */
IRREP_API irrep_status_t
irrep_pauli_new(irrep_pauli_t *out, int n);

/** @brief Free internal storage. (Does not free the struct itself.) */
IRREP_API void
irrep_pauli_free(irrep_pauli_t *p);

/** @brief Set a single qubit's Pauli letter: I, X, Y, Z. */
typedef enum {
    IRREP_PAULI_LETTER_I = 0,
    IRREP_PAULI_LETTER_X = 1,
    IRREP_PAULI_LETTER_Y = 2,
    IRREP_PAULI_LETTER_Z = 3,
} irrep_pauli_letter_t;

IRREP_API irrep_status_t
irrep_pauli_set(irrep_pauli_t *p, int qubit, irrep_pauli_letter_t letter);

/** @brief Get the Pauli letter on a single qubit. */
IRREP_API irrep_pauli_letter_t
irrep_pauli_get(const irrep_pauli_t *p, int qubit);

/** @brief Symplectic inner product of two Paulis (mod 2).
 *
 *  Returns 0 if they commute, 1 if they anti-commute. The inner product
 *  is `x_P · z_Q + z_P · x_Q (mod 2)`. This is the foundation of all
 *  stabilizer-code commutativity checks.
 *
 *  @return 0 (commute) or 1 (anti-commute), or -1 on error. */
IRREP_API int
irrep_pauli_symp_inner(const irrep_pauli_t *p, const irrep_pauli_t *q);

/** @brief Test whether two Paulis commute (ignoring signs). */
IRREP_API bool
irrep_pauli_commute(const irrep_pauli_t *p, const irrep_pauli_t *q);

/** @brief Hamming weight of a Pauli's support: number of qubits where
 *         its action is non-trivial (X, Y, or Z). */
IRREP_API int
irrep_pauli_weight(const irrep_pauli_t *p);

/* ====================================================================
 * Stabilizer group abstraction
 * ==================================================================== */

/** @brief A list of pairwise-commuting Pauli operators (a stabilizer
 *  group, possibly redundant; not necessarily independent).
 *
 *  This represents the *generators* of a stabilizer group; the actual
 *  group is closed under multiplication. For a [[n, k, d]] code with
 *  m stabilizer generators, m = n - k.
 *
 *  Toric L×L: n = 2L², k = 2, m = 2L² - 2 (one redundant vertex stabilizer
 *  + one redundant plaquette stabilizer).
 *  Surface code d: n = d² + (d-1)², k = 1, m = n - 1.
 *  BB[[144, 12, 12]]: n = 144, k = 12, m = 132. */
typedef struct {
    int n;                  /**< Number of physical qubits. */
    int n_generators;       /**< Number of stabilizer generators. */
    irrep_pauli_t *gens;    /**< Array of generators. */
} irrep_stabilizer_group_t;

/** @brief Allocate a stabilizer group with n_gens uninitialised generators
 *  on n qubits. */
IRREP_API irrep_status_t
irrep_stabilizer_group_new(irrep_stabilizer_group_t *out,
                           int n, int n_generators);

/** @brief Free internal storage. */
IRREP_API void
irrep_stabilizer_group_free(irrep_stabilizer_group_t *g);

/** @brief Verify that all generators pairwise commute.
 *
 *  Runs in O(m²) where m = n_generators. Returns IRREP_OK if the group
 *  is well-formed; IRREP_ERR_PRECONDITION if any pair anti-commutes. */
IRREP_API irrep_status_t
irrep_stabilizer_group_check_commutativity(const irrep_stabilizer_group_t *g);

/** @brief Compute the syndrome of an error operator.
 *
 *  For each generator g_i, the syndrome bit is `⟨g_i, error⟩_symp`.
 *  Returns 0 (no syndrome / commute) or 1 (syndrome / anticommute).
 *
 *  @param[out] syndrome  Bit-array of length g->n_generators
 *                         (caller-allocated, size ceil(m/64) words). */
IRREP_API irrep_status_t
irrep_stabilizer_syndrome(const irrep_stabilizer_group_t *g,
                          const irrep_pauli_t *error,
                          uint64_t *syndrome);

/** @brief Test whether a Pauli is in the F₂-linear span of the
 *  stabilizer generators (= an element of the stabilizer group,
 *  ignoring signs).
 *
 *  Concatenates the generators' symplectic vectors as rows of an
 *  `m × 2n` matrix over F₂ and runs Gaussian elimination on the
 *  augmented `(m+1) × 2n` system to test whether `p`'s symplectic
 *  vector is in the row-span. Complexity O(m · n · max(m, n) / 64).
 *
 *  Combined with `irrep_pauli_commute`, this lets callers distinguish
 *  the three cases for a Pauli p relative to the stabilizer group:
 *    - anti-commutes with some generator   → p ∉ N(<g>)         (returns 0)
 *    - commutes with all + in row-span     → p ∈ <g> (stabilizer) (returns 1)
 *    - commutes with all + not in span     → p ∈ N(<g>) \ <g>   (returns 0)
 *                                              (= a non-trivial logical)
 *
 *  @return  1 if `p` lies in the span, 0 if it does not, -1 on error.
 *           The function does NOT separately check commutativity;
 *           callers who want a full coset classification should call
 *           `irrep_pauli_commute` against every generator first. */
IRREP_API int
irrep_pauli_in_stabilizer_span(const irrep_stabilizer_group_t *g,
                               const irrep_pauli_t *p);

/** @brief Bell-pair contraction of qubits `a` and `b` in a stabilizer
 *  group. Projects the represented state onto the `|Φ+⟩` eigenspace of
 *  `(X_a X_b, Z_a Z_b)` and traces out qubits `a` and `b`, producing
 *  a stabilizer group on `(n - 2)` qubits.
 *
 *  ## Algorithm (standard stabilizer formalism)
 *
 *  Iterate over the two Pauli measurements `X_a X_b` and `Z_a Z_b`:
 *    1. Find a generator anti-commuting with the measured operator.
 *       If none, the operator is already in the centralizer; either
 *       it is in the span (skip) or adds a new logical-bound stabilizer
 *       (append).
 *    2. Multiply every other anti-commuting generator by the found
 *       one (symplectic Gaussian-elimination step), so that exactly
 *       one generator anti-commutes with the measured operator.
 *    3. Replace that generator with the measured operator.
 *  After both measurements, every generator either acts trivially on
 *  `{a, b}` or has action in the Bell-stabilizer group
 *  `{I, X⊗X, Z⊗Z, Y⊗Y}` on `{a, b}` (because it commutes with both
 *  `X_a X_b` and `Z_a Z_b`).
 *
 *  Trace-out: generators supported purely on `{a, b}` are dropped;
 *  the remaining generators have their bits on `{a, b}` zeroed and
 *  the qubits `a` and `b` removed from the index range.
 *
 *  ## Output sizing
 *
 *  `g_out` is allocated with `n - 2` qubits. The number of generators
 *  in `g_out` is at most `g_in->n_generators` (= unchanged in the
 *  typical case where neither measurement was redundant); slots beyond
 *  the actual reduced count are left as identity. Callers wanting the
 *  exact rank can post-process via `irrep_pauli_in_stabilizer_span` /
 *  Gaussian elimination.
 *
 *  ## Returns
 *
 *  `IRREP_OK` on success; `IRREP_ERR_INVALID_ARG` if inputs are NULL,
 *  if `a == b`, or if either is out of range. */
IRREP_API irrep_status_t
irrep_stabilizer_contract_bell(const irrep_stabilizer_group_t *g_in,
                               int a, int b,
                               irrep_stabilizer_group_t *g_out);

/** @brief Number of logical qubits encoded by a stabilizer group.
 *
 *  For an n-qubit stabilizer group with m generators, the code
 *  dimension is `k = n - rank(symplectic generators)`. The
 *  symplectic representation of a generator `(x, z) ∈ F₂^n × F₂^n`
 *  is the length-2n F₂ vector `(x|z)`, and the rank is over F₂.
 *
 *  Works for both CSS and non-CSS stabilizer codes. For CSS codes
 *  this returns the same value as `irrep_css_code_logical_qubits`
 *  on the equivalent CSS structure.
 *
 *  @param[in] g  Stabilizer group (read-only). All generators
 *                assumed to pairwise commute.
 *  @return  Number of logical qubits, or -1 on error. */
IRREP_API int
irrep_stabilizer_group_n_logical_qubits(const irrep_stabilizer_group_t *g);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_STABILIZER_GROUP_H */
