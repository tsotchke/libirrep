/* SPDX-License-Identifier: MIT */
/** @file point_group.h
 *  @brief Discrete-subgroup-of-O(3) character tables, projectors, and direct-sum
 *         decomposition under a point group.
 *
 *  Projector:
 *  @verbatim
 *    P_μ = (d_μ / |G|) · Σ_{g ∈ G}  χ*_μ(g) · D(g)
 *  @endverbatim
 *  where `d_μ` is the dimension of irrep μ of the point group, `|G|` is the
 *  group order, `χ_μ(g)` is the character at element g, and `D(g)` is the
 *  representation of g on the feature space specified by an
 *  @ref irrep_multiset_t. Irreps of the point group (e.g. `A₁, A₂, B₁, B₂, E`
 *  for C₄ᵥ) are distinct from SO(3) irreps (`l`) — the projector bridges them.
 *
 *  The character tables are verified against Bradley-Cracknell Vol. 1. For a
 *  pure rotation `g`, `D(g)` on an `(l, parity)` block is the real-basis
 *  Wigner-D matrix `D^l(R_g)`; improper elements pick up a `parity` factor
 *  (reflection = inversion · proper_rotation, inversion acts as `parity`
 *  on an `(l, parity)` irrep).
 *
 *  Groups supported: C₄ᵥ (square lattice, 5 irreps), D₆ (hexagonal /
 *  kagome, 6 irreps), C₃ᵥ (triangular with mirror symmetry, 3 irreps), and
 *  D₃ (purely rotational triangular, 3 irreps). C₃ᵥ and D₃ are isomorphic
 *  as abstract groups but differ in their geometric realisation — the
 *  former has improper reflections that flip parity-odd inputs, the
 *  latter is purely proper so `_reduce` on `1xLo` vs `1xLe` inputs can
 *  differ between them.
 */
#ifndef IRREP_POINT_GROUP_H
#define IRREP_POINT_GROUP_H

#include <stddef.h>

#include <irrep/export.h>
#include <irrep/multiset.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Discrete subgroup of O(3) identifier. */
typedef enum {
    IRREP_PG_C4V = 0, /**< C₄ᵥ — square-lattice symmetry.          */
    IRREP_PG_D6 = 1,  /**< D₆  — hexagonal / kagome lattice.       */
    IRREP_PG_C3V = 2, /**< C₃ᵥ — triangular lattice with σᵥ mirrors. */
    IRREP_PG_D3 = 3   /**< D₃  — triangular lattice, proper only.    */
} irrep_point_group_t;

/** @brief Opaque character table + element list for a point group. */
typedef struct irrep_pg_table irrep_pg_table_t;

/** @name Lifecycle
 *  @{ */

/** @brief Build a character table + element list for the given group.
 *  @return handle, or @c NULL on unknown / unimplemented group. */
IRREP_API irrep_pg_table_t *irrep_pg_table_build(irrep_point_group_t g);

/** @brief Release a table built by #irrep_pg_table_build. */
IRREP_API void irrep_pg_table_free(irrep_pg_table_t *t);
/** @} */

/** @name Metadata
 *  @{ */

/** @brief Number of irreducible representations of the group (5 for C₄ᵥ, 6 for D₆). */
IRREP_API int irrep_pg_num_irreps(const irrep_pg_table_t *t);

/** @brief Order of the group (8 for C₄ᵥ, 12 for D₆). */
IRREP_API int irrep_pg_order(const irrep_pg_table_t *t);

/** @brief Human-readable label for irrep index @p mu (e.g. "A1", "E").
 *         Pointer to static storage — do not free. */
IRREP_API const char *irrep_pg_irrep_label(const irrep_pg_table_t *t, int mu);
/** @} */

/** @name Projection and reduction
 *  @{ */

/** @brief Apply `P_μ` to a feature vector laid out according to @p spec.
 *  @param t       character table from @ref irrep_pg_table_build.
 *  @param mu      irrep index in `[0, num_irreps)`.
 *  @param spec    feature-space descriptor (block layout).
 *  @param in      input vector, length `spec->total_dim`.
 *  @param out     output, same length as @p in; components outside the
 *                 μ-th irreducible subspace are zeroed. In-place (`in == out`)
 *                 is not supported. */
IRREP_API void irrep_pg_project(const irrep_pg_table_t *t, int mu, const irrep_multiset_t *spec,
                                const double *in, double *out);

/** @brief Decompose @p spec as a direct sum under @p t. Writes
 *         `num_irreps(t)` multiplicities into @p out_mult — `out_mult[mu]`
 *         is the number of copies of μ in the decomposition of the
 *         representation carried by @p spec. */
IRREP_API void irrep_pg_reduce(const irrep_pg_table_t *t, const irrep_multiset_t *spec,
                               int *out_mult);
/** @} */

#ifdef __cplusplus
}
#endif

#endif /* IRREP_POINT_GROUP_H */
