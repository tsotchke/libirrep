/* SPDX-License-Identifier: MIT */
/** @file space_group.h
 *  @brief 2D wallpaper / space-group tables and site-permutation actions on
 *         a periodic lattice.
 *
 *  Supported wallpaper groups in 1.3.0-alpha:
 *    - **p1**   — translations only. Order = `num_cells`. Valid for every lattice.
 *    - **p4mm** — order 8 × `num_cells` (square lattice only).
 *    - **p6mm** — order 12 × `num_cells` (triangular / honeycomb / kagome).
 *
 *  The wallpaper group is realised as a site permutation π_g : site → site for
 *  each group element `g`. Tables are built once and cached per lattice, so
 *  downstream consumers (projected neural-quantum-state ansätze) see O(1)
 *  lookup. On a 6×6 kagome cluster (108 sites, p6mm order = 432) the full
 *  permutation table occupies 186 KB; applying one group element is a single
 *  `memcpy` of 108 integers.
 *
 *  ### PBC compatibility
 *  The space group only acts as a permutation of a *finite* lattice when the
 *  torus spanned by `(Lx·a1, Ly·a2)` is invariant under the full point group.
 *  For p4mm this requires `Lx = Ly`; for p6mm, `Lx = Ly` on triangular and
 *  kagome, and additionally that the 6-fold centre lies on the cluster. If
 *  the cluster breaks the symmetry, @ref irrep_space_group_build returns
 *  `NULL` and sets `irrep_last_error`.
 *
 *  ### Conventions
 *  For a rotation `R` and a translation `t = dx·a1 + dy·a2`, the group
 *  element `(t, R)` acts on a lattice site at cartesian position `r` as
 *
 *  \f[ g \cdot r = R \cdot (r - O) + O + t \f]
 *
 *  where `O` is the wallpaper-group origin (the lattice site for square /
 *  triangular, the hexagon centre for honeycomb / kagome). Elements are
 *  ordered as `g = t · N_point + p`, with the point-group index `p`
 *  enumerating `{E, C_n, C_n^2, …}` first, then mirrors.
 *
 *  Functions acting on configurations follow the pullback convention
 *  `(g·ψ)(r) = ψ(g⁻¹·r)`; use @ref irrep_space_group_permutation_inverse
 *  to recover `g⁻¹` as a permutation.
 */
#ifndef IRREP_SPACE_GROUP_H
#define IRREP_SPACE_GROUP_H

#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>
#include <irrep/lattice.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Enumerated wallpaper groups.
 *
 *  Lattice compatibility and the `Lx = Ly` constraint per group:
 *
 *  | Group    | Point order | Lattices                            | `Lx == Ly`? |
 *  |----------|------------:|-------------------------------------|:-----------:|
 *  | p1       | 1           | any                                 | no          |
 *  | p2       | 2           | any                                 | no          |
 *  | p4mm     | 8           | square                              | yes         |
 *  | p6mm     | 12          | triangular / honeycomb / kagome     | yes         |
 *  | p6       | 6           | triangular / honeycomb / kagome     | yes         |
 *  | p3m1     | 6           | triangular / honeycomb / kagome     | yes         |
 */
typedef enum {
    IRREP_WALLPAPER_P1 = 0,   /**< Translations only. */
    IRREP_WALLPAPER_P4MM = 1, /**< Square: C4 + 4 mirrors, 8 point ops. */
    IRREP_WALLPAPER_P6MM = 2, /**< Triangular / honeycomb / kagome: C6 + 6 mirrors, 12 point ops. */
    IRREP_WALLPAPER_P2 = 3,   /**< 180° rotation on any lattice; no `Lx = Ly` constraint. */
    IRREP_WALLPAPER_P6 = 4,   /**< Chiral hex: 6 rotations, no mirrors. */
    IRREP_WALLPAPER_P3M1 = 5, /**< Hex 3-fold + 3 mirrors through vertices. */
    IRREP_WALLPAPER_P4 = 6,   /**< Chiral square: 4 rotations, no mirrors. */
    IRREP_WALLPAPER_P31M = 7  /**< Hex 3-fold + 3 mirrors bisecting vertex pairs (30° offset from p3m1). */
} irrep_wallpaper_t;

/** @brief Opaque space-group handle. Built by @ref irrep_space_group_build
 *         and released by @ref irrep_space_group_free. */
typedef struct irrep_space_group irrep_space_group_t;

/** @brief Construct a space-group table for the given lattice and wallpaper group.
 *
 *  @return a fresh handle, or `NULL` if the cluster does not respect the
 *          full wallpaper group (e.g. `Lx != Ly` under p4mm). The reason
 *          is available via @ref irrep_last_error. */
IRREP_API irrep_space_group_t *irrep_space_group_build(const irrep_lattice_t *L,
                                                       irrep_wallpaper_t      kind);

/** @brief Release a space-group handle. */
IRREP_API void irrep_space_group_free(irrep_space_group_t *G);

/** @name Structural queries
 *  @{ */

/** @brief Total number of group elements, `point_order * num_cells`. */
IRREP_API int irrep_space_group_order(const irrep_space_group_t *G);

/** @brief Number of elements in the point subgroup (1, 8, or 12). */
IRREP_API int irrep_space_group_point_order(const irrep_space_group_t *G);

/** @brief Number of lattice sites the group acts on. */
IRREP_API int irrep_space_group_num_sites(const irrep_space_group_t *G);

/** @brief The wallpaper kind this handle was built with. */
IRREP_API irrep_wallpaper_t irrep_space_group_kind(const irrep_space_group_t *G);

/** @brief Borrow the lattice this space group was built over. The returned
 *         pointer is owned by the caller who built the lattice; do not free. */
IRREP_API const irrep_lattice_t *irrep_space_group_lattice(const irrep_space_group_t *G);
/** @} */

/** @name Group action on sites
 *  @{ */

/** @brief Apply group element @p g to @p site. Returns the image site index,
 *         or `-1` on invalid input. */
IRREP_API int irrep_space_group_apply(const irrep_space_group_t *G, int g, int site);

/** @brief Write the full permutation for @p g into @p perm[num_sites]:
 *         `perm[s] = g · s`. */
IRREP_API void irrep_space_group_permutation(const irrep_space_group_t *G, int g, int *perm);

/** @brief Write the inverse permutation `perm[g · s] = s`, i.e. `g⁻¹ · s`.
 *         Convenient for the pullback action on configurations. */
IRREP_API void irrep_space_group_permutation_inverse(const irrep_space_group_t *G, int g,
                                                     int *perm);

/** @brief Apply group element @p g to a real-valued configuration:
 *         `out[i] = in[g⁻¹ · i]` (pullback convention). `in` and `out` may
 *         alias only if `in == out` and the caller provides a scratch buffer
 *         — otherwise the permutation is written out-of-place. */
IRREP_API void irrep_space_group_apply_config(const irrep_space_group_t *G, int g, const double *in,
                                              double *out);
/** @} */

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SPACE_GROUP_H */
