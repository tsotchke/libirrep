/* SPDX-License-Identifier: MIT */
/** @file lattice.h
 *  @brief 2D Bravais-lattice primitives with periodic boundary conditions.
 *
 *  Supported lattices (@ref irrep_lattice_kind_t):
 *    - **Square** — one site per unit cell; `a₁ = (1, 0)`, `a₂ = (0, 1)`.
 *    - **Triangular** — one site per unit cell; `a₁ = (1, 0)`, `a₂ = (½, √3/2)`.
 *    - **Honeycomb** — two sites per unit cell; `a₁ = (3/2, √3/2)`, `a₂ = (3/2, −√3/2)`;
 *      sublattices A at `(0, 0)`, B at `(1, 0)`; NN bond length 1.
 *    - **Kagome** — three sites per unit cell; `a₁ = (2, 0)`, `a₂ = (1, √3)`;
 *      sublattices A at `(0, 0)`, B at `(1, 0)`, C at `(½, √3/2)`; NN bond length 1.
 *
 *  Site indices are packed `site = cell * sites_per_cell + sublattice`,
 *  `cell = ix + iy * Lx`, `ix ∈ [0, Lx)`, `iy ∈ [0, Ly)`. Periodic boundary
 *  conditions wrap all translations modulo `(Lx, Ly)`.
 *
 *  Nearest-neighbour (NN) and next-nearest-neighbour (NNN) bond lists are
 *  canonicalised (`i < j`) and deduplicated so a bond appears exactly once
 *  even on small lattices where PBC wraps a pair back onto itself.
 *
 *  The primitive-vector conventions above are chosen to match the standard
 *  references (Ashcroft-Mermin Chapter 7; Elser 1989 for kagome). They set
 *  the bond length to 1 on every NN pair of every supported lattice so that
 *  Heisenberg / Hubbard Hamiltonians expressed in units of J or t transfer
 *  verbatim.
 */
#ifndef IRREP_LATTICE_H
#define IRREP_LATTICE_H

#include <stddef.h>
#include <stdbool.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Enumerated lattice families. */
typedef enum {
    IRREP_LATTICE_SQUARE      = 0,  /**< 1 site per unit cell. */
    IRREP_LATTICE_TRIANGULAR  = 1,  /**< 1 site per unit cell. */
    IRREP_LATTICE_HONEYCOMB   = 2,  /**< 2 sites per unit cell (A, B). */
    IRREP_LATTICE_KAGOME      = 3   /**< 3 sites per unit cell (A, B, C). */
} irrep_lattice_kind_t;

/** @brief Opaque lattice handle. Built by @ref irrep_lattice_build, released
 *         by @ref irrep_lattice_free. Bond lists and positions are cached at
 *         build time so downstream queries are O(1). */
typedef struct irrep_lattice irrep_lattice_t;

/** @brief Construct a `Lx × Ly`-cell lattice of the given kind.
 *
 *  Both `Lx` and `Ly` must be ≥ 2 to avoid self-bond pathologies at PBC
 *  wrap; for kagome, `Ly` should additionally be even when the caller
 *  wants a cluster that respects the full p6mm symmetry, but the builder
 *  accepts odd `Ly`.
 *
 *  @return a fresh handle, or `NULL` on bad input / OOM; the reason is
 *          available via @ref irrep_last_error. */
IRREP_API irrep_lattice_t *irrep_lattice_build(irrep_lattice_kind_t kind,
                                               int Lx, int Ly);

/** @brief Release a lattice handle. */
IRREP_API void             irrep_lattice_free (irrep_lattice_t *L);

/** @name Structural queries
 *  @{ */

/** @brief Total site count, `Lx * Ly * sites_per_cell`. */
IRREP_API int irrep_lattice_num_sites     (const irrep_lattice_t *L);

/** @brief Unit-cell count, `Lx * Ly`. */
IRREP_API int irrep_lattice_num_cells     (const irrep_lattice_t *L);

/** @brief Sites per primitive unit cell: 1 (square/triangular), 2 (honeycomb), 3 (kagome). */
IRREP_API int irrep_lattice_sites_per_cell(const irrep_lattice_t *L);

/** @brief Number of cells along `a₁`. */
IRREP_API int irrep_lattice_Lx            (const irrep_lattice_t *L);

/** @brief Number of cells along `a₂`. */
IRREP_API int irrep_lattice_Ly            (const irrep_lattice_t *L);

/** @brief Lattice kind, as passed to @ref irrep_lattice_build. */
IRREP_API irrep_lattice_kind_t irrep_lattice_kind(const irrep_lattice_t *L);
/** @} */

/** @name Geometry
 *  @{ */

/** @brief Write primitive vectors @p a1, @p a2 in cartesian coordinates. */
IRREP_API void irrep_lattice_primitive_vectors(const irrep_lattice_t *L,
                                               double a1[2], double a2[2]);

/** @brief Write reciprocal primitive vectors @p b1, @p b2 satisfying
 *         `a_i · b_j = 2π δ_{ij}`. */
IRREP_API void irrep_lattice_reciprocal_vectors(const irrep_lattice_t *L,
                                                double b1[2], double b2[2]);

/** @brief Cartesian position of @p site. Returns #IRREP_OK or
 *         #IRREP_ERR_INVALID_ARG. */
IRREP_API irrep_status_t irrep_lattice_site_position(const irrep_lattice_t *L,
                                                     int site, double xy[2]);

/** @brief Sublattice index `[0, sites_per_cell)` for @p site; `-1` on error. */
IRREP_API int irrep_lattice_sublattice_of(const irrep_lattice_t *L, int site);

/** @brief Cell coordinates `(ix, iy)` of @p site. Returns #IRREP_OK or
 *         #IRREP_ERR_INVALID_ARG. */
IRREP_API irrep_status_t irrep_lattice_cell_of(const irrep_lattice_t *L,
                                               int site, int *ix, int *iy);

/** @brief Compose a site index from `(ix, iy, sublattice)`. PBC is applied
 *         to `ix, iy`; returns `-1` on out-of-range sublattice. */
IRREP_API int irrep_lattice_site_index(const irrep_lattice_t *L,
                                       int ix, int iy, int sublattice);
/** @} */

/** @name Translations
 *  @{ */

/** @brief Apply the primitive translation `T_{dx·a₁ + dy·a₂}` to @p site.
 *         Returns the new site index; `-1` on invalid input. PBC applied. */
IRREP_API int irrep_lattice_translate(const irrep_lattice_t *L, int site,
                                      int dx, int dy);
/** @} */

/** @name Bonds
 *  @{ */

/** @brief Number of nearest-neighbour bonds (after canonicalisation / dedup). */
IRREP_API int irrep_lattice_num_bonds_nn (const irrep_lattice_t *L);

/** @brief Number of next-nearest-neighbour bonds. */
IRREP_API int irrep_lattice_num_bonds_nnn(const irrep_lattice_t *L);

/** @brief Write the NN bond list into `[i_out, j_out]` arrays, each sized
 *         `irrep_lattice_num_bonds_nn(L)`. Either pointer may be `NULL` if
 *         only one endpoint array is needed. Each bond satisfies `i < j`. */
IRREP_API void irrep_lattice_fill_bonds_nn (const irrep_lattice_t *L,
                                            int *i_out, int *j_out);

/** @brief Write the NNN bond list; semantics match @ref irrep_lattice_fill_bonds_nn. */
IRREP_API void irrep_lattice_fill_bonds_nnn(const irrep_lattice_t *L,
                                            int *i_out, int *j_out);
/** @} */

/** @name Brillouin-zone grids
 *  @{ */

/** @brief Fill a `Lx × Ly` k-point grid, row-major, into @p kx, @p ky
 *         (each of length `num_cells`). Allowed-momentum grid `k = (n₁/Lx) b₁
 *         + (n₂/Ly) b₂` with `n₁ ∈ [0, Lx)`, `n₂ ∈ [0, Ly)`. */
IRREP_API void irrep_lattice_k_grid(const irrep_lattice_t *L,
                                    double *kx, double *ky);
/** @} */

#ifdef __cplusplus
}
#endif

#endif /* IRREP_LATTICE_H */
