/* SPDX-License-Identifier: MIT */
/** @file lattice3d.h
 *  @brief 3D Bravais-lattice primitives with periodic boundary conditions.
 *
 *  Sibling of @ref lattice.h. Conventions follow the **conventional cubic
 *  cell** with lattice constant `a = 1`:
 *
 *    - **Simple cubic (SC)** — 1 site per cell at `(0, 0, 0)`. NN distance 1,
 *      6 NN, 12 NNN at distance √2.
 *    - **Body-centered cubic (BCC)** — 2 sites per cell at `(0, 0, 0)` and
 *      `(½, ½, ½)`. NN distance √3/2, 8 NN; NNN distance 1, 6 NNN.
 *    - **Face-centered cubic (FCC)** — 4 sites per cell at `(0, 0, 0)`,
 *      `(0, ½, ½)`, `(½, 0, ½)`, `(½, ½, 0)`. NN distance √2/2, 12 NN;
 *      NNN distance 1, 6 NNN.
 *    - **Diamond** — 8 sites per cell: FCC plus the same FCC shifted by
 *      `(¼, ¼, ¼)`. NN distance √3/4, 4 NN (silicon / germanium / carbon
 *      tetrahedral coordination); NNN distance √2/2.
 *    - **Pyrochlore** — 16 sites per cell: FCC sublattices each carry the
 *      4-site tetrahedral basis `(0,0,0), (¼,¼,0), (¼,0,¼), (0,¼,¼)`.
 *      NN distance √2/4, **6 NN** (corner-sharing tetrahedra — the 3D
 *      analog of kagome's corner-sharing triangles); NNN distance √6/4.
 *      Hosts spin-ice and quantum-spin-liquid phases in materials like
 *      Tb₂Ti₂O₇, Yb₂Ti₂O₇, Dy₂Ti₂O₇.
 *
 *  Site indices are packed `site = cell · sites_per_cell + sublattice`,
 *  with `cell = ix + iy · Lx + iz · Lx · Ly`, `ix ∈ [0, Lx)`, `iy ∈ [0, Ly)`,
 *  `iz ∈ [0, Lz)`. Periodic boundary conditions wrap all translations
 *  modulo `(Lx, Ly, Lz)`.
 *
 *  Bond lists are canonicalised (`i < j`) and deduplicated, identical in
 *  semantics to the 2D variant. Both NN and NNN tables are cached at build
 *  time.
 *
 *  Bond lengths differ per lattice (1, √3/2, √2/2, √3/4) because every
 *  family is expressed in the same conventional cubic cell. This matches
 *  Ashcroft-Mermin Chapter 4 / Kittel Chapter 1; it differs from the 2D
 *  convention where every NN bond has length 1.
 */
#ifndef IRREP_LATTICE3D_H
#define IRREP_LATTICE3D_H

#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Enumerated 3D lattice families. */
typedef enum {
    IRREP_LATTICE3D_SC = 0,      /**< Simple cubic, 1 site per cell. */
    IRREP_LATTICE3D_BCC = 1,     /**< Body-centered cubic, 2 sites per cell. */
    IRREP_LATTICE3D_FCC = 2,     /**< Face-centered cubic, 4 sites per cell. */
    IRREP_LATTICE3D_DIAMOND = 3, /**< Diamond, 8 sites per cell (Si / Ge / C). */
    IRREP_LATTICE3D_PYROCHLORE = 4 /**< Corner-sharing tetrahedra, 16 sites per cell. */
} irrep_lattice3d_kind_t;

/** @brief Opaque 3D lattice handle. Built by @ref irrep_lattice3d_build,
 *         released by @ref irrep_lattice3d_free. Bond lists and positions
 *         are cached at build time so downstream queries are O(1). */
typedef struct irrep_lattice3d irrep_lattice3d_t;

/** @brief Construct an `Lx × Ly × Lz`-cell 3D lattice of the given kind.
 *
 *  All three dimensions must be ≥ 1. `L = 1` along any axis collapses all
 *  cell offsets along that axis to the same image; the bond builder dedups
 *  self-bonds, so SC at 1³ yields 0 bonds (every site is its own NN under
 *  PBC). Lattices with a non-trivial basis (BCC, FCC, Diamond, Pyrochlore)
 *  retain their intra-cell NN structure, so 1×1×1 is meaningful for those
 *  families — pyrochlore 1³ is a well-defined 16-site frustrated cluster.
 *  Consumers wanting full cubic-point-group symmetry (when 3D space groups
 *  land) will need `Lx = Ly = Lz`, which is not enforced by the builder.
 *
 *  @return a fresh handle, or `NULL` on bad input / OOM. The reason is
 *          available via @ref irrep_last_error. */
IRREP_API irrep_lattice3d_t *irrep_lattice3d_build(irrep_lattice3d_kind_t kind, int Lx, int Ly,
                                                   int Lz);

/** @brief Release a handle returned by @ref irrep_lattice3d_build. */
IRREP_API void irrep_lattice3d_free(irrep_lattice3d_t *L);

/** @name Structural queries
 *  @{ */

/** @brief Total site count, `Lx · Ly · Lz · sites_per_cell`. */
IRREP_API int irrep_lattice3d_num_sites(const irrep_lattice3d_t *L);

/** @brief Cell count, `Lx · Ly · Lz`. */
IRREP_API int irrep_lattice3d_num_cells(const irrep_lattice3d_t *L);

/** @brief Sites per primitive unit cell: 1 (SC), 2 (BCC), 4 (FCC), 8 (Diamond). */
IRREP_API int irrep_lattice3d_sites_per_cell(const irrep_lattice3d_t *L);

/** @brief Number of cells along each cubic axis. */
IRREP_API int irrep_lattice3d_Lx(const irrep_lattice3d_t *L);
IRREP_API int irrep_lattice3d_Ly(const irrep_lattice3d_t *L);
IRREP_API int irrep_lattice3d_Lz(const irrep_lattice3d_t *L);

/** @brief Lattice kind, as passed to @ref irrep_lattice3d_build. */
IRREP_API irrep_lattice3d_kind_t irrep_lattice3d_kind(const irrep_lattice3d_t *L);

/** @brief Geometric NN bond length (in conventional-cell units `a = 1`). */
IRREP_API double irrep_lattice3d_nn_distance(const irrep_lattice3d_t *L);

/** @brief Geometric NNN bond length. */
IRREP_API double irrep_lattice3d_nnn_distance(const irrep_lattice3d_t *L);
/** @} */

/** @name Geometry
 *  @{ */

/** @brief Write conventional-cell primitive vectors `a₁, a₂, a₃` (cartesian).
 *         The conventional cell is cubic with edge 1 for every supported
 *         family; the difference between SC, BCC, FCC, Diamond is the
 *         sublattice basis, not the unit cell. */
IRREP_API void irrep_lattice3d_primitive_vectors(const irrep_lattice3d_t *L, double a1[3],
                                                 double a2[3], double a3[3]);

/** @brief Write reciprocal vectors `b_i` satisfying `a_i · b_j = 2π δ_{ij}`. */
IRREP_API void irrep_lattice3d_reciprocal_vectors(const irrep_lattice3d_t *L, double b1[3],
                                                  double b2[3], double b3[3]);

/** @brief Cartesian position of @p site. Returns #IRREP_OK or
 *         #IRREP_ERR_INVALID_ARG. */
IRREP_API irrep_status_t irrep_lattice3d_site_position(const irrep_lattice3d_t *L, int site,
                                                       double xyz[3]);

/** @brief Sublattice index `[0, sites_per_cell)` for @p site; `-1` on error. */
IRREP_API int irrep_lattice3d_sublattice_of(const irrep_lattice3d_t *L, int site);

/** @brief Cell coordinates `(ix, iy, iz)` of @p site. */
IRREP_API irrep_status_t irrep_lattice3d_cell_of(const irrep_lattice3d_t *L, int site, int *ix,
                                                 int *iy, int *iz);

/** @brief Compose a site index from `(ix, iy, iz, sublattice)`. PBC is
 *         applied to `ix, iy, iz`; returns `-1` on out-of-range sublattice. */
IRREP_API int irrep_lattice3d_site_index(const irrep_lattice3d_t *L, int ix, int iy, int iz,
                                         int sublattice);
/** @} */

/** @name Translations
 *  @{ */

/** @brief Apply the primitive translation `T_{dx·a₁ + dy·a₂ + dz·a₃}`. */
IRREP_API int irrep_lattice3d_translate(const irrep_lattice3d_t *L, int site, int dx, int dy,
                                        int dz);
/** @} */

/** @name Bonds
 *  @{ */

/** @brief Number of nearest-neighbour bonds (after canonicalisation / dedup). */
IRREP_API int irrep_lattice3d_num_bonds_nn(const irrep_lattice3d_t *L);

/** @brief Number of next-nearest-neighbour bonds. */
IRREP_API int irrep_lattice3d_num_bonds_nnn(const irrep_lattice3d_t *L);

/** @brief Write the NN bond list into `[i_out, j_out]`, each sized
 *         `num_bonds_nn(L)`. Either pointer may be `NULL`. Each bond
 *         satisfies `i < j`. */
IRREP_API void irrep_lattice3d_fill_bonds_nn(const irrep_lattice3d_t *L, int *i_out, int *j_out);

/** @brief Write the NNN bond list; semantics match the NN variant. */
IRREP_API void irrep_lattice3d_fill_bonds_nnn(const irrep_lattice3d_t *L, int *i_out, int *j_out);
/** @} */

/** @name Brillouin-zone grids
 *  @{ */

/** @brief Fill an `Lx · Ly · Lz` k-point grid into @p kx, @p ky, @p kz
 *         (each of length `num_cells`). Allowed-momentum grid
 *         `k = (n₁/Lx) b₁ + (n₂/Ly) b₂ + (n₃/Lz) b₃` with
 *         `n_i ∈ [0, L_i)`. Index order `idx = n₁ + n₂ · Lx + n₃ · Lx · Ly`. */
IRREP_API void irrep_lattice3d_k_grid(const irrep_lattice3d_t *L, double *kx, double *ky,
                                      double *kz);
/** @} */

#ifdef __cplusplus
}
#endif

#endif /* IRREP_LATTICE3D_H */
