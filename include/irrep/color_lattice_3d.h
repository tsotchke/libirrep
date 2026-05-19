/* SPDX-License-Identifier: MIT */
/** @file color_lattice_3d.h
 *  @brief Generic face-and-cell-list-driven 3D color code (3-colex).
 *
 *  A **3-colex** is a 4-valent, 4-colorable 3-cell complex: every
 *  vertex has 4 incident cells, every edge is shared by 3 faces and 3
 *  cells, and the cells admit a proper 4-coloring such that no two
 *  cells sharing a face have the same colour. Bombín 2007 (PRL 98
 *  160502) constructs the topological 3D color code on any 3-colex by:
 *
 *    - **Qubits at vertices.**
 *    - **X-stabilizer per face:** `X⊗qubits-on-face`.
 *    - **Z-stabilizer per 3-cell:** `Z⊗qubits-in-cell`.
 *
 *  The face-cell incidence on a 3-colex always gives even support
 *  overlap, so CSS orthogonality `H_X · H_Z^T = 0 (mod 2)` is
 *  automatic. The number of logical qubits is `k = dim H_2(M, ℤ_2)`
 *  of the underlying 3-manifold `M`; for `M = 3-simplex` (the
 *  tetrahedron) `k = 1`, for `M = T³` (the 3-torus) `k = 3`.
 *
 *  This header exposes the **generic 3-colex framework** that lets the
 *  caller plug in any face-and-cell list and obtain the corresponding
 *  CSS code. Concrete instances:
 *
 *    - `<irrep/color_code_3d.h>` `irrep_color_3d_cube_8_3_2`
 *      ([[8, 3, 2]] cubic, distance 2) — direct CSS construction.
 *    - `<irrep/color_code_3d.h>` `irrep_color_3d_rm_15_1_3`
 *      ([[15, 1, 3]] Reed-Muller / tetrahedral, distance 3) — direct
 *      CSS construction.
 *    - `irrep_color_lattice_3d_tetrahedron` (this header) — the
 *      [[15, 1, 3]] code rebuilt as a 3-colex face-and-cell list,
 *      proving the framework subsumes the direct construction.
 *
 *  ## Rectified-cubic-honeycomb 3-colex (arbitrary L)
 *
 *  The canonical Bombín 2007 family lives on the **bitruncated cubic
 *  honeycomb** (= the rectified-cubic 3-colex). Each cell is a
 *  truncated octahedron; the lattice tessellates ℝ³ with two
 *  4-colourable orbits of cells. For an `L × L × L` periodic
 *  truncated-octahedral 3-colex on `T³`:
 *
 *    - n_qubits = 12 L³  (vertices, each shared by 4 cells)
 *    - n_faces  = 14 L³  (squares + hexagons, each shared by 2 cells)
 *    - n_cells  =  2 L³  (truncated octahedra, two-coloured)
 *
 *  Distance scales linearly with L; for L=1 the structural counts above
 *  give a non-trivial small instance. Building the full vertex/face/
 *  cell-incidence tables for the rectified-cubic honeycomb at L≥1 is
 *  the natural application of this framework; the tables themselves
 *  are unit-cell-tessellation lookups (one `.c` file per L, or one
 *  parametric generator). The tetrahedral instance shipped here serves
 *  as the framework's correctness witness.
 *
 *  ## Primary references
 *
 *  - Bombín-Martín-Delgado, *Topological computation without
 *    braiding*, PRL 98 (2007) 160502.
 *  - Bombín, *An introduction to topological quantum codes*,
 *    arXiv:1311.0277 (2013) — survey including 3-colex constructions.
 */
#ifndef IRREP_COLOR_LATTICE_3D_H
#define IRREP_COLOR_LATTICE_3D_H

#include <stdbool.h>
#include <stddef.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Specification of a 3-colex as a face list + cell list.
 *
 *  Caller owns all arrays. The struct does NOT take ownership; the
 *  build function reads from it once. */
typedef struct {
    int n_qubits;                       /**< Number of vertices (qubits). */
    int n_faces;                        /**< Number of faces (= X-stabs). */
    const int *face_sizes;              /**< face_sizes[f] = weight of face f. */
    const int * const *face_qubits;     /**< face_qubits[f] = length-face_sizes[f]
                                         *   array of qubit indices. */
    int n_cells;                        /**< Number of 3-cells (= Z-stabs). */
    const int *cell_sizes;              /**< cell_sizes[c] = weight of cell c. */
    const int * const *cell_qubits;     /**< cell_qubits[c] = length-cell_sizes[c]
                                         *   array of qubit indices. */
} irrep_color_lattice_3d_t;

/** @brief Build a CSS code from a 3-colex face-and-cell list.
 *
 *  Allocates `out` as a CSS code with `m_X = lattice->n_faces` and
 *  `m_Z = lattice->n_cells`. For each face `f`, row `f` of `H_X` is
 *  the indicator of `face_qubits[f]`; for each cell `c`, row `c` of
 *  `H_Z` is the indicator of `cell_qubits[c]`.
 *
 *  After construction, CSS orthogonality is verified via
 *  `irrep_css_code_verify`. On any inconsistency (out-of-range qubit
 *  index, duplicate qubit in a single face/cell, or odd face-cell
 *  overlap), returns `IRREP_ERR_PRECONDITION` and frees the partially
 *  constructed code. Caller must `irrep_css_code_free(out)` on success.
 *
 *  @param[in]  lattice  Face + cell list. Must be non-NULL with both
 *                       face and cell tables present.
 *  @param[out] out      Caller-allocated `irrep_css_code_t` (uninitialised
 *                       on entry; this function takes care of allocation).
 *  @return  `IRREP_OK` on success;
 *           `IRREP_ERR_INVALID_ARG` on NULL inputs;
 *           `IRREP_ERR_PRECONDITION` on structural inconsistency. */
IRREP_API irrep_status_t
irrep_color_lattice_3d_to_css(const irrep_color_lattice_3d_t *lattice,
                              irrep_css_code_t *out);

/* ====================================================================
 * Tetrahedral 3-colex — the [[15, 1, 3]] worked example
 *
 * The 3-simplex with vertices `1..15 ↔ PG(3,2)` points, indexed
 * 0-based as qubits 0..14. The face-and-cell structure:
 *
 *   - 4 X-stabilizers (weight 8): one per coordinate bit b ∈ {0..3},
 *     support = {q : (q+1) has bit b set}. These are the "tetrahedron
 *     facets" in the 3-colex picture.
 *
 *   - 10 Z-stabilizers: 4 weight-8 "single-bit" supports (same as the
 *     X-stab supports, since C_X ⊂ C_Z in the punctured-RM
 *     construction) + 6 weight-4 "pair" supports (positions where two
 *     specific bits are both set).
 *
 * Generating these via the framework demonstrates the equivalence to
 * the direct CSS construction in `irrep_color_3d_rm_15_1_3`.
 * ==================================================================== */

/** @brief Build the tetrahedral 3-colex specification for the
 *  [[15, 1, 3]] code.
 *
 *  Allocates internal storage for the face / cell / size arrays and
 *  populates `out` to point at them. Caller must free with
 *  `irrep_color_lattice_3d_tetrahedron_free` when done.
 *
 *  Same code as `irrep_color_3d_rm_15_1_3` — this is the framework
 *  re-derivation. */
IRREP_API irrep_status_t
irrep_color_lattice_3d_tetrahedron(irrep_color_lattice_3d_t *out);

/** @brief Free the internal storage allocated by
 *  `irrep_color_lattice_3d_tetrahedron`. Safe to call on a lattice that
 *  was not built by this function (it inspects the pointer chain and
 *  no-ops if the heap-owned storage marker is absent). */
IRREP_API void
irrep_color_lattice_3d_tetrahedron_free(irrep_color_lattice_3d_t *lattice);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_COLOR_LATTICE_3D_H */
