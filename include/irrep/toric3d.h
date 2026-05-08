/* SPDX-License-Identifier: MIT */
/** @file toric3d.h
 *  @brief 3D toric code on the cubic lattice with periodic boundaries.
 *
 *  The 3D toric code [Castelnovo-Chamon 2008, Hamma-Lidar-Severini 2008]
 *  is the natural extension of Kitaev's 2D toric code to three
 *  dimensions. On a 3-torus T³ = ℤ/Lx × ℤ/Ly × ℤ/Lz it encodes
 *  `k = 3` logical qubits — one per non-contractible 1-cycle. The code
 *  is the gauge theory of ℤ_2 1-form symmetry.
 *
 *  ## Lattice
 *
 *  - **Qubits live on edges** of the cubic lattice. There are 3 edges
 *    per unit cell (orientations x, y, z), so `n = 3·Lx·Ly·Lz`.
 *  - **Vertex stabilizers (X-type)**: at each vertex `v`, the product of
 *    `X` over the 6 edges incident to `v`. There are `Lx·Ly·Lz` of them,
 *    each of weight 6.
 *  - **Plaquette stabilizers (Z-type)**: on each unit-cell face, the
 *    product of `Z` over the 4 edges bounding that face. There are
 *    `3·Lx·Ly·Lz` of them (3 face orientations: xy, yz, xz), each of
 *    weight 4.
 *
 *  ## CSS orthogonality
 *
 *  Each `xy`-plaquette has 4 corner vertices; each corner contributes
 *  exactly 2 edges shared with the vertex's star. Non-corner vertices
 *  share 0 edges. Hence every (vertex, plaquette) overlap is even, and
 *  `H_X · H_Z^T = 0` in F₂.
 *
 *  ## Encoded qubits and redundancies
 *
 *  - One linear redundancy in the vertex stabilizers: the product of
 *    all of them equals the identity (each edge appears twice).
 *  - For each cube, the product of its 6 face stabilizers equals the
 *    identity (each edge appears twice). This gives `Lx·Ly·Lz` plaquette
 *    redundancies.
 *  - There are 3 logical X strings and 3 logical Z surfaces (one of each
 *    along each torus 1-cycle).
 *  - Net: `k = n - (n_X - 1) - (n_Z - Lx·Ly·Lz) = 3` logical qubits.
 *
 *  ## Edge index packing
 *
 *  An edge is `(x, y, z, orient)` with `orient ∈ {0, 1, 2}` for
 *  x-/y-/z-directed edges. Linear index:
 *
 *      e(x, y, z, orient) = orient · (Lx·Ly·Lz)
 *                          + z · (Lx·Ly) + y · Lx + x.
 *
 *  ## Primary references
 *
 *  - Castelnovo-Chamon, *Topological order in a three-dimensional toric
 *    code at finite temperature*, Phys. Rev. B 78 (2008) 155120.
 *  - Hamma-Lidar-Severini, *Entanglement and area law for the toric code
 *    on a torus*, Phys. Rev. A 79 (2009) 010102.
 *  - Bombín-Martín-Delgado, *Topological computation without braiding*,
 *    Phys. Rev. Lett. 98 (2007) 160502 — 3D color code analogue.
 */
#ifndef IRREP_TORIC3D_H
#define IRREP_TORIC3D_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief 3D toric code geometry parameters. */
typedef struct {
    int Lx, Ly, Lz;     /**< Lattice dimensions. */
    int n_qubits;       /**< 3 · Lx · Ly · Lz. */
    int n_vertex_stabs; /**< Lx · Ly · Lz. */
    int n_face_stabs;   /**< 3 · Lx · Ly · Lz. */
} irrep_toric3d_params_t;

/** @brief Initialise 3D-toric parameters for an `Lx × Ly × Lz` torus
 *  (each ≥ 1). For non-trivial homology require all ≥ 2. */
IRREP_API irrep_status_t
irrep_toric3d_init(irrep_toric3d_params_t *out, int Lx, int Ly, int Lz);

/** @brief Build the 3D toric code as a CSS code. Stabilizer ordering:
 *  H_X has `Lx·Ly·Lz` vertex stars (row-major (z, y, x)); H_Z has the
 *  3·Lx·Ly·Lz face stabilizers ordered by face-orientation (xy, then
 *  yz, then xz), each block in row-major. */
IRREP_API irrep_status_t
irrep_toric3d_build(const irrep_toric3d_params_t *p, irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_TORIC3D_H */
