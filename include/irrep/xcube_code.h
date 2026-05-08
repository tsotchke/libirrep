/* SPDX-License-Identifier: MIT */
/** @file xcube_code.h
 *  @brief X-cube fracton code (Vijay-Haah-Fu 2016).
 *
 *  The X-cube model is the canonical **type-I fracton code**: a 3D
 *  stabilizer code whose excitations are *fractons* — point-like
 *  topological defects that can only move along sub-dimensional
 *  manifolds (single planes), not freely in 3D. This restricted
 *  mobility makes the X-cube model a candidate for **self-correcting
 *  quantum memory** and a foundational example of fracton physics.
 *
 *  ## Lattice
 *
 *  Cubic lattice with periodic boundaries; **qubits live on edges**
 *  (3 edges per unit cell, total `n = 3·Lx·Ly·Lz`).
 *
 *  ## Stabilizers
 *
 *  Two families of stabilizers (both of weight 4 in the bulk):
 *
 *  ### Vertex (X-type) "planar" stabilizers
 *
 *  At each vertex `v = (vx, vy, vz)`, three stabilizers — one for each
 *  coordinate plane:
 *
 *    `A_v^{xy}` = `X X X X` on the 4 edges incident to `v` lying in
 *      the xy-plane (i.e. the 2 x-edges at `(·, vy, vz)` and the 2
 *      y-edges at `(vx, ·, vz)`).
 *    `A_v^{xz}` = same in the xz-plane.
 *    `A_v^{yz}` = same in the yz-plane.
 *
 *  These are weight-4. The product `A_v^{xy} · A_v^{xz} · A_v^{yz}`
 *  equals the standard 6-edge vertex star of the 3D toric code — so the
 *  X-cube has *more* X-stabilizers than the 3D toric, and its X-fluxes
 *  (point-like defects from anti-commuting Z) are confined to planes.
 *
 *  ### Cube (Z-type) stabilizers
 *
 *  At each elementary cube `c = (cx, cy, cz)`, the product of `Z` over
 *  the 12 edges of the cube. Weight-12.
 *
 *  ## Fracton structure
 *
 *  Z-errors create localised excitations at corners of cubes — single-
 *  qubit Z flips four cube stabilizers, creating four corner-excitations
 *  that **cannot be moved freely**. The minimal way to move them
 *  requires creating additional excitations on a *plane*. This sub-
 *  dimensional mobility is the hallmark of fracton physics.
 *
 *  ## Counts
 *
 *  - `n_qubits = 3·Lx·Ly·Lz`
 *  - `n_vertex_stabs = 3·Lx·Ly·Lz` (3 planar stabs per vertex)
 *  - `n_cube_stabs = Lx·Ly·Lz`
 *  - Logical qubits: `k = O(L)` (sub-extensive! linear in L), reflecting
 *    a sub-extensive ground-state degeneracy — also a fracton hallmark.
 *
 *  ## Primary references
 *
 *  - Vijay-Haah-Fu, *Fracton topological order, generalized lattice
 *    gauge theory, and duality*, Phys. Rev. B 94 (2016) 235157
 *    [arXiv:1603.04442].
 *  - Pretko-Chen-You, *Fracton phases of matter*, Int. J. Mod. Phys.
 *    A 35 (2020) 2030003 — review.
 *  - Slagle-Kim, *Quantum field theory of X-cube fracton topological
 *    order*, Phys. Rev. B 96 (2017) 195139.
 */
#ifndef IRREP_XCUBE_CODE_H
#define IRREP_XCUBE_CODE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief X-cube model geometry parameters. */
typedef struct {
    int Lx, Ly, Lz;     /**< Lattice dimensions. */
    int n_qubits;       /**< 3 · Lx · Ly · Lz. */
    int n_vertex_stabs; /**< 3 · Lx · Ly · Lz (3 planar stabs per vertex). */
    int n_cube_stabs;   /**< Lx · Ly · Lz. */
} irrep_xcube_params_t;

/** @brief Initialise X-cube parameters for an `Lx × Ly × Lz` torus. */
IRREP_API irrep_status_t
irrep_xcube_init(irrep_xcube_params_t *out, int Lx, int Ly, int Lz);

/** @brief Build the X-cube fracton code as a CSS code.
 *
 *  Stabilizer ordering:
 *    - H_X has 3·Lx·Ly·Lz rows: planar stabs in order xy, xz, yz, each
 *      block in row-major (vz, vy, vx). Weight 4 each.
 *    - H_Z has Lx·Ly·Lz rows: cube stabs in row-major (cz, cy, cx).
 *      Weight 12 each. */
IRREP_API irrep_status_t
irrep_xcube_build(const irrep_xcube_params_t *p, irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_XCUBE_CODE_H */
