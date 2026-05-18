/* SPDX-License-Identifier: MIT */
/** @file toric4d.h
 *  @brief 4D toric code on the hypercubic lattice with periodic boundaries.
 *
 *  The 4D toric code [Dennis-Kitaev-Landahl-Preskill 2002] is the natural
 *  extension of Kitaev's 2D toric code (qubits on edges) and the 3D toric
 *  code (qubits on edges) to four dimensions, with qubits on **2-cells
 *  (faces)** of the 4-torus. Both stabilizer types have weight 6 — a
 *  symmetric structure unavailable in lower dimensions — which makes it
 *  the smallest dimension in which self-correcting topological order is
 *  rigorously proven (Alicki-Horodecki-Horodecki-Horodecki 2010).
 *
 *  ## Lattice
 *
 *  - **Qubits on 2-cells (faces)**.  In a 4D hypercube there are
 *    `C(4, 2) = 6` face orientations `(i, j)` with `0 ≤ i < j ≤ 3`,
 *    so `n = 6·L⁴`.
 *  - **X-stabilizers on 1-cells (edges)**.  An edge of orientation `d`
 *    bounds the 3 face orientations `(d, j)` for `j ≠ d`, and for each
 *    such face orientation it sits between 2 faces (front and back),
 *    giving 6 incident faces per edge.  There are `4·L⁴` edges (4 edge
 *    orientations per vertex) so `n_X = 4·L⁴`.
 *  - **Z-stabilizers on 3-cells (cubes)**.  A 3-cube of orientation
 *    `{a, b, c}` (= 3 of the 4 axes) is bounded by `C(3, 2) = 3` face
 *    orientations within `{a, b, c}`, each contributing 2 faces
 *    (front and back), giving 6 bounding faces.  There are `4·L⁴`
 *    cubes so `n_Z = 4·L⁴`.
 *
 *  Both stabilizer types thus have weight 6 — the 4D symmetric
 *  structure that 2D and 3D toric lack.
 *
 *  ## CSS orthogonality
 *
 *  Every (edge, cube) pair shares 0 or 2 faces:
 *    - If the edge's direction `d` is not among the cube's
 *      orientation `{a, b, c}`, the edge does not touch the cube
 *      (overlap 0).
 *    - If `d ∈ {a, b, c}`, the edge bounds the 3 face orientations
 *      `(d, j)` for `j ∈ {a,b,c} \ {d}` (= 2 face orientations); for
 *      each, the cube and edge share exactly the 1 face at the cube's
 *      `d`-position (= 2 faces total). Always even.
 *
 *  Hence `H_X · H_Z^T = 0 (mod 2)` in F₂ by construction.
 *
 *  ## Encoded qubits
 *
 *  `k = dim H_2(T⁴, ℤ₂) = C(4, 2) = 6` logical qubits — one per
 *  non-trivial 2-cycle of the 4-torus (each of the 6 face orientations
 *  gives a topologically distinct 2-cycle).
 *
 *  ## Edge / face / cube index packing
 *
 *  Vertices `(x, y, z, w)` linearise as `v = w·L³ + z·L² + y·L + x`,
 *  total `V = L⁴`.  Edges and faces and cubes index as
 *
 *      e(v, d)        = d   · V + v          (d ∈ {0..3})
 *      f(v, i, j)     = o(i, j) · V + v     (o(i,j) ∈ {0..5}, i < j)
 *      c(v, {a,b,c})  = o3({a,b,c}) · V + v (o3 ∈ {0..3}, sorted triple)
 *
 *  where the face/cube orientation enumeration is in lexicographic
 *  order:
 *      o(i,j):  (0,1)→0  (0,2)→1  (0,3)→2  (1,2)→3  (1,3)→4  (2,3)→5
 *      o3:      {0,1,2}→0 {0,1,3}→1 {0,2,3}→2 {1,2,3}→3
 *
 *  ## Self-correction at finite temperature
 *
 *  Unlike 2D / 3D toric codes, the 4D toric code is provably stable
 *  to finite-temperature noise below a thermodynamic threshold — the
 *  smallest known model of self-correcting quantum memory. The
 *  Dennis-Kitaev-Landahl-Preskill argument: a logical error requires
 *  flipping all faces along a non-trivial 2-cycle, whose perimeter is
 *  topologically protected and has free-energy cost scaling with
 *  system size at low temperature.
 *
 *  ## Primary references
 *
 *  - Dennis-Kitaev-Landahl-Preskill, *Topological quantum memory*,
 *    J. Math. Phys. 43 (2002) 4452 — original 4D toric code.
 *  - Alicki-Horodecki-Horodecki-Horodecki, *On thermal stability of
 *    topological qubit in Kitaev's 4D model*, Open Sys. Inf. Dyn. 17
 *    (2010) 1 — self-correction theorem.
 *  - Hamma-Lidar-Severini, *Entanglement and area law for the toric
 *    code on a torus*, Phys. Rev. A 79 (2009) 010102 — for the lower-
 *    dimensional analogue this generalises.
 */
#ifndef IRREP_TORIC4D_H
#define IRREP_TORIC4D_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief 4D toric-code geometry parameters. */
typedef struct {
    int Lx, Ly, Lz, Lw; /**< Lattice dimensions (each ≥ 2 for non-trivial homology). */
    int V;              /**< Vertex count = Lx · Ly · Lz · Lw. */
    int n_qubits;       /**< 6 · V (faces). */
    int n_X_stabs;      /**< 4 · V (edges). */
    int n_Z_stabs;      /**< 4 · V (3-cubes). */
} irrep_toric4d_params_t;

/** @brief Initialise 4D toric-code parameters on an Lx × Ly × Lz × Lw
 *  4-torus. Returns IRREP_ERR_INVALID_ARG if any dimension < 1. */
IRREP_API irrep_status_t
irrep_toric4d_init(irrep_toric4d_params_t *out,
                   int Lx, int Ly, int Lz, int Lw);

/** @brief Build the 4D toric code as a CSS code.
 *
 *  Stabilizer ordering:
 *    - H_X: row = edge index, ordered by `(orient ∈ {0..3}, w, z, y, x)`
 *      with the four edge orientations being x, y, z, w (= dim 0..3).
 *    - H_Z: row = 3-cube index, ordered by
 *      `(orient ∈ {0..3}, w, z, y, x)` with the four 3-cube orientations
 *      being {y,z,w}=0, {x,z,w}=1, {x,y,w}=2, {x,y,z}=3 (= "drop axis d"
 *      for d = 0..3, the natural 4D Hodge dual ordering of edges). */
IRREP_API irrep_status_t
irrep_toric4d_build(const irrep_toric4d_params_t *p, irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_TORIC4D_H */
