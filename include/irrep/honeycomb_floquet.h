/* SPDX-License-Identifier: MIT */
/** @file honeycomb_floquet.h
 *  @brief Hastings-Haah honeycomb Floquet code (PRX Quantum 2 (2021) 030328).
 *
 *  The honeycomb Floquet code lives on the **brick-wall** representation
 *  of the honeycomb lattice. Edges of the lattice are 3-colored (red, green,
 *  blue) such that the three edges incident to each vertex carry the three
 *  different colors. The code is a 3-round Floquet schedule:
 *
 *    - Round 0 (red):   measure XX on every red edge.
 *    - Round 1 (green): measure YY on every green edge.
 *    - Round 2 (blue):  measure ZZ on every blue edge.
 *
 *  Two key properties:
 *    1. Within each round, measurements are on disjoint edges, so they
 *       trivially pairwise commute (they touch disjoint qubits).
 *    2. Across rounds, every consecutive-round edge pair shares exactly
 *       one vertex; the two operators (e.g. XX and YY) restrict to X·Y
 *       on that vertex, which anti-commutes — so consecutive-round
 *       measurements anti-commute on shared qubits, which is what drives
 *       the Floquet dynamics.
 *
 *  ## Brick-wall layout
 *
 *  Following Hastings-Haah, we adopt the brick-wall convention: the
 *  honeycomb lattice is laid out as alternating rows of horizontal
 *  "bricks". For an `Lx × Ly` patch with periodic boundaries, both `Lx`
 *  and `Ly` should be even. The total qubit count is `2·Lx·Ly` (one per
 *  vertex). Each round contributes `Lx·Ly` measurements (one per edge of
 *  that color); the 3-coloring means there are `3·Lx·Ly` edges in total,
 *  matching the 3-valence count.
 *
 *  Coordinate convention:
 *    Each unit cell has 2 vertices: A-sublattice at (x, y, 0) and
 *    B-sublattice at (x, y, 1) for x ∈ [0, Lx), y ∈ [0, Ly). Linear
 *    qubit index: `q(x, y, s) = 2·(y·Lx + x) + s`.
 *
 *  Edge colors at unit cell (x, y):
 *    - **red edge** (vertical within cell): A(x,y,0) ↔ B(x,y,1)
 *    - **green edge** (horizontal A→B+x): A(x,y,0) ↔ B((x-1) mod Lx, y, 1)
 *    - **blue edge** (vertical between cells): A(x,y,0) ↔ B(x, (y-1) mod Ly, 1)
 *
 *  Each B-vertex receives 3 edges, one of each color, from the 3
 *  neighbouring A-vertices. Each A-vertex emits 3 edges, one of each
 *  color, to its 3 neighbouring B-vertices.
 *
 *  ## Primary references
 *
 *  - Hastings-Haah, *Dynamically generated logical qubits*, PRX Quantum 2
 *    (2021) 030328 [arXiv:2107.02194].
 *  - Haah-Hastings, *Boundaries for the honeycomb code*, Quantum 6
 *    (2022) 693.
 *  - Gidney-Newman-Brooks-McEwen, *A Fault-Tolerant Honeycomb Memory*,
 *    Quantum 5 (2021) 605 — circuit-level threshold ~0.2-0.3%.
 */
#ifndef IRREP_HONEYCOMB_FLOQUET_H
#define IRREP_HONEYCOMB_FLOQUET_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/floquet_code.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Build the Hastings-Haah honeycomb Floquet code on an Lx×Ly
 *  brick-wall patch with periodic boundaries.
 *
 *  Both `Lx` and `Ly` must be ≥ 1; for honest boundary-free dynamics
 *  they should typically be ≥ 2. The output Floquet code has `n_qubits
 *  = 2·Lx·Ly`, `n_rounds = 3`, and each round has `Lx·Ly` 2-body Pauli
 *  measurements.
 *
 *  Caller must `irrep_floquet_code_free` when done. */
IRREP_API irrep_status_t
irrep_honeycomb_floquet_build(int Lx, int Ly, irrep_floquet_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_HONEYCOMB_FLOQUET_H */
