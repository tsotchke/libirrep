/* SPDX-License-Identifier: MIT */
/** @file css_floquet.h
 *  @brief CSS Floquet codes (Davydova-Tantivasadakarn-Balasubramanian 2023).
 *
 *  CSS Floquet codes are a generalisation of the Hastings-Haah honeycomb
 *  code that:
 *
 *    1. Have a CSS structure within each round (each round measures a
 *       set of `XX` operators AND a separate set of `ZZ` operators on
 *       specified edges; the sets within a round commute pairwise).
 *    2. Have explicit (non-gauge) logical operators surviving the cycle.
 *    3. Generalise to arbitrary 3-edge-colourable lattices.
 *
 *  ## CSS Floquet on the square lattice
 *
 *  We implement the simplest CSS Floquet variant: an `Lx × Ly` square
 *  lattice with periodic boundaries. Each unit cell has 2 qubits (one
 *  per sublattice A, B) and 2 edges (horizontal, vertical). The 4-edge
 *  schedule cycles through:
 *
 *    Round 0: measure `X_A · X_B` on every horizontal edge.
 *    Round 1: measure `Z_A · Z_B` on every horizontal edge.
 *    Round 2: measure `X_A · X_B` on every vertical edge.
 *    Round 3: measure `Z_A · Z_B` on every vertical edge.
 *
 *  At each round, all measurements within the round commute (they're on
 *  disjoint edges of one orientation). Across rounds, the `XX → ZZ`
 *  transition on the same edge anti-commutes (driving the Floquet
 *  dynamics), and the rotation horizontal → vertical anti-commutes at
 *  shared vertices.
 *
 *  ## Why CSS Floquet matters
 *
 *  The CSS structure makes the code easier to decode (X- and Z-errors
 *  decouple) and easier to compose with CSS resource-state distillation.
 *  Davydova et al. show that CSS Floquet codes can implement *all*
 *  Clifford gates via measurement schedules alone, with no two-qubit
 *  gates ever applied.
 *
 *  ## Counts (Lx × Ly periodic)
 *
 *  - `n_qubits = 2 · Lx · Ly`
 *  - `n_rounds = 4`
 *  - `n_meas` per round = `Lx · Ly` (one weight-2 measurement per edge
 *    of the appropriate orientation)
 *
 *  ## Primary references
 *
 *  - Davydova-Tantivasadakarn-Balasubramanian, *Floquet codes without
 *    parent subsystem codes*, PRX Quantum 4 (2023) 020341
 *    [arXiv:2210.02468].
 *  - Aasen-Wang-Hastings, *Adiabatic paths and CSS Floquet codes*,
 *    PRX Quantum 4 (2023) 020339.
 *  - Kesselring-Magdalena-Bauer-Tantivasadakarn-Balasubramanian-Davydova,
 *    *Anyon condensation and the color code*, PRX Quantum 5 (2024) —
 *    related construction.
 */
#ifndef IRREP_CSS_FLOQUET_H
#define IRREP_CSS_FLOQUET_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/floquet_code.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Build the 4-round square-lattice CSS Floquet code on Lx × Ly
 *  with periodic boundaries.
 *
 *  Round schedule:
 *    0: XX on horizontal edges
 *    1: ZZ on horizontal edges
 *    2: XX on vertical edges
 *    3: ZZ on vertical edges
 *
 *  Caller must `irrep_floquet_code_free` when done. */
IRREP_API irrep_status_t
irrep_css_floquet_square_build(int Lx, int Ly, irrep_floquet_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_CSS_FLOQUET_H */
