/* SPDX-License-Identifier: MIT */
/** @file happy_code.h
 *  @brief HaPPY holographic code (Pastawski-Yoshida-Harlow-Preskill 2015).
 *
 *  The HaPPY code is a tensor-network quantum error-correcting code on
 *  a hyperbolic-disk tiling, designed as an exact toy model of the
 *  AdS/CFT bulk-boundary correspondence. Logical (bulk) qubits live at
 *  the centre tensor; physical (boundary) qubits live on the outer
 *  edges of the tiling.
 *
 *  ## Construction
 *
 *  Each tile is a *perfect tensor* — a tensor `T_{i_1 ... i_n}` with
 *  the property that any partition of its indices into two halves
 *  yields an isometry from the smaller half to the larger. Concretely,
 *  the most-studied HaPPY model uses the **6-leg perfect tensor** built
 *  from the [[5, 1, 3]] **five-qubit code** (the smallest non-trivial
 *  CSS code, due to Laflamme-Miquel-Paz-Zurek 1996, and equivalent to
 *  the smallest perfect Reed-Muller / quantum-Hamming code):
 *
 *    1 logical qubit + 5 physical qubits = 6 legs in total.
 *
 *  Tiles are arranged on the {5,4} hyperbolic tiling: pentagons meeting
 *  4 at each vertex, an exponentially-growing graph in the radial
 *  direction. The interior tile has 1 leg pointing inward (the bulk
 *  logical qubit) and 5 legs pointing outward; tiles in the next layer
 *  have one leg connecting back to the interior, plus more tiles, and
 *  so on.
 *
 *  ## What this module exposes
 *
 *  This module provides the *primitive* — the **5-qubit code** — that
 *  is the workhorse perfect tensor of the HaPPY construction. The
 *  generalisation to a multi-tile network on the {5,4} tiling is
 *  conceptually straightforward but combinatorially involved (radial
 *  layers, contraction order, rate calculation), and is left to a
 *  future extension.
 *
 *  Interface:
 *    - `irrep_happy_perfect_tensor_5_1_3()` — build the canonical
 *      [[5, 1, 3]] code as a CSS code.  4 stabilizers, 5 qubits, 1
 *      logical qubit, distance 3. Each stabilizer has weight 4.
 *
 *  ## The [[5, 1, 3]] code
 *
 *  Stabilizer generators (column j is qubit j):
 *
 *      g_1 = X Z Z X I
 *      g_2 = I X Z Z X
 *      g_3 = X I X Z Z
 *      g_4 = Z X I X Z
 *
 *  These are NOT a CSS-style separation into pure-X and pure-Z stabilizers
 *  — they involve mixed Pauli letters per stabilizer. So the 5-qubit
 *  code is a stabilizer code but NOT a CSS code in the strict sense.
 *  We expose it via the abstract `irrep_stabilizer_group_t`.
 *
 *  ## Primary references
 *
 *  - Pastawski-Yoshida-Harlow-Preskill, *Holographic quantum error-
 *    correcting codes: toy models for the bulk/boundary correspondence*,
 *    JHEP 06 (2015) 149 [arXiv:1503.06237].
 *  - Laflamme-Miquel-Paz-Zurek, *Perfect quantum error correcting code*,
 *    Phys. Rev. Lett. 77 (1996) 198 — original [[5, 1, 3]] code.
 *  - Bennett-DiVincenzo-Smolin-Wootters, *Mixed-state entanglement and
 *    quantum error correction*, Phys. Rev. A 54 (1996) 3824 — early
 *    perfect-tensor structure.
 *  - Hayden-Nezami-Qi-Thomas-Walter-Yang, *Holographic duality from
 *    random tensor networks*, JHEP 11 (2016) 009 — generalised HaPPY.
 */
#ifndef IRREP_HAPPY_CODE_H
#define IRREP_HAPPY_CODE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/stabilizer_group.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Build the [[5, 1, 3]] five-qubit "perfect tensor" code as an
 *  abstract stabilizer group.
 *
 *  Generators:
 *    g_1 = X Z Z X I
 *    g_2 = I X Z Z X
 *    g_3 = X I X Z Z
 *    g_4 = Z X I X Z
 *
 *  Caller must `irrep_stabilizer_group_free` when done. */
IRREP_API irrep_status_t
irrep_happy_perfect_tensor_5_1_3(irrep_stabilizer_group_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_HAPPY_CODE_H */
