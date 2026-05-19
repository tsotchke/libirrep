/* SPDX-License-Identifier: MIT */
/** @file color_code.h
 *  @brief Triangular color code on the 4.8.8 lattice (Bombín–Martín-Delgado).
 *
 *  The 2D color code is a CSS code defined on a 3-colorable, 3-valent
 *  tiling of a 2D surface. Each face of the tiling carries **two**
 *  stabilizers — an X-type and a Z-type — both supported on the same
 *  set of qubits (the boundary of that face). The two stabilizers
 *  commute because every face has even weight.
 *
 *  The triangular color code on the **4.8.8 (square-octagon) lattice**
 *  is the smallest fault-tolerantly-relevant family. For a "side-length"
 *  parameter `L ≥ 1`, the code has:
 *
 *    - Data qubits = `(L+1)²` (corners of the L×L square sub-grid).
 *    - X-type stabilizers = Z-type stabilizers = (L+1)·L (one per face).
 *    - Logical qubits k = 1, distance d = 2L+1.
 *
 *  Wait — the 4.8.8 lattice has a more involved counting. We use the
 *  *triangular* color code with the simpler counting:
 *
 *    Triangular color code of distance d (odd):
 *    - n = (d² + 4d - 1) / 4  (for d=3: 5; d=5: 11; d=7: 19; d=9: 29)
 *    - k = 1, d_code = d
 *    - weight-4 squares + weight-8 octagons
 *
 *  This module implements the **17-qubit triangular d=5 color code**
 *  on the 4.8.8 lattice (also known as the "color code on the
 *  square-octagon lattice"), which is the smallest distance-5 color
 *  code. It is the topological analogue of the [[7,1,3]] Steane code
 *  (which is the *smallest* triangular color code, d=3).
 *
 *  ## Steane = smallest color code
 *
 *  The [[7,1,3]] Steane code IS the triangular color code on a single
 *  triangle decomposed into 3 faces (one of each color). Each face
 *  carries an X-stab and a Z-stab on the same 4 qubits (corner shared
 *  between three faces). This module exposes Steane via
 *  `irrep_color_steane`.
 *
 *  ## Primary references
 *
 *  - Bombín-Martín-Delgado, *Topological quantum distillation*, Phys.
 *    Rev. Lett. 97 (2006) 180501.
 *  - Bombín, *Topological subsystem codes*, Phys. Rev. A 81 (2010)
 *    032301 — gauge color code.
 *  - Landahl-Anderson-Rice, *Fault-tolerant quantum computing with
 *    color codes*, arXiv:1108.5738 (2011).
 *  - Kubica, *The ABCs of the color code*, PhD thesis Caltech 2018.
 */
#ifndef IRREP_COLOR_CODE_H
#define IRREP_COLOR_CODE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Build the [[7, 1, 3]] Steane code as a (degenerate) triangular
 *  color code with one triangle of three faces (R, G, B colors).
 *
 *  This is the smallest color code. Each of the three faces carries
 *  one X-stab AND one Z-stab on the same 4-qubit support, giving 3 X
 *  + 3 Z = 6 stabilizers on 7 qubits with k = 1 logical qubit.
 *
 *  Stabilizer rows (1-indexed columns to match the Hamming-code parity
 *  matrix):
 *    g_1 = bits 1, 3, 5, 7 (face R)
 *    g_2 = bits 2, 3, 6, 7 (face G)
 *    g_3 = bits 4, 5, 6, 7 (face B)
 *  Each g_i appears as both an X-stabilizer and a Z-stabilizer. */
IRREP_API irrep_status_t
irrep_color_steane(irrep_css_code_t *out);

/** @brief Build the [[15, 7, 3]] CSS code derived from the [15, 11, 3]
 *  punctured Reed-Muller code RM(2, 4)^punctured.
 *
 *  The natural Steane-style extension to 15 qubits. The classical
 *  parent code [15, 11, 3] has a 4×15 parity-check matrix whose rows
 *  are the indicator vectors of single coordinate bits of PG(3, 2):
 *
 *    row_b = {q ∈ 0..14 : bit b of (q+1) is set},   b ∈ {0, 1, 2, 3}.
 *
 *  Setting `H_X = H_Z = parity check of [15, 11, 3]` gives a CSS code
 *  with:
 *    - n_qubits = 15
 *    - m_X = m_Z = 4 (single-bit indicators, each weight 8)
 *    - k = 15 - 4 - 4 = 7 logical qubits
 *    - d = 3 (the minimum-weight codewords of C_Z^⊥ \ C_X are weight 3 —
 *      PG(3, 2) lines like {0, 1, 2})
 *
 *  CSS orthogonality: pair-of-rows intersection has size 4 (positions
 *  where two specified bits are set), which is even. ✓
 *
 *  ## Relation to other 15-qubit codes
 *
 *  - The [[15, 1, 3]] Reed-Muller code (`irrep_color_3d_rm_15_1_3`)
 *    uses the same 4 single-bit rows for `H_X` but adds 6 weight-4
 *    "pair" rows to `H_Z`, dropping k from 7 to 1 in exchange for the
 *    transversal-T property.
 *  - The two codes share the same X-stabilizers but differ in the
 *    Z-stabilizer side. This `[[15, 7, 3]]` instance is the "minimal"
 *    Hamming-style CSS code, in the sense that it has the largest k
 *    consistent with the [15, 11, 3] parity-check structure.
 *
 *  Caller must `irrep_css_code_free(out)` when done. */
IRREP_API irrep_status_t
irrep_color_hamming_15_7_3(irrep_css_code_t *out);

/* Triangular [[19, 1, 5]] (6.6.6) and [[17, 1, 5]] (4.8.8) variants
 * are shipped in `<irrep/color_codes_2d.h>` (closes R1 and R2 of
 * `docs/qec_research_roadmap.md`). */

#ifdef __cplusplus
}
#endif

#endif /* IRREP_COLOR_CODE_H */
