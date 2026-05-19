/* SPDX-License-Identifier: MIT */
/** @file color_code_3d.h
 *  @brief 3D color codes — concrete cubic-lattice instances.
 *
 *  The 3D color code [Bombín-Martín-Delgado, PRL 98 (2007) 160502]
 *  generalises the 2D color-code family to three spatial dimensions
 *  on a 4-colourable 3-cell complex (the "3-colex"). The defining
 *  features that survive to 3D:
 *
 *    - Each face has both an X- and a Z-stabilizer (the CSS-symmetric
 *      structure of color codes).
 *    - Transversal `T` gate becomes available (via the doubled-color
 *      Reed-Muller equivalence at distance 3); this is the canonical
 *      reason 3D color codes are studied for fault-tolerant logical
 *      non-Clifford gates.
 *    - Logical qubits = `dim H_2(M, ℤ_2)` of the underlying 3-manifold.
 *      For the smallest non-trivial complex — the cube — this is `k = 3`.
 *
 *  ## `[[8, 3, 2]]` cubic 3D color code (Campbell 2019)
 *
 *  Qubits at the 8 vertices of a single cube (binary positions
 *  `(x, y, z) ∈ {0, 1}^3` with linear index `4z + 2y + x`).
 *
 *  Stabilizers chosen so the X-block + Z-block has rank 5 (= n - k):
 *
 *    - 3 X-face stabilizers (one per pair of parallel cube faces):
 *        X on the `z = 1` face = `{4, 5, 6, 7}`  ("top")
 *        X on the `x = 0` face = `{0, 2, 4, 6}`  ("left")
 *        X on the `y = 0` face = `{0, 1, 4, 5}`  ("front")
 *    - 1 X-volume stabilizer: `X ⊗ X ⊗ ... ⊗ X` on all 8 qubits.
 *    - 1 Z-volume stabilizer: `Z ⊗ Z ⊗ ... ⊗ Z` on all 8 qubits.
 *
 *  Counts: `n = 8`, `m_X = 4`, `m_Z = 1`, `k = n - m_X - m_Z = 3`.
 *  Distance `d = 2` (verifiable by brute-force enumeration up to
 *  weight 2; the lowest-weight logical X is a single qubit on an
 *  X-stabilizer's complement, weight 2).
 *
 *  ## Why this matters
 *
 *  - It is **the smallest 3D color code** with non-trivial encoding
 *    (k > 0) and the simplest concrete instance of the Bombín 2007
 *    family.
 *  - It supports a transversal `CCZ` gate via the doubled-color
 *    construction at distance 3; the d=2 instance still demonstrates
 *    the X-Z-symmetric stabilizer geometry that the larger 3D-color
 *    codes inherit.
 *  - It is the natural pedagogical entry point to 3D-color-code
 *    machinery; the full rectified-3-cubic-lattice construction can
 *    be assembled from face-and-cell-list data via the generic
 *    3-colex framework in `<irrep/color_lattice_3d.h>`.
 *
 *  ## Primary references
 *
 *  - Bombín-Martín-Delgado, *Topological computation without
 *    braiding*, Phys. Rev. Lett. 98 (2007) 160502.
 *  - Campbell, *The smallest interesting colour code*,
 *    [arXiv:1907.13278] (2019) — the [[8, 3, 2]] cubic instance.
 *  - Bombín, *Topological subsystem codes*, Phys. Rev. A 81 (2010)
 *    032301 — the gauge / subsystem variant on 3-colex; tetrahedral
 *    instance shipped as `irrep_subsystem_bombin_3d_tetrahedron` in
 *    `<irrep/subsystem_code.h>`.
 */
#ifndef IRREP_COLOR_CODE_3D_H
#define IRREP_COLOR_CODE_3D_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Number of data qubits in the cubic `[[8, 3, 2]]` code. */
#define IRREP_COLOR_3D_CUBE_N           8
/** @brief Number of X-stabilizer generators (3 faces + 1 volume). */
#define IRREP_COLOR_3D_CUBE_N_X_STABS   4
/** @brief Number of Z-stabilizer generators (1 volume). */
#define IRREP_COLOR_3D_CUBE_N_Z_STABS   1
/** @brief Number of logical qubits (k = n - m_X - m_Z = 3). */
#define IRREP_COLOR_3D_CUBE_K           3
/** @brief Code distance — minimum weight of a non-trivial logical. */
#define IRREP_COLOR_3D_CUBE_DISTANCE    2

/** @brief Build the `[[8, 3, 2]]` cubic 3D color code as a CSS code.
 *
 *  Stabilizer ordering:
 *    H_X rows: top face, left face, front face, X-volume.
 *    H_Z rows: Z-volume.
 *
 *  Caller must `irrep_css_code_free(out)` when done. */
IRREP_API irrep_status_t
irrep_color_3d_cube_8_3_2(irrep_css_code_t *out);

/* ====================================================================
 * `[[15, 1, 3]]` Reed-Muller / tetrahedral 3D color code
 *
 * The smallest 3D color code with distance d = 3 (and therefore the
 * smallest Bombín 2007 3-colex instance with the transversal-T
 * property of triorthogonal codes). Qubits sit at the 15 nonzero
 * binary 4-tuples (i.e., the points of PG(3, 2)), indexed `q = 1..15`
 * with `q` interpreted as a 4-bit bitmask.
 *
 * ## Construction (punctured Reed-Muller)
 *
 *   C_X = punctured RM(1, 4)  = [15, 4, 8]   — 4 weight-8 generators
 *   C_Z = even-weight subcode of dual punctured RM(2, 4)
 *                              = [15, 10, 4] — 10 generators
 *
 * Since `C_X ⊂ C_Z^⊥` and `k = n - dim C_X - dim C_Z = 15 - 4 - 10 = 1`,
 * this is a valid CSS code; the standard analysis gives `d = 3`.
 *
 * ### X-stabilizers (4, each weight 8)
 *
 * For each bit `b ∈ {0, 1, 2, 3}`, support = `{q : (q >> b) & 1 == 1}`.
 * Each is the indicator vector of one coordinate bit.
 *
 * ### Z-stabilizers (10, mix of weight 4 and 8)
 *
 * Six "pair" generators (weight 4): for each pair of bits `(b_1, b_2)`,
 * support = `{q : (q >> b_1) & 1 == 1 AND (q >> b_2) & 1 == 1}` —
 * the points where both bits are set.
 *
 * Four "single" generators (weight 8): one per single bit. (These are
 * the same supports as the X-stabilizers, but used here as Z-type
 * generators; the inclusion `C_X ⊂ C_Z` from the construction means
 * the same supports appear in both sides.)
 *
 * ### Why the supports include the 4 weight-8 "single"s on both sides
 *
 * The CSS commutation requirement is `H_X · H_Z^T = 0` over `GF(2)`.
 * Pair × bit overlap = exactly 2 (the points with the pair's two bits
 * set AND the single's bit set, two such points). Bit × bit overlap =
 * 4 (the points with both bits set). Both even ⇒ all pairs commute.
 *
 * ## Logical operators
 *
 * Logical X (and Z) = any weight-3 codeword of `C_Z^⊥` not in `C_X`.
 * Canonical choice: the three colinear points `{1, 2, 3}` in PG(3, 2)
 * (since `001 ⊕ 010 = 011`, the three smallest indices form a
 * projective line) → qubits `{0, 1, 2}` in the 0-indexed layout.
 *
 * ## References
 *
 *  - Bombín-Martín-Delgado 2007 PRL 98 160502 — the topological 3D
 *    color code on a 3-colex.
 *  - Bravyi-Haah 2012 PRA 86 052329 — the [[15, 1, 3]] code's
 *    triorthogonal structure and transversal-T realisation.
 *  - Knill-Laflamme-Zurek 1996 quant-ph/9610011 — original Reed-Muller
 *    CSS construction.
 * ==================================================================== */

/** @brief Number of data qubits in the [[15, 1, 3]] tetrahedral code. */
#define IRREP_COLOR_3D_RM15_N           15
/** @brief Number of X-stabilizer generators (one per coordinate bit). */
#define IRREP_COLOR_3D_RM15_N_X_STABS   4
/** @brief Number of Z-stabilizer generators (4 bits + 6 pairs). */
#define IRREP_COLOR_3D_RM15_N_Z_STABS   10
/** @brief Number of logical qubits. */
#define IRREP_COLOR_3D_RM15_K           1
/** @brief Code distance (minimum weight of any logical operator). */
#define IRREP_COLOR_3D_RM15_DISTANCE    3

/** @brief Build the `[[15, 1, 3]]` Reed-Muller / tetrahedral 3D color
 *  code as a CSS code.
 *
 *  Qubits are indexed `0..14`, corresponding to nonzero binary 4-tuples
 *  `q + 1` (so qubit 0 is bitmask `0001` = point 1 of PG(3, 2), and
 *  qubit 14 is bitmask `1111` = point 15).
 *
 *  Stabilizer ordering:
 *    H_X rows (4, weight 8): one per coordinate bit `b = 0, 1, 2, 3`.
 *    H_Z rows (10): 4 weight-8 "single-bit" generators (same supports
 *      as H_X rows) + 6 weight-4 "pair" generators, lexicographic on
 *      the pair `(b_1, b_2)`.
 *
 *  The CSS code structure does not store logical operators explicitly;
 *  the canonical weight-3 representative `{0, 1, 2}` is documented in
 *  the construction notes above.
 *
 *  Caller must `irrep_css_code_free(out)` when done. */
IRREP_API irrep_status_t
irrep_color_3d_rm_15_1_3(irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_COLOR_CODE_3D_H */
