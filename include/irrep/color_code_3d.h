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
 *    machinery before the full rectified-3-cubic-lattice construction
 *    (which is documented in `docs/qec_research_roadmap.md` as a
 *    research-track item requiring a generic 3-colex generator).
 *
 *  ## Primary references
 *
 *  - Bombín-Martín-Delgado, *Topological computation without
 *    braiding*, Phys. Rev. Lett. 98 (2007) 160502.
 *  - Campbell, *The smallest interesting colour code*,
 *    [arXiv:1907.13278] (2019) — the [[8, 3, 2]] cubic instance.
 *  - Bombín, *Topological subsystem codes*, Phys. Rev. A 81 (2010)
 *    032301 — the gauge / subsystem variant on 3-colex, still a
 *    research-track item.
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

#ifdef __cplusplus
}
#endif

#endif /* IRREP_COLOR_CODE_3D_H */
