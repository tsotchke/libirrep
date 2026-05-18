/* SPDX-License-Identifier: MIT */
/** @file color_codes_2d.h
 *  @brief Concrete 2D color-code instances on the 6.6.6 honeycomb lattice.
 *
 *  Provides a ready-built constructor for the `[[19, 1, 5]]` triangular
 *  color code on the honeycomb (6.6.6) lattice — the second member of
 *  the Bombín-Martín-Delgado family `[[(3d² + 1)/4, 1, d]]` for odd `d`,
 *  the first being the `[[7, 1, 3]]` Steane code (see
 *  `<irrep/color_code.h>`).
 *
 *  ## `[[19, 1, 5]]` hexagonal layout
 *
 *  19 data qubits laid out on a triangular patch of the 6.6.6 honeycomb.
 *  Qubit indices follow the canonical color-code-stim integer grid: each
 *  data qubit lives at an even-coordinate `(x, y)` position with
 *  `y ∈ [0, 6]`, indexed row-major top-to-bottom.
 *
 *  9 hexagonal-face stabilizers split into:
 *    - 3 weight-6 fully-internal hexagons (interior of the triangle).
 *    - 6 weight-4 boundary-truncated hexagons (the 3 triangle sides each
 *      carry 2 boundary faces).
 *
 *  Sum of weights = `3·6 + 6·4 = 18 + 24 = 42 = 2·n - 4`. The 42 vertex-
 *  face incidences distribute as 3 corner vertices × 1 face + 9 edge
 *  vertices × 2 faces + 7 interior vertices × 3 faces.
 *
 *  Each face appears as both an X-stabilizer (row of `H_X`) and a
 *  Z-stabilizer (row of `H_Z`) — the defining property of CSS color
 *  codes — giving 18 stabilizer generators on 19 data qubits, hence
 *  `k = n - rank(H_X) - rank(H_Z) = 19 - 9 - 9 = 1` logical qubit.
 *
 *  The construction is verified against published references and tested
 *  in `tests/test_color_codes_2d.c` for CSS orthogonality, stabilizer
 *  commutativity, and code distance.
 *
 *  ## Primary references
 *
 *  - Bombín-Martín-Delgado, *Topological quantum distillation*, Phys.
 *    Rev. Lett. 97 (2006) 180501 [arXiv:quant-ph/0605138] — the
 *    original `[[(3d²+1)/4, 1, d]]` family.
 *  - Landahl-Anderson-Rice, *Fault-tolerant quantum computing with
 *    color codes*, arXiv:1108.5738 (2011) — `[[19, 1, 5]]` numerics
 *    and decoder benchmarks (Fig. 3).
 *  - Kubica, *The ABCs of the color code*, Caltech PhD thesis (2018)
 *    Section 2.1 — pedagogical exposition with explicit lattice
 *    diagrams.
 *  - Lee, *color-code-stim* (github.com/seokhyung-lee/color-code-stim) —
 *    open-source Python implementation; the qubit-position and face-
 *    incidence convention adopted here matches its triangular-patch
 *    geometry exactly.
 */
#ifndef IRREP_COLOR_CODES_2D_H
#define IRREP_COLOR_CODES_2D_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Number of data qubits in the `[[19, 1, 5]]` hex color code. */
#define IRREP_COLOR_HEX_19_1_5_N        19
/** @brief Number of hexagonal-face stabilizers (= X-stabs = Z-stabs). */
#define IRREP_COLOR_HEX_19_1_5_N_FACES   9
/** @brief Code distance — minimum weight of a non-trivial logical. */
#define IRREP_COLOR_HEX_19_1_5_DISTANCE  5

/** @brief Build the `[[19, 1, 5]]` triangular hexagonal color code.
 *
 *  Allocates `out` with `n = 19` qubits, `m_X = m_Z = 9` stabilizer
 *  generators (each face contributes both an X- and a Z-row). The face
 *  layout follows the convention documented at the top of this header
 *  and is verified CSS-orthogonal at construction time. Caller must
 *  `irrep_css_code_free(out)` when done.
 *
 *  @return IRREP_OK on success; IRREP_ERR_INVALID_ARG if `out` is NULL;
 *          IRREP_ERR_NOMEM if allocation fails. */
IRREP_API irrep_status_t
irrep_color_hex_19_1_5(irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_COLOR_CODES_2D_H */
