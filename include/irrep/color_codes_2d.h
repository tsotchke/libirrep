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

/** @brief Number of data qubits in the `[[17, 1, 5]]` 4.8.8 color code. */
#define IRREP_COLOR_488_17_1_5_N         17
/** @brief Number of face stabilizers (5 octagons + 3 squares). */
#define IRREP_COLOR_488_17_1_5_N_FACES    8
/** @brief Code distance for the 4.8.8 `[[17, 1, 5]]` instance. */
#define IRREP_COLOR_488_17_1_5_DISTANCE   5

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

/** @brief Build the canonical logical X̄ operator for the `[[19, 1, 5]]`
 *  hex color code.
 *
 *  `L̄_X = X` on the 5 data qubits along the bottom edge of the
 *  triangular patch (y = 0 row) — qubits `{D0, D1, D2, D3, D4}` =
 *  `{0, 1, 2, 3, 4}`. Weight 5 (= code distance).
 *
 *  Commutes with every stabilizer (every face shares 0 or 2 qubits
 *  with this support) and anti-commutes with `irrep_color_hex_19_1_5_logical_Z`
 *  on the corner qubit `D0 = 0`. The Pauli `out` is allocated by this
 *  function; caller must `irrep_pauli_free` when done. */
IRREP_API irrep_status_t
irrep_color_hex_19_1_5_logical_X(irrep_pauli_t *out);

/** @brief Build the canonical logical Z̄ operator for the `[[19, 1, 5]]`
 *  hex color code.
 *
 *  `L̄_Z = Z` on the 5 data qubits along the left edge of the triangular
 *  patch (the `x = 2y` diagonal from `D0 = 0` at the bottom-left corner
 *  to `D18 = 18` at the apex) — qubits `{D0, D9, D12, D17, D18}` =
 *  `{0, 9, 12, 17, 18}`. Weight 5.
 *
 *  Commutes with every stabilizer (every face shares 0 or 2 qubits
 *  with this support) and anti-commutes with `irrep_color_hex_19_1_5_logical_X`
 *  on the shared corner qubit `D0 = 0` (the bottom-left vertex of the
 *  triangle, the unique qubit lying on both the bottom and left edges). */
IRREP_API irrep_status_t
irrep_color_hex_19_1_5_logical_Z(irrep_pauli_t *out);

/* ====================================================================
 * `[[17, 1, 5]]` square-octagon (4.8.8) triangular color code.
 *
 * Triangular patch of the 4.8.8 (square-octagon, also called union-jack)
 * lattice with code parameters `[n, k, d] = [d²/2 + d - 1/2, 1, d]`. At
 * `d = 5`: 17 data qubits, 1 logical qubit, distance 5. Construction
 * follows the MQT QECC `SquareOctagonColorCode` algorithm (Derks, TUM,
 * github.com/munich-quantum-toolkit/qecc).
 *
 * 8 face stabilizers split into:
 *   - 1 weight-8 central octagon (fully internal, `qubits ∈
 *     {2, 3, 6, 7, 10, 11, 13, 14}`).
 *   - 4 weight-4 boundary-truncated octagons (corner / edge octagons
 *     with 4 of their 8 vertices outside the triangular patch).
 *   - 3 weight-4 squares.
 *
 * Sum of weights = `1·8 + 7·4 = 36 = 2·n + 2`. Per-qubit incidence:
 * 3 corner qubits in 1 face, 9 edge qubits in 2 faces, 5 interior
 * qubits in 3 faces.
 *
 * Each face appears as both an X- and a Z-stabilizer (color-code
 * symmetry); 16 generators on 17 qubits → k = 1 logical qubit.
 *
 * Primary references:
 *   - Bombín-Martín-Delgado, PRL 97 (2006) 180501 — original 4.8.8
 *     color-code framework.
 *   - Derks et al. (MQT QECC), `SquareOctagonColorCode` —
 *     parametric d-instance generator whose face supports we reuse.
 *   - Error Correction Zoo, "Square-octagon (4.8.8) color code",
 *     errorcorrectionzoo.org/c/488_color.
 * ==================================================================== */

/** @brief Build the `[[17, 1, 5]]` 4.8.8 (square-octagon) color code.
 *
 *  Allocates `out` with `n = 17` qubits, `m_X = m_Z = 8` stabilizer
 *  generators. CSS-orthogonal by construction. Caller must
 *  `irrep_css_code_free(out)` when done. */
IRREP_API irrep_status_t
irrep_color_488_17_1_5(irrep_css_code_t *out);

/** @brief Canonical logical X̄ for the 4.8.8 `[[17, 1, 5]]` code.
 *
 *  X-string on the 5 data qubits along the bottom edge of the
 *  triangular patch (y = 1 row) — qubits `{0, 1, 2, 3, 4}`. Weight 5. */
IRREP_API irrep_status_t
irrep_color_488_17_1_5_logical_X(irrep_pauli_t *out);

/** @brief Canonical logical Z̄ for the 4.8.8 `[[17, 1, 5]]` code.
 *
 *  Z-string on the 5 data qubits along the right edge of the
 *  triangular patch (the `x + y = 13` diagonal from `D4 = (12, 1)`
 *  to `D15 = (6, 7)`) — qubits `{4, 8, 11, 14, 15}`. Weight 5.
 *
 *  Anti-commutes with `irrep_color_488_17_1_5_logical_X` on their
 *  unique shared corner qubit `D4 = 4`. */
IRREP_API irrep_status_t
irrep_color_488_17_1_5_logical_Z(irrep_pauli_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_COLOR_CODES_2D_H */
