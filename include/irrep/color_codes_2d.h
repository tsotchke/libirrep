/* SPDX-License-Identifier: MIT */
/** @file color_codes_2d.h
 *  @brief Concrete 2D color-code instances on hex and 4.8.8 lattices.
 *
 *  Provides ready-built constructors for triangular color codes at
 *  distance 5 and beyond. Each instance is wrapped in
 *  `irrep_generic_color_build` internally; see that header for the
 *  framework.
 *
 *  ## Instances
 *
 *    - `[[19, 1, 5]]` hexagonal (6.6.6) triangular color code.
 *      19 qubits, 9 hexagonal faces.
 *
 *    - `[[17, 1, 5]]` square-octagon (4.8.8) triangular color code.
 *      17 qubits, 8 faces (1 weight-8 octagon, 4 weight-4 squares,
 *      3 weight-3 corner truncations).
 *
 *  ## References
 *
 *  - Landahl-Anderson-Rice, *Fault-tolerant quantum computing with color
 *    codes*, arXiv:1108.5738 (2011) Fig. 3 (`[[19, 1, 5]]`).
 *  - Pogorelov et al. (Quantinuum), *Experimental fault-tolerant code-
 *    switching*, arXiv:2403.13732 (2024) Fig. 2 (`[[17, 1, 5]]`).
 *  - Bombín-Martín-Delgado, PRL 97 (2006) 180501 — original 2D color
 *    code framework.
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

/** @brief Build the [[19, 1, 5]] hexagonal triangular color code.
 *
 *  Layout: 19 qubits at vertices of a triangular patch of the
 *  honeycomb lattice with side L=2 (distance d = 2L+1 = 5). The 9
 *  hexagonal faces partition into:
 *    - 3 weight-3 corner-truncated hexagons (at the 3 triangle corners)
 *    - 3 weight-5 edge hexagons (on the 3 triangle edges)
 *    - 3 weight-6 fully-internal hexagons
 *
 *  Sum of face weights = 3·3 + 3·5 + 3·6 = 9 + 15 + 18 = 42 = 2·n_qubits + 4
 *  (since the average vertex degree on the triangular hex patch is
 *  slightly above 2).
 *
 *  Caller must `irrep_css_code_free` when done. */
IRREP_API irrep_status_t
irrep_color_hex_19_1_5(irrep_css_code_t *out);

/** @brief Build the [[17, 1, 5]] square-octagon triangular color code.
 *
 *  Layout: 17 qubits on a 4.8.8 (square-octagon) lattice patch.
 *  8 faces:
 *    - 1 central weight-8 octagon
 *    - 4 weight-4 squares around the octagon
 *    - 3 weight-3 corner-truncated faces (one per triangle corner)
 *
 *  Sum of face weights = 8 + 4·4 + 3·3 = 8 + 16 + 9 = 33 = 2·17 - 1
 *  (slightly under-degree by 1 due to truncated corners).
 *
 *  Caller must `irrep_css_code_free` when done. */
IRREP_API irrep_status_t
irrep_color_488_17_1_5(irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_COLOR_CODES_2D_H */
