/* SPDX-License-Identifier: MIT */
/** @file color_codes_2d.c
 *  @brief Concrete 2D color-code instances on hex and 4.8.8 lattices.
 *
 *  This file is a research-track scaffold. The two concrete instances
 *  ([[19, 1, 5]] hex and [[17, 1, 5]] 488) require the explicit face
 *  list of a triangular patch of the honeycomb / square-octagon
 *  lattice, which is most cleanly expressed by reading a published
 *  figure (Landahl-Anderson-Rice 2011 Fig. 3 or Pogorelov 2024 Fig. 2)
 *  and translating vertex labels to qubit indices.
 *
 *  Until that translation is done, both constructors return
 *  `IRREP_ERR_NOT_IMPLEMENTED` so callers receive an explicit signal
 *  rather than silent miscompilation. Internal counts are documented
 *  here so that completing the constructor is mechanical once the
 *  face list is in hand:
 *
 *  ## [[19, 1, 5]] hex layout (LAR 2011 convention)
 *
 *  - 3 corner-truncated weight-3 faces, one per triangle corner.
 *  - 6 interior weight-6 hexagonal faces.
 *  - Vertex partition: 3 corners + 6 edge vertices (2 per triangle side)
 *    + 10 interior = 19. Each corner is in 1 face; each edge in 2 faces;
 *    each interior in 3 faces.
 *  - Self-consistency: 3·3 + 6·6 = 45 = sum of `(faces containing v)`
 *    over all v = 3·10 + 2·6 + 1·3 = 45. ✓
 *
 *  ## [[17, 1, 5]] 488 layout (Pogorelov 2024 / Quantinuum convention)
 *
 *  - 1 central weight-8 octagon.
 *  - 4 weight-4 squares, each sharing one edge (2 vertices) with the
 *    octagon.
 *  - 3 weight-3 corner-truncated boundary faces.
 *  - Total 8 faces, 17 qubits, sum of weights = 8 + 16 + 9 = 33.
 *
 *  Both are validated by the framework module `generic_color_code` once
 *  populated.
 */
#include <irrep/color_codes_2d.h>
#include <irrep/types.h>

#include <stddef.h>

irrep_status_t
irrep_color_hex_19_1_5(irrep_css_code_t *out)
{
    /* Research-track: requires explicit translation of LAR 2011 Fig. 3
     * vertex / face labels to a face_qubits[][] table. The framework is
     * in place via `irrep_generic_color_build`; only the table is
     * outstanding.
     *
     * Acceptance criteria once filled:
     *   - n_qubits = 19, n_faces = 9.
     *   - 3 weight-3 + 6 weight-6 faces.
     *   - irrep_css_code_verify returns IRREP_OK.
     *   - Brute-force distance up to weight 5 returns 5 (n=19 is in
     *     range for our brute-force enumerator).
     */
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    return IRREP_ERR_NOT_IMPLEMENTED;
}

irrep_status_t
irrep_color_488_17_1_5(irrep_css_code_t *out)
{
    /* Research-track: requires explicit translation of Pogorelov 2024
     * Fig. 2 vertex / face labels to a face_qubits[][] table. */
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    return IRREP_ERR_NOT_IMPLEMENTED;
}
