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

/* TODO(R1): Triangular hex [[19, 1, 5]] color code face list.
 *
 * Completion is mechanical once the LAR 2011 Fig. 3 / Kubica thesis
 * Fig. 2.1.5 vertex-numbering + face-membership are transcribed.
 *
 * Acceptance criteria (each verifiable in seconds):
 *   - n_qubits = 19, n_faces = 9.
 *   - Face weights: exactly 3 weight-3 corners + 6 weight-6 interior.
 *   - irrep_css_code_verify(&out) returns IRREP_OK (i.e., every pair of
 *     faces shares an even number of qubits).
 *   - Each qubit appears in either 1 face (3 corner qubits), 2 faces
 *     (6 edge qubits), or 3 faces (10 interior qubits). Self-consistency:
 *     3·3 + 6·6 = 45 = 1·3 + 2·6 + 3·10 = 45. ✓
 *   - irrep_qec_distance_brute(&g, 5) returns 5 (brute-force distance
 *     enumerator handles n=19 at weight 5 in well under a minute).
 *
 * Drop-in template (fill in the 9 face_qubits[][] arrays):
 *
 *     static const int face_R_corner[3] = { ?, ?, ? };
 *     static const int face_G_corner[3] = { ?, ?, ? };
 *     static const int face_B_corner[3] = { ?, ?, ? };
 *     static const int face_interior_0[6] = { ?, ?, ?, ?, ?, ? };
 *     ... face_interior_1 ... face_interior_5 ...
 *
 *     static const int sizes[9] = { 3, 3, 3, 6, 6, 6, 6, 6, 6 };
 *     static const int *qubits[9] = { face_R_corner, face_G_corner,
 *         face_B_corner, face_interior_0, ..., face_interior_5 };
 *     irrep_color_lattice_t L = { .n_qubits = 19, .n_faces = 9,
 *         .face_sizes = sizes, .face_qubits = qubits, .face_color = NULL };
 *     return irrep_generic_color_build(&L, out);
 *
 * Multiple prior attempts at reconstructing this face list from first
 * principles failed CSS-orthogonality (pairs of candidate faces meeting
 * in odd counts); the published figure is the load-bearing input.
 */
irrep_status_t
irrep_color_hex_19_1_5(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    return IRREP_ERR_NOT_IMPLEMENTED;
}

/* TODO(R2): Triangular [[17, 1, 5]] color code on the 4.8.8 (square-octagon)
 * lattice. Completion criteria parallel to R1:
 *   - n_qubits = 17, n_faces = 8 (1 weight-8 octagon + 4 weight-4 squares
 *     + 3 weight-3 corner-truncated boundaries).
 *   - Sum of weights = 8 + 16 + 9 = 33 = 1·3 + 2·6 + 3·8 (3 corners,
 *     6 edges, 8 interior).
 *   - Verifier: irrep_css_code_verify + brute-force distance returns 5.
 *
 * Drop-in template:
 *   static const int face_octagon[8] = { ?, ... };
 *   static const int face_sq_0[4] = { ?, ?, ?, ? };  ... face_sq_3 ...
 *   static const int face_R_corner[3] = { ?, ?, ? };  ... G, B ...
 *
 * Load-bearing input: Pogorelov 2024 (Quantinuum) Fig. 2 vertex/face labels.
 */
irrep_status_t
irrep_color_488_17_1_5(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    return IRREP_ERR_NOT_IMPLEMENTED;
}
