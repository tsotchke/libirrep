/* SPDX-License-Identifier: MIT */
/** @file color_codes_2d.c
 *  @brief Implementation of the `[[19, 1, 5]]` triangular hexagonal color
 *         code on the 6.6.6 honeycomb lattice.
 *
 *  ## Geometric construction
 *
 *  Following the integer-coordinate convention of color-code-stim
 *  (github.com/seokhyung-lee/color-code-stim), the triangular patch
 *  at distance `d = 5` is generated on the grid `y ∈ [0, L]` with
 *  `L = round(3·(d-1)/2) = 6` by enumerating all positions
 *  `(x, y)` with `2y ≤ x ≤ 4L - 2y` and `x ≡ 2y (mod 4)`. This gives
 *  `7 + 6 + 5 + 4 + 3 + 2 + 1 = 28` lattice positions, of which the
 *  9 ancilla positions are picked out by the 3-coloring rule
 *
 *      round((x/2 - y) / 2)  mod 3  ==  anc_qubit_pos[y mod 3]
 *
 *  with `anc_qubit_pos = (2, 0, 1)`. The remaining 19 positions are
 *  data qubits, indexed row-major top-to-top-then-left-to-right:
 *
 *      D0 = (0,0)    D1 = (4,0)    D2 = (12,0)   D3 = (16,0)   D4 = (24,0)
 *      D5 = (6,1)    D6 = (10,1)   D7 = (18,1)   D8 = (22,1)
 *      D9 = (4,2)    D10 = (12,2)  D11 = (16,2)
 *      D12 = (6,3)   D13 = (10,3)  D14 = (18,3)
 *      D15 = (12,4)  D16 = (16,4)
 *      D17 = (10,5)
 *      D18 = (12,6)
 *
 *  Each face (ancilla position) connects to its 6 hexagonal-neighbor
 *  positions at offsets `(-2, 1), (2, 1), (4, 0), (2, -1), (-2, -1),
 *  (-4, 0)`. Offsets that land outside the patch (boundary-truncated)
 *  drop, giving 6 weight-4 boundary faces and 3 weight-6 interior
 *  faces. All 36 face pairs intersect in 0 or 2 data qubits — i.e. the
 *  CSS orthogonality condition `H_X · H_Zᵀ = 0 (mod 2)` is satisfied
 *  by construction and verified by `irrep_css_code_verify` in the
 *  paired test.
 */
#include <irrep/color_codes_2d.h>
#include <irrep/types.h>

#include <stddef.h>

/* Hard-coded face support — derived from the geometric construction
 * above. Each row lists the data-qubit indices on one face, terminated
 * by -1.
 *
 * Layout (face | weight | type):
 *   F0  4  boundary  (top-left corner region)
 *   F1  4  boundary  (top-right corner region)
 *   F2  4  boundary  (upper-left side)
 *   F3  6  interior  (upper-middle hexagon)
 *   F4  6  interior  (middle hexagon)
 *   F5  4  boundary  (right side)
 *   F6  6  interior  (lower-middle hexagon)
 *   F7  4  boundary  (lower-left side)
 *   F8  4  boundary  (bottom corner region) */
static const int kHexFaces[9][7] = {
    {  1,  2,  5,  6, -1, -1, -1 },   /* F0 weight 4 */
    {  3,  4,  7,  8, -1, -1, -1 },   /* F1 weight 4 */
    {  0,  1,  5,  9, -1, -1, -1 },   /* F2 weight 4 */
    {  2,  3,  6,  7, 10, 11, -1 },   /* F3 weight 6 */
    {  5,  6,  9, 10, 12, 13, -1 },   /* F4 weight 6 */
    {  7,  8, 11, 14, -1, -1, -1 },   /* F5 weight 4 */
    { 10, 11, 13, 14, 15, 16, -1 },   /* F6 weight 6 */
    { 12, 13, 15, 17, -1, -1, -1 },   /* F7 weight 4 */
    { 15, 16, 17, 18, -1, -1, -1 },   /* F8 weight 4 */
};

irrep_status_t
irrep_color_hex_19_1_5(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;

    const int n_qubits = IRREP_COLOR_HEX_19_1_5_N;
    const int n_faces  = IRREP_COLOR_HEX_19_1_5_N_FACES;

    irrep_status_t s = irrep_css_code_new(out, n_qubits, n_faces, n_faces);
    if (s != IRREP_OK) return s;

    /* Color codes have H_X = H_Z = face-support matrix. */
    for (int f = 0; f < n_faces; ++f) {
        for (int i = 0; i < 7; ++i) {
            int q = kHexFaces[f][i];
            if (q < 0) break;
            irrep_parity_matrix_set(&out->H_X, f, q);
            irrep_parity_matrix_set(&out->H_Z, f, q);
        }
    }

    return IRREP_OK;
}
