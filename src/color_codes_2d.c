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

/* ====================================================================
 * Logical operators for `[[19, 1, 5]]` hex color code.
 *
 * Triangle corners (in the (x, y) coordinate system): D0 = (0, 0)
 * at bottom-left, D4 = (24, 0) at bottom-right, D18 = (12, 6) at
 * apex. Each canonical logical is a weight-5 string of identical
 * Pauli letters along one edge of the triangle:
 *
 *   bottom edge   (y = 0):           {D0, D1, D2, D3, D4}
 *   left  edge    (x = 2y diagonal): {D0, D9, D12, D17, D18}
 *   right edge    (x = 24 - 2y):     {D4, D8, D14, D16, D18}
 *
 * Any two distinct edges share exactly one corner qubit; an X-string
 * on one edge and a Z-string on a different edge anti-commute there.
 *
 * Exposed canonical choice: L̄_X = bottom edge, L̄_Z = left edge. They
 * share corner D0 only.
 * ==================================================================== */

static const int kLogicalX_qubits[5] = { 0, 1, 2, 3, 4 };
static const int kLogicalZ_qubits[5] = { 0, 9, 12, 17, 18 };

irrep_status_t
irrep_color_hex_19_1_5_logical_X(irrep_pauli_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_pauli_new(out, IRREP_COLOR_HEX_19_1_5_N);
    if (s != IRREP_OK) return s;
    for (int i = 0; i < 5; ++i) {
        irrep_pauli_set(out, kLogicalX_qubits[i], IRREP_PAULI_LETTER_X);
    }
    return IRREP_OK;
}

irrep_status_t
irrep_color_hex_19_1_5_logical_Z(irrep_pauli_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_pauli_new(out, IRREP_COLOR_HEX_19_1_5_N);
    if (s != IRREP_OK) return s;
    for (int i = 0; i < 5; ++i) {
        irrep_pauli_set(out, kLogicalZ_qubits[i], IRREP_PAULI_LETTER_Z);
    }
    return IRREP_OK;
}

/* ====================================================================
 * `[[17, 1, 5]]` 4.8.8 (square-octagon) triangular color code.
 *
 * Data-qubit layout (17 qubits at integer (x, y) coordinates produced
 * by the MQT-QECC `SquareOctagonColorCode(distance=5)` algorithm,
 * indexed row-major y-then-x):
 *
 *   D0=(0,1)   D1=(2,1)   D2=(6,1)   D3=(8,1)   D4=(12,1)
 *   D5=(3,2)   D6=(5,2)   D7=(9,2)   D8=(11,2)
 *   D9=(3,4)   D10=(5,4)  D11=(9,4)  D12=(11,4)
 *   D13=(6,5)  D14=(8,5)
 *   D15=(6,7)  D16=(8,7)
 *
 * Ancilla positions (8 faces): 5 octagons at (4,0), (10,0), (1,3),
 * (7,3), (10,6) plus 3 squares at (4,3), (10,3), (7,6). Each
 * octagon face connects to 8 neighbour offsets (±2,±1)/(±1,±2);
 * each square to 4 diagonal offsets (±1,±1). Face supports:
 *
 *   F0 (oct  4,0) = {1, 2, 5, 6}                      weight 4
 *   F1 (oct 10,0) = {3, 4, 7, 8}                      weight 4
 *   F2 (oct  1,3) = {0, 1, 5, 9}                      weight 4
 *   F3 (oct  7,3) = {2, 3, 6, 7, 10, 11, 13, 14}      weight 8
 *   F4 (oct 10,6) = {11, 12, 14, 16}                  weight 4
 *   F5 (sq   4,3) = {5, 6, 9, 10}                     weight 4
 *   F6 (sq  10,3) = {7, 8, 11, 12}                    weight 4
 *   F7 (sq   7,6) = {13, 14, 15, 16}                  weight 4
 *
 * All 28 face pairs intersect in 0 or 2 qubits (CSS-orthogonal).
 *
 * Triangle corners: D0 = (0, 1) bottom-left, D4 = (12, 1) bottom-right,
 * D15 = (6, 7) apex. Each is in exactly 1 face.
 * ==================================================================== */

static const int k488Faces[8][9] = {
    {  1,  2,  5,  6, -1, -1, -1, -1, -1 },                /* F0 oct boundary */
    {  3,  4,  7,  8, -1, -1, -1, -1, -1 },                /* F1 oct boundary */
    {  0,  1,  5,  9, -1, -1, -1, -1, -1 },                /* F2 oct boundary */
    {  2,  3,  6,  7, 10, 11, 13, 14, -1 },                /* F3 oct central  */
    { 11, 12, 14, 16, -1, -1, -1, -1, -1 },                /* F4 oct boundary */
    {  5,  6,  9, 10, -1, -1, -1, -1, -1 },                /* F5 square       */
    {  7,  8, 11, 12, -1, -1, -1, -1, -1 },                /* F6 square       */
    { 13, 14, 15, 16, -1, -1, -1, -1, -1 },                /* F7 square       */
};

irrep_status_t
irrep_color_488_17_1_5(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;

    const int n_qubits = IRREP_COLOR_488_17_1_5_N;
    const int n_faces  = IRREP_COLOR_488_17_1_5_N_FACES;

    irrep_status_t s = irrep_css_code_new(out, n_qubits, n_faces, n_faces);
    if (s != IRREP_OK) return s;

    for (int f = 0; f < n_faces; ++f) {
        for (int i = 0; i < 9; ++i) {
            int q = k488Faces[f][i];
            if (q < 0) break;
            irrep_parity_matrix_set(&out->H_X, f, q);
            irrep_parity_matrix_set(&out->H_Z, f, q);
        }
    }

    return IRREP_OK;
}

/* Logical operators for `[[17, 1, 5]]` 4.8.8.
 *
 * L̄_X = bottom-edge X-string from corner D0 to corner D4:
 *      {D0, D1, D2, D3, D4} = {0, 1, 2, 3, 4}, weight 5.
 *
 * L̄_Z = right-edge Z-string from corner D4 to corner D15 along the
 *      x + y = 13 diagonal: {D4, D8, D11, D14, D15} = {4, 8, 11, 14, 15},
 *      weight 5.
 *
 * Both verified by hand: every face has 0 or 2 overlap qubits with each
 * support. They share corner D4 only → symplectic inner product = 1
 * → anti-commute. */

static const int k488LogicalX_qubits[5] = { 0, 1, 2, 3, 4 };
static const int k488LogicalZ_qubits[5] = { 4, 8, 11, 14, 15 };

irrep_status_t
irrep_color_488_17_1_5_logical_X(irrep_pauli_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_pauli_new(out, IRREP_COLOR_488_17_1_5_N);
    if (s != IRREP_OK) return s;
    for (int i = 0; i < 5; ++i) {
        irrep_pauli_set(out, k488LogicalX_qubits[i], IRREP_PAULI_LETTER_X);
    }
    return IRREP_OK;
}

irrep_status_t
irrep_color_488_17_1_5_logical_Z(irrep_pauli_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_pauli_new(out, IRREP_COLOR_488_17_1_5_N);
    if (s != IRREP_OK) return s;
    for (int i = 0; i < 5; ++i) {
        irrep_pauli_set(out, k488LogicalZ_qubits[i], IRREP_PAULI_LETTER_Z);
    }
    return IRREP_OK;
}
