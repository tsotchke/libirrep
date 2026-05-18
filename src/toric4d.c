/* SPDX-License-Identifier: MIT */
/** @file toric4d.c
 *  @brief 4D toric code construction (qubits on 2-cells of T⁴). */
#include <irrep/toric4d.h>
#include <irrep/types.h>

#include <stddef.h>

/* ====================================================================
 * Indexing on the 4-torus.
 *
 * Vertex (x, y, z, w) ∈ [0, Lx) × ... × [0, Lw):
 *      v = ((w·Lz + z)·Ly + y)·Lx + x.
 *
 * Edge (v, d) for d ∈ {0, 1, 2, 3} = orientation:
 *      e = d · V + v.
 *
 * Face (v, orient) for orient ∈ {0..5}, ordered (i, j) lex:
 *      0 = (0,1) [xy]   1 = (0,2) [xz]   2 = (0,3) [xw]
 *      3 = (1,2) [yz]   4 = (1,3) [yw]   5 = (2,3) [zw]
 *
 * 3-cube (v, orient3) for orient3 ∈ {0..3}, ordered as "drop d":
 *      0 = {1,2,3} [yzw]   1 = {0,2,3} [xzw]
 *      2 = {0,1,3} [xyw]   3 = {0,1,2} [xyz]
 *
 * (i.e. drop axis 0..3 in order; convenient for the CSS-orthogonality
 * proof because cube `o3 = d` is dual to edge `d` under the 4D Hodge
 * star.)
 * ==================================================================== */

static const int kFaceLow[6]  = { 0, 0, 0, 1, 1, 2 };
static const int kFaceHigh[6] = { 1, 2, 3, 2, 3, 3 };

/* "Drop axis d" → 3-cube spans the remaining axes (sorted ascending). */
static const int kCubeAxes[4][3] = {
    { 1, 2, 3 },  /* drop 0 */
    { 0, 2, 3 },  /* drop 1 */
    { 0, 1, 3 },  /* drop 2 */
    { 0, 1, 2 },  /* drop 3 */
};

static inline int wrap(int v, int L) {
    int r = v % L;
    if (r < 0) r += L;
    return r;
}

static inline int vertex_idx_at(const irrep_toric4d_params_t *p,
                                int x, int y, int z, int w) {
    return ((w * p->Lz + z) * p->Ly + y) * p->Lx + x;
}

/* Step a vertex along axis `d` by `+delta` (with wraparound). */
static inline int vertex_step(const irrep_toric4d_params_t *p,
                              int x, int y, int z, int w, int d, int delta) {
    int xs = x, ys = y, zs = z, ws = w;
    if      (d == 0) xs = wrap(x + delta, p->Lx);
    else if (d == 1) ys = wrap(y + delta, p->Ly);
    else if (d == 2) zs = wrap(z + delta, p->Lz);
    else             ws = wrap(w + delta, p->Lw);
    return vertex_idx_at(p, xs, ys, zs, ws);
}

/* Face orientation index, given (i, j) with i < j. */
static int face_orient(int i, int j) {
    for (int o = 0; o < 6; ++o) {
        if (kFaceLow[o] == i && kFaceHigh[o] == j) return o;
    }
    return -1;
}

static inline int face_idx(const irrep_toric4d_params_t *p,
                           int v, int orient) {
    return orient * p->V + v;
}

static inline int edge_idx(const irrep_toric4d_params_t *p,
                           int v, int orient) {
    return orient * p->V + v;
}

static inline int cube_idx(const irrep_toric4d_params_t *p,
                           int v, int orient3) {
    return orient3 * p->V + v;
}

irrep_status_t
irrep_toric4d_init(irrep_toric4d_params_t *out,
                   int Lx, int Ly, int Lz, int Lw)
{
    if (out == NULL || Lx < 1 || Ly < 1 || Lz < 1 || Lw < 1)
        return IRREP_ERR_INVALID_ARG;
    out->Lx = Lx; out->Ly = Ly; out->Lz = Lz; out->Lw = Lw;
    out->V = Lx * Ly * Lz * Lw;
    out->n_qubits  = 6 * out->V;
    out->n_X_stabs = 4 * out->V;
    out->n_Z_stabs = 4 * out->V;
    return IRREP_OK;
}

irrep_status_t
irrep_toric4d_build(const irrep_toric4d_params_t *p, irrep_css_code_t *out)
{
    if (p == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;

    irrep_status_t s = irrep_css_code_new(out, p->n_qubits,
                                          p->n_X_stabs, p->n_Z_stabs);
    if (s != IRREP_OK) return s;

    /* -------- X-stabilizers: edge → 6 incident faces --------
     *
     * An edge of direction `d` at vertex `v_base` connects `v_base` to
     * `v_base + e_d`. The 3 face orientations containing `d` are
     * `(d, j)` for each `j ≠ d`. For each such orientation, the edge
     * bounds 2 faces: the face anchored at `v_base` and the face
     * anchored at `v_base - e_j` (one face on each side of the edge
     * in the j-direction).
     *
     * We index faces by their "lower-corner" vertex: face (i, j)
     * anchored at v_base has corners
     *      v_base, v_base + e_i, v_base + e_j, v_base + e_i + e_j.
     */
    for (int w = 0; w < p->Lw; ++w)
    for (int z = 0; z < p->Lz; ++z)
    for (int y = 0; y < p->Ly; ++y)
    for (int x = 0; x < p->Lx; ++x) {
        int v = vertex_idx_at(p, x, y, z, w);
        for (int d = 0; d < 4; ++d) {
            int row = edge_idx(p, v, d);
            for (int j = 0; j < 4; ++j) {
                if (j == d) continue;
                int i = d < j ? d : j;
                int k = d < j ? j : d;
                int orient = face_orient(i, k);
                /* Face at v_base. */
                irrep_parity_matrix_set(&out->H_X, row, face_idx(p, v, orient));
                /* Face at v_base - e_j. */
                int v_back = vertex_step(p, x, y, z, w, j, -1);
                irrep_parity_matrix_set(&out->H_X, row, face_idx(p, v_back, orient));
            }
        }
    }

    /* -------- Z-stabilizers: 3-cube → 6 bounding faces --------
     *
     * A 3-cube of orientation `o3 = "drop d"` (axes {a, b, c} = {0..3}
     * \ {d}) has vertices `v_base + ε_a e_a + ε_b e_b + ε_c e_c` for
     * ε ∈ {0, 1}^3. It is bounded by 6 faces — the C(3,2) = 3 pairs
     * of axes within {a, b, c}, each with 2 faces (ε = 0 and ε = 1
     * on the third axis):
     *
     *      Pair (a, b): faces (orient(a,b)) at v_base and v_base + e_c.
     *      Pair (a, c): faces (orient(a,c)) at v_base and v_base + e_b.
     *      Pair (b, c): faces (orient(b,c)) at v_base and v_base + e_a.
     */
    for (int w = 0; w < p->Lw; ++w)
    for (int z = 0; z < p->Lz; ++z)
    for (int y = 0; y < p->Ly; ++y)
    for (int x = 0; x < p->Lx; ++x) {
        int v = vertex_idx_at(p, x, y, z, w);
        for (int o3 = 0; o3 < 4; ++o3) {
            int row = cube_idx(p, v, o3);
            int a = kCubeAxes[o3][0];
            int b = kCubeAxes[o3][1];
            int c = kCubeAxes[o3][2];
            /* Pair (a, b), opposite axis c. */
            int orient_ab = face_orient(a, b);
            irrep_parity_matrix_set(&out->H_Z, row, face_idx(p, v, orient_ab));
            int v_c = vertex_step(p, x, y, z, w, c, +1);
            irrep_parity_matrix_set(&out->H_Z, row, face_idx(p, v_c, orient_ab));
            /* Pair (a, c), opposite axis b. */
            int orient_ac = face_orient(a, c);
            irrep_parity_matrix_set(&out->H_Z, row, face_idx(p, v, orient_ac));
            int v_b = vertex_step(p, x, y, z, w, b, +1);
            irrep_parity_matrix_set(&out->H_Z, row, face_idx(p, v_b, orient_ac));
            /* Pair (b, c), opposite axis a. */
            int orient_bc = face_orient(b, c);
            irrep_parity_matrix_set(&out->H_Z, row, face_idx(p, v, orient_bc));
            int v_a = vertex_step(p, x, y, z, w, a, +1);
            irrep_parity_matrix_set(&out->H_Z, row, face_idx(p, v_a, orient_bc));
        }
    }

    return IRREP_OK;
}
