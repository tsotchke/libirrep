/* SPDX-License-Identifier: MIT */
/** @file xcube_code.c
 *  @brief X-cube fracton code (Vijay-Haah-Fu 2016) construction. */
#include <irrep/xcube_code.h>
#include <irrep/types.h>

#include <stddef.h>

static inline int wrap(int v, int L) {
    int r = v % L;
    if (r < 0) r += L;
    return r;
}

static inline int edge_idx(const irrep_xcube_params_t *p,
                           int x, int y, int z, int orient)
{
    int per_cell = p->Lx * p->Ly * p->Lz;
    return orient * per_cell + (z * p->Lx * p->Ly + y * p->Lx + x);
}

static inline int linear_idx(const irrep_xcube_params_t *p,
                             int x, int y, int z)
{
    return z * p->Lx * p->Ly + y * p->Lx + x;
}

irrep_status_t
irrep_xcube_init(irrep_xcube_params_t *out, int Lx, int Ly, int Lz)
{
    if (out == NULL || Lx <= 0 || Ly <= 0 || Lz <= 0) return IRREP_ERR_INVALID_ARG;
    out->Lx = Lx; out->Ly = Ly; out->Lz = Lz;
    int N = Lx * Ly * Lz;
    out->n_qubits = 3 * N;
    out->n_vertex_stabs = 3 * N;
    out->n_cube_stabs = N;
    return IRREP_OK;
}

irrep_status_t
irrep_xcube_build(const irrep_xcube_params_t *p, irrep_css_code_t *out)
{
    if (p == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    int Lx = p->Lx, Ly = p->Ly, Lz = p->Lz;
    int N = Lx * Ly * Lz;

    irrep_status_t s = irrep_css_code_new(out, p->n_qubits,
                                          p->n_vertex_stabs, p->n_cube_stabs);
    if (s != IRREP_OK) return s;

    /* Planar vertex stabilizers (X-type, weight 4). Three flavours per vertex.
     * Block layout: [xy block of N rows | xz block of N rows | yz block of N rows].
     */
    for (int vz = 0; vz < Lz; ++vz) {
        for (int vy = 0; vy < Ly; ++vy) {
            for (int vx = 0; vx < Lx; ++vx) {
                int xm = wrap(vx - 1, Lx);
                int ym = wrap(vy - 1, Ly);
                int zm = wrap(vz - 1, Lz);
                int v = linear_idx(p, vx, vy, vz);

                /* xy-plane: 2 x-edges + 2 y-edges. */
                int row = v;
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, xm, vy, vz, 0));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, vz, 0));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, ym, vz, 1));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, vz, 1));

                /* xz-plane: 2 x-edges + 2 z-edges. */
                row = N + v;
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, xm, vy, vz, 0));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, vz, 0));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, zm, 2));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, vz, 2));

                /* yz-plane: 2 y-edges + 2 z-edges. */
                row = 2 * N + v;
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, ym, vz, 1));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, vz, 1));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, zm, 2));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, vz, 2));
            }
        }
    }

    /* Cube (Z-type) stabilizers, weight 12. The cube at (cx, cy, cz) has 12 edges:
     *   x-direction (4): (cx, cy+b, cz+c, 0) for b, c ∈ {0,1}
     *   y-direction (4): (cx+a, cy, cz+c, 1) for a, c ∈ {0,1}
     *   z-direction (4): (cx+a, cy+b, cz, 2) for a, b ∈ {0,1}
     */
    for (int cz = 0; cz < Lz; ++cz) {
        for (int cy = 0; cy < Ly; ++cy) {
            for (int cx = 0; cx < Lx; ++cx) {
                int row = linear_idx(p, cx, cy, cz);
                int xp1 = wrap(cx + 1, Lx);
                int yp1 = wrap(cy + 1, Ly);
                int zp1 = wrap(cz + 1, Lz);
                /* x-edges */
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, cx, cy,  cz,  0));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, cx, yp1, cz,  0));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, cx, cy,  zp1, 0));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, cx, yp1, zp1, 0));
                /* y-edges */
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, cx,  cy, cz,  1));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, xp1, cy, cz,  1));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, cx,  cy, zp1, 1));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, xp1, cy, zp1, 1));
                /* z-edges */
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, cx,  cy,  cz, 2));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, xp1, cy,  cz, 2));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, cx,  yp1, cz, 2));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, xp1, yp1, cz, 2));
            }
        }
    }
    return IRREP_OK;
}
