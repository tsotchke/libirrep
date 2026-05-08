/* SPDX-License-Identifier: MIT */
/** @file toric3d.c
 *  @brief 3D toric code construction. */
#include <irrep/toric3d.h>
#include <irrep/types.h>

#include <stddef.h>

static inline int wrap(int v, int L) {
    int r = v % L;
    if (r < 0) r += L;
    return r;
}

static inline int edge_idx(const irrep_toric3d_params_t *p,
                           int x, int y, int z, int orient)
{
    int per_cell = p->Lx * p->Ly * p->Lz;
    return orient * per_cell + (z * p->Lx * p->Ly + y * p->Lx + x);
}

static inline int vertex_idx(const irrep_toric3d_params_t *p,
                             int x, int y, int z)
{
    return z * p->Lx * p->Ly + y * p->Lx + x;
}

irrep_status_t
irrep_toric3d_init(irrep_toric3d_params_t *out, int Lx, int Ly, int Lz)
{
    if (out == NULL || Lx <= 0 || Ly <= 0 || Lz <= 0) return IRREP_ERR_INVALID_ARG;
    out->Lx = Lx; out->Ly = Ly; out->Lz = Lz;
    int N = Lx * Ly * Lz;
    out->n_qubits = 3 * N;
    out->n_vertex_stabs = N;
    out->n_face_stabs = 3 * N;
    return IRREP_OK;
}

irrep_status_t
irrep_toric3d_build(const irrep_toric3d_params_t *p, irrep_css_code_t *out)
{
    if (p == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    int Lx = p->Lx, Ly = p->Ly, Lz = p->Lz;

    irrep_status_t s = irrep_css_code_new(out, p->n_qubits,
                                          p->n_vertex_stabs, p->n_face_stabs);
    if (s != IRREP_OK) return s;

    /* Vertex stabilizers (X-type, weight 6).
     * The 6 edges incident to vertex (vx, vy, vz):
     *   x-direction: (vx-1, vy, vz, 0)  and  (vx, vy, vz, 0)
     *   y-direction: (vx, vy-1, vz, 1)  and  (vx, vy, vz, 1)
     *   z-direction: (vx, vy, vz-1, 2)  and  (vx, vy, vz, 2)
     */
    for (int vz = 0; vz < Lz; ++vz) {
        for (int vy = 0; vy < Ly; ++vy) {
            for (int vx = 0; vx < Lx; ++vx) {
                int row = vertex_idx(p, vx, vy, vz);
                int xm = wrap(vx - 1, Lx);
                int ym = wrap(vy - 1, Ly);
                int zm = wrap(vz - 1, Lz);
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, xm, vy, vz, 0));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, vz, 0));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, ym, vz, 1));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, vz, 1));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, zm, 2));
                irrep_parity_matrix_set(&out->H_X, row, edge_idx(p, vx, vy, vz, 2));
            }
        }
    }

    /* Plaquette stabilizers (Z-type, weight 4).
     * Three orientations:
     *   xy-plaquette at (px, py, pz): edges
     *     (px, py, pz, 0)        — bottom horizontal
     *     (px+1, py, pz, 1)      — right vertical
     *     (px, py+1, pz, 0)      — top horizontal
     *     (px, py, pz, 1)        — left vertical
     *   yz-plaquette at (px, py, pz): edges (orient 1 and 2)
     *     (px, py, pz, 1), (px, py+1, pz, 2), (px, py, pz+1, 1), (px, py, pz, 2)
     *   xz-plaquette at (px, py, pz): edges (orient 0 and 2)
     *     (px, py, pz, 0), (px+1, py, pz, 2), (px, py, pz+1, 0), (px, py, pz, 2)
     */
    int N = Lx * Ly * Lz;
    /* xy block: rows [0, N) */
    for (int pz = 0; pz < Lz; ++pz) {
        for (int py = 0; py < Ly; ++py) {
            for (int px = 0; px < Lx; ++px) {
                int row = vertex_idx(p, px, py, pz);
                int xp1 = wrap(px + 1, Lx);
                int yp1 = wrap(py + 1, Ly);
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, px,  py,  pz, 0));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, xp1, py,  pz, 1));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, px,  yp1, pz, 0));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, px,  py,  pz, 1));
            }
        }
    }
    /* yz block: rows [N, 2N) */
    for (int pz = 0; pz < Lz; ++pz) {
        for (int py = 0; py < Ly; ++py) {
            for (int px = 0; px < Lx; ++px) {
                int row = N + vertex_idx(p, px, py, pz);
                int yp1 = wrap(py + 1, Ly);
                int zp1 = wrap(pz + 1, Lz);
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, px, py,  pz,  1));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, px, yp1, pz,  2));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, px, py,  zp1, 1));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, px, py,  pz,  2));
            }
        }
    }
    /* xz block: rows [2N, 3N) */
    for (int pz = 0; pz < Lz; ++pz) {
        for (int py = 0; py < Ly; ++py) {
            for (int px = 0; px < Lx; ++px) {
                int row = 2 * N + vertex_idx(p, px, py, pz);
                int xp1 = wrap(px + 1, Lx);
                int zp1 = wrap(pz + 1, Lz);
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, px,  py, pz,  0));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, xp1, py, pz,  2));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, px,  py, zp1, 0));
                irrep_parity_matrix_set(&out->H_Z, row, edge_idx(p, px,  py, pz,  2));
            }
        }
    }
    return IRREP_OK;
}
