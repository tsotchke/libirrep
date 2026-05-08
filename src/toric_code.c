/* SPDX-License-Identifier: MIT */
/** @file toric_code.c
 *  @brief Implementation of toric-code stabilizer primitives.
 */
#include <irrep/toric_code.h>
#include <irrep/types.h>

#include <stdbool.h>

irrep_status_t
irrep_toric_init(irrep_toric_params_t *out, int Lx, int Ly)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    if (Lx <= 0 || Ly <= 0) return IRREP_ERR_INVALID_ARG;

    out->Lx = Lx;
    out->Ly = Ly;
    out->n_qubits = 2 * Lx * Ly;
    out->n_vertices = Lx * Ly;
    out->n_plaquettes = Lx * Ly;
    return IRREP_OK;
}

int
irrep_toric_edge_index(const irrep_toric_params_t *p,
                       int orient, int x, int y)
{
    if (p == NULL) return -1;
    if (orient != IRREP_TORIC_EDGE_HORIZONTAL &&
        orient != IRREP_TORIC_EDGE_VERTICAL) return -1;
    if (x < 0 || x >= p->Lx) return -1;
    if (y < 0 || y >= p->Ly) return -1;
    return orient * (p->Lx * p->Ly) + (y * p->Lx + x);
}

irrep_status_t
irrep_toric_edge_unpack(const irrep_toric_params_t *p,
                        int edge_index,
                        irrep_toric_edge_t *out)
{
    if (p == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    int n_edges_per_orient = p->Lx * p->Ly;
    if (edge_index < 0 || edge_index >= 2 * n_edges_per_orient)
        return IRREP_ERR_INVALID_ARG;
    out->orient = edge_index / n_edges_per_orient;
    int rem = edge_index % n_edges_per_orient;
    out->y = rem / p->Lx;
    out->x = rem % p->Lx;
    return IRREP_OK;
}

/* PBC-wrap an integer mod L. */
static int wrap_pbc(int v, int L)
{
    while (v < 0) v += L;
    while (v >= L) v -= L;
    return v;
}

irrep_status_t
irrep_toric_vertex_edges(const irrep_toric_params_t *p,
                         int vx, int vy,
                         int edges[4])
{
    if (p == NULL || edges == NULL) return IRREP_ERR_INVALID_ARG;
    if (vx < 0 || vx >= p->Lx) return IRREP_ERR_INVALID_ARG;
    if (vy < 0 || vy >= p->Ly) return IRREP_ERR_INVALID_ARG;

    int xm1 = wrap_pbc(vx - 1, p->Lx);
    int ym1 = wrap_pbc(vy - 1, p->Ly);

    /* 4 edges meeting at vertex (vx, vy):
     *  - West horizontal: (xm1, vy)
     *  - East horizontal: (vx, vy)
     *  - South vertical:  (vx, ym1)
     *  - North vertical:  (vx, vy)
     */
    edges[0] = irrep_toric_edge_index(p, IRREP_TORIC_EDGE_HORIZONTAL, xm1, vy);
    edges[1] = irrep_toric_edge_index(p, IRREP_TORIC_EDGE_HORIZONTAL, vx, vy);
    edges[2] = irrep_toric_edge_index(p, IRREP_TORIC_EDGE_VERTICAL, vx, ym1);
    edges[3] = irrep_toric_edge_index(p, IRREP_TORIC_EDGE_VERTICAL, vx, vy);

    for (int i = 0; i < 4; ++i) {
        if (edges[i] < 0) return IRREP_ERR_PRECONDITION;
    }
    return IRREP_OK;
}

irrep_status_t
irrep_toric_plaquette_edges(const irrep_toric_params_t *p,
                            int px, int py,
                            int edges[4])
{
    if (p == NULL || edges == NULL) return IRREP_ERR_INVALID_ARG;
    if (px < 0 || px >= p->Lx) return IRREP_ERR_INVALID_ARG;
    if (py < 0 || py >= p->Ly) return IRREP_ERR_INVALID_ARG;

    int xp1 = wrap_pbc(px + 1, p->Lx);
    int yp1 = wrap_pbc(py + 1, p->Ly);

    /* 4 edges bordering plaquette (px, py):
     *  - Bottom horizontal: (px, py)
     *  - Top horizontal:    (px, yp1)
     *  - Left vertical:     (px, py)
     *  - Right vertical:    (xp1, py)
     */
    edges[0] = irrep_toric_edge_index(p, IRREP_TORIC_EDGE_HORIZONTAL, px, py);
    edges[1] = irrep_toric_edge_index(p, IRREP_TORIC_EDGE_HORIZONTAL, px, yp1);
    edges[2] = irrep_toric_edge_index(p, IRREP_TORIC_EDGE_VERTICAL, px, py);
    edges[3] = irrep_toric_edge_index(p, IRREP_TORIC_EDGE_VERTICAL, xp1, py);

    for (int i = 0; i < 4; ++i) {
        if (edges[i] < 0) return IRREP_ERR_PRECONDITION;
    }
    return IRREP_OK;
}

int
irrep_toric_shared_edges(const int edges_a[4], const int edges_b[4])
{
    int count = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (edges_a[i] == edges_b[j]) {
                count++;
                break;  /* avoid double-counting if duplicates */
            }
        }
    }
    return count;
}
