/* SPDX-License-Identifier: MIT */
/** @file toric_code.h
 *  @brief Kitaev toric-code stabilizer primitives on a 2D periodic lattice.
 *
 *  The toric code is the canonical topological-stabilizer code: qubits
 *  live on the edges of an L_x × L_y square lattice, with two stabilizer
 *  types:
 *    - **Vertex (star) operators** A_v = ⊗_{e ∋ v} X_e, one per vertex.
 *    - **Plaquette operators** B_p = ⊗_{e ∈ ∂p} Z_e, one per face.
 *
 *  All stabilizers commute pairwise. The ground space is 4-fold
 *  degenerate (encoding 2 logical qubits via the homology classes of
 *  non-contractible loops on the torus).
 *
 *  This module provides:
 *    - **Edge enumeration**: 2·L_x·L_y edges (horizontal + vertical).
 *    - **Stabilizer enumeration**: A_v indexed by vertex, B_p by face.
 *    - **Edge incidence tables**: which 4 edges hit each vertex / each face.
 *    - **Code parameter check**: [[2 L_x L_y, 2, min(L_x, L_y)]] for a
 *      torus, with structural witnesses.
 *
 *  Sister to `lattice.h` (which only handles sites). We track *edges*
 *  here, the natural location for stabilizer-code qubits.
 *
 *  Primary references:
 *    - Kitaev, Ann. Phys. 303 (2003) 2 (arXiv:quant-ph/9707021)
 *    - Bombin & Martin-Delgado, Phys. Rev. A 76 (2007) 012305
 */
#ifndef IRREP_TORIC_CODE_H
#define IRREP_TORIC_CODE_H

#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Edge orientation on the square lattice. */
typedef enum {
    IRREP_TORIC_EDGE_HORIZONTAL = 0, /**< Edge from (x, y) to (x+1, y). */
    IRREP_TORIC_EDGE_VERTICAL   = 1, /**< Edge from (x, y) to (x, y+1). */
} irrep_toric_edge_orient_t;

/** @brief Edge encoded as a triple (orient, x, y).
 *
 *  The packed integer index is
 *  `e = (orient * Lx * Ly) + (y * Lx + x)`,
 *  giving 0 ≤ e < 2·Lx·Ly. */
typedef struct {
    int orient; /**< IRREP_TORIC_EDGE_HORIZONTAL or _VERTICAL. */
    int x;      /**< Cell index x ∈ [0, Lx). */
    int y;      /**< Cell index y ∈ [0, Ly). */
} irrep_toric_edge_t;

/** @brief Toric-code lattice parameters. */
typedef struct {
    int Lx;            /**< Linear size in x. */
    int Ly;            /**< Linear size in y. */
    int n_qubits;      /**< Number of qubits = 2 · Lx · Ly. */
    int n_vertices;    /**< Number of vertex stabilizers = Lx · Ly. */
    int n_plaquettes;  /**< Number of plaquette stabilizers = Lx · Ly. */
} irrep_toric_params_t;

/** @brief Initialise toric-code parameters. */
IRREP_API irrep_status_t
irrep_toric_init(irrep_toric_params_t *out, int Lx, int Ly);

/** @brief Pack edge (orient, x, y) into a linear index in [0, 2·Lx·Ly).
 *
 *  Returns -1 if any component is out of range; otherwise returns the index. */
IRREP_API int
irrep_toric_edge_index(const irrep_toric_params_t *p,
                       int orient, int x, int y);

/** @brief Unpack a linear edge index into (orient, x, y). */
IRREP_API irrep_status_t
irrep_toric_edge_unpack(const irrep_toric_params_t *p,
                        int edge_index,
                        irrep_toric_edge_t *out);

/** @brief List the 4 edges incident to vertex (vx, vy).
 *
 *  These are: horizontal edge at (vx-1, vy), horizontal at (vx, vy),
 *  vertical at (vx, vy-1), vertical at (vx, vy). All wrapped modulo
 *  (Lx, Ly) for PBC.
 *
 *  @param[out] edges  Array of 4 edge indices (caller-allocated).
 *  @return            IRREP_STATUS_OK on success. */
IRREP_API irrep_status_t
irrep_toric_vertex_edges(const irrep_toric_params_t *p,
                         int vx, int vy,
                         int edges[4]);

/** @brief List the 4 edges bounding face (px, py).
 *
 *  Plaquette p at (px, py) has corners at (px, py), (px+1, py),
 *  (px, py+1), (px+1, py+1). Its 4 boundary edges are:
 *  bottom horizontal, top horizontal, left vertical, right vertical.
 *
 *  @param[out] edges  Array of 4 edge indices (caller-allocated). */
IRREP_API irrep_status_t
irrep_toric_plaquette_edges(const irrep_toric_params_t *p,
                            int px, int py,
                            int edges[4]);

/** @brief Stabilizer commutativity witness for the toric code.
 *
 *  Counts the number of edges shared between two stabilizers (vertex,
 *  plaquette), each represented by its 4-edge boundary. If the count is
 *  even (0 or 2), the stabilizers commute; if odd (1 or 3), they
 *  anti-commute.
 *
 *  Toric-code property: every (A_v, B_p) pair shares an even number of
 *  edges, hence A_v B_p = B_p A_v.
 *
 *  @return  Number of shared edges (0..4). */
IRREP_API int
irrep_toric_shared_edges(const int edges_a[4], const int edges_b[4]);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_TORIC_CODE_H */
