/* SPDX-License-Identifier: MIT */
/** @file happy_network.h
 *  @brief Multi-tile HaPPY hyperbolic network — foundations for the
 *         Pastawski-Yoshida-Harlow-Preskill bulk/boundary holographic
 *         tensor-network construction.
 *
 *  This module extends `<irrep/happy_code.h>` (which exposes the
 *  [[5, 1, 3]] code as the single-tile depth-1 HaPPY) with the
 *  building-block needed for multi-tile networks: the **6-leg perfect
 *  tensor** form of [[5, 1, 3]], with the bulk (logical) qubit
 *  promoted to an explicit qubit leg so adjacent tiles can be contracted.
 *
 *  ## What "6-leg perfect tensor" means
 *
 *  The [[5, 1, 3]] code's encoding isometry takes 1 bulk qubit to 5
 *  boundary qubits. Promoted to a 6-leg tensor, the bulk qubit becomes
 *  qubit 0 and the 5 boundary qubits become qubits 1..5. The tensor's
 *  stabilizer group has 6 generators on 6 qubits (so `k = 0` — the
 *  tensor itself is a definite quantum state in 6-qubit Hilbert space,
 *  not a code):
 *
 *    Boundary stabilizers (4, lifted from [[5, 1, 3]] with `I` on leg 0):
 *      g_1 = I  X  Z  Z  X  I
 *      g_2 = I  I  X  Z  Z  X
 *      g_3 = I  X  I  X  Z  Z
 *      g_4 = I  Z  X  I  X  Z
 *
 *    Bulk-boundary cross-stabilizers (2):
 *      g_5 = X  X  X  X  X  X   (= bulk-X to boundary X̄)
 *      g_6 = Z  Z  Z  Z  Z  Z   (= bulk-Z to boundary Z̄)
 *
 *  This is the form used by Pastawski-Yoshida-Harlow-Preskill 2015 to
 *  build the {5, 4} hyperbolic tiling network. Two tiles connect by
 *  Bell-pair contraction of one boundary leg from each tile.
 *
 *  ## Depth-2 HaPPY network — structural data
 *
 *  Depth-2 = 1 central tile + 5 layer-1 tiles arranged around it.
 *  The central tile's 5 boundary legs (qubits 1..5) are each
 *  contracted with one boundary leg of one layer-1 tile, leaving:
 *
 *    - 6 bulk legs (1 per tile, free / encoded logical qubits).
 *    - 20 free boundary legs (5 tiles × 4 remaining per tile).
 *
 *  The resulting holographic code is `[[20, 6, ?]]`, with the 6 bulk
 *  qubits as the encoded logical qubits and the 20 free boundary
 *  qubits as the physical realisation. This module exposes both the
 *  **uncontracted joined stabilizer group** of the 6-tile network
 *  (36 qubits, 36 tile generators) AND the **Bell-pair-contracted
 *  isometry** on 26 qubits (`irrep_happy_network_depth2_contracted`),
 *  with bulk-qubit indices reported in the contracted frame.
 *
 *  ## Primary references
 *
 *  - Pastawski-Yoshida-Harlow-Preskill, *Holographic quantum error-
 *    correcting codes: toy models for the bulk/boundary correspondence*,
 *    JHEP 06 (2015) 149 [arXiv:1503.06237].
 *  - Hayden-Nezami-Qi-Thomas-Walter-Yang, *Holographic duality from
 *    random tensor networks*, JHEP 11 (2016) 009 [arXiv:1601.01694] —
 *    contraction algorithm for stabilizer tensor networks.
 */
#ifndef IRREP_HAPPY_NETWORK_H
#define IRREP_HAPPY_NETWORK_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/stabilizer_group.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Build the **6-leg perfect-tensor** form of [[5, 1, 3]] —
 *  the HaPPY tile primitive.
 *
 *  Output: a stabilizer group on 6 qubits with 6 generators. Qubit 0
 *  is the bulk (logical) leg; qubits 1..5 are the boundary legs. The
 *  group has `k = 0` (the tensor is a definite quantum state on 6
 *  qubits, not a code in itself — it becomes a code after contraction
 *  with other tiles).
 *
 *  Generators:
 *    g_1 = I X Z Z X I    (boundary [[5,1,3]] g_1, lifted)
 *    g_2 = I I X Z Z X    (boundary [[5,1,3]] g_2)
 *    g_3 = I X I X Z Z    (boundary [[5,1,3]] g_3)
 *    g_4 = I Z X I X Z    (boundary [[5,1,3]] g_4)
 *    g_5 = X X X X X X    (bulk-X to boundary X̄)
 *    g_6 = Z Z Z Z Z Z    (bulk-Z to boundary Z̄)
 *
 *  Caller must `irrep_stabilizer_group_free` when done. */
IRREP_API irrep_status_t
irrep_happy_perfect_tensor_6leg(irrep_stabilizer_group_t *out);

/** @brief Number of qubits in the depth-2 HaPPY network (1 central tile
 *  + 5 layer-1 tiles, each 6 legs). */
#define IRREP_HAPPY_DEPTH2_N_QUBITS  36

/** @brief Number of tiles in the depth-2 HaPPY network. */
#define IRREP_HAPPY_DEPTH2_N_TILES    6

/** @brief Number of contraction edges in the depth-2 HaPPY network. */
#define IRREP_HAPPY_DEPTH2_N_CONTRACTIONS  5

/** @brief Build the **uncontracted joined stabilizer group** for the
 *  depth-2 HaPPY network.
 *
 *  Layout:
 *    - 6 tiles, each 6 qubits, contiguous: tile `t` occupies qubits
 *      `[6·t, 6·t + 6)`, with the bulk leg at the lowest index.
 *    - Tile 0 = central; tiles 1..5 = layer-1.
 *    - Bulk qubits: 0, 6, 12, 18, 24, 30 (one per tile).
 *
 *  Generators: 6 per tile (`irrep_happy_perfect_tensor_6leg` lifted to
 *  the 36-qubit space) = 36 generators total. All pairwise commuting
 *  by construction (tiles are disjoint).
 *
 *  Contraction edges (NOT added to this uncontracted joined group;
 *  apply them via `irrep_happy_network_depth2_contracted` to get the
 *  [[20, 6, ?]] boundary isometry): the central tile's boundary
 *  qubits 1..5 are paired with layer-1 tile k's boundary qubit 1
 *  (= absolute index `6·(k+1) + 1` for `k = 0..4`). Adding the
 *  corresponding Bell-pair stabilizers `(X_a X_b, Z_a Z_b)` and
 *  performing symplectic elimination via the Bell-pair contraction
 *  primitive yields the [[20, 6, ?]] boundary code.
 *
 *  @param[out] out  Stabilizer group with `n = 36`, `n_generators = 36`.
 *                   Caller must `irrep_stabilizer_group_free` when done.
 *  @param[out] contraction_pairs  Optional, length `2 * 5 = 10`. If
 *                   non-NULL, written as `(a_0, b_0, a_1, b_1, ...)`
 *                   with each (a_i, b_i) the qubit indices of the i-th
 *                   contraction. Pass NULL to skip. */
IRREP_API irrep_status_t
irrep_happy_network_depth2(irrep_stabilizer_group_t *out,
                           int *contraction_pairs);

/** @brief Number of qubits after Bell-pair contraction of all 5 edges. */
#define IRREP_HAPPY_DEPTH2_N_CONTRACTED_QUBITS  26

/** @brief Build the **Bell-pair-contracted** depth-2 HaPPY network.
 *
 *  Internally:
 *    1. Build the 36-qubit joined group with `irrep_happy_network_depth2`.
 *    2. For each of the 5 contraction edges, apply
 *       `irrep_stabilizer_contract_bell` in sequence. Each call projects
 *       onto the `|Φ+⟩` eigenspace of `(X_a X_b, Z_a Z_b)` for the
 *       pair, then traces out qubits `a` and `b`, re-indexing all
 *       higher-numbered qubits.
 *    3. Track index-shift bookkeeping so later contractions reference
 *       the right qubits.
 *
 *  Output is a 26-qubit stabilizer group encoding the post-contracted
 *  HaPPY state. It has rank 26 (a pure stabilizer state, `k = 0`),
 *  but is naturally viewed as a `[[20, 6, ?]]` ENCODING ISOMETRY by
 *  splitting the 26 qubits into:
 *    - 6 bulk qubits at indices `bulk_qubits_out[0..5]` — the encoded
 *      logical qubits of the HaPPY code.
 *    - 20 boundary qubits at all other indices — the physical realisation.
 *
 *  Of the 26 stabilizer generators, the 14 acting trivially on every
 *  bulk leg are the stabilizers of the boundary code; the remaining 12
 *  cross-couple bulk and boundary, encoding the logical operators
 *  `(X̄_i, Z̄_i)` for `i = 0..5` modulo stabilizers.
 *
 *  Bulk-qubit contracted-frame indices (computed by subtracting from
 *  each original bulk index `6 t` the number of removed qubits with
 *  smaller original index): `{0, 1, 6, 11, 16, 21}`.
 *
 *  @param[out] out                 26-qubit stabilizer group; caller
 *                                  must `irrep_stabilizer_group_free`.
 *  @param[out] bulk_qubits_out     Optional length-6 array; if non-NULL,
 *                                  written with bulk-qubit indices in
 *                                  the contracted frame. */
IRREP_API irrep_status_t
irrep_happy_network_depth2_contracted(irrep_stabilizer_group_t *out,
                                      int *bulk_qubits_out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_HAPPY_NETWORK_H */
