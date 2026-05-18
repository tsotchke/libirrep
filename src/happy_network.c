/* SPDX-License-Identifier: MIT */
/** @file happy_network.c
 *  @brief Multi-tile HaPPY hyperbolic network primitives. */
#include <irrep/happy_network.h>
#include <irrep/stabilizer_group.h>
#include <irrep/types.h>

#include <stddef.h>
#include <stdlib.h>

/* [[5, 1, 3]] boundary generators, lifted to 6-qubit space with `I` on
 * leg 0. From `<irrep/happy_code.h>`:
 *   g_1 = X Z Z X I  (positions 1..5)
 *   g_2 = I X Z Z X
 *   g_3 = X I X Z Z
 *   g_4 = Z X I X Z
 *
 * The 6-leg form lifts them with leg 0 = I, then appends two cross-
 * stabilizers relating the bulk leg to the boundary X̄ and Z̄ logicals
 * (both X^⊗5 and Z^⊗5 of [[5,1,3]]):
 *   g_5 = X X X X X X
 *   g_6 = Z Z Z Z Z Z
 */
static const irrep_pauli_letter_t kPerfectTensor6Leg[6][6] = {
    /*                  leg: 0  1  2  3  4  5  */
    /* g_1 = lifted [[5,1,3]] g_1 = (XZZXI) */
    { IRREP_PAULI_LETTER_I,
      IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_Z,
      IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_X,
      IRREP_PAULI_LETTER_I },
    /* g_2 = (IXZZX) */
    { IRREP_PAULI_LETTER_I,
      IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_X,
      IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z,
      IRREP_PAULI_LETTER_X },
    /* g_3 = (XIXZZ) */
    { IRREP_PAULI_LETTER_I,
      IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_I,
      IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_Z,
      IRREP_PAULI_LETTER_Z },
    /* g_4 = (ZXIXZ) */
    { IRREP_PAULI_LETTER_I,
      IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_X,
      IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_X,
      IRREP_PAULI_LETTER_Z },
    /* g_5 = bulk-X to boundary X̄ */
    { IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_X,
      IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_X,
      IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_X },
    /* g_6 = bulk-Z to boundary Z̄ */
    { IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z,
      IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z,
      IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z },
};

irrep_status_t
irrep_happy_perfect_tensor_6leg(irrep_stabilizer_group_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_stabilizer_group_new(out, /*n=*/6, /*n_gens=*/6);
    if (s != IRREP_OK) return s;
    for (int r = 0; r < 6; ++r) {
        for (int q = 0; q < 6; ++q) {
            if (kPerfectTensor6Leg[r][q] != IRREP_PAULI_LETTER_I) {
                irrep_pauli_set(&out->gens[r], q, kPerfectTensor6Leg[r][q]);
            }
        }
    }
    return IRREP_OK;
}

irrep_status_t
irrep_happy_network_depth2(irrep_stabilizer_group_t *out,
                           int *contraction_pairs)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    const int n_qubits = IRREP_HAPPY_DEPTH2_N_QUBITS;       /* 36 */
    const int n_tiles  = IRREP_HAPPY_DEPTH2_N_TILES;        /* 6 */
    const int per_tile = 6;
    const int gens_per_tile = 6;
    const int n_gens   = n_tiles * gens_per_tile;           /* 36 */

    irrep_status_t s = irrep_stabilizer_group_new(out, n_qubits, n_gens);
    if (s != IRREP_OK) return s;

    /* Lift the 6-leg perfect tensor stabilizers onto each of the 6
     * tile slots, with tile `t` occupying qubits [6t, 6t + 6). */
    for (int t = 0; t < n_tiles; ++t) {
        int base = t * per_tile;
        for (int r = 0; r < gens_per_tile; ++r) {
            int row = t * gens_per_tile + r;
            for (int q = 0; q < per_tile; ++q) {
                if (kPerfectTensor6Leg[r][q] != IRREP_PAULI_LETTER_I) {
                    irrep_pauli_set(&out->gens[row], base + q,
                                    kPerfectTensor6Leg[r][q]);
                }
            }
        }
    }

    /* Contraction edges: central tile's boundary qubits 1..5 paired
     * with layer-1 tile k's boundary qubit 1 (= absolute index
     * 6·(k+1) + 1, for k = 0..4). */
    if (contraction_pairs != NULL) {
        for (int k = 0; k < IRREP_HAPPY_DEPTH2_N_CONTRACTIONS; ++k) {
            contraction_pairs[2 * k + 0] = 0 * per_tile + (k + 1);
            contraction_pairs[2 * k + 1] = (k + 1) * per_tile + 1;
        }
    }
    return IRREP_OK;
}

/* Map an original qubit index to its index in the contracted frame,
 * given the list of original qubit indices that have been removed by
 * previous contractions. Both `original_idx` and entries of `removed`
 * are in the ORIGINAL 36-qubit numbering. */
static int
contracted_frame_index(int original_idx, const int *removed, int n_removed)
{
    int shift = 0;
    for (int k = 0; k < n_removed; ++k) {
        if (removed[k] < original_idx) ++shift;
    }
    return original_idx - shift;
}

irrep_status_t
irrep_happy_network_depth2_contracted(irrep_stabilizer_group_t *out,
                                      int *bulk_qubits_out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;

    int pairs[2 * IRREP_HAPPY_DEPTH2_N_CONTRACTIONS];
    irrep_stabilizer_group_t g_curr;
    irrep_status_t s = irrep_happy_network_depth2(&g_curr, pairs);
    if (s != IRREP_OK) return s;

    /* Track which original indices have been removed across iterations. */
    int removed[2 * IRREP_HAPPY_DEPTH2_N_CONTRACTIONS];
    int n_removed = 0;

    for (int k = 0; k < IRREP_HAPPY_DEPTH2_N_CONTRACTIONS; ++k) {
        int orig_a = pairs[2 * k + 0];
        int orig_b = pairs[2 * k + 1];
        int curr_a = contracted_frame_index(orig_a, removed, n_removed);
        int curr_b = contracted_frame_index(orig_b, removed, n_removed);

        irrep_stabilizer_group_t g_next;
        s = irrep_stabilizer_contract_bell(&g_curr, curr_a, curr_b, &g_next);
        irrep_stabilizer_group_free(&g_curr);
        if (s != IRREP_OK) return s;
        g_curr = g_next;

        removed[n_removed++] = orig_a;
        removed[n_removed++] = orig_b;
    }

    *out = g_curr;

    if (bulk_qubits_out != NULL) {
        for (int t = 0; t < IRREP_HAPPY_DEPTH2_N_TILES; ++t) {
            int orig = 6 * t;  /* original bulk qubit of tile t */
            bulk_qubits_out[t] = contracted_frame_index(orig, removed, n_removed);
        }
    }
    return IRREP_OK;
}
