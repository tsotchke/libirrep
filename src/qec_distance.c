/* SPDX-License-Identifier: MIT */
/** @file qec_distance.c
 *  @brief Pauli multiplication + group membership + brute-force distance. */
#include <irrep/qec_distance.h>
#include <irrep/types.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

static inline int popcount_u64(uint64_t x)
{
    int c = 0;
    while (x) { x &= x - 1; ++c; }
    return c;
}

/* ====================================================================
 * Pauli multiplication
 * ==================================================================== */

irrep_status_t
irrep_pauli_multiply(irrep_pauli_t *out,
                     const irrep_pauli_t *p,
                     const irrep_pauli_t *q)
{
    if (out == NULL || p == NULL || q == NULL) return IRREP_ERR_INVALID_ARG;
    if (p->n != q->n) return IRREP_ERR_INVALID_ARG;
    if (out->n != p->n) return IRREP_ERR_INVALID_ARG;

    /* New support: out.x = p.x XOR q.x, out.z = p.z XOR q.z.
     * Read all source words into locals BEFORE writing to out, so out may
     * safely alias p or q without corrupting the carry computation. */
    int total_carry = 0;
    for (int w = 0; w < p->n_words; ++w) {
        uint64_t pxw = p->x[w], pzw = p->z[w];
        uint64_t qxw = q->x[w], qzw = q->z[w];
        out->x[w] = pxw ^ qxw;
        out->z[w] = pzw ^ qzw;
        total_carry += popcount_u64(pzw & qxw);
    }
    int sign_carry = (total_carry & 1) * 2;  /* (-1)^carry as i^(2·carry) */
    int new_sign = ((int)p->sign + (int)q->sign + sign_carry) & 3;
    out->sign = (irrep_pauli_sign_t)new_sign;
    return IRREP_OK;
}

/* ====================================================================
 * Stabilizer-group membership via F₂ Gaussian elimination
 * ==================================================================== */

static inline int bit_get(const uint64_t *row, int col)
{
    return (row[col / 64] >> (col % 64)) & 1;
}

static inline void bit_set(uint64_t *row, int col)
{
    row[col / 64] |= ((uint64_t)1 << (col % 64));
}

static inline void row_xor(uint64_t *dst, const uint64_t *src, int n_words)
{
    for (int w = 0; w < n_words; ++w) dst[w] ^= src[w];
}

int
irrep_stabilizer_group_contains(const irrep_stabilizer_group_t *g,
                                const irrep_pauli_t *p)
{
    if (g == NULL || p == NULL) return -1;
    if (g->n != p->n) return -1;
    int n = g->n;
    int n_cols = 2 * n;
    int n_words = (n_cols + 63) / 64;
    int m = g->n_generators;
    if (m == 0) {
        /* Identity is trivially in the empty group. */
        for (int w = 0; w < p->n_words; ++w) {
            if (p->x[w] != 0 || p->z[w] != 0) return 0;
        }
        return 1;
    }

    /* Allocate a working m × 2n F₂ matrix (rows = generators' symplectic
     * vectors). Layout: bits 0..n-1 = x-part, bits n..2n-1 = z-part. */
    uint64_t *rows_storage = (uint64_t *)calloc((size_t)m * n_words, sizeof(uint64_t));
    if (rows_storage == NULL) return -1;
    uint64_t **rows = (uint64_t **)malloc((size_t)m * sizeof(uint64_t *));
    if (rows == NULL) { free(rows_storage); return -1; }
    for (int i = 0; i < m; ++i) rows[i] = rows_storage + (size_t)i * n_words;

    for (int i = 0; i < m; ++i) {
        const irrep_pauli_t *gen = &g->gens[i];
        for (int q = 0; q < n; ++q) {
            if ((gen->x[q / 64] >> (q % 64)) & 1) bit_set(rows[i], q);
            if ((gen->z[q / 64] >> (q % 64)) & 1) bit_set(rows[i], n + q);
        }
    }

    /* Target vector for p. */
    uint64_t *target = (uint64_t *)calloc(n_words, sizeof(uint64_t));
    if (target == NULL) { free(rows); free(rows_storage); return -1; }
    for (int q = 0; q < n; ++q) {
        if ((p->x[q / 64] >> (q % 64)) & 1) bit_set(target, q);
        if ((p->z[q / 64] >> (q % 64)) & 1) bit_set(target, n + q);
    }

    /* Gaussian elimination: bring rows to row-echelon form. */
    int rank = 0;
    for (int col = 0; col < n_cols && rank < m; ++col) {
        int pivot = -1;
        for (int r = rank; r < m; ++r) {
            if (bit_get(rows[r], col)) { pivot = r; break; }
        }
        if (pivot < 0) continue;
        if (pivot != rank) {
            uint64_t *tmp = rows[pivot];
            rows[pivot] = rows[rank];
            rows[rank] = tmp;
        }
        for (int r = 0; r < m; ++r) {
            if (r == rank) continue;
            if (bit_get(rows[r], col)) row_xor(rows[r], rows[rank], n_words);
        }
        ++rank;
    }

    /* Reduce target by the row-reduced rows. */
    for (int r = 0; r < rank; ++r) {
        int pivot_col = -1;
        for (int c = 0; c < n_cols; ++c) {
            if (bit_get(rows[r], c)) { pivot_col = c; break; }
        }
        if (pivot_col < 0) continue;
        if (bit_get(target, pivot_col)) row_xor(target, rows[r], n_words);
    }

    int contained = 1;
    for (int w = 0; w < n_words; ++w) {
        if (target[w] != 0) { contained = 0; break; }
    }
    free(target);
    free(rows);
    free(rows_storage);
    return contained;
}

/* ====================================================================
 * Brute-force distance
 * ==================================================================== */

/* Recurse over weight-w Paulis: choose w qubit positions, then enumerate
 * 3^w letter combinations on those positions. For each, check if it's
 * a logical operator (commutes with all stabilizers, NOT in the group).
 *
 * Returns 1 on first hit (caller stops), 0 if no hit found at this weight. */
static int
search_weight_w(const irrep_stabilizer_group_t *g, int weight,
                int start_qubit, int placed,
                int *positions, irrep_pauli_t *workspace)
{
    int n = g->n;
    if (placed == weight) {
        /* Iterate over 3^weight letter assignments. Use uint64_t to avoid
         * signed overflow at weight ≥ 20 (3^20 > INT32_MAX). */
        uint64_t total = 1;
        for (int i = 0; i < weight; ++i) total *= 3;
        for (uint64_t code = 0; code < total; ++code) {
            /* Reset workspace. */
            for (int w = 0; w < workspace->n_words; ++w) {
                workspace->x[w] = 0;
                workspace->z[w] = 0;
            }
            workspace->sign = IRREP_PAULI_SIGN_POS;
            uint64_t c = code;
            for (int i = 0; i < weight; ++i) {
                int letter = (int)(c % 3) + 1; /* 1=X, 2=Y, 3=Z */
                c /= 3;
                irrep_pauli_set(workspace, positions[i],
                                (irrep_pauli_letter_t)letter);
            }
            /* Check (a) commutes with every stabilizer, (b) NOT in group. */
            bool commutes = true;
            for (int j = 0; j < g->n_generators; ++j) {
                if (!irrep_pauli_commute(&g->gens[j], workspace)) {
                    commutes = false;
                    break;
                }
            }
            if (!commutes) continue;
            int in_g = irrep_stabilizer_group_contains(g, workspace);
            if (in_g == 0) return 1;  /* found a non-trivial logical */
        }
        return 0;
    }
    /* Recurse: pick next qubit position. */
    int qubits_left = weight - placed;
    for (int q = start_qubit; q <= n - qubits_left; ++q) {
        positions[placed] = q;
        if (search_weight_w(g, weight, q + 1, placed + 1, positions, workspace))
            return 1;
    }
    return 0;
}

int
irrep_qec_distance_brute(const irrep_stabilizer_group_t *g, int max_weight)
{
    if (g == NULL || max_weight <= 0) return -1;
    if (max_weight > g->n) max_weight = g->n;

    irrep_pauli_t workspace;
    if (irrep_pauli_new(&workspace, g->n) != IRREP_OK) return -1;
    int *positions = (int *)malloc(max_weight * sizeof(int));
    if (positions == NULL) { irrep_pauli_free(&workspace); return -1; }

    int distance = max_weight + 1;
    for (int w = 1; w <= max_weight; ++w) {
        if (search_weight_w(g, w, 0, 0, positions, &workspace)) {
            distance = w;
            break;
        }
    }
    free(positions);
    irrep_pauli_free(&workspace);
    return distance;
}
