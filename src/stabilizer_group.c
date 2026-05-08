/* SPDX-License-Identifier: MIT */
/** @file stabilizer_group.c
 *  @brief Implementation of the abstract stabilizer-group machinery.
 */
#include <irrep/stabilizer_group.h>
#include <irrep/types.h>

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* ====================================================================
 * Pauli operators
 * ==================================================================== */

irrep_status_t
irrep_pauli_new(irrep_pauli_t *out, int n)
{
    if (out == NULL || n <= 0) return IRREP_ERR_INVALID_ARG;
    int n_words = (n + 63) / 64;
    out->n = n;
    out->n_words = n_words;
    out->x = (uint64_t *)calloc(n_words, sizeof(uint64_t));
    out->z = (uint64_t *)calloc(n_words, sizeof(uint64_t));
    if (out->x == NULL || out->z == NULL) {
        free(out->x);
        free(out->z);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    out->sign = IRREP_PAULI_SIGN_POS;
    return IRREP_OK;
}

void
irrep_pauli_free(irrep_pauli_t *p)
{
    if (p == NULL) return;
    free(p->x); p->x = NULL;
    free(p->z); p->z = NULL;
    p->n = 0;
    p->n_words = 0;
    p->sign = IRREP_PAULI_SIGN_POS;
}

static inline int qubit_word(int qubit) { return qubit / 64; }
static inline int qubit_bit(int qubit) { return qubit % 64; }

irrep_status_t
irrep_pauli_set(irrep_pauli_t *p, int qubit, irrep_pauli_letter_t letter)
{
    if (p == NULL || qubit < 0 || qubit >= p->n) return IRREP_ERR_INVALID_ARG;
    int w = qubit_word(qubit);
    uint64_t mask = (uint64_t)1 << qubit_bit(qubit);
    uint64_t inv_mask = ~mask;
    p->x[w] &= inv_mask;
    p->z[w] &= inv_mask;
    switch (letter) {
        case IRREP_PAULI_LETTER_I: break;
        case IRREP_PAULI_LETTER_X: p->x[w] |= mask; break;
        case IRREP_PAULI_LETTER_Y: p->x[w] |= mask; p->z[w] |= mask; break;
        case IRREP_PAULI_LETTER_Z: p->z[w] |= mask; break;
        default: return IRREP_ERR_INVALID_ARG;
    }
    return IRREP_OK;
}

irrep_pauli_letter_t
irrep_pauli_get(const irrep_pauli_t *p, int qubit)
{
    if (p == NULL || qubit < 0 || qubit >= p->n) return IRREP_PAULI_LETTER_I;
    int w = qubit_word(qubit);
    uint64_t mask = (uint64_t)1 << qubit_bit(qubit);
    bool has_x = (p->x[w] & mask) != 0;
    bool has_z = (p->z[w] & mask) != 0;
    if (has_x && has_z) return IRREP_PAULI_LETTER_Y;
    if (has_x) return IRREP_PAULI_LETTER_X;
    if (has_z) return IRREP_PAULI_LETTER_Z;
    return IRREP_PAULI_LETTER_I;
}

/* Population count of a single uint64_t. */
static inline int popcount_u64(uint64_t x)
{
    int c = 0;
    while (x) { x &= x - 1; ++c; }
    return c;
}

int
irrep_pauli_symp_inner(const irrep_pauli_t *p, const irrep_pauli_t *q)
{
    if (p == NULL || q == NULL) return -1;
    if (p->n != q->n) return -1;
    /* ⟨P, Q⟩_symp = sum_i (x_P[i] z_Q[i] + z_P[i] x_Q[i]) (mod 2)
     * Implement word-wise: count parity of (x_P AND z_Q) XOR (z_P AND x_Q). */
    uint64_t parity = 0;
    for (int w = 0; w < p->n_words; ++w) {
        parity ^= p->x[w] & q->z[w];
        parity ^= p->z[w] & q->x[w];
    }
    return popcount_u64(parity) & 1;
}

bool
irrep_pauli_commute(const irrep_pauli_t *p, const irrep_pauli_t *q)
{
    int s = irrep_pauli_symp_inner(p, q);
    return s == 0;
}

int
irrep_pauli_weight(const irrep_pauli_t *p)
{
    if (p == NULL) return -1;
    int w_total = 0;
    for (int w = 0; w < p->n_words; ++w) {
        /* A qubit has non-trivial action iff x[i] | z[i]. */
        w_total += popcount_u64(p->x[w] | p->z[w]);
    }
    return w_total;
}

/* ====================================================================
 * Stabilizer group
 * ==================================================================== */

irrep_status_t
irrep_stabilizer_group_new(irrep_stabilizer_group_t *out,
                           int n, int n_generators)
{
    if (out == NULL || n <= 0 || n_generators < 0) return IRREP_ERR_INVALID_ARG;
    out->n = n;
    out->n_generators = n_generators;
    if (n_generators == 0) {
        out->gens = NULL;
        return IRREP_OK;
    }
    out->gens = (irrep_pauli_t *)calloc(n_generators, sizeof(irrep_pauli_t));
    if (out->gens == NULL) return IRREP_ERR_OUT_OF_MEMORY;
    for (int i = 0; i < n_generators; ++i) {
        irrep_status_t s = irrep_pauli_new(&out->gens[i], n);
        if (s != IRREP_OK) {
            for (int j = 0; j < i; ++j) irrep_pauli_free(&out->gens[j]);
            free(out->gens);
            out->gens = NULL;
            return s;
        }
    }
    return IRREP_OK;
}

void
irrep_stabilizer_group_free(irrep_stabilizer_group_t *g)
{
    if (g == NULL) return;
    if (g->gens) {
        for (int i = 0; i < g->n_generators; ++i) {
            irrep_pauli_free(&g->gens[i]);
        }
        free(g->gens);
        g->gens = NULL;
    }
    g->n = 0;
    g->n_generators = 0;
}

irrep_status_t
irrep_stabilizer_group_check_commutativity(const irrep_stabilizer_group_t *g)
{
    if (g == NULL) return IRREP_ERR_INVALID_ARG;
    for (int i = 0; i < g->n_generators; ++i) {
        for (int j = i + 1; j < g->n_generators; ++j) {
            if (!irrep_pauli_commute(&g->gens[i], &g->gens[j])) {
                return IRREP_ERR_PRECONDITION;
            }
        }
    }
    return IRREP_OK;
}

irrep_status_t
irrep_stabilizer_syndrome(const irrep_stabilizer_group_t *g,
                          const irrep_pauli_t *error,
                          uint64_t *syndrome)
{
    if (g == NULL || error == NULL || syndrome == NULL) {
        return IRREP_ERR_INVALID_ARG;
    }
    if (error->n != g->n) return IRREP_ERR_INVALID_ARG;
    int n_syndrome_words = (g->n_generators + 63) / 64;
    memset(syndrome, 0, n_syndrome_words * sizeof(uint64_t));
    for (int i = 0; i < g->n_generators; ++i) {
        int s = irrep_pauli_symp_inner(&g->gens[i], error);
        if (s == 1) {
            int w = i / 64;
            int b = i % 64;
            syndrome[w] |= ((uint64_t)1 << b);
        }
    }
    return IRREP_OK;
}
