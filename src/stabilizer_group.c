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

/* ====================================================================
 * F₂ row-span test for stabilizer membership.
 *
 * Concatenates each generator's symplectic vector (x | z) into a row
 * of an `m × 2n` matrix M over F₂, augments with `p`'s symplectic
 * vector as the bottom row, and runs Gaussian elimination on the
 * (m+1) × 2n system. `p` is in the row-span of M iff elimination
 * reduces the augmented row to zero.
 *
 * Symplectic columns are packed as `n_words = ceil(2n / 64)` uint64_t
 * — but we re-use the per-Pauli `n_words` (= ceil(n/64)) twice, once
 * for x and once for z, so each row's pivot search runs over
 * 2·n_words words.
 * ==================================================================== */

/* Find the lowest set bit in two concatenated bit-arrays {x[], z[]},
 * each of `n_words` uint64_t (representing 2n bits in symplectic
 * convention, with x ordered before z). Returns the global bit index
 * in [0, 2n), or -1 if both arrays are zero. */
static int
symp_lsb(const uint64_t *x, const uint64_t *z, int n_words, int n_bits)
{
    for (int w = 0; w < n_words; ++w) {
        if (x[w] != 0) {
            uint64_t v = x[w];
            int b = 0;
            while ((v & 1) == 0) { v >>= 1; ++b; }
            int idx = 64 * w + b;
            return idx < n_bits ? idx : -1;
        }
    }
    /* z half, starting at column n_bits. */
    for (int w = 0; w < n_words; ++w) {
        if (z[w] != 0) {
            uint64_t v = z[w];
            int b = 0;
            while ((v & 1) == 0) { v >>= 1; ++b; }
            int idx = 64 * w + b;
            if (idx < n_bits) return n_bits + idx;
        }
    }
    return -1;
}

/* XOR row src into row dst (both have x and z components). */
static void
symp_row_xor(uint64_t *dst_x, uint64_t *dst_z,
             const uint64_t *src_x, const uint64_t *src_z,
             int n_words)
{
    for (int w = 0; w < n_words; ++w) {
        dst_x[w] ^= src_x[w];
        dst_z[w] ^= src_z[w];
    }
}

/* Test whether the bit at symplectic column `col` (0..2n-1) is set
 * in row (x[], z[]). Columns [0, n) live in x; [n, 2n) in z (offset
 * by n). */
static int
symp_bit(const uint64_t *x, const uint64_t *z, int n_bits, int col)
{
    const uint64_t *arr = (col < n_bits) ? x : z;
    int local = (col < n_bits) ? col : (col - n_bits);
    int w = local / 64, b = local % 64;
    return (int)((arr[w] >> b) & 1ULL);
}

int
irrep_pauli_in_stabilizer_span(const irrep_stabilizer_group_t *g,
                               const irrep_pauli_t *p)
{
    if (g == NULL || p == NULL) return -1;
    if (g->n != p->n) return -1;
    if (g->n_generators == 0) {
        /* Span is the identity-only group; p ∈ span iff p == I. */
        for (int w = 0; w < p->n_words; ++w) {
            if (p->x[w] != 0 || p->z[w] != 0) return 0;
        }
        return 1;
    }

    const int n = g->n;
    const int n_words = p->n_words;

    /* Copy the m generators + p into a working matrix of (m+1) rows,
     * each carrying its own (x, z) word arrays. Row m is the target p. */
    const int n_rows = g->n_generators + 1;
    uint64_t **X = (uint64_t **)calloc((size_t)n_rows, sizeof(uint64_t *));
    uint64_t **Z = (uint64_t **)calloc((size_t)n_rows, sizeof(uint64_t *));
    if (!X || !Z) { free(X); free(Z); return -1; }
    for (int r = 0; r < n_rows; ++r) {
        X[r] = (uint64_t *)calloc((size_t)n_words, sizeof(uint64_t));
        Z[r] = (uint64_t *)calloc((size_t)n_words, sizeof(uint64_t));
        if (!X[r] || !Z[r]) {
            for (int rr = 0; rr <= r; ++rr) { free(X[rr]); free(Z[rr]); }
            free(X); free(Z); return -1;
        }
        if (r < g->n_generators) {
            memcpy(X[r], g->gens[r].x, (size_t)n_words * sizeof(uint64_t));
            memcpy(Z[r], g->gens[r].z, (size_t)n_words * sizeof(uint64_t));
        } else {
            memcpy(X[r], p->x, (size_t)n_words * sizeof(uint64_t));
            memcpy(Z[r], p->z, (size_t)n_words * sizeof(uint64_t));
        }
    }

    /* Gaussian elimination on the first m rows. For each pivot column
     * (lowest-bit-first), find a row with that bit set in {pivot_row..m-1},
     * swap to position pivot_row, then XOR it into all other rows that
     * have that bit set (including row m = p). */
    int pivot_row = 0;
    for (int r = 0; r < g->n_generators && pivot_row < g->n_generators; ++r) {
        /* Find pivot column of row pivot_row (lowest set bit). */
        int piv_col = symp_lsb(X[pivot_row], Z[pivot_row], n_words, n);
        if (piv_col < 0) {
            /* This row is zero — skip; continue with next unique row. */
            /* Swap it to the back to keep pivot_row contiguous. */
            for (int r2 = pivot_row + 1; r2 < g->n_generators; ++r2) {
                if (symp_lsb(X[r2], Z[r2], n_words, n) >= 0) {
                    uint64_t *tx = X[pivot_row], *tz = Z[pivot_row];
                    X[pivot_row] = X[r2]; Z[pivot_row] = Z[r2];
                    X[r2] = tx; Z[r2] = tz;
                    break;
                }
            }
            piv_col = symp_lsb(X[pivot_row], Z[pivot_row], n_words, n);
            if (piv_col < 0) break;  /* all remaining rows zero */
        }
        /* Eliminate piv_col from every other row (including the target). */
        for (int r2 = 0; r2 < n_rows; ++r2) {
            if (r2 == pivot_row) continue;
            if (symp_bit(X[r2], Z[r2], n, piv_col)) {
                symp_row_xor(X[r2], Z[r2], X[pivot_row], Z[pivot_row], n_words);
            }
        }
        ++pivot_row;
        (void)r;  /* loop variable unused; pivot_row drives progress */
    }

    /* p is in the span iff its symplectic vector reduces to zero. */
    int residual = 0;
    for (int w = 0; w < n_words; ++w) {
        if (X[g->n_generators][w] != 0 || Z[g->n_generators][w] != 0) {
            residual = 1;
            break;
        }
    }

    for (int r = 0; r < n_rows; ++r) { free(X[r]); free(Z[r]); }
    free(X); free(Z);
    return residual == 0 ? 1 : 0;
}
