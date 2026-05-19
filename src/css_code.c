/* SPDX-License-Identifier: MIT */
/** @file css_code.c
 *  @brief Implementation of CSS-code primitives. */
#include <irrep/css_code.h>
#include <irrep/types.h>

#include <stdlib.h>
#include <string.h>

/* ====================================================================
 * Bit-packed parity-check matrix
 * ==================================================================== */

irrep_status_t
irrep_parity_matrix_new(irrep_parity_matrix_t *out, int n_rows, int n_cols)
{
    if (out == NULL || n_rows <= 0 || n_cols <= 0) return IRREP_ERR_INVALID_ARG;
    out->n_rows = n_rows;
    out->n_cols = n_cols;
    out->n_words = (n_cols + 63) / 64;
    size_t total = (size_t)n_rows * (size_t)out->n_words;
    out->data = (uint64_t *)calloc(total, sizeof(uint64_t));
    if (out->data == NULL) return IRREP_ERR_OUT_OF_MEMORY;
    return IRREP_OK;
}

void
irrep_parity_matrix_free(irrep_parity_matrix_t *m)
{
    if (m == NULL) return;
    free(m->data); m->data = NULL;
    m->n_rows = 0;
    m->n_cols = 0;
    m->n_words = 0;
}

static inline uint64_t *
row_ptr(irrep_parity_matrix_t *m, int row)
{
    return m->data + (size_t)row * (size_t)m->n_words;
}

static inline const uint64_t *
row_ptr_const(const irrep_parity_matrix_t *m, int row)
{
    return m->data + (size_t)row * (size_t)m->n_words;
}

irrep_status_t
irrep_parity_matrix_set(irrep_parity_matrix_t *m, int row, int col)
{
    if (m == NULL) return IRREP_ERR_INVALID_ARG;
    if (row < 0 || row >= m->n_rows) return IRREP_ERR_INVALID_ARG;
    if (col < 0 || col >= m->n_cols) return IRREP_ERR_INVALID_ARG;
    int w = col / 64;
    int b = col % 64;
    row_ptr(m, row)[w] |= ((uint64_t)1 << b);
    return IRREP_OK;
}

int
irrep_parity_matrix_get(const irrep_parity_matrix_t *m, int row, int col)
{
    if (m == NULL) return -1;
    if (row < 0 || row >= m->n_rows) return -1;
    if (col < 0 || col >= m->n_cols) return -1;
    int w = col / 64;
    int b = col % 64;
    return (row_ptr_const(m, row)[w] >> b) & 1;
}

static inline int popcount_u64(uint64_t x)
{
    int c = 0;
    while (x) { x &= x - 1; ++c; }
    return c;
}

int
irrep_parity_matrix_row_inner(const irrep_parity_matrix_t *a, int row_a,
                              const irrep_parity_matrix_t *b, int row_b)
{
    if (a == NULL || b == NULL) return -1;
    if (row_a < 0 || row_a >= a->n_rows) return -1;
    if (row_b < 0 || row_b >= b->n_rows) return -1;
    if (a->n_cols != b->n_cols) return -1;
    const uint64_t *ra = row_ptr_const(a, row_a);
    const uint64_t *rb = row_ptr_const(b, row_b);
    uint64_t parity = 0;
    for (int w = 0; w < a->n_words; ++w) {
        parity ^= ra[w] & rb[w];
    }
    return popcount_u64(parity) & 1;
}

/* ====================================================================
 * CSS code
 * ==================================================================== */

irrep_status_t
irrep_css_code_new(irrep_css_code_t *out, int n, int m_X, int m_Z)
{
    if (out == NULL || n <= 0 || m_X < 0 || m_Z < 0) return IRREP_ERR_INVALID_ARG;
    out->n = n;
    /* Zero-row sides are legitimate for degenerate CSS codes (e.g. a
     * classical-only X-code embedded as a CSS instance). Skip the
     * allocation but keep the struct valid. */
    if (m_X > 0) {
        irrep_status_t s = irrep_parity_matrix_new(&out->H_X, m_X, n);
        if (s != IRREP_OK) return s;
    } else {
        out->H_X.n_rows = 0;
        out->H_X.n_cols = n;
        out->H_X.n_words = (n + 63) / 64;
        out->H_X.data = NULL;
    }
    if (m_Z > 0) {
        irrep_status_t s = irrep_parity_matrix_new(&out->H_Z, m_Z, n);
        if (s != IRREP_OK) {
            if (m_X > 0) irrep_parity_matrix_free(&out->H_X);
            return s;
        }
    } else {
        out->H_Z.n_rows = 0;
        out->H_Z.n_cols = n;
        out->H_Z.n_words = (n + 63) / 64;
        out->H_Z.data = NULL;
    }
    return IRREP_OK;
}

void
irrep_css_code_free(irrep_css_code_t *c)
{
    if (c == NULL) return;
    irrep_parity_matrix_free(&c->H_X);
    irrep_parity_matrix_free(&c->H_Z);
    c->n = 0;
}

irrep_status_t
irrep_css_code_verify(const irrep_css_code_t *c)
{
    if (c == NULL) return IRREP_ERR_INVALID_ARG;
    /* CSS orthogonality: every row of H_X is F₂-orthogonal to every row of H_Z. */
    for (int i = 0; i < c->H_X.n_rows; ++i) {
        for (int j = 0; j < c->H_Z.n_rows; ++j) {
            if (irrep_parity_matrix_row_inner(&c->H_X, i, &c->H_Z, j) != 0) {
                return IRREP_ERR_PRECONDITION;
            }
        }
    }
    return IRREP_OK;
}

irrep_status_t
irrep_css_code_to_stabilizer_group(const irrep_css_code_t *c,
                                   irrep_stabilizer_group_t *out)
{
    if (c == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    int m_X = c->H_X.n_rows;
    int m_Z = c->H_Z.n_rows;
    int total = m_X + m_Z;
    irrep_status_t s = irrep_stabilizer_group_new(out, c->n, total);
    if (s != IRREP_OK) return s;
    /* X-stabilizers: x = row of H_X, z = 0. */
    for (int i = 0; i < m_X; ++i) {
        const uint64_t *src = c->H_X.data + (size_t)i * (size_t)c->H_X.n_words;
        memcpy(out->gens[i].x, src,
               (size_t)c->H_X.n_words * sizeof(uint64_t));
        /* z is already zero from calloc. */
    }
    /* Z-stabilizers: x = 0, z = row of H_Z. */
    for (int i = 0; i < m_Z; ++i) {
        const uint64_t *src = c->H_Z.data + (size_t)i * (size_t)c->H_Z.n_words;
        memcpy(out->gens[m_X + i].z, src,
               (size_t)c->H_Z.n_words * sizeof(uint64_t));
    }
    return IRREP_OK;
}


/* ====================================================================
 * F₂-rank of a parity matrix and logical-qubit count of a CSS code.
 * ==================================================================== */

int
irrep_parity_matrix_rank(const irrep_parity_matrix_t *H)
{
    if (H == NULL) return -1;
    int m = H->n_rows;
    int n = H->n_cols;
    if (m == 0 || n == 0) return 0;

    /* Work-copy of H. */
    irrep_parity_matrix_t W;
    if (irrep_parity_matrix_new(&W, m, n) != IRREP_OK) return -1;
    for (int r = 0; r < m; ++r) {
        for (int w = 0; w < H->n_words; ++w) {
            W.data[(size_t)r * (size_t)W.n_words + w] =
                H->data[(size_t)r * (size_t)H->n_words + w];
        }
    }

    int pivot_row = 0;
    for (int c = 0; c < n && pivot_row < m; ++c) {
        int r = -1;
        for (int i = pivot_row; i < m; ++i) {
            if (irrep_parity_matrix_get(&W, i, c) == 1) { r = i; break; }
        }
        if (r == -1) continue;
        if (r != pivot_row) {
            uint64_t *a = &W.data[(size_t)pivot_row * (size_t)W.n_words];
            uint64_t *b = &W.data[(size_t)r          * (size_t)W.n_words];
            for (int w = 0; w < W.n_words; ++w) {
                uint64_t t = a[w]; a[w] = b[w]; b[w] = t;
            }
        }
        for (int i = 0; i < m; ++i) {
            if (i == pivot_row) continue;
            if (irrep_parity_matrix_get(&W, i, c) == 1) {
                uint64_t *dst = &W.data[(size_t)i         * (size_t)W.n_words];
                uint64_t *src = &W.data[(size_t)pivot_row * (size_t)W.n_words];
                for (int w = 0; w < W.n_words; ++w) dst[w] ^= src[w];
            }
        }
        ++pivot_row;
    }
    irrep_parity_matrix_free(&W);
    return pivot_row;
}

int
irrep_css_code_logical_qubits(const irrep_css_code_t *c)
{
    if (c == NULL) return -1;
    int rX = irrep_parity_matrix_rank(&c->H_X);
    int rZ = irrep_parity_matrix_rank(&c->H_Z);
    if (rX < 0 || rZ < 0) return -1;
    int k = c->n - rX - rZ;
    return k < 0 ? -1 : k;
}

#include <irrep/qec_distance.h>
#include <irrep/stabilizer_group.h>

int
irrep_css_code_distance(const irrep_css_code_t *c, int max_weight)
{
    if (c == NULL || max_weight <= 0) return -1;
    irrep_stabilizer_group_t g;
    irrep_status_t s = irrep_css_code_to_stabilizer_group(c, &g);
    if (s != IRREP_OK) return -1;
    int d = irrep_qec_distance_brute(&g, max_weight);
    irrep_stabilizer_group_free(&g);
    return d;
}

/* ====================================================================
 * Automatic minimum-weight logical X̄ discovery via brute search.
 *
 * Algorithm:
 *   For each w = 1..max_weight:
 *     For each C(n, w) support choice:
 *       Form the X^v Pauli on those qubits.
 *       Check: commutes with every Z-stab (v · z = 0 for every z row).
 *       Check: not in rowspace(H_X).
 *       Return the first hit.
 *   Return max_weight + 1 if no logical found in range.
 *
 * Helper: pauli_in_X_rowspace tests if v ∈ rowspan(H_X) via the same
 * F₂ Gaussian elimination as irrep_pauli_in_stabilizer_span.
 * ==================================================================== */

/* Test whether the F_2 vector v of length n is in the row-space of H. */
static int
parity_vec_in_rowspace(const irrep_parity_matrix_t *H, const int *v_support, int sz)
{
    int m = H->n_rows;
    int n = H->n_cols;
    if (m == 0) return 0; /* Empty rowspace contains only zero. */

    /* Build augmented matrix: [H | v]^T then row-reduce. We just need
     * to check if v is a linear combination of rows of H. */
    irrep_parity_matrix_t W;
    if (irrep_parity_matrix_new(&W, m + 1, n) != IRREP_OK) return -1;
    /* Copy H rows. */
    for (int r = 0; r < m; ++r) {
        for (int w = 0; w < H->n_words; ++w) {
            W.data[(size_t)r * (size_t)W.n_words + w] =
                H->data[(size_t)r * (size_t)H->n_words + w];
        }
    }
    /* Append v as the last row. */
    for (int i = 0; i < sz; ++i) {
        irrep_parity_matrix_set(&W, m, v_support[i]);
    }
    /* Row-reduce. v is in rowspan(H) iff rank(W) == rank(H), i.e., the
     * augmented matrix has the same rank as H alone. */
    int rank_H = irrep_parity_matrix_rank(H);
    int rank_W = irrep_parity_matrix_rank(&W);
    irrep_parity_matrix_free(&W);
    return (rank_W == rank_H) ? 1 : 0;
}

/* Test whether x-support v commutes with all rows of H_Z (i.e., the
 * F_2 inner product v · z_row = 0 for every Z-stabilizer row). */
static int
support_commutes_all(const irrep_parity_matrix_t *H_Z, const int *v_support, int sz)
{
    for (int r = 0; r < H_Z->n_rows; ++r) {
        int parity = 0;
        for (int i = 0; i < sz; ++i) {
            parity ^= irrep_parity_matrix_get(H_Z, r, v_support[i]);
        }
        if (parity != 0) return 0;
    }
    return 1;
}

/* Iterate over weight-w supports of {0..n-1}. Lex order via a
 * combination-counter. Returns 0 if no next combination, 1 otherwise. */
static int
combo_next(int *combo, int w, int n)
{
    /* Find rightmost element that can be incremented. */
    int i = w - 1;
    while (i >= 0 && combo[i] == n - w + i) --i;
    if (i < 0) return 0;
    ++combo[i];
    for (int j = i + 1; j < w; ++j) combo[j] = combo[j - 1] + 1;
    return 1;
}

int
irrep_css_code_compute_logical_X(const irrep_css_code_t *c, int max_weight,
                                  irrep_pauli_t *out_pauli)
{
    if (c == NULL || out_pauli == NULL || max_weight <= 0) return -1;
    int n = c->n;
    if (n <= 0) return -1;

    int *combo = (int *)malloc((size_t)max_weight * sizeof(int));
    if (combo == NULL) return -1;

    for (int w = 1; w <= max_weight; ++w) {
        /* Initialise to {0, 1, ..., w-1}. */
        for (int i = 0; i < w; ++i) combo[i] = i;
        if (w > n) break;
        do {
            /* Test this support: commutes with all Z-stabs AND not in
             * the X-stab rowspace. */
            if (support_commutes_all(&c->H_Z, combo, w) &&
                parity_vec_in_rowspace(&c->H_X, combo, w) == 0) {
                /* Found minimum-weight logical X. */
                irrep_status_t s = irrep_pauli_new(out_pauli, n);
                if (s != IRREP_OK) { free(combo); return -1; }
                for (int i = 0; i < w; ++i) {
                    irrep_pauli_set(out_pauli, combo[i], IRREP_PAULI_LETTER_X);
                }
                free(combo);
                return w;
            }
        } while (combo_next(combo, w, n));
    }
    free(combo);
    return max_weight + 1; /* No logical found in range. */
}

int
irrep_css_code_compute_logical_Z(const irrep_css_code_t *c, int max_weight,
                                  irrep_pauli_t *out_pauli)
{
    if (c == NULL || out_pauli == NULL || max_weight <= 0) return -1;
    int n = c->n;
    if (n <= 0) return -1;

    int *combo = (int *)malloc((size_t)max_weight * sizeof(int));
    if (combo == NULL) return -1;

    for (int w = 1; w <= max_weight; ++w) {
        for (int i = 0; i < w; ++i) combo[i] = i;
        if (w > n) break;
        do {
            /* Z-type: commutes with H_X (v · x = 0 for every x row of H_X)
             * AND not in rowspace(H_Z). */
            if (support_commutes_all(&c->H_X, combo, w) &&
                parity_vec_in_rowspace(&c->H_Z, combo, w) == 0) {
                irrep_status_t s = irrep_pauli_new(out_pauli, n);
                if (s != IRREP_OK) { free(combo); return -1; }
                for (int i = 0; i < w; ++i) {
                    irrep_pauli_set(out_pauli, combo[i], IRREP_PAULI_LETTER_Z);
                }
                free(combo);
                return w;
            }
        } while (combo_next(combo, w, n));
    }
    free(combo);
    return max_weight + 1;
}

irrep_status_t
irrep_css_code_audit(const irrep_css_code_t *c, int max_distance_check,
                     irrep_css_code_params_t *out)
{
    if (c == NULL || out == NULL || max_distance_check <= 0)
        return IRREP_ERR_INVALID_ARG;

    out->n      = c->n;
    out->m_X    = c->H_X.n_rows;
    out->m_Z    = c->H_Z.n_rows;
    out->rank_X = irrep_parity_matrix_rank(&c->H_X);
    out->rank_Z = irrep_parity_matrix_rank(&c->H_Z);
    if (out->rank_X < 0 || out->rank_Z < 0) return IRREP_ERR_OUT_OF_MEMORY;
    out->k = out->n - out->rank_X - out->rank_Z;
    if (out->k < 0) out->k = 0;

    /* Compute d_X and d_Z without keeping the Pauli supports. */
    irrep_pauli_t tmp;
    out->d_X = irrep_css_code_compute_logical_X(c, max_distance_check, &tmp);
    if (out->d_X >= 0 && out->d_X <= max_distance_check) {
        out->d_X_bounded = 1;
        irrep_pauli_free(&tmp);
    } else {
        out->d_X_bounded = 0;
    }
    out->d_Z = irrep_css_code_compute_logical_Z(c, max_distance_check, &tmp);
    if (out->d_Z >= 0 && out->d_Z <= max_distance_check) {
        out->d_Z_bounded = 1;
        irrep_pauli_free(&tmp);
    } else {
        out->d_Z_bounded = 0;
    }
    out->d = (out->d_X < out->d_Z) ? out->d_X : out->d_Z;

    return IRREP_OK;
}
