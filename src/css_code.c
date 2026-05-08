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
