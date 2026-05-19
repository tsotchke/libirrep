/* SPDX-License-Identifier: MIT */
/** @file single_shot.c
 *  @brief Single-shot QEC framework: meta-checks on stabilizer syndromes. */
#include <irrep/single_shot.h>
#include <irrep/types.h>

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

irrep_status_t
irrep_single_shot_code_new(irrep_single_shot_code_t *out,
                           irrep_css_code_t *css,
                           int m_X_meta, int m_Z_meta)
{
    if (out == NULL || css == NULL) return IRREP_ERR_INVALID_ARG;
    if (m_X_meta <= 0 || m_Z_meta <= 0) return IRREP_ERR_INVALID_ARG;

    /* Take ownership of css by shallow copy, then fully zero the source so
     * the caller's free is a no-op AND any subsequent stray accesses
     * (e.g. parity_matrix_set, parity_matrix_get) hit guarded paths
     * rather than dereferencing NULL with non-zero n_rows. */
    out->css = *css;
    css->n = 0;
    css->H_X.data = NULL;
    css->H_X.n_rows = 0;
    css->H_X.n_cols = 0;
    css->H_X.n_words = 0;
    css->H_Z.data = NULL;
    css->H_Z.n_rows = 0;
    css->H_Z.n_cols = 0;
    css->H_Z.n_words = 0;

    irrep_status_t s = irrep_parity_matrix_new(&out->M_X, m_X_meta, out->css.H_X.n_rows);
    if (s != IRREP_OK) {
        irrep_css_code_free(&out->css);
        return s;
    }
    s = irrep_parity_matrix_new(&out->M_Z, m_Z_meta, out->css.H_Z.n_rows);
    if (s != IRREP_OK) {
        irrep_parity_matrix_free(&out->M_X);
        irrep_css_code_free(&out->css);
        return s;
    }
    return IRREP_OK;
}

void
irrep_single_shot_code_free(irrep_single_shot_code_t *c)
{
    if (c == NULL) return;
    irrep_parity_matrix_free(&c->M_X);
    irrep_parity_matrix_free(&c->M_Z);
    irrep_css_code_free(&c->css);
}

/* M · H = 0 in F₂: every (M-row, H-col) pair has even inner product. */
static int
matrix_product_zero(const irrep_parity_matrix_t *M,
                    const irrep_parity_matrix_t *H)
{
    /* M is m_meta × m_check; H is m_check × n_qubits. */
    if (M->n_cols != H->n_rows) return 0;
    for (int i = 0; i < M->n_rows; ++i) {
        for (int j = 0; j < H->n_cols; ++j) {
            int parity = 0;
            for (int k = 0; k < M->n_cols; ++k) {
                parity ^= (irrep_parity_matrix_get(M, i, k) &
                           irrep_parity_matrix_get(H, k, j));
            }
            if (parity != 0) return 0;
        }
    }
    return 1;
}

irrep_status_t
irrep_single_shot_verify_meta(const irrep_single_shot_code_t *c)
{
    if (c == NULL) return IRREP_ERR_INVALID_ARG;
    if (!matrix_product_zero(&c->M_X, &c->css.H_X)) return IRREP_ERR_PRECONDITION;
    if (!matrix_product_zero(&c->M_Z, &c->css.H_Z)) return IRREP_ERR_PRECONDITION;
    return IRREP_OK;
}

/* Apply meta matrix M to syndrome bit-vector s: out_i = (M · s)_i mod 2. */
static irrep_status_t
apply_meta(const irrep_parity_matrix_t *M,
           const uint64_t *syndrome, uint64_t *out)
{
    int n_out_words = (M->n_rows + 63) / 64;
    memset(out, 0, n_out_words * sizeof(uint64_t));
    for (int i = 0; i < M->n_rows; ++i) {
        int parity = 0;
        for (int j = 0; j < M->n_cols; ++j) {
            int s_bit = (syndrome[j / 64] >> (j % 64)) & 1;
            parity ^= (irrep_parity_matrix_get(M, i, j) & s_bit);
        }
        if (parity) out[i / 64] |= ((uint64_t)1 << (i % 64));
    }
    return IRREP_OK;
}

irrep_status_t
irrep_single_shot_meta_syndrome_X(const irrep_single_shot_code_t *c,
                                  const uint64_t *syndrome,
                                  uint64_t *meta_syndrome)
{
    if (c == NULL || syndrome == NULL || meta_syndrome == NULL)
        return IRREP_ERR_INVALID_ARG;
    return apply_meta(&c->M_X, syndrome, meta_syndrome);
}

irrep_status_t
irrep_single_shot_meta_syndrome_Z(const irrep_single_shot_code_t *c,
                                  const uint64_t *syndrome,
                                  uint64_t *meta_syndrome)
{
    if (c == NULL || syndrome == NULL || meta_syndrome == NULL)
        return IRREP_ERR_INVALID_ARG;
    return apply_meta(&c->M_Z, syndrome, meta_syndrome);
}

/* ====================================================================
 * Generic single-shot lifter: compute meta-check matrices from any CSS.
 *
 * Algorithm: F₂ left-nullspace of H via Gaussian elimination on the
 * augmented matrix `[H | I_m]`. After row reduction:
 *   - Rows whose H-portion has a pivot span the row-space of H.
 *   - Rows whose H-portion is all zero have their I-portion as a
 *     basis vector of the left-nullspace `{v : v · H = 0}`.
 *
 * Bit-packed XOR keeps the inner loop tight; the asymptotic cost is
 * `O(m · n · m / 64)`, fine for n, m up to a few thousand.
 * ==================================================================== */

static void
parity_row_xor(irrep_parity_matrix_t *m, int dst_row, int src_row)
{
    uint64_t *dst = &m->data[(size_t)dst_row * (size_t)m->n_words];
    const uint64_t *src = &m->data[(size_t)src_row * (size_t)m->n_words];
    for (int w = 0; w < m->n_words; ++w) dst[w] ^= src[w];
}

static void
parity_row_swap(irrep_parity_matrix_t *m, int row_a, int row_b)
{
    if (row_a == row_b) return;
    uint64_t *a = &m->data[(size_t)row_a * (size_t)m->n_words];
    uint64_t *b = &m->data[(size_t)row_b * (size_t)m->n_words];
    for (int w = 0; w < m->n_words; ++w) {
        uint64_t t = a[w]; a[w] = b[w]; b[w] = t;
    }
}

/* Compute basis of the left-nullspace of H (m × n matrix).
 *
 * Output `nullspace_out` is a `nullity × m` parity matrix; each row v
 * satisfies `v · H = 0` in F₂. Returns IRREP_OK and fills `nullspace_out`
 * even when nullity = 0 (degenerate empty matrix). */
static irrep_status_t
parity_matrix_left_nullspace(const irrep_parity_matrix_t *H,
                             irrep_parity_matrix_t *nullspace_out)
{
    int m = H->n_rows;
    int n = H->n_cols;

    /* Work-copy of H and an m × m identity that tracks row operations. */
    irrep_parity_matrix_t W, I;
    irrep_status_t s = irrep_parity_matrix_new(&W, m, n > 0 ? n : 1);
    if (s != IRREP_OK) return s;
    s = irrep_parity_matrix_new(&I, m, m > 0 ? m : 1);
    if (s != IRREP_OK) { irrep_parity_matrix_free(&W); return s; }

    /* Copy H into W bit-for-bit, then I = identity. */
    if (n > 0 && m > 0) {
        for (int r = 0; r < m; ++r) {
            for (int w = 0; w < H->n_words; ++w) {
                W.data[(size_t)r * (size_t)W.n_words + w] =
                    H->data[(size_t)r * (size_t)H->n_words + w];
            }
        }
    }
    for (int i = 0; i < m; ++i) {
        irrep_parity_matrix_set(&I, i, i);
    }

    /* Gaussian elimination: walk columns of W left-to-right; for each
     * column c, find a row r ≥ pivot_row with W[r][c] = 1, swap that
     * row into pivot_row, and XOR W[i] (and I[i]) by W[pivot_row] for
     * every other row i with W[i][c] = 1. */
    int pivot_row = 0;
    for (int c = 0; c < n && pivot_row < m; ++c) {
        int r = -1;
        for (int i = pivot_row; i < m; ++i) {
            if (irrep_parity_matrix_get(&W, i, c) == 1) { r = i; break; }
        }
        if (r == -1) continue;
        if (r != pivot_row) {
            parity_row_swap(&W, pivot_row, r);
            parity_row_swap(&I, pivot_row, r);
        }
        for (int i = 0; i < m; ++i) {
            if (i == pivot_row) continue;
            if (irrep_parity_matrix_get(&W, i, c) == 1) {
                parity_row_xor(&W, i, pivot_row);
                parity_row_xor(&I, i, pivot_row);
            }
        }
        ++pivot_row;
    }

    /* pivot_row = rank(H). Rows pivot_row..m-1 of I are the basis of
     * the left-nullspace of H. */
    int nullity = m - pivot_row;
    s = irrep_parity_matrix_new(nullspace_out, nullity > 0 ? nullity : 1, m);
    if (s != IRREP_OK) {
        irrep_parity_matrix_free(&W);
        irrep_parity_matrix_free(&I);
        return s;
    }
    /* If nullity == 0 the matrix was allocated with 1 row to satisfy the
     * non-zero-row precondition of irrep_parity_matrix_new; we still set
     * n_rows = 0 to honour the caller-visible row count. */
    nullspace_out->n_rows = nullity;
    for (int k = 0; k < nullity; ++k) {
        for (int w = 0; w < I.n_words; ++w) {
            nullspace_out->data[(size_t)k * (size_t)nullspace_out->n_words + w] =
                I.data[(size_t)(pivot_row + k) * (size_t)I.n_words + w];
        }
    }

    irrep_parity_matrix_free(&W);
    irrep_parity_matrix_free(&I);
    return IRREP_OK;
}

/* Copy `src` into a freshly-allocated parity matrix `dst`. */
static irrep_status_t
parity_matrix_copy(const irrep_parity_matrix_t *src, irrep_parity_matrix_t *dst)
{
    irrep_status_t s = irrep_parity_matrix_new(dst, src->n_rows > 0 ? src->n_rows : 1,
                                               src->n_cols > 0 ? src->n_cols : 1);
    if (s != IRREP_OK) return s;
    dst->n_rows = src->n_rows;
    dst->n_cols = src->n_cols;
    if (src->n_rows > 0 && src->n_cols > 0) {
        size_t bytes = (size_t)src->n_rows * (size_t)src->n_words * sizeof(uint64_t);
        memcpy(dst->data, src->data, bytes);
    }
    return IRREP_OK;
}

irrep_status_t
irrep_single_shot_lift(const irrep_css_code_t *css,
                       irrep_single_shot_code_t *out)
{
    if (css == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;

    /* Deep-copy the CSS code so the caller retains ownership of `css`. */
    out->css.n = css->n;
    irrep_status_t s = parity_matrix_copy(&css->H_X, &out->css.H_X);
    if (s != IRREP_OK) return s;
    s = parity_matrix_copy(&css->H_Z, &out->css.H_Z);
    if (s != IRREP_OK) {
        irrep_parity_matrix_free(&out->css.H_X);
        return s;
    }

    /* Compute meta-checks. */
    s = parity_matrix_left_nullspace(&css->H_X, &out->M_X);
    if (s != IRREP_OK) {
        irrep_css_code_free(&out->css);
        return s;
    }
    s = parity_matrix_left_nullspace(&css->H_Z, &out->M_Z);
    if (s != IRREP_OK) {
        irrep_parity_matrix_free(&out->M_X);
        irrep_css_code_free(&out->css);
        return s;
    }
    return IRREP_OK;
}
