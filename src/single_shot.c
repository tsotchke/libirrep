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
