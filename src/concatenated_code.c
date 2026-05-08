/* SPDX-License-Identifier: MIT */
/** @file concatenated_code.c
 *  @brief Concatenated CSS code construction. */
#include <irrep/concatenated_code.h>
#include <irrep/types.h>

#include <stddef.h>

irrep_status_t
irrep_css_concatenate(const irrep_css_code_t *inner,
                      const irrep_pauli_t *inner_x_bar,
                      const irrep_pauli_t *inner_z_bar,
                      const irrep_css_code_t *outer,
                      irrep_css_code_t *out)
{
    if (inner == NULL || inner_x_bar == NULL || inner_z_bar == NULL ||
        outer == NULL || out == NULL) {
        return IRREP_ERR_INVALID_ARG;
    }
    if (inner->n != inner_x_bar->n || inner->n != inner_z_bar->n) {
        return IRREP_ERR_INVALID_ARG;
    }

    /* Pure-X / pure-Z preconditions on the inner logicals. */
    for (int w = 0; w < inner_x_bar->n_words; ++w) {
        if (inner_x_bar->z[w] != 0) return IRREP_ERR_PRECONDITION;
    }
    for (int w = 0; w < inner_z_bar->n_words; ++w) {
        if (inner_z_bar->x[w] != 0) return IRREP_ERR_PRECONDITION;
    }

    int n_in = inner->n;
    int n_out_q = outer->n;
    int n_total = n_in * n_out_q;

    int m_X_in = inner->H_X.n_rows;
    int m_Z_in = inner->H_Z.n_rows;
    int m_X_out = outer->H_X.n_rows;
    int m_Z_out = outer->H_Z.n_rows;
    int total_X = n_out_q * m_X_in + m_X_out;
    int total_Z = n_out_q * m_Z_in + m_Z_out;

    irrep_status_t s = irrep_css_code_new(out, n_total, total_X, total_Z);
    if (s != IRREP_OK) return s;

    /* H_X: blockwise inner X-stabs first, then outer X-stabs lifted via X̄. */
    int row = 0;
    for (int b = 0; b < n_out_q; ++b) {
        int offset = b * n_in;
        for (int i = 0; i < m_X_in; ++i, ++row) {
            for (int j = 0; j < n_in; ++j) {
                if (irrep_parity_matrix_get(&inner->H_X, i, j)) {
                    irrep_parity_matrix_set(&out->H_X, row, offset + j);
                }
            }
        }
    }
    for (int i = 0; i < m_X_out; ++i, ++row) {
        for (int q = 0; q < n_out_q; ++q) {
            if (!irrep_parity_matrix_get(&outer->H_X, i, q)) continue;
            int offset = q * n_in;
            for (int j = 0; j < n_in; ++j) {
                if ((inner_x_bar->x[j / 64] >> (j % 64)) & 1) {
                    irrep_parity_matrix_set(&out->H_X, row, offset + j);
                }
            }
        }
    }

    /* H_Z: blockwise inner Z-stabs, then outer Z-stabs lifted via Z̄. */
    row = 0;
    for (int b = 0; b < n_out_q; ++b) {
        int offset = b * n_in;
        for (int i = 0; i < m_Z_in; ++i, ++row) {
            for (int j = 0; j < n_in; ++j) {
                if (irrep_parity_matrix_get(&inner->H_Z, i, j)) {
                    irrep_parity_matrix_set(&out->H_Z, row, offset + j);
                }
            }
        }
    }
    for (int i = 0; i < m_Z_out; ++i, ++row) {
        for (int q = 0; q < n_out_q; ++q) {
            if (!irrep_parity_matrix_get(&outer->H_Z, i, q)) continue;
            int offset = q * n_in;
            for (int j = 0; j < n_in; ++j) {
                if ((inner_z_bar->z[j / 64] >> (j % 64)) & 1) {
                    irrep_parity_matrix_set(&out->H_Z, row, offset + j);
                }
            }
        }
    }

    return IRREP_OK;
}
