/* SPDX-License-Identifier: MIT */
/** @file hypergraph_product.c
 *  @brief Implementation of the Tillich-Zemor hypergraph-product CSS code. */
#include <irrep/hypergraph_product.h>
#include <irrep/types.h>

#include <stdlib.h>

/* For Kronecker-product encoding into the CSS check matrices, recall
 *
 *   (H_a ⊗ I_{n_2})[i, j] = H_a[i / n_2, j / n_2] · [i mod n_2 == j mod n_2]
 *   (I_{m_1} ⊗ H_b^T)[i, j] = [i / n_2 == j / m_2] · H_b^T[i mod n_2, j mod m_2]
 *                            = [i / n_2 == j / m_2] · H_b[j mod m_2, i mod n_2]
 *
 * etc. The two left blocks have row-index space `i ∈ [0, m_1·n_2)`. */

irrep_status_t
irrep_hypergraph_product_build(const irrep_parity_matrix_t *H_a,
                               const irrep_parity_matrix_t *H_b,
                               irrep_css_code_t *out)
{
    if (H_a == NULL || H_b == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    int m1 = H_a->n_rows, n1 = H_a->n_cols;
    int m2 = H_b->n_rows, n2 = H_b->n_cols;
    int n_q = n1 * n2 + m1 * m2;
    int m_X = m1 * n2;
    int m_Z = n1 * m2;

    irrep_status_t s = irrep_css_code_new(out, n_q, m_X, m_Z);
    if (s != IRREP_OK) return s;

    /* H_X = [ H_a ⊗ I_{n_2}  |  I_{m_1} ⊗ H_b^T ]
     *
     * Row index i ∈ [0, m_1 n_2): write i = i_a · n_2 + i_b with i_a ∈ [0, m_1), i_b ∈ [0, n_2).
     * Column index j ∈ [0, n_q):
     *   - j ∈ [0, n_1 n_2): write j = j_a · n_2 + j_b. Set bit if H_a[i_a, j_a] AND i_b == j_b.
     *   - j ∈ [n_1 n_2, n_q): k = j - n_1 n_2; write k = k_a · m_2 + k_b with k_a ∈ [0, m_1), k_b ∈ [0, m_2).
     *     Set bit if i_a == k_a AND H_b[k_b, i_b]. */
    for (int i_a = 0; i_a < m1; ++i_a) {
        for (int i_b = 0; i_b < n2; ++i_b) {
            int row = i_a * n2 + i_b;
            /* Left block: (H_a ⊗ I_{n_2}) — fix j_b = i_b; vary j_a where H_a[i_a, j_a] = 1. */
            for (int j_a = 0; j_a < n1; ++j_a) {
                if (irrep_parity_matrix_get(H_a, i_a, j_a)) {
                    int col = j_a * n2 + i_b;
                    irrep_parity_matrix_set(&out->H_X, row, col);
                }
            }
            /* Right block: (I_{m_1} ⊗ H_b^T) — fix k_a = i_a; vary k_b where H_b[k_b, i_b] = 1. */
            for (int k_b = 0; k_b < m2; ++k_b) {
                if (irrep_parity_matrix_get(H_b, k_b, i_b)) {
                    int col = (n1 * n2) + (i_a * m2 + k_b);
                    irrep_parity_matrix_set(&out->H_X, row, col);
                }
            }
        }
    }

    /* H_Z = [ I_{n_1} ⊗ H_b  |  H_a^T ⊗ I_{m_2} ]
     *
     * Row index i ∈ [0, n_1 m_2): write i = i_a · m_2 + i_b with i_a ∈ [0, n_1), i_b ∈ [0, m_2).
     * Column index j:
     *   - j ∈ [0, n_1 n_2): j = j_a · n_2 + j_b. Set bit if i_a == j_a AND H_b[i_b, j_b].
     *   - j ∈ [n_1 n_2, n_q): k = j - n_1 n_2 = k_a · m_2 + k_b with k_a ∈ [0, m_1), k_b ∈ [0, m_2).
     *     Set bit if H_a^T[i_a, k_a] (= H_a[k_a, i_a]) AND k_b == i_b. */
    for (int i_a = 0; i_a < n1; ++i_a) {
        for (int i_b = 0; i_b < m2; ++i_b) {
            int row = i_a * m2 + i_b;
            for (int j_b = 0; j_b < n2; ++j_b) {
                if (irrep_parity_matrix_get(H_b, i_b, j_b)) {
                    int col = i_a * n2 + j_b;
                    irrep_parity_matrix_set(&out->H_Z, row, col);
                }
            }
            for (int k_a = 0; k_a < m1; ++k_a) {
                if (irrep_parity_matrix_get(H_a, k_a, i_a)) {
                    int col = (n1 * n2) + (k_a * m2 + i_b);
                    irrep_parity_matrix_set(&out->H_Z, row, col);
                }
            }
        }
    }
    return IRREP_OK;
}

/* ====================================================================
 * Named instances
 * ==================================================================== */

irrep_status_t
irrep_hgp_repetition_3_13_1_3(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;

    /* Build the 2×3 parity-check of the [3, 1, 3] repetition code:
     *   H = [[1, 1, 0],
     *        [0, 1, 1]]                                              */
    irrep_parity_matrix_t H;
    irrep_status_t s = irrep_parity_matrix_new(&H, /*rows=*/2, /*cols=*/3);
    if (s != IRREP_OK) return s;
    irrep_parity_matrix_set(&H, 0, 0);
    irrep_parity_matrix_set(&H, 0, 1);
    irrep_parity_matrix_set(&H, 1, 1);
    irrep_parity_matrix_set(&H, 1, 2);

    /* HGP(H, H). */
    s = irrep_hypergraph_product_build(&H, &H, out);
    irrep_parity_matrix_free(&H);
    return s;
}

irrep_status_t
irrep_hgp_repetition_4_25_1_4(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;

    /* 3×4 parity check of the [4, 1, 4] repetition code:
     *   H = [[1, 1, 0, 0],
     *        [0, 1, 1, 0],
     *        [0, 0, 1, 1]]                                            */
    irrep_parity_matrix_t H;
    irrep_status_t s = irrep_parity_matrix_new(&H, /*rows=*/3, /*cols=*/4);
    if (s != IRREP_OK) return s;
    irrep_parity_matrix_set(&H, 0, 0); irrep_parity_matrix_set(&H, 0, 1);
    irrep_parity_matrix_set(&H, 1, 1); irrep_parity_matrix_set(&H, 1, 2);
    irrep_parity_matrix_set(&H, 2, 2); irrep_parity_matrix_set(&H, 2, 3);

    s = irrep_hypergraph_product_build(&H, &H, out);
    irrep_parity_matrix_free(&H);
    return s;
}

irrep_status_t
irrep_hgp_repetition_5_41_1_5(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_parity_matrix_t H;
    irrep_status_t s = irrep_parity_matrix_new(&H, /*rows=*/4, /*cols=*/5);
    if (s != IRREP_OK) return s;
    irrep_parity_matrix_set(&H, 0, 0); irrep_parity_matrix_set(&H, 0, 1);
    irrep_parity_matrix_set(&H, 1, 1); irrep_parity_matrix_set(&H, 1, 2);
    irrep_parity_matrix_set(&H, 2, 2); irrep_parity_matrix_set(&H, 2, 3);
    irrep_parity_matrix_set(&H, 3, 3); irrep_parity_matrix_set(&H, 3, 4);
    s = irrep_hypergraph_product_build(&H, &H, out);
    irrep_parity_matrix_free(&H);
    return s;
}
