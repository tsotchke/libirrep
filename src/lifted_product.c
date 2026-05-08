/* SPDX-License-Identifier: MIT */
/** @file lifted_product.c
 *  @brief Implementation of Panteleev-Kalachev lifted-product CSS qLDPC codes.
 *
 *  Construction (abelian, ring R = F₂[ℤ/ℓ] = F₂[x]/(xˡ-1)):
 *
 *      H_X^R = [ A ⊗ I_{n_b}    |  I_{m_a} ⊗ B^T ]
 *      H_Z^R = [ I_{n_a} ⊗ σ(B) |  σ(A^T) ⊗ I_{m_b} ]
 *
 *  where σ is the antipode of the cyclic group, σ(x^a) = x^{-a}; on
 *  coefficient vectors, σ(p)[a] = p[(ℓ - a) mod ℓ]. The σ in `H_Z^R`
 *  is what makes CSS orthogonality at the F₂ level hold:
 *
 *      L(H_X^R) · L(H_Z^R)^T  =  H_X^R · σ(H_Z^R)^T (lifted)
 *                              =  A⊗B^T + A⊗B^T
 *                              =  0  (mod 2).
 *
 *  ## Why σ
 *
 *  The lift `R → F₂^{ℓ × ℓ}` sends `p(x) ↦ M_p` where
 *  `M_p[r, c] = p[(r - c) mod ℓ]`. Importantly, the F₂ transpose of
 *  the lift equals the lift of the σ-conjugate: `(M_p)^T = M_{σ(p)}`.
 *  Therefore taking F₂-transpose at the lifted level corresponds to
 *  applying σ at the polynomial level. The CSS condition
 *  `L(H_X) · L(H_Z)^T = 0` translates into `Σ_j H_X[i, j] · σ(H_Z[i', j]) = 0`
 *  in R, which the formula above satisfies.
 */
#include <irrep/lifted_product.h>
#include <irrep/types.h>

#include <stdlib.h>
#include <string.h>

/* ====================================================================
 * Polynomial-matrix bookkeeping
 * ==================================================================== */

irrep_status_t
irrep_poly_matrix_new(irrep_poly_matrix_t *out,
                      int n_rows, int n_cols, int ell)
{
    if (out == NULL || n_rows <= 0 || n_cols <= 0 || ell <= 0) return IRREP_ERR_INVALID_ARG;
    out->n_rows = n_rows;
    out->n_cols = n_cols;
    out->ell = ell;
    out->n_words = (ell + 63) / 64;
    size_t total = (size_t)n_rows * (size_t)n_cols * (size_t)out->n_words;
    out->data = (uint64_t *)calloc(total, sizeof(uint64_t));
    if (out->data == NULL) return IRREP_ERR_OUT_OF_MEMORY;
    return IRREP_OK;
}

void
irrep_poly_matrix_free(irrep_poly_matrix_t *m)
{
    if (m == NULL) return;
    free(m->data); m->data = NULL;
    m->n_rows = m->n_cols = m->ell = m->n_words = 0;
}

static inline size_t pm_offset(const irrep_poly_matrix_t *m, int row, int col)
{
    /* size_t avoids overflow at research-scale lifts (n_rows · n_cols · n_words
     * easily exceeds INT32_MAX for ℓ ≳ 100, n_a ≳ 100). */
    return ((size_t)row * (size_t)m->n_cols + (size_t)col) * (size_t)m->n_words;
}

irrep_status_t
irrep_poly_matrix_add_monomial(irrep_poly_matrix_t *m,
                               int row, int col, int a)
{
    if (m == NULL) return IRREP_ERR_INVALID_ARG;
    if (row < 0 || row >= m->n_rows) return IRREP_ERR_INVALID_ARG;
    if (col < 0 || col >= m->n_cols) return IRREP_ERR_INVALID_ARG;
    if (a < 0 || a >= m->ell) return IRREP_ERR_INVALID_ARG;
    int base = pm_offset(m, row, col);
    m->data[base + a / 64] ^= ((uint64_t)1 << (a % 64));
    return IRREP_OK;
}

int
irrep_poly_matrix_get(const irrep_poly_matrix_t *m,
                      int row, int col, int a)
{
    if (m == NULL) return -1;
    if (row < 0 || row >= m->n_rows) return -1;
    if (col < 0 || col >= m->n_cols) return -1;
    if (a < 0 || a >= m->ell) return -1;
    int base = pm_offset(m, row, col);
    return (m->data[base + a / 64] >> (a % 64)) & 1;
}

/* ====================================================================
 * Lifted-product CSS construction
 * ==================================================================== */

/* Set the (ℓ × ℓ) circulant block for polynomial p at row-block i_R, col-block j_R.
 * If sigma == 0: M_p[r, c] = p[(r - c) mod ℓ].
 * If sigma == 1: M_{σ(p)}[r, c] = σ(p)[(r - c) mod ℓ] = p[((c - r) mod ℓ)]. */
static void
write_circulant_block(irrep_parity_matrix_t *H,
                      const irrep_poly_matrix_t *P,
                      int p_row, int p_col,
                      int i_R, int j_R,
                      int sigma)
{
    int ell = P->ell;
    for (int r = 0; r < ell; ++r) {
        for (int c = 0; c < ell; ++c) {
            int a = sigma ? ((c - r) % ell + ell) % ell
                          : ((r - c) % ell + ell) % ell;
            if (irrep_poly_matrix_get(P, p_row, p_col, a)) {
                irrep_parity_matrix_set(H, ell * i_R + r, ell * j_R + c);
            }
        }
    }
}

irrep_status_t
irrep_lifted_product_build(const irrep_poly_matrix_t *A,
                           const irrep_poly_matrix_t *B,
                           irrep_css_code_t *out)
{
    if (A == NULL || B == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    if (A->ell != B->ell) return IRREP_ERR_INVALID_ARG;

    int ell = A->ell;
    int m_a = A->n_rows, n_a = A->n_cols;
    int m_b = B->n_rows, n_b = B->n_cols;
    int n_qubits = ell * (n_a * n_b + m_a * m_b);
    int m_X = ell * m_a * n_b;
    int m_Z = ell * n_a * m_b;

    irrep_status_t s = irrep_css_code_new(out, n_qubits, m_X, m_Z);
    if (s != IRREP_OK) return s;

    /* H_X^R = [A ⊗ I_{n_b} | I_{m_a} ⊗ B^T]
     *
     * Row index in H_X^R: i = i_a · n_b + i_b.
     * Column in left block: j = j_a · n_b + j_b. Entry = A[i_a, j_a] if i_b == j_b.
     * Column in right block: k = k_a · m_b + k_b. Entry = (B^T)[i_b, k_b] = B[k_b, i_b]
     *                                                    if i_a == k_a.
     */
    for (int i_a = 0; i_a < m_a; ++i_a) {
        for (int i_b = 0; i_b < n_b; ++i_b) {
            int i_R = i_a * n_b + i_b;
            for (int j_a = 0; j_a < n_a; ++j_a) {
                int j_R = j_a * n_b + i_b;
                write_circulant_block(&out->H_X, A, i_a, j_a, i_R, j_R, /*sigma=*/0);
            }
            for (int k_b = 0; k_b < m_b; ++k_b) {
                int k_R = (n_a * n_b) + (i_a * m_b + k_b);
                /* B[k_b, i_b] (no σ; the σ goes into H_Z). */
                write_circulant_block(&out->H_X, B, k_b, i_b, i_R, k_R, /*sigma=*/0);
            }
        }
    }

    /* H_Z^R = [I_{n_a} ⊗ σ(B) | σ(A^T) ⊗ I_{m_b}]
     *
     * Row index in H_Z^R: i = i_a · m_b + i_b.
     * Column in left block: j = j_a · n_b + j_b. Entry = σ(B[i_b, j_b]) if i_a == j_a.
     * Column in right block: k = k_a · m_b + k_b. Entry = σ(A^T)[i_a, k_a] = σ(A[k_a, i_a])
     *                                                    if i_b == k_b.
     */
    for (int i_a = 0; i_a < n_a; ++i_a) {
        for (int i_b = 0; i_b < m_b; ++i_b) {
            int i_R = i_a * m_b + i_b;
            for (int j_b = 0; j_b < n_b; ++j_b) {
                int j_R = i_a * n_b + j_b;
                /* σ(B[i_b, j_b]). */
                write_circulant_block(&out->H_Z, B, i_b, j_b, i_R, j_R, /*sigma=*/1);
            }
            for (int k_a = 0; k_a < m_a; ++k_a) {
                int k_R = (n_a * n_b) + (k_a * m_b + i_b);
                /* σ(A[k_a, i_a]). */
                write_circulant_block(&out->H_Z, A, k_a, i_a, i_R, k_R, /*sigma=*/1);
            }
        }
    }
    return IRREP_OK;
}
