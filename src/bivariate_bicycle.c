/* SPDX-License-Identifier: MIT */
/** @file bivariate_bicycle.c
 *  @brief Implementation of bivariate-bicycle (BB) code construction. */
#include <irrep/bivariate_bicycle.h>
#include <irrep/types.h>

#include <stdlib.h>
#include <string.h>

/* ====================================================================
 * Polynomials in F₂[x, y] / (xˡ - 1, yᵐ - 1)
 * ==================================================================== */

irrep_status_t
irrep_bb_poly_new(irrep_bb_poly_t *out, int ell, int m)
{
    if (out == NULL || ell <= 0 || m <= 0) return IRREP_ERR_INVALID_ARG;
    out->ell = ell;
    out->m = m;
    out->dim = ell * m;
    out->n_words = (out->dim + 63) / 64;
    out->coeffs = (uint64_t *)calloc(out->n_words, sizeof(uint64_t));
    if (out->coeffs == NULL) return IRREP_ERR_OUT_OF_MEMORY;
    return IRREP_OK;
}

void
irrep_bb_poly_free(irrep_bb_poly_t *p)
{
    if (p == NULL) return;
    free(p->coeffs); p->coeffs = NULL;
    p->ell = 0;
    p->m = 0;
    p->dim = 0;
    p->n_words = 0;
}

static inline int poly_index(const irrep_bb_poly_t *p, int a, int b)
{
    return a + p->ell * b;
}

irrep_status_t
irrep_bb_poly_add_monomial(irrep_bb_poly_t *p, int a, int b)
{
    if (p == NULL) return IRREP_ERR_INVALID_ARG;
    if (a < 0 || a >= p->ell) return IRREP_ERR_INVALID_ARG;
    if (b < 0 || b >= p->m)   return IRREP_ERR_INVALID_ARG;
    int idx = poly_index(p, a, b);
    int w = idx / 64;
    int bit = idx % 64;
    p->coeffs[w] ^= ((uint64_t)1 << bit);
    return IRREP_OK;
}

int
irrep_bb_poly_get(const irrep_bb_poly_t *p, int a, int b)
{
    if (p == NULL) return -1;
    if (a < 0 || a >= p->ell) return -1;
    if (b < 0 || b >= p->m)   return -1;
    int idx = poly_index(p, a, b);
    return (p->coeffs[idx / 64] >> (idx % 64)) & 1;
}

/* ====================================================================
 * BB code construction
 * ==================================================================== */

/* Set H[r, c] to the matrix-element of M_P at row r, column c, where
 * M_P encodes "left multiplication by P" in the basis {x^a y^b}.
 *
 * Identifying r = a_r + ℓ b_r and c = a_c + ℓ b_c, we have
 *   M_P[r, c] = P[(a_r - a_c) mod ℓ, (b_r - b_c) mod m].
 */
static int
matrix_element_M_P(const irrep_bb_poly_t *P, int r, int c)
{
    int ell = P->ell, m = P->m;
    int a_r = r % ell, b_r = r / ell;
    int a_c = c % ell, b_c = c / ell;
    int da = (a_r - a_c) % ell; if (da < 0) da += ell;
    int db = (b_r - b_c) % m;   if (db < 0) db += m;
    return irrep_bb_poly_get(P, da, db);
}

/* Set H[r, c] to the matrix-element of M_P^T = M_{Pᵀ} at row r, column c.
 *
 *   M_P^T[r, c] = M_P[c, r] = P[(a_c - a_r) mod ℓ, (b_c - b_r) mod m].
 */
static int
matrix_element_M_P_T(const irrep_bb_poly_t *P, int r, int c)
{
    return matrix_element_M_P(P, c, r);
}

irrep_status_t
irrep_bb_code_build(const irrep_bb_poly_t *A,
                    const irrep_bb_poly_t *B,
                    irrep_css_code_t *out)
{
    if (A == NULL || B == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    if (A->ell != B->ell || A->m != B->m) return IRREP_ERR_INVALID_ARG;

    int ell = A->ell;
    int m = A->m;
    int dim = ell * m;
    int n_qubits = 2 * dim;

    irrep_status_t s = irrep_css_code_new(out, n_qubits, dim, dim);
    if (s != IRREP_OK) return s;

    /* H_X = [A | B]: m_X = dim rows, n_qubits cols.
     * For row r: cols 0..dim-1 ← M_A row r; cols dim..2dim-1 ← M_B row r. */
    for (int r = 0; r < dim; ++r) {
        for (int c = 0; c < dim; ++c) {
            if (matrix_element_M_P(A, r, c)) {
                irrep_parity_matrix_set(&out->H_X, r, c);
            }
            if (matrix_element_M_P(B, r, c)) {
                irrep_parity_matrix_set(&out->H_X, r, dim + c);
            }
        }
    }

    /* H_Z = [B^T | A^T]: m_Z = dim rows, n_qubits cols.
     * For row r: cols 0..dim-1 ← M_B^T row r; cols dim..2dim-1 ← M_A^T row r. */
    for (int r = 0; r < dim; ++r) {
        for (int c = 0; c < dim; ++c) {
            if (matrix_element_M_P_T(B, r, c)) {
                irrep_parity_matrix_set(&out->H_Z, r, c);
            }
            if (matrix_element_M_P_T(A, r, c)) {
                irrep_parity_matrix_set(&out->H_Z, r, dim + c);
            }
        }
    }
    return IRREP_OK;
}

/* ====================================================================
 * Named IBM bivariate-bicycle instances (Bravyi et al. 2024, Table 3)
 * ==================================================================== */

irrep_status_t
irrep_bb_code_ibm_144_12_12(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_bb_poly_t A, B;
    irrep_status_t s = irrep_bb_poly_new(&A, /*ell=*/12, /*m=*/6);
    if (s != IRREP_OK) return s;
    s = irrep_bb_poly_new(&B, /*ell=*/12, /*m=*/6);
    if (s != IRREP_OK) { irrep_bb_poly_free(&A); return s; }

    /* A = x³ + y + y². */
    irrep_bb_poly_add_monomial(&A, /*a=*/3, /*b=*/0);
    irrep_bb_poly_add_monomial(&A, /*a=*/0, /*b=*/1);
    irrep_bb_poly_add_monomial(&A, /*a=*/0, /*b=*/2);
    /* B = y³ + x + x². */
    irrep_bb_poly_add_monomial(&B, /*a=*/0, /*b=*/3);
    irrep_bb_poly_add_monomial(&B, /*a=*/1, /*b=*/0);
    irrep_bb_poly_add_monomial(&B, /*a=*/2, /*b=*/0);

    s = irrep_bb_code_build(&A, &B, out);
    irrep_bb_poly_free(&A);
    irrep_bb_poly_free(&B);
    return s;
}

irrep_status_t
irrep_bb_code_ibm_288_12_18(irrep_css_code_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    irrep_bb_poly_t A, B;
    irrep_status_t s = irrep_bb_poly_new(&A, /*ell=*/12, /*m=*/12);
    if (s != IRREP_OK) return s;
    s = irrep_bb_poly_new(&B, /*ell=*/12, /*m=*/12);
    if (s != IRREP_OK) { irrep_bb_poly_free(&A); return s; }

    /* A = x³ + y² + y⁷. */
    irrep_bb_poly_add_monomial(&A, /*a=*/3, /*b=*/0);
    irrep_bb_poly_add_monomial(&A, /*a=*/0, /*b=*/2);
    irrep_bb_poly_add_monomial(&A, /*a=*/0, /*b=*/7);
    /* B = y³ + x + x². */
    irrep_bb_poly_add_monomial(&B, /*a=*/0, /*b=*/3);
    irrep_bb_poly_add_monomial(&B, /*a=*/1, /*b=*/0);
    irrep_bb_poly_add_monomial(&B, /*a=*/2, /*b=*/0);

    s = irrep_bb_code_build(&A, &B, out);
    irrep_bb_poly_free(&A);
    irrep_bb_poly_free(&B);
    return s;
}
