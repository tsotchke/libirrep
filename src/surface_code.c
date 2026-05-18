/* SPDX-License-Identifier: MIT */
/** @file surface_code.c
 *  @brief Implementation of the rotated planar surface code [[d², 1, d]]. */
#include <irrep/surface_code.h>
#include <irrep/types.h>

#include <stddef.h>
#include <stdlib.h>

static inline int qcoord(int r, int c, int d) { return r * d + c; }

/* Count helpers — see header for the layout rules. */
static int count_bulk_X(int d) {
    int n = 0;
    for (int a = 0; a < d - 1; ++a)
        for (int b = 0; b < d - 1; ++b)
            if (((a + b) & 1) == 0) ++n;
    return n;
}
static int count_bulk_Z(int d) { return (d - 1) * (d - 1) - count_bulk_X(d); }
static int count_top_X(int d) {
    int n = 0;
    for (int b = 0; b < d - 1; ++b) if ((b & 1) == 1) ++n;
    return n;
}
static int count_bottom_X(int d) {
    int n = 0;
    for (int b = 0; b < d - 1; ++b) if (((d - 1 + b) & 1) == 0) ++n;
    return n;
}
static int count_left_Z(int d) {
    int n = 0;
    for (int a = 0; a < d - 1; ++a) if ((a & 1) == 0) ++n;
    return n;
}
static int count_right_Z(int d) {
    int n = 0;
    for (int a = 0; a < d - 1; ++a) if (((a + d - 1) & 1) == 1) ++n;
    return n;
}

irrep_status_t
irrep_surface_init(irrep_surface_params_t *out, int distance)
{
    if (out == NULL || distance < 2) return IRREP_ERR_INVALID_ARG;
    out->distance = distance;
    out->n_qubits = distance * distance;
    out->n_X_stabs = count_bulk_X(distance) + count_top_X(distance) + count_bottom_X(distance);
    out->n_Z_stabs = count_bulk_Z(distance) + count_left_Z(distance) + count_right_Z(distance);
    return IRREP_OK;
}

irrep_status_t
irrep_surface_build(const irrep_surface_params_t *p, irrep_css_code_t *out)
{
    if (p == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    int d = p->distance;
    if (d < 2) return IRREP_ERR_INVALID_ARG;

    irrep_status_t s = irrep_css_code_new(out, p->n_qubits, p->n_X_stabs, p->n_Z_stabs);
    if (s != IRREP_OK) return s;

    /* H_X rows: bulk-X (row-major (a,b)), then top-X, then bottom-X. */
    int xrow = 0;
    for (int a = 0; a < d - 1; ++a) {
        for (int b = 0; b < d - 1; ++b) {
            if (((a + b) & 1) == 0) {
                irrep_parity_matrix_set(&out->H_X, xrow, qcoord(a,     b,     d));
                irrep_parity_matrix_set(&out->H_X, xrow, qcoord(a + 1, b,     d));
                irrep_parity_matrix_set(&out->H_X, xrow, qcoord(a,     b + 1, d));
                irrep_parity_matrix_set(&out->H_X, xrow, qcoord(a + 1, b + 1, d));
                ++xrow;
            }
        }
    }
    for (int b = 0; b < d - 1; ++b) {
        if ((b & 1) == 1) {
            irrep_parity_matrix_set(&out->H_X, xrow, qcoord(0, b,     d));
            irrep_parity_matrix_set(&out->H_X, xrow, qcoord(0, b + 1, d));
            ++xrow;
        }
    }
    for (int b = 0; b < d - 1; ++b) {
        if (((d - 1 + b) & 1) == 0) {
            irrep_parity_matrix_set(&out->H_X, xrow, qcoord(d - 1, b,     d));
            irrep_parity_matrix_set(&out->H_X, xrow, qcoord(d - 1, b + 1, d));
            ++xrow;
        }
    }

    /* H_Z rows: bulk-Z, then left-Z, then right-Z. */
    int zrow = 0;
    for (int a = 0; a < d - 1; ++a) {
        for (int b = 0; b < d - 1; ++b) {
            if (((a + b) & 1) == 1) {
                irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a,     b,     d));
                irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a + 1, b,     d));
                irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a,     b + 1, d));
                irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a + 1, b + 1, d));
                ++zrow;
            }
        }
    }
    for (int a = 0; a < d - 1; ++a) {
        if ((a & 1) == 0) {
            irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a,     0, d));
            irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a + 1, 0, d));
            ++zrow;
        }
    }
    for (int a = 0; a < d - 1; ++a) {
        if (((a + d - 1) & 1) == 1) {
            irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a,     d - 1, d));
            irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a + 1, d - 1, d));
            ++zrow;
        }
    }

    return IRREP_OK;
}

/* ====================================================================
 * Logical operators
 * ==================================================================== */

irrep_status_t
irrep_surface_logical_X(const irrep_surface_params_t *p, irrep_pauli_t *out)
{
    if (p == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    int d = p->distance;
    if (d < 2) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_pauli_new(out, p->n_qubits);
    if (s != IRREP_OK) return s;
    /* L_X = X-string on column 0: qubits q(r, 0) = r*d for r ∈ [0, d). */
    for (int r = 0; r < d; ++r) {
        irrep_pauli_set(out, qcoord(r, 0, d), IRREP_PAULI_LETTER_X);
    }
    return IRREP_OK;
}

irrep_status_t
irrep_surface_logical_Z(const irrep_surface_params_t *p, irrep_pauli_t *out)
{
    if (p == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    int d = p->distance;
    if (d < 2) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_pauli_new(out, p->n_qubits);
    if (s != IRREP_OK) return s;
    /* L_Z = Z-string on row 0: qubits q(0, c) = c for c ∈ [0, d). */
    for (int c = 0; c < d; ++c) {
        irrep_pauli_set(out, qcoord(0, c, d), IRREP_PAULI_LETTER_Z);
    }
    return IRREP_OK;
}
