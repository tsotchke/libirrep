/* SPDX-License-Identifier: MIT */
/** @file lattice_surgery.c
 *  @brief Implementation of smooth-merge lattice surgery on rotated
 *         surface codes. The merged geometry is a `d × 2d` rectangular
 *         rotated surface code; the stabilizer rules generalise those
 *         of `<irrep/surface_code.h>` from square `d × d` to rectangular
 *         `n_rows × n_cols` without changing the parity assignment. */
#include <irrep/lattice_surgery.h>
#include <irrep/types.h>

#include <stddef.h>

static inline int qcoord(int r, int c, int n_cols) { return r * n_cols + c; }

/* Counts on a rectangular n_rows × n_cols rotated surface code with the
 * standard alternation rule (bulk X if a+b even, Z if odd; boundary
 * stabs at the four edges). */
static int count_bulk_X(int n_rows, int n_cols) {
    int n = 0;
    for (int a = 0; a < n_rows - 1; ++a)
        for (int b = 0; b < n_cols - 1; ++b)
            if (((a + b) & 1) == 0) ++n;
    return n;
}
static int count_bulk_Z(int n_rows, int n_cols) {
    return (n_rows - 1) * (n_cols - 1) - count_bulk_X(n_rows, n_cols);
}
static int count_top_X(int n_cols) {
    int n = 0;
    for (int b = 0; b < n_cols - 1; ++b) if ((b & 1) == 1) ++n;
    return n;
}
static int count_bottom_X(int n_rows, int n_cols) {
    int n = 0;
    for (int b = 0; b < n_cols - 1; ++b)
        if (((n_rows - 1 + b) & 1) == 0) ++n;
    return n;
}
static int count_left_Z(int n_rows) {
    int n = 0;
    for (int a = 0; a < n_rows - 1; ++a) if ((a & 1) == 0) ++n;
    return n;
}
static int count_right_Z(int n_rows, int n_cols) {
    int n = 0;
    for (int a = 0; a < n_rows - 1; ++a)
        if (((a + n_cols - 1) & 1) == 1) ++n;
    return n;
}

irrep_status_t
irrep_lattice_surgery_smooth_init(irrep_lattice_surgery_smooth_t *out, int d)
{
    if (out == NULL || d < 2) return IRREP_ERR_INVALID_ARG;
    out->d = d;
    out->rows = d;
    out->cols = 2 * d;
    out->n_qubits = 2 * d * d;
    out->n_X_stabs = count_bulk_X(d, 2 * d)
                   + count_top_X(2 * d)
                   + count_bottom_X(d, 2 * d);
    out->n_Z_stabs = count_bulk_Z(d, 2 * d)
                   + count_left_Z(d)
                   + count_right_Z(d, 2 * d);
    return IRREP_OK;
}

irrep_status_t
irrep_lattice_surgery_smooth_build(const irrep_lattice_surgery_smooth_t *p,
                                   irrep_css_code_t *out)
{
    if (p == NULL || out == NULL || p->d < 2) return IRREP_ERR_INVALID_ARG;
    const int n_rows = p->rows;
    const int n_cols = p->cols;

    irrep_status_t s = irrep_css_code_new(out, p->n_qubits,
                                          p->n_X_stabs, p->n_Z_stabs);
    if (s != IRREP_OK) return s;

    /* H_X rows: bulk-X (row-major (a,b)), then top-X (r=0), then bottom-X. */
    int xrow = 0;
    for (int a = 0; a < n_rows - 1; ++a) {
        for (int b = 0; b < n_cols - 1; ++b) {
            if (((a + b) & 1) == 0) {
                irrep_parity_matrix_set(&out->H_X, xrow, qcoord(a,     b,     n_cols));
                irrep_parity_matrix_set(&out->H_X, xrow, qcoord(a + 1, b,     n_cols));
                irrep_parity_matrix_set(&out->H_X, xrow, qcoord(a,     b + 1, n_cols));
                irrep_parity_matrix_set(&out->H_X, xrow, qcoord(a + 1, b + 1, n_cols));
                ++xrow;
            }
        }
    }
    for (int b = 0; b < n_cols - 1; ++b) {
        if ((b & 1) == 1) {
            irrep_parity_matrix_set(&out->H_X, xrow, qcoord(0, b,     n_cols));
            irrep_parity_matrix_set(&out->H_X, xrow, qcoord(0, b + 1, n_cols));
            ++xrow;
        }
    }
    for (int b = 0; b < n_cols - 1; ++b) {
        if (((n_rows - 1 + b) & 1) == 0) {
            irrep_parity_matrix_set(&out->H_X, xrow, qcoord(n_rows - 1, b,     n_cols));
            irrep_parity_matrix_set(&out->H_X, xrow, qcoord(n_rows - 1, b + 1, n_cols));
            ++xrow;
        }
    }

    /* H_Z rows: bulk-Z, then left-Z, then right-Z. */
    int zrow = 0;
    for (int a = 0; a < n_rows - 1; ++a) {
        for (int b = 0; b < n_cols - 1; ++b) {
            if (((a + b) & 1) == 1) {
                irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a,     b,     n_cols));
                irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a + 1, b,     n_cols));
                irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a,     b + 1, n_cols));
                irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a + 1, b + 1, n_cols));
                ++zrow;
            }
        }
    }
    for (int a = 0; a < n_rows - 1; ++a) {
        if ((a & 1) == 0) {
            irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a,     0, n_cols));
            irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a + 1, 0, n_cols));
            ++zrow;
        }
    }
    for (int a = 0; a < n_rows - 1; ++a) {
        if (((a + n_cols - 1) & 1) == 1) {
            irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a,     n_cols - 1, n_cols));
            irrep_parity_matrix_set(&out->H_Z, zrow, qcoord(a + 1, n_cols - 1, n_cols));
            ++zrow;
        }
    }

    return IRREP_OK;
}

irrep_status_t
irrep_lattice_surgery_smooth_logical_X(const irrep_lattice_surgery_smooth_t *p,
                                       irrep_pauli_t *out)
{
    if (p == NULL || out == NULL || p->d < 2) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_pauli_new(out, p->n_qubits);
    if (s != IRREP_OK) return s;
    /* L̄_X(M) = X-string on column 0 spanning all d rows (between the
     * top and bottom rough X-boundaries). Matches `irrep_surface_logical_X`. */
    for (int r = 0; r < p->rows; ++r) {
        irrep_pauli_set(out, qcoord(r, 0, p->cols), IRREP_PAULI_LETTER_X);
    }
    return IRREP_OK;
}

irrep_status_t
irrep_lattice_surgery_smooth_logical_Z(const irrep_lattice_surgery_smooth_t *p,
                                       irrep_pauli_t *out)
{
    if (p == NULL || out == NULL || p->d < 2) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_pauli_new(out, p->n_qubits);
    if (s != IRREP_OK) return s;
    /* L̄_Z(M) = Z-string on row 0 spanning all 2d columns (between the
     * left and right smooth Z-boundaries — now `2d` apart after merge).
     * This equals `Z̄_A · Z̄_B`: the smooth merge collapses both logical-Z
     * operators into the merged L̄_Z. */
    for (int c = 0; c < p->cols; ++c) {
        irrep_pauli_set(out, qcoord(0, c, p->cols), IRREP_PAULI_LETTER_Z);
    }
    return IRREP_OK;
}

irrep_status_t
irrep_lattice_surgery_smooth_joint_X_parity(
    const irrep_lattice_surgery_smooth_t *p, irrep_pauli_t *out)
{
    if (p == NULL || out == NULL || p->d < 2) return IRREP_ERR_INVALID_ARG;
    irrep_status_t s = irrep_pauli_new(out, p->n_qubits);
    if (s != IRREP_OK) return s;
    /* P_joint = X̄_A · X̄_B = X on column 0 AND column d of the merged frame.
     * Pre-merge: column 0 is A's logical-X column (X̄_A), column d is B's
     * logical-X column (X̄_B). After smooth merge their product becomes a
     * stabilizer — the joint-X parity measurement. */
    const int d = p->d;
    const int n_cols = p->cols;
    for (int r = 0; r < d; ++r) {
        irrep_pauli_set(out, qcoord(r, 0, n_cols), IRREP_PAULI_LETTER_X);
        irrep_pauli_set(out, qcoord(r, d, n_cols), IRREP_PAULI_LETTER_X);
    }
    return IRREP_OK;
}
