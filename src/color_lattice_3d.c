/* SPDX-License-Identifier: MIT */
/** @file color_lattice_3d.c
 *  @brief Generic 3-colex framework + tetrahedral instance. */
#include <irrep/color_lattice_3d.h>
#include <irrep/css_code.h>
#include <irrep/types.h>

#include <stddef.h>

/* ====================================================================
 * Generic builder: 3-colex face-and-cell list → CSS code
 * ==================================================================== */

irrep_status_t
irrep_color_lattice_3d_to_css(const irrep_color_lattice_3d_t *lattice,
                              irrep_css_code_t *out)
{
    if (lattice == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    if (lattice->n_qubits <= 0 || lattice->n_faces < 0 || lattice->n_cells < 0)
        return IRREP_ERR_INVALID_ARG;
    if ((lattice->n_faces > 0 &&
         (lattice->face_sizes == NULL || lattice->face_qubits == NULL)) ||
        (lattice->n_cells > 0 &&
         (lattice->cell_sizes == NULL || lattice->cell_qubits == NULL)))
        return IRREP_ERR_INVALID_ARG;

    irrep_status_t s = irrep_css_code_new(out, lattice->n_qubits,
                                          lattice->n_faces, lattice->n_cells);
    if (s != IRREP_OK) return s;

    /* H_X rows: one per face. */
    for (int f = 0; f < lattice->n_faces; ++f) {
        int sz = lattice->face_sizes[f];
        if (sz <= 0 || lattice->face_qubits[f] == NULL) {
            irrep_css_code_free(out);
            return IRREP_ERR_PRECONDITION;
        }
        for (int k = 0; k < sz; ++k) {
            int q = lattice->face_qubits[f][k];
            if (q < 0 || q >= lattice->n_qubits) {
                irrep_css_code_free(out);
                return IRREP_ERR_PRECONDITION;
            }
            irrep_parity_matrix_set(&out->H_X, f, q);
        }
    }
    /* H_Z rows: one per cell. */
    for (int c = 0; c < lattice->n_cells; ++c) {
        int sz = lattice->cell_sizes[c];
        if (sz <= 0 || lattice->cell_qubits[c] == NULL) {
            irrep_css_code_free(out);
            return IRREP_ERR_PRECONDITION;
        }
        for (int k = 0; k < sz; ++k) {
            int q = lattice->cell_qubits[c][k];
            if (q < 0 || q >= lattice->n_qubits) {
                irrep_css_code_free(out);
                return IRREP_ERR_PRECONDITION;
            }
            irrep_parity_matrix_set(&out->H_Z, c, q);
        }
    }
    /* Verify CSS orthogonality: H_X · H_Z^T = 0 (mod 2). */
    if (irrep_css_code_verify(out) != IRREP_OK) {
        irrep_css_code_free(out);
        return IRREP_ERR_PRECONDITION;
    }
    return IRREP_OK;
}

/* ====================================================================
 * Tetrahedral 3-colex — the [[15, 1, 3]] worked example
 *
 * Static face-and-cell tables. Qubit indexing: qubit q = 0..14
 * corresponds to PG(3, 2) point (q + 1), i.e. binary 4-tuple (q + 1).
 *
 * The face list contains the 4 "single-bit" supports (positions where
 * a specified bit of (q + 1) is set). The cell list contains the same
 * 4 single-bit supports (C_X ⊂ C_Z in the punctured-RM construction)
 * plus 6 "pair" supports (positions where two specified bits of
 * (q + 1) are both set).
 * ==================================================================== */

/* Bit-0 set in (q+1): (q+1) odd → q even, 8 qubits. */
static const int kFaceBit0[8] = { 0, 2, 4, 6, 8, 10, 12, 14 };
/* Bit-1 set in (q+1): 8 qubits. */
static const int kFaceBit1[8] = { 1, 2, 5, 6, 9, 10, 13, 14 };
/* Bit-2 set in (q+1): 8 qubits. */
static const int kFaceBit2[8] = { 3, 4, 5, 6, 11, 12, 13, 14 };
/* Bit-3 set in (q+1): 8 qubits. */
static const int kFaceBit3[8] = { 7, 8, 9, 10, 11, 12, 13, 14 };

/* Bits (0, 1) both set in (q+1): 4 qubits. */
static const int kCellPair01[4] = { 2, 6, 10, 14 };
/* Bits (0, 2): 4 qubits. */
static const int kCellPair02[4] = { 4, 6, 12, 14 };
/* Bits (0, 3): 4 qubits. */
static const int kCellPair03[4] = { 8, 10, 12, 14 };
/* Bits (1, 2): 4 qubits. */
static const int kCellPair12[4] = { 5, 6, 13, 14 };
/* Bits (1, 3): 4 qubits. */
static const int kCellPair13[4] = { 9, 10, 13, 14 };
/* Bits (2, 3): 4 qubits. */
static const int kCellPair23[4] = { 11, 12, 13, 14 };

static const int kTetraFaceSizes[4] = { 8, 8, 8, 8 };
static const int * const kTetraFaceQubits[4] = {
    kFaceBit0, kFaceBit1, kFaceBit2, kFaceBit3,
};

/* 10 cells: 4 single-bit (weight 8, same supports as the faces) + 6
 * pair (weight 4). */
static const int kTetraCellSizes[10] = { 8, 8, 8, 8, 4, 4, 4, 4, 4, 4 };
static const int * const kTetraCellQubits[10] = {
    kFaceBit0, kFaceBit1, kFaceBit2, kFaceBit3,
    kCellPair01, kCellPair02, kCellPair03,
    kCellPair12, kCellPair13, kCellPair23,
};

irrep_status_t
irrep_color_lattice_3d_tetrahedron(irrep_color_lattice_3d_t *out)
{
    if (out == NULL) return IRREP_ERR_INVALID_ARG;
    out->n_qubits = 15;
    out->n_faces  = 4;
    out->face_sizes  = kTetraFaceSizes;
    out->face_qubits = kTetraFaceQubits;
    out->n_cells  = 10;
    out->cell_sizes  = kTetraCellSizes;
    out->cell_qubits = kTetraCellQubits;
    return IRREP_OK;
}

void
irrep_color_lattice_3d_tetrahedron_free(irrep_color_lattice_3d_t *lattice)
{
    /* All storage is static-const; nothing to free. Zero the struct as
     * a courtesy so use-after-free becomes a NULL-deref. */
    if (lattice == NULL) return;
    lattice->n_qubits = 0;
    lattice->n_faces = 0;
    lattice->n_cells = 0;
    lattice->face_sizes = NULL;
    lattice->face_qubits = NULL;
    lattice->cell_sizes = NULL;
    lattice->cell_qubits = NULL;
}
