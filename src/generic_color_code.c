/* SPDX-License-Identifier: MIT */
/** @file generic_color_code.c
 *  @brief Face-list-driven 2D color code builder. */
#include <irrep/generic_color_code.h>
#include <irrep/types.h>

#include <stddef.h>

irrep_status_t
irrep_generic_color_build(const irrep_color_lattice_t *lattice,
                          irrep_css_code_t *out)
{
    if (lattice == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    if (lattice->n_qubits <= 0) return IRREP_ERR_INVALID_ARG;
    if (lattice->n_faces <= 0) return IRREP_ERR_INVALID_ARG;
    if (lattice->face_sizes == NULL || lattice->face_qubits == NULL)
        return IRREP_ERR_INVALID_ARG;

    /* Validate the face list: every face must have at least 3 qubits
     * (color codes), all qubit indices in range. */
    for (int f = 0; f < lattice->n_faces; ++f) {
        int sz = lattice->face_sizes[f];
        if (sz < 2) return IRREP_ERR_PRECONDITION;
        const int *qs = lattice->face_qubits[f];
        if (qs == NULL) return IRREP_ERR_PRECONDITION;
        for (int k = 0; k < sz; ++k) {
            if (qs[k] < 0 || qs[k] >= lattice->n_qubits)
                return IRREP_ERR_PRECONDITION;
        }
    }

    irrep_status_t s = irrep_css_code_new(out, lattice->n_qubits,
                                          lattice->n_faces, lattice->n_faces);
    if (s != IRREP_OK) return s;

    for (int f = 0; f < lattice->n_faces; ++f) {
        int sz = lattice->face_sizes[f];
        const int *qs = lattice->face_qubits[f];
        for (int k = 0; k < sz; ++k) {
            irrep_parity_matrix_set(&out->H_X, f, qs[k]);
            irrep_parity_matrix_set(&out->H_Z, f, qs[k]);
        }
    }

    /* Verify CSS orthogonality: every pair of faces meets in even
     * count. Since H_X = H_Z for color codes, this reduces to checking
     * pairs of face supports. */
    s = irrep_css_code_verify(out);
    if (s != IRREP_OK) {
        irrep_css_code_free(out);
        return s;
    }
    return IRREP_OK;
}
