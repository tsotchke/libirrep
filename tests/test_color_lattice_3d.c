/* SPDX-License-Identifier: MIT */
/* Tests for `<irrep/color_lattice_3d.h>` — generic 3-colex framework
 * + tetrahedral instance.
 *
 * Verifies:
 *   - The framework reconstructs the [[15, 1, 3]] code from the
 *     tetrahedral face-and-cell list, with the same H_X and H_Z
 *     supports as the direct construction in
 *     `irrep_color_3d_rm_15_1_3`.
 *   - CSS orthogonality and full pairwise commutativity hold on the
 *     framework-built code.
 *   - Brute-force distance verification finds d = 3.
 *   - The framework rejects malformed inputs (out-of-range qubit
 *     indices, NULL face/cell lists).
 */
#include "harness.h"
#include <irrep/color_code_3d.h>
#include <irrep/color_lattice_3d.h>
#include <irrep/css_code.h>
#include <irrep/qec_distance.h>
#include <irrep/stabilizer_group.h>

#include <stdio.h>

/* Compare two parity matrices' supports cell-by-cell. */
static int parity_matrices_equal(const irrep_parity_matrix_t *a,
                                 const irrep_parity_matrix_t *b,
                                 int n_rows, int n_cols)
{
    for (int r = 0; r < n_rows; ++r) {
        for (int q = 0; q < n_cols; ++q) {
            if (irrep_parity_matrix_get(a, r, q) !=
                irrep_parity_matrix_get(b, r, q)) {
                return 0;
            }
        }
    }
    return 1;
}

static int test_tetrahedron_matches_direct_construction(void) {
    IRREP_TEST_START("color_lattice_3d_tetrahedron_matches_direct");

    /* Build via the framework. */
    irrep_color_lattice_3d_t lat;
    IRREP_ASSERT(irrep_color_lattice_3d_tetrahedron(&lat) == IRREP_OK);
    IRREP_ASSERT(lat.n_qubits == 15);
    IRREP_ASSERT(lat.n_faces == 4);
    IRREP_ASSERT(lat.n_cells == 10);

    irrep_css_code_t cs_framework;
    IRREP_ASSERT(irrep_color_lattice_3d_to_css(&lat, &cs_framework) == IRREP_OK);
    IRREP_ASSERT(cs_framework.n == 15);
    IRREP_ASSERT(cs_framework.H_X.n_rows == 4);
    IRREP_ASSERT(cs_framework.H_Z.n_rows == 10);

    /* Build via the direct construction. */
    irrep_css_code_t cs_direct;
    IRREP_ASSERT(irrep_color_3d_rm_15_1_3(&cs_direct) == IRREP_OK);

    /* The two H_X matrices must have identical supports. */
    IRREP_ASSERT(parity_matrices_equal(&cs_framework.H_X, &cs_direct.H_X,
                                       4, 15));
    /* The two H_Z matrices must have identical supports. */
    IRREP_ASSERT(parity_matrices_equal(&cs_framework.H_Z, &cs_direct.H_Z,
                                       10, 15));

    /* CSS orthogonality + commutativity on the framework code. */
    IRREP_ASSERT(irrep_css_code_verify(&cs_framework) == IRREP_OK);
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs_framework, &g)
                 == IRREP_OK);
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);

    /* Brute-force distance d = 3. */
    int d = irrep_qec_distance_brute(&g, /*max_weight*/ 3);
    IRREP_ASSERT(d == 3);

    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs_framework);
    irrep_css_code_free(&cs_direct);
    irrep_color_lattice_3d_tetrahedron_free(&lat);
    return IRREP_TEST_END();
}

static int test_rejects_malformed_input(void) {
    IRREP_TEST_START("color_lattice_3d_rejects_malformed");

    /* NULL inputs. */
    IRREP_ASSERT(irrep_color_lattice_3d_to_css(NULL, NULL)
                 == IRREP_ERR_INVALID_ARG);
    {
        irrep_css_code_t cs;
        IRREP_ASSERT(irrep_color_lattice_3d_to_css(NULL, &cs)
                     == IRREP_ERR_INVALID_ARG);
    }
    /* Bad qubit index. */
    {
        const int bad_face[2] = { 0, 99 };
        const int sizes[1] = { 2 };
        const int *const faces[1] = { bad_face };
        const int dummy_cell[1] = { 0 };
        const int cell_sizes[1] = { 1 };
        const int *const cells[1] = { dummy_cell };
        irrep_color_lattice_3d_t lat = {
            .n_qubits = 4,
            .n_faces = 1, .face_sizes = sizes, .face_qubits = faces,
            .n_cells = 1, .cell_sizes = cell_sizes, .cell_qubits = cells,
        };
        irrep_css_code_t cs;
        IRREP_ASSERT(irrep_color_lattice_3d_to_css(&lat, &cs)
                     == IRREP_ERR_PRECONDITION);
    }
    /* Odd face-cell overlap → CSS orthogonality fails. */
    {
        /* 3 qubits; one face {0, 1}; one cell {1, 2}. Overlap = {1},
         * weight 1 (odd) — should fail CSS. */
        const int face0[2] = { 0, 1 };
        const int cell0[2] = { 1, 2 };
        const int fs[1] = { 2 };
        const int cs_arr[1] = { 2 };
        const int *const faces[1] = { face0 };
        const int *const cells[1] = { cell0 };
        irrep_color_lattice_3d_t lat = {
            .n_qubits = 3,
            .n_faces = 1, .face_sizes = fs, .face_qubits = faces,
            .n_cells = 1, .cell_sizes = cs_arr, .cell_qubits = cells,
        };
        irrep_css_code_t cs;
        IRREP_ASSERT(irrep_color_lattice_3d_to_css(&lat, &cs)
                     == IRREP_ERR_PRECONDITION);
    }
    return IRREP_TEST_END();
}

int main(void) {
    int rc = 0;
    rc |= test_tetrahedron_matches_direct_construction();
    rc |= test_rejects_malformed_input();
    return rc;
}
