/* SPDX-License-Identifier: MIT */
/* Tests for the subsystem-code (gauge-group) framework — Bacon-Shor [[9,1,3]].
 *
 * Verifies:
 *  - 12 gauge generators built (6 X-gauges + 6 Z-gauges).
 *  - The gauge group is NON-abelian: at least one pair of generators
 *    anti-commutes (e.g. X_0 X_1 and Z_0 Z_3 share qubit 0).
 *  - The 4 weight-6 Bacon-Shor stabilizers (S_X1, S_X2, S_Z1, S_Z2) are
 *    all in the centraliser of the gauge group.
 *  - The logical operators L_X = X column 0, L_Z = Z row 0 are also
 *    in the centraliser.
 *  - L_X and L_Z anti-commute on qubit 0 (the defining logical algebra).
 *  - A typical gauge generator (X_0 X_1) is NOT in the centraliser
 *    (it anti-commutes with Z_0 Z_3).
 */
#include "harness.h"
#include <irrep/stabilizer_group.h>
#include <irrep/subsystem_code.h>
#include <stdio.h>

static void set_pauli(irrep_pauli_t *p, const irrep_pauli_letter_t *letters, int n) {
    for (int q = 0; q < n; ++q) {
        if (letters[q] != IRREP_PAULI_LETTER_I) {
            irrep_pauli_set(p, q, letters[q]);
        }
    }
}

static int test_bacon_shor_shape(void) {
    irrep_subsystem_code_t c;
    if (irrep_subsystem_bacon_shor_9_1_3(&c) != IRREP_OK) return 1;
    int rc = (c.n == 9 && c.n_gauge == 12) ? 0 : 1;
    /* Each gauge generator has weight 2. */
    for (int i = 0; i < 12 && rc == 0; ++i) {
        if (irrep_pauli_weight(&c.gauge[i]) != 2) rc = 1;
    }
    irrep_subsystem_code_free(&c);
    return rc;
}

static int test_gauge_non_abelian(void) {
    irrep_subsystem_code_t c;
    if (irrep_subsystem_bacon_shor_9_1_3(&c) != IRREP_OK) return 1;
    int found_anticomm = 0;
    for (int i = 0; i < c.n_gauge; ++i) {
        for (int j = i + 1; j < c.n_gauge; ++j) {
            if (!irrep_pauli_commute(&c.gauge[i], &c.gauge[j])) {
                found_anticomm = 1;
                break;
            }
        }
        if (found_anticomm) break;
    }
    irrep_subsystem_code_free(&c);
    return found_anticomm ? 0 : 1;
}

/* Bacon-Shor stabilizers: all 4 must be in the centraliser. */
static int test_bacon_shor_stabilizers_central(void) {
    irrep_subsystem_code_t c;
    if (irrep_subsystem_bacon_shor_9_1_3(&c) != IRREP_OK) return 1;
    /* S_X1 = X_0 X_1 X_3 X_4 X_6 X_7  (cols 0+1) */
    irrep_pauli_letter_t S_X1[9] = {
        IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_I,
        IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_I,
        IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_I,
    };
    /* S_X2 = X_1 X_2 X_4 X_5 X_7 X_8  (cols 1+2) */
    irrep_pauli_letter_t S_X2[9] = {
        IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_X,
        IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_X,
        IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_X, IRREP_PAULI_LETTER_X,
    };
    /* S_Z1 = Z_0 Z_1 Z_2 Z_3 Z_4 Z_5  (rows 0+1) */
    irrep_pauli_letter_t S_Z1[9] = {
        IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z,
        IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z,
        IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_I,
    };
    /* S_Z2 = Z_3 Z_4 Z_5 Z_6 Z_7 Z_8  (rows 1+2) */
    irrep_pauli_letter_t S_Z2[9] = {
        IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_I, IRREP_PAULI_LETTER_I,
        IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z,
        IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z, IRREP_PAULI_LETTER_Z,
    };
    irrep_pauli_t p;
    int rc = 0;
    irrep_pauli_letter_t *stabs[] = { S_X1, S_X2, S_Z1, S_Z2 };
    for (int k = 0; k < 4; ++k) {
        irrep_pauli_new(&p, 9);
        set_pauli(&p, stabs[k], 9);
        if (!irrep_subsystem_in_centraliser(&c, &p)) rc = 1;
        if (irrep_pauli_weight(&p) != 6) rc = 1;
        irrep_pauli_free(&p);
    }
    irrep_subsystem_code_free(&c);
    return rc;
}

/* Logical operators must be in the centraliser, and L_X · L_Z must
 * anti-commute (defining property of the logical X-Z pair). */
static int test_bacon_shor_logical(void) {
    irrep_subsystem_code_t c;
    if (irrep_subsystem_bacon_shor_9_1_3(&c) != IRREP_OK) return 1;
    /* L_X = X_0 X_3 X_6 (column 0). */
    irrep_pauli_t Lx, Lz;
    irrep_pauli_new(&Lx, 9);
    irrep_pauli_set(&Lx, 0, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&Lx, 3, IRREP_PAULI_LETTER_X);
    irrep_pauli_set(&Lx, 6, IRREP_PAULI_LETTER_X);
    /* L_Z = Z_0 Z_1 Z_2 (row 0). */
    irrep_pauli_new(&Lz, 9);
    irrep_pauli_set(&Lz, 0, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&Lz, 1, IRREP_PAULI_LETTER_Z);
    irrep_pauli_set(&Lz, 2, IRREP_PAULI_LETTER_Z);

    int rc = 0;
    if (!irrep_subsystem_in_centraliser(&c, &Lx)) rc = 1;
    if (!irrep_subsystem_in_centraliser(&c, &Lz)) rc = 1;
    /* L_X and L_Z share exactly qubit 0 → anti-commute. */
    if (irrep_pauli_commute(&Lx, &Lz)) rc = 1;
    irrep_pauli_free(&Lx);
    irrep_pauli_free(&Lz);
    irrep_subsystem_code_free(&c);
    return rc;
}

/* A typical gauge generator (X_0 X_1) must NOT be in the centraliser. */
static int test_gauge_not_central(void) {
    irrep_subsystem_code_t c;
    if (irrep_subsystem_bacon_shor_9_1_3(&c) != IRREP_OK) return 1;
    int rc = irrep_subsystem_in_centraliser(&c, &c.gauge[0]) ? 1 : 0;
    irrep_subsystem_code_free(&c);
    return rc;
}

/* Tetrahedral 3D Bombín gauge color code substrate.
 *
 *   - shape: 4 qubits, 12 gauge generators (6 XX-edges + 6 ZZ-edges).
 *   - non-abelian: XX-edge and ZZ-edge sharing one vertex anti-commute.
 *   - center contains X⊗4 and Z⊗4 (commute with every weight-2 edge).
 *   - center does NOT contain a single ZZ-edge (anti-commutes with the
 *     5 other ZZ-edges that share one vertex with it).
 */
static int test_bombin_3d_tetrahedron_shape(void) {
    irrep_subsystem_code_t c;
    if (irrep_subsystem_bombin_3d_tetrahedron(&c) != IRREP_OK) return 1;
    int rc = (c.n == 4 && c.n_gauge == 12) ? 0 : 1;
    /* Each gauge generator has weight exactly 2. */
    for (int i = 0; i < c.n_gauge; ++i) {
        if (irrep_pauli_weight(&c.gauge[i]) != 2) rc = 1;
    }
    irrep_subsystem_code_free(&c);
    return rc;
}

static int test_bombin_3d_tetrahedron_non_abelian(void) {
    irrep_subsystem_code_t c;
    if (irrep_subsystem_bombin_3d_tetrahedron(&c) != IRREP_OK) return 1;
    /* gauge[0] = X_0 X_1, gauge[6] = Z_0 Z_1. Sharing both qubits,
     * XX and ZZ commute (two Pauli sign flips cancel). */
    int rc = irrep_pauli_commute(&c.gauge[0], &c.gauge[6]) ? 0 : 1;
    /* gauge[0] = X_0 X_1, gauge[7] = Z_0 Z_2. Sharing one qubit (q=0),
     * XX and ZZ anti-commute. */
    if (irrep_pauli_commute(&c.gauge[0], &c.gauge[7])) rc = 1;
    irrep_subsystem_code_free(&c);
    return rc;
}

static int test_bombin_3d_tetrahedron_center(void) {
    irrep_subsystem_code_t c;
    if (irrep_subsystem_bombin_3d_tetrahedron(&c) != IRREP_OK) return 1;
    int rc = 0;
    /* X⊗4 commutes with every edge generator → in centraliser. */
    irrep_pauli_t Xall;
    irrep_pauli_new(&Xall, 4);
    for (int q = 0; q < 4; ++q) irrep_pauli_set(&Xall, q, IRREP_PAULI_LETTER_X);
    if (!irrep_subsystem_in_centraliser(&c, &Xall)) rc = 1;
    irrep_pauli_free(&Xall);
    /* Z⊗4 likewise. */
    irrep_pauli_t Zall;
    irrep_pauli_new(&Zall, 4);
    for (int q = 0; q < 4; ++q) irrep_pauli_set(&Zall, q, IRREP_PAULI_LETTER_Z);
    if (!irrep_subsystem_in_centraliser(&c, &Zall)) rc = 1;
    irrep_pauli_free(&Zall);
    /* A single ZZ-edge (gauge[6] = Z_0 Z_1) anti-commutes with the
     * other ZZ-edges sharing one vertex (e.g. gauge[7] = Z_0 Z_2 — wait,
     * ZZ and ZZ always commute). Use XX-edge as the test: gauge[1] =
     * X_0 X_2 anti-commutes with Z_0 Z_1 (shared q=0, lone Pauli flip),
     * so Z_0 Z_1 NOT in centraliser. */
    if (irrep_subsystem_in_centraliser(&c, &c.gauge[6])) rc = 1;
    irrep_subsystem_code_free(&c);
    return rc;
}

int main(void) {
    int rc = 0;
    if (test_bacon_shor_shape())             { fprintf(stderr, "FAIL test_bacon_shor_shape\n"); rc = 1; }
    if (test_gauge_non_abelian())            { fprintf(stderr, "FAIL test_gauge_non_abelian\n"); rc = 1; }
    if (test_bacon_shor_stabilizers_central()) { fprintf(stderr, "FAIL test_bacon_shor_stabilizers_central\n"); rc = 1; }
    if (test_bacon_shor_logical())           { fprintf(stderr, "FAIL test_bacon_shor_logical\n"); rc = 1; }
    if (test_gauge_not_central())            { fprintf(stderr, "FAIL test_gauge_not_central\n"); rc = 1; }
    if (test_bombin_3d_tetrahedron_shape())       { fprintf(stderr, "FAIL test_bombin_3d_tetrahedron_shape\n"); rc = 1; }
    if (test_bombin_3d_tetrahedron_non_abelian()) { fprintf(stderr, "FAIL test_bombin_3d_tetrahedron_non_abelian\n"); rc = 1; }
    if (test_bombin_3d_tetrahedron_center())      { fprintf(stderr, "FAIL test_bombin_3d_tetrahedron_center\n"); rc = 1; }
    return rc;
}
