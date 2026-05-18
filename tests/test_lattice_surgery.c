/* SPDX-License-Identifier: MIT */
/* Tests for smooth-merge lattice surgery on rotated surface codes.
 *
 * Verifies that the merged d × 2d rectangular code is a well-formed CSS
 * code (orthogonality + stabilizer commutativity), that its logical X̄
 * (column 0, weight d) and Z̄ (row 0, weight 2d) commute with every
 * stabilizer and anti-commute pairwise, and — most importantly — that
 * the joint-X parity operator X̄_A · X̄_B (= X on columns 0 and d) is
 * in the merged stabilizer group:
 *
 *   P_joint commutes with every merged stabilizer AND with L̄_X(M)
 *   AND with L̄_Z(M)  →  in a [[2d², 1, d]] code, this forces
 *   P_joint ∈ S_M (a stabilizer, not a logical).
 *
 * Stabilizer counts at d ∈ {2, 3, 5, 7} are also pinned against the
 * closed-form prediction n - k = 2d² - 1.
 */
#include "harness.h"
#include <irrep/css_code.h>
#include <irrep/lattice_surgery.h>
#include <irrep/stabilizer_group.h>

#include <stdio.h>

static int run_for_d(int d) {
    char banner[64];
    snprintf(banner, sizeof banner, "lattice_surgery_smooth_d%d", d);
    IRREP_TEST_START(banner);

    irrep_lattice_surgery_smooth_t p;
    IRREP_ASSERT(irrep_lattice_surgery_smooth_init(&p, d) == IRREP_OK);
    IRREP_ASSERT(p.d == d);
    IRREP_ASSERT(p.rows == d);
    IRREP_ASSERT(p.cols == 2 * d);
    IRREP_ASSERT(p.n_qubits == 2 * d * d);

    /* n - k = 2d² - 1 since merged code is [[2d², 1, d]]. */
    IRREP_ASSERT(p.n_X_stabs + p.n_Z_stabs == 2 * d * d - 1);

    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_lattice_surgery_smooth_build(&p, &cs) == IRREP_OK);
    IRREP_ASSERT(cs.n == 2 * d * d);
    IRREP_ASSERT(cs.H_X.n_rows == p.n_X_stabs);
    IRREP_ASSERT(cs.H_Z.n_rows == p.n_Z_stabs);

    /* CSS orthogonality: H_X · H_Z^T = 0 (mod 2). */
    IRREP_ASSERT(irrep_css_code_verify(&cs) == IRREP_OK);

    /* Full stabilizer-group commutativity (subsumes CSS but also catches
     * any internal anti-commutation we might have introduced). */
    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);
    IRREP_ASSERT(g.n_generators == 2 * d * d - 1);

    /* Logical operators. */
    irrep_pauli_t Lx, Lz, Pj;
    IRREP_ASSERT(irrep_lattice_surgery_smooth_logical_X(&p, &Lx) == IRREP_OK);
    IRREP_ASSERT(irrep_lattice_surgery_smooth_logical_Z(&p, &Lz) == IRREP_OK);
    IRREP_ASSERT(irrep_lattice_surgery_smooth_joint_X_parity(&p, &Pj) == IRREP_OK);

    /* L̄_X = X on column 0 → weight d.  L̄_Z = Z on row 0 spanning 2d
     * columns → weight 2d.  P_joint = X on columns 0 and d → weight 2d. */
    IRREP_ASSERT(irrep_pauli_weight(&Lx) == d);
    IRREP_ASSERT(irrep_pauli_weight(&Lz) == 2 * d);
    IRREP_ASSERT(irrep_pauli_weight(&Pj) == 2 * d);

    /* L̄_X and L̄_Z share exactly one qubit q(0,0) where one is X and the
     * other Z — they anti-commute there. */
    IRREP_ASSERT(irrep_pauli_symp_inner(&Lx, &Lz) == 1);

    /* Logicals commute with every stabilizer. */
    int anti_Lx = 0, anti_Lz = 0, anti_Pj = 0;
    int anti_Pj_Lx = irrep_pauli_symp_inner(&Pj, &Lx);
    int anti_Pj_Lz = irrep_pauli_symp_inner(&Pj, &Lz);
    for (int i = 0; i < g.n_generators; ++i) {
        if (!irrep_pauli_commute(&g.gens[i], &Lx)) ++anti_Lx;
        if (!irrep_pauli_commute(&g.gens[i], &Lz)) ++anti_Lz;
        if (!irrep_pauli_commute(&g.gens[i], &Pj)) ++anti_Pj;
    }
    IRREP_ASSERT(anti_Lx == 0);
    IRREP_ASSERT(anti_Lz == 0);

    /* P_joint commutes with every stabilizer of the merged code. */
    IRREP_ASSERT(anti_Pj == 0);

    /* P_joint commutes with BOTH merged logicals → forces P_joint ∈ S_M.
     * (In a [[n, 1, d]] code, the normalizer / stabilizer quotient is
     * Z_2 × Z_2 = {I, L̄_X, L̄_Z, L̄_X L̄_Z}; commuting with both L̄_X
     * and L̄_Z is the I coset, i.e. a stabilizer.) */
    IRREP_ASSERT(anti_Pj_Lx == 0);
    IRREP_ASSERT(anti_Pj_Lz == 0);

    irrep_pauli_free(&Lx);
    irrep_pauli_free(&Lz);
    irrep_pauli_free(&Pj);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);

    return IRREP_TEST_END();
}

static int test_invalid_args(void) {
    IRREP_TEST_START("lattice_surgery_invalid_args");
    irrep_lattice_surgery_smooth_t p;
    IRREP_ASSERT(irrep_lattice_surgery_smooth_init(NULL, 3) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_lattice_surgery_smooth_init(&p, 1)   == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_lattice_surgery_smooth_init(&p, 0)   == IRREP_ERR_INVALID_ARG);

    IRREP_ASSERT(irrep_lattice_surgery_smooth_init(&p, 3) == IRREP_OK);
    IRREP_ASSERT(irrep_lattice_surgery_smooth_build(NULL, NULL) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_lattice_surgery_smooth_logical_X(NULL, NULL) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_lattice_surgery_smooth_logical_Z(NULL, NULL) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_lattice_surgery_smooth_joint_X_parity(NULL, NULL) == IRREP_ERR_INVALID_ARG);
    return IRREP_TEST_END();
}

/* ====================================================================
 * Rough-merge tests: 2d × d rectangular code, joint-Z parity stabilizer.
 * ==================================================================== */

static int run_rough_for_d(int d) {
    char banner[64];
    snprintf(banner, sizeof banner, "lattice_surgery_rough_d%d", d);
    IRREP_TEST_START(banner);

    irrep_lattice_surgery_rough_t p;
    IRREP_ASSERT(irrep_lattice_surgery_rough_init(&p, d) == IRREP_OK);
    IRREP_ASSERT(p.d == d);
    IRREP_ASSERT(p.rows == 2 * d);
    IRREP_ASSERT(p.cols == d);
    IRREP_ASSERT(p.n_qubits == 2 * d * d);
    IRREP_ASSERT(p.n_X_stabs + p.n_Z_stabs == 2 * d * d - 1);

    irrep_css_code_t cs;
    IRREP_ASSERT(irrep_lattice_surgery_rough_build(&p, &cs) == IRREP_OK);
    IRREP_ASSERT(cs.n == 2 * d * d);
    IRREP_ASSERT(cs.H_X.n_rows == p.n_X_stabs);
    IRREP_ASSERT(cs.H_Z.n_rows == p.n_Z_stabs);
    IRREP_ASSERT(irrep_css_code_verify(&cs) == IRREP_OK);

    irrep_stabilizer_group_t g;
    IRREP_ASSERT(irrep_css_code_to_stabilizer_group(&cs, &g) == IRREP_OK);
    IRREP_ASSERT(irrep_stabilizer_group_check_commutativity(&g) == IRREP_OK);
    IRREP_ASSERT(g.n_generators == 2 * d * d - 1);

    irrep_pauli_t Lx, Lz, Pj;
    IRREP_ASSERT(irrep_lattice_surgery_rough_logical_X(&p, &Lx) == IRREP_OK);
    IRREP_ASSERT(irrep_lattice_surgery_rough_logical_Z(&p, &Lz) == IRREP_OK);
    IRREP_ASSERT(irrep_lattice_surgery_rough_joint_Z_parity(&p, &Pj) == IRREP_OK);

    /* Rough merge: L̄_X spans 2d rows (column 0), L̄_Z spans d cols (row 0);
     * joint-Z parity = Z on row 0 AND row d = weight 2d. */
    IRREP_ASSERT(irrep_pauli_weight(&Lx) == 2 * d);
    IRREP_ASSERT(irrep_pauli_weight(&Lz) == d);
    IRREP_ASSERT(irrep_pauli_weight(&Pj) == 2 * d);

    /* L̄_X and L̄_Z anti-commute on corner q(0, 0). */
    IRREP_ASSERT(irrep_pauli_symp_inner(&Lx, &Lz) == 1);

    int anti_Lx = 0, anti_Lz = 0, anti_Pj = 0;
    for (int i = 0; i < g.n_generators; ++i) {
        if (!irrep_pauli_commute(&g.gens[i], &Lx)) ++anti_Lx;
        if (!irrep_pauli_commute(&g.gens[i], &Lz)) ++anti_Lz;
        if (!irrep_pauli_commute(&g.gens[i], &Pj)) ++anti_Pj;
    }
    IRREP_ASSERT(anti_Lx == 0);
    IRREP_ASSERT(anti_Lz == 0);
    IRREP_ASSERT(anti_Pj == 0);

    /* P_joint commutes with BOTH merged logicals → stabilizer of M. */
    IRREP_ASSERT(irrep_pauli_symp_inner(&Pj, &Lx) == 0);
    IRREP_ASSERT(irrep_pauli_symp_inner(&Pj, &Lz) == 0);

    irrep_pauli_free(&Lx);
    irrep_pauli_free(&Lz);
    irrep_pauli_free(&Pj);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);

    return IRREP_TEST_END();
}

int main(void) {
    int rc = 0;
    rc |= run_for_d(2);
    rc |= run_for_d(3);
    rc |= run_for_d(5);
    rc |= run_for_d(7);
    rc |= run_rough_for_d(2);
    rc |= run_rough_for_d(3);
    rc |= run_rough_for_d(5);
    rc |= run_rough_for_d(7);
    rc |= test_invalid_args();
    return rc;
}
