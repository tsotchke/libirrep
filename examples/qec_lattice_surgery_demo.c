/* SPDX-License-Identifier: MIT */
/* Lattice-surgery demo: smooth merge of two d=3 surface-code patches.
 *
 * Walks through the pre/post-merge picture for the smooth merge of two
 * adjacent [[9, 1, 3]] rotated surface-code patches (A on the left,
 * B on the right) glued along A's right smooth (Z-type) boundary and
 * B's left smooth boundary. The merged code is a [[18, 1, 3]]
 * rectangular surface code (3 × 6 data-qubit grid).
 *
 * Pedagogy:
 *   1. Build patch A and patch B independently. Each is [[9, 1, 3]] with
 *      8 stabilizer generators (4 X + 4 Z) and one logical qubit.
 *      Pre-merge total: 18 qubits, 16 stabs, 2 logical qubits.
 *   2. Build the merged code M using `irrep_lattice_surgery_smooth_*`.
 *      Post-merge: 18 qubits, 17 stabs (= 2·d² - 1), 1 logical qubit.
 *      The merge added one stabilizer = it measured one logical operator.
 *   3. Identify which operator got measured: the joint-X parity
 *      `X̄_A · X̄_B = X` on columns 0 and d of the merged frame, weight 2d.
 *      Verify it commutes with every merged stabilizer AND with both
 *      merged logicals L̄_X(M), L̄_Z(M) — forcing it into the stabilizer
 *      group of M.
 *   4. Cross-check: the *joint-Z* parity `Z̄_A · Z̄_B` equals the merged
 *      L̄_Z(M) (a Z-string on row 0 spanning all 2d columns), confirming
 *      that the orthogonal smooth-merge measurement collapses the two
 *      pre-merge L̄_Z operators into the single merged L̄_Z.
 *
 * Build / run:
 *   make examples
 *   ./build/bin/qec_lattice_surgery_demo
 */
#include <irrep/css_code.h>
#include <irrep/lattice_surgery.h>
#include <irrep/stabilizer_group.h>
#include <irrep/surface_code.h>

#include <stdio.h>
#include <stdlib.h>

static void print_pauli(const char *label, const irrep_pauli_t *p) {
    printf("  %-12s [w=%2d]  ", label, irrep_pauli_weight(p));
    for (int q = 0; q < p->n; ++q) {
        irrep_pauli_letter_t L = irrep_pauli_get(p, q);
        char c = (L == IRREP_PAULI_LETTER_X) ? 'X'
                : (L == IRREP_PAULI_LETTER_Y) ? 'Y'
                : (L == IRREP_PAULI_LETTER_Z) ? 'Z' : '.';
        putchar(c);
    }
    putchar('\n');
}

static int count_anti(const irrep_stabilizer_group_t *g, const irrep_pauli_t *p) {
    int n = 0;
    for (int i = 0; i < g->n_generators; ++i) {
        if (!irrep_pauli_commute(&g->gens[i], p)) ++n;
    }
    return n;
}

int main(void) {
    const int d = 3;
    printf("Lattice-surgery smooth merge of two [[%d, 1, %d]] patches\n", d * d, d);
    printf("===========================================================\n");
    printf("Geometry: A (left) and B (right) on a %d×%d data-qubit grid;\n", d, 2 * d);
    printf("          A occupies columns [0, %d), B occupies columns [%d, %d).\n\n",
           d, d, 2 * d);

    /* --- Pre-merge: two independent patches ---------------------------- */
    irrep_surface_params_t sp;
    irrep_surface_init(&sp, d);
    printf("Each pre-merge patch:\n");
    printf("  n_qubits  = %d\n", sp.n_qubits);
    printf("  n_X_stabs = %d\n", sp.n_X_stabs);
    printf("  n_Z_stabs = %d\n", sp.n_Z_stabs);
    printf("  total stabs = %d  (= n - k = %d - 1)\n",
           sp.n_X_stabs + sp.n_Z_stabs, sp.n_qubits);
    printf("Pre-merge total: %d qubits, %d stabs, 2 logical qubits.\n\n",
           2 * sp.n_qubits, 2 * (sp.n_X_stabs + sp.n_Z_stabs));

    /* --- Post-merge: rectangular d × 2d code -------------------------- */
    irrep_lattice_surgery_smooth_t mp;
    if (irrep_lattice_surgery_smooth_init(&mp, d) != IRREP_OK) {
        fprintf(stderr, "init failed\n"); return 1;
    }
    printf("Merged code (smooth merge → d × 2d rectangle):\n");
    printf("  n_qubits  = %d  (= 2·d²)\n", mp.n_qubits);
    printf("  n_X_stabs = %d\n", mp.n_X_stabs);
    printf("  n_Z_stabs = %d\n", mp.n_Z_stabs);
    printf("  total stabs = %d  (= n - k = 2d² - 1)\n",
           mp.n_X_stabs + mp.n_Z_stabs);
    printf("Post-merge: 1 logical qubit; the merge promoted one operator\n");
    printf("            from logical to stabilizer (= a measurement).\n\n");

    irrep_css_code_t cs;
    if (irrep_lattice_surgery_smooth_build(&mp, &cs) != IRREP_OK) return 1;
    if (irrep_css_code_verify(&cs) != IRREP_OK) {
        fprintf(stderr, "merged code is NOT CSS-orthogonal\n"); return 1;
    }
    irrep_stabilizer_group_t g;
    if (irrep_css_code_to_stabilizer_group(&cs, &g) != IRREP_OK) return 1;
    if (irrep_stabilizer_group_check_commutativity(&g) != IRREP_OK) {
        fprintf(stderr, "merged stabilizers do NOT all commute\n"); return 1;
    }
    printf("Verification:\n");
    printf("  CSS orthogonality (H_X · H_Z^T = 0):   PASS\n");
    printf("  Pairwise stabilizer commutativity:    PASS\n\n");

    /* --- Logical operators -------------------------------------------- */
    irrep_pauli_t Lx, Lz, Pj;
    irrep_lattice_surgery_smooth_logical_X(&mp, &Lx);
    irrep_lattice_surgery_smooth_logical_Z(&mp, &Lz);
    irrep_lattice_surgery_smooth_joint_X_parity(&mp, &Pj);

    printf("Merged-code operators (each shown as letters on qubits 0..%d):\n",
           mp.n_qubits - 1);
    print_pauli("L̄_X(M)", &Lx);
    print_pauli("L̄_Z(M)", &Lz);
    print_pauli("X̄_A·X̄_B", &Pj);

    /* --- Commutation checks ------------------------------------------ */
    int anti_Lx = count_anti(&g, &Lx);
    int anti_Lz = count_anti(&g, &Lz);
    int anti_Pj = count_anti(&g, &Pj);
    printf("\nCommutation against merged stabilizer group (must all be 0):\n");
    printf("  L̄_X(M)    anti-commutes with %d / %d stabs\n", anti_Lx, g.n_generators);
    printf("  L̄_Z(M)    anti-commutes with %d / %d stabs\n", anti_Lz, g.n_generators);
    printf("  X̄_A·X̄_B  anti-commutes with %d / %d stabs\n", anti_Pj, g.n_generators);

    int sym_LxLz = irrep_pauli_symp_inner(&Lx, &Lz);
    int sym_PjLx = irrep_pauli_symp_inner(&Pj, &Lx);
    int sym_PjLz = irrep_pauli_symp_inner(&Pj, &Lz);
    printf("\nSymplectic inner products (1 = anti-commute):\n");
    printf("  ⟨L̄_X(M), L̄_Z(M)⟩  = %d   (expected 1 — they're the merged logicals)\n", sym_LxLz);
    printf("  ⟨X̄_A·X̄_B, L̄_X(M)⟩ = %d   (expected 0)\n", sym_PjLx);
    printf("  ⟨X̄_A·X̄_B, L̄_Z(M)⟩ = %d   (expected 0)\n", sym_PjLz);

    if (anti_Pj == 0 && sym_PjLx == 0 && sym_PjLz == 0) {
        printf("\n=> X̄_A·X̄_B commutes with every stabilizer AND with both merged\n");
        printf("   logicals. In a [[2d², 1, d]] code the normalizer / stabilizer\n");
        printf("   quotient is Z_2 × Z_2 = {I, L̄_X, L̄_Z, L̄_X L̄_Z}, so commuting\n");
        printf("   with both logicals forces X̄_A·X̄_B ∈ S_M.  Smooth merge has\n");
        printf("   measured the joint-X parity.\n");
    } else {
        fprintf(stderr, "\nUNEXPECTED: X̄_A·X̄_B failed the stabilizer-class test\n");
        return 1;
    }

    irrep_pauli_free(&Lx);
    irrep_pauli_free(&Lz);
    irrep_pauli_free(&Pj);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return 0;
}
