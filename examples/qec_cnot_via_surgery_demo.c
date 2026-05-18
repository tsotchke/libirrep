/* SPDX-License-Identifier: MIT */
/* CNOT-via-lattice-surgery: the canonical fault-tolerant logical-CNOT
 * protocol on surface-code patches (Horsman-Fowler-Devitt-Van Meter 2012).
 *
 * # Protocol
 *
 * Three rotated surface-code patches:
 *   - C (control) at distance d
 *   - T (target)  at distance d
 *   - A (ancilla) at distance d, initialised in |+⟩_L (the +1 eigenstate
 *     of X̄_A — produced by measuring all Z-stabilizers and applying
 *     correction to fix any -1 eigenvalues; this initialisation step
 *     is standard in any surface-code implementation and we treat it
 *     here as a precondition).
 *
 * Steps:
 *   1. SMOOTH-MERGE C with A (along their facing smooth boundaries).
 *      This measures the joint-X parity X̄_C · X̄_A; outcome ∈ {+1, -1}
 *      becomes the syndrome bit m_X.
 *   2. SMOOTH-SPLIT the merged code back into C and A patches (turn
 *      off the seam stabilizers).
 *   3. ROUGH-MERGE A with T (along their facing rough boundaries).
 *      This measures the joint-Z parity Z̄_A · Z̄_T; outcome = m_Z.
 *   4. ROUGH-SPLIT back into A and T patches.
 *   5. Measure A in the Z basis (= measure Z̄_A) → outcome m_A.
 *   6. Apply CONDITIONAL CORRECTIONS:
 *         X̄_T  if  m_X = -1  (= +1 in Z₂-parity convention)
 *         Z̄_C  if  m_Z · m_A = -1
 *      (signs are convention-dependent; the structural fact is that
 *      the gate is achieved up to a known Pauli-frame correction).
 *   7. The encoded state is now |CNOT⟩(|ψ⟩_C ⊗ |φ⟩_T).
 *
 * # Pedagogical scope of this demo
 *
 * This demo prints the structural data for steps 1–6 at d = 3:
 *   - The 3 patches' geometries (each [[9, 1, 3]]).
 *   - The smooth-merge stabilizer count for C ⊕ A → 2d × d rectangle.
 *   - The rough-merge stabilizer count for A ⊕ T → d × 2d rectangle.
 *   - The joint-parity operators that get measured at each merge.
 *   - The correction table relating measurement outcomes (m_X, m_Z, m_A)
 *     to Pauli-frame corrections on the C/T logical qubits.
 * This makes the gadget concrete without needing a full simulator;
 * the actual stabilizer measurements + corrections are deferred to
 * a quantum-circuit simulator (which would consume libirrep's
 * stabilizer-group output).
 *
 * # Primary references
 *
 *   - Horsman-Fowler-Devitt-Van Meter, *Surface code quantum computing
 *     by lattice surgery*, New J. Phys. 14 (2012) 123011
 *     [arXiv:1111.4022].
 *   - Litinski, *A game of surface codes*, Quantum 3 (2019) 128
 *     [arXiv:1808.02892].
 *
 * Build / run:
 *   make examples
 *   ./build/bin/qec_cnot_via_surgery_demo
 */
#include <irrep/css_code.h>
#include <irrep/lattice_surgery.h>
#include <irrep/stabilizer_group.h>
#include <irrep/surface_code.h>

#include <stdio.h>
#include <stdlib.h>

static void print_pauli(const char *label, const irrep_pauli_t *p) {
    printf("    %-12s [w=%2d]  ", label, irrep_pauli_weight(p));
    for (int q = 0; q < p->n; ++q) {
        irrep_pauli_letter_t L = irrep_pauli_get(p, q);
        char c = (L == IRREP_PAULI_LETTER_X) ? 'X'
                : (L == IRREP_PAULI_LETTER_Y) ? 'Y'
                : (L == IRREP_PAULI_LETTER_Z) ? 'Z' : '.';
        putchar(c);
    }
    putchar('\n');
}

int main(void) {
    const int d = 3;
    printf("CNOT-via-lattice-surgery on three [[%d, 1, %d]] patches\n", d * d, d);
    printf("====================================================================\n\n");

    /* Per-patch geometry. */
    irrep_surface_params_t patch;
    irrep_surface_init(&patch, d);
    printf("Each patch: n = %d data qubits, %d X-stabs + %d Z-stabs, k = 1.\n",
           patch.n_qubits, patch.n_X_stabs, patch.n_Z_stabs);
    printf("Three patches total: %d data qubits across C, A, T.\n\n",
           3 * patch.n_qubits);

    /* --- Step 1: Smooth-merge C with A ---------------------------- */
    printf("Step 1: Smooth-merge C with A → measures X̄_C · X̄_A.\n");
    irrep_lattice_surgery_smooth_t smooth;
    irrep_lattice_surgery_smooth_init(&smooth, d);
    printf("  Merged geometry: %d × %d rectangular surface code.\n",
           smooth.rows, smooth.cols);
    printf("  Merged-code stabilizers: %d X + %d Z = %d total.\n",
           smooth.n_X_stabs, smooth.n_Z_stabs,
           smooth.n_X_stabs + smooth.n_Z_stabs);

    irrep_pauli_t joint_X;
    irrep_lattice_surgery_smooth_joint_X_parity(&smooth, &joint_X);
    printf("  Joint-X parity operator (the measurement outcome m_X):\n");
    print_pauli("X̄_C · X̄_A", &joint_X);
    irrep_pauli_free(&joint_X);
    printf("\n");

    /* --- Step 2: Smooth-split (= reverse of merge) ---------------- */
    printf("Step 2: Smooth-split → restores separate C and A patches.\n");
    printf("  Operationally: turn off the seam-X stabilizers; the merged\n");
    printf("  code's L̄_Z (row-0 Z-string of length 2d) splits into two\n");
    printf("  half-strings — one per patch — each = the patch's L̄_Z.\n\n");

    /* --- Step 3: Rough-merge A with T ----------------------------- */
    printf("Step 3: Rough-merge A with T → measures Z̄_A · Z̄_T.\n");
    irrep_lattice_surgery_rough_t rough;
    irrep_lattice_surgery_rough_init(&rough, d);
    printf("  Merged geometry: %d × %d rectangular surface code.\n",
           rough.rows, rough.cols);
    printf("  Merged-code stabilizers: %d X + %d Z = %d total.\n",
           rough.n_X_stabs, rough.n_Z_stabs,
           rough.n_X_stabs + rough.n_Z_stabs);

    irrep_pauli_t joint_Z;
    irrep_lattice_surgery_rough_joint_Z_parity(&rough, &joint_Z);
    printf("  Joint-Z parity operator (the measurement outcome m_Z):\n");
    print_pauli("Z̄_A · Z̄_T", &joint_Z);
    irrep_pauli_free(&joint_Z);
    printf("\n");

    /* --- Step 4: Rough-split -------------------------------------- */
    printf("Step 4: Rough-split → restores separate A and T patches.\n\n");

    /* --- Steps 5, 6: ancilla measurement + corrections ------------ */
    printf("Step 5: Measure A in the Z basis → outcome m_A.\n");
    printf("\n");
    printf("Step 6: Conditional Pauli-frame corrections\n");
    printf("        based on the three measurement bits (m_X, m_Z, m_A):\n");
    printf("\n");
    printf("        +-------+-------+-------+--------------------------+\n");
    printf("        | m_X   | m_Z   | m_A   | Correction               |\n");
    printf("        +-------+-------+-------+--------------------------+\n");
    printf("        | +1    | +1    | +1    | (none)                   |\n");
    printf("        | -1    | +1    | +1    | X̄_T                      |\n");
    printf("        | +1    | -1    | +1    | (none — Z̄_C cancels)     |\n");
    printf("        | +1    | +1    | -1    | Z̄_C                      |\n");
    printf("        | -1    | -1    | +1    | X̄_T                      |\n");
    printf("        | -1    | +1    | -1    | X̄_T · Z̄_C                |\n");
    printf("        | +1    | -1    | -1    | (none)                   |\n");
    printf("        | -1    | -1    | -1    | X̄_T · Z̄_C                |\n");
    printf("        +-------+-------+-------+--------------------------+\n");
    printf("\n");
    printf("        (Conventions match Horsman et al. 2012 Eq. (10).)\n");
    printf("\n");

    /* --- Net result ---------------------------------------------- */
    printf("Result: The protocol implements the logical CNOT(C → T) up to\n");
    printf("a Pauli-frame correction determined by the three measurement\n");
    printf("outcomes. Total physical-qubit cost: %d data + ~%d ancilla\n",
           3 * patch.n_qubits, 3 * (patch.n_X_stabs + patch.n_Z_stabs));
    printf("(measurement ancillas, ~1 per stabilizer per round). Wall-\n");
    printf("clock cost: O(d) rounds of stabilizer measurement for the\n");
    printf("merge (= the seam needs to be measured d times for fault\n");
    printf("tolerance), so the gadget has logical-time cost d.\n");
    return 0;
}
