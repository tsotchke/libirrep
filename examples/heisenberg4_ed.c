/* SPDX-License-Identifier: MIT */
/* Exact diagonalisation of the 4-site periodic Heisenberg antiferromagnet on
 * a 2×2 square cluster, end-to-end through three 1.3 modules:
 *
 *   lattice.h       → the 4-site lattice and its 4 NN bonds (under PBC)
 *   rdm.h           → Hermitian eigendecomposition and partial-trace
 *                     entanglement entropy
 *   spin_project.h  → total-J = 0 singlet projector verifying the ground
 *                     state carries zero net spin
 *
 * The Hamiltonian is
 *
 *      H = J Σ_{⟨i,j⟩}  S_i · S_j ,       J = 1, S = ½ on each site,
 *
 * with S_i · S_j = ¼ σ_i · σ_j on spin-½. For the 2×2 torus the NN-bond
 * set wraps to 4 unique bonds (each site has 2 distinct neighbours under
 * PBC on a 2×2 square). Analytical ground state:
 *
 *      E_0  = −2 J       (verify)
 *      S²   = 0          (singlet, verify via spin_project)
 *      S_VN ≈ 0.837 nats at a 2-vs-2 bipartition on this torus cluster
 *                        — the reduced density matrix is mixed across the
 *                        two-site sub-Hilbert space because the ground
 *                        state is a non-factorising superposition of the
 *                        two distinct singlets compatible with the torus
 *
 * Compile with `make examples`.
 */

#include <irrep/lattice.h>
#include <irrep/rdm.h>
#include <irrep/spin_project.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N           4
#define DIM         (1 << N)    /* 16 */

/* Single-qubit Pauli-Z, Pauli-X, Pauli-Y as action on |ψ⟩ stored as length-DIM
 * complex vector. bit 0 = site 0 (lowest), bit N-1 = site N-1 (highest). */

static void apply_ZZ(int i, int j, double J, const double _Complex *psi, double _Complex *out) {
    for (int k = 0; k < DIM; ++k) {
        int zi = (k >> i) & 1;
        int zj = (k >> j) & 1;
        double s = (zi ^ zj) ? -1.0 : +1.0;
        out[k] += (J * 0.25 * s) * psi[k];        /* (¼) S_z^i S_z^j factor */
    }
}

static void apply_XX_plus_YY(int i, int j, double J,
                             const double _Complex *psi, double _Complex *out) {
    /* σ_x^i σ_x^j + σ_y^i σ_y^j = 2 (σ_+^i σ_-^j + σ_-^i σ_+^j).
     * Contribution to H_ij: (¼) * (σ_x σ_x + σ_y σ_y) / 2 = (¼) · (σ_+ σ_- + σ_- σ_+).
     * Wait: S_i · S_j = ¼ σ_i · σ_j = ¼ (σ_x^i σ_x^j + σ_y^i σ_y^j + σ_z^i σ_z^j).
     * σ_x σ_x + σ_y σ_y = 2 (|01⟩⟨10| + |10⟩⟨01|) on the (i,j) pair.
     * So the contribution is (¼) · 2 · (pair-flip). */
    for (int k = 0; k < DIM; ++k) {
        int bi = (k >> i) & 1, bj = (k >> j) & 1;
        if (bi == bj) continue;
        int k2 = k ^ (1 << i) ^ (1 << j);
        out[k2] += (J * 0.5) * psi[k];
    }
}

static void build_H(const int *bi, const int *bj, int nb, double J,
                    double _Complex *H) {
    memset(H, 0, sizeof(double _Complex) * DIM * DIM);
    double _Complex e[DIM], h[DIM];
    for (int col = 0; col < DIM; ++col) {
        for (int k = 0; k < DIM; ++k) { e[k] = 0.0; h[k] = 0.0; }
        e[col] = 1.0;
        for (int b = 0; b < nb; ++b) {
            apply_ZZ       (bi[b], bj[b], J, e, h);
            apply_XX_plus_YY(bi[b], bj[b], J, e, h);
        }
        for (int row = 0; row < DIM; ++row) H[row * DIM + col] = h[row];
    }
}

int main(void) {
    /* 1. Lattice: 2×2 periodic square. Emits 4 NN bonds (2 unique under PBC). */
    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 2, 2);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    printf("2×2 square cluster: %d NN bonds (under PBC)\n", nb);

    /* 2. Hamiltonian. */
    double _Complex *H    = malloc(sizeof(double _Complex) * DIM * DIM);
    double _Complex *Hcpy = malloc(sizeof(double _Complex) * DIM * DIM);
    double *ev            = malloc(sizeof(double) * DIM);
    build_H(bi, bj, nb, /*J=*/1.0, H);
    memcpy(Hcpy, H, sizeof(double _Complex) * DIM * DIM);

    /* 3. Eigendecomposition — Jacobi destroys its input matrix, so we keep a
     * copy of H to reuse for the ground-state vector (obtained next via
     * shifted power iteration). */
    irrep_hermitian_eigvals(DIM, Hcpy, ev);

    /* rdm.h sorts eigenvalues descending (density-matrix convention); for a
     * Hamiltonian the ground state lives at the *bottom* of that list. */
    double E0_ed = ev[DIM - 1];
    double Emax  = ev[0];
    printf("  spectrum span:  E_max = %+.6f,  E_0 = %+.6f  (expect %+.6f)\n",
           Emax, E0_ed, -2.0);

    /* 4. Recover the ground-state vector by shifted power iteration on the
     * operator  M = Emax·I − H, whose dominant eigenvalue corresponds to
     * the smallest eigenvalue of H. Seed with a uniform Sz = 0 mix. */
    double _Complex gs[DIM] = {0};
    static const int sz0_states[6] = { 0x3, 0x5, 0x6, 0x9, 0xA, 0xC };
    for (int k = 0; k < 6; ++k) gs[sz0_states[k]] = 1.0;

    double shift = Emax + 1.0;
    for (int it = 0; it < 600; ++it) {
        double _Complex out[DIM] = {0};
        for (int r = 0; r < DIM; ++r) {
            double _Complex acc = shift * gs[r];
            for (int c = 0; c < DIM; ++c) acc -= H[r * DIM + c] * gs[c];
            out[r] = acc;
        }
        double norm = 0.0;
        for (int k = 0; k < DIM; ++k) norm += creal(out[k])*creal(out[k]) + cimag(out[k])*cimag(out[k]);
        norm = sqrt(norm);
        if (norm < 1e-300) break;
        for (int k = 0; k < DIM; ++k) gs[k] = out[k] / norm;
    }
    double E_check = 0.0;
    for (int r = 0; r < DIM; ++r) {
        double _Complex acc = 0.0;
        for (int c = 0; c < DIM; ++c) acc += H[r * DIM + c] * gs[c];
        E_check += creal(conj(gs[r]) * acc);
    }
    printf("⟨gs|H|gs⟩ = %+.6f  (should match E_0 = %+.6f)\n", E_check, E0_ed);

    /* 5. Total-J projection: project the ground state onto J=0, J=1, J=2.
     * Expected singlet weight ≈ 1. */
    double _Complex *pgs = malloc(sizeof(double _Complex) * DIM);
    double sumsq[3] = {0};
    for (int two_J = 0; two_J <= 4; two_J += 2) {
        irrep_spin_project_spin_half(two_J, N, 8, 6, 8, gs, pgs);
        double w = 0.0;
        for (int k = 0; k < DIM; ++k) w += creal(pgs[k])*creal(pgs[k]) + cimag(pgs[k])*cimag(pgs[k]);
        sumsq[two_J / 2] = w;
        printf("‖P_{J=%d} |gs⟩‖² = %.6f\n", two_J / 2, w);
    }
    printf("Σ_J ‖P_J |gs⟩‖² = %.6f   (expected 1.0 by completeness)\n",
           sumsq[0] + sumsq[1] + sumsq[2]);

    /* 6. Partial trace to a 2-site subsystem (sites 0, 1) and entropy. */
    int sites_A[2] = { 0, 1 };
    double _Complex rho_A[16];
    irrep_partial_trace(N, 2, gs, sites_A, 2, rho_A);
    double S = irrep_entropy_vonneumann(rho_A, 4);
    printf("S_VN at the {0,1} vs {2,3} bipartition = %.6f nats  (≈ %.4f bits)\n",
           S, S / log(2.0));

    /* Tr(ρ_A) = 1 sanity */
    double tr = 0.0;
    for (int k = 0; k < 4; ++k) tr += creal(rho_A[k * 4 + k]);
    printf("Tr(ρ_A) = %.6f   (expected 1.0)\n", tr);

    free(pgs);
    free(H); free(Hcpy); free(ev);
    free(bi); free(bj);
    irrep_lattice_free(L);
    return 0;
}
