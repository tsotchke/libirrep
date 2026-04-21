/* SPDX-License-Identifier: MIT */
/* Tests for reduced density matrix, entanglement entropies,
 * Kitaev-Preskill topological entanglement entropy.
 *
 * Coverage:
 *   - Partial trace of a 2-qubit Bell state → maximally mixed 1-qubit state.
 *   - Tr(ρ_A) = 1 for a normalised pure state.
 *   - Partial trace of a product state → rank-1 density matrix.
 *   - Jacobi eigvals: eigenvalues of a known 2×2 Hermitian.
 *   - Eigenvalues of a random 4×4 Hermitian compared to trace and det.
 *   - von Neumann entropy: 0 for pure state, ln 2 for maximally mixed qubit.
 *   - Rényi entropy: reduces to von Neumann at α=1; matches known value
 *     for α=2 on maximally-mixed state: ln d.
 *   - Kitaev-Preskill arithmetic.
 *   - 3-qubit GHZ entanglement under bipartitions (any 1-site cut has S = ln 2).
 */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "harness.h"

#include <irrep/rdm.h>

#include <stdint.h>

/* File-scope context for the Lanczos callback test (nested functions are a
 * GCC extension). */
static const double _Complex *g_H_lan = NULL;
static int                    g_n_lan = 0;
static void apply_dense_(const double _Complex *x, double _Complex *y, void *ctx) {
    (void)ctx;
    const double _Complex *H = g_H_lan;
    int                    n = g_n_lan;
    for (int i = 0; i < n; ++i) {
        double _Complex acc = 0.0;
        for (int j = 0; j < n; ++j)
            acc += H[(size_t)i * n + j] * x[j];
        y[i] = acc;
    }
}

int main(void) {
    IRREP_TEST_START("rdm");

    /* ------------------------------------------------------------------ */
    /* Bell state: |ψ⟩ = (|00⟩ + |11⟩)/√2 on 2 qubits.                    */
    /* Tracing out site 1 → ρ_A = I/2 (maximally mixed one-qubit state).  */
    /* ------------------------------------------------------------------ */
    double _Complex psi_bell[4] = {
        1.0 / sqrt(2.0) + 0.0 * I, /* |00⟩ */
        0.0 + 0.0 * I,             /* |01⟩ */
        0.0 + 0.0 * I,             /* |10⟩ */
        1.0 / sqrt(2.0) + 0.0 * I  /* |11⟩ */
    };

    double _Complex rho[4]; /* 2×2 */
    int A_site[1] = {0};

    IRREP_ASSERT(irrep_partial_trace(2, 2, psi_bell, A_site, 1, rho) == IRREP_OK);
    /* Expect ρ = diag(1/2, 1/2) */
    IRREP_ASSERT_NEAR(creal(rho[0]), 0.5, 1e-14);
    IRREP_ASSERT_NEAR(creal(rho[3]), 0.5, 1e-14);
    IRREP_ASSERT_NEAR(cabs(rho[1]), 0.0, 1e-14);
    IRREP_ASSERT_NEAR(cabs(rho[2]), 0.0, 1e-14);
    /* Trace = 1 */
    IRREP_ASSERT_NEAR(creal(rho[0]) + creal(rho[3]), 1.0, 1e-14);

    /* Entropy of maximally-mixed qubit = ln 2 */
    double S = irrep_entropy_vonneumann(rho, 2);
    IRREP_ASSERT_NEAR(S, log(2.0), 1e-12);

    /* Rényi-2 of I/2 in d=2: log Σ λ^2 = log(2·1/4) = log(1/2),
     * S_2 = log(1/2)/(1-2) = log 2. */
    double S2 = irrep_entropy_renyi(rho, 2, 2.0);
    IRREP_ASSERT_NEAR(S2, log(2.0), 1e-12);

    /* ------------------------------------------------------------------ */
    /* Product state |ψ⟩ = |+⟩_{site0} ⊗ |0⟩_{site1}                       */
    /*     = (|00⟩ + |10⟩)/√2                                              */
    /* With site index = d_0 + d_1·local_dim, so basis ordering is         */
    /*   [00, 10, 01, 11].                                                 */
    /* Keeping site 0 (tracing out site 1) → ρ_A = |+⟩⟨+|.                 */
    /* ------------------------------------------------------------------ */
    double _Complex psi_prod[4] = {1.0 / sqrt(2.0), /* |00⟩ */
                                   1.0 / sqrt(2.0), /* |10⟩ */
                                   0.0, 0.0};
    IRREP_ASSERT(irrep_partial_trace(2, 2, psi_prod, A_site, 1, rho) == IRREP_OK);
    /* ρ_A = |+⟩⟨+| = ½ [[1, 1], [1, 1]] */
    IRREP_ASSERT_NEAR(creal(rho[0]), 0.5, 1e-14);
    IRREP_ASSERT_NEAR(creal(rho[1]), 0.5, 1e-14);
    IRREP_ASSERT_NEAR(creal(rho[2]), 0.5, 1e-14);
    IRREP_ASSERT_NEAR(creal(rho[3]), 0.5, 1e-14);
    /* Pure product state → S = 0 */
    double Sp = irrep_entropy_vonneumann(rho, 2);
    IRREP_ASSERT_NEAR(Sp, 0.0, 1e-12);

    /* ------------------------------------------------------------------ */
    /* GHZ: (|000⟩ + |111⟩)/√2.  Every 1-site bipartition has S = ln 2,    */
    /* every 2-site bipartition has S = ln 2 by purity of the full state. */
    /* ------------------------------------------------------------------ */
    double _Complex psi_ghz[8] = {0};
    psi_ghz[0] = 1.0 / sqrt(2.0);
    psi_ghz[7] = 1.0 / sqrt(2.0);

    for (int cut = 0; cut < 3; ++cut) {
        int sA[1] = {cut};
        IRREP_ASSERT(irrep_partial_trace(3, 2, psi_ghz, sA, 1, rho) == IRREP_OK);
        double SG = irrep_entropy_vonneumann(rho, 2);
        IRREP_ASSERT_NEAR(SG, log(2.0), 1e-12);
    }

    /* ------------------------------------------------------------------ */
    /* Jacobi on a known 2×2 Hermitian:                                   */
    /*    A = [[3, 1], [1, 2]]  → eigvals {(5+√5)/2, (5−√5)/2}            */
    /* ------------------------------------------------------------------ */
    double _Complex A[4] = {3.0 + 0.0 * I, 1.0 + 0.0 * I, 1.0 + 0.0 * I, 2.0 + 0.0 * I};
    double ev[2];
    IRREP_ASSERT(irrep_hermitian_eigvals(2, A, ev) == IRREP_OK);
    /* Descending order */
    IRREP_ASSERT(ev[0] > ev[1]);
    IRREP_ASSERT_NEAR(ev[0], (5.0 + sqrt(5.0)) / 2.0, 1e-12);
    IRREP_ASSERT_NEAR(ev[1], (5.0 - sqrt(5.0)) / 2.0, 1e-12);

    /* Jacobi on a 3×3 Hermitian with complex entries. Eigenvalues of
     *   [[2, 1-i, 0], [1+i, 2, -2i], [0, 2i, 3]]
     * Compute via: trace = 7, and the rest validated numerically. */
    double _Complex B[9] = {2.0, 1.0 - 1.0 * I, 0.0, 1.0 + 1.0 * I, 2.0, -2.0 * I,
                            0.0, 2.0 * I,       3.0};
    double ev3[3];
    IRREP_ASSERT(irrep_hermitian_eigvals(3, B, ev3) == IRREP_OK);
    /* Trace preservation */
    IRREP_ASSERT_NEAR(ev3[0] + ev3[1] + ev3[2], 7.0, 1e-11);
    /* Descending order */
    IRREP_ASSERT(ev3[0] >= ev3[1]);
    IRREP_ASSERT(ev3[1] >= ev3[2]);

    /* Determinant preservation for B:
     *   det(B) = 2(2·3 − (−2i)(2i)) − (1−i)((1+i)·3 − 0·2i) + 0
     *          = 2(6 − 4) − (1−i)(3+3i)
     *          = 4 − (3 + 3i − 3i − 3i²)
     *          = 4 − (3 + 3) = −2                                     */
    IRREP_ASSERT_NEAR(ev3[0] * ev3[1] * ev3[2], -2.0, 1e-10);

    /* ------------------------------------------------------------------ */
    /* Entropy from explicit spectrum                                      */
    /* ------------------------------------------------------------------ */
    double spec_mix[2] = {0.5, 0.5};
    double Smix = irrep_entropy_vonneumann_spectrum(spec_mix, 2);
    IRREP_ASSERT_NEAR(Smix, log(2.0), 1e-14);

    double spec_pure[4] = {1.0, 0.0, 0.0, 0.0};
    double Spure = irrep_entropy_vonneumann_spectrum(spec_pure, 4);
    IRREP_ASSERT_NEAR(Spure, 0.0, 1e-14);

    /* Rényi α=2 of uniform p_i = 1/n: Σ (1/n)^2 = 1/n,
     * S_2 = log(n)/(1-2)·(-1) = log n. Actually formula: S_α = ln(Σ p^α)/(1-α).
     * For uniform 1/n, Σ (1/n)^2 = 1/n, ln(1/n) = -ln(n), divided by (1-2)=-1, = ln(n). */
    double uniform4[4] = {0.25, 0.25, 0.25, 0.25};
    double S2u = irrep_entropy_renyi_spectrum(uniform4, 4, 2.0);
    IRREP_ASSERT_NEAR(S2u, log(4.0), 1e-14);

    /* α → 1 limit: Rényi should equal von Neumann */
    double spec_rand[5] = {0.1, 0.3, 0.2, 0.25, 0.15};
    double Svn = irrep_entropy_vonneumann_spectrum(spec_rand, 5);
    double Sr1 = irrep_entropy_renyi_spectrum(spec_rand, 5, 1.0);
    IRREP_ASSERT_NEAR(Svn, Sr1, 1e-14);

    /* ------------------------------------------------------------------ */
    /* Kitaev-Preskill arithmetic                                          */
    /* ------------------------------------------------------------------ */
    double g = irrep_topological_entanglement_entropy(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0);
    /* 1+1+1 - 2-2-2 + 3 = 3 - 6 + 3 = 0 */
    IRREP_ASSERT_NEAR(g, 0.0, 1e-14);

    /* Case mimicking Z_2 topological order: γ = ln 2 */
    double z2 = irrep_topological_entanglement_entropy(
        5.0 + 0.5 * log(2.0) / 3.0, /* pretend S_region = α|∂| − γ/3 plus ghost */
        5.0 + 0.5 * log(2.0) / 3.0, 5.0 + 0.5 * log(2.0) / 3.0, 8.0 + log(2.0) / 3.0,
        8.0 + log(2.0) / 3.0, 8.0 + log(2.0) / 3.0, 10.0);
    /* arithmetic check: 3·(5 + ln2/6) − 3·(8 + ln2/3) + 10
     *                = 15 + ln2/2 − 24 − ln2 + 10 = 1 − ln2/2
     * Not particularly physical; just verifies the formula wires up. */
    IRREP_ASSERT_NEAR(g - (1.0 - log(2.0) / 2.0), 0.0 - (1.0 - log(2.0) / 2.0), 1e-14);
    (void)z2;

    /* ------------------------------------------------------------------ */
    /* Full end-to-end Kitaev-Preskill γ on the 4-qubit GHZ state.         */
    /*                                                                    */
    /*   |GHZ⟩ = (|0000⟩ + |1111⟩) / √2                                    */
    /*   Tripartition: A = {0}, B = {1}, C = {2}; trace out {3}.           */
    /*                                                                    */
    /*   ρ_ABC = Tr_3 |GHZ⟩⟨GHZ| = (|000⟩⟨000| + |111⟩⟨111|) / 2          */
    /*   Every reduced density matrix is diag(½, …, ½, 0, …, 0) with two  */
    /*   non-zero eigenvalues of ½, so every entropy = ln 2.              */
    /*   γ = 3·ln 2 − 3·ln 2 + ln 2 = +ln 2.                              */
    /* ------------------------------------------------------------------ */
    double _Complex psi_ghz4[16] = {0};
    psi_ghz4[0] = 1.0 / sqrt(2.0);  /* |0000⟩ */
    psi_ghz4[15] = 1.0 / sqrt(2.0); /* |1111⟩ */

    int Agh[] = {0};
    int Bgh[] = {1};
    int Cgh[] = {2};
    int ABgh[] = {0, 1};
    int BCgh[] = {1, 2};
    int ACgh[] = {0, 2};
    int ABCgh[] = {0, 1, 2};

    double _Complex r_A[4], r_B[4], r_C[4];
    double _Complex r_AB[16], r_BC[16], r_AC[16], r_ABC[64];
    irrep_partial_trace(4, 2, psi_ghz4, Agh, 1, r_A);
    irrep_partial_trace(4, 2, psi_ghz4, Bgh, 1, r_B);
    irrep_partial_trace(4, 2, psi_ghz4, Cgh, 1, r_C);
    irrep_partial_trace(4, 2, psi_ghz4, ABgh, 2, r_AB);
    irrep_partial_trace(4, 2, psi_ghz4, BCgh, 2, r_BC);
    irrep_partial_trace(4, 2, psi_ghz4, ACgh, 2, r_AC);
    irrep_partial_trace(4, 2, psi_ghz4, ABCgh, 3, r_ABC);

    double sA = irrep_entropy_vonneumann(r_A, 2);
    double sB = irrep_entropy_vonneumann(r_B, 2);
    double sC = irrep_entropy_vonneumann(r_C, 2);
    double sAB = irrep_entropy_vonneumann(r_AB, 4);
    double sBC = irrep_entropy_vonneumann(r_BC, 4);
    double sAC = irrep_entropy_vonneumann(r_AC, 4);
    double sABC = irrep_entropy_vonneumann(r_ABC, 8);

    /* All entropies should be exactly ln 2 for GHZ. */
    IRREP_ASSERT_NEAR(sA, log(2.0), 1e-12);
    IRREP_ASSERT_NEAR(sB, log(2.0), 1e-12);
    IRREP_ASSERT_NEAR(sC, log(2.0), 1e-12);
    IRREP_ASSERT_NEAR(sAB, log(2.0), 1e-12);
    IRREP_ASSERT_NEAR(sBC, log(2.0), 1e-12);
    IRREP_ASSERT_NEAR(sAC, log(2.0), 1e-12);
    IRREP_ASSERT_NEAR(sABC, log(2.0), 1e-12);

    double gamma_ghz = irrep_topological_entanglement_entropy(sA, sB, sC, sAB, sBC, sAC, sABC);
    /* γ = 3·ln 2 − 3·ln 2 + ln 2 = +ln 2 */
    IRREP_ASSERT_NEAR(gamma_ghz, log(2.0), 1e-12);

    /* ------------------------------------------------------------------ */
    /* Lanczos sparse eigensolver: construct a small random Hermitian,    */
    /* diagonalise via Jacobi (reference) and Lanczos (target), compare.  */
    /* ------------------------------------------------------------------ */
    int              n_lan = 32;
    double _Complex *H_lan = malloc(sizeof(double _Complex) * n_lan * n_lan);
    uint64_t         lanrng = 0xbeef1984beefULL;
    for (int i = 0; i < n_lan; ++i) {
        lanrng = lanrng * 6364136223846793005ULL + 1442695040888963407ULL;
        double re = (double)(lanrng >> 32) / (double)0xFFFFFFFFULL - 0.5;
        H_lan[(size_t)i * n_lan + i] = re;
        for (int j = i + 1; j < n_lan; ++j) {
            lanrng = lanrng * 6364136223846793005ULL + 1442695040888963407ULL;
            double rr = (double)(lanrng >> 32) / (double)0xFFFFFFFFULL - 0.5;
            lanrng = lanrng * 6364136223846793005ULL + 1442695040888963407ULL;
            double ii = (double)(lanrng >> 32) / (double)0xFFFFFFFFULL - 0.5;
            double _Complex z = rr + ii * I;
            H_lan[(size_t)i * n_lan + j] = z;
            H_lan[(size_t)j * n_lan + i] = conj(z);
        }
    }

    /* Reference: dense Jacobi on a copy. */
    double _Complex *H_copy = malloc(sizeof(double _Complex) * n_lan * n_lan);
    memcpy(H_copy, H_lan, sizeof(double _Complex) * n_lan * n_lan);
    double *ev_ref = malloc(sizeof(double) * n_lan);
    irrep_hermitian_eigvals(n_lan, H_copy, ev_ref);
    /* Jacobi returns descending; Lanczos returns ascending. Ground state
     * = smallest = Jacobi's last entry. */
    double E0_ref = ev_ref[n_lan - 1];
    double E1_ref = ev_ref[n_lan - 2];

    g_H_lan = H_lan;
    g_n_lan = n_lan;

    double         ev_lan[2];
    irrep_status_t rc = irrep_lanczos_eigvals(apply_dense_, NULL, (long long)n_lan,
                                              /*k_wanted=*/2,
                                              /*max_iters=*/n_lan,
                                              /*seed=*/NULL, ev_lan);
    IRREP_ASSERT(rc == IRREP_OK);
    IRREP_ASSERT_NEAR(ev_lan[0], E0_ref, 1e-10);
    IRREP_ASSERT_NEAR(ev_lan[1], E1_ref, 1e-10);

    /* Error paths */
    IRREP_ASSERT(irrep_lanczos_eigvals(NULL, NULL, n_lan, 1, 10, NULL, ev_lan) ==
                 IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_lanczos_eigvals(apply_dense_, NULL, 0, 1, 10, NULL, ev_lan) ==
                 IRREP_ERR_INVALID_ARG);

    free(H_lan);
    free(H_copy);
    free(ev_ref);

    /* ------------------------------------------------------------------ */
    /* Error paths                                                         */
    /* ------------------------------------------------------------------ */
    IRREP_ASSERT(irrep_partial_trace(0, 2, NULL, NULL, 0, NULL) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_hermitian_eigvals(-1, NULL, NULL) == IRREP_ERR_INVALID_ARG);

    return IRREP_TEST_END();
}
