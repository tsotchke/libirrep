/* SPDX-License-Identifier: MIT */
/* Tests for total-J projection on spin-½ systems.
 *
 * Coverage:
 *   - Singlet state (|↑↓⟩ - |↓↑⟩)/√2 on 2 sites: J=0 projection is the
 *     identity on this state; J=1 projection is zero.
 *   - Triplet state (|↑↓⟩ + |↓↑⟩)/√2 on 2 sites: J=1 projection is
 *     identity (on this m=0 triplet component); J=0 projection is zero.
 *   - Trivial 1-qubit case: β=0 rotation leaves the state invariant.
 *   - Idempotence: P(P ψ) = P ψ for J=0 on 2 sites.
 *   - Completeness: sum over J={0,1} for 2-site state = original ψ.
 *   - Spin-½ rotation correctness: applying D^{½}(α=π, β=0, γ=0) to |↑⟩
 *     produces a phase on the state (−i |↑⟩).
 */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "harness.h"

#include <irrep/spin_project.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double norm2_(const double _Complex *psi, long long dim) {
    double s = 0.0;
    for (long long k = 0; k < dim; ++k) {
        double a = creal(psi[k]);
        double b = cimag(psi[k]);
        s += a * a + b * b;
    }
    return s;
}

static void normalize_(double _Complex *psi, long long dim) {
    double n = sqrt(norm2_(psi, dim));
    if (n > 1e-300)
        for (long long k = 0; k < dim; ++k)
            psi[k] /= n;
}

int main(void) {
    IRREP_TEST_START("spin_project");

    /* ------------------------------------------------------------------ */
    /* Single-qubit rotations: D^{½}(π, 0, 0) should produce e^{-iπ/2}|↑⟩ */
    /* = -i |↑⟩ (up to the ZYZ sign convention).                          */
    /* ------------------------------------------------------------------ */
    double _Complex up[2] = {1.0 + 0.0 * I, 0.0 + 0.0 * I};
    double _Complex out[2] = {0};
    irrep_spin_half_apply_rotation(1, M_PI, 0.0, 0.0, up, out);
    /* D^{½}(α, 0, 0) = diag(e^{-iα/2}, e^{+iα/2}).  At α=π: (-i, +i) */
    IRREP_ASSERT_NEAR(creal(out[0]), 0.0, 1e-14);
    IRREP_ASSERT_NEAR(cimag(out[0]), -1.0, 1e-14);
    IRREP_ASSERT_NEAR(cabs(out[1]), 0.0, 1e-14);

    /* β = 0 rotation preserves state up to an α+γ phase — so its norm is preserved. */
    double _Complex rand_state[4] = {0.3 + 0.1 * I, 0.5 - 0.2 * I, -0.1 + 0.4 * I, 0.2 + 0.3 * I};
    double n0 = norm2_(rand_state, 4);
    double _Complex rand_out[4];
    irrep_spin_half_apply_rotation(2, 0.7, 0.0, 0.3, rand_state, rand_out);
    IRREP_ASSERT_NEAR(norm2_(rand_out, 4), n0, 1e-12);

    /* Arbitrary rotation still preserves the norm of the full wavefunction. */
    irrep_spin_half_apply_rotation(2, 0.7, 1.1, 0.3, rand_state, rand_out);
    IRREP_ASSERT_NEAR(norm2_(rand_out, 4), n0, 1e-12);

    /* ------------------------------------------------------------------ */
    /* 2-site singlet: ψ_s = (|↑↓⟩ − |↓↑⟩)/√2                               */
    /* Basis ordering: |↑↑⟩=0, |↓↑⟩=1, |↑↓⟩=2, |↓↓⟩=3  (qubit 0 low bit,   */
    /* |↑⟩ = 0, |↓⟩ = 1 in computational basis).                           */
    /* P_{J=0} ψ_s = ψ_s; P_{J=1} ψ_s = 0.                                 */
    /* ------------------------------------------------------------------ */
    double _Complex singlet[4] = {0};
    singlet[1] = -1.0 / sqrt(2.0); /* |↓↑⟩ */
    singlet[2] = 1.0 / sqrt(2.0);  /* |↑↓⟩ */

    double _Complex proj[4] = {0};

    /* J=0 projection */
    IRREP_ASSERT(irrep_spin_project_spin_half(0, 2, 8, 6, 8, singlet, proj) == IRREP_OK);
    /* Should be ≈ singlet */
    for (int k = 0; k < 4; ++k) {
        IRREP_ASSERT_NEAR(creal(proj[k]), creal(singlet[k]), 1e-10);
        IRREP_ASSERT_NEAR(cimag(proj[k]), cimag(singlet[k]), 1e-10);
    }
    /* ‖P₀ ψ_s‖² = 1 */
    IRREP_ASSERT_NEAR(norm2_(proj, 4), 1.0, 1e-10);

    /* J=1 projection of singlet → zero */
    IRREP_ASSERT(irrep_spin_project_spin_half(2, 2, 8, 6, 8, singlet, proj) == IRREP_OK);
    IRREP_ASSERT_NEAR(norm2_(proj, 4), 0.0, 1e-10);

    /* ------------------------------------------------------------------ */
    /* 2-site triplet m=0: ψ_t = (|↑↓⟩ + |↓↑⟩)/√2                          */
    /* P_{J=1} ψ_t = ψ_t; P_{J=0} ψ_t = 0                                  */
    /* ------------------------------------------------------------------ */
    double _Complex triplet[4] = {0};
    triplet[1] = 1.0 / sqrt(2.0);
    triplet[2] = 1.0 / sqrt(2.0);

    IRREP_ASSERT(irrep_spin_project_spin_half(2, 2, 8, 6, 8, triplet, proj) == IRREP_OK);
    IRREP_ASSERT_NEAR(norm2_(proj, 4), 1.0, 1e-10);
    /* Should be ≈ triplet */
    for (int k = 0; k < 4; ++k) {
        IRREP_ASSERT_NEAR(creal(proj[k]), creal(triplet[k]), 1e-10);
        IRREP_ASSERT_NEAR(cimag(proj[k]), cimag(triplet[k]), 1e-10);
    }

    IRREP_ASSERT(irrep_spin_project_spin_half(0, 2, 8, 6, 8, triplet, proj) == IRREP_OK);
    IRREP_ASSERT_NEAR(norm2_(proj, 4), 0.0, 1e-10);

    /* ------------------------------------------------------------------ */
    /* |↑↑⟩ is a triplet m=+1 eigenstate of J² with J=1; no J=0 component. */
    /* ------------------------------------------------------------------ */
    double _Complex uu[4] = {1.0, 0.0, 0.0, 0.0};
    IRREP_ASSERT(irrep_spin_project_spin_half(0, 2, 8, 6, 8, uu, proj) == IRREP_OK);
    IRREP_ASSERT_NEAR(norm2_(proj, 4), 0.0, 1e-10);
    IRREP_ASSERT(irrep_spin_project_spin_half(2, 2, 8, 6, 8, uu, proj) == IRREP_OK);
    IRREP_ASSERT_NEAR(norm2_(proj, 4), 1.0, 1e-10);

    /* ------------------------------------------------------------------ */
    /* Idempotence: P_J (P_J ψ) = P_J ψ                                    */
    /* ------------------------------------------------------------------ */
    double _Complex mix[4] = {0};
    mix[0] = 0.3;
    mix[1] = 0.5;
    mix[2] = 0.4;
    mix[3] = 0.2;
    normalize_(mix, 4);
    double _Complex proj1[4] = {0}, proj2[4] = {0};
    IRREP_ASSERT(irrep_spin_project_spin_half(0, 2, 8, 6, 8, mix, proj1) == IRREP_OK);
    IRREP_ASSERT(irrep_spin_project_spin_half(0, 2, 8, 6, 8, proj1, proj2) == IRREP_OK);
    for (int k = 0; k < 4; ++k) {
        IRREP_ASSERT_NEAR(creal(proj1[k]), creal(proj2[k]), 1e-9);
        IRREP_ASSERT_NEAR(cimag(proj1[k]), cimag(proj2[k]), 1e-9);
    }

    /* ------------------------------------------------------------------ */
    /* Completeness on 2 sites: J ∈ {0, 1}. P_0 ψ + P_1 ψ = ψ.             */
    /* ------------------------------------------------------------------ */
    double _Complex pJ0[4] = {0}, pJ1[4] = {0};
    IRREP_ASSERT(irrep_spin_project_spin_half(0, 2, 8, 6, 8, mix, pJ0) == IRREP_OK);
    IRREP_ASSERT(irrep_spin_project_spin_half(2, 2, 8, 6, 8, mix, pJ1) == IRREP_OK);
    for (int k = 0; k < 4; ++k) {
        double _Complex sum = pJ0[k] + pJ1[k];
        IRREP_ASSERT_NEAR(creal(sum), creal(mix[k]), 1e-9);
        IRREP_ASSERT_NEAR(cimag(sum), cimag(mix[k]), 1e-9);
    }

    /* ------------------------------------------------------------------ */
    /* Error paths                                                          */
    /* ------------------------------------------------------------------ */
    IRREP_ASSERT(irrep_spin_project_spin_half(-1, 2, 1, 1, 1, singlet, proj) ==
                 IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_spin_project_spin_half(0, 0, 1, 1, 1, singlet, proj) ==
                 IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_spin_project_spin_half(0, 2, 0, 1, 1, singlet, proj) ==
                 IRREP_ERR_INVALID_ARG);

    return IRREP_TEST_END();
}
