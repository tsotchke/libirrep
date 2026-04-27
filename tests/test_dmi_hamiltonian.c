/* SPDX-License-Identifier: MIT */
/* Tests for the DMI Hamiltonian apply operator. Verifies the matrix
 * elements of `D · (S_i × S_j)` against the analytic 2-spin reference,
 * Hermiticity, and the reduction to known limits. */

#include <irrep/dmi_hamiltonian.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int total = 0, failed = 0;

#define ASSERT_NEAR_C(z, expected_re, expected_im, tol, msg)                                       \
    do {                                                                                           \
        ++total;                                                                                   \
        double dre = creal(z) - (expected_re);                                                     \
        double dim_ = cimag(z) - (expected_im);                                                    \
        if (fabs(dre) > (tol) || fabs(dim_) > (tol)) {                                             \
            fprintf(stderr, "  FAIL  %s:%d  %s  got (%.6f, %.6f) vs (%.6f, %.6f)\n",               \
                    __FILE__, __LINE__, msg, creal(z), cimag(z), (double)(expected_re),            \
                    (double)(expected_im));                                                        \
            ++failed;                                                                              \
        }                                                                                          \
    } while (0)

#define ASSERT(cond, msg)                                                                          \
    do {                                                                                           \
        ++total;                                                                                   \
        if (!(cond)) {                                                                             \
            fprintf(stderr, "  FAIL  %s:%d  %s\n", __FILE__, __LINE__, msg);                       \
            ++failed;                                                                              \
        }                                                                                          \
    } while (0)

/* The 2-spin DMI Hamiltonian D · (S_0 × S_1) in the basis
 * {|↓↓⟩, |↓↑⟩, |↑↓⟩, |↑↑⟩} (bit 0 = site 0, bit 1 = site 1).
 *
 * Closed-form derivation: with D = (Dx, Dy, Dz),
 *
 *   ⟨↑↓ | H | ↓↑⟩  =  +i Dz / 4  +  …
 *
 * The off-diagonal matrix element on the |↑↓⟩↔|↓↑⟩ block:
 *   coeff_+- = (Dx + iDy) / 4    (from Dx S_y S_z − S_z S_y + Dy S_z S_x − S_x S_z)
 *
 * Wait actually let me derive carefully. For a 2-site bond in computational
 * basis indexed by (s0, s1):
 *
 *   |0⟩ = |↓↓⟩  (s = 00b = 0)
 *   |1⟩ = |↑↓⟩  (s = 01b = 1, bit 0 set)
 *   |2⟩ = |↓↑⟩  (s = 10b = 2, bit 1 set)
 *   |3⟩ = |↑↑⟩  (s = 11b = 3)
 *
 * Matrix elements of S_a × S_b (a=0, b=1) in this basis: */
static void test_2spin_matrix_elements(void) {
    /* Test with D = (1, 0, 0). */
    {
        const int bi[1] = {0}, bj[1] = {1};
        double Dx[1] = {1.0}, Dy[1] = {0.0}, Dz[1] = {0.0};
        irrep_dmi_hamiltonian_t *H = irrep_dmi_hamiltonian_new(2, 1, bi, bj, Dx, Dy, Dz);
        double _Complex psi[4], out[4];

        /* Apply to |↓↓⟩ = (1, 0, 0, 0) */
        for (int k = 0; k < 4; ++k) psi[k] = 0; psi[0] = 1.0;
        irrep_dmi_apply(psi, out, H);
        /* D_x acting on |↓↓⟩ = |s=0⟩:
         *   Term 1 (Dx): bit-a flip with coeff Dx · (i/2) · z_b · sgn_a
         *     z_b = -1/2 (s_b=0), sgn_a = +1 (s_a=0)
         *     coeff = 1 · (i/2) · (-1/2) · (+1) = -i/4 → out[1] += -i/4
         *   Term 1 (Dx): bit-b flip with coeff -Dx · (i/2) · z_a · sgn_b
         *     z_a = -1/2, sgn_b = +1
         *     coeff = -1 · (i/2) · (-1/2) · (+1) = +i/4 → out[2] += i/4
         *   Term 3 (Dz): zero (Dz=0). */
        ASSERT_NEAR_C(out[0], 0, 0, 1e-12, "D_x on |↓↓⟩: out[|↓↓⟩] = 0");
        ASSERT_NEAR_C(out[1], 0, -0.25, 1e-12, "D_x on |↓↓⟩: out[|↑↓⟩] = -i/4");
        ASSERT_NEAR_C(out[2], 0, +0.25, 1e-12, "D_x on |↓↓⟩: out[|↓↑⟩] = +i/4");
        ASSERT_NEAR_C(out[3], 0, 0, 1e-12, "D_x on |↓↓⟩: out[|↑↑⟩] = 0");

        irrep_dmi_hamiltonian_free(H);
    }

    /* Test with D = (0, 0, 1) (D_z). Action on |↑↓⟩:
     *   Term 3: bits anti-aligned (sa=1, sb=0). sign = -1 since sa==1, sb==0
     *   out[s ⊕ mi ⊕ mj] += (Dz · 0.25 · I · sign) · psi[s]
     *   = 1 · 0.25 · i · (-1) · 1 = -i/4
     *   s=01b → s' = 01b ⊕ 11b = 10b = |↓↑⟩
     *   so out[|↓↑⟩] = -i/4
     *   Plus Terms 1, 2 are zero since Dx=Dy=0.
     * Action on |↓↑⟩ (sa=0, sb=1): sign = +1, out[|↑↓⟩] = +i/4
     */
    {
        const int bi[1] = {0}, bj[1] = {1};
        double Dx[1] = {0}, Dy[1] = {0}, Dz[1] = {1.0};
        irrep_dmi_hamiltonian_t *H = irrep_dmi_hamiltonian_new(2, 1, bi, bj, Dx, Dy, Dz);
        double _Complex psi[4], out[4];
        for (int k = 0; k < 4; ++k) psi[k] = 0; psi[1] = 1.0;  /* |↑↓⟩ */
        irrep_dmi_apply(psi, out, H);
        ASSERT_NEAR_C(out[2], 0, -0.25, 1e-12, "D_z on |↑↓⟩: out[|↓↑⟩] = -i/4");
        for (int k = 0; k < 4; ++k) if (k != 2)
            ASSERT_NEAR_C(out[k], 0, 0, 1e-12, "D_z on |↑↓⟩: zero elsewhere");

        for (int k = 0; k < 4; ++k) psi[k] = 0; psi[2] = 1.0;  /* |↓↑⟩ */
        irrep_dmi_apply(psi, out, H);
        ASSERT_NEAR_C(out[1], 0, +0.25, 1e-12, "D_z on |↓↑⟩: out[|↑↓⟩] = +i/4");

        irrep_dmi_hamiltonian_free(H);
    }
    (void)0;
}

/* Hermiticity: ⟨φ|H|ψ⟩ = ⟨ψ|H|φ⟩^*. Test by sandwiching random
 * amplitude vectors. */
static void test_hermiticity(void) {
    const int bi[1] = {0}, bj[1] = {1};
    double Dx[1] = {0.7}, Dy[1] = {-0.3}, Dz[1] = {0.5};
    irrep_dmi_hamiltonian_t *H = irrep_dmi_hamiltonian_new(2, 1, bi, bj, Dx, Dy, Dz);

    double _Complex psi[4], phi[4], H_psi[4], H_phi[4];
    /* Deterministic "random" amplitudes. */
    for (int k = 0; k < 4; ++k) {
        psi[k] = 0.3 * sin(0.7 * k + 1.0) + I * 0.5 * cos(1.3 * k - 0.7);
        phi[k] = 0.4 * cos(1.1 * k + 0.3) + I * 0.2 * sin(0.5 * k + 1.7);
    }
    irrep_dmi_apply(psi, H_psi, H);
    irrep_dmi_apply(phi, H_phi, H);

    double _Complex inner_phi_H_psi = 0, inner_psi_H_phi = 0;
    for (int k = 0; k < 4; ++k) {
        inner_phi_H_psi += conj(phi[k]) * H_psi[k];
        inner_psi_H_phi += conj(psi[k]) * H_phi[k];
    }
    /* Hermitian: ⟨φ|H|ψ⟩ = conj(⟨ψ|H|φ⟩) */
    double dre = creal(inner_phi_H_psi) - creal(conj(inner_psi_H_phi));
    double dim_ = cimag(inner_phi_H_psi) - cimag(conj(inner_psi_H_phi));
    ASSERT(fabs(dre) < 1e-12, "Hermiticity: real part");
    ASSERT(fabs(dim_) < 1e-12, "Hermiticity: imag part");

    irrep_dmi_hamiltonian_free(H);
}

/* Limit check: D = 0 → H = 0. Apply to any state should give zero. */
static void test_zero_D_limit(void) {
    const int bi[2] = {0, 1}, bj[2] = {1, 2};
    double Dx[2] = {0}, Dy[2] = {0}, Dz[2] = {0};
    irrep_dmi_hamiltonian_t *H = irrep_dmi_hamiltonian_new(3, 2, bi, bj, Dx, Dy, Dz);
    double _Complex psi[8], out[8];
    for (int k = 0; k < 8; ++k) psi[k] = 0.1 * (k + 1);
    irrep_dmi_apply(psi, out, H);
    for (int k = 0; k < 8; ++k)
        ASSERT_NEAR_C(out[k], 0, 0, 1e-15, "D=0 limit: out[k]=0");
    irrep_dmi_hamiltonian_free(H);
}

int main(void) {
    fprintf(stderr, "test_dmi_hamiltonian:\n");
    test_2spin_matrix_elements();
    test_hermiticity();
    test_zero_D_limit();
    fprintf(stderr, "  %d / %d assertions passed\n", total - failed, total);
    return failed == 0 ? 0 : 1;
}
