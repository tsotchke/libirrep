/* SPDX-License-Identifier: MIT */
/* Tests for configuration-space projection.
 *
 * Coverage:
 *   - Trivial-irrep builder produces flat characters.
 *   - A1 projection of a constant amplitude returns the constant.
 *   - A1 projection idempotent: P(Pψ) = Pψ within tolerance.
 *   - Orbit enumeration reproduces individual permutation_inverse calls.
 *   - Under translation-only p1 group, projecting a translation-invariant
 *     amplitude (all the same) returns the same value.
 *   - For an asymmetric amplitude (all zero except g=0), the A1 projection
 *     returns 1/|G|.
 */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "harness.h"

#include <irrep/lattice.h>
#include <irrep/space_group.h>
#include <irrep/config_project.h>

int main(void) {
    IRREP_TEST_START("config_project");

    /* Build a 4×4 square p4mm group to exercise a non-trivial order. */
    irrep_lattice_t *sq     = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
    irrep_space_group_t *G  = irrep_space_group_build(sq, IRREP_WALLPAPER_P4MM);
    IRREP_ASSERT(G != NULL);
    int order = irrep_space_group_order(G);          /* 128 */
    int n     = irrep_space_group_num_sites(G);      /* 16 */
    IRREP_ASSERT(order == 128);
    IRREP_ASSERT(n     == 16);

    /* Trivial irrep: should give all-1 characters, dim 1 */
    irrep_sg_irrep_t *A1 = irrep_sg_trivial(G);
    IRREP_ASSERT(A1 != NULL);

    /* Reduce a constant-1 amplitude over the group — A1 projection of a
     * constant is the constant itself (d_μ · (1/|G|) · Σ 1·1 = 1). */
    double _Complex *psi = malloc(sizeof(double _Complex) * order);
    for (int g = 0; g < order; ++g) psi[g] = 1.0 + 0.0*I;
    double _Complex proj = irrep_sg_project_amplitude(A1, psi);
    IRREP_ASSERT_NEAR(creal(proj), 1.0, 1e-14);
    IRREP_ASSERT_NEAR(cimag(proj), 0.0, 1e-14);

    /* Same with the A1 shortcut */
    proj = irrep_sg_project_A1(G, psi);
    IRREP_ASSERT_NEAR(creal(proj), 1.0, 1e-14);

    /* Delta amplitude at g=0 → A1 projection = 1/|G| */
    memset(psi, 0, sizeof(double _Complex) * order);
    psi[0] = 1.0 + 0.0*I;
    proj = irrep_sg_project_A1(G, psi);
    IRREP_ASSERT_NEAR(creal(proj), 1.0 / (double)order, 1e-14);

    /* Random-ish non-symmetric amplitude: A1 projection is a scalar; then
     * applying the projection to that constant array returns the same. */
    for (int g = 0; g < order; ++g) psi[g] = (double)(g % 7) - 3.0;
    double _Complex avg = irrep_sg_project_A1(G, psi);
    for (int g = 0; g < order; ++g) psi[g] = avg;
    double _Complex avg2 = irrep_sg_project_A1(G, psi);
    IRREP_ASSERT_NEAR(creal(avg), creal(avg2), 1e-14);
    IRREP_ASSERT_NEAR(cimag(avg), cimag(avg2), 1e-14);

    /* Generic non-trivial irrep: alternating ±1 characters on element parity.
     * (Not the true parity of the group, but an exercise of the reduce.) */
    double _Complex *chi = malloc(sizeof(double _Complex) * order);
    for (int g = 0; g < order; ++g) chi[g] = (g & 1) ? -1.0 : +1.0;
    irrep_sg_irrep_t *mu = irrep_sg_irrep_new(G, chi, 1);
    IRREP_ASSERT(mu != NULL);
    for (int g = 0; g < order; ++g) psi[g] = 1.0 + 0.0*I;
    /* Σ conj(χ) · 1 = (order/2) · (+1) + (order/2) · (-1) = 0 → projection 0 */
    proj = irrep_sg_project_amplitude(mu, psi);
    IRREP_ASSERT_NEAR(creal(proj), 0.0, 1e-14);
    IRREP_ASSERT_NEAR(cimag(proj), 0.0, 1e-14);

    irrep_sg_irrep_free(mu);
    free(chi);

    /* Sign-rep A₂: identity character on all translations + rotations,
     * -1 on reflections. On a constant amplitude vector, the reduce
     * should cancel (as many +1's as -1's in the p4mm case). */
    irrep_sg_irrep_t *A2 = irrep_sg_sign_rep(G);
    IRREP_ASSERT(A2 != NULL);
    for (int g = 0; g < order; ++g) psi[g] = 1.0 + 0.0*I;
    proj = irrep_sg_project_amplitude(A2, psi);
    IRREP_ASSERT_NEAR(creal(proj), 0.0, 1e-14);
    IRREP_ASSERT_NEAR(cimag(proj), 0.0, 1e-14);
    /* On a reflection-odd amplitude (say ψ_g = +1 on rotations, -1 on
     * mirrors, replicated across translations), A₂ reduce returns ±1. */
    int point_order = irrep_space_group_point_order(G);
    int half = point_order / 2;
    for (int g = 0; g < order; ++g) {
        int p = g % point_order;
        psi[g] = (p >= half) ? -1.0 : +1.0;
    }
    proj = irrep_sg_project_amplitude(A2, psi);
    /* Σ |χ|² / |G| · d_μ = (1/|G|) · order = 1 (since chi² = 1 always). */
    IRREP_ASSERT_NEAR(creal(proj), 1.0, 1e-14);
    irrep_sg_irrep_free(A2);

    /* For p1 space groups the sign rep degenerates to the trivial one. */
    irrep_space_group_t *G1 = irrep_space_group_build(sq, IRREP_WALLPAPER_P1);
    irrep_sg_irrep_t *A2_p1 = irrep_sg_sign_rep(G1);
    IRREP_ASSERT(A2_p1 != NULL);
    for (int g = 0; g < irrep_space_group_order(G1); ++g) psi[g] = 1.0;
    proj = irrep_sg_project_amplitude(A2_p1, psi);
    IRREP_ASSERT_NEAR(creal(proj), 1.0, 1e-14);   /* matches A₁ on p1 */
    irrep_sg_irrep_free(A2_p1);
    irrep_space_group_free(G1);

    /* Symmetry-adapted basis builder test on a 2×2 square / p4mm cluster
     * (4 sites, Hilbert-space dim 16). The Γ-invariant (A₁) sector
     * dim equals the number of T-orbits = (16 + 3·4)/4 = 7 by Burnside
     * on the order-4 translation subgroup. Builder returns an orthonormal
     * basis of this sector. */
    {
        irrep_lattice_t     *sq2   = irrep_lattice_build(IRREP_LATTICE_SQUARE, 2, 2);
        irrep_space_group_t *G4 = irrep_space_group_build(sq2, IRREP_WALLPAPER_P4MM);
        IRREP_ASSERT(G4 != NULL);

        long long D = 1LL << 4;                 /* local_dim = 2, N = 4 */
        double _Complex *basis = malloc(sizeof(double _Complex) * (size_t)D * (size_t)D);
        IRREP_ASSERT(basis != NULL);

        irrep_sg_irrep_t *A1_sq = irrep_sg_trivial(G4);
        int n_basis = irrep_sg_adapted_basis(G4, A1_sq, 4, 2, basis, (int)D);
        IRREP_ASSERT(n_basis > 0);
        IRREP_ASSERT(n_basis <= (int)D);

        /* Orthonormality: ⟨b_i|b_j⟩ = δ_ij */
        for (int i = 0; i < n_basis; ++i) {
            for (int j = i; j < n_basis; ++j) {
                double _Complex overlap = 0.0;
                for (long long t = 0; t < D; ++t) {
                    overlap += conj(basis[(size_t)i * D + t]) * basis[(size_t)j * D + t];
                }
                double expected = (i == j) ? 1.0 : 0.0;
                IRREP_ASSERT_NEAR(creal(overlap), expected, 1e-10);
                IRREP_ASSERT_NEAR(cimag(overlap), 0.0,      1e-10);
            }
        }

        /* Error paths */
        IRREP_ASSERT(irrep_sg_adapted_basis(NULL, A1_sq, 4, 2, basis, (int)D) == -1);
        IRREP_ASSERT(irrep_sg_adapted_basis(G4, NULL,    4, 2, basis, (int)D) == -1);
        IRREP_ASSERT(irrep_sg_adapted_basis(G4, A1_sq,   4, 1, basis, (int)D) == -1);
        /* Wrong num_sites */
        IRREP_ASSERT(irrep_sg_adapted_basis(G4, A1_sq,   8, 2, basis, (int)D) == -1);

        irrep_sg_irrep_free(A1_sq);
        irrep_space_group_free(G4);
        irrep_lattice_free(sq2);
        free(basis);
    }

    /* Orbit enumeration: check that orbit[g·n + s] = sigma[g^{-1}·s]. */
    double *sigma = malloc(sizeof(double) * n);
    for (int s = 0; s < n; ++s) sigma[s] = (double)(s + 1);
    double *orbit = malloc(sizeof(double) * (size_t)order * n);
    irrep_sg_enumerate_orbit(G, sigma, orbit);

    int *inv = malloc(sizeof(int) * n);
    for (int g = 0; g < order; g += 13) {
        irrep_space_group_permutation_inverse(G, g, inv);
        int ok = 1;
        for (int s = 0; s < n && ok; ++s) {
            if (orbit[g*n + s] != sigma[inv[s]]) ok = 0;
        }
        IRREP_ASSERT(ok);
    }
    free(inv);
    free(orbit);
    free(sigma);

    free(psi);
    irrep_sg_irrep_free(A1);
    irrep_space_group_free(G);
    irrep_lattice_free(sq);

    /* Error paths */
    IRREP_ASSERT(irrep_sg_irrep_new(NULL, NULL, 1) == NULL);
    IRREP_ASSERT(irrep_sg_trivial(NULL)            == NULL);
    irrep_sg_irrep_free(NULL);                                    /* no-op */

    return IRREP_TEST_END();
}
