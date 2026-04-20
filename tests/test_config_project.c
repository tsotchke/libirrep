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
#include <stdint.h>
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

    /* -------------------------------------------------------------------- *
     * Bloch (non-Γ) momentum projection and basis tests.                   *
     *                                                                      *
     * On a 3×3 square cluster under p1 (translations only, 9 elements),   *
     * sum of sector dimensions across all 9 k-points must equal 2^9 = 512 *
     * (full Hilbert-space dimension). Basis vectors across different k    *
     * sectors must be orthogonal. At Γ (k=0), Bloch reduces to A₁.        *
     * -------------------------------------------------------------------- */
    {
        irrep_lattice_t     *L3  = irrep_lattice_build(IRREP_LATTICE_SQUARE, 3, 3);
        irrep_space_group_t *G3  = irrep_space_group_build(L3, IRREP_WALLPAPER_P1);
        IRREP_ASSERT(G3 != NULL);
        IRREP_ASSERT(irrep_space_group_lattice(G3) == L3);

        int Nsites = 9;
        long long D = 1LL << Nsites;

        /* Γ-sector via Bloch must match A₁ via irrep_sg_project_A1 on any
         * amplitude array. */
        int order3 = irrep_space_group_order(G3);
        double _Complex *psi9 = malloc(sizeof(double _Complex) * order3);
        uint64_t rng = 0xC0FFEEULL;
        for (int g = 0; g < order3; ++g) {
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            double re = ((rng >> 16) & 0xFFFF) / 65536.0 - 0.5;
            double im = ((rng >> 32) & 0xFFFF) / 65536.0 - 0.5;
            psi9[g] = re + I * im;
        }
        double _Complex a1 = irrep_sg_project_A1(G3, psi9);
        double _Complex b0 = irrep_sg_bloch_amplitude(G3, 0, 0, psi9);
        IRREP_ASSERT_NEAR(creal(a1), creal(b0), 1e-12);
        IRREP_ASSERT_NEAR(cimag(a1), cimag(b0), 1e-12);

        /* Sum of Bloch projections over all k equals ψ at g=0 (Fourier
         * inversion on the translation group). */
        double _Complex acc = 0.0;
        for (int kx = 0; kx < 3; ++kx)
            for (int ky = 0; ky < 3; ++ky)
                acc += irrep_sg_bloch_amplitude(G3, kx, ky, psi9);
        IRREP_ASSERT_NEAR(creal(acc), creal(psi9[0]), 1e-12);
        IRREP_ASSERT_NEAR(cimag(acc), cimag(psi9[0]), 1e-12);
        free(psi9);

        /* Sector-dimension sum must equal full Hilbert dim. */
        double _Complex *basis = malloc(sizeof(double _Complex) * (size_t)D * (size_t)D);
        IRREP_ASSERT(basis != NULL);
        int total_dim = 0;
        int sector_dim[9] = {0};
        for (int kx = 0; kx < 3; ++kx) {
            for (int ky = 0; ky < 3; ++ky) {
                int nb = irrep_sg_bloch_basis(G3, kx, ky, Nsites, 2, basis, (int)D);
                IRREP_ASSERT(nb > 0);
                sector_dim[ky * 3 + kx] = nb;
                total_dim += nb;

                /* Orthonormality within the sector. */
                for (int i = 0; i < nb && i < 10; ++i) {
                    for (int j = i; j < nb && j < 10; ++j) {
                        double _Complex ov = 0.0;
                        for (long long t = 0; t < D; ++t)
                            ov += conj(basis[(size_t)i * D + t]) * basis[(size_t)j * D + t];
                        double ex = (i == j) ? 1.0 : 0.0;
                        IRREP_ASSERT_NEAR(creal(ov), ex,  1e-10);
                        IRREP_ASSERT_NEAR(cimag(ov), 0.0, 1e-10);
                    }
                }
            }
        }
        IRREP_ASSERT(total_dim == (int)D);

        /* Cross-sector orthogonality: a basis vector from k=(1,0) should be
         * orthogonal to every basis vector from k=(0,0). Build k=(0,0) first,
         * keep a copy of its basis, then k=(1,0), test. */
        double _Complex *gamma = malloc(sizeof(double _Complex) * (size_t)D * (size_t)D);
        double _Complex *mpoint = malloc(sizeof(double _Complex) * (size_t)D * (size_t)D);
        int nG = irrep_sg_bloch_basis(G3, 0, 0, Nsites, 2, gamma, (int)D);
        int nM = irrep_sg_bloch_basis(G3, 1, 0, Nsites, 2, mpoint, (int)D);
        IRREP_ASSERT(nG > 0 && nM > 0);
        for (int i = 0; i < nG && i < 5; ++i) {
            for (int j = 0; j < nM && j < 5; ++j) {
                double _Complex ov = 0.0;
                for (long long t = 0; t < D; ++t)
                    ov += conj(gamma[(size_t)i * D + t]) * mpoint[(size_t)j * D + t];
                IRREP_ASSERT_NEAR(creal(ov), 0.0, 1e-10);
                IRREP_ASSERT_NEAR(cimag(ov), 0.0, 1e-10);
            }
        }
        free(gamma);
        free(mpoint);
        free(basis);

        /* Canonicalisation: kx = -1 must match kx = Lx - 1; kx = 2·Lx + 1 must
         * match kx = 1. Same for ky. Test on a small random amplitude. */
        double _Complex *psi_bound = malloc(sizeof(double _Complex) * order3);
        for (int g = 0; g < order3; ++g) psi_bound[g] = (double)(g % 5) - 2.0;
        double _Complex k_pos = irrep_sg_bloch_amplitude(G3,  1,  2, psi_bound);
        double _Complex k_neg = irrep_sg_bloch_amplitude(G3,  1, -1, psi_bound);   /* ky = -1 ≡ 2 */
        double _Complex k_far = irrep_sg_bloch_amplitude(G3, -2,  2, psi_bound);   /* kx = -2 ≡ 1 */
        IRREP_ASSERT_NEAR(creal(k_pos), creal(k_neg), 1e-12);
        IRREP_ASSERT_NEAR(cimag(k_pos), cimag(k_neg), 1e-12);
        IRREP_ASSERT_NEAR(creal(k_pos), creal(k_far), 1e-12);
        IRREP_ASSERT_NEAR(cimag(k_pos), cimag(k_far), 1e-12);
        free(psi_bound);

        /* Error paths */
        IRREP_ASSERT(irrep_sg_bloch_basis(NULL, 0, 0, Nsites, 2, NULL, 1) == -1);
        IRREP_ASSERT(irrep_sg_bloch_basis(G3, 0, 0, 8, 2, (double _Complex*)&total_dim, 1) == -1);
        /* NULL inputs: NaN result (IEEE-754 way to distinguish error from a
         * legitimate zero projection). */
        double _Complex bad = irrep_sg_bloch_amplitude(NULL, 0, 0, NULL);
        IRREP_ASSERT(isnan(creal(bad)) && isnan(cimag(bad)));
        bad = irrep_sg_project_A1(NULL, NULL);
        IRREP_ASSERT(isnan(creal(bad)) && isnan(cimag(bad)));
        bad = irrep_sg_project_amplitude(NULL, NULL);
        IRREP_ASSERT(isnan(creal(bad)) && isnan(cimag(bad)));

        irrep_space_group_free(G3);
        irrep_lattice_free(L3);
    }

    /* Error paths */
    IRREP_ASSERT(irrep_sg_irrep_new(NULL, NULL, 1) == NULL);
    IRREP_ASSERT(irrep_sg_trivial(NULL)            == NULL);
    irrep_sg_irrep_free(NULL);                                    /* no-op */

    return IRREP_TEST_END();
}
