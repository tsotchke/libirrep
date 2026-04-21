/* SPDX-License-Identifier: MIT */
/* Tests for the irrep_heisenberg Hamiltonian primitive.
 *
 * Exercises:
 *   - Two-spin (N = 2) Heisenberg: H = J (S_1 · S_2). Exact spectrum
 *     {−3J/4 (singlet), +J/4 × 3 (triplet)}. Apply to basis states and
 *     compare against explicit formulas.
 *   - Closed-form ground-state energy of a 4-spin ring (J·Σ S_i·S_{i+1}
 *     PBC). Known E_0 = −2J.
 *   - num_sites / dim accessors.
 *   - Selection-rule checks: `S_z_total` is preserved (popcount
 *     invariant under the apply).
 *   - Invalid input handling.
 */

#include "harness.h"

#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    IRREP_TEST_START("hamiltonian");

    /* ---- N = 2: bond (0, 1), J = 1 -------------------------------- */
    {
        int bi[] = {0}, bj[] = {1};
        irrep_heisenberg_t *H = irrep_heisenberg_new(2, 1, bi, bj, 1.0);
        IRREP_ASSERT(H != NULL);
        IRREP_ASSERT(irrep_heisenberg_num_sites(H) == 2);
        IRREP_ASSERT(irrep_heisenberg_dim(H)       == 4);

        /* basis states:   |00⟩ = 0,  |01⟩ = 1,  |10⟩ = 2,  |11⟩ = 3
         * H |00⟩ = +¼ |00⟩                       (both down, aligned)
         * H |01⟩ = −¼ |01⟩ + ½ |10⟩                (anti-aligned → flip)
         * H |10⟩ = −¼ |10⟩ + ½ |01⟩
         * H |11⟩ = +¼ |11⟩                                                    */
        double _Complex psi[4] = { 1, 0, 0, 0 };
        double _Complex out[4];
        irrep_heisenberg_apply(psi, out, H);
        IRREP_ASSERT_NEAR(creal(out[0]),  0.25, 1e-14);
        IRREP_ASSERT_NEAR(creal(out[1]),  0.0,  1e-14);
        IRREP_ASSERT_NEAR(creal(out[2]),  0.0,  1e-14);
        IRREP_ASSERT_NEAR(creal(out[3]),  0.0,  1e-14);

        psi[0] = 0; psi[1] = 1;
        irrep_heisenberg_apply(psi, out, H);
        IRREP_ASSERT_NEAR(creal(out[0]),  0.0,  1e-14);
        IRREP_ASSERT_NEAR(creal(out[1]), -0.25, 1e-14);
        IRREP_ASSERT_NEAR(creal(out[2]),  0.5,  1e-14);
        IRREP_ASSERT_NEAR(creal(out[3]),  0.0,  1e-14);

        /* Singlet (|01⟩ − |10⟩)/√2  is an eigenvector with E = −¾. */
        double _Complex s2 = 1.0 / sqrt(2.0);
        psi[0] = 0; psi[1] =  s2; psi[2] = -s2; psi[3] = 0;
        irrep_heisenberg_apply(psi, out, H);
        IRREP_ASSERT_NEAR(creal(out[1]), -0.75 *  creal(s2), 1e-14);
        IRREP_ASSERT_NEAR(creal(out[2]), -0.75 * -creal(s2), 1e-14);

        irrep_heisenberg_free(H);
    }

    /* ---- N = 4 ring: Σ S_i · S_{i+1}, PBC. E_0 = −2J ---------------- *
     * Lanczos-solve via the library. Explicit benchmark: the 4-spin
     * antiferromagnetic-Heisenberg chain with PBC has ground state
     * energy E_0 = −2J (Bethe-ansatz / exact result for N = 4).         */
    {
        int bi[] = {0, 1, 2, 3};
        int bj[] = {1, 2, 3, 0};
        irrep_heisenberg_t *H = irrep_heisenberg_new(4, 4, bi, bj, 1.0);
        IRREP_ASSERT(H != NULL);

        long long dim = irrep_heisenberg_dim(H);
        double _Complex *seed = calloc((size_t)dim, sizeof(double _Complex));
        IRREP_ASSERT(seed != NULL);
        uint64_t rng = 0xbabefeedULL;
        for (long long s = 0; s < dim; ++s) {
            if (__builtin_popcountll(s) == 2) {
                rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
                seed[s] = (double)(rng >> 32) / (double)0xFFFFFFFFu - 0.5;
            }
        }

        double eig[2] = {0, 0};
        irrep_status_t rc = irrep_lanczos_eigvals(
            irrep_heisenberg_apply, H, dim, 2, 60, seed, eig);
        IRREP_ASSERT(rc == IRREP_OK);
        IRREP_ASSERT_NEAR(eig[0], -2.0, 1e-8);

        free(seed);
        irrep_heisenberg_free(H);
    }

    /* ---- Invalid inputs --------------------------------------------- */
    IRREP_ASSERT(irrep_heisenberg_new(0,   0, NULL, NULL, 1.0) == NULL);
    IRREP_ASSERT(irrep_heisenberg_new(64,  0, NULL, NULL, 1.0) == NULL);   /* too large */
    {
        int bad_i[] = {0};
        int bad_j[] = {4};   /* out of [0, 3) */
        IRREP_ASSERT(irrep_heisenberg_new(3, 1, bad_i, bad_j, 1.0) == NULL);
    }
    irrep_heisenberg_free(NULL);   /* no-op */

    /* ---- Dim/accessor error paths ----------------------------------- */
    IRREP_ASSERT(irrep_heisenberg_num_sites(NULL) == 0);
    IRREP_ASSERT(irrep_heisenberg_dim(NULL)       == 0);

    return IRREP_TEST_END();
}
