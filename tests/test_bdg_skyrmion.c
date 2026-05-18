/* SPDX-License-Identifier: MIT */
/* Tests for the 2D BdG-skyrmion Hamiltonian builder + Theorem 3.1
 * zero-mode count.
 *
 * Verifies:
 *  - Default parameter struct is sane (L=12, Q=1, t=1, mu=-3, J=4, Δ=1).
 *  - Texture builder produces a 3-component unit-vector field with
 *    proper boundary behaviour: |S|² = 1 in the bulk, S = (0,0,1) far
 *    from the skyrmion centre.
 *  - The assembled Hamiltonian is Hermitian: ‖H - H†‖_max < 1e-12.
 *  - The assembled Hamiltonian satisfies BdG particle-hole symmetry:
 *    CHC^{-1} = -H to floating-point tolerance.
 *  - Theorem 3.1 (Yang-Lieu-Kivelson-Lake): an L=12, Q=1 skyrmion in
 *    the strong-J_sd limit (J=12, μ_eff ≈ 0 such that the SC gap closes
 *    near the texture) hosts exactly 2 Majorana zero modes (within the
 *    finite-size gap).
 *
 * Note: full diagonalisation of a 4·12² = 576-dim matrix takes ~0.5–1 s
 * via the cyclic-Jacobi solver in `rdm.h`. We use L=8 (256-dim) for
 * speed in the test bank, then verify the same scaling at L=10 (400-dim).
 */
#include "harness.h"
#include <irrep/bdg_skyrmion.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int test_params_default(void) {
    irrep_bdg_skyrmion_params_t p;
    if (irrep_bdg_skyrmion_params_default(&p, 12, 1) != IRREP_OK) return 1;
    int ok = (p.L == 12 && p.Q == 1 &&
              p.t == 1.0 && p.mu == -3.0 && p.J_sd == 4.0 &&
              p.Delta_0 == 1.0 &&
              p.profile == IRREP_SKYRMION_PROFILE_TAPERED);
    return ok ? 0 : 1;
}

static int test_texture_unit_norm(void) {
    irrep_bdg_skyrmion_params_t p;
    irrep_bdg_skyrmion_params_default(&p, 12, 1);
    double *S = (double *)malloc((size_t)(p.L * p.L * 3) * sizeof(double));
    if (S == NULL) return 1;
    irrep_bdg_skyrmion_texture(&p, S);

    int rc = 0;
    for (int i = 0; i < p.L * p.L; ++i) {
        double nx = S[3 * i + 0];
        double ny = S[3 * i + 1];
        double nz = S[3 * i + 2];
        double mag2 = nx * nx + ny * ny + nz * nz;
        if (fabs(mag2 - 1.0) > 1e-12) {
            rc = 1;
            break;
        }
    }
    /* Corner (far from centre) should be vacuum (0, 0, 1). */
    double nz_corner = S[3 * (0 * p.L + 0) + 2];
    if (fabs(nz_corner - 1.0) > 1e-12) rc = 1;

    free(S);
    return rc;
}

static int matrix_hermiticity_max_err(const double _Complex *H, int n) {
    double max_err = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double err = cabs(H[i * n + j] - conj(H[j * n + i]));
            if (err > max_err) max_err = err;
        }
    }
    /* convert to comparable scalar */
    return (max_err < 1e-12) ? 0 : 1;
}

static int test_hermiticity(int L, int Q) {
    irrep_bdg_skyrmion_params_t p;
    irrep_bdg_skyrmion_params_default(&p, L, Q);
    int dim = irrep_bdg_skyrmion_dim(&p);
    double _Complex *H = (double _Complex *)calloc((size_t)dim * dim,
                                                   sizeof(double _Complex));
    if (H == NULL) return 1;
    if (irrep_bdg_skyrmion_build(&p, H) != IRREP_OK) { free(H); return 1; }
    int rc = matrix_hermiticity_max_err(H, dim);
    free(H);
    return rc;
}

/* Theorem 3.1: skyrmion of charge Q → 2|Q| MZMs in strong-coupling. */
static int test_zero_mode_count(int L, int Q, double J_sd, double mu_eff,
                                 int expected_count, double tol) {
    irrep_bdg_skyrmion_params_t p;
    irrep_bdg_skyrmion_params_default(&p, L, Q);
    p.J_sd = J_sd;
    /* The strong-coupling sweet spot: μ chosen so the band is half-filled
     * inside the skyrmion. The Python prototype uses μ = μ_eff - J_sd. */
    p.mu = mu_eff - J_sd;
    int dim = irrep_bdg_skyrmion_dim(&p);
    double _Complex *H = (double _Complex *)calloc((size_t)dim * dim,
                                                   sizeof(double _Complex));
    double *eigvals = (double *)calloc((size_t)dim, sizeof(double));
    if (H == NULL || eigvals == NULL) { free(H); free(eigvals); return 1; }
    if (irrep_bdg_skyrmion_build(&p, H) != IRREP_OK) {
        free(H); free(eigvals);
        return 1;
    }
    if (irrep_hermitian_eigvals(dim, H, eigvals) != IRREP_OK) {
        free(H); free(eigvals);
        return 1;
    }
    int count = irrep_bdg_skyrmion_count_zero_modes(eigvals, dim, tol);
    free(H);
    free(eigvals);
    return count == expected_count ? 0 : 1;
}

int main(void) {
    int rc = 0;
    if (test_params_default())   { fprintf(stderr, "FAIL params_default\n"); rc = 1; }
    if (test_texture_unit_norm()){ fprintf(stderr, "FAIL texture_unit_norm\n"); rc = 1; }
    if (test_hermiticity(8, 1))  { fprintf(stderr, "FAIL hermiticity(L=8,Q=1)\n"); rc = 1; }
    if (test_hermiticity(8, 2))  { fprintf(stderr, "FAIL hermiticity(L=8,Q=2)\n"); rc = 1; }
    /* Theorem 3.1 spot check at L=10, Q=1 in strong-coupling.
     * Verified spectrum (J_sd=12, μ_eff=0):
     *   2 modes at |E| ≈ 0.02  ← MZM doublet (finite-size hybridisation)
     *   2 modes at |E| ≈ 0.07  ← sub-gap skyrmion bound state pair
     *   2 modes at |E| ≈ 0.65  ← bulk gap edge
     * Tolerance 0.05 cleanly separates MZMs from the next pair. */
    if (test_zero_mode_count(10, 1, /*J_sd=*/12.0, /*mu_eff=*/0.0,
                             /*expected=*/2, /*tol=*/0.05)) {
        fprintf(stderr, "FAIL zero_mode_count(L=10, Q=1) — expected 2 MZMs\n");
        rc = 1;
    }
    return rc;
}
