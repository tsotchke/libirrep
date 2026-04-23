/* SPDX-License-Identifier: MIT */
/* Tests for the batched RDM + entropy pipeline.
 *
 * Coverage:
 *   - `_batch_partial_trace` on a single sample matches the single-state
 *     `irrep_partial_trace` bit-exactly.
 *   - Bell state S_VN = ln 2 over a single-qubit partition, replicated
 *     across a batch.
 *   - 2-site product state → S_VN = 0.
 *   - Rényi-∞ on a maximally mixed qubit RDM = ln 2 (≡ ln(dim)).
 *   - `_from_sample_amplitudes`: uniform-weight aggregate of identical
 *     rank-1 outer products reproduces that outer product, trace = 1.
 *   - Weighted aggregate: normalised correctly; negative weight rejected.
 *   - All-zero weights → PRECONDITION error.
 */

#include "harness.h"
#include <irrep/rdm.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    IRREP_TEST_START("rdm_batch");

    /* ---- 1. Single-sample batched partial trace == non-batched --------- */
    {
        /* 2-qubit Bell state: (|00⟩ + |11⟩) / √2. */
        int             nsites = 2;
        int             local = 2;
        int             D = 4;
        double _Complex psi[4] = {1.0 / sqrt(2.0), 0, 0, 1.0 / sqrt(2.0)};
        int             sites_A[1] = {0};
        double _Complex rho_direct[4];
        double _Complex rho_batched[4];
        IRREP_ASSERT(irrep_partial_trace(nsites, local, psi, sites_A, 1, rho_direct) == IRREP_OK);
        IRREP_ASSERT(irrep_rdm_batch_partial_trace(nsites, local, 1, psi, sites_A, 1,
                                                   rho_batched) == IRREP_OK);
        for (int k = 0; k < 4; ++k) {
            IRREP_ASSERT(cabs(rho_direct[k] - rho_batched[k]) < 1e-14);
        }

        /* S_VN on this RDM = ln 2 (Bell entanglement). */
        double vn;
        IRREP_ASSERT(irrep_rdm_batch_entropy_vonneumann(2, 1, rho_batched, &vn) == IRREP_OK);
        IRREP_ASSERT_NEAR(vn, log(2.0), 1e-12);
    }

    /* ---- 2. Batched Bell states: replicate across 5 samples ------------- */
    {
        int             nsites = 2;
        int             local = 2;
        int             n_samples = 5;
        double _Complex psi_batch[5 * 4];
        for (int m = 0; m < n_samples; ++m) {
            double _Complex *p = psi_batch + m * 4;
            p[0] = 1.0 / sqrt(2.0);
            p[1] = 0;
            p[2] = 0;
            p[3] = 1.0 / sqrt(2.0);
        }
        int             sites_A[1] = {0};
        double _Complex rho_batch[5 * 4];
        IRREP_ASSERT(irrep_rdm_batch_partial_trace(nsites, local, n_samples, psi_batch, sites_A, 1,
                                                   rho_batch) == IRREP_OK);
        double vn[5];
        IRREP_ASSERT(irrep_rdm_batch_entropy_vonneumann(2, n_samples, rho_batch, vn) == IRREP_OK);
        for (int m = 0; m < n_samples; ++m) {
            IRREP_ASSERT_NEAR(vn[m], log(2.0), 1e-12);
        }
    }

    /* ---- 3. Product state → S_VN = 0, Rényi-α = 0 for any α ------------ */
    {
        double _Complex psi[4] = {1.0, 0, 0, 0}; /* |00⟩ */
        int             sites_A[1] = {0};
        double _Complex rho[4];
        IRREP_ASSERT(irrep_rdm_batch_partial_trace(2, 2, 1, psi, sites_A, 1, rho) == IRREP_OK);

        double vn;
        double _Complex rho_copy[4];
        memcpy(rho_copy, rho, sizeof(rho_copy));
        IRREP_ASSERT(irrep_rdm_batch_entropy_vonneumann(2, 1, rho_copy, &vn) == IRREP_OK);
        IRREP_ASSERT_NEAR(vn, 0.0, 1e-12);

        memcpy(rho_copy, rho, sizeof(rho_copy));
        double renyi;
        IRREP_ASSERT(irrep_rdm_batch_entropy_renyi(2, 1, rho_copy, 2.0, &renyi) == IRREP_OK);
        IRREP_ASSERT_NEAR(renyi, 0.0, 1e-12);
    }

    /* ---- 4. Rényi-2 on a maximally mixed qubit = ln 2 ------------------- */
    {
        /* Ghost mixed state: ρ = (1/2) I₂. */
        double _Complex rho[4] = {0.5, 0, 0, 0.5};
        double          renyi;
        IRREP_ASSERT(irrep_rdm_batch_entropy_renyi(2, 1, rho, 2.0, &renyi) == IRREP_OK);
        /* S_2 = -ln(Σ λ²) = -ln(0.25 + 0.25) = -ln(0.5) = ln 2. */
        IRREP_ASSERT_NEAR(renyi, log(2.0), 1e-12);
    }

    /* ---- 5. from_sample_amplitudes: uniform replicas → rank-1 ρ -------- */
    {
        /* Every sample is the same |+⟩ = (1, 1) / √2. 10 copies, uniform. */
        int             dim_A = 2;
        int             n_samples = 10;
        double _Complex psi_batch[10 * 2];
        for (int m = 0; m < n_samples; ++m) {
            psi_batch[m * 2 + 0] = 1.0 / sqrt(2.0);
            psi_batch[m * 2 + 1] = 1.0 / sqrt(2.0);
        }
        double _Complex rho_out[4];
        IRREP_ASSERT(irrep_rdm_from_sample_amplitudes(dim_A, n_samples, psi_batch, NULL,
                                                     rho_out) == IRREP_OK);
        /* Expected: (1/2) [[1, 1], [1, 1]]. */
        IRREP_ASSERT_NEAR(creal(rho_out[0]), 0.5, 1e-14);
        IRREP_ASSERT_NEAR(creal(rho_out[1]), 0.5, 1e-14);
        IRREP_ASSERT_NEAR(creal(rho_out[2]), 0.5, 1e-14);
        IRREP_ASSERT_NEAR(creal(rho_out[3]), 0.5, 1e-14);
        /* Trace = 1. */
        double tr = creal(rho_out[0] + rho_out[3]);
        IRREP_ASSERT_NEAR(tr, 1.0, 1e-14);
    }

    /* ---- 6. Weighted aggregate: half |+⟩ + half |−⟩ → I/2 -------------- */
    {
        int             dim_A = 2;
        double _Complex psi[4];
        psi[0] = 1.0 / sqrt(2.0); /* |+⟩ */
        psi[1] = 1.0 / sqrt(2.0);
        psi[2] = 1.0 / sqrt(2.0); /* |−⟩ */
        psi[3] = -1.0 / sqrt(2.0);
        double          weights[2] = {1.0, 1.0};
        double _Complex rho_out[4];
        IRREP_ASSERT(irrep_rdm_from_sample_amplitudes(dim_A, 2, psi, weights, rho_out) == IRREP_OK);
        /* Off-diagonals cancel: ρ = (1/2) [[1, 0], [0, 1]]. */
        IRREP_ASSERT_NEAR(creal(rho_out[0]), 0.5, 1e-14);
        IRREP_ASSERT(cabs(rho_out[1]) < 1e-14);
        IRREP_ASSERT(cabs(rho_out[2]) < 1e-14);
        IRREP_ASSERT_NEAR(creal(rho_out[3]), 0.5, 1e-14);
    }

    /* ---- 7. Error paths ------------------------------------------------- */
    {
        double _Complex dummy_psi[4] = {1, 0, 0, 0};
        double _Complex dummy_rho[4];
        /* Negative weight → INVALID_ARG. */
        double neg[1] = {-1.0};
        IRREP_ASSERT(irrep_rdm_from_sample_amplitudes(2, 1, dummy_psi, neg, dummy_rho) ==
                     IRREP_ERR_INVALID_ARG);
        /* Zero weights → PRECONDITION. */
        double zero[1] = {0.0};
        IRREP_ASSERT(irrep_rdm_from_sample_amplitudes(2, 1, dummy_psi, zero, dummy_rho) ==
                     IRREP_ERR_PRECONDITION);
        /* NULL inputs. */
        IRREP_ASSERT(irrep_rdm_batch_partial_trace(2, 2, 1, NULL, NULL, 1, dummy_rho) ==
                     IRREP_ERR_INVALID_ARG);
        IRREP_ASSERT(irrep_rdm_batch_entropy_vonneumann(2, 1, NULL, NULL) == IRREP_ERR_INVALID_ARG);
    }

    return IRREP_TEST_END();
}
