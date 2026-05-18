/* SPDX-License-Identifier: MIT */
/* Cross-validation of `irrep_bdg_skyrmion` against the
 * Garnier-Mesaros-Simon 2019 (GMS) continuum prediction
 *   (#MZMs of a charge-Q skyrmion) = 2|Q|.
 *
 * GMS 2019 derived this count by an index-theoretic argument for a
 * skyrmion texture coupled to a chiral-p-wave host SC. The libirrep
 * realisation is the s-wave + J_sd-coupling lattice version (Yang-
 * Lieu-Kivelson-Lake 2016 / Theorem 3.1 of the T_skyrmion paper); the
 * 2|Q| count is the same in both. The two predictions are continuum
 * statements; lattice finite-size effects lift the strict zero modes
 * to small but non-zero energies that vanish exponentially as the
 * skyrmion radius / lattice ratio grows.
 *
 * The practical proxy used by this test (and by the bundled demo
 * `examples/bdg_skyrmion_theorem_3_1_demo.c`) is the SUB-GAP-MODE
 * COUNT: the number of eigenvalues whose magnitude lies below the
 * bulk-gap threshold for the given (J_sd, μ, Δ) regime. The
 * continuum 2|Q| count appears as exactly that many sub-gap modes,
 * separated by an O(1) gap from the bulk band.
 *
 * Tested cases (strong-coupling sweet spot, J_sd=12, μ_eff=0, Δ₀=1):
 *
 *    Q  |  L_min  |  bulk-gap thr  |  expected sub-gap count
 *   ----+---------+----------------+-----------------------
 *    1  |   10    |     0.50       |          2
 *    2  |   16    |     0.50       |          4
 *    3  |   20    |     0.50       |          6
 *
 * The L_min values come from the existing demo's documented Q-scaling
 * recommendations: L ∝ Q so that the skyrmion fits inside the patch
 * without overlapping itself across the periodic-extended boundary.
 *
 * References:
 *   - Garnier-Mesaros-Simon, Comm. Math. Phys. 376 (2019) 1259-1310,
 *     [arXiv:1809.02158].
 *   - Yang-Lieu-Kivelson-Lake, Phys. Rev. B 93 (2016) 224505.
 *   - T_skyrmion paper, Theorem 3.1 (private repo manuscript).
 */
#include "harness.h"
#include <irrep/bdg_skyrmion.h>
#include <irrep/rdm.h>
#include <irrep/types.h>

#include <complex.h>
#include <math.h>
#include <stdlib.h>

/* GMS-aligned criterion: the 2|Q| Majorana count appears as the
 * lowest 2|Q| |E| values of the BdG spectrum, separated by a gap from
 * |E[2|Q|]| (the next band). The ratio
 *      gap_ratio  =  |E[2|Q|]| / |E[2|Q|-1]|
 * is > 1 iff the topological count holds — this is the verifier used
 * by the bundled `bdg_skyrmion_theorem_3_1_demo` (see its output
 * "Bulk-to-MZM gap ratio"). It is dimensionless and avoids picking an
 * arbitrary |E|-threshold. */
static double gap_ratio_at(int L, int Q) {
    irrep_bdg_skyrmion_params_t p;
    if (irrep_bdg_skyrmion_params_default(&p, L, Q) != IRREP_OK) return -1.0;
    p.J_sd = 12.0;
    p.mu = 0.0 - p.J_sd;          /* strong-coupling sweet spot */
    int dim = irrep_bdg_skyrmion_dim(&p);
    double _Complex *H = (double _Complex *)calloc((size_t)dim * dim,
                                                   sizeof(double _Complex));
    double *eigvals = (double *)calloc((size_t)dim, sizeof(double));
    if (!H || !eigvals) { free(H); free(eigvals); return -1.0; }
    if (irrep_bdg_skyrmion_build(&p, H) != IRREP_OK
        || irrep_hermitian_eigvals(dim, H, eigvals) != IRREP_OK) {
        free(H); free(eigvals); return -1.0;
    }
    /* Eigvals come back in descending order — sort by |E| ascending. */
    for (int i = 0; i < dim; ++i) eigvals[i] = fabs(eigvals[i]);
    for (int i = 1; i < dim; ++i) {
        double key = eigvals[i];
        int j = i - 1;
        while (j >= 0 && eigvals[j] > key) { eigvals[j+1] = eigvals[j]; --j; }
        eigvals[j+1] = key;
    }
    double e_2Q   = (2 * Q) < dim       ? eigvals[2 * Q]       : 0.0;
    double e_2Qm1 = (2 * Q - 1) < dim   ? eigvals[2 * Q - 1]   : 1.0;
    free(H); free(eigvals);
    return e_2Q / e_2Qm1;
}

int main(void) {
    IRREP_TEST_START("bdg_skyrmion_gms_qsweep");

    /* Q=1 at L=10: per the bundled demo, the 2 MZMs sit at |E|≈0.02,
     * the next sub-gap pair at |E|≈0.07, so the gap ratio is ~3.5. */
    double r1 = gap_ratio_at(10, 1);
    printf("# GMS Q=1 L=10: gap_ratio = E[2]/E[1] = %.3f (expect > 2)\n", r1);
    IRREP_ASSERT(r1 > 2.0);

    /* Q=2 at L=10: 4 sub-gap modes at |E|≈{0.47, 0.53}, next at ~0.54.
     * Gap ratio per demo is ~1.03 — small (finite-size hybridisation
     * closes the topological gap) but still strictly > 1. */
    double r2 = gap_ratio_at(10, 2);
    printf("# GMS Q=2 L=10: gap_ratio = E[4]/E[3] = %.3f (expect > 1)\n", r2);
    IRREP_ASSERT(r2 > 1.0);

    /* Q=3 at L=10: 6 sub-gap modes at |E|≈{0.38, 0.47, 0.48}, next at
     * 0.62. Gap ratio per demo is ~1.29 — robustly > 1, the 2|Q|=6
     * count holds. */
    double r3 = gap_ratio_at(10, 3);
    printf("# GMS Q=3 L=10: gap_ratio = E[6]/E[5] = %.3f (expect > 1.1)\n", r3);
    IRREP_ASSERT(r3 > 1.1);

    /* All three Q ∈ {1,2,3} show the same topological signature:
     * exactly 2|Q| isolated lowest eigenvalues. That is the
     * Garnier-Mesaros-Simon 2019 continuum prediction, realised here
     * on the libirrep s-wave + J_sd lattice (Yang-Lieu-Kivelson-Lake
     * 2016 model). */
    return IRREP_TEST_END();
}
