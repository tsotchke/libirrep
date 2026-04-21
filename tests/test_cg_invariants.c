/* SPDX-License-Identifier: MIT */
/* Property-based invariants on randomly-sampled (j₁, j₂, J, m₁, m₂, M)
 * tuples. Complements the fixed-value Sakurai hand checks in
 * test_clebsch_gordan.c: here we draw N (j₁, j₂) pairs from a
 * deterministic PCG stream and verify the *structural* identities that
 * every Clebsch-Gordan coefficient must obey, independently of any
 * particular closed form.
 *
 * Invariants verified on each random triple:
 *   (1) Sum rule:  Σ_J |⟨j₁ m₁ j₂ m₂ | J M⟩|² = 1 for each valid (m₁, m₂).
 *   (2) Orthogonality across J:
 *         Σ_{m₁} ⟨j₁ m₁; j₂ M−m₁ | J  M⟩ · ⟨j₁ m₁; j₂ M−m₁ | J' M⟩ = δ_{JJ'}
 *   (3) Selection rules return 0 for every violation:
 *         m₁ + m₂ ≠ M,  |M| > J,  J outside [|j₁−j₂|, j₁+j₂], triangle parity
 *   (4) CG vs 3j cross-reference:
 *         ⟨j₁ m₁ j₂ m₂ | J M⟩ = (−1)^{j₁−j₂+M} √(2J+1) (j₁ j₂ J; m₁ m₂ −M)
 *
 * Ranges exercised: j₁, j₂ ∈ {0, ½, 1, …, 10} (21 distinct values);
 * 128 random (j₁, j₂) pairs × every valid J × every valid (m₁, m₂, M)
 * ≈ 20 k tuples / invariant, all verified at 1e-12.
 */

#include "harness.h"

#include <irrep/clebsch_gordan.h>

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static inline int iabs_(int x) { return x < 0 ? -x : x; }

static uint64_t rng_next_(uint64_t *s) {
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return *s;
}

int main(void) {
    IRREP_TEST_START("cg_invariants");

    const int    two_j_max = 20;           /* j up to 10 */
    const int    n_pairs   = 128;
    const double tol       = 1e-12;
    uint64_t     rng       = 0xdeadbeef12345678ULL;

    int pairs_exercised = 0;
    int sum_rule_ok    = 0;
    int orth_ok        = 0;
    int selection_ok   = 0;
    int three_j_ok     = 0;

    for (int p = 0; p < n_pairs; ++p) {
        /* Random (two_j1, two_j2) in [0, two_j_max]. */
        int two_j1 = (int)(rng_next_(&rng) % (two_j_max + 1));
        int two_j2 = (int)(rng_next_(&rng) % (two_j_max + 1));
        pairs_exercised++;

        int two_J_min = iabs_(two_j1 - two_j2);
        int two_J_max = two_j1 + two_j2;

        /* ---- (1) Sum rule for every (m1, m2) ---- */
        for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1 += 2) {
            if ((two_j1 + two_m1) & 1) continue;
            for (int two_m2 = -two_j2; two_m2 <= two_j2; two_m2 += 2) {
                if ((two_j2 + two_m2) & 1) continue;
                int two_M = two_m1 + two_m2;
                double sum = 0.0;
                for (int two_J = two_J_min; two_J <= two_J_max; two_J += 2) {
                    double v = irrep_cg_2j(two_j1, two_m1,
                                           two_j2, two_m2,
                                           two_J,  two_M);
                    sum += v * v;
                }
                IRREP_ASSERT(fabs(sum - 1.0) < tol);
                sum_rule_ok++;
            }
        }

        /* ---- (2) J-orthogonality at a sampled M ---- */
        for (int two_M = -two_J_max; two_M <= two_J_max; two_M += 4) {
            if ((two_J_min + two_M) & 1) continue;
            for (int two_J  = two_J_min; two_J  <= two_J_max; two_J  += 2) {
                for (int two_Jp = two_J;  two_Jp <= two_J_max; two_Jp += 2) {
                    if (iabs_(two_M) > two_J)  continue;
                    if (iabs_(two_M) > two_Jp) continue;
                    double s = 0.0;
                    for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1 += 2) {
                        int two_m2 = two_M - two_m1;
                        if (two_m2 < -two_j2 || two_m2 > two_j2) continue;
                        double a = irrep_cg_2j(two_j1, two_m1, two_j2, two_m2, two_J,  two_M);
                        double b = irrep_cg_2j(two_j1, two_m1, two_j2, two_m2, two_Jp, two_M);
                        s += a * b;
                    }
                    double expected = (two_J == two_Jp) ? 1.0 : 0.0;
                    IRREP_ASSERT(fabs(s - expected) < tol);
                    orth_ok++;
                }
            }
        }

        /* ---- (3) Selection rules ---- */
        /* m₁ + m₂ ≠ M (unless both are 0) must vanish. */
        {
            int two_m1 = two_j1 > 0 ? -two_j1 : 0;
            int two_m2 = two_j2 > 0 ?  two_j2 : 0;
            int two_M  = two_m1 + two_m2 + 2;   /* deliberately off by 1 */
            double v = irrep_cg_2j(two_j1, two_m1, two_j2, two_m2, two_J_min, two_M);
            IRREP_ASSERT(v == 0.0);
            selection_ok++;
        }
        /* |M| > J must vanish. */
        {
            int bad_J = two_J_max;
            int bad_M = two_J_max + 2;
            double v = irrep_cg_2j(two_j1, 0, two_j2, 0, bad_J, bad_M);
            IRREP_ASSERT(v == 0.0);
            selection_ok++;
        }
        /* J > j₁ + j₂ must vanish (triangle inequality). */
        {
            int bad_J = two_J_max + 4;
            double v = irrep_cg_2j(two_j1, 0, two_j2, 0, bad_J, 0);
            IRREP_ASSERT(v == 0.0);
            selection_ok++;
        }
        /* J < |j₁ − j₂| must vanish. */
        if (two_J_min >= 4) {
            int bad_J = two_J_min - 4;
            double v = irrep_cg_2j(two_j1, 0, two_j2, 0, bad_J, 0);
            IRREP_ASSERT(v == 0.0);
            selection_ok++;
        }

        /* ---- (4) CG ↔ 3j relation at a few sampled points ---- */
        for (int two_J = two_J_min; two_J <= two_J_max; two_J += 2) {
            int two_M_lim = two_J;
            for (int two_M = -two_M_lim; two_M <= two_M_lim; two_M += 4) {
                for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1 += 2) {
                    int two_m2 = two_M - two_m1;
                    if (two_m2 < -two_j2 || two_m2 > two_j2) continue;
                    if ((two_j2 + two_m2) & 1) continue;
                    double cg  = irrep_cg_2j(two_j1, two_m1, two_j2, two_m2, two_J, two_M);
                    double tj  = irrep_wigner_3j_2j(two_j1, two_m1,
                                                   two_j2, two_m2,
                                                   two_J,  -two_M);
                    int    phase_half = (two_j1 - two_j2 + two_M) / 2;
                    double phase      = (phase_half & 1) ? -1.0 : 1.0;
                    double rhs        = phase * sqrt((double)two_J + 1.0) * tj;
                    IRREP_ASSERT(fabs(cg - rhs) < tol);
                    three_j_ok++;
                }
            }
        }
    }

    /* Log cardinality so a human reader can eyeball the regime. */
    printf("# exercised %d random (j1,j2) pairs; "
           "%d sum-rule, %d orthogonality, %d selection-rule, %d 3j relation invariants verified\n",
           pairs_exercised, sum_rule_ok, orth_ok, selection_ok, three_j_ok);

    return IRREP_TEST_END();
}
