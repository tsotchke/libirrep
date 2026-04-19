/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/clebsch_gordan.h>
#include "reference_data/cg_reference.h"
#include <math.h>

static int iabs_(int x) { return x < 0 ? -x : x; }

int main(void) {
    IRREP_TEST_START("clebsch_gordan");

    const double tol = 1e-12;

    /* -------- hand-tabulated values (Sakurai Appendix A) -------- */

    /* ⟨½ ½; ½ −½ | 1 0⟩ = 1/√2 */
    IRREP_ASSERT_NEAR(irrep_cg_2j(1, 1, 1, -1, 2, 0),  1.0 / sqrt(2.0), tol);
    /* ⟨½ −½; ½ ½ | 1 0⟩ = 1/√2 */
    IRREP_ASSERT_NEAR(irrep_cg_2j(1, -1, 1, 1, 2, 0),  1.0 / sqrt(2.0), tol);
    /* ⟨½ ½; ½ −½ | 0 0⟩ = 1/√2 */
    IRREP_ASSERT_NEAR(irrep_cg_2j(1, 1, 1, -1, 0, 0),  1.0 / sqrt(2.0), tol);
    /* ⟨½ −½; ½ ½ | 0 0⟩ = −1/√2 */
    IRREP_ASSERT_NEAR(irrep_cg_2j(1, -1, 1, 1, 0, 0), -1.0 / sqrt(2.0), tol);
    /* ⟨½ ½; ½ ½ | 1 1⟩ = 1 */
    IRREP_ASSERT_NEAR(irrep_cg_2j(1, 1, 1, 1, 2, 2),   1.0,             tol);

    /* ⟨1 0; 1 0 | 0 0⟩ = −1/√3;  ⟨1 1; 1 −1 | 0 0⟩ = 1/√3 */
    IRREP_ASSERT_NEAR(irrep_cg(1, 0, 1,  0, 0, 0), -1.0 / sqrt(3.0), tol);
    IRREP_ASSERT_NEAR(irrep_cg(1, 1, 1, -1, 0, 0),  1.0 / sqrt(3.0), tol);
    /* ⟨1 0; 1 0 | 2 0⟩ = √(2/3) */
    IRREP_ASSERT_NEAR(irrep_cg(1, 0, 1,  0, 2, 0),  sqrt(2.0/3.0),   tol);
    /* ⟨1 1; 1 −1 | 2 0⟩ = √(1/6) */
    IRREP_ASSERT_NEAR(irrep_cg(1, 1, 1, -1, 2, 0),  sqrt(1.0/6.0),   tol);
    /* ⟨1 1; 1 0 | 2 1⟩ = 1/√2 */
    IRREP_ASSERT_NEAR(irrep_cg(1, 1, 1,  0, 2, 1),  1.0 / sqrt(2.0), tol);

    /* ⟨3/2 1/2; 1 0 | 3/2 1/2⟩ = 2/3 · √(1/5) ... let's use a simpler half-int case */
    /* ⟨1 1; ½ −½ | 3/2 ½⟩ = √(1/3)  [Varshalovich 8.12 / Sakurai] */
    IRREP_ASSERT_NEAR(irrep_cg_2j(2, 2, 1, -1, 3, 1), sqrt(1.0/3.0), tol);
    /* ⟨1 0; ½ ½ | 3/2 ½⟩ = √(2/3) */
    IRREP_ASSERT_NEAR(irrep_cg_2j(2, 0, 1,  1, 3, 1), sqrt(2.0/3.0), tol);
    /* ⟨1 1; ½ ½ | 3/2 3/2⟩ = 1 */
    IRREP_ASSERT_NEAR(irrep_cg_2j(2, 2, 1,  1, 3, 3), 1.0,           tol);

    /* -------- selection rules return 0 -------- */
    IRREP_ASSERT_NEAR(irrep_cg(1, 0, 1, 0, 5, 0), 0.0, tol);   /* J > j1+j2 */
    IRREP_ASSERT_NEAR(irrep_cg(1, 0, 1, 0, 1, 0), 0.0, tol);   /* wait, this is 0 due to selection — (1,0,1,0|1,0) actually vanishes by parity of (j1+j2+J)=3 */
    IRREP_ASSERT_NEAR(irrep_cg(1, 1, 1, 0, 2, 0), 0.0, tol);   /* m1+m2≠M */
    IRREP_ASSERT_NEAR(irrep_cg(1, 2, 1, 0, 2, 0), 0.0, tol);   /* |m1|>j1 */

    /* -------- completeness: Σ_J |⟨j1 m1; j2 m2 | J M⟩|² = 1 -------- */
    for (int two_j1 = 0; two_j1 <= 6; ++two_j1) {
        for (int two_j2 = 0; two_j2 <= 6; ++two_j2) {
            for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1 += 2) {
                for (int two_m2 = -two_j2; two_m2 <= two_j2; two_m2 += 2) {
                    int two_M = two_m1 + two_m2;
                    double s = 0.0;
                    int two_J_min = iabs_(two_j1 - two_j2);
                    int two_J_max = two_j1 + two_j2;
                    for (int two_J = two_J_min; two_J <= two_J_max; two_J += 2) {
                        double v = irrep_cg_2j(two_j1, two_m1,
                                               two_j2, two_m2,
                                               two_J,  two_M);
                        s += v * v;
                    }
                    IRREP_ASSERT_NEAR(s, 1.0, 1e-10);
                }
            }
        }
    }

    /* -------- orthogonality: Σ_{m1} CG(..|JM) · CG(..|J'M) = δ_JJ' -------- */
    {
        int two_j1 = 2, two_j2 = 2;
        int two_M = 0;
        int two_J_a = 0, two_J_b = 2;
        double s_aa = 0.0, s_ab = 0.0, s_bb = 0.0;
        for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1 += 2) {
            int two_m2 = two_M - two_m1;
            if (two_m2 < -two_j2 || two_m2 > two_j2) continue;
            double va = irrep_cg_2j(two_j1, two_m1, two_j2, two_m2, two_J_a, two_M);
            double vb = irrep_cg_2j(two_j1, two_m1, two_j2, two_m2, two_J_b, two_M);
            s_aa += va * va;
            s_ab += va * vb;
            s_bb += vb * vb;
        }
        IRREP_ASSERT_NEAR(s_aa, 1.0, 1e-10);
        IRREP_ASSERT_NEAR(s_bb, 1.0, 1e-10);
        IRREP_ASSERT_NEAR(s_ab, 0.0, 1e-10);
    }

    /* -------- 3j symbol matches CG via the standard phase/normalization -------- */
    for (int two_j1 = 0; two_j1 <= 4; ++two_j1) {
        for (int two_j2 = 0; two_j2 <= 4; ++two_j2) {
            int two_J_min = iabs_(two_j1 - two_j2);
            int two_J_max = two_j1 + two_j2;
            for (int two_j3 = two_J_min; two_j3 <= two_J_max; two_j3 += 2) {
                for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1 += 2) {
                    for (int two_m2 = -two_j2; two_m2 <= two_j2; two_m2 += 2) {
                        int two_m3 = -(two_m1 + two_m2);
                        if (two_m3 < -two_j3 || two_m3 > two_j3) continue;
                        double w  = irrep_wigner_3j_2j(two_j1, two_m1,
                                                       two_j2, two_m2,
                                                       two_j3, two_m3);
                        double cg = irrep_cg_2j(two_j1, two_m1,
                                                two_j2, two_m2,
                                                two_j3, -two_m3);
                        int ph = (two_j1 - two_j2 - two_m3) / 2;
                        double phase = (ph & 1) ? -1.0 : 1.0;
                        double expected = phase * cg / sqrt((double)two_j3 + 1.0);
                        IRREP_ASSERT_NEAR(w, expected, 1e-12);
                    }
                }
            }
        }
    }

    /* -------- 3j cyclic column symmetry: (j1 j2 j3; m1 m2 m3) = (j2 j3 j1; m2 m3 m1) -------- */
    {
        double a = irrep_wigner_3j(1, 1, 2, -1, 1, 0);
        double b = irrep_wigner_3j(2, -1, 1, 0, 1, 1);
        IRREP_ASSERT_NEAR(a, b, 1e-12);
    }

    /* -------- table lookup matches direct computation -------- */
    {
        cg_table_t *t = irrep_cg_table_build(3, 3);
        IRREP_ASSERT(t != NULL);
        int mismatches = 0;
        for (int j1 = 0; j1 <= 3; ++j1) {
            for (int j2 = 0; j2 <= 3; ++j2) {
                for (int m1 = -j1; m1 <= j1; ++m1) {
                    for (int m2 = -j2; m2 <= j2; ++m2) {
                        int M = m1 + m2;
                        for (int J = iabs_(j1 - j2); J <= j1 + j2; ++J) {
                            if (iabs_(M) > J) continue;
                            double direct = irrep_cg(j1, m1, j2, m2, J, M);
                            double cached = irrep_cg_lookup(t, j1, m1, j2, m2, J, M);
                            if (fabs(direct - cached) > 1e-14) mismatches++;
                        }
                    }
                }
            }
        }
        IRREP_ASSERT(mismatches == 0);
        irrep_cg_table_free(t);
    }

    /* -------- cg_block fills the correct (m1, m2, M) layout -------- */
    {
        int j1 = 1, j2 = 1, J = 2;
        int d1 = 2*j1 + 1, d2 = 2*j2 + 1, dJ = 2*J + 1;
        double buf[3 * 3 * 5];
        irrep_cg_block(j1, j2, J, buf);
        for (int i1 = 0; i1 < d1; ++i1) {
            int m1 = i1 - j1;
            for (int i2 = 0; i2 < d2; ++i2) {
                int m2 = i2 - j2;
                for (int iJ = 0; iJ < dJ; ++iJ) {
                    int M = iJ - J;
                    double v = buf[(i1 * d2 + i2) * dJ + iJ];
                    double expected = irrep_cg(j1, m1, j2, m2, J, M);
                    IRREP_ASSERT(fabs(v - expected) < 1e-14);
                }
            }
        }
    }

    /* -------- large j numerical stability (j = 10, various m) -------- */
    {
        /* Not against a closed form, just verify finite values and sum rule. */
        int two_j = 20;  /* j = 10 */
        double s = 0.0;
        for (int two_J = 0; two_J <= 2 * two_j; two_J += 2) {
            double v = irrep_cg_2j(two_j, 0, two_j, 0, two_J, 0);
            IRREP_ASSERT(isfinite(v));
            s += v * v;
        }
        IRREP_ASSERT_NEAR(s, 1.0, 1e-10);
    }

    /* -------- Hand-tabulated reference sweep (Sakurai / Varshalovich) -------- *
     * Check every entry in tests/reference_data/cg_reference.h to 1e-12. This
     * catches silent drift — e.g., a convention flip — that sum rules wouldn't. */
    for (int i = 0; i < IRREP_CG_NUM_REFERENCES; ++i) {
        const irrep_cg_reference_t *ref = &IRREP_CG_REFERENCES[i];
        double expected = (ref->sign == 0)
            ? 0.0
            : (double)ref->sign * sqrt((double)ref->num / (double)ref->den);
        double got = irrep_cg_2j(ref->two_j1, ref->two_m1,
                                 ref->two_j2, ref->two_m2,
                                 ref->two_J,  ref->two_M);
        IRREP_ASSERT(fabs(got - expected) < 1e-12);
    }

    return IRREP_TEST_END();
}
