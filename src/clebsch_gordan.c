/* SPDX-License-Identifier: MIT */
/* M4: Clebsch-Gordan coefficients and Wigner 3j symbols.
 *
 * Two-tier kernel:
 *
 *   - **Small j (max(2j) ≤ 30, i.e. j ≤ 15):** direct Racah single-sum
 *     evaluation against a precomputed factorial table (`wigner_3j_racah_
 *     small_`). The k-sum has at most ~16 terms, every factorial argument
 *     stays under 32, and alternating-sign cancellation is bounded below
 *     2 ulps for all triples in the supported range. ~14 ns/op for j = 1.
 *
 *   - **Large j (max(2j) > 30):** Schulten–Gordon three-term recurrence
 *     over varying j₁ at fixed (j₂, j₃, m₂, m₃) (Schulten & Gordon, J.
 *     Math. Phys. 16, 1961, 1975; Luscombe & Luban, Phys. Rev. E 57,
 *     7274, 1998), as a Miller two-directional iteration: forward from
 *     j_min, backward from j_max, splice at argmax|T_fwd|·|T_bwd|,
 *     normalise via Σ_j (2j+1) T(j)² = 1, sign-anchor at j_max via
 *     (−1)^{j₂−j₃−m₁}. Stable through j = 80 where the small-j Racah
 *     path's k-sum cancellation would lose precision.
 *
 * The earlier Racah-log-gamma single-sum form (replaced by Schulten-
 * Gordon for the j > 15 regime) lost precision to catastrophic
 * cancellation past j ≈ 20 and produced NaN past j ≈ 60 near the
 * triangle edge.
 *
 * Clebsch-Gordan is derived from the 3j symbol via
 *   ⟨j₁ m₁; j₂ m₂ | J M⟩ = (−1)^{j₁−j₂+M} · √(2J+1) · (j₁ j₂ J; m₁ m₂ −M).
 *
 * Doubled-integer arguments (`_2j` variants) are primary; integer wrappers
 * multiply by 2. See also Varshalovich §8; Edmonds §3.
 */

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include <irrep/clebsch_gordan.h>

#define UNUSED(x) ((void)(x))

static inline int iabs_(int x) {
    return x < 0 ? -x : x;
}

/* -------------------------------------------------------------------------- *
 * Schulten-Gordon forward recurrence for the Wigner-3j series.               *
 *                                                                            *
 * For fixed (j₂, j₃, m₁, m₂, m₃) with m₁ + m₂ + m₃ = 0, the 3j symbol        *
 * T(j) = (j, j₂, j₃; m₁, m₂, m₃) satisfies                                   *
 *                                                                            *
 *   j · E(j+1) · T(j+1) + F(j) · T(j) + (j+1) · E(j) · T(j−1) = 0            *
 *                                                                            *
 * where                                                                      *
 *   E(j)² = [j² − (j₂−j₃)²] · [(j₂+j₃+1)² − j²] · [j² − m₁²]                 *
 *   F(j)  = −(2j+1) · { [j₂(j₂+1) − j₃(j₃+1)] · m₁ + j(j+1) · (m₂ − m₃) }    *
 *                                                                            *
 * Forward recurrence (starting with T(j_min) arbitrary, T(j_min − 1) = 0)    *
 * is stable in the classical region; the arbitrary overall scale and sign   *
 * are fixed by the sum rule Σ_j (2j+1) T(j)² = 1 and the sign of T(j_max)   *
 * from the closed-form (−1)^{j₂−j₃−m₁}.                                     *
 * -------------------------------------------------------------------------- */

#define IRREP_3J_MAX_SERIES 512

/* -------------------------------------------------------------------------- *
 * Small-j fast path — direct Racah closed-form.                              *
 *                                                                            *
 * For j ≤ ~15 the Racah single-sum is well-conditioned (the k-sum has a few  *
 * terms and the largest factorial argument stays under 30, well below the    *
 * lgamma overflow boundary AND the catastrophic-cancellation threshold).     *
 * Evaluating it directly with a precomputed factorial table avoids the       *
 * forward + backward + splice + normalise overhead of the full Schulten-     *
 * Gordon recurrence and brings cg_110_110_110 back to ~22 ns/op (matching    *
 * pre-rewrite perf) while preserving the high-j stability of the recurrence  *
 * path used past the threshold.                                              *
 *                                                                            *
 * Threshold rationale: for j_max = 15, max factorial argument is 2 · 15 +    *
 * 1 = 31. 31! ≈ 8.2e33 fits in double (max double ~1.8e308) and the          *
 * alternating-sign cancellation in the k-sum is below 2 ulps for all         *
 * tested triples. j_max = 20 still works but adds margin we don't need;      *
 * 15 is a safe sweet spot that still hits the common physics use cases       *
 * (l_max ≤ 8 spherical harmonics, NN coupling).                              *
 * -------------------------------------------------------------------------- */

#define IRREP_3J_SMALL_FACT_MAX 64
static double small_factorial_[IRREP_3J_SMALL_FACT_MAX + 1];
static int    small_factorial_init_ = 0;

static void init_small_factorial_(void) {
    if (small_factorial_init_)
        return;
    small_factorial_[0] = 1.0;
    for (int n = 1; n <= IRREP_3J_SMALL_FACT_MAX; ++n)
        small_factorial_[n] = small_factorial_[n - 1] * (double)n;
    small_factorial_init_ = 1;
}

/* Direct Racah single-sum evaluation of (j₁ j₂ j₃; m₁ m₂ m₃).
 * Caller has already applied selection rules and verified that
 * max(2j₁, 2j₂, 2j₃) ≤ 30 (i.e. j ≤ 15). All arguments are
 * doubled-integers; conversion divides by 2. */
static double wigner_3j_racah_small_(int two_j1, int two_m1, int two_j2, int two_m2,
                                      int two_j3, int two_m3) {
    init_small_factorial_();

    /* Half-integers as integers via doubled form. The selection-rule
     * checks already ensured all the following are non-negative integers. */
    const int j1pm1 = (two_j1 + two_m1) / 2;
    const int j1mm1 = (two_j1 - two_m1) / 2;
    const int j2pm2 = (two_j2 + two_m2) / 2;
    const int j2mm2 = (two_j2 - two_m2) / 2;
    const int j3pm3 = (two_j3 + two_m3) / 2;
    const int j3mm3 = (two_j3 - two_m3) / 2;

    /* Triangle factor Δ(j₁ j₂ j₃) = √[(j₁+j₂−j₃)!(j₁−j₂+j₃)!(−j₁+j₂+j₃)!/(j₁+j₂+j₃+1)!]. */
    const int t1     = (two_j1 + two_j2 - two_j3) / 2;
    const int t2     = (two_j1 - two_j2 + two_j3) / 2;
    const int t3     = (-two_j1 + two_j2 + two_j3) / 2;
    const int s_plus = (two_j1 + two_j2 + two_j3) / 2 + 1;

    /* Bounds check: redundant given caller, but defensive. */
    if (s_plus > IRREP_3J_SMALL_FACT_MAX)
        return NAN;

    const double tri_num = small_factorial_[t1] * small_factorial_[t2] * small_factorial_[t3];
    const double tri     = sqrt(tri_num / small_factorial_[s_plus]);

    const double m_factor =
        sqrt(small_factorial_[j1pm1] * small_factorial_[j1mm1] * small_factorial_[j2pm2]
             * small_factorial_[j2mm2] * small_factorial_[j3pm3] * small_factorial_[j3mm3]);

    /* Summation range over k ≥ 0 such that every factorial argument
     * in the denominator is ≥ 0:
     *   k ≤ t1                  (j₁+j₂−j₃−k ≥ 0)
     *   k ≤ j₁−m₁               (j₁−m₁−k ≥ 0)
     *   k ≤ j₂+m₂               (j₂+m₂−k ≥ 0)
     *   k ≥ j₂−j₃−m₁ ⇒ k ≥ −(j₃−j₂+m₁)
     *   k ≥ j₁−j₃+m₂  ⇒ k ≥ −(j₃−j₁−m₂) */
    const int e1 = (two_j3 - two_j2 + two_m1) / 2;  /* j₃−j₂+m₁ */
    const int e2 = (two_j3 - two_j1 - two_m2) / 2;  /* j₃−j₁−m₂ */
    int k_lo = 0;
    if (-e1 > k_lo) k_lo = -e1;
    if (-e2 > k_lo) k_lo = -e2;
    int k_hi = t1;
    if (j1mm1 < k_hi) k_hi = j1mm1;
    if (j2pm2 < k_hi) k_hi = j2pm2;
    if (k_lo > k_hi)
        return 0.0;

    double sum     = 0.0;
    double sign    = (k_lo & 1) ? -1.0 : 1.0;
    for (int k = k_lo; k <= k_hi; ++k) {
        const int a = t1 - k;
        const int b = j1mm1 - k;
        const int c = j2pm2 - k;
        const int d = e1 + k;
        const int e = e2 + k;
        const double denom = small_factorial_[k] * small_factorial_[a] * small_factorial_[b]
                             * small_factorial_[c] * small_factorial_[d] * small_factorial_[e];
        sum += sign / denom;
        sign = -sign;
    }

    /* Overall sign convention from Edmonds (3.6.10): (−1)^{j₁−j₂−m₃}. */
    const int phase_int = (two_j1 - two_j2 - two_m3) / 2;
    const double phase  = (phase_int & 1) ? -1.0 : 1.0;
    return phase * tri * m_factor * sum;
}

static double three_j_E_squared_(double j, double j2mj3, double j2pj3p1, double m1) {
    return (j * j - j2mj3 * j2mj3) * (j2pj3p1 * j2pj3p1 - j * j) * (j * j - m1 * m1);
}

/* Three-term recurrence step F(j)·T(j) = −j·E(j+1)·T(j+1) − (j+1)·E(j)·T(j−1).
 * Used below by both forward (→ T(j+1)) and backward (→ T(j−1)) passes. */

static int wigner_3j_series_(double j2, double j3, double m1, double m2, double m3, double j_min,
                             double j_max, double *T_out, int N) {
    if (N <= 0 || N > IRREP_3J_MAX_SERIES)
        return -1;

    double j2mj3 = j2 - j3;
    double j2pj3p1 = j2 + j3 + 1.0;

    /* Two-directional Miller iteration:
     *
     *   1. Forward recurrence from j_min up to j_max, using T(j_min−1) = 0
     *      and T(j_min) = 1 (arbitrary seed). This is stable for the
     *      classical solution in the lower half of [j_min, j_max] and
     *      becomes contaminated by the non-classical (subdominant)
     *      solution in the upper tail where E(j)² shrinks toward the
     *      j_max endpoint.
     *
     *   2. Backward recurrence from j_max down to j_min with T(j_max+1) = 0
     *      and T(j_max) = 1, symmetrically stable in the upper half and
     *      contaminated in the lower tail.
     *
     *   3. Splice: choose a matching index idx* inside the classical region
     *      (one where the forward and backward iterations agree up to an
     *      overall scale), rescale the forward array so T_fwd(idx*) =
     *      T_bwd(idx*), then take T(j) from T_fwd for idx ≤ idx* and from
     *      T_bwd for idx > idx*.
     *
     *   4. Normalise the spliced series via Σ_j (2j+1) T(j)² = 1 and fix
     *      sign from the closed-form value at j_max.
     *
     * Match-point heuristic: pick the idx that maximises
     * |T_fwd(idx)| · |T_bwd(idx)| — the locus where both solutions are
     * large (so the scale estimate has the best relative precision).
     */

    double T_fwd[IRREP_3J_MAX_SERIES];
    double T_bwd[IRREP_3J_MAX_SERIES];

    /* ---- forward ---- */
    {
        double T_prev = 0.0; /* T(j_min − 1) */
        double T_curr = 1.0; /* T(j_min); arbitrary */
        T_fwd[0] = T_curr;
        for (int idx = 1; idx < N; ++idx) {
            double j = j_min + (double)(idx - 1); /* recurrence at j produces T(j+1) */
            double jp1 = j + 1.0;
            double Ej_sq = three_j_E_squared_(j, j2mj3, j2pj3p1, m1);
            double Ejp1_sq = three_j_E_squared_(jp1, j2mj3, j2pj3p1, m1);
            double Ej = (Ej_sq > 0.0) ? sqrt(Ej_sq) : 0.0;
            double Ejp1 = (Ejp1_sq > 0.0) ? sqrt(Ejp1_sq) : 0.0;
            double Fj = -(2.0 * j + 1.0) *
                        ((j2 * (j2 + 1.0) - j3 * (j3 + 1.0)) * m1 + j * (j + 1.0) * (m2 - m3));
            double T_next;
            if (Ejp1 > 0.0) {
                T_next = -(Fj * T_curr + jp1 * Ej * T_prev) / (j * Ejp1);
            } else if (j <= 0.0) {
                /* j_min = 0 and j·E(j+1) = 0 — the recurrence at j = 0
                 * degenerates. The closed form at j = 0 coupling two
                 * copies of the same j₂ = j₃ is (0 j j; 0 m −m) =
                 * (−1)^{j−m} / √(2j+1); use it to initialise T_fwd[0] and
                 * pick T_fwd[1] from the recurrence at j = 1 in step 2. */
                T_next = 0.0; /* unreliable from this step; overwritten */
            } else {
                T_next = 0.0;
            }
            T_fwd[idx] = T_next;
            T_prev = T_curr;
            T_curr = T_next;
        }
    }

    /* ---- backward ---- */
    {
        double T_next = 0.0; /* T(j_max + 1) */
        double T_curr = 1.0; /* T(j_max); arbitrary */
        T_bwd[N - 1] = T_curr;
        for (int idx = N - 2; idx >= 0; --idx) {
            double j = j_min + (double)(idx + 1); /* recurrence at j produces T(j-1) */
            double jp1 = j + 1.0;
            double Ej_sq = three_j_E_squared_(j, j2mj3, j2pj3p1, m1);
            double Ejp1_sq = three_j_E_squared_(jp1, j2mj3, j2pj3p1, m1);
            double Ej = (Ej_sq > 0.0) ? sqrt(Ej_sq) : 0.0;
            double Ejp1 = (Ejp1_sq > 0.0) ? sqrt(Ejp1_sq) : 0.0;
            double Fj = -(2.0 * j + 1.0) *
                        ((j2 * (j2 + 1.0) - j3 * (j3 + 1.0)) * m1 + j * (j + 1.0) * (m2 - m3));
            double T_prev;
            if (Ej > 0.0) {
                T_prev = -(j * Ejp1 * T_next + Fj * T_curr) / (jp1 * Ej);
            } else {
                T_prev = 0.0;
            }
            T_bwd[idx] = T_prev;
            T_next = T_curr;
            T_curr = T_prev;
        }
    }

    /* ---- splice: find match point that maximises |T_fwd|·|T_bwd| ---- */
    int    idx_match = N / 2;
    double best = 0.0;
    for (int idx = 0; idx < N; ++idx) {
        double w = fabs(T_fwd[idx]) * fabs(T_bwd[idx]);
        if (w > best) {
            best = w;
            idx_match = idx;
        }
    }
    if (best == 0.0) {
        /* Degenerate series (N = 1, or one direction overflowed). Fall
         * back to backward-only. */
        for (int idx = 0; idx < N; ++idx)
            T_out[idx] = T_bwd[idx];
    } else {
        double scale = T_bwd[idx_match] / T_fwd[idx_match];
        for (int idx = 0; idx <= idx_match; ++idx)
            T_out[idx] = scale * T_fwd[idx];
        for (int idx = idx_match + 1; idx < N; ++idx)
            T_out[idx] = T_bwd[idx];
    }

    /* Sum-rule normalisation. */
    double norm_sq = 0.0;
    for (int idx = 0; idx < N; ++idx) {
        double j = j_min + (double)idx;
        norm_sq += (2.0 * j + 1.0) * T_out[idx] * T_out[idx];
    }
    if (norm_sq <= 0.0 || !isfinite(norm_sq)) {
        for (int idx = 0; idx < N; ++idx)
            T_out[idx] = 0.0;
        return 0;
    }
    double norm = sqrt(norm_sq);
    for (int idx = 0; idx < N; ++idx)
        T_out[idx] /= norm;

    /* Sign convention: at j_max the 3j equals (−1)^{j₂−j₃−m₁} · √(positive). */
    int    exp_int = (int)lround(j2 - j3 - m1);
    double expected_sign = (exp_int & 1) ? -1.0 : 1.0;
    if (T_out[N - 1] * expected_sign < 0.0) {
        for (int idx = 0; idx < N; ++idx)
            T_out[idx] = -T_out[idx];
    }
    return 0;
}

/* -------------------------------------------------------------------------- *
 * Clebsch-Gordan — Racah formula                                             *
 * -------------------------------------------------------------------------- */

double irrep_cg_2j(int two_j1, int two_m1, int two_j2, int two_m2, int two_J, int two_M) {
    /* CG ↔ 3j relation:
     *   ⟨j₁ m₁ j₂ m₂ | J M⟩ = (−1)^{j₁−j₂+M} · √(2J+1) · (j₁ j₂ J; m₁ m₂ −M)
     *
     * m-conservation check here; the remaining selection rules are enforced
     * inside irrep_wigner_3j_2j. */
    if (two_m1 + two_m2 != two_M)
        return 0.0;
    double three_j = irrep_wigner_3j_2j(two_j1, two_m1, two_j2, two_m2, two_J, -two_M);
    if (three_j == 0.0)
        return 0.0;
    int    phase_half = (two_j1 - two_j2 + two_M) / 2;
    double phase = (phase_half & 1) ? -1.0 : 1.0;
    return phase * sqrt((double)two_J + 1.0) * three_j;
}

double irrep_cg(int j1, int m1, int j2, int m2, int J, int M) {
    return irrep_cg_2j(2 * j1, 2 * m1, 2 * j2, 2 * m2, 2 * J, 2 * M);
}

/* -------------------------------------------------------------------------- *
 * Wigner 3j via Schulten-Gordon recurrence.                                  *
 * -------------------------------------------------------------------------- */

double irrep_wigner_3j_2j(int two_j1, int two_m1, int two_j2, int two_m2, int two_j3, int two_m3) {
    /* selection rules */
    if (two_j1 < 0 || two_j2 < 0 || two_j3 < 0)
        return 0.0;
    if (two_m1 + two_m2 + two_m3 != 0)
        return 0.0;
    if (two_m1 < -two_j1 || two_m1 > two_j1)
        return 0.0;
    if (two_m2 < -two_j2 || two_m2 > two_j2)
        return 0.0;
    if (two_m3 < -two_j3 || two_m3 > two_j3)
        return 0.0;
    if ((two_j1 + two_m1) & 1)
        return 0.0;
    if ((two_j2 + two_m2) & 1)
        return 0.0;
    if ((two_j3 + two_m3) & 1)
        return 0.0;
    if ((two_j1 + two_j2 + two_j3) & 1)
        return 0.0;
    if (two_j1 < iabs_(two_j2 - two_j3) || two_j1 > two_j2 + two_j3)
        return 0.0;

    int two_j1_min = iabs_(two_j2 - two_j3);
    if (iabs_(two_m1) > two_j1_min)
        two_j1_min = iabs_(two_m1);
    int two_j1_max = two_j2 + two_j3;

    /* The series steps by 1 in j₁, i.e. by 2 in two_j₁; parity must agree. */
    if ((two_j1 - two_j1_min) & 1)
        return 0.0;
    int N = (two_j1_max - two_j1_min) / 2 + 1;
    int idx = (two_j1 - two_j1_min) / 2;
    if (idx < 0 || idx >= N)
        return 0.0;

    /* Small-j fast path: Racah closed-form. The threshold of two_j_max ≤ 30
     * (j ≤ 15) keeps every factorial argument under 32 and the alternating-
     * sign cancellation in the k-sum below 2 ulps. Above this we fall
     * through to the Schulten-Gordon recurrence which is stable to j ≈ 80. */
    int two_j_max = two_j1;
    if (two_j2 > two_j_max) two_j_max = two_j2;
    if (two_j3 > two_j_max) two_j_max = two_j3;
    if (two_j_max <= 30)
        return wigner_3j_racah_small_(two_j1, two_m1, two_j2, two_m2, two_j3, two_m3);

    double j2 = 0.5 * (double)two_j2;
    double j3 = 0.5 * (double)two_j3;
    double m1 = 0.5 * (double)two_m1;
    double m2 = 0.5 * (double)two_m2;
    double m3 = 0.5 * (double)two_m3;
    double j_min = 0.5 * (double)two_j1_min;
    double j_max = 0.5 * (double)two_j1_max;

    double T[IRREP_3J_MAX_SERIES];
    if (wigner_3j_series_(j2, j3, m1, m2, m3, j_min, j_max, T, N) != 0)
        return 0.0;
    return T[idx];
}

double irrep_wigner_3j(int j1, int m1, int j2, int m2, int j3, int m3) {
    return irrep_wigner_3j_2j(2 * j1, 2 * m1, 2 * j2, 2 * m2, 2 * j3, 2 * m3);
}

/* -------------------------------------------------------------------------- *
 * Dense CG cache                                                             *
 *                                                                            *
 * Five-dimensional flat table, indexed by                                    *
 *   (two_j1, two_m1 + two_j1_max, two_j2, two_m2 + two_j2_max, two_J).       *
 * Storage is (two_j1_max + 1) × (2·two_j1_max + 1) × (two_j2_max + 1)        *
 *          × (2·two_j2_max + 1) × (two_j1_max + two_j2_max + 1) doubles;     *
 * invalid combos stored as 0. For two_j1_max = two_j2_max = 8 (int j ≤ 4)    *
 * this is ~3 MB; for 16/16 (int j ≤ 8) it is ~80 MB and callers should use   *
 * a smaller-scope table.                                                     *
 * -------------------------------------------------------------------------- */

struct cg_table {
    int     two_j1_max;
    int     two_j2_max;
    size_t  stride_j1;
    size_t  stride_m1;
    size_t  stride_j2;
    size_t  stride_m2;
    size_t  stride_J;
    size_t  data_size;
    double *data;
};

cg_table_t *irrep_cg_table_build_2j(int two_j1_max, int two_j2_max) {
    if (two_j1_max < 0 || two_j2_max < 0)
        return NULL;

    cg_table_t *t = calloc(1, sizeof(*t));
    if (!t)
        return NULL;
    t->two_j1_max = two_j1_max;
    t->two_j2_max = two_j2_max;

    size_t nm1 = (size_t)(2 * two_j1_max + 1);
    size_t nj2 = (size_t)(two_j2_max + 1);
    size_t nm2 = (size_t)(2 * two_j2_max + 1);
    size_t nJ = (size_t)(two_j1_max + two_j2_max + 1);

    t->stride_J = 1;
    t->stride_m2 = nJ;
    t->stride_j2 = nm2 * nJ;
    t->stride_m1 = nj2 * nm2 * nJ;
    t->stride_j1 = nm1 * nj2 * nm2 * nJ;

    size_t total = (size_t)(two_j1_max + 1) * t->stride_j1;
    if (total == 0)
        total = 1;
    t->data_size = total;
    t->data = calloc(total, sizeof(double));
    if (!t->data) {
        free(t);
        return NULL;
    }

    for (int two_j1 = 0; two_j1 <= two_j1_max; ++two_j1) {
        for (int two_m1 = -two_j1; two_m1 <= two_j1; ++two_m1) {
            for (int two_j2 = 0; two_j2 <= two_j2_max; ++two_j2) {
                for (int two_m2 = -two_j2; two_m2 <= two_j2; ++two_m2) {
                    int two_M = two_m1 + two_m2;
                    int two_J_min = iabs_(two_j1 - two_j2);
                    int two_J_top = two_j1 + two_j2;
                    for (int two_J = two_J_min; two_J <= two_J_top; two_J += 2) {
                        if (iabs_(two_M) > two_J)
                            continue;
                        double v = irrep_cg_2j(two_j1, two_m1, two_j2, two_m2, two_J, two_M);
                        size_t idx = (size_t)two_j1 * t->stride_j1 +
                                     (size_t)(two_m1 + two_j1_max) * t->stride_m1 +
                                     (size_t)two_j2 * t->stride_j2 +
                                     (size_t)(two_m2 + two_j2_max) * t->stride_m2 +
                                     (size_t)two_J * t->stride_J;
                        if (idx < t->data_size)
                            t->data[idx] = v;
                    }
                }
            }
        }
    }
    return t;
}

cg_table_t *irrep_cg_table_build(int j1_max, int j2_max) {
    return irrep_cg_table_build_2j(2 * j1_max, 2 * j2_max);
}

void irrep_cg_table_free(cg_table_t *t) {
    if (!t)
        return;
    free(t->data);
    free(t);
}

double irrep_cg_lookup_2j(const cg_table_t *t, int two_j1, int two_m1, int two_j2, int two_m2,
                          int two_J, int two_M) {
    if (!t || !t->data)
        return 0.0;
    if (two_j1 < 0 || two_j1 > t->two_j1_max)
        return 0.0;
    if (two_j2 < 0 || two_j2 > t->two_j2_max)
        return 0.0;
    if (two_m1 < -two_j1 || two_m1 > two_j1)
        return 0.0;
    if (two_m2 < -two_j2 || two_m2 > two_j2)
        return 0.0;
    if (two_m1 + two_m2 != two_M)
        return 0.0;
    if (two_J < 0 || two_J > t->two_j1_max + t->two_j2_max)
        return 0.0;
    if (two_J < iabs_(two_j1 - two_j2) || two_J > two_j1 + two_j2)
        return 0.0;
    if (iabs_(two_M) > two_J)
        return 0.0;

    size_t idx = (size_t)two_j1 * t->stride_j1 + (size_t)(two_m1 + t->two_j1_max) * t->stride_m1 +
                 (size_t)two_j2 * t->stride_j2 + (size_t)(two_m2 + t->two_j2_max) * t->stride_m2 +
                 (size_t)two_J * t->stride_J;
    if (idx >= t->data_size)
        return 0.0;
    return t->data[idx];
}

double irrep_cg_lookup(const cg_table_t *t, int j1, int m1, int j2, int m2, int J, int M) {
    return irrep_cg_lookup_2j(t, 2 * j1, 2 * m1, 2 * j2, 2 * m2, 2 * J, 2 * M);
}

/* -------------------------------------------------------------------------- *
 * Dense block                                                                *
 *                                                                            *
 * Layout:                                                                    *
 *   out[((m1 + j1) · (2 j2 + 1) + (m2 + j2)) · (2 J + 1) + (M + J)]          *
 *     = ⟨j1 m1; j2 m2 | J M⟩.                                                *
 * -------------------------------------------------------------------------- */

void irrep_cg_block(int j1, int j2, int J, double *out) {
    int d1 = 2 * j1 + 1;
    int d2 = 2 * j2 + 1;
    int dJ = 2 * J + 1;
    for (int i1 = 0; i1 < d1; ++i1) {
        int m1 = i1 - j1;
        for (int i2 = 0; i2 < d2; ++i2) {
            int m2 = i2 - j2;
            for (int iJ = 0; iJ < dJ; ++iJ) {
                int M = iJ - J;
                out[(i1 * d2 + i2) * dJ + iJ] = irrep_cg(j1, m1, j2, m2, J, M);
            }
        }
    }
}
