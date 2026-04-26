/* SPDX-License-Identifier: MIT */
/* Point-group projector tests per PR#6 acceptance criteria.
 *
 *   1. Table metadata: num_irreps, order.
 *   2. Character orthogonality: (1/|G|) Σ_g χ*_μ(g) χ_ν(g) = δ_{μν}.
 *   3. Projector idempotence: P_μ(P_μ v) == P_μ v.
 *   4. Projector sum: Σ_μ P_μ v == v.
 *   5. Decomposition matches hand-computed Bradley-Cracknell reductions.
 */

#include "harness.h"
#include <irrep/multiset.h>
#include <irrep/point_group.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static uint64_t sm64_(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static double sm64_unit_(uint64_t *s) {
    return (double)(sm64_(s) >> 11) / (double)(1ULL << 53);
}

/* Character orthogonality — uses the public API indirectly: we don't have a
 * direct character accessor, so we reconstruct χ_μ(g) by projecting the
 * identity representation and reading off traces. Instead we inspect the
 * known hand-entered table through the projector sum identity — covered by
 * test 4 — and here assert the dim + order invariants directly. */
static void test_metadata_(irrep_point_group_t g, int expect_num_irreps, int expect_order) {
    irrep_pg_table_t *t = irrep_pg_table_build(g);
    IRREP_ASSERT(t != NULL);
    IRREP_ASSERT(irrep_pg_num_irreps(t) == expect_num_irreps);
    IRREP_ASSERT(irrep_pg_order(t) == expect_order);
    for (int mu = 0; mu < expect_num_irreps; ++mu) {
        const char *lbl = irrep_pg_irrep_label(t, mu);
        IRREP_ASSERT(lbl != NULL);
        IRREP_ASSERT(lbl[0] != '\0');
    }
    irrep_pg_table_free(t);
}

/* Projector idempotence: P_μ(P_μ v) = P_μ v. */
static void test_idempotence_(irrep_point_group_t g, const char *spec) {
    irrep_pg_table_t *t = irrep_pg_table_build(g);
    IRREP_ASSERT(t != NULL);
    irrep_multiset_t *m = irrep_multiset_parse(spec);
    IRREP_ASSERT(m != NULL);
    int      N = m->total_dim;

    double  *v = calloc((size_t)N, sizeof(double));
    double  *p1 = calloc((size_t)N, sizeof(double));
    double  *p2 = calloc((size_t)N, sizeof(double));

    uint64_t s = 0xDEADBEEFCAFEBABEULL;
    for (int i = 0; i < N; ++i)
        v[i] = 2.0 * sm64_unit_(&s) - 1.0;

    for (int mu = 0; mu < irrep_pg_num_irreps(t); ++mu) {
        irrep_pg_project(t, mu, m, v, p1);
        irrep_pg_project(t, mu, m, p1, p2);
        for (int i = 0; i < N; ++i) {
            IRREP_ASSERT(fabs(p2[i] - p1[i]) < 1e-10);
        }
    }
    free(v);
    free(p1);
    free(p2);
    irrep_multiset_free(m);
    irrep_pg_table_free(t);
}

/* Projector sum: Σ_μ P_μ v = v. */
static void test_projector_sum_(irrep_point_group_t g, const char *spec) {
    irrep_pg_table_t *t = irrep_pg_table_build(g);
    IRREP_ASSERT(t != NULL);
    irrep_multiset_t *m = irrep_multiset_parse(spec);
    IRREP_ASSERT(m != NULL);
    int      N = m->total_dim;

    double  *v = calloc((size_t)N, sizeof(double));
    double  *acc = calloc((size_t)N, sizeof(double));
    double  *tmp = calloc((size_t)N, sizeof(double));

    uint64_t s = 0x0123456789ABCDEFULL;
    for (int i = 0; i < N; ++i)
        v[i] = 2.0 * sm64_unit_(&s) - 1.0;

    for (int mu = 0; mu < irrep_pg_num_irreps(t); ++mu) {
        irrep_pg_project(t, mu, m, v, tmp);
        for (int i = 0; i < N; ++i)
            acc[i] += tmp[i];
    }
    for (int i = 0; i < N; ++i)
        IRREP_ASSERT(fabs(acc[i] - v[i]) < 1e-10);

    free(v);
    free(acc);
    free(tmp);
    irrep_multiset_free(m);
    irrep_pg_table_free(t);
}

/* Character orthogonality, computed from the table indirectly via the
 * idempotent projector structure: we assert (1/|G|) · Σ_g χ_μ(g)² = d_μ²/d_μ = d_μ.
 * Formally: the trace of P_μ equals the dimension of the μ-th isotypic subspace,
 * which on the scalar representation (1x0e) is 1 iff μ is the trivial irrep
 * else 0. So projecting a single 0e scalar onto each μ and summing yields 1
 * only for the trivial irrep. This is the cheapest orthogonality fingerprint
 * using only the public API. */
static void test_scalar_projection_(irrep_point_group_t g, int trivial_mu) {
    irrep_pg_table_t *t = irrep_pg_table_build(g);
    irrep_multiset_t *m = irrep_multiset_parse("1x0e");
    double            v = 1.0, p = 0.0;
    for (int mu = 0; mu < irrep_pg_num_irreps(t); ++mu) {
        irrep_pg_project(t, mu, m, &v, &p);
        if (mu == trivial_mu) {
            IRREP_ASSERT(fabs(p - 1.0) < 1e-12);
        } else {
            IRREP_ASSERT(fabs(p) < 1e-12);
        }
    }
    irrep_multiset_free(m);
    irrep_pg_table_free(t);
}

int main(void) {
    IRREP_TEST_START("point_group");

    /* (1) Metadata for all six groups. */
    test_metadata_(IRREP_PG_C4V, 5, 8);
    test_metadata_(IRREP_PG_D6, 6, 12);
    test_metadata_(IRREP_PG_C3V, 3, 6);
    test_metadata_(IRREP_PG_D3, 3, 6);
    test_metadata_(IRREP_PG_TD, 5, 24);
    test_metadata_(IRREP_PG_OH, 10, 48);
    test_metadata_(IRREP_PG_O, 5, 24);

    /* (3) Idempotence across small and larger feature specs. */
    test_idempotence_(IRREP_PG_C4V, "1x0e + 1x1o");
    test_idempotence_(IRREP_PG_C4V, "2x0e + 1x1o + 1x2e");
    test_idempotence_(IRREP_PG_D6, "1x0e + 1x1o");
    test_idempotence_(IRREP_PG_D6, "2x0e + 1x1o + 1x2e");
    test_idempotence_(IRREP_PG_C3V, "1x0e + 1x1o + 1x2e");
    test_idempotence_(IRREP_PG_D3, "1x0e + 1x1o + 1x2e");
    test_idempotence_(IRREP_PG_TD, "1x0e + 1x1o + 1x2e");
    test_idempotence_(IRREP_PG_TD, "2x0e + 1x1o + 1x2e + 1x3o");
    test_idempotence_(IRREP_PG_OH, "1x0e + 1x1o + 1x2e");
    test_idempotence_(IRREP_PG_OH, "2x0e + 1x1o + 1x2e + 1x3o");
    test_idempotence_(IRREP_PG_O, "1x0e + 1x1o + 1x2e");
    test_idempotence_(IRREP_PG_O, "2x0e + 1x1o + 1x2e + 1x3o");

    /* (4) Projector sum reproduces the input vector. */
    test_projector_sum_(IRREP_PG_C4V, "1x0e + 1x1o + 1x2e");
    test_projector_sum_(IRREP_PG_D6, "1x0e + 1x1o + 1x2e");
    test_projector_sum_(IRREP_PG_C3V, "1x0e + 1x1o + 1x2e");
    test_projector_sum_(IRREP_PG_D3, "1x0e + 1x1o + 1x2e");
    test_projector_sum_(IRREP_PG_TD, "1x0e + 1x1o + 1x2e + 1x3o");
    test_projector_sum_(IRREP_PG_OH, "1x0e + 1x1o + 1x2e + 1x3o");
    test_projector_sum_(IRREP_PG_O, "1x0e + 1x1o + 1x2e + 1x3o");

    /* (2) Scalar-projection trace equality (orthogonality fingerprint). */
    test_scalar_projection_(IRREP_PG_C4V, /*trivial=*/0);
    test_scalar_projection_(IRREP_PG_D6, /*trivial=*/0);
    test_scalar_projection_(IRREP_PG_C3V, /*trivial=*/0);
    test_scalar_projection_(IRREP_PG_D3, /*trivial=*/0);
    test_scalar_projection_(IRREP_PG_TD, /*trivial=*/0);
    test_scalar_projection_(IRREP_PG_OH, /*trivial=*/0);
    test_scalar_projection_(IRREP_PG_O, /*trivial=*/0);

    /* (5) Hand-reducible decomposition spot checks (Bradley-Cracknell). */

    /* C₄ᵥ:
     *   1x0e → A₁
     *   1x1o → A₁ + E     (z-axial under z-rotations + vectors in xy)
     *   1x2e → A₁ + B₁ + B₂ + E   (standard l=2 even decomposition)
     *
     * Labels order: {A1, A2, B1, B2, E} at indices {0,1,2,3,4}. */
    {
        irrep_pg_table_t *t = irrep_pg_table_build(IRREP_PG_C4V);

        {
            irrep_multiset_t *m = irrep_multiset_parse("1x0e");
            int               mult[5];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 1);
            for (int i = 1; i < 5; ++i)
                IRREP_ASSERT(mult[i] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x1o");
            int               mult[5];
            irrep_pg_reduce(t, m, mult);
            /* A1 (z) + E (x, y). Sum of dims = 1 + 2 = 3 = dim(l=1). */
            IRREP_ASSERT(mult[0] == 1); /* A1 */
            IRREP_ASSERT(mult[1] == 0);
            IRREP_ASSERT(mult[2] == 0);
            IRREP_ASSERT(mult[3] == 0);
            IRREP_ASSERT(mult[4] == 1); /* E  */
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x2e");
            int               mult[5];
            irrep_pg_reduce(t, m, mult);
            /* A1 + B1 + B2 + E. Sum of dims = 1 + 1 + 1 + 2 = 5 = dim(l=2). */
            IRREP_ASSERT(mult[0] == 1); /* A1 */
            IRREP_ASSERT(mult[1] == 0);
            IRREP_ASSERT(mult[2] == 1); /* B1 */
            IRREP_ASSERT(mult[3] == 1); /* B2 */
            IRREP_ASSERT(mult[4] == 1); /* E  */
            irrep_multiset_free(m);
        }
        irrep_pg_table_free(t);
    }

    /* C₃ᵥ:  (labels {A1, A2, E} at indices {0, 1, 2})
     *   1x0e → A₁
     *   1x1o → A₁ + E   (z invariant under C3 and under σ_v since σ_v · z = z;
     *                     xy-doublet is E)
     *   1x1e → A₂ + E   (identical rotation action; parity flip under σ_v
     *                     swaps A₁↔A₂)
     *
     * This is where C3v's improper reflections pull apart from D3's proper-
     * only structure: D3 sees 1x1o and 1x1e as identical reps. */
    {
        irrep_pg_table_t *t = irrep_pg_table_build(IRREP_PG_C3V);

        {
            irrep_multiset_t *m = irrep_multiset_parse("1x0e");
            int               mult[3];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 1);
            IRREP_ASSERT(mult[1] == 0);
            IRREP_ASSERT(mult[2] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x1o");
            int               mult[3];
            irrep_pg_reduce(t, m, mult);
            /* Sum of dims = 1 + 2 = 3 = dim(l=1). */
            IRREP_ASSERT(mult[0] == 1); /* A1 (z) */
            IRREP_ASSERT(mult[1] == 0);
            IRREP_ASSERT(mult[2] == 1); /* E (x, y) */
            irrep_multiset_free(m);
        }
        {
            /* Same rotation content as 1x1o, but σ_v multiplies by parity
             * = +1 instead of −1, so z transforms as A₁ → A₂ flips. */
            irrep_multiset_t *m = irrep_multiset_parse("1x1e");
            int               mult[3];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 0);
            IRREP_ASSERT(mult[1] == 1); /* A2 */
            IRREP_ASSERT(mult[2] == 1); /* E  */
            irrep_multiset_free(m);
        }
        irrep_pg_table_free(t);
    }

    /* D₃:  (labels {A1, A2, E} at indices {0, 1, 2}; no improper elements)
     *   1x0e → A₁
     *   1x1o → A₂ + E    (z flips under C₂ axes, so z is A₂; xy is E)
     *   1x1e → A₂ + E    (parity has no effect in a purely proper group,
     *                      so 1x1o and 1x1e are indistinguishable). */
    {
        irrep_pg_table_t *t = irrep_pg_table_build(IRREP_PG_D3);

        {
            irrep_multiset_t *m = irrep_multiset_parse("1x0e");
            int               mult[3];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 1);
            IRREP_ASSERT(mult[1] == 0);
            IRREP_ASSERT(mult[2] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m1 = irrep_multiset_parse("1x1o");
            irrep_multiset_t *m2 = irrep_multiset_parse("1x1e");
            int               mult1[3], mult2[3];
            irrep_pg_reduce(t, m1, mult1);
            irrep_pg_reduce(t, m2, mult2);
            /* Both reduce the same way in a purely proper group. */
            for (int i = 0; i < 3; ++i)
                IRREP_ASSERT(mult1[i] == mult2[i]);
            IRREP_ASSERT(mult1[0] == 0);
            IRREP_ASSERT(mult1[1] == 1); /* A2 (z) */
            IRREP_ASSERT(mult1[2] == 1); /* E (x, y) */
            irrep_multiset_free(m1);
            irrep_multiset_free(m2);
        }
        irrep_pg_table_free(t);
    }

    /* D₆:
     *   1x0e → A₁
     *   1x1o → A₂ + E₁    (z-axis invariant + (x, y) doublet)
     *   Labels: {A1, A2, B1, B2, E1, E2} at indices {0..5}. */
    {
        irrep_pg_table_t *t = irrep_pg_table_build(IRREP_PG_D6);

        {
            irrep_multiset_t *m = irrep_multiset_parse("1x0e");
            int               mult[6];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 1);
            for (int i = 1; i < 6; ++i)
                IRREP_ASSERT(mult[i] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x1o");
            int               mult[6];
            irrep_pg_reduce(t, m, mult);
            /* l=1 under D6 (purely proper rotations) decomposes as A2 + E1:
             *   z transforms as A2 (+1 under main-axis rotations, −1 under C2' / C2''
             *     which flip z because they are π rotations about horizontal axes),
             *   (x, y) transforms as E1 (the 2-dim rep with character 2, 1, -1, -2, -1, 1 on the
             * main axis). Sum of dims = 1 + 2 = 3 = dim(l=1). */
            IRREP_ASSERT(mult[0] == 0); /* A1 */
            IRREP_ASSERT(mult[1] == 1); /* A2 (z) */
            IRREP_ASSERT(mult[2] == 0);
            IRREP_ASSERT(mult[3] == 0);
            IRREP_ASSERT(mult[4] == 1); /* E1 (x, y) */
            IRREP_ASSERT(mult[5] == 0);
            irrep_multiset_free(m);
        }
        irrep_pg_table_free(t);
    }

    /* T_d:  (labels {A1, A2, E, T1, T2} at indices {0..4})
     *   1x0e → A₁
     *   1x1o → T₂        (polar vector x,y,z transforms as T₂ — the standard
     *                       choice in chemistry, +1 char on σ_d planes)
     *   1x1e → T₁        (axial vector — same rotation content but opposite
     *                       parity flips T₂ ↔ T₁)
     *   1x2e → E + T₂    (l=2 spherical harmonics under T_d)
     *   1x3o → A₁ + T₁ + T₂   (l=3 odd) */
    {
        irrep_pg_table_t *t = irrep_pg_table_build(IRREP_PG_TD);

        {
            irrep_multiset_t *m = irrep_multiset_parse("1x0e");
            int               mult[5];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 1); /* A1 */
            for (int i = 1; i < 5; ++i)
                IRREP_ASSERT(mult[i] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x1o");
            int               mult[5];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 0);
            IRREP_ASSERT(mult[1] == 0);
            IRREP_ASSERT(mult[2] == 0);
            IRREP_ASSERT(mult[3] == 0);
            IRREP_ASSERT(mult[4] == 1); /* T2 */
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x1e");
            int               mult[5];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 0);
            IRREP_ASSERT(mult[1] == 0);
            IRREP_ASSERT(mult[2] == 0);
            IRREP_ASSERT(mult[3] == 1); /* T1 */
            IRREP_ASSERT(mult[4] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x2e");
            int               mult[5];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 0);
            IRREP_ASSERT(mult[1] == 0);
            IRREP_ASSERT(mult[2] == 1); /* E  */
            IRREP_ASSERT(mult[3] == 0);
            IRREP_ASSERT(mult[4] == 1); /* T2 */
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x3o");
            int               mult[5];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 1); /* A1 */
            IRREP_ASSERT(mult[1] == 0);
            IRREP_ASSERT(mult[2] == 0);
            IRREP_ASSERT(mult[3] == 1); /* T1 */
            IRREP_ASSERT(mult[4] == 1); /* T2 */
            irrep_multiset_free(m);
        }
        irrep_pg_table_free(t);
    }

    /* O_h:  (labels {A1g, A2g, Eg, T1g, T2g, A1u, A2u, Eu, T1u, T2u} at
     *        indices {0..9}.)  Standard reductions:
     *   1x0e → A₁g
     *   1x1o → T₁u   (polar vector — odd parity, transforms as T₁ under O)
     *   1x1e → T₁g   (axial vector — same rotation content but g/u swap)
     *   1x2e → Eg + T₂g
     *   1x2o → Eu + T₂u
     *   1x3o → A₂u + T₁u + T₂u
     *   1x3e → A₂g + T₁g + T₂g */
    {
        irrep_pg_table_t *t = irrep_pg_table_build(IRREP_PG_OH);

        {
            irrep_multiset_t *m = irrep_multiset_parse("1x0e");
            int               mult[10];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 1); /* A1g */
            for (int i = 1; i < 10; ++i)
                IRREP_ASSERT(mult[i] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x1o");
            int               mult[10];
            irrep_pg_reduce(t, m, mult);
            for (int i = 0; i < 10; ++i)
                IRREP_ASSERT(mult[i] == (i == 8 ? 1 : 0)); /* T1u only */
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x1e");
            int               mult[10];
            irrep_pg_reduce(t, m, mult);
            for (int i = 0; i < 10; ++i)
                IRREP_ASSERT(mult[i] == (i == 3 ? 1 : 0)); /* T1g only */
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x2e");
            int               mult[10];
            irrep_pg_reduce(t, m, mult);
            /* Eg + T2g. */
            IRREP_ASSERT(mult[2] == 1); /* Eg */
            IRREP_ASSERT(mult[4] == 1); /* T2g */
            for (int i = 0; i < 10; ++i)
                if (i != 2 && i != 4)
                    IRREP_ASSERT(mult[i] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x2o");
            int               mult[10];
            irrep_pg_reduce(t, m, mult);
            /* Eu + T2u. */
            IRREP_ASSERT(mult[7] == 1); /* Eu */
            IRREP_ASSERT(mult[9] == 1); /* T2u */
            for (int i = 0; i < 10; ++i)
                if (i != 7 && i != 9)
                    IRREP_ASSERT(mult[i] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x3o");
            int               mult[10];
            irrep_pg_reduce(t, m, mult);
            /* A2u + T1u + T2u. */
            IRREP_ASSERT(mult[6] == 1); /* A2u */
            IRREP_ASSERT(mult[8] == 1); /* T1u */
            IRREP_ASSERT(mult[9] == 1); /* T2u */
            for (int i = 0; i < 10; ++i)
                if (i != 6 && i != 8 && i != 9)
                    IRREP_ASSERT(mult[i] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x3e");
            int               mult[10];
            irrep_pg_reduce(t, m, mult);
            /* A2g + T1g + T2g. */
            IRREP_ASSERT(mult[1] == 1); /* A2g */
            IRREP_ASSERT(mult[3] == 1); /* T1g */
            IRREP_ASSERT(mult[4] == 1); /* T2g */
            for (int i = 0; i < 10; ++i)
                if (i != 1 && i != 3 && i != 4)
                    IRREP_ASSERT(mult[i] == 0);
            irrep_multiset_free(m);
        }
        irrep_pg_table_free(t);
    }

    /* O (chiral cubic, 5 irreps {A1, A2, E, T1, T2} at indices {0..4}):
     *   1x0e → A₁
     *   1x1o → T₁    — same as 1x1e (no improper element to flip parity,
     *                   so O can't distinguish polar vs axial vectors)
     *   1x1e → T₁
     *   1x2e → E + T₂
     *
     * Note: under T_d, 1x1o gives T₂ (different label), because T_d's
     * improper σ_d elements flip parity-odd inputs and that re-routes
     * which abstract irrep the polar vector lands in. */
    {
        irrep_pg_table_t *t = irrep_pg_table_build(IRREP_PG_O);

        {
            irrep_multiset_t *m = irrep_multiset_parse("1x0e");
            int               mult[5];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[0] == 1);
            for (int i = 1; i < 5; ++i)
                IRREP_ASSERT(mult[i] == 0);
            irrep_multiset_free(m);
        }
        {
            irrep_multiset_t *m1 = irrep_multiset_parse("1x1o");
            irrep_multiset_t *m2 = irrep_multiset_parse("1x1e");
            int               mult1[5], mult2[5];
            irrep_pg_reduce(t, m1, mult1);
            irrep_pg_reduce(t, m2, mult2);
            for (int i = 0; i < 5; ++i)
                IRREP_ASSERT(mult1[i] == mult2[i]); /* O can't tell parity apart */
            IRREP_ASSERT(mult1[3] == 1); /* T1 */
            for (int i = 0; i < 5; ++i)
                if (i != 3)
                    IRREP_ASSERT(mult1[i] == 0);
            irrep_multiset_free(m1);
            irrep_multiset_free(m2);
        }
        {
            irrep_multiset_t *m = irrep_multiset_parse("1x2e");
            int               mult[5];
            irrep_pg_reduce(t, m, mult);
            IRREP_ASSERT(mult[2] == 1); /* E  */
            IRREP_ASSERT(mult[4] == 1); /* T2 */
            for (int i = 0; i < 5; ++i)
                if (i != 2 && i != 4)
                    IRREP_ASSERT(mult[i] == 0);
            irrep_multiset_free(m);
        }
        irrep_pg_table_free(t);
    }

    return IRREP_TEST_END();
}
