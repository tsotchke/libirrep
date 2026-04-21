/* SPDX-License-Identifier: MIT */
/* Tests for parity (inversion) handling in the O(3) extension of the
 * library: scalar `parity()` of an irrep label, parity products, and
 * parity-aware path filtering used when building tensor products in O(3).
 *
 * Coverage:
 *   - `parity()` and `parity_product()` on scalar labels (even × odd = odd).
 *   - `parity_filter_paths` keeps only paths with p_a · p_b = p_c and
 *     drops parity-violating ones.
 *   - Out-of-range indices are silently dropped rather than asserting.
 */
#include "harness.h"
#include <irrep/parity.h>
#include <irrep/multiset.h>

int main(void) {
    IRREP_TEST_START("parity");

    /* -------- parity() and parity_product() on scalar labels -------- */
    {
        irrep_label_t e = {.l = 1, .parity = IRREP_EVEN};
        irrep_label_t o = {.l = 1, .parity = IRREP_ODD};
        IRREP_ASSERT(irrep_parity(e) == +1);
        IRREP_ASSERT(irrep_parity(o) == -1);
        IRREP_ASSERT(irrep_parity_product(e, e) == +1);
        IRREP_ASSERT(irrep_parity_product(e, o) == -1);
        IRREP_ASSERT(irrep_parity_product(o, e) == -1);
        IRREP_ASSERT(irrep_parity_product(o, o) == +1);
    }

    /* -------- parity_filter_paths keeps valid, drops parity-violating paths -------- */
    {
        irrep_multiset_t *a = irrep_multiset_parse("1x1o + 1x1e"); /* [0]=1o, [1]=1e */
        irrep_multiset_t *b = irrep_multiset_parse("1x1o + 1x1e");
        irrep_multiset_t *c = irrep_multiset_parse("1x0e + 1x0o"); /* [0]=0e, [1]=0o */

        /* Paths triplets (i_a, i_b, i_c). Enumerate all 2×2×2 = 8, expect 4 kept:
         *   (o, o, e) = +1  ✓
         *   (o, e, o) = -1  ✓
         *   (e, o, o) = -1  ✓
         *   (e, e, e) = +1  ✓
         *   (o, o, o) = +1  ✗   (want parity_c = -1 but pa·pb = +1)
         *   (o, e, e) = -1  ✗
         *   (e, o, e) = -1  ✗
         *   (e, e, o) = +1  ✗ */
        int paths[8 * 3] = {
            0,
            0,
            0, /* o o · e → 1·1=1 ≠ pc=+1 wait */
            /* Let me be careful: a[0]=1o so parity=-1, a[1]=1e so parity=+1, etc. */
            0,
            0,
            1, /* o·o·(0o): (-1)(-1)=+1, pc=-1 → fail */
            0,
            1,
            0, /* o·e·(0e): (-1)(+1)=-1, pc=+1 → fail */
            0,
            1,
            1, /* o·e·(0o): (-1)(+1)=-1, pc=-1 → ok */
            1,
            0,
            0, /* e·o·(0e): (+1)(-1)=-1, pc=+1 → fail */
            1,
            0,
            1, /* e·o·(0o): -1 == -1 → ok */
            1,
            1,
            0, /* e·e·(0e): +1 == +1 → ok */
            1,
            1,
            1, /* e·e·(0o): +1 != -1 → fail */
        };
        /* Row 0 duplicated by accident above — let me count carefully. The
         * 8 valid triples are the 8 distinct (ia, ib, ic) tuples; I'll just
         * build them explicitly. */
        int all[8 * 3] = {
            0, 0, 0, /* o o 0e: +1 vs +1 ok */
            0, 0, 1, /* o o 0o: +1 vs -1 no */
            0, 1, 0, /* o e 0e: -1 vs +1 no */
            0, 1, 1, /* o e 0o: -1 vs -1 ok */
            1, 0, 0, /* e o 0e: -1 vs +1 no */
            1, 0, 1, /* e o 0o: -1 vs -1 ok */
            1, 1, 0, /* e e 0e: +1 vs +1 ok */
            1, 1, 1, /* e e 0o: +1 vs -1 no */
        };
        (void)paths;
        int kept = irrep_parity_filter_paths(a, b, c, all, 8);
        IRREP_ASSERT(kept == 4);
        /* The kept triplets should be (0,0,0), (0,1,1), (1,0,1), (1,1,0). */
        IRREP_ASSERT(all[0] == 0 && all[1] == 0 && all[2] == 0);
        IRREP_ASSERT(all[3] == 0 && all[4] == 1 && all[5] == 1);
        IRREP_ASSERT(all[6] == 1 && all[7] == 0 && all[8] == 1);
        IRREP_ASSERT(all[9] == 1 && all[10] == 1 && all[11] == 0);

        irrep_multiset_free(a);
        irrep_multiset_free(b);
        irrep_multiset_free(c);
    }

    /* -------- out-of-range indices are silently dropped -------- */
    {
        irrep_multiset_t *a = irrep_multiset_parse("1x0e");
        irrep_multiset_t *b = irrep_multiset_parse("1x0e");
        irrep_multiset_t *c = irrep_multiset_parse("1x0e");
        int               paths[2 * 3] = {0, 0, 0, 5, 5, 5};
        int               kept = irrep_parity_filter_paths(a, b, c, paths, 2);
        IRREP_ASSERT(kept == 1);

        irrep_multiset_free(a);
        irrep_multiset_free(b);
        irrep_multiset_free(c);
    }

    return IRREP_TEST_END();
}
