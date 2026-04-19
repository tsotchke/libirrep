/* SPDX-License-Identifier: MIT */
/* Tests for the doubled-integer multiset type. */

#include "harness.h"
#include <irrep/multiset_2j.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    IRREP_TEST_START("multiset_2j");

    /* ---- Parser: integer-only form backwards-compatible with multiset.h. ---- */
    {
        irrep_multiset_2j_t *m = irrep_multiset_2j_parse("1x0e + 2x1o + 1x2e");
        IRREP_ASSERT(m != NULL);
        IRREP_ASSERT(m->num_terms == 3);
        IRREP_ASSERT(m->labels[0].two_j == 0 && m->labels[0].parity == +1);
        IRREP_ASSERT(m->labels[1].two_j == 2 && m->labels[1].parity == -1);
        IRREP_ASSERT(m->labels[2].two_j == 4 && m->labels[2].parity == +1);
        IRREP_ASSERT(m->multiplicities[0] == 1);
        IRREP_ASSERT(m->multiplicities[1] == 2);
        IRREP_ASSERT(m->multiplicities[2] == 1);
        /* total_dim = 1·(0+1) + 2·(2+1) + 1·(4+1) = 1 + 6 + 5 = 12 */
        IRREP_ASSERT(m->total_dim == 12);
        IRREP_ASSERT(irrep_multiset_2j_has_half_integer(m) == 0);
        IRREP_ASSERT(irrep_time_reversal_square_sign_2j(m) == +1);
        irrep_multiset_2j_free(m);
    }

    /* ---- Parser: half-integer form. ---- */
    {
        irrep_multiset_2j_t *m = irrep_multiset_2j_parse("1x1/2o + 2x3/2e");
        IRREP_ASSERT(m != NULL);
        IRREP_ASSERT(m->num_terms == 2);
        /* two_j = 1 (spin 1/2), parity odd */
        IRREP_ASSERT(m->labels[0].two_j == 1 && m->labels[0].parity == -1);
        /* two_j = 3 (spin 3/2), parity even */
        IRREP_ASSERT(m->labels[1].two_j == 3 && m->labels[1].parity == +1);
        /* total_dim = 1·(1+1) + 2·(3+1) = 2 + 8 = 10 */
        IRREP_ASSERT(m->total_dim == 10);
        IRREP_ASSERT(irrep_multiset_2j_has_half_integer(m) == 1);
        IRREP_ASSERT(irrep_time_reversal_square_sign_2j(m) == -1);
        irrep_multiset_2j_free(m);
    }

    /* ---- Mixed integer + half-integer. ---- */
    {
        irrep_multiset_2j_t *m = irrep_multiset_2j_parse("1x0e + 1x1/2o + 1x1e");
        IRREP_ASSERT(m != NULL);
        IRREP_ASSERT(m->num_terms == 3);
        /* total_dim = 1·1 + 1·2 + 1·3 = 6 */
        IRREP_ASSERT(m->total_dim == 6);
        IRREP_ASSERT(irrep_multiset_2j_has_half_integer(m) == 1);
        /* Any half-integer block → Kramers sign −1 */
        IRREP_ASSERT(irrep_time_reversal_square_sign_2j(m) == -1);
        irrep_multiset_2j_free(m);
    }

    /* ---- Empty spec → empty multiset. ---- */
    {
        irrep_multiset_2j_t *m = irrep_multiset_2j_parse("");
        IRREP_ASSERT(m != NULL);
        IRREP_ASSERT(m->num_terms == 0);
        IRREP_ASSERT(m->total_dim == 0);
        IRREP_ASSERT(irrep_time_reversal_square_sign_2j(m) == 0);
        irrep_multiset_2j_free(m);
    }

    /* ---- Malformed inputs return NULL. ---- */
    IRREP_ASSERT(irrep_multiset_2j_parse(NULL) == NULL);
    IRREP_ASSERT(irrep_multiset_2j_parse("foo") == NULL);
    IRREP_ASSERT(irrep_multiset_2j_parse("-1x0e") == NULL);
    IRREP_ASSERT(irrep_multiset_2j_parse("1x2/2e") == NULL);    /* even numerator */
    IRREP_ASSERT(irrep_multiset_2j_parse("1x0/2e") == NULL);    /* zero numerator */
    IRREP_ASSERT(irrep_multiset_2j_parse("1x1/3e") == NULL);    /* denominator ≠ 2 */
    IRREP_ASSERT(irrep_multiset_2j_parse("1x0e + ") == NULL);   /* trailing + */
    IRREP_ASSERT(irrep_multiset_2j_parse("1x0x") == NULL);      /* bad parity */

    /* ---- Format round trip. ---- */
    {
        irrep_multiset_2j_t *m = irrep_multiset_2j_parse("2x0e + 1x1/2o + 3x1e");
        IRREP_ASSERT(m != NULL);
        char buf[128];
        int n = irrep_multiset_2j_format(m, buf, sizeof(buf));
        IRREP_ASSERT(n > 0);
        IRREP_ASSERT(strstr(buf, "2x0e") != NULL);
        IRREP_ASSERT(strstr(buf, "1x1/2o") != NULL);
        IRREP_ASSERT(strstr(buf, "3x1e") != NULL);
        irrep_multiset_2j_free(m);
    }

    /* ---- append directly with a half-integer label. ---- */
    {
        irrep_multiset_2j_t *m = irrep_multiset_2j_new(4);
        IRREP_ASSERT(m != NULL);
        IRREP_ASSERT(irrep_multiset_2j_append(m, (irrep_label_2j_t){3, +1}, 2) == IRREP_OK);
        IRREP_ASSERT(m->num_terms == 1);
        IRREP_ASSERT(m->labels[0].two_j == 3);
        IRREP_ASSERT(m->total_dim == 2 * 4);
        IRREP_ASSERT(irrep_time_reversal_square_sign_2j(m) == -1);
        /* Invalid append: negative two_j */
        IRREP_ASSERT(irrep_multiset_2j_append(m, (irrep_label_2j_t){-1, +1}, 1) == IRREP_ERR_INVALID_ARG);
        /* Invalid append: zero multiplicity */
        IRREP_ASSERT(irrep_multiset_2j_append(m, (irrep_label_2j_t){0, +1}, 0) == IRREP_ERR_INVALID_ARG);
        irrep_multiset_2j_free(m);
    }

    return IRREP_TEST_END();
}
