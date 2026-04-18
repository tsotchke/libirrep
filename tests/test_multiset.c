/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/multiset.h>

int main(void) {
    IRREP_TEST_START("multiset");

    /* -------- empty multiset via new / direct-sum of empty -------- */
    {
        irrep_multiset_t *m = irrep_multiset_new(4);
        IRREP_ASSERT(m != NULL);
        IRREP_ASSERT(m->capacity  == 4);
        IRREP_ASSERT(m->num_terms == 0);
        IRREP_ASSERT(m->total_dim == 0);
        irrep_multiset_free(m);
    }

    /* -------- append grows capacity and tracks total_dim -------- */
    {
        irrep_multiset_t *m = irrep_multiset_new(0);
        IRREP_ASSERT(m != NULL);
        irrep_status_t s = irrep_multiset_append(m,
            (irrep_label_t){ .l = 0, .parity = IRREP_EVEN }, 1);
        IRREP_ASSERT(s == IRREP_OK);
        IRREP_ASSERT(m->num_terms == 1);
        IRREP_ASSERT(m->total_dim == 1);
        IRREP_ASSERT(m->capacity  >= 1);

        irrep_multiset_append(m,
            (irrep_label_t){ .l = 1, .parity = IRREP_ODD }, 2);
        IRREP_ASSERT(m->num_terms == 2);
        IRREP_ASSERT(m->total_dim == 1 + 2 * 3);     /* 0e + 2·1o = 1 + 6 */

        irrep_multiset_append(m,
            (irrep_label_t){ .l = 2, .parity = IRREP_EVEN }, 1);
        IRREP_ASSERT(m->total_dim == 1 + 6 + 5);

        /* reject invalid args */
        IRREP_ASSERT(irrep_multiset_append(m,
            (irrep_label_t){ .l = -1, .parity = IRREP_EVEN }, 1)
            == IRREP_ERR_INVALID_ARG);
        IRREP_ASSERT(irrep_multiset_append(m,
            (irrep_label_t){ .l = 0, .parity = 7 }, 1)
            == IRREP_ERR_INVALID_ARG);
        IRREP_ASSERT(irrep_multiset_append(m,
            (irrep_label_t){ .l = 0, .parity = IRREP_EVEN }, 0)
            == IRREP_ERR_INVALID_ARG);

        irrep_multiset_free(m);
    }

    /* -------- parse simple case -------- */
    {
        irrep_multiset_t *m = irrep_multiset_parse("1x0e + 2x1o + 1x2e");
        IRREP_ASSERT(m != NULL);
        IRREP_ASSERT(m->num_terms == 3);
        IRREP_ASSERT(m->labels[0].l == 0 && m->labels[0].parity == IRREP_EVEN);
        IRREP_ASSERT(m->multiplicities[0] == 1);
        IRREP_ASSERT(m->labels[1].l == 1 && m->labels[1].parity == IRREP_ODD);
        IRREP_ASSERT(m->multiplicities[1] == 2);
        IRREP_ASSERT(m->labels[2].l == 2 && m->labels[2].parity == IRREP_EVEN);
        IRREP_ASSERT(m->multiplicities[2] == 1);
        IRREP_ASSERT(m->total_dim == 1 + 2 * 3 + 5);
        irrep_multiset_free(m);
    }

    /* -------- parse tolerates whitespace, no '+' between single term -------- */
    {
        irrep_multiset_t *m = irrep_multiset_parse("   3 x 4 e  ");
        IRREP_ASSERT(m != NULL);
        IRREP_ASSERT(m->num_terms == 1);
        IRREP_ASSERT(m->labels[0].l == 4 && m->labels[0].parity == IRREP_EVEN);
        IRREP_ASSERT(m->multiplicities[0] == 3);
        IRREP_ASSERT(m->total_dim == 3 * 9);
        irrep_multiset_free(m);
    }

    /* -------- parse empty string → empty multiset -------- */
    {
        irrep_multiset_t *m = irrep_multiset_parse("");
        IRREP_ASSERT(m != NULL);
        IRREP_ASSERT(m->num_terms == 0);
        irrep_multiset_free(m);
    }

    /* -------- parse error cases -------- */
    IRREP_ASSERT(irrep_multiset_parse(NULL)          == NULL);
    IRREP_ASSERT(irrep_multiset_parse("0x0e")        == NULL);  /* zero mult */
    IRREP_ASSERT(irrep_multiset_parse("1x0q")        == NULL);  /* bad parity */
    IRREP_ASSERT(irrep_multiset_parse("1x-1e")       == NULL);  /* negative l */
    IRREP_ASSERT(irrep_multiset_parse("1y0e")        == NULL);  /* missing x */
    IRREP_ASSERT(irrep_multiset_parse("1x0e + ")     == NULL);  /* trailing + */

    /* -------- format round-trip -------- */
    {
        const char *src = "1x0e + 2x1o + 1x2e";
        irrep_multiset_t *m = irrep_multiset_parse(src);
        IRREP_ASSERT(m != NULL);
        char buf[64];
        int n = irrep_multiset_format(m, buf, sizeof(buf));
        IRREP_ASSERT(n > 0);
        IRREP_ASSERT(strcmp(buf, src) == 0);
        irrep_multiset_free(m);
    }

    /* -------- format returns required length when buf is NULL or too small -------- */
    {
        irrep_multiset_t *m = irrep_multiset_parse("1x0e + 2x1o");
        char small[4];
        int  full_len = irrep_multiset_format(m, NULL, 0);
        int  short_len = irrep_multiset_format(m, small, sizeof(small));
        IRREP_ASSERT(full_len == short_len);   /* required length */
        IRREP_ASSERT(small[sizeof(small) - 1] == '\0');
        irrep_multiset_free(m);
    }

    /* -------- simplify merges + sorts -------- */
    {
        irrep_multiset_t *m = irrep_multiset_parse("2x1o + 1x0e + 3x1o + 1x2e + 2x0e");
        IRREP_ASSERT(m != NULL);
        int pre_dim = m->total_dim;
        irrep_multiset_simplify(m);
        /* Expected canonical form: 3x0e + 5x1o + 1x2e; dim preserved */
        IRREP_ASSERT(m->num_terms == 3);
        IRREP_ASSERT(m->labels[0].l == 0 && m->labels[0].parity == IRREP_EVEN);
        IRREP_ASSERT(m->multiplicities[0] == 3);
        IRREP_ASSERT(m->labels[1].l == 1 && m->labels[1].parity == IRREP_ODD);
        IRREP_ASSERT(m->multiplicities[1] == 5);
        IRREP_ASSERT(m->labels[2].l == 2 && m->labels[2].parity == IRREP_EVEN);
        IRREP_ASSERT(m->multiplicities[2] == 1);
        IRREP_ASSERT(m->total_dim == pre_dim);

        /* Idempotent */
        irrep_multiset_simplify(m);
        IRREP_ASSERT(m->num_terms == 3);
        IRREP_ASSERT(m->total_dim == pre_dim);

        irrep_multiset_free(m);
    }

    /* -------- direct sum -------- */
    {
        irrep_multiset_t *a = irrep_multiset_parse("1x0e + 2x1o");
        irrep_multiset_t *b = irrep_multiset_parse("3x0e + 1x2e");
        irrep_multiset_t *c = irrep_multiset_direct_sum(a, b);
        IRREP_ASSERT(c != NULL);
        /* After simplify: 4x0e + 2x1o + 1x2e; total = 4 + 6 + 5 = 15 */
        IRREP_ASSERT(c->num_terms == 3);
        IRREP_ASSERT(c->total_dim == 4 + 6 + 5);
        IRREP_ASSERT(c->labels[0].l == 0 && c->multiplicities[0] == 4);
        IRREP_ASSERT(c->labels[1].l == 1 && c->multiplicities[1] == 2);
        IRREP_ASSERT(c->labels[2].l == 2 && c->multiplicities[2] == 1);
        irrep_multiset_free(a);
        irrep_multiset_free(b);
        irrep_multiset_free(c);
    }

    /* -------- block_offset running sum -------- */
    {
        irrep_multiset_t *m = irrep_multiset_parse("1x0e + 2x1o + 1x2e");
        IRREP_ASSERT(irrep_multiset_block_offset(m, 0) == 0);
        IRREP_ASSERT(irrep_multiset_block_offset(m, 1) == 1);        /* + 1·(1) */
        IRREP_ASSERT(irrep_multiset_block_offset(m, 2) == 1 + 6);    /* + 2·(3) */
        IRREP_ASSERT(irrep_multiset_block_offset(m, 3) == 1 + 6 + 5);
        IRREP_ASSERT(irrep_multiset_dim(m)              == 1 + 6 + 5);
        irrep_multiset_free(m);
    }

    return IRREP_TEST_END();
}
