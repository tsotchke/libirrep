/* SPDX-License-Identifier: MIT */
/* M6: irrep-multiset algebra.
 *
 * Parser grammar (matches e3nn canonical form):
 *
 *     multiset := term ( '+' term )*      whitespace permitted anywhere
 *     term     := mult 'x' l parity
 *     mult     := positive integer
 *     l        := non-negative integer  (≤ IRREP_L_MAX)
 *     parity   := 'e' | 'o'
 *
 * Canonical simplified form: sorted by (l ascending, parity even-before-odd),
 * with like terms merged. Empty string parses to an empty multiset.
 */

#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/multiset.h>

#define UNUSED(x) ((void)(x))

/* forward: irrep_set_error_ lives in src/error.c */
extern void irrep_set_error_(const char *fmt, ...);

/* -------------------------------------------------------------------------- *
 * Construction / destruction                                                 *
 * -------------------------------------------------------------------------- */

irrep_multiset_t *irrep_multiset_new(int capacity) {
    if (capacity < 0)
        return NULL;
    irrep_multiset_t *m = calloc(1, sizeof(*m));
    if (!m)
        return NULL;
    if (capacity > 0) {
        m->labels = calloc((size_t)capacity, sizeof(*m->labels));
        m->multiplicities = calloc((size_t)capacity, sizeof(*m->multiplicities));
        if (!m->labels || !m->multiplicities) {
            free(m->labels);
            free(m->multiplicities);
            free(m);
            return NULL;
        }
    }
    m->capacity = capacity;
    m->num_terms = 0;
    m->total_dim = 0;
    return m;
}

void irrep_multiset_free(irrep_multiset_t *m) {
    if (!m)
        return;
    free(m->labels);
    free(m->multiplicities);
    free(m);
}

/* -------------------------------------------------------------------------- *
 * Append with amortized growth                                               *
 * -------------------------------------------------------------------------- */

irrep_status_t irrep_multiset_append(irrep_multiset_t *m, irrep_label_t label, int multiplicity) {
    if (!m)
        return IRREP_ERR_INVALID_ARG;
    if (multiplicity <= 0)
        return IRREP_ERR_INVALID_ARG;
    if (label.l < 0)
        return IRREP_ERR_INVALID_ARG;
    if (label.parity != IRREP_EVEN && label.parity != IRREP_ODD)
        return IRREP_ERR_INVALID_ARG;

    if (m->num_terms >= m->capacity) {
        int            new_cap = m->capacity > 0 ? m->capacity * 2 : 4;
        irrep_label_t *nl = realloc(m->labels, (size_t)new_cap * sizeof(*m->labels));
        int *nm = realloc(m->multiplicities, (size_t)new_cap * sizeof(*m->multiplicities));
        if (!nl || !nm) {
            /* best-effort: partial realloc may have freed old buffer */
            if (nl)
                m->labels = nl;
            if (nm)
                m->multiplicities = nm;
            return IRREP_ERR_OUT_OF_MEMORY;
        }
        m->labels = nl;
        m->multiplicities = nm;
        m->capacity = new_cap;
    }
    m->labels[m->num_terms] = label;
    m->multiplicities[m->num_terms] = multiplicity;
    m->num_terms++;
    m->total_dim += multiplicity * (2 * label.l + 1);
    return IRREP_OK;
}

/* -------------------------------------------------------------------------- *
 * Parser                                                                     *
 * -------------------------------------------------------------------------- */

static void skip_ws_(const char **p) {
    while (**p && isspace((unsigned char)**p))
        (*p)++;
}

irrep_multiset_t *irrep_multiset_parse(const char *spec) {
    if (!spec) {
        irrep_set_error_("irrep_multiset_parse: NULL spec");
        return NULL;
    }
    irrep_multiset_t *m = irrep_multiset_new(4);
    if (!m)
        return NULL;

    const char *p = spec;
    skip_ws_(&p);
    if (!*p)
        return m; /* empty spec → empty multiset */

    while (*p) {
        /* multiplicity */
        char *end;
        long  mult = strtol(p, &end, 10);
        if (end == p || mult <= 0 || mult > INT_MAX) {
            irrep_set_error_("irrep_multiset_parse: invalid multiplicity near '%s'", p);
            irrep_multiset_free(m);
            return NULL;
        }
        p = end;

        skip_ws_(&p);
        if (*p != 'x') {
            irrep_set_error_("irrep_multiset_parse: expected 'x' near '%s'", p);
            irrep_multiset_free(m);
            return NULL;
        }
        p++;
        skip_ws_(&p);

        /* l */
        long l = strtol(p, &end, 10);
        if (end == p || l < 0 || l > IRREP_L_MAX) {
            irrep_set_error_("irrep_multiset_parse: invalid l near '%s'", p);
            irrep_multiset_free(m);
            return NULL;
        }
        p = end;
        skip_ws_(&p);

        /* parity */
        int parity;
        if (*p == 'e')
            parity = IRREP_EVEN;
        else if (*p == 'o')
            parity = IRREP_ODD;
        else {
            irrep_set_error_("irrep_multiset_parse: expected 'e' or 'o' near '%s'", p);
            irrep_multiset_free(m);
            return NULL;
        }
        p++;

        irrep_label_t lbl = {.l = (int)l, .parity = parity};
        if (irrep_multiset_append(m, lbl, (int)mult) != IRREP_OK) {
            irrep_multiset_free(m);
            return NULL;
        }

        skip_ws_(&p);
        if (*p == '+') {
            p++;
            skip_ws_(&p);
            if (!*p) {
                irrep_set_error_("irrep_multiset_parse: trailing '+' with no following term");
                irrep_multiset_free(m);
                return NULL;
            }
            continue;
        }
        if (*p == '\0')
            break;

        irrep_set_error_("irrep_multiset_parse: expected '+' or end near '%s'", p);
        irrep_multiset_free(m);
        return NULL;
    }
    return m;
}

/* -------------------------------------------------------------------------- *
 * Formatter                                                                  *
 * -------------------------------------------------------------------------- */

int irrep_multiset_format(const irrep_multiset_t *m, char *buf, size_t buflen) {
    if (!m || m->num_terms == 0) {
        if (buf && buflen > 0)
            buf[0] = '\0';
        return 0;
    }
    size_t written = 0;
    size_t total = 0;
    for (int i = 0; i < m->num_terms; ++i) {
        char chunk[64];
        int  n =
            snprintf(chunk, sizeof(chunk), "%s%dx%d%c", i == 0 ? "" : " + ", m->multiplicities[i],
                     m->labels[i].l, m->labels[i].parity == IRREP_ODD ? 'o' : 'e');
        if (n < 0)
            return -1;
        if (buf && written + 1 < buflen) {
            size_t copy = (size_t)n < buflen - 1 - written ? (size_t)n : buflen - 1 - written;
            memcpy(buf + written, chunk, copy);
            written += copy;
        }
        total += (size_t)n;
    }
    if (buf && buflen > 0)
        buf[written < buflen ? written : buflen - 1] = '\0';
    return (int)total;
}

/* -------------------------------------------------------------------------- *
 * Simplify: canonical sort + merge                                           *
 * -------------------------------------------------------------------------- */

void irrep_multiset_simplify(irrep_multiset_t *m) {
    if (!m || m->num_terms <= 0)
        return;

    int n = m->num_terms;
    /* Bubble sort — n is small in practice. Stable enough, and we're about to
     * merge anyway. Order: l ascending, then parity even (+1) before odd (-1). */
    for (int i = 0; i < n - 1; ++i) {
        for (int k = 0; k < n - 1 - i; ++k) {
            int lk = m->labels[k].l;
            int lkp = m->labels[k + 1].l;
            int pk = m->labels[k].parity;
            int pkp = m->labels[k + 1].parity;
            int swap = 0;
            if (lk > lkp)
                swap = 1;
            else if (lk == lkp && pk < pkp)
                swap = 1; /* −1 before +1 → swap */
            if (swap) {
                irrep_label_t lt = m->labels[k];
                int           mt = m->multiplicities[k];
                m->labels[k] = m->labels[k + 1];
                m->multiplicities[k] = m->multiplicities[k + 1];
                m->labels[k + 1] = lt;
                m->multiplicities[k + 1] = mt;
            }
        }
    }

    int write = 0;
    for (int r = 0; r < n; ++r) {
        if (m->multiplicities[r] == 0)
            continue;
        if (write > 0 && m->labels[write - 1].l == m->labels[r].l &&
            m->labels[write - 1].parity == m->labels[r].parity) {
            m->multiplicities[write - 1] += m->multiplicities[r];
        } else {
            m->labels[write] = m->labels[r];
            m->multiplicities[write] = m->multiplicities[r];
            write++;
        }
    }
    m->num_terms = write;

    int total = 0;
    for (int i = 0; i < write; ++i) {
        total += m->multiplicities[i] * (2 * m->labels[i].l + 1);
    }
    m->total_dim = total;
}

/* -------------------------------------------------------------------------- *
 * Direct sum                                                                 *
 * -------------------------------------------------------------------------- */

irrep_multiset_t *irrep_multiset_direct_sum(const irrep_multiset_t *m1,
                                            const irrep_multiset_t *m2) {
    int               cap = (m1 ? m1->num_terms : 0) + (m2 ? m2->num_terms : 0);
    irrep_multiset_t *out = irrep_multiset_new(cap > 0 ? cap : 1);
    if (!out)
        return NULL;
    if (m1) {
        for (int i = 0; i < m1->num_terms; ++i) {
            if (irrep_multiset_append(out, m1->labels[i], m1->multiplicities[i]) != IRREP_OK) {
                irrep_multiset_free(out);
                return NULL;
            }
        }
    }
    if (m2) {
        for (int i = 0; i < m2->num_terms; ++i) {
            if (irrep_multiset_append(out, m2->labels[i], m2->multiplicities[i]) != IRREP_OK) {
                irrep_multiset_free(out);
                return NULL;
            }
        }
    }
    irrep_multiset_simplify(out);
    return out;
}

/* -------------------------------------------------------------------------- *
 * Accessors                                                                  *
 * -------------------------------------------------------------------------- */

int irrep_multiset_dim(const irrep_multiset_t *m) {
    return m ? m->total_dim : 0;
}

int irrep_multiset_block_offset(const irrep_multiset_t *m, int term_idx) {
    if (!m || term_idx <= 0)
        return 0;
    if (term_idx > m->num_terms)
        term_idx = m->num_terms;
    int offset = 0;
    for (int i = 0; i < term_idx; ++i) {
        offset += m->multiplicities[i] * (2 * m->labels[i].l + 1);
    }
    return offset;
}
