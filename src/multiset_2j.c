/* SPDX-License-Identifier: MIT */
/* Doubled-integer multiset implementation. See include/irrep/multiset_2j.h.
 *
 * Parser is structurally parallel to `src/multiset.c` but extended with an
 * optional `/2` suffix on the label for half-integer content. All storage
 * uses `two_j` (doubled) consistently — spin-½ is `two_j = 1`, block
 * dimension is `two_j + 1`.
 *
 * Kept minimal and additive. The integer-only `irrep_multiset_t` remains
 * the canonical type for the existing ML-pipeline APIs; this type is
 * specifically for consumers carrying half-integer spin content (spinor
 * wavefunctions, magnetic-moment representations with spin-orbit
 * coupling, Kramers-paired states). */

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/multiset_2j.h>

extern void irrep_set_error_(const char *fmt, ...);

/* -------------------------------------------------------------------------- *
 * Construction                                                               *
 * -------------------------------------------------------------------------- */

irrep_multiset_2j_t *irrep_multiset_2j_new(int capacity) {
    if (capacity < 0)
        return NULL;
    irrep_multiset_2j_t *m = calloc(1, sizeof(*m));
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
    return m;
}

void irrep_multiset_2j_free(irrep_multiset_2j_t *m) {
    if (!m)
        return;
    free(m->labels);
    free(m->multiplicities);
    free(m);
}

static int grow_(irrep_multiset_2j_t *m, int needed) {
    if (needed <= m->capacity)
        return 0;
    int new_cap = m->capacity ? m->capacity * 2 : 4;
    while (new_cap < needed)
        new_cap *= 2;
    irrep_label_2j_t *nl = realloc(m->labels, (size_t)new_cap * sizeof(*m->labels));
    int              *nm = realloc(m->multiplicities, (size_t)new_cap * sizeof(*m->multiplicities));
    if (!nl || !nm) {
        free(nl);
        free(nm);
        return -1;
    }
    m->labels = nl;
    m->multiplicities = nm;
    m->capacity = new_cap;
    return 0;
}

irrep_status_t irrep_multiset_2j_append(irrep_multiset_2j_t *m, irrep_label_2j_t label,
                                        int multiplicity) {
    if (!m)
        return IRREP_ERR_INVALID_ARG;
    if (label.two_j < 0)
        return IRREP_ERR_INVALID_ARG;
    if (label.parity != +1 && label.parity != -1)
        return IRREP_ERR_INVALID_ARG;
    if (multiplicity <= 0)
        return IRREP_ERR_INVALID_ARG;

    /* Guard the total_dim accumulator against `int` overflow. The block
     * contribution is `multiplicity * (two_j + 1)`; the running sum
     * `m->total_dim + block` must stay within INT_MAX so that every
     * downstream consumer indexing a feature buffer sees a non-negative
     * dimension. At ML-pipeline scales this is never tight, but the
     * guard is cheap and the silent wrap is a long-tail correctness
     * hazard we do not want to ship. */
    long block_dim = (long)label.two_j + 1L;
    if (multiplicity > INT_MAX / block_dim)
        return IRREP_ERR_INVALID_ARG;
    long added = (long)multiplicity * block_dim;
    if (m->total_dim > INT_MAX - added)
        return IRREP_ERR_INVALID_ARG;

    if (grow_(m, m->num_terms + 1) != 0)
        return IRREP_ERR_OUT_OF_MEMORY;
    m->labels[m->num_terms] = label;
    m->multiplicities[m->num_terms] = multiplicity;
    m->total_dim += (int)added;
    ++m->num_terms;
    return IRREP_OK;
}

/* -------------------------------------------------------------------------- *
 * Parser                                                                     *
 * -------------------------------------------------------------------------- */

static void skip_ws(const char **p) {
    while (**p == ' ' || **p == '\t')
        ++(*p);
}

static int parse_uint_(const char **pp, long *out) {
    const char *p = *pp;
    if (!isdigit((unsigned char)*p))
        return -1;
    char *end = NULL;
    errno = 0;
    long v = strtol(p, &end, 10);
    if (end == p)
        return -1;
    *pp = end;
    *out = v;
    return 0;
}

irrep_multiset_2j_t *irrep_multiset_2j_parse(const char *spec) {
    if (!spec) {
        irrep_set_error_("irrep_multiset_2j_parse: NULL input");
        return NULL;
    }
    irrep_multiset_2j_t *m = irrep_multiset_2j_new(0);
    if (!m)
        return NULL;

    const char *p = spec;
    skip_ws(&p);
    if (*p == '\0')
        return m; /* empty spec → empty multiset */

    for (;;) {
        skip_ws(&p);

        long mult;
        if (parse_uint_(&p, &mult) != 0 || mult <= 0 || mult > INT_MAX) {
            irrep_set_error_("irrep_multiset_2j_parse: expected positive "
                             "multiplicity at column %ld",
                             (long)(p - spec + 1));
            irrep_multiset_2j_free(m);
            return NULL;
        }
        if (*p != 'x') {
            irrep_set_error_("irrep_multiset_2j_parse: expected 'x' after "
                             "multiplicity at column %ld",
                             (long)(p - spec + 1));
            irrep_multiset_2j_free(m);
            return NULL;
        }
        ++p;

        long num;
        if (parse_uint_(&p, &num) != 0 || num < 0 || num > INT_MAX) {
            irrep_set_error_("irrep_multiset_2j_parse: expected non-negative "
                             "integer for l or 2j numerator at column %ld",
                             (long)(p - spec + 1));
            irrep_multiset_2j_free(m);
            return NULL;
        }
        int two_j;
        if (*p == '/') {
            /* Half-integer form: `N/2`, N must be odd positive. */
            ++p;
            if (*p != '2' || (num & 1) == 0 || num == 0) {
                irrep_set_error_("irrep_multiset_2j_parse: half-integer "
                                 "denominator must be 2 with odd numerator "
                                 "at column %ld",
                                 (long)(p - spec + 1));
                irrep_multiset_2j_free(m);
                return NULL;
            }
            ++p;
            two_j = (int)num;
        } else {
            two_j = 2 * (int)num;
        }

        int parity;
        if (*p == 'e')
            parity = +1;
        else if (*p == 'o')
            parity = -1;
        else {
            irrep_set_error_("irrep_multiset_2j_parse: expected 'e' or 'o' "
                             "parity sigil at column %ld",
                             (long)(p - spec + 1));
            irrep_multiset_2j_free(m);
            return NULL;
        }
        ++p;

        irrep_label_2j_t lbl = {.two_j = two_j, .parity = parity};
        if (irrep_multiset_2j_append(m, lbl, (int)mult) != IRREP_OK) {
            irrep_set_error_("irrep_multiset_2j_parse: append failed");
            irrep_multiset_2j_free(m);
            return NULL;
        }

        skip_ws(&p);
        if (*p == '\0')
            break;
        if (*p != '+') {
            irrep_set_error_("irrep_multiset_2j_parse: expected '+' or end "
                             "at column %ld",
                             (long)(p - spec + 1));
            irrep_multiset_2j_free(m);
            return NULL;
        }
        ++p;
    }
    return m;
}

int irrep_multiset_2j_format(const irrep_multiset_2j_t *m, char *buf, size_t buflen) {
    if (!m || !buf || buflen == 0)
        return 0;
    size_t written = 0;
    for (int i = 0; i < m->num_terms; ++i) {
        const irrep_label_2j_t *lbl = &m->labels[i];
        char                    label_buf[32];
        int                     n;
        if (lbl->two_j & 1) {
            n = snprintf(label_buf, sizeof(label_buf), "%dx%d/2%c", m->multiplicities[i],
                         lbl->two_j, lbl->parity > 0 ? 'e' : 'o');
        } else {
            n = snprintf(label_buf, sizeof(label_buf), "%dx%d%c", m->multiplicities[i],
                         lbl->two_j / 2, lbl->parity > 0 ? 'e' : 'o');
        }
        if (n < 0)
            return (int)written;
        int sep = (i > 0);
        if (sep) {
            if (written + 3 < buflen) {
                buf[written++] = ' ';
                buf[written++] = '+';
                buf[written++] = ' ';
            } else
                written += 3;
        }
        if (written + (size_t)n < buflen) {
            memcpy(buf + written, label_buf, (size_t)n);
        }
        written += (size_t)n;
    }
    if (written < buflen)
        buf[written] = '\0';
    else
        buf[buflen - 1] = '\0';
    return (int)written;
}

int irrep_multiset_2j_dim(const irrep_multiset_2j_t *m) {
    return m ? m->total_dim : 0;
}

int irrep_multiset_2j_has_half_integer(const irrep_multiset_2j_t *m) {
    if (!m)
        return 0;
    for (int i = 0; i < m->num_terms; ++i) {
        if (m->labels[i].two_j & 1)
            return 1;
    }
    return 0;
}

int irrep_time_reversal_square_sign_2j(const irrep_multiset_2j_t *m) {
    if (!m || m->num_terms == 0)
        return 0;
    return irrep_multiset_2j_has_half_integer(m) ? -1 : +1;
}
