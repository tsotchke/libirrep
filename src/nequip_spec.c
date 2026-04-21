/* SPDX-License-Identifier: MIT */
/* Spec-string constructor for NequIP layers. Implements the grammar in
 * include/irrep/nequip.h — hand-coded recursive descent, reuses the
 * existing `irrep_multiset_parse` for the two multiset chunks, and parses
 * the option list token by token.
 *
 * Error messages go through the same thread-local last-error channel as the
 * rest of the library; every failure path sets a message naming the column
 * so downstream debuggers can locate the offending token. */

#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/nequip.h>
#include <irrep/multiset.h>

extern void irrep_set_error_(const char *fmt, ...);

/* Defaults per coord doc; kept as macros so the CHANGELOG / header and the
 * implementation stay in lock-step if one ever changes. */
#define IRREP_NEQUIP_SPEC_DEFAULT_SH 2
#define IRREP_NEQUIP_SPEC_DEFAULT_RADIAL 8
#define IRREP_NEQUIP_SPEC_DEFAULT_R_CUT 1.0
#define IRREP_NEQUIP_SPEC_DEFAULT_CUTOFF IRREP_NEQUIP_CUTOFF_POLYNOMIAL
#define IRREP_NEQUIP_SPEC_DEFAULT_POLY_P 6

/* Position cursor that also tracks 1-based column from the original spec. */
typedef struct {
    const char *base;
    const char *p;
} cursor_t;

static int col_of_(const cursor_t *c) {
    return (int)(c->p - c->base) + 1;
}

static void skip_ws_(cursor_t *c) {
    while (*c->p == ' ' || *c->p == '\t')
        ++c->p;
}

/* Extract a NUL-terminated substring between c->p and `end` (exclusive),
 * trimmed of leading/trailing spaces. Caller owns the returned buffer. */
static char *extract_trimmed_(const char *begin, const char *end) {
    while (begin < end && (*begin == ' ' || *begin == '\t'))
        ++begin;
    while (end > begin && (end[-1] == ' ' || end[-1] == '\t'))
        --end;
    size_t n = (size_t)(end - begin);
    char  *out = malloc(n + 1);
    if (!out)
        return NULL;
    memcpy(out, begin, n);
    out[n] = '\0';
    return out;
}

/* Advance past an identifier (lowercase alphabetics + underscore). Writes
 * the identifier into `buf` (at most buflen-1 chars). Returns chars written,
 * or -1 if the identifier is empty or overflows. */
static int consume_ident_(cursor_t *c, char *buf, size_t buflen) {
    size_t n = 0;
    while ((*c->p >= 'a' && *c->p <= 'z') || *c->p == '_') {
        if (n + 1 >= buflen)
            return -1;
        buf[n++] = *c->p++;
    }
    buf[n] = '\0';
    return (int)n == 0 ? -1 : (int)n;
}

/* Consume a non-negative integer. Returns -1 if the cursor doesn't start
 * on a digit. */
static long consume_int_(cursor_t *c) {
    if (!isdigit((unsigned char)*c->p))
        return -1;
    char *end = NULL;
    long  v = strtol(c->p, &end, 10);
    c->p = end ? end : c->p;
    return v;
}

/* Consume a positive real number (int or float, optional exponent). */
static int consume_real_(cursor_t *c, double *out) {
    const char *start = c->p;
    char       *end = NULL;
    double      v = strtod(start, &end);
    if (!end || end == start)
        return -1;
    c->p = end;
    *out = v;
    return 0;
}

static int expect_char_(cursor_t *c, char ch) {
    skip_ws_(c);
    if (*c->p != ch)
        return -1;
    ++c->p;
    return 0;
}

/* Find the "->" separator outside of any brackets. Returns pointer to its
 * first char, or NULL if absent. */
static const char *find_arrow_outside_brackets_(const char *s) {
    int depth = 0;
    for (const char *q = s; *q; ++q) {
        if (*q == '[')
            ++depth;
        else if (*q == ']')
            --depth;
        else if (depth == 0 && q[0] == '-' && q[1] == '>')
            return q;
    }
    return NULL;
}

/* Parse the [option, option, ...] block. On entry `c->p` is at `[`. */
static int parse_options_(cursor_t *c, int *out_sh, int *out_radial, double *out_r_cut,
                          irrep_nequip_cutoff_t *out_cutoff_kind, int *out_cutoff_p) {
    if (expect_char_(c, '[') != 0)
        return -1;

    for (;;) {
        skip_ws_(c);
        if (*c->p == ']') {
            ++c->p;
            break;
        }

        char name[24];
        int  nk = consume_ident_(c, name, sizeof(name));
        if (nk < 0) {
            irrep_set_error_("nequip spec: expected option name at column %d", col_of_(c));
            return -1;
        }
        if (expect_char_(c, '=') != 0) {
            irrep_set_error_("nequip spec: expected '=' after option name "
                             "'%s' at column %d",
                             name, col_of_(c));
            return -1;
        }
        skip_ws_(c);

        if (strcmp(name, "sh") == 0) {
            long v = consume_int_(c);
            if (v < 0 || v > INT_MAX) {
                irrep_set_error_("nequip spec: 'sh' expects non-negative "
                                 "integer ≤ INT_MAX at column %d",
                                 col_of_(c));
                return -1;
            }
            *out_sh = (int)v;
        } else if (strcmp(name, "radial") == 0) {
            long v = consume_int_(c);
            if (v < 1 || v > INT_MAX) {
                irrep_set_error_("nequip spec: 'radial' expects positive "
                                 "integer ≤ INT_MAX at column %d",
                                 col_of_(c));
                return -1;
            }
            *out_radial = (int)v;
        } else if (strcmp(name, "r_cut") == 0) {
            double v;
            if (consume_real_(c, &v) != 0 || v <= 0.0) {
                irrep_set_error_("nequip spec: 'r_cut' expects positive real "
                                 "at column %d",
                                 col_of_(c));
                return -1;
            }
            *out_r_cut = v;
        } else if (strcmp(name, "cutoff") == 0) {
            char kind[24];
            int  kn = consume_ident_(c, kind, sizeof(kind));
            if (kn < 0) {
                irrep_set_error_("nequip spec: expected cutoff kind "
                                 "at column %d",
                                 col_of_(c));
                return -1;
            }
            if (strcmp(kind, "cosine") == 0) {
                *out_cutoff_kind = IRREP_NEQUIP_CUTOFF_COSINE;
            } else if (strcmp(kind, "polynomial") == 0) {
                if (expect_char_(c, '(') != 0) {
                    irrep_set_error_("nequip spec: 'polynomial' expects "
                                     "'(N)' at column %d",
                                     col_of_(c));
                    return -1;
                }
                skip_ws_(c);
                long v = consume_int_(c);
                if (v < 1 || v > INT_MAX) {
                    irrep_set_error_("nequip spec: polynomial order must be "
                                     "in [1, INT_MAX] at column %d",
                                     col_of_(c));
                    return -1;
                }
                if (expect_char_(c, ')') != 0) {
                    irrep_set_error_("nequip spec: expected ')' after "
                                     "polynomial order at column %d",
                                     col_of_(c));
                    return -1;
                }
                *out_cutoff_kind = IRREP_NEQUIP_CUTOFF_POLYNOMIAL;
                *out_cutoff_p = (int)v;
            } else {
                irrep_set_error_("nequip spec: unknown cutoff '%s' (expected "
                                 "'cosine' or 'polynomial') at column %d",
                                 kind, col_of_(c));
                return -1;
            }
        } else {
            irrep_set_error_("nequip spec: unknown option '%s' at column %d", name, col_of_(c));
            return -1;
        }

        skip_ws_(c);
        if (*c->p == ',') {
            ++c->p;
            continue;
        }
        if (*c->p == ']') {
            ++c->p;
            break;
        }
        irrep_set_error_("nequip spec: expected ',' or ']' at column %d", col_of_(c));
        return -1;
    }

    /* Everything after the option block must be whitespace only. */
    skip_ws_(c);
    if (*c->p != '\0') {
        irrep_set_error_("nequip spec: unexpected trailing input at column %d", col_of_(c));
        return -1;
    }
    return 0;
}

irrep_nequip_layer_t *irrep_nequip_layer_from_spec(const char *spec) {
    if (!spec) {
        irrep_set_error_("nequip spec: NULL input");
        return NULL;
    }

    /* Locate the `->` separator — must be outside any bracket pair so a
     * stray `->` inside the (never-legal) option list can't fool us. */
    const char *arrow = find_arrow_outside_brackets_(spec);
    if (!arrow) {
        irrep_set_error_("nequip spec: missing '->' between hidden_in and "
                         "hidden_out multisets");
        return NULL;
    }

    /* LHS = everything before the arrow. */
    char *lhs = extract_trimmed_(spec, arrow);
    if (!lhs || *lhs == '\0') {
        free(lhs);
        irrep_set_error_("nequip spec: hidden_in multiset is empty");
        return NULL;
    }

    /* RHS = everything between `->` and either `[` or end. */
    const char *after_arrow = arrow + 2;
    const char *rhs_end = after_arrow;
    while (*rhs_end && *rhs_end != '[')
        ++rhs_end;
    char *rhs = extract_trimmed_(after_arrow, rhs_end);
    if (!rhs || *rhs == '\0') {
        free(lhs);
        free(rhs);
        irrep_set_error_("nequip spec: hidden_out multiset is empty");
        return NULL;
    }

    /* Parse the two multisets via the existing canonical parser. */
    irrep_multiset_t *h_in = irrep_multiset_parse(lhs);
    irrep_multiset_t *h_out = irrep_multiset_parse(rhs);
    free(lhs);
    free(rhs);
    if (!h_in || !h_out) {
        irrep_multiset_free(h_in);
        irrep_multiset_free(h_out);
        /* irrep_multiset_parse already set a message — prepend context. */
        irrep_set_error_("nequip spec: failed to parse multiset — %s", irrep_last_error());
        return NULL;
    }

    /* Options, or defaults. */
    int                   sh_max = IRREP_NEQUIP_SPEC_DEFAULT_SH;
    int                   n_radial = IRREP_NEQUIP_SPEC_DEFAULT_RADIAL;
    double                r_cut = IRREP_NEQUIP_SPEC_DEFAULT_R_CUT;
    irrep_nequip_cutoff_t cutoff_kind = IRREP_NEQUIP_SPEC_DEFAULT_CUTOFF;
    int                   cutoff_poly = IRREP_NEQUIP_SPEC_DEFAULT_POLY_P;

    if (*rhs_end == '[') {
        cursor_t c = {.base = spec, .p = rhs_end};
        if (parse_options_(&c, &sh_max, &n_radial, &r_cut, &cutoff_kind, &cutoff_poly) != 0) {
            irrep_multiset_free(h_in);
            irrep_multiset_free(h_out);
            return NULL;
        }
    } else if (*rhs_end != '\0') {
        /* Shouldn't reach — the while loop above stops at `[` or `\0`. */
        irrep_multiset_free(h_in);
        irrep_multiset_free(h_out);
        irrep_set_error_("nequip spec: stray characters after hidden_out "
                         "multiset");
        return NULL;
    }

    irrep_nequip_layer_t *layer =
        irrep_nequip_layer_build(h_in, sh_max, n_radial, r_cut, cutoff_kind, cutoff_poly, h_out);

    /* `irrep_nequip_layer_build` copies the multisets internally, so we
     * release our originals regardless of build success. */
    irrep_multiset_free(h_in);
    irrep_multiset_free(h_out);
    return layer;
}
