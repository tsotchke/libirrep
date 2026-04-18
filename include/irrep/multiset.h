/* SPDX-License-Identifier: MIT */
#ifndef IRREP_MULTISET_H
#define IRREP_MULTISET_H

#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

IRREP_API irrep_multiset_t *irrep_multiset_new  (int capacity);
IRREP_API void              irrep_multiset_free (irrep_multiset_t *m);
IRREP_API irrep_status_t    irrep_multiset_append(irrep_multiset_t *m,
                                                  irrep_label_t label,
                                                  int multiplicity);

/* Parse e3nn-style spec: "1x0e + 2x1o + 1x2e". Returns NULL on malformed input
 * (the reason is available via irrep_last_error()). */
IRREP_API irrep_multiset_t *irrep_multiset_parse(const char *spec);

/* Format into `buf`; returns required length (excluding terminator). Writes at
 * most buflen-1 characters plus a terminator. */
IRREP_API int  irrep_multiset_format(const irrep_multiset_t *m, char *buf, size_t buflen);

/* Merge like terms and sort canonically (by (l, parity)). Idempotent. */
IRREP_API void irrep_multiset_simplify(irrep_multiset_t *m);

/* Direct sum m1 ⊕ m2. Returns a newly-allocated multiset. */
IRREP_API irrep_multiset_t *irrep_multiset_direct_sum(const irrep_multiset_t *m1,
                                                      const irrep_multiset_t *m2);

/* Convenience accessors. */
IRREP_API int irrep_multiset_dim         (const irrep_multiset_t *m);
IRREP_API int irrep_multiset_block_offset(const irrep_multiset_t *m, int term_idx);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_MULTISET_H */
