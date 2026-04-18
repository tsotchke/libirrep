/* SPDX-License-Identifier: MIT */
#ifndef IRREP_PARITY_H
#define IRREP_PARITY_H

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Parity of a single irrep label (+1 or -1). */
IRREP_API int irrep_parity(irrep_label_t label);

/* Parity of the direct-product irrep a ⊗ b (parity_a * parity_b). */
IRREP_API int irrep_parity_product(irrep_label_t a, irrep_label_t b);

/* Filter a flat list of path triplets (i_a, i_b, i_c) in-place, removing any
 * path whose parities do not satisfy parity_a * parity_b == parity_c. Returns
 * the new path count. */
IRREP_API int irrep_parity_filter_paths(const irrep_multiset_t *a,
                                        const irrep_multiset_t *b,
                                        const irrep_multiset_t *c,
                                        int *paths,
                                        int num_paths);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_PARITY_H */
