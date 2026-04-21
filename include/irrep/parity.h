/* SPDX-License-Identifier: MIT */
/** @file parity.h
 *  @brief O(3) parity helpers: labels, products, and path filtering.
 */
#ifndef IRREP_PARITY_H
#define IRREP_PARITY_H

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief `+1` or `−1` — parity component of an irrep label. */
IRREP_API int irrep_parity(irrep_label_t label);

/** @brief Parity of the tensor-product irrep `a ⊗ b` (product of parities). */
IRREP_API int irrep_parity_product(irrep_label_t a, irrep_label_t b);

/** @brief Filter a flat list of path triplets `(i_a, i_b, i_c)` in place,
 *         removing any path whose parities violate `p_a · p_b = p_c`.
 *  @return surviving path count. */
IRREP_API int irrep_parity_filter_paths(const irrep_multiset_t *a, const irrep_multiset_t *b,
                                        const irrep_multiset_t *c, int *paths, int num_paths);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_PARITY_H */
