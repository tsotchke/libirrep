/* SPDX-License-Identifier: MIT */
#ifndef IRREP_TENSOR_PRODUCT_H
#define IRREP_TENSOR_PRODUCT_H

#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tp_descriptor tp_descriptor_t;

/* Build a tensor-product descriptor from irrep multisets a, b, c and a list of
 * selected paths as triplets (idx_a, idx_b, idx_c) flattened into `selected_paths`.
 * Returns NULL on failure (reason via irrep_last_error()). */
IRREP_API tp_descriptor_t *irrep_tp_build(const irrep_multiset_t *a,
                                          const irrep_multiset_t *b,
                                          const irrep_multiset_t *c,
                                          const int *selected_paths,
                                          int num_selected_paths);

IRREP_API void irrep_tp_free(tp_descriptor_t *desc);

/* Enumerate valid paths (i.e. triples (i_a, i_b, i_c) with
 * |l_a - l_b| <= l_c <= l_a + l_b and parity_a * parity_b == parity_c).
 * If `out_paths` is non-NULL, writes up to `max_paths` triplets flattened;
 * returns the total count. */
IRREP_API int irrep_tp_enumerate_paths(const irrep_multiset_t *a,
                                       const irrep_multiset_t *b,
                                       const irrep_multiset_t *c,
                                       int *out_paths,
                                       int  max_paths);

/* Forward: c = a ⊗ b along the descriptor's selected paths. */
IRREP_API void irrep_tp_apply(const tp_descriptor_t *desc,
                              const double *a_in,
                              const double *b_in,
                              double       *c_out);

IRREP_API void irrep_tp_apply_weighted(const tp_descriptor_t *desc,
                                       const double *weights,  /* num_selected_paths */
                                       const double *a_in,
                                       const double *b_in,
                                       double       *c_out);

/* Backward: grad_a and grad_b accumulate; caller pre-zeros them. */
IRREP_API void irrep_tp_apply_backward(const tp_descriptor_t *desc,
                                       const double *a_in,
                                       const double *b_in,
                                       const double *grad_c_out,
                                       double       *grad_a,
                                       double       *grad_b);

IRREP_API void irrep_tp_apply_backward_weighted(const tp_descriptor_t *desc,
                                                const double *weights,
                                                const double *a_in,
                                                const double *b_in,
                                                const double *grad_c_out,
                                                double       *grad_a,
                                                double       *grad_b,
                                                double       *grad_w);

/* Batched variants. Each input and output is stride-separated by its total_dim. */
IRREP_API void irrep_tp_apply_batch                (const tp_descriptor_t *desc,
                                                    size_t batch,
                                                    const double *a_in,
                                                    const double *b_in,
                                                    double       *c_out);
IRREP_API void irrep_tp_apply_weighted_batch       (const tp_descriptor_t *desc,
                                                    size_t batch,
                                                    const double *weights,
                                                    const double *a_in,
                                                    const double *b_in,
                                                    double       *c_out);
IRREP_API void irrep_tp_apply_backward_batch       (const tp_descriptor_t *desc,
                                                    size_t batch,
                                                    const double *a_in,
                                                    const double *b_in,
                                                    const double *grad_c_out,
                                                    double       *grad_a,
                                                    double       *grad_b);
IRREP_API void irrep_tp_apply_backward_weighted_batch(const tp_descriptor_t *desc,
                                                     size_t batch,
                                                     const double *weights,
                                                     const double *a_in,
                                                     const double *b_in,
                                                     const double *grad_c_out,
                                                     double       *grad_a,
                                                     double       *grad_b,
                                                     double       *grad_w);

/* Descriptor introspection. */
IRREP_API int irrep_tp_output_dim(const tp_descriptor_t *desc);
IRREP_API int irrep_tp_num_paths (const tp_descriptor_t *desc);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_TENSOR_PRODUCT_H */
