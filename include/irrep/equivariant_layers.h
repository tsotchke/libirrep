/* SPDX-License-Identifier: MIT */
#ifndef IRREP_EQUIVARIANT_LAYERS_H
#define IRREP_EQUIVARIANT_LAYERS_H

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---------- linear-on-irreps ---------- *
 * Per-irrep-block dense mixing across channels; preserves (l, parity). */
typedef struct irrep_linear irrep_linear_t;

IRREP_API irrep_linear_t *irrep_linear_build(const irrep_multiset_t *in,
                                             const irrep_multiset_t *out,
                                             int in_channels,
                                             int out_channels);
IRREP_API void irrep_linear_free(irrep_linear_t *lin);

IRREP_API int  irrep_linear_num_weights(const irrep_linear_t *lin);

IRREP_API void irrep_linear_apply(const irrep_linear_t *lin,
                                  const double *weights,
                                  const double *in,
                                  double       *out);

IRREP_API void irrep_linear_backward(const irrep_linear_t *lin,
                                     const double *weights,
                                     const double *in,
                                     const double *grad_out,
                                     double       *grad_weights,
                                     double       *grad_in);

/* ---------- RMS norm per irrep block ---------- *
 * One learnable scale per (term_idx, channel). */
IRREP_API void irrep_norm_rms(const irrep_multiset_t *m,
                              int channels,
                              const double *scales,
                              const double *in,
                              double       *out);

IRREP_API void irrep_norm_rms_backward(const irrep_multiset_t *m,
                                       int channels,
                                       const double *scales,
                                       const double *in,
                                       const double *grad_out,
                                       double       *grad_scales,
                                       double       *grad_in);

/* ---------- gate activation ---------- *
 * Caller provides per-(term, channel) scalar gates (after sigmoid, tanh, etc.),
 * which multiply the corresponding irrep block. */
IRREP_API void irrep_gate_apply(const irrep_multiset_t *m,
                                int channels,
                                const double *scalar_gates,
                                const double *in,
                                double       *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_EQUIVARIANT_LAYERS_H */
