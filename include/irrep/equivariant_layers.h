/* SPDX-License-Identifier: MIT */
/** @file equivariant_layers.h
 *  @brief Minimal equivariant-MLP building blocks: per-irrep linear mixing,
 *         RMS normalisation per block, and multiplicative gates.
 *
 *  Each primitive commutes with the real @ref irrep_wigner_D_multiset action
 *  by construction, so stacking them preserves SO(3) equivariance.
 */
#ifndef IRREP_EQUIVARIANT_LAYERS_H
#define IRREP_EQUIVARIANT_LAYERS_H

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Opaque descriptor for a linear-on-irreps layer. */
typedef struct irrep_linear irrep_linear_t;

/** @brief Compile a linear descriptor that mixes channels within each
 *         `(l, parity)` block from @p in to @p out. */
IRREP_API irrep_linear_t *irrep_linear_build(const irrep_multiset_t *in,
                                             const irrep_multiset_t *out, int in_channels,
                                             int out_channels);
/** @brief Release a linear descriptor. */
IRREP_API void irrep_linear_free(irrep_linear_t *lin);

/** @brief Total number of learnable scalars. */
IRREP_API int irrep_linear_num_weights(const irrep_linear_t *lin);

/** @brief Forward: `out = W · in` per irrep block. */
IRREP_API void irrep_linear_apply(const irrep_linear_t *lin, const double *weights,
                                  const double *in, double *out);

/** @brief Backward pass — accumulates into @p grad_weights and @p grad_in
 *         (caller zero-initialises). */
IRREP_API void irrep_linear_backward(const irrep_linear_t *lin, const double *weights,
                                     const double *in, const double *grad_out, double *grad_weights,
                                     double *grad_in);

/** @brief RMS-norm per irrep block. `scales[t · channels + c]` rescales the
 *         output of term @p t, channel @p c; equivariant by construction. */
IRREP_API void irrep_norm_rms(const irrep_multiset_t *m, int channels, const double *scales,
                              const double *in, double *out);

/** @brief RMS-norm backward — accumulates into @p grad_scales and @p grad_in. */
IRREP_API void irrep_norm_rms_backward(const irrep_multiset_t *m, int channels,
                                       const double *scales, const double *in,
                                       const double *grad_out, double *grad_scales,
                                       double *grad_in);

/** @brief Multiplicative gate: `out = scalar_gate · in` per `(term, channel)`.
 *         Caller supplies the activated gate values (sigmoid, tanh, …). */
IRREP_API void irrep_gate_apply(const irrep_multiset_t *m, int channels, const double *scalar_gates,
                                const double *in, double *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_EQUIVARIANT_LAYERS_H */
