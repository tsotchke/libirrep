/* SPDX-License-Identifier: MIT */
/** @file tensor_product.h
 *  @brief e3nn-style path-indexed tensor products on irrep multisets, forward
 *         and backward, scalar and batched.
 *
 *  Given features `a ∈ A = ⊕ mul_a · irrep_a` and `b ∈ B = ⊕ mul_b · irrep_b`,
 *  produce `c ∈ C = ⊕ mul_c · irrep_c` by contracting along CG coefficients
 *  for a selected set of paths `(l_a, l_b, l_c)`:
 *
 *  @verbatim
 *    c[l_c, m_c, ch_c] = Σ_path  w[path] · Σ_{m_a, m_b}  ⟨l_a m_a; l_b m_b | l_c m_c⟩
 *                                        · a[l_a, m_a, ch_a(path)]
 *                                        · b[l_b, m_b, ch_b(path)]
 *  @endverbatim
 *
 *  The real-basis coefficients are computed from complex-basis CGs with an
 *  `i^{l_a + l_b - l_c}` phase correction and the real↔complex change of
 *  basis, so equivariance holds under the real `D^l(R)` on each space.
 *
 *  Two channel modes are provided:
 *
 *  - #IRREP_TP_MODE_UUU — one scalar weight per path. Requires matching
 *    multiplicities across `a, b, c` for that path.
 *  - #IRREP_TP_MODE_UVW — independent multiplicities per space. Weight shape
 *    per path is `[W, V, U]`; the apply/backward kernels mix channels.
 */
#ifndef IRREP_TENSOR_PRODUCT_H
#define IRREP_TENSOR_PRODUCT_H

#include <stddef.h>

#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>
#include <irrep/multiset_2j.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Opaque handle describing a compiled tensor-product plan. */
typedef struct tp_descriptor tp_descriptor_t;

/** @brief Channel-weight layout for tensor-product apply/backward. */
typedef enum {
    IRREP_TP_MODE_UUU = 0,  /**< one scalar per path; mult_a = mult_b = mult_c. */
    IRREP_TP_MODE_UVW = 1   /**< full [W, V, U] weight tensor per path. */
} irrep_tp_mode_t;

/** @brief Compile a UUU tensor-product descriptor from @p a, @p b, @p c and
 *         selected paths.
 *  @param a                   input multiset A.
 *  @param b                   input multiset B.
 *  @param c                   output multiset C.
 *  @param selected_paths      flat array of triplets (idx_a, idx_b, idx_c);
 *                             pass @c NULL to auto-enumerate all valid paths.
 *  @param num_selected_paths  number of triplets (not array length).
 *  @return opaque descriptor, or @c NULL on failure (see @c irrep_last_error()). */
IRREP_API tp_descriptor_t *irrep_tp_build(const irrep_multiset_t *a,
                                          const irrep_multiset_t *b,
                                          const irrep_multiset_t *c,
                                          const int *selected_paths,
                                          int num_selected_paths);

/** @brief Release resources owned by a descriptor. */
IRREP_API void irrep_tp_free(tp_descriptor_t *desc);

/** @brief Enumerate all CG-valid paths obeying `|l_a − l_b| ≤ l_c ≤ l_a + l_b`
 *         and parity consistency.
 *  @param a, b, c   multisets.
 *  @param out_paths buffer for flattened triplets; may be @c NULL to just count.
 *  @param max_paths capacity in triplets.
 *  @return total number of valid paths (may exceed @p max_paths). */
IRREP_API int irrep_tp_enumerate_paths(const irrep_multiset_t *a,
                                       const irrep_multiset_t *b,
                                       const irrep_multiset_t *c,
                                       int *out_paths,
                                       int  max_paths);

/** @brief Forward: `c = a ⊗ b` with unit weights per path. */
IRREP_API void irrep_tp_apply(const tp_descriptor_t *desc,
                              const double *a_in,
                              const double *b_in,
                              double       *c_out);

/** @brief Forward with one learnable scalar weight per selected path. */
IRREP_API void irrep_tp_apply_weighted(const tp_descriptor_t *desc,
                                       const double *weights,
                                       const double *a_in,
                                       const double *b_in,
                                       double       *c_out);

/** @brief Unweighted backward. Accumulates into @p grad_a, @p grad_b (caller
 *         zero-initialises). */
IRREP_API void irrep_tp_apply_backward(const tp_descriptor_t *desc,
                                       const double *a_in,
                                       const double *b_in,
                                       const double *grad_c_out,
                                       double       *grad_a,
                                       double       *grad_b);

/** @brief Weighted backward: also populates @p grad_w, one scalar per path. */
IRREP_API void irrep_tp_apply_backward_weighted(const tp_descriptor_t *desc,
                                                const double *weights,
                                                const double *a_in,
                                                const double *b_in,
                                                const double *grad_c_out,
                                                double       *grad_a,
                                                double       *grad_b,
                                                double       *grad_w);

/** @name Batched variants
 *  Each of @p a_in, @p b_in, @p c_out is a buffer of @p batch contiguous
 *  blocks of the respective multiset's total dimension.
 *  @{ */
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
/** @} */

/** @brief Total dimension of the descriptor's output space. */
IRREP_API int               irrep_tp_output_dim(const tp_descriptor_t *desc);
/** @brief Number of selected paths in the descriptor. */
IRREP_API int               irrep_tp_num_paths (const tp_descriptor_t *desc);
/** @brief Channel mode the descriptor was built in. */
IRREP_API irrep_tp_mode_t   irrep_tp_mode      (const tp_descriptor_t *desc);

/** @name Half-integer (spinor) tensor products (1.3)
 *  These operate on doubled-integer multisets @ref irrep_multiset_2j_t and on
 *  **complex** amplitudes (half-integer irreps have no natural real basis;
 *  the tensor product is genuinely complex). Selection rules: triangle on
 *  `two_j`, parity multiplication, and `two_j_a + two_j_b + two_j_c ≡ 0 (mod 2)`
 *  (automatic under the triangle). Clebsch-Gordan coefficients come from
 *  @ref irrep_cg_2j and are real (Condon-Shortley).
 *  @{ */

/** @brief Opaque half-integer tensor-product descriptor. */
typedef struct tp_2j_descriptor tp_2j_descriptor_t;

/** @brief Enumerate valid `(i_a, i_b, i_c)` paths on doubled-integer multisets. */
IRREP_API int irrep_tp_2j_enumerate_paths(const irrep_multiset_2j_t *a,
                                          const irrep_multiset_2j_t *b,
                                          const irrep_multiset_2j_t *c,
                                          int *out_paths, int max_paths);

/** @brief Compile a UVW-mode half-integer tensor-product descriptor.
 *
 *  @param selected_paths      triplet list `(i_a, i_b, i_c)`; pass `NULL`
 *                             to auto-enumerate all valid paths.
 *  @param num_selected_paths  number of triplets.
 *  @return opaque descriptor, or `NULL` on failure. */
IRREP_API tp_2j_descriptor_t *irrep_tp_2j_build(
    const irrep_multiset_2j_t *a,
    const irrep_multiset_2j_t *b,
    const irrep_multiset_2j_t *c,
    const int *selected_paths, int num_selected_paths);

/** @brief Release resources owned by a half-integer tp descriptor. */
IRREP_API void irrep_tp_2j_free(tp_2j_descriptor_t *desc);

/** @brief Unweighted forward for half-integer TP: complex `c = a ⊗ b`. */
IRREP_API void irrep_tp_2j_apply(const tp_2j_descriptor_t *desc,
                                 const double _Complex *a_in,
                                 const double _Complex *b_in,
                                 double _Complex *c_out);

/** @brief Weighted forward: one complex scalar per path.
 *         `c = Σ_path w_p · (a_p ⊗ b_p)`. */
IRREP_API void irrep_tp_2j_apply_weighted(const tp_2j_descriptor_t *desc,
                                          const double _Complex *weights,
                                          const double _Complex *a_in,
                                          const double _Complex *b_in,
                                          double _Complex *c_out);

/** @brief Backward of @ref irrep_tp_2j_apply_weighted. Accumulates into
 *         `grad_a, grad_b, grad_w` (caller zero-initialises). */
IRREP_API void irrep_tp_2j_apply_backward(const tp_2j_descriptor_t *desc,
                                          const double _Complex *weights,
                                          const double _Complex *a_in,
                                          const double _Complex *b_in,
                                          const double _Complex *grad_c_out,
                                          double _Complex *grad_a,
                                          double _Complex *grad_b,
                                          double _Complex *grad_w);

/** @brief Output dimension of a half-integer tp descriptor. */
IRREP_API int irrep_tp_2j_output_dim(const tp_2j_descriptor_t *desc);

/** @brief Number of selected paths in a half-integer tp descriptor. */
IRREP_API int irrep_tp_2j_num_paths(const tp_2j_descriptor_t *desc);
/** @} */

/** @name UVW mode
 *  Per-path weight tensor shape is `[W, V, U]`; full channel mixing between
 *  the three multisets.
 *  @{ */

/** @brief Compile a UVW-mode descriptor. Parameters match #irrep_tp_build. */
IRREP_API tp_descriptor_t *irrep_tp_build_uvw(
    const irrep_multiset_t *a, const irrep_multiset_t *b,
    const irrep_multiset_t *c,
    const int *selected_paths, int num_selected_paths);

/** @brief Total weight count for a UVW descriptor: Σ_path (W_p · V_p · U_p). */
IRREP_API int irrep_tp_num_weights_uvw(const tp_descriptor_t *desc);

/** @brief UVW forward. @p weights has `irrep_tp_num_weights_uvw(desc)` entries,
 *         packed path-major, then (w, v, u). */
IRREP_API void irrep_tp_apply_uvw(
    const tp_descriptor_t *desc,
    const double *weights,
    const double *a_in, const double *b_in, double *c_out);

/** @brief UVW backward. Accumulates gradients into @p grad_a, @p grad_b, @p grad_w. */
IRREP_API void irrep_tp_apply_uvw_backward(
    const tp_descriptor_t *desc,
    const double *weights,
    const double *a_in, const double *b_in,
    const double *grad_c_out,
    double *grad_a, double *grad_b, double *grad_w);

/** @brief Per-path L2 norm of UVW weights:
 *         `out[k] = Σ_{w, v, u} weights[k, w, v, u]²`.
 *  @p out has `irrep_tp_num_paths(desc)` entries, path-major. Useful for
 *  SR-style holomorphic training where different paths have different
 *  curvatures on the variational manifold and deserve path-dependent
 *  regularisation strengths — combine as `Σ_k λ_k · out[k]`. */
IRREP_API void irrep_tp_weight_l2_per_path_uvw(
    const tp_descriptor_t *desc,
    const double *weights,
    double       *out);

/** @brief Backward of #irrep_tp_weight_l2_per_path_uvw.
 *         `grad_weights[k, w, v, u] += grad_per_path[k] · 2 · weights[k, w, v, u]`.
 *  Accumulates — caller zero-initialises @p grad_weights. */
IRREP_API void irrep_tp_weight_l2_per_path_uvw_backward(
    const tp_descriptor_t *desc,
    const double *weights,
    const double *grad_per_path,
    double       *grad_weights);
/** @} */

#ifdef __cplusplus
}
#endif

#endif /* IRREP_TENSOR_PRODUCT_H */
