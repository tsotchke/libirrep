/* SPDX-License-Identifier: MIT */
/** @file nequip.h
 *  @brief NequIP-style E(3)-equivariant message-passing layer.
 *
 *  Composes the following into one opaque object:
 *
 *   1. Edge embedding — `Y(r̂_ij)` cartesian real spherical harmonics up to
 *      `l_sh_max` on each edge unit-vector.
 *   2. Radial basis — `φ_n(r_ij)` Bessel-RBF expansion of the edge length,
 *      `n_radial` channels.
 *   3. Cutoff — polynomial or cosine cutoff smoothly vanishing at `r_cut`.
 *   4. Tensor product — weighted UVW tensor product of each neighbour's
 *      hidden features `h_j` with `Y(r̂_ij)`, producing an edge message
 *      `m_ij`.
 *   5. Aggregation — `Σ_j m_ij` onto node `i`, scaled by `φ_n · cutoff`.
 *
 *  One call (`irrep_nequip_layer_apply`) consumes a small graph plus per-edge
 *  geometry and produces new node features. Backward, geometry-force, and
 *  weight-regulariser paths are wired through the underlying TP primitives.
 *
 *  Weight layout: `tp_weights` matches `irrep_tp_apply_uvw` — flat,
 *  path-indexed, `irrep_nequip_layer_num_weights(layer)` scalars total. No
 *  separate learnable radial head here; wrap this layer with additional
 *  logic if you want one. */

#ifndef IRREP_NEQUIP_H
#define IRREP_NEQUIP_H

#include <stddef.h>

#include <irrep/export.h>
#include <irrep/multiset.h>
#include <irrep/tensor_product.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Opaque compiled NequIP layer descriptor. */
typedef struct irrep_nequip_layer irrep_nequip_layer_t;

/** @brief Cutoff-function family. */
typedef enum {
    IRREP_NEQUIP_CUTOFF_COSINE     = 0,  /**< `½ (1 + cos(π r / r_cut))`. */
    IRREP_NEQUIP_CUTOFF_POLYNOMIAL = 1   /**< NequIP polynomial cutoff (smooth to order p). */
} irrep_nequip_cutoff_t;

/** @brief Compile a NequIP layer from a single spec string.
 *
 *  Shorthand for #irrep_nequip_layer_build that bundles the input multiset,
 *  output multiset, and edge-side parameters into one e3nn-style string.
 *  The spec grammar is
 *
 *  @verbatim
 *    layer_spec ::= multiset "->" multiset option_list?
 *    option_list ::= "[" option ( "," option )* "]"
 *    option      ::= "sh"     "=" int
 *                 | "radial" "=" int
 *                 | "r_cut"  "=" float
 *                 | "cutoff" "=" ( "cosine" | "polynomial(" int ")" )
 *  @endverbatim
 *
 *  Whitespace is lenient between tokens, strict within them (`"1x0e"`, not
 *  `"1 x 0 e"`). Omitted options take defaults `sh = 2`, `radial = 8`,
 *  `r_cut = 1.0`, `cutoff = polynomial(6)`.
 *
 *  Example: `"4x0e + 2x1o + 1x2e -> 2x0e + 1x1o [sh=3, radial=8, r_cut=1.5]"`.
 *
 *  @return opaque layer, or @c NULL on parse / build failure (see
 *          @ref irrep_last_error). */
IRREP_API irrep_nequip_layer_t *irrep_nequip_layer_from_spec(const char *spec);

/** @brief Compile a NequIP layer with the requested hidden-space shapes and
 *         edge-side parameters.
 *  @param hidden_in      input multiset on each node.
 *  @param l_sh_max       maximum SH order used on the edge vector.
 *  @param n_radial       number of Bessel-RBF channels.
 *  @param r_cut          cutoff radius (edges beyond are ignored).
 *  @param cutoff_kind    cosine or polynomial.
 *  @param cutoff_poly_p  polynomial order (≥ 1) when @p cutoff_kind is polynomial.
 *  @param hidden_out     output multiset.
 *  @return opaque descriptor, or @c NULL on failure (see @ref irrep_last_error). */
IRREP_API irrep_nequip_layer_t *irrep_nequip_layer_build(
    const irrep_multiset_t *hidden_in,
    int                     l_sh_max,
    int                     n_radial,
    double                  r_cut,
    irrep_nequip_cutoff_t   cutoff_kind,
    int                     cutoff_poly_p,
    const irrep_multiset_t *hidden_out);

/** @brief Release a descriptor built by #irrep_nequip_layer_build. */
IRREP_API void irrep_nequip_layer_free(irrep_nequip_layer_t *layer);

/** @brief Total number of learnable TP weights expected by
 *         #irrep_nequip_layer_apply. */
IRREP_API int  irrep_nequip_layer_num_weights(const irrep_nequip_layer_t *layer);

/** @brief Forward pass.
 *
 *  Computes per-node `h_out[i] = Σ_{j ∈ N(i)} tp(h_in[j], Y(r̂_ij)) · rbf(r_ij) · cutoff(r_ij)`.
 *
 *  @param layer        descriptor from #irrep_nequip_layer_build.
 *  @param tp_weights   `num_weights(layer)` scalars.
 *  @param n_nodes      node count.
 *  @param n_edges      edge count.
 *  @param edge_src     length @p n_edges; source node index per edge (j).
 *  @param edge_dst     length @p n_edges; target node index per edge (i).
 *  @param edge_vec     length `n_edges * 3`; real 3-vector `r_ij = r_j − r_i`.
 *                      Not required to be unit-norm; normalisation happens inside.
 *  @param h_in         length `n_nodes * hidden_in->total_dim`; input node features.
 *  @param h_out        length `n_nodes * hidden_out->total_dim`; **zeroed internally**
 *                      before message accumulation. */
IRREP_API void irrep_nequip_layer_apply(
    const irrep_nequip_layer_t *layer,
    const double               *tp_weights,
    int                         n_nodes,
    int                         n_edges,
    const int                  *edge_src,
    const int                  *edge_dst,
    const double               *edge_vec,
    const double               *h_in,
    double                     *h_out);

/** @brief Backward pass through hidden features and weights.
 *
 *  Convention (match PyTorch / e3nn): @p grad_h_in and @p grad_tp_weights
 *  are **accumulated** (`+=`). Caller MUST zero them before the first call of
 *  a training iteration; uninitialised input gives garbage output. This lets
 *  you share buffers across edges or micro-batches without repeated alloc.
 *
 *  Edge-geometry gradients are deliberately not produced here — use
 *  #irrep_nequip_layer_apply_forces if you need them.
 *
 *  @param grad_h_out    pullback of the loss w.r.t. the layer output.
 *  @param grad_h_in     out: accumulated ∂L/∂h_in.
 *  @param grad_tp_weights  out: accumulated ∂L/∂tp_weights. */
IRREP_API void irrep_nequip_layer_apply_backward(
    const irrep_nequip_layer_t *layer,
    const double               *tp_weights,
    int                         n_nodes,
    int                         n_edges,
    const int                  *edge_src,
    const int                  *edge_dst,
    const double               *edge_vec,
    const double               *h_in,
    const double               *grad_h_out,
    double                     *grad_h_in,
    double                     *grad_tp_weights);

/** @brief Edge-geometry gradient `∂L/∂edge_vec[e, axis]`.
 *
 *  Chain rule:
 *  @verbatim
 *    h_out[dst] += scale(r) · tp(h_src, Y(r̂))
 *    ∂h_out/∂edge_vec[e] = (∂scale/∂r · r̂) · tp_out
 *                        + (scale / r) · (∂tp/∂Y_m) · ∂Y_m/∂r̂
 *  @endverbatim
 *
 *  No `∂/∂r` term on the TP output because `Y` depends only on `r̂`; no
 *  radial component of `∂Y/∂r̂` because the SH gradient is already tangent
 *  to `S²` (asserted to 1e-8 for `l ≤ 3` in the SH test suite).
 *
 *  @param grad_edge_vec  out: accumulated `∂L/∂edge_vec`, length `n_edges * 3`.
 *                        Caller pre-zeros. */
IRREP_API void irrep_nequip_layer_apply_forces(
    const irrep_nequip_layer_t *layer,
    const double               *tp_weights,
    int                         n_nodes,
    int                         n_edges,
    const int                  *edge_src,
    const int                  *edge_dst,
    const double               *edge_vec,
    const double               *h_in,
    const double               *grad_h_out,
    double                     *grad_edge_vec);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_NEQUIP_H */
