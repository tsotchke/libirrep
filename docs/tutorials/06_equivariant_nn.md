# Tutorial 06 — Equivariant neural networks

This tutorial assembles the libirrep primitives into an equivariant
neural-network stack. The goal is not to reproduce the state of the
art — that is what NequIP (Batzner et al. 2022), MACE (Batatia et al.
2022), and Allegro (Musaelian et al. 2023) already do — but to make the
relationship between the library's building blocks and a working
equivariant GNN concrete and end-to-end verifiable.

## 1. Equivariance: the formal requirement

Let `V` and `W` be vector spaces carrying representations `D_V` and
`D_W` of a group `G`. A function `f : V → W` is **equivariant** if

```
f(D_V(g) · v) = D_W(g) · f(v) for all g ∈ G, v ∈ V.
```

For equivariant graph neural networks targeting molecular / materials
systems, the group of interest is E(3) = SO(3) ⋉ ℝ³ (rotations and
translations of `ℝ³`). Translations are trivialised by using relative
vectors `r_ij = r_j − r_i` between atoms, so the operative symmetry
reduces to SO(3). If the system has inversion symmetry — magnetic
Hamiltonians, atomic electronic structure — we extend to O(3), the
representations acquire a parity index, and libirrep's `(l, p)` irrep
labels become material.

The practical engineering consequence: every operation in the network
must be either a linear combination of equivariant operations
(preserving equivariance by closure) or an explicit projection onto an
irrep subspace (preserving equivariance by construction). Libirrep
exposes exactly these:

- **Linear-on-irreps** (`irrep_linear_*`) — channel-mixing within each
 `(l, p)` block, preserving the irrep structure.
- **RMS-norm per block** (`irrep_norm_rms`) — `O(l)`-scalar-gated
 rescaling of each block, preserving equivariance because the norm of
 a rotation-covariant vector is a scalar.
- **Gate activation** (`irrep_gate_apply`) — multiplicative gate driven
 by a separate scalar feature; non-linear, equivariant by construction.
- **Tensor product** (`irrep_tp_apply_*`) — the main non-linear
 operation, combining features in different irreps via
 Clebsch-Gordan-weighted contractions.
- **Point-group projection** (`irrep_pg_project`) — projection onto
 the subspace invariant under a discrete subgroup of O(3), for
 symmetric neural quantum state (NQS) ansätze that live in a
 specific irrep of the point group.

## 2. The NequIP message-passing layer

NequIP (Batzner et al. 2022) composes five operations into one
message-passing layer:

1. **Edge embedding.** For each directed edge `(i, j)`, compute the
 cartesian spherical harmonics `Y(r̂_ij)` up to some `l_sh_max`,
 producing an irrep-valued feature on each edge.
2. **Radial basis.** Expand the edge length `r_ij` in a set of
 `n_radial` Bessel functions `φ_n(r_ij)`.
3. **Cutoff.** Multiply by a smooth envelope `u(r_ij)` that vanishes at
 and beyond `r_cut` (cosine or NequIP-polynomial form), ensuring
 the interaction has compact support.
4. **Tensor product.** Contract each neighbour's hidden feature
 `h_j ∈ A` with `Y(r̂_ij) ∈ B` via a UVW tensor product
 `m_ij ∈ C`, parameterised by learnable weights.
5. **Aggregation.** Sum `m_ij` over neighbours to produce the updated
 node feature `h_i' ∈ C`.

Libirrep ships this composition as the first-class
`irrep_nequip_layer_t`:

```c
#include <irrep/multiset.h>
#include <irrep/nequip.h>

irrep_multiset_t *h_in = irrep_multiset_parse("2x0e + 1x1o");
irrep_multiset_t *h_out = irrep_multiset_parse("2x0e + 1x1o");

irrep_nequip_layer_t *layer = irrep_nequip_layer_build(
 h_in,
 /*l_sh_max=*/ 2,
 /*n_radial=*/ 4,
 /*r_cut=*/ 3.0,
 IRREP_NEQUIP_CUTOFF_POLYNOMIAL,
 /*cutoff_poly_p=*/ 6,
 h_out);

int nw = irrep_nequip_layer_num_weights(layer);
```

Or via the spec-string constructor (1.2+):

```c
irrep_nequip_layer_t *layer = irrep_nequip_layer_from_spec(
 "2x0e + 1x1o -> 2x0e + 1x1o [sh=2, radial=4, r_cut=3.0, cutoff=polynomial(6)]");
```

The string grammar, with defaults, is in
[`include/irrep/nequip.h`](../../include/irrep/nequip.h); consult the
Doxygen block on `irrep_nequip_layer_from_spec`.

## 3. Running a forward pass

```c
int n_nodes = 3;
int n_edges = 4;
int edge_src[4] = { 0, 1, 1, 2 };
int edge_dst[4] = { 1, 0, 2, 1 };
double edge_vec[4 * 3] = {
 0.5, -0.3, 0.8,
 -0.5, 0.3, -0.8,
 0.2, 0.7, -0.4,
 -0.2, -0.7, 0.4,
};

double *w = calloc(nw, sizeof(double));
double *h_in_buf = calloc(n_nodes * h_in->total_dim, sizeof(double));
double *h_out_buf = calloc(n_nodes * h_out->total_dim, sizeof(double));
/* ... fill w from your optimiser, h_in_buf from the upstream layer ... */

irrep_nequip_layer_apply(layer, w, n_nodes, n_edges,
 edge_src, edge_dst, edge_vec,
 h_in_buf, h_out_buf);
/* h_out_buf is zeroed internally before message accumulation. */
```

Weights are laid out UVW-style: per selected `(i_a, i_b, i_c)` path, a
`[W, V, U]` tensor. The total count is
`irrep_nequip_layer_num_weights(layer)`.

## 4. Backward passes

Two distinct gradient paths:

### 4.1. Through hidden features and weights

`_apply_backward` accumulates gradients into `grad_h_in` and
`grad_tp_weights` given the pullback `grad_h_out`:

```c
double *grad_h_out = calloc(n_nodes * h_out->total_dim, sizeof(double));
double *grad_h_in = calloc(n_nodes * h_in->total_dim, sizeof(double));
double *grad_w = calloc(nw, sizeof(double));

/* Pre-zero grad_h_in and grad_w at the start of each training iteration;
 * _apply_backward accumulates (+=). */

irrep_nequip_layer_apply_backward(layer, w, n_nodes, n_edges,
 edge_src, edge_dst, edge_vec,
 h_in_buf, grad_h_out,
 grad_h_in, grad_w);
```

### 4.2. Through edge geometry

Needed for force evaluation in equivariant MD (LLG integrators,
atomistic dynamics). `_apply_forces` accumulates
`∂L / ∂edge_vec[e, axis]`:

```c
double *grad_edge_vec = calloc(n_edges * 3, sizeof(double));

irrep_nequip_layer_apply_forces(layer, w, n_nodes, n_edges,
 edge_src, edge_dst, edge_vec,
 h_in_buf, grad_h_out,
 grad_edge_vec);
```

The chain rule flows through three pieces: the spherical-harmonic
gradient (tangent to `S²`; see
[`tutorials/02_spherical_harmonics.md`](02_spherical_harmonics.md)),
the radial-basis derivative (with a Taylor fallback at `|ar| < 10⁻³`
to avoid catastrophic cancellation), and the cutoff derivative. All
three pieces are implemented and tested; the composite path is
finite-difference cross-checked to `10⁻⁶` on every edge × axis
combination in `tests/test_nequip.c`.

## 5. Equivariance check

The single most important property: if you rotate every input by the
same `R`, the output rotates by the same `R`. Concretely:

```c
/* 1. Rotate each node's hidden feature by D_{h_in}(R). */
/* 2. Rotate every edge_vec by R (cartesian matrix apply). */
/* 3. Apply the layer to the rotated inputs → h_out_rotated_in. */
/* 4. Apply the layer to the originals, then rotate h_out → h_out_direct. */
/* 5. Assert ‖ h_out_rotated_in - h_out_direct ‖_∞ < 1e-10. */
```

This is the pattern used in `tests/test_nequip.c` against random
`(α, β, γ)` triples; if this test fails, one of:

- your tensor-product path has the wrong `i^{l_a + l_b − l_c}` phase,
- a real↔complex basis change has drifted (rare — the library's
 `U` matrices are cached and tested),
- the cutoff is discontinuous at `r_cut` (which it isn't for the two
 shipped families),

or, most likely, a downstream re-implementation doesn't match. Libirrep
runs this test automatically on every PR.

## 6. Stacking layers

Between NequIP message-passing layers, you typically normalise,
mix channels linearly, and gate:

```c
#include <irrep/equivariant_layers.h>

/* RMS-norm per (l, parity) block — learnable per-(term, channel) scale. */
int channels = /* from your MLP scheme */;
double *scales = calloc(h_out->num_terms * channels, sizeof(double));
/* ... fill scales ... */
irrep_norm_rms(h_out, channels, scales, h_buf, h_normed);

/* Channel mixing within each irrep block, preserving (l, parity). */
irrep_linear_t *lin = irrep_linear_build(h_out, h_out, channels, channels);
int nw_lin = irrep_linear_num_weights(lin);
double *w_lin = calloc(nw_lin, sizeof(double));
irrep_linear_apply(lin, w_lin, h_normed, h_next);

/* Multiplicative gate — a separate scalar feature gates each block. */
double *gate_scalars = /* activated per-(term, channel) gate values */;
irrep_gate_apply(h_out, channels, gate_scalars, h_next, h_gated);
```

Each primitive commutes with the block-diagonal Wigner-D by
construction — stacking them preserves equivariance automatically, and
the standard decomposition

```
message-pass → normalise → linear-mix → gate → next block
```

is what NequIP, MACE, and Allegro share architecturally.

## 7. Symmetric-NQS projection (1.2)

For neural quantum states on lattices with a discrete point symmetry
(square C₄ᵥ, hexagonal / kagome D₆, triangular C₃ᵥ / D₃), projecting
the network output onto a specific irrep of the point group produces an
ansatz that lives in a well-defined symmetry sector. Libirrep ships
this as `irrep_pg_project`:

```c
#include <irrep/point_group.h>

irrep_pg_table_t *pg = irrep_pg_table_build(IRREP_PG_C4V);
irrep_multiset_t *spec = irrep_multiset_parse("2x0e + 1x1o + 1x2e");

double *h_projected = calloc(spec->total_dim, sizeof(double));
irrep_pg_project(pg, /*mu=A1*/ 0, spec, h_feature, h_projected);
/* h_projected lives in the totally-symmetric (A₁) subspace. */
```

The projector `P_μ = (d_μ / |G|) Σ_g χ_μ*(g) D(g)` satisfies
`P_μ² = P_μ`, `P_μ P_ν = 0 for μ ≠ ν`, and `Σ_μ P_μ = I` — all
verified numerically to `10⁻¹⁰`. As of 1.2 the `D^l(g)` matrices are
pre-computed at table-build time, so `_project` costs only an
element-sum of matrix-vector multiplies — ~80× faster than the previous
rebuild-per-call form.

See `examples/symmetric_nqs_projection.c` for an end-to-end
demonstration on a 4-node square lattice.

## 8. Training-loop topology

Libirrep does not ship an autograd engine. Every primitive exposes a
forward and a backward routine; the caller owns the topological order
in which they are called. The price is ~20 lines of glue per layer;
the benefit is framework-agnosticism — the same library binds into
PyTorch C++ extensions, JAX's FFI, a hand-rolled C trainer, or a
research prototype in Rust.

The canonical topological order for one training step:

```
(1) forward pass:
 for each layer L:
 call L.apply [writes cache of intermediates if needed]
(2) compute loss from final layer's output
(3) backward pass:
 zero grad_w, grad_h for every layer
 for each layer L in reverse:
 call L.apply_backward [accumulates into grad_w[L], grad_h[L-1]]
(4) optimizer step on accumulated grad_w
```

For force-aware training (molecular dynamics, LLG spin dynamics), add
`L.apply_forces` calls at the output layer and propagate the
`grad_edge_vec` through whatever geometry graph produced `edge_vec`
(typically an automatic-differentiation hook on the position tensor).

## 9. When to not use `irrep_nequip_layer_t`

The shipped layer covers the common case — single scalar-weight scheme,
Bessel radial basis, polynomial / cosine cutoff, fixed aggregation. Roll
your own when:

- You want a **learnable** radial head — e.g., the per-neighbour
 Chebyshev polynomial in DimeNet (Klicpera et al. 2020) or a small
 MLP on `r_ij`. Build the layer manually from
 `irrep_sph_harm_cart_all` + your radial net +
 `irrep_tp_apply_uvw`.
- You need an **edge-output** rather than aggregated-node-output
 layer — useful for edge attention or message-embedding tasks.
- You want **many-body correlation features** (MACE) — build the
 second-order tensor product yourself using `irrep_tp_apply_uvw` and
 the 6j / 9j recoupling symbols.
- You need to expose `D(R)` at a particular layer for integrating into
 an O(3)-invariant pooling — the block-diagonal
 `irrep_wigner_D_multiset` gives you the action directly.

In each case, libirrep provides the low-level machinery; the layer
composition is your orchestration.

## 10. References

- **Thomas, N., Smidt, T., Kearnes, S., et al.** "Tensor Field
 Networks," arXiv:1802.08219 (2018) — the foundational paper
 introducing the path-indexed SO(3) tensor product as the basis for
 equivariant neural networks.
- **Batzner, S., Musaelian, A., Sun, L., et al.** "E(3)-Equivariant
 Graph Neural Networks for Data-Efficient and Accurate Interatomic
 Potentials," *Nature Communications* **13**, 2453 (2022).
 DOI: [10.1038/s41467-022-29939-5](https://doi.org/10.1038/s41467-022-29939-5).
- **Batatia, I., Kovacs, D. P., Simm, G. N. C., et al.** "MACE:
 Higher Order Equivariant Message Passing Neural Networks for Fast
 and Accurate Force Fields," NeurIPS **35**, 11423 (2022).
- **Musaelian, A., Batzner, S., Johansson, A., et al.** "Learning
 local equivariant representations for large-scale atomistic
 dynamics," *Nature Communications* **14**, 579 (2023).
- **Klicpera, J., Groß, J., & Günnemann, S.** "Directional Message
 Passing for Molecular Graphs," ICLR (2020), arXiv:2003.03123.
- **Geiger, M., Smidt, T. E., Miller, B. K., et al.** "e3nn:
 Euclidean Neural Networks," arXiv:2207.09453 (2022).

Full bibliography in [`REFERENCES.md`](../REFERENCES.md).
