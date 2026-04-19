# Tutorial 05 — Tensor Products

This is where irreps earn their keep. Two irrep-structured feature vectors
get combined into a third, respecting every symmetry (rotation, parity,
time reversal) that we started with.

## 1. What a tensor product actually computes

Given two input feature vectors

- `a` living in irrep space `A = ⊕_i mult^a_i · (l^a_i, p^a_i)`
- `b` living in irrep space `B = ⊕_j mult^b_j · (l^b_j, p^b_j)`

and a target space

- `c` living in `C = ⊕_k mult^c_k · (l^c_k, p^c_k)`

a "path" `(i, j, k)` is an allowed coupling of the i-th term of `a` with
the j-th term of `b` into the k-th term of `c`. Allowed means:

1. Triangle: `|l^a_i − l^b_j| ≤ l^c_k ≤ l^a_i + l^b_j`.
2. Parity: `p^a_i · p^b_j = p^c_k`.

Both even and odd `l_a + l_b + l_c` are supported; odd-sum paths (e.g. the
cross-product `(1, 1, 1)`) carry an `i^{l_a + l_b − l_c}` phase factor
internally so the real-basis output is guaranteed real.

For each allowed path and each channel index `u`, the tensor product
computes

```
c[k, u, m_c] += w_path · Σ_{m_a, m_b} W[n_a, n_b, n_c]
 · a[i, u, n_a] · b[j, u, n_b]
```

where `W` is a pre-computed real-basis coupling tensor (the real-basis
analogue of Clebsch-Gordan), and `m_a, m_b, m_c` are the m-component
indices of each irrep block.

## 2. Build / apply API

```c
#include <irrep/multiset.h>
#include <irrep/tensor_product.h>

irrep_multiset_t *a = irrep_multiset_parse("1x1o");
irrep_multiset_t *b = irrep_multiset_parse("1x1o");
irrep_multiset_t *c = irrep_multiset_parse("1x0e + 1x2e");

int n_paths = irrep_tp_enumerate_paths(a, b, c, NULL, 0);
int *paths = malloc(n_paths * 3 * sizeof(int));
irrep_tp_enumerate_paths(a, b, c, paths, n_paths);

tp_descriptor_t *desc = irrep_tp_build(a, b, c, paths, n_paths);
/* n_paths == 2: (1o ⊗ 1o → 0e) and (1o ⊗ 1o → 2e) */

double a_vec[3] = { 1, 0, 0 }; /* m = -1, 0, +1 (y, z, x) */
double b_vec[3] = { 0, 1, 0 };
double c_vec[6]; /* 0e dim 1 + 2e dim 5 */

irrep_tp_apply(desc, a_vec, b_vec, c_vec);
```

The real-basis scalar (`l=0`) result for two l=1 vectors is
`c[0] = (a · b) / √3`. For `a = (1, 0, 0)` and `b = (0, 1, 0)` that's 0 —
they are orthogonal. (Sign empirically verified against cartesian
`a · b` via `examples/torque_net_tp_paths.c` to machine precision.)

The `l=2` output (entries 1..5) carries the traceless symmetric part
of the outer product, again normalised to the real-basis spherical
harmonics convention.

## 3. Weighted tensor products

One scalar weight per path lets you parameterise the TP:

```c
double weights[2] = { 1.0, -0.5 };
irrep_tp_apply_weighted(desc, weights, a_vec, b_vec, c_vec);
```

Weights multiply path-by-path and therefore preserve equivariance — they
cannot mix paths of different `(l, parity)`.

## 4. Gradients for training

```c
double grad_c[6] = { ... }; /* pulled back from the loss */
double grad_a[3] = { 0 };
double grad_b[3] = { 0 };
double grad_w[2] = { 0 };

irrep_tp_apply_backward_weighted(desc, weights, a_vec, b_vec, grad_c,
 grad_a, grad_b, grad_w);
```

All three output buffers accumulate; zero them before the call if you
are not combining with another accumulation.

## 5. Batching

The batched variants stride by `total_dim` per batch slot:

```c
size_t B = 1024;
double *a_batch = malloc(B * 3 * sizeof(double)); /* a.total_dim = 3 */
double *b_batch = malloc(B * 3 * sizeof(double));
double *c_batch = malloc(B * 6 * sizeof(double));
/* ... fill a_batch, b_batch ... */
irrep_tp_apply_batch(desc, B, a_batch, b_batch, c_batch);
```

The batch dimension is outside the multiset, so every feature vector has
the same irrep structure. (Heterogeneous batches require multiple
descriptors.)

## 6. Verifying equivariance

You can (and should) test a new descriptor with a random rotation:

```c
/* 1. Build the real-basis rotation D_R = U · D_complex(α, β, γ) · U†
 * for each l appearing in a, b, c.
 * 2. Check that tp(D_A · a, D_B · b) == D_C · tp(a, b). */
```

The M7 test suite in `tests/test_tensor_product.c` has a worked example.

## 7. Worked example — mapping a hand-rolled Cartesian basis to TP paths

Spin-based torque networks, many effective-Hamiltonian terms in magnetism,
and several attention mechanisms in geometric GNNs are often written using
an explicit Cartesian basis — combinations of dot products, cross products,
and scalar triples of the inputs. Every one of those constructs factorises
into a specific `(l_a, l_b, l_c)` tensor-product path. Here is the mapping
for the five basis terms that show up in a typical equivariant torque net:

```
 T1 = (m_j · r̂) · m_i — scalar projection times m_i (l=1)
 T2 = m_j × r̂ — axial cross product (l=1)
 T3 = m_i × m_j — pair axial cross product (l=1)
 T4 = (m_i · m_j) · m_i — pair scalar × m_i (l=1)
 T5 = m_j — passthrough (l=1)
```

The rule of thumb:
- **dot product** of two l=1 vectors → the `l=0` output path of `1o ⊗ 1o`
 (parity `odd · odd = even`, so the target is `0e`)
- **cross product** of two l=1 vectors → the `l=1, even-parity` output
 path of `1o ⊗ 1o` (target `1e`, the axial pseudovector). `1o ⊗ 1o → 1o`
 is parity-forbidden and has zero paths.
- **scalar × vector** → the `l=1` output path of `0e ⊗ 1o`

Concretely, each term corresponds to one TP path with the inputs shown:

Parity is multiplicative under the tensor product: `p_a · p_b = p_c`. Two
polar (odd-parity) vectors produce an even output — their cross product is
an axial pseudovector (`1e`), not a polar vector. Libirrep's parity filter
enforces this, so `1o ⊗ 1o → 1o` has zero paths; the correct mapping uses
`1e` for any axial intermediate.

| Term | `A` (irreps) | `B` (irreps) | `C` (irreps) | Path `(l_a, l_b, l_c)` |
| ---- | ------------ | ------------ | ------------ | ---------------------- |
| `T1` | `m_j: 1x1o` | `r̂: 1x1o` | `tmp: 1x0e` × `m_i: 1x1o` → `1x1o` | two-step: `(1,1,0)` then `(0,1,1)` |
| `T2` | `m_j: 1x1o` | `r̂: 1x1o` | `1x1e` (axial, parity = odd·odd = even) | `(1,1,1)` |
| `T3` | `m_i: 1x1o` | `m_j: 1x1o` | `1x1e` (axial) | `(1,1,1)` |
| `T4` | `m_i: 1x1o` | `m_j: 1x1o` | `tmp: 1x0e` × `m_i: 1x1o` → `1x1o` | two-step: `(1,1,0)` then `(0,1,1)` |
| `T5` | `m_j: 1x1o` | — | `1x1o` | identity (no TP) |

If your physics has `m` as axial (a magnetic moment / pseudovector), swap
the inputs to `1x1e` and adjust the output parities accordingly —
axial · polar = parity-odd scalar, axial × polar = polar, etc. See
`examples/torque_net_tp_paths.c` for a worked end-to-end check.

A full UVW descriptor for `T2` (cross product `m_j × r̂` landing in an `l=1`
channel) builds like this:

```c
#include <irrep/multiset.h>
#include <irrep/tensor_product.h>

irrep_multiset_t *mj = irrep_multiset_parse("1x1o"); /* spin-1, odd parity */
irrep_multiset_t *rhat= irrep_multiset_parse("1x1o");
irrep_multiset_t *out = irrep_multiset_parse("1x1o");

int path[3] = { 0, 0, 0 }; /* the one (l_a, l_b, l_c) = (1, 1, 1) option */
tp_descriptor_t *tp = irrep_tp_build_uvw(mj, rhat, out, path, /*num_paths=*/1);

int nw = irrep_tp_num_weights_uvw(tp); /* 1·1·1 = 1 scalar */
double w[1] = { 1.0 }; /* bare cross product */
double c[3];
irrep_tp_apply_uvw(tp, w, mj_buf, rhat_buf, c);
/* c is m_j × r̂ up to an overall sign / √3 — see note below. */
```

### Sign and normalisation note

The real-basis `(1, 1, 1)` path in libirrep, with the `i^{l_a + l_b - l_c}`
phase convention, the real↔complex SH basis change, and the Clebsch-Gordan
`⟨1 m_a; 1 m_b | 1 m_c⟩` coefficients, produces `+(1/√2) · (a × b)` in the
`(y, z, x)` real-basis layout. Concretely:

```
irrep_tp_apply_uvw( {1o ⊗ 1o → 1e}, a, b ) = (1 / √2) · (a ×_real b)
```

where `a_real = (a_y, a_z, a_x)` is the cartesian-to-real-SH permutation
libirrep uses internally for l=1. To recover the bare cartesian cross
product, multiply by `√2` and then permute `(y, z, x) → (x, y, z)` on the
way out. See `examples/torque_net_tp_paths.c` for the full round-trip,
bit-exact to machine precision (residual ~2e-16).

The `l = 0` path likewise picks up the `1 / √3` from
`⟨1 m_a; 1 m_b | 0 0⟩`, so:

```
irrep_tp_apply_uvw( {1o ⊗ 1o → 0e}, a, b ) = (1 / √3) · (a · b)
```

Multiply by `√3` for the bare dot. No basis permutation needed for l=0
(it's a scalar).

### Why this matters

Once the five Cartesian terms are expressed as TP paths, you can replace
the hand-rolled basis with a small `irrep_multiset_t + tp_descriptor_t`
pair and let the equivariance check in `tests/test_nequip.c` cover the
whole thing — no more manual sign hunting for new terms. The torque-net's
existing `test_torque_net_irrep.c` cross-check (SH addition theorem plus
explicit Y₂⁰ spot check) already pins the SH side; porting the
cross/dot primitives to `irrep_tp_apply_uvw` closes the last hand-rolled
loop in that pipeline.

## 8. What's next

- Higher-order couplings via Wigner 6j / 9j (`irrep/recoupling.h`) —
 lets you short-circuit three-body tensor products.
- Gating, RMS-norm, and linear mixing before the next TP
 (`irrep/equivariant_layers.h`) — see tutorial 06.
- Building a full NequIP-style message pass in pure C — see
 `examples/nequip_message.c`.
