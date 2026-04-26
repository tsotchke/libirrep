# Migrating from Python `e3nn` to `libirrep`

A side-by-side cheat sheet for users fluent in the Python `e3nn` stack.
Where conventions differ, the difference is flagged in a **"Note"**.

## Irreps and multisets

| `e3nn` (Python) | `libirrep` (C) |
|---|---|
| `o3.Irreps("1x0e + 2x1o + 1x2e")` | `irrep_multiset_parse("1x0e + 2x1o + 1x2e")` |
| `irreps.dim` | `irrep_multiset_dim(m)` |
| `irreps.simplify()` | `irrep_multiset_simplify(m)` |
| `irreps1 + irreps2` | `irrep_multiset_direct_sum(m1, m2)` |
| string repr (`str(irreps)`) | `irrep_multiset_format(m, buf, n)` |

**Note.** `irrep_multiset_parse` rejects zero multiplicity, negative `l`,
and trailing `+` with an explicit parse error (available via
`irrep_last_error()`).

## Rotations

| `e3nn` | `libirrep` |
|---|---|
| `o3.rand_matrix()` | `irrep_quat_random(&state)` + `irrep_rot_from_quat(q)` |
| `o3.angles_to_matrix(α, β, γ)` | `irrep_rot_from_euler_zyz((irrep_euler_zyz_t){α, β, γ})` |
| `o3.matrix_to_angles(R)` | `irrep_euler_zyz_from_rot(R)` |
| `o3.compose(R1, R2)` | `irrep_rot_compose(R1, R2)` |

**Convention.** Both conventions use **ZYZ Euler angles** and **active**
rotations. Quaternion layout in `libirrep` is `{x, y, z, w}` with `w` last
(not `wxyz`); Shoemake-sampled quaternions are canonicalised to `w ≥ 0`.

## Spherical harmonics

| `e3nn` | `libirrep` |
|---|---|
| `o3.spherical_harmonics(l, x, normalize=True)` | `irrep_sph_harm_cart_all(l_max, out, r_hat)` |
| `o3.spherical_harmonics(ls, x, normalize=True)` | same — batched over `l = 0..l_max` |
| `o3.spherical_harmonics_alpha_beta(...)` | `irrep_sph_harm_real(l, m, theta, phi)` |

**Convention.** `libirrep` uses the Condon–Shortley phase and the e3nn sign
choice for real SH (`Y_{1,+1} ∝ +x`). Inputs to `_cart_all` must be unit
vectors — no normalisation is done internally.

## Wigner-D

| `e3nn` | `libirrep` |
|---|---|
| `o3.wigner_D(l, α, β, γ)` | `irrep_wigner_D_matrix(l, out, α, β, γ)` |
| `o3.wigner_D(l, R).to_real()` (rough analog) | `U_l D(R) U_l†` using `irrep_sph_harm_complex_to_real` |

`libirrep` exposes both the scalar `D^j_{m'm}(α,β,γ)` and the full
`(2j+1)²` matrix. Small-d `d^j(β)` is available as
`irrep_wigner_d_small_2j` (doubled-integer API, stable past `j = 50`).

## Tensor products

| `e3nn` | `libirrep` |
|---|---|
| `o3.TensorProduct(in1, in2, out, instructions)` | `irrep_tp_build(a, b, c, selected_paths, n)` |
| `tp(a, b)` | `irrep_tp_apply(desc, a, b, c_out)` |
| `tp(a, b, weights)` | `irrep_tp_apply_weighted(desc, w, a, b, c_out)` |

**v1.2 scope.** `libirrep` supports both the UUU channel mode (matched
multiplicities across `a`, `b`, `c`, one scalar per path — call
`irrep_tp_build` / `irrep_tp_apply_weighted`) and the UVW mode (independent
multiplicities, full `[W, V, U]` weight tensor per path — call
`irrep_tp_build_uvw` / `irrep_tp_apply_uvw`). Paths must satisfy the
parity rule `parity_a · parity_b = parity_c`, and **all** valid
`(l_a, l_b, l_c)` triples are supported — including odd-l-sum paths like
the cross-product `1o ⊗ 1o → 1e`. Those carry an
`i^{l_a + l_b - l_c}` phase correction internally so the real-basis
output is guaranteed real. See `examples/torque_net_tp_paths.c` for the
cross-product sign convention, verified bit-exact against cartesian
`a × b` to ~2e-16.

Differences from `e3nn` at parity: libirrep's `(1, 1, 1)` real-basis path
produces `+(1/√2) · (a × b)` in the `(y, z, x)` SH layout; e3nn's sign
and layout may differ depending on their version and build. Consumers
porting between the two should use the included
`tests/test_downstream_compat/` golden vectors to pin conventions.

## Equivariant neural-network layers

| `e3nn` | `libirrep` |
|---|---|
| `nn.Linear(in_irreps, out_irreps)` | `irrep_linear_build(in, out, 1, 1)` |
| `nn.BatchNorm(irreps)` (RMS-like) | `irrep_norm_rms(m, ch, scales, in, out)` |
| `nn.Gate(irreps_scalars, ...)` | `irrep_gate_apply(m, ch, gates, in, out)` |

`libirrep`'s layers are the minimal set needed to assemble a two-layer
equivariant MLP in pure C (see `examples/equivariant_mlp.c`). Channel
dimension parameters are accepted for future "uvw"-mode use but are
currently treated as `1`.

## Radial basis functions

| `e3nn` / NequIP | `libirrep` |
|---|---|
| `o3.spherical_bessel(n, x)` | `irrep_rbf_bessel(n, r, r_cut)` |
| `nequip.cutoffs.poly` | `irrep_cutoff_polynomial(r, r_cut, p)` |
| `nequip.cutoffs.cosine` | `irrep_cutoff_cosine(r, r_cut)` |

Analytic derivatives are provided for all cutoffs (`_d` suffix).

## Quadrature on the sphere

| `e3nn` | `libirrep` |
|---|---|
| (external `quadpy` / `lebedev_laikov`) | `irrep_lebedev_size`, `irrep_lebedev_fill` |
| (not provided) | `irrep_gauss_legendre(n, nodes, weights)` |

`libirrep` ships Lebedev rules for orders 3, 5, 7 (6/14/26 points).
Additional orders arrive as the Lebedev–Laikov 1999 data is imported.

## Autograd

`libirrep` does not embed an autograd engine. Every kernel that has a
natural backward pass exposes one directly:
`irrep_tp_apply_backward{,_weighted,_batch}`, `irrep_linear_backward`,
`irrep_norm_rms_backward`. Callers stitch them into whatever autograd
framework they target — plain C function pointers into PyTorch custom
ops, MLIR, or a hand-rolled reverse-mode engine.

## Precision

- `e3nn` defaults to float32; `libirrep` defaults to float64.
- `_f32` suffixed variants (currently SH cartesian only) wrap the double
 kernel and cast. A fully-SIMD float path is scheduled for 1.3; the 1.2
 NEON kernels are double-precision, bit-exact against the scalar path.

## What's new in 1.1 / 1.2

### NequIP message-passing layer (1.1 →)

e3nn callers stitch `o3.spherical_harmonics`, `o3.FullyConnectedTensorProduct`,
a radial MLP, and a scatter-add themselves. libirrep ships this as
`irrep_nequip_layer_t`:

| e3nn | libirrep |
| --------------------------------------------------------------------------- | -------------------------------------------------------------- |
| `e3nn.nn.models.gate_points_2101.Network` (forward pass equivalent) | `irrep_nequip_layer_build` + `irrep_nequip_layer_apply` |
| `torch.autograd.backward` on the layer output | `irrep_nequip_layer_apply_backward` (accumulates into `grad_h_in`, `grad_tp_weights`) |
| `edge_vec.grad` after backward | `irrep_nequip_layer_apply_forces` (separate call, accumulates into `grad_edge_vec`) |
| `e3nn.o3.Irreps("2x0e + 1x1o").randn(...)` | `irrep_multiset_parse("2x0e + 1x1o")` + manual fill |
| `IrrepsFromSpec("... -> ... [sh=2]")`-like convenience | **1.2:** `irrep_nequip_layer_from_spec("2x0e + 1x1o -> 1x1o [sh=2, radial=8, r_cut=1.5]")` |

### Per-path L2 regulariser (1.2)

For SR-style or natural-gradient training where different tensor-product
paths sit on differently-curved sections of the variational manifold and
want per-path regularisation strengths:

```c
double *per_path = malloc(irrep_tp_num_paths(desc) * sizeof(double));
irrep_tp_weight_l2_per_path_uvw(desc, w, per_path);

/* backward: accumulates grad_w += 2·λ_k·w[k,w,v,u] with λ_k per path */
irrep_tp_weight_l2_per_path_uvw_backward(desc, w, lambdas, grad_w);
```

e3nn has no direct analogue; the usual pattern is `sum(w**2)` in Python
with a scalar coefficient.

### Point-group projection (1.2)

e3nn users build permutation-symmetric ansätze via the `group_projection`
helpers in the experimental tree. libirrep ships:

| e3nn experimental | libirrep |
| ---------------------------------------------------------------- | ----------------------------------------------------------- |
| `from e3nn.nn.so3 import SymmetricTensorProduct` | `irrep_pg_project(table, mu, spec, in, out)` |
| `group_projection(..., "A1", ...)` | `irrep_pg_project(table, /*mu=A1*/ 0, spec, in, out)` |
| Decomposition of a representation under the group | `irrep_pg_reduce(table, spec, out_mult)` (Bradley–Cracknell) |

Supported groups in 1.2: **C₄ᵥ** (square lattice, 5 irreps) and **D₆**
(hexagonal / kagome, 6 irreps), C₃ᵥ (triangular, 3 irreps), and D₃
(triangular, proper only, 3 irreps).

### Analytic radial derivatives (1.2)

`irrep_rbf_bessel_d` (plus `_all` and `_batch` variants) closes the chain
rule for force evaluation. Taylor-expanded at small r so the naïve
`(a cos(ar)/r − sin(ar)/r²)` cancellation doesn't cost precision. e3nn
equivalent: manual autograd through `torch.sin / r`.

### Runtime SIMD dispatch (1.2)

Batched radial / SH entry points (`_batch` suffix) route through a
runtime-selected kernel — NEON on aarch64 today, AVX2 / NEON-SH in 1.3.
No API change; the dispatch is transparent.

### Half-integer spin — `multiset_2j.h`

e3nn's `Irreps` type permits labels like `"1/2e"` (spin-½, even
parity). Libirrep's integer-l `irrep_multiset_t` cannot express these;
use `irrep_multiset_2j_t` from `irrep/multiset_2j.h` instead:

| e3nn | libirrep |
| ---------------------------------------------------------------- | -------------------------------------------------------------- |
| `o3.Irreps("1x1/2e")` | `irrep_multiset_2j_parse("1x1/2e")` |
| `o3.Irreps("1x0e + 1x1/2o")` | `irrep_multiset_2j_parse("1x0e + 1x1/2o")` |
| `irreps.dim` | `.total_dim` (equivalent) |
| n/a — must iterate `irreps.ls` | `irrep_multiset_2j_has_half_integer(m)` |
| Kramers-degeneracy check on the representation | `irrep_time_reversal_square_sign_2j(m)` returns `−1` iff any block has odd `two_j` |

Internally, labels store `two_j = 2·j`, so spin-½ has `two_j = 1` —
consistent with the `_2j`-suffixed function family throughout the
library (`irrep_cg_2j`, `irrep_wigner_d_small_2j`, etc.).

## What's new in 1.3 (no e3nn equivalent)

The 1.3.0-alpha cycle added entire categories of functionality that
have no analog in `e3nn`. e3nn is fundamentally an equivariant-neural-
network library; libirrep's 1.3 substrate adds **physics primitives**
(condensed matter, exact diagonalisation) and **materials-search
analyzers** (symmetry-allowed exchange tensor decomposition) that
sit alongside the e3nn-style core.

### Spin Hamiltonians, ED, and entanglement

| capability | libirrep entry point | e3nn equivalent |
|---|---|---|
| Spin-½ Heisenberg apply operator | `irrep_heisenberg_new` + `_apply` | none |
| J₁-J₂ Heisenberg apply | `irrep_heisenberg_j1j2_new` | none |
| Sparse Lanczos eigensolver (3-vector / reorth / +eigvecs) | `irrep_lanczos_eigvals(_reorth)(_eigvecs_reorth)` in `rdm.h` | none — would call out to scipy |
| Reduced density matrix of a sub-system | `irrep_partial_trace` | none |
| Von Neumann / Rényi entropies | `irrep_entropy_vonneumann`, `_renyi` | none |
| Kitaev–Preskill γ topological diagnostic | `irrep_topological_entanglement_entropy` | none |
| Total-J projection on spin-½ chains | `irrep_spin_project_spin_half` | none |

The use case here is **substrate validation** for symmetric NQS and
neural-quantum-state ansätze — verifying that a small-cluster ground
state matches expected symmetry properties before trusting a 108-site
NQS prediction.

### 2D and 3D Bravais lattices

| capability | libirrep | e3nn |
|---|---|---|
| 2D Bravais (square, triangular, honeycomb, kagome) | `irrep_lattice_build` | none — manual graph construction |
| 3D Bravais (SC, BCC, FCC, Diamond, Pyrochlore) | `irrep_lattice3d_build` | none |
| Periodic boundary conditions, NN/NNN bond enumeration | `_num_bonds_nn(_nnn)` + `_fill_*` | none |
| Brillouin zone k-point grid | `irrep_lattice_k_grid`, `irrep_lattice3d_k_grid` | none |
| 2D wallpaper-group site permutations | `irrep_space_group_build` | none |

For molecule-style equivariant-NN consumers (no PBC, no extended
crystal), this layer is irrelevant. For **crystal-property prediction**
(magnetism, phonons, transport), it's the foundation.

### Bond-exchange-tensor + chirality symmetry analyzer (`dmi.h`)

The cleanest libirrep-unique 1.3 deliverable. Given a candidate
crystal symmetry (with optional time-reversal augmentation for
magnetic point groups), return the complete symmetry-allowed
exchange-tensor structure of the bilinear and trilinear spin
couplings:

| capability | libirrep | e3nn / equivalent |
|---|---|---|
| Moriya's five DMI rules in one projector | `irrep_dmi_allowed_basis(_from_pg)` | none — historically hand-derived from International Tables vol. A |
| Symmetric exchange (Heisenberg + Kitaev-Γ-type anisotropy) | `irrep_exchange_symmetric_basis(_from_pg)` | none |
| Three-spin scalar chirality `S_i · (S_j × S_k)` symmetry verdict | `irrep_chirality_allowed(_from_pg)` | none |
| Magnetic-point-group antiunitary `T·g` operations | `irrep_dmi_sym_op_t.antiunitary` flag | none |

Materials-search workflow: propose a candidate space group + magnetic
ordering → run the analyzer → get the parameter scaffolding for a
downstream DFT or micromagnetic simulation. The pyrochlore "up
tetrahedron" catalog example (`pyrochlore_tetra_complete_catalog.c`)
produces 72 group-theoretic verdicts in a single program — the
parameter scaffold for Yb₂Ti₂O₇-style quantum-spin-ice modelling.
See `docs/tutorials/08_materials_search.md`.

### Cubic point groups (T_d, O_h, O)

The 1.2 point-group support was 2D-lattice-oriented (C₄ᵥ, D₆,
C₃ᵥ, D₃). 1.3 adds the three cubic groups needed for diamond /
zincblende / SC / BCC / FCC / perovskite / chiral-cubic site
symmetry. e3nn handles point-group projection via its tensor-product
+ representation-theory machinery, so the closest equivalent is
"compose the projector as a sum of representation matrices weighted
by characters" — but the cubic character tables and conjugacy-class
structure aren't pre-tabulated in e3nn. libirrep ships them as
data plus the `irrep_pg_project` / `_reduce` API.

## ABI stability

`libirrep` tracks a stable C ABI within a major version. Every release
publishes `ABI_HASH = SHA-256(sorted exported symbols)` under
`release/<version>/`. If you vendor the static library, you can guard
against binary drift with

```c
#include <irrep/version.h>
/* assert irrep_abi_hash() matches the hash your build was linked against */
```
