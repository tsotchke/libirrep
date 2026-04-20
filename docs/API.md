# API reference

This file is the narrative entry point into the public C API. Per-symbol
documentation — argument semantics, preconditions, algorithmic choices,
numerical tolerances — is in the Doxygen blocks of the individual
headers under `include/irrep/`; `make docs` renders them to HTML under
`build/docs/`.

Consumers who just want a quick index should scan the module table
below. Consumers who want to understand the semantics without opening
the headers can read the per-module notes. For the mathematical
definitions of every returned quantity, see
[`PHYSICS_APPENDIX.md`](PHYSICS_APPENDIX.md); for the primary sources of
every formula, [`REFERENCES.md`](REFERENCES.md).

---

## Module index

| Header | Module | Primary types | Entry points |
| ------ | ------ | ------------- | ------------ |
| `<irrep/types.h>` | core types | `irrep_quaternion_t`, `irrep_rot_matrix_t`, `irrep_euler_zyz_t`, `irrep_axis_angle_t`, `irrep_label_t`, `irrep_multiset_t`, `irrep_status_t` | `irrep_strerror`, `irrep_last_error` |
| `<irrep/version.h>` | versioning | — | `IRREP_VERSION_*` macros, `irrep_version_*()`, `irrep_abi_hash()` |
| `<irrep/export.h>` | symbol visibility | — | `IRREP_API` macro |
| `<irrep/simd.h>` | CPU features | `irrep_cpu_features_t` | `irrep_cpu_features()`, `irrep_cpu_has_neon()`, `…_has_avx2()` |
| `<irrep/so3.h>` | rotation math | (types only) | conversions, `_compose`, `_inverse`, `_rot_exp`, `_rot_log`, `_quat_slerp`, `_quat_random`, `_quat_frechet_mean` |
| `<irrep/su2.h>` | SU(2) ops | `irrep_su2_t` | Pauli matrices, SU(2) ↔ quaternion ↔ rotation matrix |
| `<irrep/multiset.h>` | integer-l multiset algebra | — | `_new`, `_free`, `_parse`, `_format`, `_simplify`, `_direct_sum`, `_append` |
| `<irrep/multiset_2j.h>` | doubled-integer multiset (half-integer-spin) | `irrep_multiset_2j_t`, `irrep_label_2j_t` | `_2j_new`, `_2j_free`, `_2j_parse`, `_2j_format`, `_2j_append`, `_2j_dim`, `_2j_has_half_integer`, `irrep_time_reversal_square_sign_2j` |
| `<irrep/parity.h>` | parity | — | `_parity`, `_parity_product`, `_parity_filter_paths` |
| `<irrep/time_reversal.h>` | time reversal | — | `_time_reversal_integer`, `_time_reversal_half_integer`, `_time_reversal_multiset`, `_time_reversal_square_sign` |
| `<irrep/clebsch_gordan.h>` | CG coefficients | `cg_table_t` | `irrep_cg`, `irrep_cg_2j`, `irrep_wigner_3j`, `irrep_wigner_3j_2j`, `irrep_cg_table_build`, `irrep_cg_lookup`, `irrep_cg_block` |
| `<irrep/recoupling.h>` | 6j / 9j / Racah W | — | `irrep_wigner_6j`, `irrep_wigner_9j`, `irrep_racah_w` |
| `<irrep/spherical_harmonics.h>` | spherical harmonics | — | `irrep_sph_harm` (complex), `_real`, `_cart`, `_cart_all`, `_cart_grad`, `_cart_all_batch`, `_cart_all_grad_batch`, `_complex_to_real` |
| `<irrep/solid_harmonics.h>` | solid harmonics | — | `irrep_solid_harm_cart`, `_cart_grad`, `_cart_all` |
| `<irrep/wigner_d.h>` | Wigner-D | — | `irrep_wigner_d_small`, `_D`, `_D_matrix`, `_d_matrix`, `_D_multiset`, `_d_small_dbeta` |
| `<irrep/tensor_product.h>` | e3nn-style tensor products | `tp_descriptor_t`, `irrep_tp_mode_t` | `irrep_tp_build`, `_build_uvw`, `_apply`, `_apply_weighted`, `_apply_uvw`, `_apply_backward`, `_apply_uvw_backward`, `_weight_l2_per_path_uvw` |
| `<irrep/radial.h>` | radial basis | — | `irrep_rbf_bessel`, `_bessel_all`, `_bessel_d`, `_bessel_d_all`, `_gaussian`, `_gaussian_grid`, `irrep_cutoff_cosine`, `_cosine_d`, `_polynomial`, `_polynomial_d`, plus `_batch` variants |
| `<irrep/quadrature.h>` | quadrature | — | `irrep_lebedev_size`, `_fill`, `_gauss_legendre`, `_quadrature_sphere_size`, `_fill` |
| `<irrep/equivariant_layers.h>` | NN primitives | `irrep_linear_t` | `irrep_linear_build`, `_apply`, `_backward`, `irrep_norm_rms`, `_backward`, `irrep_gate_apply` |
| `<irrep/nequip.h>` | NequIP message-passing layer | `irrep_nequip_layer_t`, `irrep_nequip_cutoff_t` | `_build`, `_from_spec`, `_apply`, `_apply_backward`, `_apply_forces`, `_num_weights`, `_free` |
| `<irrep/point_group.h>` | point-group projection | `irrep_pg_table_t`, `irrep_point_group_t` | `_table_build`, `_table_free`, `_num_irreps`, `_order`, `_irrep_label`, `_project`, `_reduce` |
| `<irrep/lattice.h>` (1.3) | 2D lattices with PBC | `irrep_lattice_t`, `irrep_lattice_kind_t` | `_build`, `_free`, `_num_sites`, `_num_cells`, `_primitive_vectors`, `_reciprocal_vectors`, `_site_position`, `_sublattice_of`, `_cell_of`, `_site_index`, `_translate`, `_num_bonds_nn`, `_num_bonds_nnn`, `_fill_bonds_nn`, `_fill_bonds_nnn`, `_k_grid` |
| `<irrep/space_group.h>` (1.3) | 2D wallpaper-group tables | `irrep_space_group_t`, `irrep_wallpaper_t` | `_build`, `_free`, `_order`, `_point_order`, `_num_sites`, `_kind`, `_lattice`, `_apply`, `_permutation`, `_permutation_inverse`, `_apply_config` |
| `<irrep/config_project.h>` (1.3) | configuration-space projection | `irrep_sg_irrep_t` | `irrep_sg_irrep_new`, `_free`, `irrep_sg_trivial`, `irrep_sg_sign_rep`, `irrep_sg_project_amplitude`, `irrep_sg_project_A1`, `irrep_sg_enumerate_orbit`, `irrep_sg_adapted_basis`, `irrep_sg_bloch_amplitude`, `irrep_sg_bloch_basis` |
| `<irrep/rdm.h>` (1.3) | reduced density matrix / entropies | — | `irrep_partial_trace`, `irrep_hermitian_eigvals`, `irrep_entropy_vonneumann_spectrum`, `_renyi_spectrum`, `_vonneumann`, `_renyi`, `irrep_topological_entanglement_entropy` |
| `<irrep/sym_group.h>` (1.3) | S_N, Young tableaux, (anti)sym projectors | — | `irrep_factorial`, `_permutation_sign`, `_permutations_all`, `_young_dim`, `irrep_sym_group_antisymmetrize`, `_symmetrize` |
| `<irrep/spin_project.h>` (1.3) | total-J projection on spin-½ chains | — | `irrep_spin_half_apply_rotation`, `irrep_spin_project_spin_half` |
| `<irrep/tensor_product.h>` (half-int path, 1.3) | spinor tensor products | `tp_2j_descriptor_t` | `irrep_tp_2j_enumerate_paths`, `_build`, `_free`, `_apply`, `_apply_weighted`, `_apply_backward`, `_output_dim`, `_num_paths` |
| `<irrep/irrep.h>` | umbrella | — | all of the above |

---

## Calling-convention conventions

**Suffix grammar** — uniform across the API.

- `_2j` — the function takes doubled-integer arguments. Every label
 (`j`, `m`, `J`, `M`) is passed as `2 ·` its physical value, so spin-½
 enters as `two_j = 1` and there is no ambiguity about what a
 half-integer means in C's `int` type.
- `_f32` — single-precision wrapper. Internally casts to `double` and
 back; a true float SIMD path is a 1.3 deliverable.
- `_batch` — batched variant. Accepts an extra leading `size_t N`
 argument and writes `N` copies of the single-call output to a
 caller-provided contiguous buffer. Routes through the runtime SIMD
 dispatch table.
- `_backward` — backward pass for a forward routine. Accumulates into
 the output-gradient buffers; the caller must pre-zero them at the
 start of each training iteration.

**Error handling** — one of three patterns per function.

1. **Pure math functions** (`irrep_cg`, `irrep_sph_harm_real`,
 `irrep_wigner_d_small`, `irrep_rot_apply`, …) return the
 mathematically correct value on valid input. On invalid input —
 including every selection-rule violation — the return is `0.0`
 for scalar-valued functions and a zero-filled buffer for matrix
 and array outputs, which coincides with the mathematical
 definition in both cases. Consequently, callers may sum over
 full index ranges without pre-filtering forbidden terms.
2. **Builder functions** (`irrep_multiset_parse`, `irrep_tp_build`,
 `irrep_nequip_layer_build`, `irrep_pg_table_build`, …) return
 `NULL` on failure. A thread-local error message is available via
 `irrep_last_error()`; it is stable until the next failing call on
 the same thread.
3. **Validating functions** (currently `irrep_multiset_append`, and
 the reserved `_from_spec_ex` slot in 1.3) return an
 `irrep_status_t` enum from `<irrep/types.h>`.

**Memory ownership** — every allocating function is paired with a
`_free` routine that releases its output. No public API hands back a
`malloc`-ed pointer that the caller must `free()` with libc. Opaque
handle types (`tp_descriptor_t`, `cg_table_t`, `irrep_linear_t`,
`irrep_nequip_layer_t`, `irrep_pg_table_t`) are released via their
respective `_free` calls.

**Thread safety** — see [`DESIGN.md`](DESIGN.md) §4. The one-sentence
summary: pure math is reentrant, built tables are concurrent-read-safe,
`irrep_last_error()` is thread-local, and the SIMD dispatch table is
race-free under the C11 memory model.

---

## Per-module notes

### `types.h` — the cross-cutting types

Rotation types, the `(l, p)` irrep label, and the `irrep_multiset_t`
direct-sum container appear here. Every other header includes this one.
Two library-wide constants:

- `IRREP_L_MAX = 16` — the maximum `l` supported by the spherical- and
 solid-harmonic routines, and the limit against which internal stack
 buffers are sized.
- `IRREP_TWO_J_MAX = 32` — the maximum `2j` for Clebsch-Gordan and
 Wigner-D (so `j_max = 16`).

### `so3.h` — SO(3) rotation math

Every pairwise conversion between `irrep_rot_matrix_t`,
`irrep_quaternion_t`, `irrep_euler_zyz_t`, and `irrep_axis_angle_t`,
tested to round-trip to `10⁻¹²`. Rodrigues exponential with a first-
order Taylor fallback for `|ω| < 10⁻¹²`; Markley branch-switching
logarithm stable near `θ = π`; Shoemake uniform sampling with `w ≥ 0`
canonicalisation; SLERP; iterated Karcher-Fréchet mean on a weighted
set of unit quaternions. `irrep_rot_apply` and `irrep_quat_apply` act a
rotation on a 3-vector in-place.

Numerical gotchas are named with macros and documented inline:
`IRREP_EULER_GIMBAL_EPS = 10⁻¹²`, `IRREP_ROT_EXP_SMALL = 10⁻¹²`. See
[`PHYSICS_APPENDIX.md`](PHYSICS_APPENDIX.md) §12.3 for the derivation
of each threshold.

### `su2.h` — SU(2), spin-½, double cover

Pauli matrices in `irrep_pauli_x`, `_y`, `_z`; SU(2) ↔ quaternion
conversion (`irrep_su2_from_quat`, `_quat_from_su2`); the explicit
SU(2) → SO(3) double cover in `irrep_rot_from_su2`. Compose, inverse,
matrix exponential of a 2×2 complex generator, and apply-to-state are
the group operations.

### `multiset.h` — irrep algebra

The string grammar `"1x0e + 2x1o + 1x2e"` is parsed by
`irrep_multiset_parse`; formatted back via `_format`. `_simplify` merges
like terms and sorts canonically (by `l` ascending, even-before-odd).
`_direct_sum` concatenates. The underlying type is a pair of dynamic
arrays `(labels, multiplicities)` with a running `total_dim =
Σ mult_i · (2 l_i + 1)` maintained by the API.

### `parity.h` — O(3) parity

`_parity` extracts the parity bit of a single `irrep_label_t`;
`_parity_product` computes `p_c = p_a · p_b` for a tensor-product
output; `_parity_filter_paths` drops parity-forbidden entries from a
flat path list in place.

### `time_reversal.h` — antiunitary T operator

Integer-l and half-integer-j forms of the linear part of
`T = U_T · K`. The multiset-level `_time_reversal_multiset` writes a
`total_dim × total_dim` block-diagonal `U_T`. The sign check
`_time_reversal_square_sign` returns `+1` if every block is integer-l
and `−1` if any is half-integer — i.e., whether the input space has
Kramers degeneracy.

### `clebsch_gordan.h` — the CG implementation

Integer signatures `irrep_cg(j1, m1, j2, m2, J, M)`; half-integer
`_cg_2j(two_j1, two_m1, …)`. Wigner 3j symbols follow the relation
derived in [`PHYSICS_APPENDIX.md`](PHYSICS_APPENDIX.md) §7.3. A cached
table `cg_table_t` is available for batch lookups; the internal
representation is a flat array indexed by
`(j1, m1, j2, m2, J, M)` offsets. `_cg_block` fills a dense
`(2j1+1)(2j2+1)(2J+1)` block for kernel-style use.

### `recoupling.h` — 6j, 9j, Racah W

Direct implementations of the single-sum formulae in
[`PHYSICS_APPENDIX.md`](PHYSICS_APPENDIX.md) §8, with log-gamma
stabilisation. Selection-rule violations return `0.0`.

### `spherical_harmonics.h` — three API surfaces

- Polar complex: `irrep_sph_harm(l, m, theta, phi)`.
- Polar real: `irrep_sph_harm_real(l, m, theta, phi)`.
- Cartesian real (the hot path): `irrep_sph_harm_cart(l, out, r_hat)`
 writes `2l + 1` values; `_cart_all` writes `(l_max+1)²` values for
 the whole triangle `0 ≤ l ≤ l_max`; `_cart_all_batch` is the
 SIMD-dispatched batched form for many edges at once.

Gradients: `_cart_grad` per-l, `_cart_all_grad_batch` over the full
triangle and many edges. The complex-to-real basis change matrix is
available via `_complex_to_real`.

### `solid_harmonics.h` — homogeneous-polynomial form

`R_{l, m}(r) = |r|^l · Y^{real}_{l, m}(r / |r|)`. Explicit polynomial
forms for `l ≤ 4`, a Limpanuparb-Milthorpe-style recurrence with
co-evolved gradients for larger `l`.

### `wigner_d.h` — rotations on irreps

Scalar small-d (`_wigner_d_small`, half-integer `_2j`), full-D scalar
(`_wigner_D`, `_2j`), full real-small-d matrix (`_d_matrix`), full
complex-D matrix (`_D_matrix`), block-diagonal complex-D on a multiset
(`_D_multiset`), and the `∂d/∂β` derivative for forces.

### `tensor_product.h` — the main ML kernel

Two modes: UUU (matched multiplicities, one scalar weight per path,
`irrep_tp_build` / `_apply` / `_apply_weighted`) and UVW (independent
multiplicities, full `[W, V, U]` weight tensor per path,
`_build_uvw` / `_apply_uvw`). Both have backward routines. Batched
variants. A per-path L2 regulariser (`_weight_l2_per_path_uvw`) for SR-
style training.

### `radial.h` — radial basis and cutoffs

Bessel RBF `irrep_rbf_bessel(n, r, r_cut)` with batched and derivative
forms. Gaussian RBF with a regular-grid placement helper. Cosine and
polynomial cutoffs with derivatives. Every `_d` routine is tested
against central differences.

### `quadrature.h` — numerical integration

Lebedev rules up to order 131 (131st-algebraic-order accurate on the
sphere; Lebedev & Laikov 1999); tabulated order list starts at 3 and
runs through 131. Gauss-Legendre on `[−1, 1]` for generic 1-D
integrals. A tensor-product `(Gauss-Legendre, uniform-φ)` fill on the
sphere for cases where a specific Lebedev order isn't tabulated.

### `equivariant_layers.h` — NN building blocks

`irrep_linear_t` does per-irrep-block channel mixing, preserving
`(l, p)`. `irrep_norm_rms` is per-block RMS normalisation with
learnable scales. `irrep_gate_apply` multiplies each irrep block by a
caller-supplied scalar gate (typically sigmoid or tanh of a separate
scalar feature). All three commute with `_wigner_D_multiset`.

### `nequip.h` — first-class NequIP message-passing layer

`_build` takes `(hidden_in, l_sh_max, n_radial, r_cut, cutoff_kind,
cutoff_poly_p, hidden_out)` and returns an opaque descriptor.
`_from_spec` is the spec-string shortcut. `_apply` runs the forward
pass on a small graph; `_apply_backward` accumulates gradients through
hidden features and weights; `_apply_forces` accumulates gradients
through the edge geometry (for force evaluation in equivariant MD).

### `point_group.h` — discrete-subgroup projection

`_table_build` assembles the character table, element list, and the
cached real-basis Wigner-D matrices for each `(element, l)` with
`l ≤ IRREP_PG_CACHED_L_MAX`. `_project` applies
`P_μ = (d_μ / |G|) Σ_g χ_μ*(g) D(g)`; `_reduce` runs
`m_μ = (1 / |G|) Σ_g χ_μ*(g) χ_V(g)` via the closed-form
`χ_l(θ) = sin((2l+1) θ / 2) / sin(θ / 2)` rotation character. Four
groups supported in 1.2: C₄ᵥ (5 irreps), D₆ (6 irreps), C₃ᵥ (3 irreps),
D₃ (3 irreps).

---

## Where to go next

- [`DESIGN.md`](DESIGN.md) — the architectural document, including
 layering, threading, ABI policy, and testing strategy.
- [`PHYSICS_APPENDIX.md`](PHYSICS_APPENDIX.md) — mathematical
 conventions, with derivations and a cross-reference to
 `REFERENCES.md`.
- [`tutorials/`](tutorials/) — step-by-step walkthroughs, one per
 major topic.
- [`MIGRATION_FROM_E3NN.md`](MIGRATION_FROM_E3NN.md) — side-by-side
 mapping from e3nn (Python) to libirrep (C), with sign conventions
 and gotchas.
