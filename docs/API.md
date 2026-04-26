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
| `<irrep/lattice3d.h>` (1.3) | 3D lattices with PBC (SC, BCC, FCC, Diamond, Pyrochlore) | `irrep_lattice3d_t`, `irrep_lattice3d_kind_t` | `_build`, `_free`, `_num_sites`, `_num_cells`, `_sites_per_cell`, `_Lx` / `_Ly` / `_Lz`, `_kind`, `_nn_distance`, `_nnn_distance`, `_primitive_vectors`, `_reciprocal_vectors`, `_site_position`, `_sublattice_of`, `_cell_of`, `_site_index`, `_translate`, `_num_bonds_nn`, `_num_bonds_nnn`, `_fill_bonds_nn`, `_fill_bonds_nnn`, `_k_grid` |
| `<irrep/space_group.h>` (1.3) | 2D wallpaper-group tables | `irrep_space_group_t`, `irrep_wallpaper_t` | `_build`, `_free`, `_order`, `_point_order`, `_num_sites`, `_kind`, `_lattice`, `_apply`, `_permutation`, `_permutation_inverse`, `_apply_config` |
| `<irrep/config_project.h>` (1.3) | configuration-space projection | `irrep_sg_irrep_t` | `irrep_sg_irrep_new`, `_free`, `irrep_sg_trivial`, `irrep_sg_sign_rep`, `irrep_sg_project_amplitude`, `irrep_sg_project_A1`, `irrep_sg_enumerate_orbit`, `irrep_sg_adapted_basis`, `irrep_sg_bloch_amplitude`, `irrep_sg_bloch_basis` |
| `<irrep/rdm.h>` (1.3) | reduced density matrix / entropies | — | `irrep_partial_trace`, `irrep_hermitian_eigvals`, `irrep_entropy_vonneumann_spectrum`, `_renyi_spectrum`, `_vonneumann`, `_renyi`, `irrep_topological_entanglement_entropy` |
| `<irrep/hamiltonian.h>` (1.3) | spin-½ Hamiltonian apply operators | `irrep_heisenberg_t` | `irrep_heisenberg_new`, `irrep_heisenberg_j1j2_new`, `irrep_xy_new`, `_free`, `_apply`, `_num_sites`, `_dim` |
| `<irrep/sym_group.h>` (1.3) | S_N, Young tableaux, (anti)sym projectors | — | `irrep_factorial`, `_permutation_sign`, `_permutations_all`, `_young_dim`, `irrep_sym_group_antisymmetrize`, `_symmetrize` |
| `<irrep/spin_project.h>` (1.3) | total-J projection on spin-½ chains | — | `irrep_spin_half_apply_rotation`, `irrep_spin_project_spin_half` |
| `<irrep/tensor_product.h>` (half-int path, 1.3) | spinor tensor products | `tp_2j_descriptor_t` | `irrep_tp_2j_enumerate_paths`, `_build`, `_free`, `_apply`, `_apply_weighted`, `_apply_backward`, `_output_dim`, `_num_paths` |
| `<irrep/dmi.h>` (1.3) | Bond + triangle exchange-tensor symmetry analyzers (DMI + symmetric exchange + scalar chirality + magnetic-point-group antiunitary) | `irrep_dmi_sym_op_t` | `irrep_dmi_allowed_basis`, `_from_pg`, `irrep_exchange_symmetric_basis`, `_from_pg`, `irrep_chirality_allowed`, `_from_pg`, `irrep_pg_element` |
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

Lebedev rules at orders 3, 5, 7 (octahedral-orbit form, polynomial
exactness `l = 1, 2, 3`; Lebedev & Laikov 1999). Higher-order Lebedev
tables (9..41) are on the roadmap; for arbitrary exactness, use
the tensor-product `irrep_quadrature_sphere_fill` (Gauss-Legendre in
`cos θ` × uniform φ) at ≈ 2× the point count. Gauss-Legendre on
`[−1, 1]` is also exposed for generic 1-D integrals.

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
`χ_l(θ) = sin((2l+1) θ / 2) / sin(θ / 2)` rotation character.

Seven groups supported as of 1.3-alpha:

- **2D lattice point groups** (4): C₄ᵥ (5 irreps, square), D₆ (6
  irreps, hexagonal / kagome), C₃ᵥ (3 irreps, triangular with σ_v),
  D₃ (3 irreps, triangular proper-only).
- **3D cubic point groups** (3): T_d (5 irreps A₁/A₂/E/T₁/T₂, 24
  elements; diamond / zincblende site sym), O_h (10 irreps with g/u
  parity split, 48 elements; SC / BCC / FCC / perovskite site sym),
  O (5 irreps proper-only, 24 elements; chiral cubic, MnSi /
  Cu₂OSeO₃ skyrmion-host magnets).

### `lattice3d.h` — 3D Bravais lattices

Sibling of `lattice.h`. Five families parameterised by
`irrep_lattice3d_kind_t`:

- `IRREP_LATTICE3D_SC`        — 1 site/cell, NN=1
- `IRREP_LATTICE3D_BCC`       — 2 sites/cell, NN=√3/2
- `IRREP_LATTICE3D_FCC`       — 4 sites/cell, NN=√2/2
- `IRREP_LATTICE3D_DIAMOND`   — 8 sites/cell, NN=√3/4
- `IRREP_LATTICE3D_PYROCHLORE` — 16 sites/cell, NN=√2/4

All five families use the conventional cubic cell with `a = 1`. NN
and NNN bond directions are auto-discovered per family by scanning
sublattice pairs over cell offsets in `[-1, +1]³` and retaining
those at the family's NN (or NNN) distance with a single canonical
orientation per undirected bond. Bond lists feed directly into
`irrep_heisenberg_new` for 3D ED — see
`examples/lattice3d_heisenberg_ed.c`,
`examples/lattice3d_sector_ed.c`,
`examples/lattice3d_kspace_ed.c`,
`examples/pyrochlore16_heisenberg.c`.

The `Lx ≥ 1` constraint (relaxed from L≥2) lets pyrochlore 1×1×1 =
16 sites build cleanly. Self-bonds from PBC wrap dedup automatically;
the non-trivial sublattice basis preserves NN coordination at L=1
for all families except SC (which becomes a single-site graph with
0 bonds at 1×1×1).

### `lattice.h` — 2D Bravais lattices

Sibling to `lattice3d.h` for 2D systems. Four lattice families
parameterised by `irrep_lattice_kind_t`: `IRREP_LATTICE_SQUARE`
(1 site/cell), `IRREP_LATTICE_TRIANGULAR` (1), `IRREP_LATTICE_HONEYCOMB`
(2: A, B), `IRREP_LATTICE_KAGOME` (3: A, B, C). Bond lists are
canonicalised (`i < j`) and deduplicated; primitive vectors,
reciprocal vectors, and BZ k-grids are pre-computed at build time.

Bond directions are tabulated per family rather than auto-discovered —
they're hand-tuned for the four 2D conventions to match Ashcroft–Mermin
and Elser-1989-for-kagome. The 2D bond length is **always 1** on every
NN pair so Heisenberg `J` and Hubbard `t` transfer between lattices
without rescaling. (The 3D analog uses the conventional cubic cell with
`a = 1`, so NN bond lengths differ between SC, BCC, FCC, etc. — see the
`lattice3d.h` notes.)

The kagome cluster at `Lx = 6, Ly = 6` lands at the N = 108 size that
drives the kagome NQS substrate (`docs/PHYSICS_RESULTS.md` §5).

### `space_group.h` — 2D wallpaper groups

Site-permutation tables for the 2D wallpaper groups
`IRREP_WALLPAPER_P1` (translations only), `_P2`, `_P4`, `_P4MM`, `_P6`,
`_P6MM`, `_P3M1`, `_P31M`, plus `_P4GM` scaffolded but rejected at
build time on every currently-shipped lattice (the glide ½(a₁+a₂)
maps a single-sublattice square site off the integer lattice; needs
a two-basis square lattice to land).

The group is realised as a permutation `π_g : site → site` per group
element `g`, cached per lattice. On a 6×6 kagome (108 sites,
p6mm = 432 elements), the full permutation table is 186 KB and a
single-element apply is one `memcpy(108 ints)`. This is the kernel
that drives orbit canonicalisation in
`irrep_sg_heisenberg_sector_apply` (see `hamiltonian.h`).

PBC compatibility is checked at build time: p6mm requires `Lx = Ly`
on triangular / honeycomb / kagome; p4mm requires `Lx = Ly` on square.
Mismatch returns `NULL` with `irrep_last_error()` populated.

### `config_project.h` — configuration-space projection

Build symmetry-adapted bases and project amplitudes:

- `irrep_sg_project_amplitude(G, mu, psi_orbit, ...)` — character-
  weighted projector `P_μ ψ(σ) = (d_μ / |G|) Σ_g χ_μ(g)* ψ(g·σ)`.
- `irrep_sg_adapted_basis(G, mu, ...)` — orbit-representative basis
  for sector ED in the μ-irrep block.
- `irrep_sg_bloch_amplitude(G, kx, ky, ...)` — non-Γ Bloch-momentum
  projection on the translation subgroup only.
- `irrep_sg_little_group_build(G, kx, ky)` — stabiliser at k.
- `irrep_sg_project_at_k(lg, mu_k, ...)` — composite (k, μ_k) projector.

The (k, μ_k)-block-diagonal basis from `irrep_sg_adapted_basis_at_k`
is the input to per-block dense ED for symmetry-resolved spectra
on small clusters.

### `rdm.h` — reduced density matrices, entropies, Lanczos

Three families of routines:

1. **Reduced density matrices**: `irrep_partial_trace(N, n_keep, psi,
   sites_keep, n_sites_keep, rho_out)` — trace out the complement of
   `sites_keep` from `psi` ∈ `C^{2^N}`, returning ρ_A ∈ `C^{2^{n_keep}
   × 2^{n_keep}}`.

2. **Entropies**: `irrep_entropy_vonneumann(rho, dA)`,
   `irrep_entropy_renyi(rho, dA, alpha)`, plus the spectrum
   variants that take a pre-computed eigenspectrum and the
   Kitaev–Preskill γ formula
   `irrep_topological_entanglement_entropy(S_A, S_B, S_C, S_AB, S_BC,
   S_AC, S_ABC)`.

3. **Sparse eigensolver**: three Lanczos variants —
   `irrep_lanczos_eigvals` (3-vector recurrence, ghost-prone past 100
   iterations on near-degenerate spectra), `_eigvals_reorth` (full
   Gram–Schmidt reorth, ghost-suppressed), `_eigvecs_reorth` (also
   returns Ritz eigenvectors). All take a callback
   `void apply(const double _Complex *x, double _Complex *y, void *ctx)`
   so the same kernel works on Heisenberg apply, sector apply, or any
   user-defined operator.

Hermitian eigendecomposition uses cyclic-Jacobi (`irrep_hermitian_eigvals`)
with no LAPACK dependency — sufficient for the 2D / 3D ED-validation
problem sizes (matrices up to 4096 × 4096 are routine; 65 536 × 65 536
takes minutes but works).

### `hamiltonian.h` — spin-½ apply operators

On-the-fly Heisenberg (and J₁-J₂, XY) apply operators for ED. The
canonical workflow:

```c
irrep_lattice_t   *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 4, 4);
int n = irrep_lattice_num_sites(L), nb = irrep_lattice_num_bonds_nn(L);
int *bi = malloc(nb*sizeof(int)), *bj = malloc(nb*sizeof(int));
irrep_lattice_fill_bonds_nn(L, bi, bj);
irrep_heisenberg_t *H = irrep_heisenberg_new(n, nb, bi, bj, /*J=*/1.0);
double E0;
irrep_lanczos_eigvals_reorth(irrep_heisenberg_apply, H,
                             irrep_heisenberg_dim(H), 1, 200, NULL, &E0);
```

For larger N, the **sector apply** path uses orbit-canonicalised
representatives: `irrep_sg_heisenberg_sector_build(H, T)` with
`T = irrep_sg_rep_table_build(G, ...)` precomputes
`(rep, bond) → target_rep + coefficient` transitions and then
`irrep_sg_heisenberg_sector_apply` is pure memory-bound arithmetic.
On kagome 3×3 (N = 27, 186616 reps, 27 bonds) this is ~60× faster
than the uncached path.

The **off-diagonal coefficient** for the sector apply is
`y_v += c_pm · √(N_u/N_v) · x_u` with **source-orbit-size in numerator**
— the docstring used to claim the inverse, which was a documentation
bug discovered while implementing `examples/lattice3d_sector_ed.c`
(see `CHANGELOG.md`'s "Fixed" entries).

### `dmi.h` — bond + triangle exchange-tensor symmetry analyzers

Three companion analyzers covering the bilinear and trilinear spin-
spin couplings:

1. **Antisymmetric bond exchange (DMI vector)** —
   `irrep_dmi_allowed_basis(_from_pg)`. Takes a bond `{r_a, r_b}`
   and returns the orthonormal basis of allowed `D` vectors
   (0–3 dim). Implements all five Moriya rules (1960) by direct
   projector construction.

2. **Symmetric bond exchange** —
   `irrep_exchange_symmetric_basis(_from_pg)`. Takes a bond and
   returns the orthonormal basis of allowed 3×3 symmetric matrices
   (0–6 dim). Heisenberg + Kitaev-Γ-type anisotropy.

3. **Scalar three-spin chirality** —
   `irrep_chirality_allowed(_from_pg)`. Takes a triangle
   `{r_a, r_b, r_c}` and decides whether `χ_ijk = S_i · (S_j × S_k)`
   is symmetry-allowed (returns 0 / 1 — no basis since χ is a scalar).
   The pseudoscalar selection rule is `det(g) · σ_perm(g) · σ_T(g)
   = +1` for every preserving operation (where `σ_perm = ±1` is
   the permutation parity, `σ_T = -1` for antiunitary).

The shared input type `irrep_dmi_sym_op_t = {R_proper[9], det,
antiunitary}` carries the magnetic-point-group augmentation: set
`antiunitary = 1` for `T·g` operations to handle Shubnikov-type
magnetic groups. The DMI and chirality analyzers consume this flag;
the symmetric-exchange analyzer is invariant to it (rank-2 axial-
axial → polar tensors are invariant under `T`).

Verified against textbook results:
- **B20 chiral magnet** pattern `D ∥ bond` from chiral cubic O —
  reproduced from group theory alone (`D · bond̂ = +1.000000`).
- **Pyrochlore NN** under O_h: `D = 0`, J^s 3-dim Curnoe-Ross-Kao
  parametrisation.
- **Kagome NN** under D_6: `D ∥ bond` (chiral hexagonal).
- **Kagome triangle** chirality under D_3h: forbidden; under T·σ_h
  ("magnetic mirror" of Mn₃Sn): allowed → topological Hall.

Worked examples: `examples/dmi_pyrochlore_pattern.c`,
`examples/dmi_kagome_pattern.c`,
`examples/kagome_triangle_chirality.c`. Mathematical derivations
in `docs/PHYSICS_APPENDIX.md` §13 (bond bilinear) and §14 (magnetic
point groups + antiunitary).

This module is the foundation of libirrep's **materials-search
pipeline**: given a candidate space group + magnetic ordering, the
analyzers return the complete symmetry-allowed exchange-tensor
structure (DMI, J^s, χ) that downstream codes (DFT, mumax, OOMMF)
need as parameter scaffolding. It automates the crystallographer's
hand-derivation from International Tables vol. A.

### `sym_group.h` — symmetric group, Young tableaux

Permutations and (anti)symmetrisation utilities for fermion / boson
Fock-space construction. `irrep_factorial`, `_permutation_sign`,
`_permutations_all` enumerate `S_N` for small N. `_young_dim` gives
the dimension of the Young-tableaux irrep. `irrep_sym_group_antisymmetrize`
and `_symmetrize` apply the corresponding projector to a vector in
`C^{N!}`-equivalent indexing.

Used by Hubbard / fermion-Heisenberg consumers; the spin-½ Heisenberg
in `hamiltonian.h` doesn't need it (Sz_total conservation handles
the sector decomposition cheaper).

### `spin_project.h` — total-J projection

`irrep_spin_project_spin_half(two_J_target, N, l_max_quad, ...)`
projects an `N`-site spin-½ amplitude vector onto the total-J
sector with `J = two_J_target/2`. Implementation: integrate
`χ_J(R)·U(R)|ψ⟩` over SO(3) using a Lebedev × Gauss-Legendre
quadrature on Euler angles, with `χ_J(R)` the SO(3) character.
The Marshall sign-rule structure of `1/2`-singlets falls out of
`two_J_target = 0`.

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
