# Design

An architecture-level overview of libirrep. For the mathematical
conventions and their primary sources, see
[`PHYSICS_APPENDIX.md`](PHYSICS_APPENDIX.md) and
[`REFERENCES.md`](REFERENCES.md). For per-module API, the individual
headers under `include/irrep/` carry Doxygen blocks that `make docs`
renders into HTML. For concrete physics results produced with the
library, see [`PHYSICS_RESULTS.md`](PHYSICS_RESULTS.md).

---

## 0. Feature → module mapping

Every 1.3-cycle module exists because a specific open physics problem
needs a specific primitive. The table below is the load-bearing
connection between the 1.3 scope (headline: Kagome Heisenberg
ground-state nature) and the library's surface area:

| research need | module | concrete action |
|------------------------------------------------------------------|-------------------------------------|--------------------------------------------------------------------------------------|
| Build lattices, enumerate bonds, access the Brillouin zone | `lattice.h` | `irrep_lattice_build` + NN/NNN bond lists + `_k_grid` |
| Project NQS ansätze onto a specific p6mm irrep | `space_group.h` + `config_project.h` | site-permutation tables + `irrep_sg_project_amplitude` / `_project_A1` |
| Decompose an ED ground state across symmetry sectors | `config_project.h` | `irrep_sg_adapted_basis` + standard Jacobi on each block |
| Diagnose gapped-Z₂-vs-trivial via topological entanglement | `rdm.h` | `irrep_partial_trace` + `_vonneumann` + `irrep_topological_entanglement_entropy` |
| Restrict to singlet sector (Marshall / Lieb-Mattis) | `spin_project.h` | `irrep_spin_project_spin_half(two_J_target = 0)` |
| ED at clusters too big for dense storage | `rdm.h` | `irrep_lanczos_eigvals` with callback H-apply |
| Fermion antisymmetrisation for 2D Hubbard target | `sym_group.h` | `irrep_sym_group_antisymmetrize` |
| Spin-orbit-coupled / spinor equivariant NN layers | `tensor_product.h` half-integer | `irrep_tp_2j_build` + `_apply_weighted` + `_apply_backward` |

The 1.2 core (spherical harmonics, Clebsch-Gordan, Wigner-D, NequIP
layer, point-group projection on feature vectors) sits below this
physics substrate and is reused intact; see §2 for its layering.

---

## 1. Layering

The library is stratified by dependency direction. Every module depends
only on modules strictly below it in the stack, never laterally or
upward. This lets a consumer include only the headers they need without
pulling in transitive implementation details.

```
 ┌────────────────────────────────┐
 │ types.h version.h export.h │ atomic primitives
 │ simd.h │ (no math)
 └───────────────┬────────────────┘
 │
 ┌──────────────────────────────┼──────────────────────────────┐
 ▼ ▼ ▼
 ┌──────────┐ ┌────────────────┐ ┌──────────────┐
 │ so3.h │ │ multiset.h │ │ parity.h │
 │ su2.h │ │ multiset_2j.h │ │ time_reversal│
 │ │ │ (irrep algebra)│ └──────────────┘
 └────┬─────┘ └────────┬───────┘
 │ │
 │ ▼
 │ ┌──────────────────────┐
 ▼ │ clebsch_gordan.h │
┌─────────────────┐ │ recoupling.h │
│ spherical_ │ │ (CG, 3j, 6j, 9j) │
│ harmonics.h │ └──────────┬───────────┘
│ solid_harmonics│ │
│ radial.h │ ▼
│ quadrature.h │ ┌───────────────────────┐
└────────┬────────┘ │ wigner_d.h │
 │ │ tensor_product.h │
 └─────────┬───────┴────────┬──────────────┘
 │ │
 ▼ ▼
 ┌────────────────┐ ┌──────────────────────┐
 │ equivariant_ │ │ nequip.h │
 │ layers.h │ │ point_group.h │
 │ (linear / norm │ │ (composite layers, │
 │ / gate) │ │ projectors) │
 └────────────────┘ └──────────────────────┘
 │
 ▼
 ┌───────────────────────────────────────────────┐
 │ 1.3 physics substrate () │
 │ lattice.h (2D Bravais, PBC) │
 │ space_group.h (p1, p4mm, p6mm tables) │
 │ config_project.h (P_μ ψ(σ) reducer) │
 │ rdm.h (partial trace, S_VN, γ) │
 │ sym_group.h (S_N, hook-length, A/S) │
 │ spin_project.h (total-J projection) │
 │ tensor_product.h [_2j] (half-integer UVW TP) │
 └───────────────────────┬───────────────────────┘
 │
 ▼
 ┌────────────┐
 │ irrep.h │ umbrella (pulls all)
 └────────────┘
```

The layering has three practical consequences:

- Consumers who want only `irrep_quat_slerp` include `<irrep/so3.h>`
 alone and pay for nothing else at link time.
- The test suite can exercise each module in isolation (one
 `tests/test_<module>.c` per module), so a failure in tensor products
 cannot mask a bug in Clebsch-Gordan.
- New modules land at an explicit layer — `point_group.h` in 1.2 slotted
 above `wigner_d.h` because it consumes the real-basis D matrix, and
 below the umbrella because it has no dependents in core.

---

## 2. Module responsibilities

### 2.1. Core primitives (`types.h`, `version.h`, `simd.h`, `export.h`)

- **`types.h`**: `irrep_quaternion_t`, `irrep_rot_matrix_t`,
 `irrep_euler_zyz_t`, `irrep_axis_angle_t`, `irrep_label_t`,
 `irrep_multiset_t`, status enum, and library-wide constants
 (`IRREP_L_MAX = 16`, `IRREP_TWO_J_MAX = 32`).
- **`version.h`**: compile-time `IRREP_VERSION_*` macros, runtime
 accessors, and the `irrep_abi_hash()` accessor whose value is baked in
 at release time via `-DIRREP_BAKED_ABI_HASH`.
- **`simd.h`**: CPU feature detection surface (`irrep_cpu_has_neon()`,
 `…has_avx2()`, etc.). Populated lazily on first call via a C11
 `atomic_flag`-guarded double-checked lock in `src/simd_runtime.c`.
- **`export.h`**: the `IRREP_API` visibility macro
 (`__attribute__((visibility("default")))` on gcc/clang, `__declspec`
 on Windows).

### 2.2. Rotation math (`so3.h`, `su2.h`)

Every pair of conversions between quaternion / rotation matrix / Euler
ZYZ / axis-angle, round-tripping to `10⁻¹²`; composition, inverse,
Rodrigues exponential with a small-angle Taylor guard, Markley
branch-switching logarithm stable near π; Shoemake uniform sampling;
SLERP; Karcher-Fréchet mean. SU(2) provides the Pauli matrices, SU(2) ↔
quaternion conversion, and the explicit double-cover map
`π: SU(2) → SO(3)`.

### 2.3. Irrep algebra (`multiset.h`, `multiset_2j.h`, `parity.h`, `time_reversal.h`)

The integer-l multiset type `irrep_multiset_t` encodes
`A = ⊕ mult_i · (l_i, p_i)` as a `(labels, multiplicities)` pair with
builders: `irrep_multiset_parse` accepts the e3nn-style string grammar
`"1x0e + 2x1o + 1x2e"`, `_append` grows the backing storage,
`_simplify` merges like terms and sorts canonically, `_direct_sum`
concatenates two multisets. The parallel `irrep_multiset_2j_t` type
(`multiset_2j.h`) carries doubled-integer labels so half-integer spin
content (`"1x1/2o"`, `"2x3/2e"`) is expressible; it has its own
parser and a Kramers-sign query
`irrep_time_reversal_square_sign_2j` that returns `−1` exactly when
any block has odd `two_j`. Parity helpers enforce `p_c = p_a · p_b`
on tensor-product paths. Time-reversal exposes the antiunitary T
operator on integer-l, half-integer-j, and the multiset-level
composition.

### 2.4. Coupling coefficients (`clebsch_gordan.h`, `recoupling.h`)

Wigner 3j via Miller two-directional iteration over the Schulten–Gordon
three-term recurrence in `j₁` (Luscombe–Luban 1998): forward pass from
`j_min` and backward pass from `j_max` are spliced at the index that
maximises `|T_fwd|·|T_bwd|` (deep in each iteration's classical region);
the result is normalised by `Σ_j (2j+1)·T(j)² = 1` and sign-anchored to
`(−1)^{j₂−j₃−m₁}` at `j_max`. Clebsch–Gordan derived from 3j via
`⟨j₁ m₁ j₂ m₂ | J M⟩ = (−1)^{j₁−j₂+M} · √(2J+1) · (j₁ j₂ J; m₁ m₂ −M)`.
Machine precision through at least `j ≈ 200` in regression tests; the
prior backward-only pass degraded past `j ≈ 80` from subdominant-solution
contamination. Selection rules return `0.0`; a small `cg_table_t` cache
supports batch lookups. Wigner 3j, 6j, 9j, and Racah W are in
`recoupling.h` and use the same log-gamma kernel.

### 2.5. Spherical harmonics and geometry (`spherical_harmonics.h`,
 `solid_harmonics.h`, `radial.h`, `quadrature.h`)

Three API surfaces for SH: polar complex `Y_l^m(θ, φ)`, polar real
`Y_{l, m}^{real}`, and cartesian `Y^{real}(r̂)`. The cartesian form is
the hot path for graph neural networks — every NequIP / MACE message
layer calls it per edge. Associated Legendre via the stable three-term
recurrence; complex↔real basis change. Solid harmonics
`R_{l, m}(r) = |r|^l · Y^{real}_{l, m}(r̂)` for analytic-gradient force
evaluation. Radial primitives (Bessel, Gaussian) with derivatives and
cutoffs (cosine, polynomial). Lebedev and tensor-product-of-
Gauss-Legendre quadrature on the sphere.

### 2.6. Rotations on irreps (`wigner_d.h`)

Scalar and matrix `D^j_{m', m}(α, β, γ)`, `d^j(β)`, half-integer forms,
the `∂d/∂β` derivative, and the block-diagonal Wigner-D matrix on a
whole multiset. Real-basis Wigner-D is built by sandwiching the complex
form between the complex-to-real basis change (`U · D · U†`) — used
internally by both `tensor_product.c` (for the `i^{l_a + l_b − l_c}`
phase-corrected real-basis CG tensor) and `point_group.c` (for the
projection operator).

### 2.7. Tensor products (`tensor_product.h`)

The heart of the library for ML consumers. Path-indexed e3nn-style
tensor products, UUU mode (matched multiplicities, one weight per path)
and UVW mode (independent multiplicities, full `[W, V, U]` weight
tensor). Forward, backward, batched. The parity filter enforces
`p_a · p_b = p_c`; the `i^{l_a + l_b − l_c}` phase factor is applied
internally so that odd-l-sum paths (cross-product and similar) produce
guaranteed-real output.

### 2.8. Composite layers (`equivariant_layers.h`, `nequip.h`,
 `point_group.h`)

- **`equivariant_layers.h`** — the minimum viable set of primitives
 needed to build an equivariant MLP from libirrep alone:
 linear-on-irreps (channel mixing within each block), RMS norm per
 block, and multiplicative gate. All commute with the block-diagonal
 Wigner-D by construction.
- **`nequip.h`** — first-class `irrep_nequip_layer_t` composing edge SH +
 radial basis + cutoff + UVW tensor product + aggregation. Forward
 (`_apply`), backward through weights and hidden state
 (`_apply_backward`), and backward through edge geometry
 (`_apply_forces`, wired for the Landau–Lifshitz–Gilbert (LLG)
 torque-network / force-evaluation
 consumer). Spec-string constructor (`_from_spec`).
- **`point_group.h`** — character tables, projection operator, and
 direct-sum decomposition under the supported point groups (C₄ᵥ, D₆,
 C₃ᵥ, D₃). D^l matrices for the elements are pre-computed at
 table-build time for small l, giving an ~80× speedup on the
 projection hot path.

### 2.9. 1.3 physics substrate (`lattice.h`, `space_group.h`, `config_project.h`, `rdm.h`, `sym_group.h`, `spin_project.h`, `tensor_product.h` half-integer path)

The 1.3 cycle adds a substrate for neural-quantum-state work.
Its scope is pinned in
the 1.3 CHANGELOG
and each module here is sized against the target: diagnosing the Kagome
Heisenberg S = ½ ground-state nature on a 6×6 × 3 = 108-site cluster
with genuine p6mm space-group projection.

- **`lattice.h`** — 2D Bravais lattices (square, triangular, honeycomb,
 kagome) with periodic boundary conditions. Site coordinates, sublattice
 lookup, NN / NNN bond enumeration (canonicalised and deduplicated under
 PBC wrap), primitive / reciprocal vectors, cell translations, and
 Brillouin-zone k-grids. Primitive-vector conventions are chosen so every
 NN bond has length 1 on every lattice — Heisenberg J and Hubbard t
 transfer verbatim.
- **`space_group.h`** — p1, p4mm, p6mm wallpaper-group actions. At build
 time each element is realised as a site permutation via cartesian
 matrix × quantised-inverse-basis lookup; the builder rejects clusters
 whose torus breaks the symmetry. Queries are O(1): `apply` is an array
 read, `permutation` is a `memcpy`. The 6×6 kagome cluster costs 186 KB
 per direction (forward + inverse), and the full 432-element orbit sum
 benchmarks at 17.6 µs end-to-end — subordinate to any realistic
 wavefunction evaluation.
- **`config_project.h`** — character-weighted reducer for
 `P_μ ψ(σ) = (d_μ / |G|) Σ_g χ_μ*(g) ψ(g·σ)`. The library holds the
 group-theoretic bookkeeping (characters, orbit pullbacks); the caller
 supplies the amplitudes. Distinct from `point_group.h`, which projects
 *internal irrep-space feature vectors* rather than *configurations*.
- **`rdm.h`** — the entanglement primitives the kagome diagnosis needs.
 Partial trace for `N ≤ 30` sites of arbitrary local dimension;
 phase-reduction + real-Givens cyclic-Jacobi Hermitian
 eigendecomposition (stable, no LAPACK dependency); von Neumann and
 Rényi entropies; the Kitaev–Preskill topological-entanglement-entropy
 formula `γ = S_A + S_B + S_C − S_{AB} − S_{BC} − S_{AC} + S_{ABC}`
 (γ = ln 2 signals Z₂ topological order; γ = 0 signals no topological
 order).
- **`sym_group.h`** — `S_N` infrastructure for the fermion programme:
 factorial, permutation sign, lexicographic full-permutation enumeration
 (`N ≤ 10`), hook-length irrep dimension formula, and antisymmetric
 (fermion) and totally-symmetric (boson) projectors on tensor-factored
 states.
- **`spin_project.h`** — total-J projection on `N` spin-½ sites via the
 Wigner-D character-weighted SU(2) integral. The N-fold tensor-product
 rotation factors as sequential single-qubit updates (O(N · 2^N) per
 Euler grid point). Used to restrict NQS ansätze to the singlet
 (`J = 0`) sector at half filling — where the even-site antiferromagnet
 ground state lives.
- **`tensor_product.h` (half-integer path)** — The half-integer path of `tensor_product.h` extends
 `tensor_product.h` to `irrep_multiset_2j_t`, on complex amplitudes
 (half-integer irreps carry no natural real basis). Enumeration, build,
 forward, weighted forward, and backward all available; spin-½ ⊗ spin-½
 = spin-0 ⊕ spin-1 matches Clebsch–Gordan to machine precision.

Each 1.3 module has a downstream deliverable listed in
`CHANGELOG.md §Coordination`.

### 2.10. SIMD and dispatch (`src/internal/dispatch.h`, `simd_runtime.c`)

Not exposed publicly. An internal function-pointer table (the
`irrep_dispatch_t` struct in `src/internal/dispatch.h`) holds per-kernel
pointers; at first use, `irrep_dispatch_get()` populates the table from
the widest available SIMD variant. Currently populated:
- `cutoff_polynomial_batch` → AVX2+FMA on x86_64 / NEON on aarch64 / scalar
- `cutoff_polynomial_d_batch` → AVX2+FMA / NEON / scalar
- `sph_harm_cart_all_batch` → AVX2+FMA / NEON / scalar

Dispatch overhead is one indirect call per batched entry point;
measurably unchanged against direct-call baseline in the benchmark
suite.

---

## 3. Numerical conventions (authoritative)

See [`PHYSICS_APPENDIX.md`](PHYSICS_APPENDIX.md) for derivations. The
one-page summary:

- **Precision**: `double` default; `_f32` wrappers cast from double.
- **Euler ordering**: ZYZ physics convention.
- **Phase**: Condon–Shortley, applied once.
- **Rotations**: active, right-handed.
- **Quaternion layout**: `{x, y, z, w}` with scalar `w` last.
- **Half-integer**: `_2j` suffix takes doubled-integer arguments.
- **FP contraction**: off in the SH scalar and NEON paths (so the two
 stay bit-identical); default elsewhere.

Error-handling convention:

- **Pure math functions** that cannot meaningfully fail return the
 mathematically-correct value on valid input and `0.0` (or the natural
 zero) on invalid input. A caller summing over a large index range
 doesn't need to filter selection-rule violations first — they sum to
 the right answer because forbidden terms are exact zero.
- **Builder functions** (`irrep_tp_build`, `irrep_multiset_parse`,
 `irrep_nequip_layer_build`, `irrep_pg_table_build`) return `NULL` on
 failure and set the thread-local `irrep_last_error()` buffer.
- **Validating functions** that could report structured errors return an
 `irrep_status_t` enum.
- No `errno`, no `setjmp`/`longjmp`, no asserts that fire in release
 builds.

---

## 4. Threading model

- **Pure math functions** are reentrant and thread-safe by construction
 (no shared state).
- **Built tables** (`cg_table_t`, `tp_descriptor_t`, `irrep_multiset_t`,
 `irrep_linear_t`, `irrep_nequip_layer_t`, `irrep_pg_table_t`) are
 safe to **read** concurrently from multiple threads once fully built,
 and unsafe to build / free concurrently with any read on the same
 handle.
- **`irrep_last_error()`** is thread-local (`_Thread_local` with
 `__thread` fallback).
- **SIMD feature table and dispatch pointers** are populated once via an
 atomic_flag-guarded critical section in `src/simd_runtime.c`; readers
 see a fully-initialised struct after the first caller completes the
 populate step. The implementation is formally race-free under the C11
 memory model (see the commentary in `src/simd_runtime.c` — an earlier
 revision with a plain `atomic_bool` flag was upgraded after the
 1.2.0-rc1 audit identified the race).

---

## 5. ABI policy

The ABI surface is defined as:

1. The sorted set of exported symbol names (those tagged `IRREP_API` in
 the public headers).
2. The memory layout of every publicly-visible struct (field names,
 types, offsets, and total size).

`scripts/generate_abi_hash.sh` produces SHA-256 over the canonicalised
form of that set. At release time, `scripts/build_release.sh` performs a
two-pass build: pass 1 computes the hash from the placeholder binary,
pass 2 recompiles with `-DIRREP_BAKED_ABI_HASH="<hash>"` so that
`irrep_abi_hash()` at runtime returns the same string that
`scripts/generate_abi_hash.sh` produces, and that the installed
pkg-config file exposes as the `abi_hash` variable.

Versioning follows SemVer 2.0:

- **MAJOR**: ABI hash changes in a way that breaks consumers — renamed
 or removed symbols, changed struct layout, changed calling
 convention.
- **MINOR**: additive API surface. ABI hash changes but existing
 consumers are unaffected.
- **PATCH**: internal fixes. ABI hash should not change.

The CI workflow (`.github/workflows/ci.yml`) has an `abi-check` job that
diffs the built hash against `abi/baseline/ABI_HASH.<triple>` and fails
the build if the hash changed without a corresponding major bump.

---

## 6. Release engineering

Release artefacts are produced by `scripts/build_release.sh` and laid
out under `release/<version>/`:

```
release/<version>/
├── include/irrep/ # all 21 public headers
├── lib/<triple>/ # static .a + shared .so|.dylib|.dll
├── pkgconfig/libirrep.pc # with baked abi_hash variable + embedded rpath
├── VERSION # plain text version string
├── ABI_HASH # sha256 of the canonicalised export set
├── CHECKSUMS # sha256 over every file in the release tree
└── SBOM.spdx.json # SPDX 2.3 with per-file hash and license
```

Per-triple tarballs `liblibirrep-<version>-<triple>.tar.gz` are produced
alongside. Cross-compilation to multiple triples from a single host is a
roadmap item (see
the 1.3 CHANGELOG);
the 1.2 release produces the current host's triple only.

---

## 7. Testing strategy

- **Per-module pass / fail tests** in `tests/test_<module>.c`, one per
 public header. TAP-style output via `tests/harness.h` — plain C with
 no external framework.
- **Reference-data cross-checks** in `tests/reference_data/`. CG values
 are stored as `(sign, num, den)` rationals so the test reconstructs
 `sign · sqrt(num / den)` via libm `sqrt` and comparisons are bit-
 aligned with the library's rounding.
- **Error-path tests** in `tests/test_error_paths.c` exercise NULL
 pointers, out-of-range indices, malformed multiset strings, invalid
 NequIP parameters, etc.
- **Downstream-compat tests** in `tests/test_downstream_compat.c` — a
 placeholder that lights up bit-exact assertions once the downstream
 vendored golden-vector fixtures land.
- **Fuzz targets** in `tests/fuzz/` (currently six: CG, multiset parse,
 tensor product, NequIP apply_forces, NequIP from_spec,
 point-group projection), run in CI for 60 s per target.
- **Sanitizer CI** — every PR rebuilds and tests under ASan and UBSan.
 As of 1.2.0 the suite passes clean under both.
- **Benchmark baseline** — `benchmarks/results/baseline/<triple>.json`
 is checked into the tree; `scripts/perf_compare.sh` diffs a fresh run
 against it with a 5%-plus-1ns threshold and prints CPU-model
 metadata on mismatch so a different CI machine doesn't false-positive.

---

## 8. Coordination with downstream consumers

libirrep is developed in tandem with its downstream consumers (primarily
`spin_based_neural_network` and `quantum_geometric_tensor`). Each
release cycle has a coordination document that enumerates:

1. Deliverables owed by libirrep (new APIs, performance targets).
2. Deliverables owed by the consumer (fixtures, parser specs, PRs).
3. Pinned compatibility (VERSION ranges both trees expect).
4. Schedule with owner per item.

The current cycle's doc is
the 1.3 CHANGELOG.
Mirror-doc on the downstream side lives in
`spin_based_neural_network/docs/`. Both trees commit acceptance
criteria for cross-tree changes before implementation, so a regression
in either codebase surfaces in the other's CI before reaching main.

---

## 9. Open items

Tracked per-release in `TODO.md` and, more formally, in the coordination
document for the current cycle. Architectural items actively deferred
past 1.3 are called out as non-goals in the README.
