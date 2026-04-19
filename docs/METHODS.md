# Methods

A methods-paper-style writeup of the primitives that libirrep
implements and the algorithmic choices behind them. Intended as the
technical backbone of a future software-description publication (JOSS
or CPC) paired with a physics-results publication (see
[`PHYSICS_RESULTS.md`](PHYSICS_RESULTS.md)).

Every equation here cites a primary source; a full bibliography is in
[`REFERENCES.md`](REFERENCES.md).

---

## 1. Statement of the problem

Equivariant machine learning on atomic-scale physics (materials,
molecules, spin systems) and numerical studies of frustrated quantum
magnetism share a common mathematical substrate: the irreducible
representations of the rotation group SO(3) (and its covers SU(2)
and O(3)). Every production ML code in the e3nn family (NequIP,
MACE, Allegro, …), every serious study of a 2D frustrated Heisenberg
model, and every t-matrix calculation in atomic physics requires a
library that can evaluate:

- Spherical harmonics `Y_l^m(r̂)` at moderate `l` (up to 16),
- Clebsch-Gordan coefficients and Wigner 3j / 6j / 9j symbols for
 both integer and half-integer spin,
- Wigner small-d and full-D rotation matrices, stable past `j = 50`,
- Tensor products `a ⊗ b` with irrep-preserving contraction,
- Point-group and space-group projection operators on feature or
 configuration vectors,
- Entanglement entropies and partial traces on small-system state
 vectors,
- Total-J and symmetry-sector restrictions on spin wavefunctions.

Existing public implementations are fragmented: the numerical machinery
lives in the private source of individual ML codebases (with
mutually incompatible sign conventions), in specialist ED codes (with
their own coupling geometries), or in general-purpose libraries
(SymPy, sympy.physics) at prohibitive cost. **libirrep** is a single
pure-C11 library, MIT-licensed, with a stable ABI, that shipped these
primitives at `10⁻¹⁰` relative accuracy under one canonical convention
catalogue (Condon-Shortley phase, ZYZ Euler, right-handed active
rotations, xyzw quaternion layout).

## 2. Conventions

Every formula in libirrep follows the following conventions, matching
Sakurai (1994), Varshalovich *et al.* (1988), and Geiger *et al.* 2022
(the e3nn reference):

- **Angles** — radians throughout.
- **Rotations** — active, right-handed. Positive rotation about `ẑ`
 sends `x̂ → ŷ`.
- **Euler** — ZYZ (Sakurai §3.3), with ranges `α ∈ [0, 2π)`,
 `β ∈ [0, π]`, `γ ∈ [0, 2π)`.
- **Phase** — Condon-Shortley `(-1)^m` applied once, in the
 associated Legendre polynomial.
- **Real spherical harmonics** — e3nn sign convention,
 `Y_{1,+1}^{real} ∝ +x` at the equator.
- **Quaternions** — `{x, y, z, w}` layout with scalar `w` last
 (SIMD-ergonomic; matches Eigen and glTF 2.0); Shoemake-sampled
 quaternions are canonicalised to `w ≥ 0`.
- **Half-integer** — `_2j` suffix, so spin-½ enters as `two_j = 1`
 and the type system distinguishes doubled-integer arguments from
 integer ones by name.

Full derivations and edge-case handling (gimbal lock, rot_log near
π, Shoemake's uniformity proof) are in
[`PHYSICS_APPENDIX.md`](PHYSICS_APPENDIX.md).

## 3. Algorithmic choices

### 3.1. Clebsch-Gordan via Racah's single-sum formula in log-gamma form

Direct factorial evaluation of Racah's formula overflows at modest
`j`. libirrep uses the log-gamma form throughout, summing signed
`exp(log Γ(…))` terms; intermediate values stay within double-
precision range past `j = 100`. Sign tracking avoids the catastrophic
cancellation that plagues naive exp-after-log chains.

Reference: Racah 1942, Phys. Rev. 62, 438. Implemented in
`src/clebsch_gordan.c`; cross-validated against Sakurai Appendix A
hand-tabulated values at
`(j₁, j₂) ∈ {(½,½), (1,½), (1,1), (3/2,1), (2,2), (3,2)}` to `10⁻¹²`.

### 3.2. Wigner small-d via Sakurai's direct sum

Jacobi-polynomial recurrences for the small-d function
`d^j_{m'm}(β)` are standard but can lose precision near `β = π`.
libirrep uses Sakurai's direct sum (Sakurai 3.8.33) summed in
log-gamma form; verified stable past `j = 50` in our unit tests,
and up to `j = 100` by cross-check against a SymPy-generated
reference.

### 3.3. Cartesian spherical harmonics via stable three-term recurrence

Associated Legendre polynomials `P_l^m(cos θ)` are computed by the
stable three-term recurrence in `z = cos θ`, paired with a
`cos(mφ), sin(mφ)` recurrence computed directly from `(x, y)` —
all in Cartesian, no transcendental calls. Near the poles
(`r_{xy} < 10⁻⁹`), a degenerate branch handles the `sin(mφ) = 0,
cos(mφ) = 1` limit smoothly. NEON batched kernel is bit-exact
against the scalar path by way of `#pragma STDC FP_CONTRACT OFF` in
both translation units.

Reference: Numerical Recipes 3e §6.7; Limpanuparb & Milthorpe 2014.

### 3.4. Tensor products with real-basis output

e3nn-style path-indexed tensor products are implemented in real-basis
arithmetic for efficiency. The complex-basis CG coefficients connect
via an `i^(l_a + l_b − l_c)` phase factor on odd-l-sum paths, applied
internally so that the `(1, 1, 1)` cross-product path — `a × b =
√2 · (1o ⊗ 1o → 1e)` in our layout — produces guaranteed-real output.
Cross-checked bit-exactly against a Cartesian reference in
`examples/torque_net_tp_paths.c`.

### 3.5. Point-group projection with pre-computed real-basis D matrices

The character-weighted projector on irrep-space feature vectors is
`P_μ = (d_μ / |G|) Σ_g χ_μ^*(g) D(g)`, where `D(g)` is the real-basis
Wigner-D matrix for the element `g`. For `l ≤ 4` (covering every
common NequIP / MACE SH degree), all `D(g)` matrices are computed
once at table-build time; each subsequent projection becomes a pure
matrix-vector multiply. This delivers a sustained ~80× speedup on
the projection hot path at moderate cache cost (~2–16 KB per group
table).

Groups supported in 1.3: C₄ᵥ, D₆, C₃ᵥ, D₃ (Bradley-Cracknell 1972,
cross-checked against Altmann-Herzig 1994).

### 3.6. Space-group site permutations with inverse-basis lookup

A space-group element `g = (t, R)` acts on a lattice site at
cartesian position `r` as `r' = R · (r − O) + O + t` where `O` is the
wallpaper-group origin. libirrep materialises every element as a
site permutation `π_g : [0, N) → [0, N)` at table-build time; each
subsequent application is a pointer read.

Computing the permutation requires mapping a cartesian image `r'`
back to a site index. We decompose `r' = Σ n_i a_i + δ_s` in
the lattice basis via inverse-basis multiplication; integer-ness of
`(n₁, n₂)` within tolerance `10⁻⁸` selects the site. This is the
stable alternative to quantised-position hashing, which breaks on
rotation-induced ULP differences between `√3/2` and `cos(π/3)`.

On the 6×6 × 3 = 108-site kagome cluster, the full p6mm space group
has 432 elements; applying one element is 1.35 ns, applying to a
full configuration is 40.6 ns. See
[`PHYSICS_RESULTS.md`](PHYSICS_RESULTS.md) §3.1.

### 3.7. Symmetry-adapted basis via orbit-representative projection

The standard character-weighted basis builder — project every
computational basis state `|s⟩` with `P_μ|s⟩`, Gram-Schmidt
orthogonalise — is O(D · |G|) in the apply phase, where `D` is the
total Hilbert-space dimension and `|G|` is the space-group order.
Because `P_μ|s⟩` and `P_μ|g · s⟩` are proportional (for 1D irreps)
or span the same invariant subspace (for higher-dim irreps), a large
fraction of seeds produces redundant basis vectors.

libirrep's `irrep_sg_adapted_basis` filters to **orbit
representatives**: for each seed `s`, check whether some `g · s < s`
(where comparison is the natural `<` on integer indices); if so, skip
`s`. Each orbit is processed exactly once. On 12-site kagome p6mm
this is a **31× speedup** (from 33 s to 1.07 s) with bit-identical
output. The correctness relies on Schur's lemma for 1D irreps and
degrades gracefully for 2D irreps (extra vectors are caught by
Gram-Schmidt).

### 3.8. Sparse Hermitian eigensolver

`rdm.h` ships two complementary eigensolvers:

- `irrep_hermitian_eigvals` — cyclic-Jacobi with phase-reduction +
 real-Givens rotations. Converges in ~10 sweeps at `10⁻¹⁴` for a
 dense Hermitian; O(n³) work per call. Used for block ED of
 symmetry-adapted 10²–10³ dim blocks, and for small partial-trace
 density-matrix spectra.
- `irrep_lanczos_eigvals` — 3-term-recurrence Lanczos with a
 callback-based `apply_op`. Extracts the lowest few eigenvalues of
 a sparse Hermitian given only matrix-free multiplication. No
 Lanczos reorthogonalisation (storage: 3 state vectors); ~50
 iterations converge the ground state at high accuracy.
 Cross-validated against Jacobi on 32×32 random Hermitian matrices
 at `10⁻¹⁰` relative accuracy (see `tests/test_rdm.c`).

For the 24-site kagome cluster (Hilbert dim 2²⁴ = 16 777 216),
Lanczos with 80 iterations converges the ground state in ~74 s on
an Apple M2 Ultra; 3 state vectors consume ~768 MB. The same
callback infrastructure supports spin gap extraction via S_z = 1
seeding, deflated first-excited-state extraction, and arbitrary
Hermitian operators other than Heisenberg.

### 3.9. Kitaev-Preskill topological entanglement entropy

Given a tripartition `A, B, C` of a quantum region and a density
matrix `ρ` on the full region, the topological entanglement entropy
is (Kitaev & Preskill 2006):

 γ = S_A + S_B + S_C − S_AB − S_BC − S_AC + S_ABC

For a gapped state with a proper annular tripartition geometry (each
pair of regions shares a boundary; `A ∪ B ∪ C` is topologically an
annulus), `γ` is a topological invariant. For Z₂ topological order
γ = ln 2 ≈ 0.693; for a trivially-ordered state γ = 0.

libirrep computes `γ` via seven calls to `irrep_partial_trace`
followed by seven calls to `irrep_hermitian_eigvals`. Validated to
machine precision on the 4-qubit GHZ state (known γ = +ln 2; see
`tests/test_rdm.c`). On the 12-site kagome Heisenberg ground state
with an annular tripartition (trace out 6 sites, split the kept 6
into 2+2+2), the pipeline returns a finite-size-dominated value of
−0.33 nats — not a physical γ (2×2 torus is too small for a proper
annulus) but a proof the pipeline runs end-to-end on real ED data.

### 3.10. Total-J projection via Wigner-D integration

For spin-½ on `N` sites, the total-J = `j_target` projector is

 𝒫_J = (2J + 1)/(8π²) ∫ dΩ χ_J*(Ω) · R(Ω)

where `χ_J(Ω) = Tr D^J(α, β, γ)` is the character and `R(Ω) =
[D^{½}(α, β, γ)]^{⊗N}` is the tensor-product rotation. libirrep
discretises the SU(2) integral as a tensor-product quadrature —
uniform in α and γ (trigonometric polynomials, exact for sufficient
n), Gauss-Legendre in cos β. Each rotation is applied as a sequence
of N single-qubit updates (O(N · 2^N) per grid point). Validated on
2-spin singlet / triplet states at machine precision; used at 12,
18, and 24-site kagome to confirm the Heisenberg ground state is a
pure J = 0 singlet (`‖P_{J=0}|gs⟩‖² = 1.000000`).

## 4. Performance characteristics

libirrep is written for a single modern CPU core first. Key
performance targets for the research workflow:

- **Space-group application** on a 108-site cluster: the 432-element
 orbit sum must be < 0.5 ms (subdominant to any realistic NN
 forward pass). libirrep delivers **17.6 µs**, 28× under budget.
- **Symmetry-adapted basis** for a 12-site cluster across all six
 p6mm Γ-irreps: 1.07 s end-to-end.
- **Sparse Lanczos** at 24-site kagome: ~74 s per eigensolve (2 min
 including memory allocation). Memory footprint under 1 GB.

Runtime SIMD dispatch (NEON on aarch64, AVX2 reserved for x86_64) is
wired via a function-pointer table populated at first use. Currently
the polynomial cutoff and the cartesian SH batched kernel are NEON-
accelerated; Wigner-D and tensor-product hot paths are scheduled for
NEON in a 1.3.x point release (see [`../TODO.md`](../TODO.md) § M10).

## 5. Validation

### 5.1. Correctness

Every primitive has a unit test; every unit test cross-validates
against a primary-source value where one exists. Summary as of
1.3.0-alpha:

- 28 test suites, ~11 000 assertions total.
- All suites pass under normal, ASan, and UBSan builds.
- 28 public headers are self-contained (each compiles standalone
 under `-Wall -Wextra -Wpedantic -Werror -std=c11`).
- 9 libFuzzer targets run 60 s each in CI.
- Bit-exactness tests for every SIMD kernel against the scalar path.

### 5.2. ABI stability

libirrep tracks a stable C ABI within a major version. Every release
embeds an `irrep_abi_hash()` — SHA-256 over the sorted exported-
symbol set — into the binary, and exposes the same hash via
pkg-config. Consumers can guard against binary drift via a one-line
`-DIRREP_ABI_HASH_EXPECTED` compile-time check. Full policy in
[`DESIGN.md`](DESIGN.md) §5.

### 5.3. Reproducibility

The physics numbers in [`PHYSICS_RESULTS.md`](PHYSICS_RESULTS.md)
were produced by running the `examples/` programs on a single Apple
M2 Ultra. They are deterministic (all random seeds are explicit; no
wall-clock dependence in any kernel), bit-reproducible across macOS
arm64 and Linux x86_64 at the scalar path, and within `10⁻¹⁴`
relative at the NEON / AVX2 paths.

## 6. Relationship to prior art

- **e3nn** (Geiger *et al.* 2022) — libirrep implements the same
 mathematical primitives as the Python `e3nn` package, with
 matched sign conventions (Condon-Shortley + e3nn real-SH sign).
 The mapping is documented side-by-side in
 [`MIGRATION_FROM_E3NN.md`](MIGRATION_FROM_E3NN.md).
- **NetKet** (Carleo *et al.* 2019) — Python / JAX NQS framework.
 libirrep is the C/C++ companion for consumers who need to embed
 symmetry + entropy primitives in performance-critical code; it
 does not ship an NQS ansatz or MCMC sampler, which NetKet does.
- **HPhi** (Kawamura *et al.* 2017) — specialised ED code for quantum
 lattice models. libirrep's symmetry and Lanczos infrastructure is
 a subset of HPhi's capability, but ships as a reusable library
 rather than a monolithic executable.
- **Alps / ALPSCore** — widely-used CMake-based condensed-matter
 framework. Overlap is minimal; libirrep is narrower-purpose and
 has no dependency on Boost or other C++ infrastructure.

The library's differentiating contribution is its combination of
(a) genuine space-group projection, not just translation-symmetry,
(b) a true character-weighted reducer for arbitrary irreps, and
(c) a pure-C ABI suitable for consumption by C, C++, Rust, Julia,
JAX (via `ctypes`), and any other language with a C FFI.

## 7. Release cadence

The 1.3 cycle is structured around the 1.3 release gates:

- **1.3.0-alpha** (current) — full research-substrate code + ED
 validation at 12, 18, 24 sites.
- **1.3.0-beta** — non-Γ Bloch-wave projections + k-point-resolved
 adapted basis; enables full-Hilbert-space block ED.
- **1.3.0-rc1** — downstream integration with
 `spin_based_neural_network` at 108-site kagome NQS scale.
- **1.3.0 final** — paired publication: methods paper (this
 document as the technical backbone) + physics paper (the kagome
 result that actually settles the ground-state-nature question).

## 8. Citation

Please cite libirrep via the BibTeX entry in [`../README.md`](../README.md)
§ Citation.
