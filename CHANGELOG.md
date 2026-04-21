# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`irrep/hamiltonian.h` — on-the-fly Hamiltonian apply operators.**
 Every ED example previously re-implemented the same spin-½ Heisenberg
 `apply_op` callback by hand (~40 LOC each). Promoted to a library
 primitive with three constructors sharing one data-driven apply kernel:
   - `irrep_heisenberg_new` — NN Heisenberg `J·Σ S_i·S_j`
   - `irrep_heisenberg_j1j2_new` — J₁–J₂ on separate NN + NNN bond sets
   - `irrep_xy_new` — XY-only `J·Σ (S^x S^x + S^y S^y)`
 All three produce an `irrep_heisenberg_t` whose `_apply` function has
 the signature `irrep_lanczos_eigvals` expects, so they plug directly
 into the sparse eigensolver. Closes the audit item "physics substrate
 claim is scaffolding, not API." `examples/kagome24_ed.c` migrated to
 the new primitive; other examples will follow in a 1.3.x point
 release. Tested: `tests/test_hamiltonian.c` exercises N=2 closed-form
 eigenvalues, N=4 Heisenberg-ring E₀ = −2J (Bethe-ansatz reference) via
 Lanczos, XY flip-only semantics, J₁–J₂ with J₂=0 / J₂=J₁ boundary
 conditions (E₀ = −¾J at the 3-site equilateral triangle), and input
 validation — 49 assertions.

- **`examples/EXPECTED_OUTPUT.md` — reproducibility reference.** Every
 non-toy example's RNG seed is documented inline in the source; this
 file catalogues the expected numerical output per example so a reader
 building from a fresh clone can verify the library matches this tree.
 Tolerances per-example reflect Lanczos convergence bound, not
 floating-point precision (the latter is machine-epsilon through every
 documented `j` regime).

### Changed

- **Wigner 3j / CG now use Miller two-directional iteration**
 (`src/clebsch_gordan.c`). The Round-2 backward-only Schulten–Gordon
 recurrence was machine-precision through `j ≈ 50` but leaked
 subdominant-solution contamination in the lower tail past `j ≈ 80`
 (sum-rule err `6 × 10⁻⁴` at `j = 80`, diverging at `j = 200`).
 Upgrade to full Miller iteration: forward pass from `j_min`, backward
 pass from `j_max`, splice at `argmax |T_fwd|·|T_bwd|`, rescale,
 normalise by the sum rule, sign-anchor at `j_max`. Measured sum-rule
 precision now machine-level across the entire tested range:

       j =  20 : 0            j =  80 : 4 × 10⁻¹⁶
       j =  50 : 2 × 10⁻¹⁶    j = 120 : 0
       j = 200 : 0

 vs. the prior Racah log-gamma implementation (2 × 10⁻⁹ at `j = 50`,
 NaN past `j ≈ 60`). Sakurai hand values unchanged. No public-API
 changes; no ABI change.

- **Wigner-d rewritten to Jacobi-polynomial form** (`src/wigner_d.c`).
 Replaced the Sakurai (3.8.33) direct-sum implementation — which lost
 precision to catastrophic cancellation past `j ≈ 20` (unitarity
 `≈ 2 × 10⁻³` at `j = 50`, divergent past `j = 60`) — with the Edmonds
 (4.1.23) Jacobi-polynomial form via the NIST DLMF §18.9.1 forward
 three-term recurrence and symmetry reduction to the canonical
 `m ≥ |m'|` region. Measured unitarity is now `≤ 1 × 10⁻¹²` for every
 `j ≤ 80` tested; bounded only by the IEEE-754 `lgamma` overflow limit
 (`j ≈ 170`) past that. Analytic ∂d/∂β updated to match, using the
 Jacobi derivative identity `(d/dx) P_n^{(α,β)} = (n+α+β+1)/2 ·
 P_{n−1}^{(α+1,β+1)}`. No public-API changes. Test suite pins the new
 stability regime across `j ∈ {20, 30, 50, 80}` at 1e-12 / 1e-11.

### Added

- **Non-Γ Bloch-momentum projection** (`irrep/config_project.h`,
 `src/config_project.c`). `irrep_sg_bloch_amplitude` and
 `irrep_sg_bloch_basis` project amplitudes and build symmetry-adapted
 bases at arbitrary Bloch momentum `k = (kx/Lx) b1 + (ky/Ly) b2`,
 using the translation subgroup only. Indices are canonicalised mod
 `(Lx, Ly)` so negative or out-of-range `(kx, ky)` are accepted.
 Enables k-resolved exact diagonalisation; Γ-sector matches A₁
 projection exactly, and the dimensions across all `Lx·Ly` k-sectors
 sum to the full Hilbert-space dimension. Fourier inversion and
 cross-sector orthogonality verified in the test suite.
- **Space-group lattice accessor** (`irrep/space_group.h`).
 `irrep_space_group_lattice(G)` returns the borrowed lattice handle the
 space group was built over; saves callers from threading the lattice
 pointer alongside the space-group pointer when both are needed.

## [1.3.0-alpha] — 2026-04-19

First public tag of the 1.3 cycle. The tested 1.2 core
(spherical harmonics, Clebsch-Gordan, Wigner-D, tensor products, NequIP
message-passing layer, point-group projection, equivariant-NN building
blocks) is preserved intact; the cycle adds a physics substrate pinned
to the Kagome Heisenberg S = ½ ground-state-nature problem (open since
Yan–Huse–White 2011). Seven new headers have landed (`lattice.h`,
`space_group.h`, `config_project.h`, `rdm.h`, `sym_group.h`,
`spin_project.h`, plus a half-integer path in `tensor_product.h`);
end-to-end ED at 12/18/24 sites reproduces published values.

### Added — 1.3 modules

- **2D lattice primitives** (`irrep/lattice.h`, `src/lattice.c`).
 Square, triangular, honeycomb, kagome lattices under PBC. Site
 coordinates, sublattice lookup, NN / NNN bond enumeration (canonicalised
 and deduplicated), primitive / reciprocal vectors, cell translations,
 Brillouin-zone k-grids. The 6×6 kagome target cluster resolves
 to 108 sites, 216 NN bonds, 216 NNN bonds.
- **2D wallpaper-group tables** (`irrep/space_group.h`,
 `src/space_group.c`). p1, p4mm, p6mm site-permutation actions; full
 + inverse permutation caches; `memcpy`-speed application; rejects clusters
 whose torus breaks the full point-group symmetry. Kagome p6mm on 6×6:
 432 elements × 108 sites, ~186 KB per direction.
- **Configuration-space projection** (`irrep/config_project.h`,
 `src/config_project.c`). `P_μ ψ(σ) = (d_μ/|G|) Σ_g χ_μ*(g) ψ(g·σ)` as
 a character-weighted reducer; generic irrep handle plus totally-symmetric
 (A₁) and sign-representation (A₂) shortcuts; orbit enumeration helper;
 `irrep_sg_adapted_basis` — the full symmetry-adapted basis builder
 (character-weighted projection + Gram-Schmidt) that enables symmetry-
 reduced block ED of any Hamiltonian commuting with `G`. Uses an
 orbit-representative filter (skip seeds that are not the minimum of
 their G-orbit) for a ~30× speedup on typical clusters; 12-site kagome
 p6mm symmetry-block ED drops from 33 s to 1.07 s with identical output.
- **Reduced density matrix / entanglement entropy** (`irrep/rdm.h`,
 `src/rdm.c`). Partial trace on `N ≤ 30` sites of arbitrary local
 dimension; cyclic-Jacobi Hermitian eigendecomposition (phase reduction +
 real Givens, convergent in ~10 sweeps at `10⁻¹⁴`); **sparse Lanczos
 eigensolver** (`irrep_lanczos_eigvals`) with callback-based `apply_op`,
 3-term recurrence, ~50-iteration convergence to ground-state energy on
 Heisenberg clusters; von Neumann and Rényi entropies; Kitaev–Preskill
 topological-entanglement-entropy `γ = S_A + S_B + S_C − S_{AB} − S_{BC}
 − S_{AC} + S_{ABC}` — the kagome Z₂-vs-trivial diagnostic (`γ = ln 2`
 vs. `γ = 0`).
- **Symmetric group / Young tableaux** (`irrep/sym_group.h`,
 `src/sym_group.c`). Factorial, permutation sign, lexicographic permutation
 enumeration, hook-length dimension formula, antisymmetric (fermion)
 and totally-symmetric (boson) projectors on tensor-factored states.
 Plancherel identity `Σ dim(λ)² = N!` verified on S₅.
- **Total-J projection** (`irrep/spin_project.h`,
 `src/spin_project.c`). Character-weighted SU(2) integral that projects an
 `N`-spin-½ wavefunction onto a fixed total-J subspace. Single-qubit
 sequential rotation (O(N·2^N) per Euler point); character formula via
 `χ_J = sin((2J+1)ω/2)/sin(ω/2)` with `ω` the total rotation angle.
 Singlet / triplet bipartitions on 2 and 3 sites validated.
- **Half-integer (spinor) tensor products** (extension to
 `irrep/tensor_product.h` + `src/tensor_product_2j.c`). Complex-basis UVW
 tensor product on `irrep_multiset_2j_t`: enumeration, build, forward,
 weighted forward, backward. Spin-½ ⊗ spin-½ = spin-0 ⊕ spin-1 verified
 against Clebsch–Gordan coefficients on singlet / triplet states.
- `CHANGELOG.md` — formal research
 agenda pinning the cycle to real physics work.

### Coverage
Seven new test suites add ~3100 assertions (`test_lattice`: 2341,
`test_space_group`: 330, `test_config_project`: 40, `test_rdm`: 37,
`test_sym_group`: 169, `test_spin_project`: 56, `test_tensor_product_2j`:
29). All 28 test suites green under normal / ASan / UBSan builds. 28
public headers self-contained.

### End-to-end physics regressions
- `examples/kagome_a1_projection.c` — 6×6 × 3 = 108-site kagome orbit sum
 under the full 432-element p6mm space group. Drifts from the
 translation-invariant reference by `1e-14`.
- `examples/heisenberg4_ed.c` — 2×2 square Heisenberg on a 4-site torus.
 E_0 = −2 J (analytical), ground state is a pure J = 0 singlet
 (`‖P_{J=0}|gs⟩‖² = 1.00`), S_VN at the 2-vs-2 cut ≈ 0.837 nats.
- `examples/kagome12_ed.c` — 2×2 × 3 = 12-site kagome Heisenberg on a
 torus. 4096-dim Hilbert space, on-the-fly H-apply, shifted + deflated
 power iterations. Runtime ~2 s end-to-end on M2 Ultra for the full
 ED + projections. Results:
 - E_0 / N = −0.4537 J (ground-state energy, matching Elser 1989 /
 Lecheminant et al. 1997).
 - Ground state is a pure J = 0 singlet (`‖P_{J=0}‖² = 1.00000000`).
 - Ground state lives in a unique 1D p6mm irrep; the library's
 character-weighted projection decomposes |gs⟩ across all six
 irreps (A₁, A₂, B₁, B₂, E₁, E₂) and reports weights summing to
 1.000000000 by completeness.
 - First excited state found by deflating against |gs⟩: also J = 0,
 in a different p6mm sector. Singlet-singlet gap Δ_ss = 0.117 J.
 - Lowest triplet found by seeding in the S_z = 1 sector: confirmed
 J = 1, giving the cluster spin gap Δ_S = 0.383 J.
 - Bipartite S_VN at an up-triangle (3-site) cut ≈ 1.574 nats.
 - Kitaev-Preskill γ on the 12-site ground state, via an annular
 tripartition (trace out sites {6..11}, split kept {0..5} into
 A={0,1}, B={2,3}, C={4,5}): γ ≈ −0.330 nats. Not clean +ln 2 or 0
 — 2×2 torus is too small for a proper annular KP geometry — but
 the full pipeline (7 partial traces, 7 Jacobi eigendecompositions,
 final formula) works end-to-end, which is the machinery-level
 claim the 1.3 scope relies on.
 First libirrep-based kagome ED physics regression, and evidence that
 the 1.3 substrate is fit for purpose at the research
 target's scale.

 Further exercising `kagome12_ed.c`:
 - Static structure factor S(k) over the 4 allowed k-points of the
 2×2 BZ. S(Γ) = 0.000000 to machine precision (a hard consistency
 check: S(Γ) vanishes iff the state is a total-S = 0 singlet).
 S(M_a) = S(M_b) = S(K) = 0.395014 — equal across the three
 non-Γ k-points related by C_6, a non-trivial numerical check that
 the ground state respects the full p6mm point symmetry.
 - Per-bond energy ⟨S_i · S_j⟩_NN = −0.226870 J, consistent with
 E_0 / nb = −0.226870 J to 6 digits.

- `examples/kagome24_ed.c` — 24-site kagome Heisenberg on a 2×4 torus
 via the sparse Lanczos solver. Hilbert space 2^24 = 16 777 216
 (≈ 256 MB per state vector); 3 Lanczos vectors + workspace fit under
 1 GB. Runtime 2.5 min on M2 Ultra (2 × 74 s for ground + triplet
 eigensolves). Results:
 - E_0 = −10.75989 J, **E_0/N = −0.44833 J** (24-site kagome
 literature: −0.438 to −0.443; 2×4 torus published ≈ −0.441).
 - Singlet-singlet gap in Sz=0: Δ_ss = +0.0798 J — the classic
 small-gap VBC / low-lying-singlet-tower signature on kagome.
 - Spin gap Δ_S = +0.26391 J.
 - **Finite-size-scaling** across the 12/18/24 libirrep ED series:
 linear 1/N extrapolation gives **Δ_S(N→∞) ≈ +0.132 J**, which is
 remarkably close to the Yan-Huse-White 2011 DMRG value of 0.13 J
 for the gapped Z₂-spin-liquid picture and distinctly **away from
 zero** (the gapless Dirac SL hypothesis would predict Δ_∞ = 0).
 (E_0/N)(N→∞) ≈ −0.441 J matches the published thermodynamic limit.
 Three data points with no error bars — not definitive, but
 infrastructure-wise this is a concrete end-to-end chain from a
 stable C library to a physics observable relevant to the current
 kagome ground-state-nature debate.

- `examples/kagome18_ed.c` — 18-site kagome Heisenberg on a 2×3 torus,
 the next non-square cluster up from the 12-site one. Hilbert space
 2^18 = 262144 (≈ 4 MB per state vector); the cluster breaks C_6 so
 only p1 (translations) applies. Power iteration on on-the-fly sparse
 H-apply (36 NN bonds). Runtime 27 s on M2 Ultra. Results:
 - E_0 = −8.04719493 J, E_0/N = −0.44706639 J/site, within the
 published 18-site kagome range of −0.43 to −0.44 J/site
 (Waldtmann et al. 1998, Läuchli 2011).
 - Pure J = 0 singlet (‖P_{J=0}|gs⟩‖² = 1.00000000).
 - Bipartite S_VN at an up-triangle cut ≈ 1.587 nats.
 - Per-bond consistency ⟨S_i·S_j⟩_NN = −0.2235 J matches E_0/nb to
 six digits.
 Extended physics on the same cluster:
 - Spin gap via S_z = 1 power iteration: Δ_S = +0.28349 J
 (12-site value: 0.3827 J). The gap **dropping** with cluster
 size is the observation at the heart of the Kagome ground-state-
 nature debate (consistent with a gapless spin-liquid picture
 or slowly-converging gapped picture; not yet resolvable from
 two data points).
 - Static structure factor S(k) across the 6 k-points of the 2×3
 BZ. S(Γ) = 0 to machine precision (J_tot = 0 consistency check);
 maximum at (1,2) with S(k) = 0.6314.
 Gives a second data point for finite-size scaling (12-site: −0.4537;
 18-site: −0.4471); the downward trend toward the thermodynamic limit
 is correct for kagome AFM.

- `examples/kagome12_symmetry_ed.c` — full symmetry-adapted block ED on
 12-site kagome. For each of the six p6mm Γ-irreps (A₁, A₂, B₁, B₂,
 E₁, E₂), build an orthonormal basis via character-weighted projection
 + Gram-Schmidt, then diagonalise H within that sector via `rdm.h`'s
 Jacobi. Sector dims 76–420 (vs the 4096 full Hilbert space, 10–50×
 reduction). Sum = 1072 = dim V_Γ by Burnside on the translation
 subgroup. Global minimum reproduces E_0 = −5.44487522 J in the B₁
 sector, matching the power-iteration value bit-for-bit. This is the
 "beyond naive ED" path the 1.3 scope commits to for 18- and
 24-site kagome clusters where full ED is intractable but per-sector
 block sizes remain small.

### Added (post-1.2.0 — lands in the first 1.2.x point release or folds into 1.3)
- `irrep/multiset_2j.h` + `src/multiset_2j.c` — doubled-integer multiset
 type for mixed integer + half-integer spin content. Parser accepts
 both `"1x0e"` (integer-l, e3nn-compatible) and `"1x1/2o"`
 (half-integer) forms; stores `two_j` throughout so spin-½ has
 `two_j = 1`. The companion `irrep_time_reversal_square_sign_2j(m)`
 correctly returns `−1` on any multiset carrying half-integer spin
 (Kramers degeneracy; Sakurai §4.4), closing the scope gap noted in
 1.2.0's audit. `irrep_multiset_2j_append` guards the `total_dim`
 accumulator against `int` overflow via an `INT_MAX / (two_j + 1)`
 check before multiplication.
- `tests/test_multiset_2j.c` — 48 assertions across nine scenarios:
 integer-only parsing, half-integer parsing, mixed content, empty
 input, eight malformed-input rejections, format round-trip,
 direct `_append` validation with error codes.
- `tests/test_downstream_compat/torque_net_tp_paths/vectors.json` +
 `generate_golden.c` — five fixed `(m_i, m_j, r̂)` configurations
 with pre-computed cartesian reference values for the five torque-net
 basis terms (T1 through T5). The upgraded `test_downstream_compat`
 now asserts bit-exact agreement through the UVW tensor-product
 paths with the documented `√3 · (1o⊗1o→0e)` and `√2 · (1o⊗1o→1e)`
 prefactors, to `10⁻¹²` residual. 37 golden-vector assertions active.
- Windows / MinGW-w64 path in the Makefile (DLL build with
 `--out-implib`, `IRREP_BUILDING_DLL` preprocessor define) and a
 Windows smoke job in `.github/workflows/ci.yml`. Cross-compilation
 for the release tarballs is handled by CI fan-out across native
 runners (macOS arm64/x86_64, Linux x86_64, Linux arm64, Windows
 x86_64); see `.github/workflows/release.yml`.

---

> **Note on git tags.** Releases `1.0.0`, `1.1.0`, and `1.2.0` were internal
> development cycles prior to the public repository. Their entries below
> are retained as narrative history so the scope and precision trajectory
> of the project is transparent, but they have no corresponding git tag
> in this repository. Only `v1.3.0-alpha` and later carry tags.

---

## [1.2.0] — 2026-04-19

Superset of 1.2.0-rc1 plus the post-rc audit fixes, the four-point-group
projection module (C₄ᵥ, D₆, C₃ᵥ, D₃), the 80× `irrep_pg_project` speedup,
the Legendre-grid SH batch speedup, two integration examples, and the 1.3
roadmap doc.

### 1.2.0-rc1 audit fixes + C3v/D3 point groups

Internal hardening from the 1.2.0-rc1 audit sweep plus the promised C3v/D3
point-group support. No public API changes on the audit side; two new
values of `irrep_point_group_t` become build-able.

### Performance
- `irrep_pg_project` is now **≈80× faster** on the standard
 `2x0e + 1x1o + 1x2e` multiset across all four supported groups. At table-
 build time we pre-compute the real-basis `D^l(g) = U · D_complex · U†`
 matrices for every element and every `l ≤ IRREP_PG_CACHED_L_MAX = 4`
 (covers every common NequIP / MACE `l_sh` value); projection becomes a
 pure matrix-vector multiply with a cache hit. Larger `l` transparently
 falls back to per-call recomputation. Cache footprint is ~2–16 KB per
 table depending on group order — negligible for long-lived tables.
 Measured on Apple M2 Ultra:
 - C₃ᵥ: 19 µs → 224 ns (83×)
 - D₃: 19 µs → 224 ns (85×)
 - C₄ᵥ: 26 µs → 296 ns (88×)
 - D₆: 38 µs → 443 ns (85×)

### Added
- `examples/symmetric_nqs_projection.c` — end-to-end integration: parse a
 NequIP layer from a spec string, forward-apply on a 4-node square lattice,
 project the output onto the C₄ᵥ A₁ (totally symmetric) irrep per node.
 Demonstrates the full 1.2 surface (spec parser, NequIP layer, point-group
 projector) in ~100 LOC. Verifies projector idempotence to ~4e-18.
- `examples/torque_net_tp_paths.c` — worked coord-doc example mapping the
 torque-net's cartesian `a · b` and `a × b` to `irrep_tp_apply_uvw` paths
 via `1o ⊗ 1o → 0e` and `1o ⊗ 1o → 1e`. Bit-exact round trip to ~2e-16.

### Fixed (docs)
- `docs/tutorials/05_tensor_products.md` §7 — corrected two errors in the
 worked-example table: the cross-product path lands in `1e` (axial), not
 `1o` (parity is multiplicative, `odd × odd = even`), and the `(1,1,1)`
 path's sign convention in libirrep is `+(1/√2) · (a × b)` (verified
 empirically against `examples/torque_net_tp_paths.c`), not `−(1/√2)`.

- `tests/fuzz/fuzz_nequip_from_spec.c` and `tests/fuzz/fuzz_pg_project.c` —
 libFuzzer harnesses for the spec-string parser and for all four point-
 group projectors. Both wired into the CI fuzz-smoke job for 60-second
 runs alongside `fuzz_cg`, `fuzz_multiset_parse`, `fuzz_tp_apply`, and
 `fuzz_nequip_apply_forces`.
- `benchmarks/bench_point_group.c` — per-group projection throughput on
 the canonical `2x0e + 1x1o + 1x2e` multiset. Numbers scale as
 `|G| · Σ_l (2l+1)²` as expected: C₃ᵥ/D₃ ~19 µs, C₄ᵥ ~26 µs, D₆ ~38 µs
 per call on Apple M2 Ultra. Lands in the baseline as the first
 per-group cost reference.
- `IRREP_PG_C3V` and `IRREP_PG_D3` — triangular-lattice symmetry groups
 (order 6, 3 irreps: A₁, A₂, E). Character tables per Bradley-Cracknell.
 C₃ᵥ has three improper σᵥ reflections at 120°; D₃ has three proper C₂
 axes at 120° in the xy-plane. The two are isomorphic as abstract groups
 but differ on parity-odd inputs: `_reduce` on `1x1o` under C₃ᵥ gives
 `A₁ + E` (z survives the σᵥ parity flip onto A₁), under D₃ gives
 `A₂ + E` (z sees only the C₂ sign flip). `1x1e` flips this for C₃ᵥ
 (→ A₂ + E) but not D₃ (parity is irrelevant in purely-proper groups).

### Fixed
- `irrep_nequip_layer_from_spec` — integer options (`sh`, `radial`,
 `polynomial(p)`) now reject values that would wrap on the `long → int`
 downcast (`v > INT_MAX` guard). Previously `sh=999999999999999` silently
 became a negative `int`.
- `libirrep.pc` now emits `-Wl,-rpath,${libdir}` in the `Libs` line.
 Consumers compiling dynamic binaries via
 `cc $(pkg-config --cflags --libs libirrep)` no longer hit
 "Library not loaded: @rpath/liblibirrep.1.dylib" at runtime. Static
 builds (`pkg-config --static`) use `Libs.private`.

### Hardened
- Parser tests extended to 65 assertions: explicit coverage for
 `polynomial(0)`, `polynomial(-3)`, `sh=-1`, `radial=0`, `r_cut=0.0`,
 INT_MAX overflow, empty hidden_in/out, empty `[]` block, and the
 duplicate-option last-wins behaviour.
- Point-group reduction tests extended to 252 assertions with hand-
 computed Bradley-Cracknell decompositions: `1x1o → A₁+E` and
 `1x2e → A₁+B₁+B₂+E` under C₄ᵥ; `1x0e → A₁` and `1x1o → A₂+E₁`
 under D₆.
- NEON SH bit-exactness sweep extended to 4099 assertions across
 `l_max ∈ {0..6}`, exact poles, near-poles on both sides of the 1e-14
 `r_xy` threshold, and equator inputs.
- Benchmark results JSON now carries a `_meta` sentinel with
 `cpu_model`, `triple`, and `stamp`. `scripts/perf_compare.sh` prints
 both CPU models and warns on mismatch — suppresses false-positive
 regressions when CI runs on a different-generation core than the
 frozen baseline.
- `fuzz_nequip_apply_forces` wired into `.github/workflows/ci.yml` for
 60-second runs (matches the other three libFuzzer targets).

### Docs
- README quickstart extended with a 1.2-era `irrep_nequip_layer_from_spec`
 example and an ABI-drift check snippet using `pkg-config --variable=abi_hash`.

## [1.1.0] — 2026-04-18

Cut for downstream consumer feedback (spin-based neural network vendoring).

### Added
- Initial project skeleton (M1): Makefile, CMake, public headers, stub sources,
 empty TAP-style test harness, documentation scaffolding, release scripts.
- **M2 — SO(3) / SU(2)**: all rotation conversions (quat ↔ rot matrix ↔ Euler
 ZYZ ↔ axis-angle), Shepperd's quaternion-from-matrix, Rodrigues exp,
 quaternion-routed log (π-safe), Hamilton-product composition, Shoemake
 uniform quaternion sampling, SLERP, geodesic distance, chordal Fréchet mean,
 shortest-arc from two vectors, validity check. SU(2) Pauli matrices, SU(2)
 ↔ quaternion ↔ SO(3) conversions (double cover), SU(2) apply, compose,
 inverse, 2×2 matrix exponential via scaling-and-squaring. Tests pass at
 `1e-12` (277 SO(3) + 32 SU(2) = 309 assertions).
- **M3 — Spherical harmonics**: associated Legendre via stable three-term
 recurrence with Condon-Shortley phase; complex `irrep_sph_harm`, real
 `irrep_sph_harm_real`, and cartesian `irrep_sph_harm_cart` /
 `irrep_sph_harm_cart_all` up to `l = IRREP_L_MAX = 16`. Real-SH convention
 matches e3nn / Wikipedia (`Y_{1,+1} ∝ +x`). Complex↔real basis-change matrix
 `irrep_sph_harm_complex_to_real` (unitary to `1e-12`). Finite-difference
 gradient `irrep_sph_harm_cart_grad` (analytic solid-harmonic form deferred
 to M10). `_f32` single-precision wrappers. Tests (386 assertions) cover
 Condon-Shortley identity `Y_l^{-m} = (-1)^m Y_l^{m*}`, sum rule, addition
 theorem against Legendre, pole behaviour, cartesian-vs-polar equivalence,
 and radial-derivative-is-zero tangency.
- **M4 — Clebsch-Gordan + Wigner 3j**: `irrep_cg_2j` via the Racah single-sum
 in log-gamma space; integer wrapper `irrep_cg`. Wigner 3j symbols follow
 the standard phase / normalization from CG. Selection-rule / parity checks
 zero out invalid arguments exactly. Dense five-dimensional `cg_table_t`
 cache (`irrep_cg_table_build` + `irrep_cg_lookup`) populated at build time;
 space is ~3 MB for `j_max = 4`, ~80 MB for `j_max = 8`. `irrep_cg_block`
 fills a flat `(2j1+1)(2j2+1)(2J+1)` buffer. Tests (1391 assertions)
 cover Sakurai hand values, completeness `Σ_J |CG|² = 1`, orthogonality
 between distinct J columns, 3j↔CG relation, 3j cyclic column symmetry,
 cache/direct agreement, and large-j finite-value stability through j = 10.
- **M5 — Wigner-D matrices**: small-d `irrep_wigner_d_small_2j` via Sakurai's
 direct-sum formula evaluated in log-gamma form. Full D factorized as
 `D = e^{-i m' α} · d(β) · e^{-i m γ}`. Scalar, per-matrix
 (`irrep_wigner_D_matrix`, `irrep_wigner_d_matrix`), and multiset
 block-diagonal (`irrep_wigner_D_multiset`) builders. Analytic `∂d/∂β` via
 term-wise differentiation. Tests (1065 assertions) cover closed-form values
 at j=½ and j=1 (full 3×3 d-matrix), `d(0) = I` and `d(π)` reflection with
 sign `(−1)^{j−m'}`, unitarity `D D† = I` at j ≤ 4 and j = 10, composition
 `D(R₁) D(R₂) = D(R₁·R₂)` against SO(3), analytic ∂d/∂β vs centered FD, and
 multiset block-diagonal placement with cross-block zero entries.
- **M12 / M13 / M14 partial** — sanitizer-clean, release artifacts, working
 examples.
 - `make asan` and `make ubsan` pass the full 14-suite test matrix with zero
 findings (after fixing a stack-buffer overflow in the Wigner-D test).
 - `scripts/build_release.sh 1.0.0` produces `release/1.0.0/{include, lib,
 pkgconfig, VERSION, ABI_HASH, CHECKSUMS}` plus a per-triple tarball.
 Baseline ABI hash committed.
 - Three examples promoted from stubs to real demos: `spin_half_rotation`
 verifies Berry phase `D^{1/2}(2π, ẑ)|↑⟩ = −|↑⟩`; `equivariant_mlp`
 stacks Linear → RMS-norm → Gate → Linear and verifies
 `f(R·x) = R·f(x)` at `3e-17`; `nequip_message` composes cartesian SH,
 Bessel RBF, polynomial cutoff, and a weighted tensor product to produce
 one message.
 - Tensor product now rejects `l_a + l_b + l_c` odd paths in both
 `irrep_tp_build` and `_enumerate_paths` (they produce pure-imaginary
 real-basis couplings not supported in v1.0).
 - Benchmarks run and write JSON: `sph_harm_cart_all_l4` ~2.8 ns/output
 element, `wigner_D_matrix_j4` ~134 ns/element, on macos-arm64.
- **M9 — Equivariant-NN primitives**:
 `irrep_linear_build` matches `(l, parity)` pairs across input and output
 multisets and allocates one `mult_out × mult_in` weight block per matched
 pair (unmatched labels contribute zero). `irrep_linear_apply` /
 `_backward` apply per-block channel mixing that preserves equivariance.
 `irrep_norm_rms` / `_backward` normalize each `(term, copy)` block by its
 own RMS and scale by a learnable scalar; analytic backward derived
 term-by-term. `irrep_gate_apply` multiplies each block by a caller-supplied
 scalar (typically the sigmoid of a scalar feature). Tests (39 assertions):
 weight-count formula, scalar-only mixing on mismatched labels, full
 real-basis SO(3) equivariance through `l=1`, analytic backward matches
 centered FD for linear and RMS-norm at `1e-6`.
- **M8 — Recoupling, radial basis, quadrature**:
 - `irrep_wigner_6j_2j` via the Racah single-sum in log-gamma form; 9j
 via Edmonds' sum over k of three 6j products; Racah W is 6j with the
 standard `(−1)^{j1+j2+j3+J}` phase.
 - `irrep_rbf_bessel` / `_all` (DimeNet form with r→0 sinc limit),
 `irrep_rbf_gaussian` / `_grid`, `irrep_cutoff_cosine` + analytic
 derivative, `irrep_cutoff_polynomial` (NequIP form) + analytic
 derivative `−p(p+1)(p+2)/(2 r_cut) · u^{p−1} · (1−u)²`.
 - `irrep_gauss_legendre(n, nodes, weights)` via Newton iteration on
 `P_n(x)` starting from a Chebyshev initial guess (any n).
 - `irrep_lebedev_size` / `_fill` for orders 3, 5, 7 (6, 14, 26 points)
 using the a₁/a₂/a₃ generators; higher orders deferred to the full
 Lebedev-Laikov 1999 data import.
 - Also closes the M3 orthonormality TODO via a Gauss-Legendre × uniform-φ
 tensor-product quadrature of `∫ Y_l^m Y_{l'}^{m'}* dΩ = δδ`.
 - Tests (13 + 54 + 168 = 235 new assertions, plus 256 added to the SH
 orthonormality block). Covers: Edmonds 6j hand values incl. the
 `j_6 = 0` reduction, 6j column-permutation symmetry, 9j reduction when
 the last column is zero, Bessel orthonormality on `[0, r_cut]` to `1e-8`,
 cutoff continuity at `r_cut`, analytic-vs-FD derivatives, Gauss-Legendre
 exactness on polynomials through `2n − 1`, Lebedev weight-sum = 1, point
 unit-sphere check, `⟨x²⟩_{S²} = 1/3` exact, `⟨x⁴⟩ = 1/5`, `⟨x²y²⟩ = 1/15`.
- **M7 — Tensor products (real-basis, uuu channel mode)**:
 `irrep_tp_enumerate_paths` scans all (i_a, i_b, i_c) triples for triangle +
 parity validity. `irrep_tp_build` requires matching multiplicities per
 selected path and precomputes a real-basis coupling tensor
 `W[n_a,n_b,n_c] = Σ U_c[n_c,m_c]·CG[m_a,m_b,m_c]·conj(U_a[n_a,m_a])·conj(U_b[n_b,m_b])`
 from the complex CG block and the M3 complex-to-real basis matrices —
 so apply works directly on real-basis SH feature vectors. Forward + weighted
 (one scalar per path) + backward for a, b, and w; batched variants stride
 by per-multiset total_dim. Tests (34 assertions) cover path enumeration
 with parity filtering, the `1×1→0e` real-basis scalar
 `c = −(a·b)/√3`, multiplicity-2 copy pairing, SO(3) equivariance with
 D_real = U·D_complex·U† up to l=2, finite-difference backward matching for
 grad_a, grad_b, and grad_w, batched ≡ per-sample, and build-time rejection
 of mismatched multiplicities or parity-violating paths.
- **M6 — Multiset algebra + parity + time reversal**:
 `irrep_multiset_parse` accepts the e3nn canonical form (`"1x0e + 2x1o + 1x2e"`),
 `irrep_multiset_format` round-trips, `irrep_multiset_simplify` sorts by
 (l, parity) and merges like terms, `irrep_multiset_direct_sum` concatenates
 then simplifies, `irrep_multiset_append` grows amortized, `block_offset`
 is O(num_terms). `irrep_parity_filter_paths` compacts path lists in place
 by the `parity_a · parity_b = parity_c` rule. Time-reversal as U·K: exposes
 the linear factor U with `U_{m',m} = (−1)^{j+m'} δ_{m,−m'}`, confirmed
 T² = +I for integer-l and T² = −I for half-integer (Kramers). Tests
 (64 + 12 + 261 = 337 assertions) cover parser error paths, format
 round-trip, simplify idempotency, parity filter truth table, explicit
 spin-½ matrix matching i σ_y, and multiset block-diagonal placement.

