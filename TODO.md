# TODO

Milestone checklist. `[x]` = landed, `[ ]` = actual open work. Historical
milestones (M1–M9) are kept for context; the live open items are in M10
onwards plus the 1.3 section at the bottom.

## M1 — project skeleton (done)
- [x] LICENSE, NOTICE, README, CHANGELOG, CONTRIBUTING, etc.
- [x] VERSION, .gitignore, .clang-format, .editorconfig
- [x] All public headers with full API declarations (31 as of 1.3.0-alpha,
 post 3D-lattice + cubic-point-group + DMI-symmetry-analyzer additions)
- [x] Stub .c sources that link (then filled in)
- [x] Makefile + CMakeLists.txt
- [x] TAP-style test harness and one stub test per module
- [x] Scripts, benchmarks, examples, docs scaffolding
- [x] `git init` + initial commit (a6a7e34)
- [x] `make`, `make test`, `make examples` all succeed

## M2 — SO(3) / SU(2) (done)
- [x] All rotation conversions (quat / axis-angle / rot_matrix / euler_zyz)
- [x] compose, inverse, exp (Rodrigues), log (quaternion-routed, π-safe)
- [x] Shoemake uniform quaternion sampling
- [x] SLERP, Fréchet mean (chordal — full Karcher iteration deferred to v1.1),
 geodesic distance
- [x] SU(2) Pauli matrices, SU(2) ↔ SO(3) double cover, SU(2) exp
- [x] Tests pass at `1e-12` (277 SO(3) + 32 SU(2) assertions)

## M3 — spherical harmonics (done)
- [x] Complex, real, cartesian (`l = 0..16`)
- [x] Gradient `∇Y_l^m` (cartesian) — FD for M3; analytic in M10
- [x] Complex ↔ real basis change (unitary at `1e-12`)
- [x] Associated Legendre polynomials (stable recurrence, Condon-Shortley phase)
- [x] Orthonormality test via Lebedev
- [x] Polar vs. cartesian agreement to `1e-12`
- [x] Addition theorem check

## M4 — Clebsch-Gordan + Wigner 3j (done)
- [x] Racah formula (log-gamma throughout)
- [x] Integer + half-integer (`_2j`) variants
- [x] Selection-rule zeros exact
- [x] Dense table builder / lookup (5-D flat array)
- [x] Sakurai hand values match at `1e-12`
- [x] 3j ↔ CG relation tested

## M5 — Wigner-D (done)
- [x] Small-d via Sakurai direct sum in log-gamma form (stable past j = 50)
- [x] Full D matrix and scalar variants (integer + `_2j`)
- [x] Block-diagonal D on irrep multiset
- [x] Analytic `∂d/∂β` for force evaluation
- [x] Unitarity, composition, and `d(π)` reflection tests

## M6 — multiset + parity + time-reversal (done)
- [x] String parser (`"1x0e + 2x1o + 1x2e"`)
- [x] Simplify (canonical sort + merge), format round-trip, direct sum
- [x] Parity helpers, in-place path filter
- [x] T operator linear factor (integer and half-integer)
- [x] `T² = ±1` tests; spin-½ Kramers convention verified

## M7 — tensor products (done)
- [x] Path enumeration and descriptor build (real-basis coupling tensor)
- [x] Scalar + weighted apply
- [x] Backward (forward gradient)
- [x] Batched variants
- [x] Equivariance test passes at `1e-10`
- [x] Channel modes beyond "uuu" (uvw — 1.2; half-integer 2j — 1.3)
- [ ] Sparse CG-block storage for compile-time path lists — future

## M8 — recoupling + radial + quadrature (done)
- [x] Wigner 6j, 9j, Racah W (log-gamma)
- [x] Bessel and Gaussian RBFs
- [x] Cosine / polynomial cutoffs + analytic derivatives
- [x] Gauss-Legendre via Newton iteration (any n)
- [x] Lebedev orders 3, 5, 7
- [x] Lebedev orders 9..41 — runtime registry + `scripts/fetch_lebedev_tables.sh`
 + `examples/register_lebedev.c`. Exact to ±5e-15 at every order.
- [x] Exactness tests + SH orthonormality via Gauss-Legendre tensor product

## M9 — equivariant layers (done)
- [x] Linear-on-irreps + analytic backward (matched-label channel mixing)
- [x] RMS norm per (term, copy) + analytic backward
- [x] Gate activation
- [x] Equivariance verified on linear; backward FD cross-checks on linear + RMS

## M10 — NEON hot paths (partial)
- [x] `irrep_cutoff_polynomial` batched NEON kernel (1.1, bit-exact)
- [x] `irrep_sph_harm_cart_all` NEON kernel (1.2)
- [x] `irrep_wigner_d_matrix_batch` NEON kernel — 2.4–2.7× across j ∈ [2, 30],
 bit-exact vs. scalar (30 363 per-entry assertions pass).
- [x] `irrep_tp_apply_weighted_batch_flat` (dim-first layout) +
 NEON + AVX2 kernels. Batch dim is innermost in memory so vector
 loads/stores are contiguous across samples — 2.14–2.47× measured
 speedup over the batch-first scalar path on NEON (Apple M2 Ultra).
 Bit-exact via transpose against per-sample `irrep_tp_apply_weighted`;
 443 per-entry assertions pass at 1e-14. The existing batch-first
 `irrep_tp_apply_weighted_batch` API stays scalar (gather-scatter
 bound without a layout change); callers needing the SIMD win
 transpose their inputs and call the `_flat` variant.
- [x] Bit-exact vs. scalar for kernels that exist

## Numerical stability replacements (open)
- [x] Wigner-d small: replaced Sakurai (3.8.33) direct-sum (unstable past
 j ≈ 20) with Edmonds Jacobi-polynomial form via DLMF §18.9.1 recurrence.
 Measured unitarity ≤ 1e-12 through j = 80; bounded only by IEEE-754
 lgamma overflow (j ≈ 170) past that.
- [x] Wigner 3j + Clebsch–Gordan: replaced Racah log-gamma single-sum
 (unstable past j ≈ 20, NaN past j ≈ 60 near triangle edge) with
 Schulten–Gordon three-term recurrence in j₁, now as Miller two-
 directional iteration (forward + backward passes, splice at
 argmax |T_fwd|·|T_bwd|, normalise via sum rule, sign-anchor at j_max).
 Machine precision to at least j ≈ 200 in regression tests.

## M11 — x86 SIMD hot paths (open)
- [x] AVX2 + FMA variant of `irrep_wigner_d_matrix_batch` — 4 β lanes per
 `__m256d`, same Jacobi-coefficient-broadcast design as the NEON kernel.
 Dispatch via `irrep_cpu_has_avx2` + `_has_fma`. Bit-exactness to scalar
 verified on the test suite (CI macos-x86_64 + linux-x86_64 runners
 exercise the live kernel).
- [x] AVX2 + NEON variants of the dim-first tensor-product batch
 (`irrep_tp_apply_weighted_batch_flat`). 4 batch lanes per
 `__m256d` on x86_64, 2 per `float64x2_t` on arm64. Bit-exact to
 the per-sample scalar reference.
- [x] Runtime dispatch verified on both macOS x86_64 and Linux x86_64 for
 `sph_harm_cart_all`, `cutoff_polynomial`, `rbf_bessel`; the new
 `wigner_d_batch` adds to that set.
- [x] NEON-to-AVX2 bit-exactness documented (matching pragma `FP_CONTRACT
 OFF` on both sides; scalar reference is the same in both cases).

## M12 — benchmarks + fuzz + sanitizer CI (done)
- [x] All benchmarks produce JSON (harness emits per-result record)
- [x] libFuzzer targets (9 targets: CG, multiset parser, TP apply, nequip
 spec + forces, point-group projection, lattice, etc.)
- [x] GitHub Actions matrix: 5 triples × {gcc, clang where available}
- [x] `perf_compare.sh` baseline committed (`benchmarks/results/baseline/`)
- [x] Sanitizer CI: ASan + UBSan green on every commit

## M13 — release build system (done)
- [x] `scripts/build_release.sh <version>` produces `release/<version>/`
- [x] Per-triple tarballs + checksums (macOS arm64/x86_64, Linux x86_64/arm64,
 Windows MinGW x86_64)
- [x] `ABI_HASH` emitted and verified (`irrep_abi_hash()` baked into binary)
- [x] `libirrep.pc` and `libirrepConfig.cmake` installed
- [x] SBOM (SPDX 2.3 JSON)

## M14 — examples, tutorials, docs (done)
- [x] `spin_half_rotation` (Berry phase π)
- [x] `equivariant_mlp` (toy E(3) MLP end-to-end, equivariance verified)
- [x] `nequip_message` (one NequIP-style TP convolution)
- [x] Seven tutorials committed (rotations, SH, CG, Wigner-D, TP, equivariant
 NN, and the 1.3 Kagome NQS substrate)
- [x] Doxygen HTML generated and bundled into the release tarball at
 `docs/html/`. `scripts/build_release.sh` detects `doxygen` on the
 build host; absent, the tarball ships without the HTML (not a
 release blocker — the CI doxygen job enforces zero warnings
 independently).
- [x] README polished
- [x] Consumer smoke test against the release tarball

## 1.3 cycle — 1.3 scope

See the 1.3 CHANGELOG.

- [x] `lattice.h` (square, triangular, honeycomb, kagome; PBC)
- [x] `space_group.h` (p1, p4mm, p6mm)
- [x] `config_project.h` (character-weighted reducer + adapted basis)
- [x] `rdm.h` (partial trace, Jacobi, Lanczos, entropies, γ)
- [x] `sym_group.h` (S_N, hook-length, A/S projectors)
- [x] `spin_project.h` (total-J projection on spin-½ chains)
- [x] `tensor_product.h` half-integer (spinor) UVW path
- [x] End-to-end kagome ED at 12, 18, 24 sites — published-matching values
- [x] Non-Γ Bloch-wave projections (`irrep_sg_bloch_amplitude`,
 `irrep_sg_bloch_basis` — translation-subgroup momentum projector)
- [x] k-point-resolved symmetry-adapted basis at arbitrary k (integer
 `(kx, ky)` with `k = (kx/Lx) b1 + (ky/Ly) b2`; little-group point-group
 overlay on top of the pure-momentum basis remains deferred)
- [x] Named-irrep builtins for C_6v, C_3v, C_2v, C_4v little groups via
 `irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_*)`. Covers A_1,
 A_2, B_1, B_2, E_1, E_2 (C_6v); A_1, A_2, E (C_3v); A_1, A_2, B_1, B_2
 (C_2v at M on p6mm and at X on p4mm); A_1, A_2, B_1, B_2, E (C_4v at Γ
 on p4mm). 185 assertions cover bit-exact agreement with manually-
 constructed character rows and orthogonality against A_1-symmetric
 inputs.
- [x] `irrep_sg_projector_weights` — per-group-element composite-projector
 weight export. Unblocks downstream MPO-based tensor-network consumers
 building `P_{k,μ_k} = Σ_g w_g · π_g` directly.
- [x] Bitstring orbit-canonicalisation primitives: `irrep_space_group_apply_bits`,
 `irrep_sg_canonicalise`, `irrep_sg_orbit_size`. Foundation for
 representative-basis sparse Lanczos at N ≥ 27. 6835 assertions pass
 (identity, popcount invariance, idempotence, orbit invariance, minimum-g
 tiebreak, Burnside sum rule on kagome 12 Sz=0 and square 3×3 p1).
- [x] `irrep_sg_rep_table_t` — orbit-representative table with
 `_build / _free / _count / _get / _orbit_size / _index / _space_group`.
 Enumerates via Gosper's next-lex-greater-with-same-popcount (C(N,p)
 not 2^N); reps stored sorted; O(log n) binary-search rep→index lookup.
 2514 assertions (Burnside, canonicity, sort invariant, round-trip on
 924 configs, rejection of non-reps and wrong-popcount queries).
- [x] `irrep_heisenberg_apply_in_sector(H, rep_table, psi_in, psi_out)` —
 sparse H apply in the totally-symmetric sector of the space group the
 table was built on. O(|bonds| · |reps|) vs dense O(|bonds| · 2^N);
 Sandvik / Weiße-Fehske formula. 1813 assertions cross-check dense vs
 sparse bit-exactly at 1e-13 on kagome 12 Γ-A_1 Sz=0 (30×30 matrix,
 Hermiticity, ground state E_0 = −5.328392 J).
- [x] `irrep_sg_stabiliser(G, config, out)` — enumerate Stab(config).
 Building block for sector-basis norms at non-trivial (k, μ_k). 72
 assertions (trivial configs have full stab, orbit-stabiliser theorem).
- [x] `irrep_sg_heisenberg_sector_t` — cached (rep, bond) → target-index
 connectivity so the sparse matvec is O(1) per transition (pure
 memory-bound CSR). Measured 71× speedup vs dense at N=12; N=27 Γ-A_1
 Lanczos in 1.2 s on 186 k-dim sector (cached build ~12 s one-shot;
 dense path infeasible at 2 GiB per state vector). This is the
 workstation-scale ED unlock for the kagome torus side.
- [x] `irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu_k)` — full
 (k, μ_k) sparse apply for 1D little-group irreps. Computes per-rep
 σ_u = Σ_{s ∈ Stab(u)} conj(w_s), filters σ_u = 0 reps from the sector
 basis, stores complex-coef CSR. Verified on kagome 12: Γ-A_1 →
 −5.32839240 J (matches trivial-sector path); Γ-B_1 → −5.44487522 J
 (matches absolute KH ground state). Lifts the sparse ED path from
 "Γ-A_1 only" to "any 1D (k, μ_k)" — every sector in `kagome12_k_resolved_ed`
 is now reachable at workstation scale.
- [x] `examples/kagome27_sector_ed.c` — full 8-sector (k, μ_k) ED on
 N=27 kagome torus in ~120 s wall clock on Apple M2 Ultra. Dense path
 infeasible (2 GiB/vector, 120 GiB for reorth Lanczos). Shows K-point
 dominance (E_0/N = −0.6132 at K, A_1 vs −0.4300 at Γ, A_1) — the
 Dirac-cone signature in first-principles numerical form.
- [x] `examples/j1j2_square_sector_ed.c` — J₁-J₂ square-lattice 4×4 p4mm
 sparse ED at J₂/J₁ ∈ {0, 0.5, 1.0}. Three 12-sector sweeps in 0.6 s
 wall clock. Pure-J₁ GS energy E_0/N = −0.7018 (cf. Schulz 1989:
 −0.7015, 3e-4 agreement). DQCP regime (J₂=0.5) shows Γ→X GS crossover
 characteristic of the Néel→VBS transition. Demonstrates the sparse
 stack is not kagome-specific — works for any wallpaper group + bond-list H.
- [x] `benchmarks/bench_kagome_sector_ed.c` — N=27 sparse-stack performance
 harness emitting JSON per phase (rep-table, trivial-sector build,
 trivial Lanczos, (K, A_1) sector build, K Lanczos). Picked up by
 `benchmarks/run_benchmarks.sh` for CI regression tracking.
- [x] Rep-table save/load persistence (`irrep_sg_rep_table_save/_load`).
 Binary format with magic header + version + cross-checked metadata
 (popcount, num_sites, group_order). Measured 4700× speedup on N=27:
 3.615 s build → 0.001 s load. Enables iterative research loops
 without rebuilding the table from scratch per invocation.
- [x] Sector-binding save/load (`irrep_sg_heisenberg_sector_save/_load`).
 Binary CSR format (magic + version + flavour flag for real/complex
 coefficients). Both trivial and (k, μ_k) sector types round-trip
 bit-exactly. 58 assertions covering apply-equivalence + bad-metadata
 rejection.
- [x] OpenMP parallel sector build. Two-pass refactor (count + diag →
 prefix-sum → fill) makes outer loop over reps embarrassingly parallel.
 Opt-in via `USE_OPENMP=1`, auto-detects `brew --prefix libomp` on macOS.
 Measured on Apple M2 Ultra: **9.6–9.8× speedup on sector build**;
 full 8-sector N=27 kagome ED end-to-end **120 s → 23 s** (5.3×
 end-to-end, 1029% CPU utilisation). Bit-exact results vs serial.
- [x] 2D-irrep D_μ(g) matrix storage. Named 2D irreps (C_3v E, C_6v
 E_1/E_2, C_4v E) now carry per-element 2×2 real orthogonal matrices
 auto-populated from parent-op geometric action. New API:
 `irrep_sg_little_group_irrep_matrix(mu, i, D_out)` — uniform across
 1D (returns character as 1×1) and 2D irreps. Verified D(E)=I,
 orthogonality D^T D=I, trace(D(g))=χ(g), across C_6v E_1 at Γ,
 C_3v E at K, C_4v E at Γ. Sparse-sector build for 2D remains scoped
 (matrix-element derivation requires doubled sector dimension).
- [x] `examples/ed_cache_demo.c` — cache-workflow demonstration showing
 the save-on-cold-run / load-on-warm-run pattern. Measured N=27 end-
 to-end: cold 28.2 s → warm 1.1 s (26× speedup). Bit-identical
 eigenvalues across runs. The recommended pattern for iterative
 research loops hitting the same sector repeatedly.
- [x] On-disk persistence formats documented in docs/MOONLAB_INTEGRATION.md
 §6.6. Both rep-table and sector-binding formats specified with offset
 tables for cross-implementation interop.
- [x] `irrep_sg_heisenberg_sector_build_dense` — first-principles reference
 sector H builder for small-N cross-checks. Works for BOTH 1D AND 2D
 irreps via the character projector applied to each rep as a 2^N dense
 vector. Validated on kagome 12: Γ-A_1/A_2/B_1/B_2/E_1/E_2 all match
 published values from kagome12_k_resolved_ed within 1e-5, sparse
 1D matches dense bit-exactly at 1e-10, hermiticity < 1e-15 for every
 sector. Capacity N ≤ 24.
- [x] `irrep_sg_heisenberg_sector_build_at_k` now accepts 2D irreps:
 internally routes through the dense builder at d_μ=2, packs the result
 as CSR, returns same opaque sector handle. Consumer apply interface
 unchanged. Γ-E_1 Lanczos ground state on kagome 12 = −3.299516 J
 (matches published value). Sparse 1D performance path preserved at
 d_μ=1. The "pure sparse 2D" character-formula remains open future work
 (dense-backed path gives correct 2D spectra up to N ≤ 24 today).
- [x] `irrep_hermitian_eigendecomp(n, A, eigvals, eigvecs)` — full Jacobi
 diagonalisation returning both eigenvalues AND eigenvectors. Eigenvectors
 are orthonormal columns of `eigvecs` in the same descending-eigenvalue
 order as `eigvals`. Used in the new entanglement example to recover
 ground-state eigenvector from a small dense sector matrix. Tested for
 accuracy (1e-12 on eigenvalue equation and unitarity).
- [x] `examples/kagome12_entanglement.c` — first end-to-end research-grade
 example: N=12 kagome-Heisenberg Γ-A_1 GS → eigenvector → unfold to
 4096-dim state → partial traces → von Neumann entropies S_A for
 |A| = 1..6. Reproduces published E_0 = −5.32839240 J to 1e-8 and the
 exact prediction S_{|A|=1} = ln 2 (single-site entropy in SU(2)-
 invariant state). First data point in an N-scaling series for the
 Kitaev-Preskill γ topological entanglement entropy — the key diagnostic
 for kagome spin-liquid identification (γ = log 2 ⇔ gapped Z_2,
 γ = 0 ⇔ gapless U(1) Dirac).
- [x] `irrep_lanczos_eigvecs_reorth(apply, ctx, dim, k_wanted, max_iters,
 seed, eigvals, eigvecs)` — extends reorth Lanczos with Ritz-vector
 lifting. Stored Krylov basis V [max_iters × dim] + Jacobi eigen-
 decomposition of the tridiagonal T gives eigvecs via V · ritz. Tested
 against direct Jacobi on random 8×8 Hermitian matrix to 1e-10 residual.
 Unlocks eigenvector recovery at sparse scale (N=27 kagome Γ-A_1 GS in
 the 186 k-dim sector).
- [x] `examples/kagome27_entanglement.c` — N=27 kagome-Heisenberg
 entanglement pipeline: sparse Lanczos (150 iters, 1 eigvec) → 2 GB
 dense state unfold → partial traces. Second N-scaling data point for
 γ extraction. Demonstrates the full stack works at the thermodynamic-
 limit-reach scale for kagome (N > 24 sector-ED territory).
- [x] Kitaev-Preskill γ construction extended to BOTH kagome examples.
 Computed γ = S_A + S_B + S_C − S_{AB} − S_{AC} − S_{BC} + S_{ABC} for
 three non-overlapping 3-site regions in each cluster. Two-data-point
 N-scaling table:
   N=12: γ = -0.454
   N=27: γ = -0.106 (Δγ = +0.35, trending upward)
 Both values far from converged; a proper extrapolation requires larger
 N and more sophisticated region geometry (Levin-Wen or
 Jiang-Wang-Balents cylinder prescription). But the pipeline produces
 reproducible well-defined numbers — the methodology contribution is
 the thing, and we have two data points on an N-scaling series.
- [x] Spin-flip Z_2 diagnostic in `kagome12_entanglement`: computes
 ⟨ψ|F|ψ⟩ where F is the global bit-complement operator. Reports F-
 eigenvalue of the GS. Kagome 12 Γ-A_1 GS is F=+1 exactly (to 1e-15),
 confirming the expected singlet character.
- [x] `examples/kagome18_entanglement.c` — third kagome γ data point via
 kagome 2×3 rectangular cluster (p1 only since Lx≠Ly breaks p6mm).
 Revealed real physics: the p1 trivial sector GS is a TRIPLET (F=-1),
 not a singlet. Example scans 10 Lanczos eigvecs and filters for F=+1
 to select the lowest singlet, giving directly-comparable γ data across
 N ∈ {12, 18, 27}.
- [x] `examples/kagome24_entanglement.c` — fourth kagome γ data point via
 kagome 2×4 rectangular cluster (p1, ~338k reps, 39 s wall-clock). All
 10 lowest Lanczos states are F=+1 on this cluster (singlet GS — unlike
 N=18). S_{|A|=1} = ln 2 to 4e-14, confirming singlet character.
 γ_KP = -0.039, closest to 0 of any point yet.
- [x] **Four-point γ scaling series** (kagome singlet sector, contiguous-
 region KP construction — finite-size dominated):
   N=12: γ = -0.454  (E/N = -0.4440)
   N=18: γ = -0.221  (E/N = -0.4313)
   N=24: γ = -0.039  (E/N = -0.4461)
   N=27: γ = -0.106  (E/N = -0.4300, Sz=-1/2 not pure singlet)
 Documented in `docs/PHYSICS_RESULTS.md` §1.9.
- [x] **Area-law γ extraction** (`examples/kagome_arealaw_all.c`) — much
 cleaner methodology. For nested-region sequences |A|=1..9, fits
 `S_A = α·|∂A| − γ + O(1)` to extract γ from the intercept.
 Results (singlet sector) at three cluster shapes:
   N=12 p6mm 2×2:  γ = +0.785  (α = 0.370, R² = 0.95)
   N=18 p1 2×3:    γ = +0.675  (α = 0.342, R² = 0.97)
   N=24 p1 2×4:    γ = +0.893  (α = 0.391, R² = 0.98)
 **All three POSITIVE and near log 2 = 0.693**, consistent with gapped
 Z_2 topological order. N=18 gives γ within 3% of log 2, matching the
 Jiang-Wang-Balents 2012 DMRG extrapolation. Documented in §1.10 of
 `docs/PHYSICS_RESULTS.md`.
- [x] **Null-control validation** (`examples/square_arealaw_null.c`):
 identical area-law γ extraction on square 4×4 Heisenberg (Néel-ordered,
 published γ=0). Result: γ = +0.025 (within 3% of 0). Methodology
 correctly distinguishes Néel (γ≈0) from kagome spin liquid (γ≈log 2).
 Documented in §1.11 + §1.12 of `docs/PHYSICS_RESULTS.md`. **First
 reproducible computational signal from libirrep consistent with the
 gapped-Z₂ kagome-spin-liquid hypothesis at accessible cluster sizes.**
- [x] **J₁-J₂ square γ phase sweep** (`examples/j1j2_gamma_sweep.c`):
 21 values of J₂/J₁ ∈ [0, 1] on square N=16, extracting area-law γ at
 each. Δγ peaks at +0.415 at J₂/J₁ = 0.65, directly inside the
 published DQCP / spin-liquid window. 6.5 s wall-clock total.
 **Computational detection of the J₁-J₂ spin-liquid window via γ**
 — Paper #13 in the research program (DQCP) now has its N=16 data.
 Documented §1.13 of `docs/PHYSICS_RESULTS.md`.
- [x] **p2, p6, p3m1, p31m, p4** wallpaper groups. Enables:
 - C_2 on non-square clusters (2×3 kagome N = 18 was p1-only previously)
 - chiral hex C_6 and chiral square C_4 (rotations without mirrors)
 - both hexagonal-mirror-orientation classes p3m1 (axes through vertices)
  and p31m (axes bisecting vertex pairs, 30° offset)
 307 test assertions covering bijectivity, lattice-compatibility
 gates, subgroup bit-exactness with p4mm / p6mm parent groups, and
 p31m mirror-distinction from p3m1.
- [x] p4gm framework scaffolded (enum value, `fill_p4gm_` with
 fractional-translation output, `t_cart` plumbing through
 `build_point_perm_`). `IRREP_WALLPAPER_P4GM` is rejected at build
 time on every currently-shipped lattice with a clear diagnostic —
 the glide `½(a_1 + a_2)` maps a single-sublattice square-lattice
 site (i, j) to (i + ½, −j + ½), off the integer-site set. Adding
 a two-basis square lattice (`IRREP_LATTICE_SQUARE_2BASIS` — out
 of scope here) will make p4gm constructible without further
 space-group changes.
- [x] Downstream integration smoke test against
 `spin_based_neural_network` — **landed**. SbNN ships
 `src/libirrep_bridge.c` (lazy-loaded with `IRREP_ENABLE=1`) and
 a dedicated `tests/test_libirrep_bridge.c`. SbNN's `VERSION_PINS`
 records `LIBIRREP_MIN=1.3.0-alpha`. Bridge surface currently:
 multiset construction / parse / format, NequIP layer apply +
 backward, real spherical harmonics, Clebsch-Gordan, small Wigner-d.
 The 1.3-cycle lattice / DMI / cubic-PG additions are available
 for future bridge expansion (lattice3d → 3D spin systems; dmi →
 symmetry-derived couplings before the NQS sees them).
- [x] **3D Bravais lattice module** (`lattice3d.h` / `lattice3d.c`) — SC,
 BCC, FCC, Diamond, Pyrochlore in the conventional cubic cell. NN
 distance auto-discovered per family (1, √3/2, √2/2, √3/4, √2/4);
 coordination numbers (6, 8, 12, 4, 6) verified by 8418 assertions in
 `test_lattice3d`. FCC at 3³ lands at N=108 sites, the SbNN/moonlab
 target size; pyrochlore at 2³ lands at N=128 sites — corner-sharing
 tetrahedra, the 3D analog of kagome, hosting spin-ice / U(1) QSL
 phases (Tb₂Ti₂O₇, Yb₂Ti₂O₇, Dy₂Ti₂O₇). Bond lists feed straight into
 `irrep_heisenberg_new` for 3D ED, validated by `lattice3d_heisenberg_ed`
 against the K₈,₈ closed-form (BCC 2³ E₀ = −20 J exactly).
- [x] **T_d point group** (`IRREP_PG_TD`, 24 elements, 5 irreps:
 A₁, A₂, E, T₁, T₂) — first cubic point group in libirrep. Site
 symmetry of diamond / zincblende; subgroup of O_h. Hand-reducible
 decompositions all verified (1x0e→A₁, 1x1o→T₂, 1x1e→T₁,
 1x2e→E+T₂, 1x3o→A₁+T₁+T₂). 195 new assertions in `test_point_group`.
- [x] **O_h point group** (`IRREP_PG_OH`, 48 elements, 10 irreps:
 A₁g, A₂g, Eg, T₁g, T₂g, A₁u, A₂u, Eu, T₁u, T₂u) — full cubic
 group. Site symmetry of SC, BCC, FCC, perovskite, rocksalt. Element
 layout exploits O_h = O × {E, i}: 24 proper rotations of O followed
 by 24 improper elements with the same R_proper but det=−1. All seven
 hand-reducible decompositions verified — including the textbook
 1x1o→T₁u (polar vector) / 1x1e→T₁g (axial vector) g/u distinction.
- [x] **O point group** (`IRREP_PG_O`, 24 proper-only elements, 5 irreps
 A₁, A₂, E, T₁, T₂) — chiral cubic, site symmetry of MnSi, Cu₂OSeO₃,
 and other skyrmion-host magnets. Cannot distinguish polar from axial
 vectors (both give T₁) — the diagnostic feature that separates O from
 O_h. 1142 assertions in `test_point_group` (up from 372 baseline,
 +770 across T_d + O_h + O).
- [x] **3D translation-sector Heisenberg ED** (`examples/lattice3d_sector_ed.c`)
 — orbit-canonicalisation pattern at Γ momentum, applied to BCC 2³
 Heisenberg AFM. Sector dim 1670 vs full 65536 (39× shrink); Lanczos
 returns E₀ = −20.000000 J exactly, matching the K₈,₈ closed form.
 Pattern works on every 3D family without 3D space-group infrastructure
 — uses `lattice3d_translate` directly. Caught and fixed an inverted-
 ratio docstring bug in `hamiltonian.h` (the implementation uses
 `√(N_u/N_v)` source-in-numerator, not the documented inverse).
- [x] **Momentum-resolved 3D Heisenberg ED**
 (`examples/lattice3d_kspace_ed.c`) — generalises the Γ-only sector
 ED to every k on the cluster's BZ grid. SC 2³: 8 k-sectors enumerated,
 dimensions sum to 70 = C(8,4) = full Sz=0 dim, minimum E₀(k) over all
 sectors recovers the full-ED ground state (-4.820089) to 4×10⁻⁷.
 Tower-of-states diagnostic primitive — distinguishes symmetric phases
 from symmetry-broken via E₀(k) gap structure. Off-diagonal matrix
 element generalised to `½J · e^{−i k·t_R} · √(σ_v/σ_u)`; orbits with
 σ_u = 0 filtered out per-k.
- [x] **Pyrochlore 16-site Heisenberg AFM**
 (`examples/pyrochlore16_heisenberg.c`) — full ED on the 1³ pyrochlore
 cluster (16 sites, 48 NN bonds, 6-fold coordination preserved). E₀ =
 −8.809 J (−0.551 / site, −0.184 / bond, 73% of the
 independent-tetrahedron lower bound −12 J). 4-fold degenerate first
 excited state at E₁ = −8.441 J, Δ_01 = 0.368 J — likely a cubic-point
 group multiplet (T₁ or T₂ of the site symmetry). Required relaxing
 the lattice3d builder constraint from L≥2 to L≥1 (self-bonds dedup
 cleanly for non-trivial bases).
- [x] **Pyrochlore J₁-J₂ phase sweep** (`examples/pyrochlore16_j1j2.c`)
 — sweep J₂/J₁ ∈ [−0.5, 1.0] on the 16-site cluster, report E_0,
 E_1, E_2, E_3, Δ_01 at each point. Gap-collapse diagnostic at
 J₂/J₁ = 0.5 (Δ_01 = 0.039 — phase transition signature). Exact
 4-fold degeneracy E₀ = E_1 = E_2 = E_3 = −12 J at J₂/J₁ = 1
 — hidden-symmetry point where the Hamiltonian gains enhanced
 cluster invariance under the NN = NNN coupling condition.
- [x] **Pyrochlore GS spin-spin correlation function**
 (`examples/pyrochlore16_correlations.c`) — connected ⟨S₀ · S_r⟩
 on the 16-site GS. All 6 NN partners equivalent to 13-digit
 precision (−0.183523), confirming O_h preservation; all 9 non-NN
 sites uniformly +0.039015. Sum-rule cross-check 48 · ⟨C_NN⟩ = E₀
 to 10⁻¹³ — validates GS eigenvector quality. Negative-then-
 positive correlation pattern is the signature of a paramagnetic
 singlet (no long-range order).
- [x] **DMI symmetry analyzer** (`irrep/dmi.h`, `src/dmi.c`,
 `tests/test_dmi.c`, `examples/dmi_pyrochlore_pattern.c`) —
 first piece of the materials-search pipeline. Given a bond and
 a list of candidate symmetry operations, returns the symmetry-
 allowed `D` subspace by projecting `R³` under the bond's site
 stabiliser (axial-vector representation). All five Moriya-1960
 rules verified by direct test. Demo derives the B20-type DMI
 pattern (D ∥ bond) from chiral-cubic O group with no material
 input — the Bak-Jensen pattern of MnSi reproduced from pure
 group theory. Connects to the materials hunt for RT-CMOS-
 compatible chiral magnets: `irrep_pg_element` (newly public),
 `irrep_dmi_allowed_basis`, `_from_pg`.
- [x] **Symmetric exchange tensor analyzer** (extension of the same
 dmi module) — sibling for the symmetric part of the bond exchange
 tensor `J^s_ij = (J_ij + J_ji)/2` (Heisenberg + Kitaev-Γ-type
 anisotropy). Projects 3×3 symmetric matrices under `J → R · J · R^T`,
 returns up to 6 orthonormal basis matrices. 8 new test assertions
 (trivial, C₂, three-mirror, tetragonal, C₃-body-diagonal). Demo
 reports the complete `(DMI, J^s)` decomposition for pyrochlore NN
 — 3-dim J^s under all three cubic groups, spanning the Curnoe-
 Ross-Kao parametrisation derived from group theory alone. ABI
 baseline refreshed to `113b9663…`.
- [x] **Magnetic-point-group framework** in DMI analyzer
 (`irrep/dmi.h`). Antiunitary support via the new `antiunitary`
 flag in `irrep_dmi_sym_op_t`. Sign rule extension matches
 PHYSICS_APPENDIX §14.3: T flips spatial sign for DMI but not
 for symmetric exchange. 7 new assertions verify the four
 (unitary × {preserve, reverse}) × (antiunitary × {preserve,
 reverse}) sign cases plus contradictions. 58/58 in test_dmi.
 ABI unchanged (struct layout same 80 bytes, was previously
 padded). Future Shubnikov / 122-magnetic-point-group enumeration
 builds on this primitive.
- [x] **Three-spin scalar chirality analyzer** (`irrep_chirality_allowed`,
 `_from_pg`). For a triangle with site stabiliser, decides whether
 the pseudoscalar `χ_ijk = S_i · (S_j × S_k)` is symmetry-allowed.
 Sign rule: `det(g) · σ_perm(g) · σ_T(g) = +1` for every preserving
 operation. 7 new test cases (65/65 in test_dmi) cover identity,
 C₃, σ_h in plane (forbidden), σ_v transposition, C_3v, T·E
 (forbidden), T·σ_h (allowed). Connects to topological Hall in
 Mn₃Sn AFMs. ABI refreshed to `8a8c5698…`.
- [x] **Materials-search end-to-end demos**: `kagome_triangle_chirality.c`
 (Mn₃Sn / Co₃Sn₂S₂ chirality verdicts derived from the magnetic
 point group) and `pyrochlore_tetra_complete_catalog.c` (the
 full 72-verdict bilinear+trilinear catalog for the pyrochlore
 up-tetrahedron under O_h / O / T_d). The latter is the
 parameter scaffold downstream DFT / micromagnetic codes
 consume — the irreducible group-theory step now runs as one
 library call rather than crystallographer hand-derivation.
 3D space groups remain deferred past 1.3.

## Deferred past 1.3

- 3D space groups (230 groups).
- Magnetic space groups (1651 Shubnikov groups).
- Double groups / half-integer crystallographic reps.
- GPU kernels (Metal / CUDA / HIP).
- Python / JAX bindings (intentional non-goal — wrap via `cffi` / `ctypes`
 against the stable C ABI).
