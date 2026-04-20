# TODO

Milestone checklist. `[x]` = landed, `[ ]` = actual open work. Historical
milestones (M1–M9) are kept for context; the live open items are in M10
onwards plus the 1.3 section at the bottom.

## M1 — project skeleton (done)
- [x] LICENSE, NOTICE, README, CHANGELOG, CONTRIBUTING, etc.
- [x] VERSION, .gitignore, .clang-format, .editorconfig
- [x] All public headers with full API declarations (28 as of 1.3.0-alpha)
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
- [ ] Lebedev orders 9..41 — pending data import
- [x] Exactness tests + SH orthonormality via Gauss-Legendre tensor product

## M9 — equivariant layers (done)
- [x] Linear-on-irreps + analytic backward (matched-label channel mixing)
- [x] RMS norm per (term, copy) + analytic backward
- [x] Gate activation
- [x] Equivariance verified on linear; backward FD cross-checks on linear + RMS

## M10 — NEON hot paths (partial)
- [x] `irrep_cutoff_polynomial` batched NEON kernel (1.1, bit-exact)
- [x] `irrep_sph_harm_cart_all` NEON kernel (1.2)
- [ ] `irrep_wigner_d_matrix` NEON kernel
- [ ] `irrep_tp_apply_weighted_batch` NEON kernel
- [x] Bit-exact vs. scalar for kernels that exist

## Numerical stability replacements (open)
- [x] Wigner-d small: replaced Sakurai (3.8.33) direct-sum (unstable past
 j ≈ 20) with Edmonds Jacobi-polynomial form via DLMF §18.9.1 recurrence.
 Measured unitarity ≤ 1e-12 through j = 80; bounded only by IEEE-754
 lgamma overflow (j ≈ 170) past that.
- [x] Wigner 3j + Clebsch–Gordan: replaced Racah log-gamma single-sum
 (unstable past j ≈ 20, NaN past j ≈ 60 near triangle edge) with
 Schulten–Gordon backward three-term recurrence in j₁, normalised by the
 sum rule and sign-anchored at j_max. Machine precision through j ≈ 50.
- [ ] Miller two-directional iteration for 3j past j ≈ 80. Backward-only
 SG leaks the non-classical solution into the subdominant regime; Miller
 (forward + backward, match + renormalise) recovers machine precision
 through j ≈ 150+. Current regime is documented in tests/test_clebsch_
 gordan.c and the METHODS appendix.

## M11 — x86 SIMD hot paths (open)
- [ ] SSE4.2 + AVX2 variants of the above
- [ ] Runtime dispatch verified on both macOS x86_64 and Linux x86_64
- [ ] NEON-to-AVX2 bit-exactness acceptable tolerance documented

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
- [ ] Doxygen HTML generated and checked into release artifacts
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
- [ ] Full point-group irrep tables for p3m1, p31m, p4, p4gm, p2 — enables
 non-square kagome and triangular tilings under their symmetry-reduced
 point groups
- [ ] Downstream integration smoke test against
 `spin_based_neural_network` once that project consumes the 1.3 headers

## Deferred past 1.3

- 3D space groups (230 groups).
- Magnetic space groups (1651 Shubnikov groups).
- Double groups / half-integer crystallographic reps.
- GPU kernels (Metal / CUDA / HIP).
- Python / JAX bindings (intentional non-goal — wrap via `cffi` / `ctypes`
 against the stable C ABI).
