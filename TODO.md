# TODO

Milestones per the authoritative implementation plan.

## M1 — project skeleton (in progress)
- [x] LICENSE, NOTICE, README, CHANGELOG, CONTRIBUTING, etc.
- [x] VERSION, .gitignore, .clang-format, .editorconfig
- [ ] All public headers with full API declarations
- [ ] Stub .c sources that link
- [ ] Makefile + CMakeLists.txt
- [ ] TAP-style test harness and one stub test per module
- [ ] Scripts, benchmarks, examples, docs scaffolding
- [ ] `git init` + initial commit
- [ ] `make`, `make test`, `make examples` all succeed

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
- [ ] Orthonormality test via Lebedev — deferred to M8 when tables ship
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
- [ ] Channel modes beyond "uuu" (uvw, uvu) — future
- [ ] Sparse CG-block storage for compile-time path lists — future

## M8 — recoupling + radial + quadrature (done)
- [x] Wigner 6j, 9j, Racah W (log-gamma)
- [x] Bessel and Gaussian RBFs
- [x] Cosine / polynomial cutoffs + analytic derivatives
- [x] Gauss-Legendre via Newton iteration (any n)
- [x] Lebedev orders 3, 5, 7 — remaining orders (9..41) pending data import
- [x] Exactness tests + SH orthonormality via Gauss-Legendre tensor product

## M9 — equivariant layers (done)
- [x] Linear-on-irreps + analytic backward (matched-label channel mixing)
- [x] RMS norm per (term, copy) + analytic backward
- [x] Gate activation
- [x] Equivariance verified on linear; backward FD cross-checks on linear + RMS

## M10 — NEON hot paths
- [ ] `irrep_sph_harm_cart_all` NEON kernel
- [ ] `irrep_wigner_d_matrix` NEON kernel
- [ ] `irrep_tp_apply_weighted_batch` NEON kernel
- [ ] Bit-exact vs. scalar for representative inputs

## M11 — x86 SIMD hot paths
- [ ] SSE4.2 + AVX2 variants of the above
- [ ] Runtime dispatch verified on both macOS x86_64 and Linux x86_64

## M12 — benchmarks + fuzz + sanitizer CI
- [ ] All benchmarks produce JSON
- [ ] libFuzzer targets for CG, multiset parser, TP apply
- [ ] GitHub Actions matrix: 5 triples × {gcc, clang} × {Debug, Release} × sanitizers
- [ ] `perf_compare.sh` baseline committed

## M13 — release build system
- [ ] `scripts/build_release.sh 1.0.0` produces `release/1.0.0/`
- [ ] Per-triple tarballs + checksums
- [ ] `ABI_HASH` emitted and verified
- [ ] `libirrep.pc` and `libirrepConfig.cmake` installed
- [ ] SBOM (SPDX 2.3 JSON)

## M14 — examples, tutorials, docs finalization
- [x] `spin_half_rotation` (Berry phase π)
- [x] `equivariant_mlp` (toy E(3) MLP end-to-end, equivariance verified)
- [x] `nequip_message` (one NequIP-style TP convolution)
- [ ] All six tutorials committed
- [ ] Doxygen HTML generated
- [ ] README polished
- [ ] Consumer smoke test: 20-line C file evaluating `Y_1^0(π/4)` passes
