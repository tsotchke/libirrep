# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

## [1.0.0] - TBD

First public release. Targeted scope per the authoritative implementation plan.
