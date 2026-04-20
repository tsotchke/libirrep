# libirrep

> A pure-C11 representation-theory library, built to support
> work on open problems in frustrated quantum magnetism —
> with a tested SO(3)/SU(2)/O(3) / Clebsch-Gordan / Wigner-D /
> e3nn-style tensor-product core underneath.

## What this is for

libirrep exists to serve a concrete 1.3 cycle. The 1.3 cycle is
pinned to the **Kagome Heisenberg S = ½ ground-state-nature question** —
the flagship open problem in 2D frustrated magnetism, unresolved since
Yan, Huse & White, *Science* **332**, 1173 (2011). Two candidate
phases (gapped Z₂ spin liquid vs. gapless Dirac spin liquid) differ in
an extrapolated spin gap `Δ_S(N→∞) = 0.13 J` vs `0`. Settling it
requires symmetric neural-quantum-state (NQS) calculations at 108 sites
with genuine p6mm space-group projection — a combination no existing
public tool ships cleanly.

The library provides the group-theoretic, projection, and entanglement
machinery that NQS consumers need; it does not itself embed a
neural-network or an MCMC sampler (those live in the downstream
`spin_based_neural_network` project). Scope is pinned in
the 1.3 CHANGELOG;
physics results already demonstrated with the stack are catalogued in
[`docs/PHYSICS_RESULTS.md`](docs/PHYSICS_RESULTS.md).

## Results to date

End-to-end exact-diagonalisation on the 2×L_y kagome Heisenberg series,
produced entirely through libirrep's public C API (no LAPACK, no
external solver):

| N | E_0/N (J) | Δ_S (J) | runtime | literature cross-ref |
|----|---------------|----------|---------|-------------------------------|
| 12 | −0.45374 | 0.3827 | 0.4 s | Elser 1989, Lecheminant 1997 |
| 18 | −0.44707 | 0.2835 | 27 s | Waldtmann 1998, Läuchli 2011 |
| 24 | −0.44833 | 0.2639 | 2.5 min | 2×4 torus published ≈ −0.441 |

A three-point linear 1/N extrapolation of Δ_S across (12, 18, 24)
gives **≈ 0.132 J**, consistent with the Yan–Huse–White 2011 DMRG
value of 0.13 J for the gapped Z₂ picture. This is a sanity check,
not a scaling study: the three clusters have different 2×L_y
geometries, the extrapolation is linear-in-1/N, and reaching a regime
where a gapless Dirac alternative could be ruled out would require
larger sizes and more points. The 12-site cluster's ground state is
additionally decomposed across all six p6mm Γ-irreps (sector
dimensions summing to 1072 = dim V_Γ by Burnside); the ground-state
energy in the resolved B₁ block matches the power-iteration answer to
eight digits. Full details and per-result runtime / memory / code
pointers in [`docs/PHYSICS_RESULTS.md`](docs/PHYSICS_RESULTS.md).

## Library foundation

Underneath the physics substrate sits a tested representation-theory
core that stands on its own:

1. **Equivariant neural networks** for chemistry, materials, and physics —
 NequIP, MACE, Allegro, TENN, and the broader e3nn lineage. These networks
 factor their hidden state as direct sums of SO(3) / O(3) irreducible
 representations, and every non-trivial layer reduces to evaluating
 spherical harmonics on edge vectors, applying Wigner-D matrices to rotate
 features, or contracting features through Clebsch-Gordan-weighted tensor
 products.

2. **Quantum and classical spin systems** — magnetism simulations, density-
 matrix rotations, spin-coupling of multi-particle states, time-reversal
 symmetry on Kramers pairs, and group-theoretic classification of magnetic
 Hamiltonians.

3. **Disciplined rotation mathematics** for scientific computing and
 graphics — quaternions, Euler-ZYZ conversions, SO(3) exp/log stable near
 π, Shoemake-uniform sampling on S³, SLERP, and Karcher-Fréchet means.

All three rely on the same primitives; historically each re-implements them
with mutually incompatible sign and normalisation conventions (Condon-
Shortley vs. quantum-chemistry phase, ZYZ vs. ZXZ Euler, `{w, x, y, z}` vs.
`{x, y, z, w}` quaternions). libirrep picks one convention catalogue — the
physics one, matching Sakurai and Varshalovich and the e3nn reference
implementation — documents it exhaustively in
[`docs/PHYSICS_APPENDIX.md`](docs/PHYSICS_APPENDIX.md), cites every formula
to a primary source in [`docs/REFERENCES.md`](docs/REFERENCES.md), and
implements it to `1e-10` relative accuracy in double precision with a
stable ABI.

## Scope

Production features in libirrep 1.2:

**Core mathematical kernels.** Complex, real, and cartesian spherical
harmonics up to `l = 16`; Clebsch-Gordan coefficients for integer and
half-integer spin evaluated by the Racah single-sum formula in log-gamma
form; Wigner 3j, 6j, 9j, and Racah-W symbols; Wigner small-d and full-D
matrices stable past `j = 50`; all pairwise conversions between rotation
representations (matrix, quaternion, Euler ZYZ, axis-angle) round-trip to
`10⁻¹²`; Markley's branch-switching SO(3) logarithm; Karcher-Fréchet mean
on quaternions; Shoemake uniform sampling; Lebedev and Gauss-Legendre
quadrature.

**Irrep algebra.** An e3nn-compatible integer-l multiset type parsed
from strings like `"2x0e + 1x1o + 1x2e"`, with direct-sum, simplify, and
block-offset helpers. A parallel doubled-integer type `irrep_multiset_2j_t`
accepts half-integer labels (`"1x1/2o"`, `"2x3/2e"`) for mixed
integer-plus-half-integer spin content. Parity, time-reversal (integer
and half-integer), and the `T² = ±1` Kramers-degeneracy dichotomy
(detectable directly via `irrep_time_reversal_square_sign_2j`).

**Tensor products.** Path-indexed e3nn-style tensor products in both UUU
(matched-multiplicity) and UVW (independent-multiplicity) channel modes,
forward and backward, batched, with per-path scalar or full-tensor
weights. All valid `(l_a, l_b, l_c)` triples are supported including
odd-l-sum paths (cross-product-like), with the `i^{l_a + l_b − l_c}` phase
factor applied internally so real-basis output is guaranteed real.

**NequIP-style message-passing layer.** First-class `irrep_nequip_layer_t`
composing edge SH + Bessel radial basis + smooth cutoff + UVW tensor
product + neighbour aggregation; gradients through hidden features,
weights, and the edge geometry (`irrep_nequip_layer_apply_forces`);
spec-string constructor accepting e3nn-style layer descriptions
(`"2x0e + 1x1o -> 1x1o [sh=2, radial=8, r_cut=1.5]"`); per-path L2
regulariser for SR-style training.

**Point-group projection.** Character tables and projection operators
under C₄ᵥ, D₆, C₃ᵥ, and D₃ (Bradley-Cracknell 1972, cross-checked against
Altmann-Herzig 1994). Real-basis Wigner-D matrices are pre-computed at
table build time for `l ≤ 4`, yielding ~80× speedup on the projection hot
path. Direct-sum decomposition (`irrep_pg_reduce`) via the character
orthogonality formula.

**Equivariant-NN building blocks.** Linear-on-irreps layer with per-block
weights; RMS normalisation per irrep block; multiplicative gate activation.
All commute with the block-diagonal Wigner-D action by construction.

**Performance.** Runtime SIMD dispatch through a function-pointer table.
NEON kernels (2 lanes, aarch64) and AVX2+FMA kernels (4 lanes, x86_64)
ship today for the cartesian spherical harmonics batch, the polynomial
cutoff, and its derivative. Each vector kernel is bit-identical to the
scalar reference on representative inputs (FP_CONTRACT discipline kept
consistent between scalar and vector translation units). Hot-path
algorithmic wins that precede SIMD: `wigner_d_matrix` is 3.5× via
symmetry reduction, `tp_apply` 3.8× via sparse CG storage; both are
now candidates for a further 2–4× via SIMD in a 1.4 point release.

**Release engineering.** Stable ABI tracked by SHA-256 over the exported-
symbol set, baked into the binary (`irrep_abi_hash()`) and the installed
pkg-config variable (`pkg-config --variable=abi_hash libirrep`). SPDX 2.3
SBOM enumerating every source file with its license and checksum.
Continuous integration across macOS (arm64, x86_64) and Linux x86_64, with
sanitizer builds (ASan, UBSan), six libFuzzer targets running 60 s each,
ABI-drift gate, and regression-aware benchmark comparison.

## In development — 1.3 physics substrate

The 1.3 cycle is a 1.3 cycle: the candidate modules,
all already in-tree behind unreleased headers, are sized against the
Kagome Heisenberg S = ½ ground-state-nature problem (open since
Yan–Huse–White 2011, open problem). Scope and coordination with downstream
consumers are pinned in
the 1.3 CHANGELOG.

- **`irrep/lattice.h`** — square, triangular, honeycomb, kagome Bravais
 lattices with periodic boundary conditions; NN / NNN bond enumeration,
 primitive and reciprocal vectors, Brillouin-zone k-grids.
- **`irrep/space_group.h`** — `p1`, `p4mm`, `p6mm` wallpaper-group
 site-permutation tables. Full 6×6 kagome p6mm orbit sum (432 elements ×
 108 sites) benches at **17.6 µs** end-to-end on M2 Ultra — well under
 the 0.5 ms ceiling the 1.3 scope commits to.
- **`irrep/config_project.h`** — character-weighted reducer
 `P_μ ψ(σ) = (d_μ / |G|) Σ_g χ_μ*(g) ψ(g·σ)` for symmetric
 neural-quantum-state ansätze; orbit-enumeration helper.
- **`irrep/rdm.h`** — partial trace, cyclic-Jacobi Hermitian
 eigendecomposition (no LAPACK dependency), von Neumann / Rényi
 entropies, Kitaev–Preskill topological-entanglement-entropy
 `γ = S_A + S_B + S_C − S_{AB} − S_{BC} − S_{AC} + S_{ABC}`. The
 γ = ln 2 signal for Z₂ topological order is the primary physics
 observable driving the cycle.
- **`irrep/sym_group.h`** — `S_N`, permutation sign, hook-length
 irrep-dimension formula, antisymmetric (fermion) and totally-symmetric
 (boson) projectors on tensor-factored states.
- **`irrep/spin_project.h`** — total-`J` projection on `N` spin-½ sites
 via the Wigner-D character-weighted SU(2) integral. Restricts NQS
 ansätze to the singlet sector at half filling (Marshall sign).
- **`irrep/tensor_product.h` (half-integer path)** — complex-basis UVW
 tensor product on `irrep_multiset_2j_t` for spin-orbit-coupled
 equivariant layers and fermion spinors.

End-to-end demos exercising the substrate: `examples/kagome_a1_projection.c`
(full 108-site orbit + A₁ projection), `examples/heisenberg4_ed.c`
(4-site Heisenberg ED with J=0 singlet verification and entanglement
entropy at a 2-site cut).

## Deliberate non-goals

- **GPU kernels.** Consumers targeting Metal / CUDA / HIP / Vulkan compose
 their own on the stable C ABI. The internal function-pointer dispatch is
 ready for a GPU backend but we do not ship one.
- **Lie groups beyond SO(3) / SU(2) / O(3).** No SE(3), SE(2), SU(3), or
 higher-rank groups; these have their own natural libraries.
- **Autograd engine.** Every backward-capable primitive exposes a forward
 and a backward routine; composing them into a training loop is the
 caller's responsibility.
- **Python bindings.** Consumers wrap via the stable C ABI. A fast
 `cffi` or `ctypes` wrapper is mechanically straightforward.

## Quickstart

### Rotation math

```c
#include <stdio.h>
#include <irrep/wigner_d.h>

int main(void) {
 /* Wigner-D matrix for j = 1 at (α, β, γ) = (0.1, 0.2, 0.3). */
 double _Complex D[9];
 irrep_wigner_D_matrix(1, D, 0.1, 0.2, 0.3);
 for (int mp = 0; mp < 3; ++mp) {
 for (int m = 0; m < 3; ++m) {
 double _Complex z = D[mp * 3 + m];
 printf("D[%+d][%+d] = % .6f %+ .6fi\n", mp - 1, m - 1,
 (double)__real__ z, (double)__imag__ z);
 }
 }
 return 0;
}
```

### NequIP layer from a spec string (1.2+)

```c
#include <irrep/nequip.h>

int main(void) {
 irrep_nequip_layer_t *layer = irrep_nequip_layer_from_spec(
 "2x0e + 1x1o -> 1x1o [sh=2, radial=8, r_cut=1.5, cutoff=polynomial(6)]");
 if (!layer) return 1;

 int nw = irrep_nequip_layer_num_weights(layer);
 double w[nw]; /* from your optimiser */
 /* ... irrep_nequip_layer_apply(layer, w, n_nodes, n_edges, ...) */
 /* ... irrep_nequip_layer_apply_backward(...) + _apply_forces(...) */

 irrep_nequip_layer_free(layer);
 return 0;
}
```

### Point-group projection

```c
#include <irrep/point_group.h>
#include <irrep/multiset.h>

irrep_pg_table_t *c4v = irrep_pg_table_build(IRREP_PG_C4V);
irrep_multiset_t *spec = irrep_multiset_parse("2x0e + 1x1o + 1x2e");

/* Project a feature vector onto the totally-symmetric A₁ subspace. */
double in[10], out[10];
/* ... fill in ... */
irrep_pg_project(c4v, /*mu=A1*/ 0, spec, in, out);

irrep_multiset_free(spec);
irrep_pg_table_free(c4v);
```

### ABI drift check

```c
#include <irrep/version.h>
#include <string.h>

/* Built-against hash comes from pkg-config; pass via CFLAGS as
 * -DIRREP_ABI_HASH_EXPECTED="$(pkg-config --variable=abi_hash libirrep)"
 */
int libirrep_abi_ok(void) {
 return strcmp(IRREP_ABI_HASH_EXPECTED, irrep_abi_hash()) == 0;
}
```

Compile via:

```sh
cc -o demo demo.c $(pkg-config --cflags --libs libirrep)
```

or by directly pointing at the extracted release tarball:

```sh
cc -o demo demo.c \
 -I/path/to/libirrep/include \
 -L/path/to/libirrep/lib \
 -Wl,-rpath,/path/to/libirrep/lib \
 -llibirrep -lm
```

## Conventions (authoritative)

- **Angles** — radians throughout.
- **Rotations** — active, right-handed.
- **Euler** — ZYZ (Sakurai §3.3), `α ∈ [0, 2π)`, `β ∈ [0, π]`, `γ ∈ [0, 2π)`.
- **Phase** — Condon-Shortley `(−1)^m` applied once, in the associated
 Legendre polynomial.
- **Real SH** — e3nn sign convention, `Y_{1, +1}^{real} ∝ +x` at the equator.
- **Quaternions** — `{x, y, z, w}` layout with scalar `w` last (SIMD-
 ergonomic; matches Eigen / glTF 2.0).
- **Precision** — `double` by default; `_f32` wrappers for mixed pipelines.
- **Half-integer** — `_2j` suffix means every integer argument is
 doubled; spin-½ enters as `two_j = 1`.

Full derivations, edge-case handling, and citations in
[`docs/PHYSICS_APPENDIX.md`](docs/PHYSICS_APPENDIX.md). Migration
side-by-side with e3nn in
[`docs/MIGRATION_FROM_E3NN.md`](docs/MIGRATION_FROM_E3NN.md).

## Build and install

```
make # default: lib + test + examples + check-headers
make test # run 28 test suites
make bench # benchmark JSON under benchmarks/results/<utc>/
make asan # rebuild and test with -fsanitize=address
make ubsan # rebuild and test with -fsanitize=undefined
make docs # Doxygen HTML (if doxygen is installed)
make install # install to $(PREFIX) (default /usr/local)
make release # produce release/<VERSION>/ with tarball
```

CMake is also supported:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Documentation

| Document | Purpose |
| -------- | ------- |
| [`docs/PHYSICS_RESULTS.md`](docs/PHYSICS_RESULTS.md) | Concrete numerical results produced with libirrep — kagome ED at 12/18/24 sites, finite-size extrapolation, symmetry-sector decompositions, structure factors. |
| [`docs/METHODS.md`](docs/METHODS.md) | Methods-paper-style writeup of algorithmic choices (Racah log-gamma CG, Sakurai direct-sum Wigner-D, real-basis TP with `i^{l_a+l_b−l_c}` phase, orbit-rep basis builder, callback-based Lanczos, KP γ pipeline). |
| [`docs/DESIGN.md`](docs/DESIGN.md) | Architecture, module graph, feature → module mapping, numerical conventions, ABI policy, threading model. |
| [`docs/PHYSICS_APPENDIX.md`](docs/PHYSICS_APPENDIX.md) | Mathematical foundations, convention derivations, threshold rationales. |
| [`docs/REFERENCES.md`](docs/REFERENCES.md) | Annotated bibliography with DOIs and edition specificity. |
| [`docs/API.md`](docs/API.md) | Module-by-module API reference, linking to individual headers. |
| [`docs/MIGRATION_FROM_E3NN.md`](docs/MIGRATION_FROM_E3NN.md) | e3nn-to-libirrep mapping with sign conventions and caveats. |
| the 1.3 CHANGELOG | Research-agenda scope lock for the 1.3 cycle (the Kagome Heisenberg ground-state-nature problem, secondary targets, module table). |
| [`docs/tutorials/`](docs/tutorials/) | Seven tutorials: rotations, spherical harmonics, Clebsch-Gordan, Wigner-D, tensor products, equivariant NNs, and the 1.3 Kagome NQS substrate. |

Per-symbol Doxygen blocks are in the headers under `include/irrep/`;
`make docs` renders them to HTML.

## Citation

If you use libirrep in academic work, please cite:

```bibtex
@software{libirrep2026,
 author = {tsotchke},
 title = {libirrep: a high-performance SO(3)/SU(2)/O(3)
 irreducible-representation library for equivariant
 neural networks and spin systems},
 year = {2026},
 version = {1.3.0-alpha},
 url = {https://github.com/tsotchke/libirrep},
 note = {Pure C11; stable ABI; includes a physics substrate
 for symmetric neural-quantum-state work on the kagome
 Heisenberg S = 1/2 ground-state-nature problem}
}
```

## License

MIT. See [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE). Bundled Lebedev
quadrature tables (Lebedev & Laikov 1999) are public domain; attribution
in [`docs/REFERENCES.md`](docs/REFERENCES.md).
