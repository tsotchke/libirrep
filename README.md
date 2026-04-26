# libirrep

> **Status: pre-release — `v1.3.0-alpha`.** Public API guarded by an ABI
> hash and a regression gate. Numerical kernels at machine precision across
> every documented `j` regime. **31 public headers, 42 test suites, 0
> Doxygen warnings, ASan + UBSan clean.** Active downstream consumer:
> `spin_based_neural_network` (SbNN) pins `LIBIRREP_MIN=1.3.0-alpha`
> and links via `src/libirrep_bridge.c` against the multiset / NequIP /
> spherical-harmonic / CG / Wigner-d API surfaces.

A pure-C11 library of SO(3) / SU(2) / O(3) / SE(3) representation-theory
primitives — spherical harmonics, Wigner-D rotation matrices, Clebsch–
Gordan coupling, e3nn-style tensor products — with a stable ABI, no
runtime dependencies beyond `libc` + `libm`, and runtime-dispatched
NEON / AVX2 kernels on the hot paths.

### What's validated, what's experimental

| layer | maturity | notes |
|---|---|---|
| Rotation conversions, SH, CG, Wigner-D | **production** | machine precision through `j ≥ 80`, SIMD-accelerated, tested against analytical references |
| Multisets, tensor products, NequIP layer | **production** | e3nn parity; 2488× faster than e3nn on per-call CG (`scripts/bench_vs_e3nn.py`) |
| 2D lattices, wallpaper groups, sector ED | **production** | kagome 24-site ED matches Läuchli 2011 to 4 decimals |
| 3D Bravais lattices (5 families) | **validated** | BCC 2³ Heisenberg matches K₈,₈ closed-form `−n(n+2)/4` exactly |
| Cubic point groups (T_d, O_h, O) | **validated** | textbook crystal-field decompositions (eg + t_2g) reproduced |
| DMI / J^s / scalar-chirality analyzer | **validated against textbook** | reproduces Bak-Jensen B20 pattern, Curnoe-Ross-Kao pyrochlore, Mn₃Sn-via-T·σ_h |
| Magnetic point groups via `antiunitary` flag | **framework only** | sign rules verified per Moriya/Shubnikov; no pre-tabulated 122-group set yet |
| 3D space groups (P2_1 3, Fd-3m, etc.) | **deferred past 1.3** | only 2D wallpaper + 3D point groups in 1.3.0-alpha |
| Pyrochlore N≥32 sector ED | **out of reach** | requires the deferred 3D space-group machinery |
| Downstream consumer linkage | **live (SbNN)** | `spin_based_neural_network` pins `LIBIRREP_MIN=1.3.0-alpha`, integrates via `src/libirrep_bridge.c` and `tests/test_libirrep_bridge.c`. Currently consumes the multiset, NequIP layer, real spherical harmonics, CG, and small-d kernels; the 1.3 lattice / DMI / cubic-PG additions are available for future bridge expansion. |

## Quickstart

```c
/* demo.c — compile: cc demo.c $(pkg-config --cflags --libs libirrep) */
#include <stdio.h>
#include <irrep/wigner_d.h>

int main(void) {
    double _Complex D[9];                           /* j = 1, 3×3 */
    irrep_wigner_D_matrix(1, D, 0.1, 0.2, 0.3);     /* Euler ZYZ (α, β, γ) */
    for (int mp = 0; mp < 3; ++mp) {
        for (int m = 0; m < 3; ++m) {
            double _Complex z = D[mp * 3 + m];
            printf("D[%+d][%+d] = % .6f %+.6fi\n", mp - 1, m - 1,
                   __real__ z, __imag__ z);
        }
    }
    return 0;
}
```

```
D[-1][-1] =  0.955336 −0.041976i
D[-1][ 0] = −0.099003 −0.009903i
D[-1][+1] = −0.009950 −0.003324i
...
```

Longer quickstarts — NequIP layer from a spec string, point-group
projection, ABI-drift check — further down; see also the seven
[tutorials](docs/tutorials/).

## Who it's for

Three audiences share the same primitives. Pick the entry point that
matches your use case; the library is the same underneath.

### 1. Equivariant-NN inference engines

NequIP, MACE, Allegro, TENN, and the broader e3nn lineage all factor
hidden state as direct sums of SO(3) / O(3) irreps. Every non-trivial
layer reduces to (a) evaluating spherical harmonics on edge vectors,
(b) applying Wigner-D matrices to rotate features, and (c) contracting
features through Clebsch–Gordan-weighted tensor products. libirrep
ships all three as ABI-stable C primitives — with an e3nn-compatible
convention set (Condon–Shortley phase, ZYZ Euler) — so you can vendor
it into a C++/CUDA inference engine without pulling in Python.

Key entry points: `irrep_sph_harm_cart_all`, `irrep_wigner_d_matrix`,
`irrep_tp_build` / `_apply` / `_apply_weighted`, plus a first-class
`irrep_nequip_layer_t` composing edge SH + radial basis + smooth
cutoff + UVW tensor product + aggregation with full forward and
backward passes.

### 2. Physics / quantum-simulation code

Exact diagonalisation, symmetric-NQS ansätze, and entanglement
diagnostics on finite spin systems. The library carries explicit
lattice / space-group / reduced-density-matrix / Lanczos / Kitaev–
Preskill γ primitives, plus
[`irrep/hamiltonian.h`](include/irrep/hamiltonian.h) so the "Heisenberg
H-apply" that every ED example re-implements is one library call, not
forty lines of bit-twiddling per project. See
the "Physics application" section below.

### 3. Scientific and graphics code that needs disciplined rotations

Quaternions with `{x, y, z, w}` layout, Shoemake-uniform sampling on
S³, SLERP, Karcher–Fréchet mean, Markley's branch-switching SO(3)
logarithm stable near π, Euler-ZYZ with gimbal-lock guards,
axis-angle conversions round-trip to `10⁻¹²`.

## Measured performance

All numbers on Apple M2 Ultra. `ns/op` is wall-clock per call with
O2 optimization, `-ffp-contract=on`, FP_CONTRACT parity kept between
scalar and vector paths (bit-exact output). The Rosetta column is
`x86_64` cross-compiled with `-mavx2 -mfma` and run under Apple's
Rosetta 2 translator — a correctness oracle for the AVX2 path; native
x86_64 timings land when CI fires post-push.

| Kernel                                   | arm64 (NEON) | x86_64 (Rosetta, AVX2) |
|------------------------------------------|-------------:|-----------------------:|
| `sph_harm_cart_all` @ l=4                |   3.0 ns/op  |   6.8 ns/op            |
| `sph_harm_cart_all_batch` @ l=4, 4096 edges | 38.0 ns/edge | 175.2 ns/edge         |
| `tp_apply` on NequIP-like descriptor     | 138.3 ns/call| 229.6 ns/call          |
| `cg` @ j=1⊗1→1                           |  22.1 ns/call|  94.7 ns/call          |
| `wigner_d_matrix` @ j=8 (17×17)          |   4.8 µs/call|          — [1]         |
| `rbf_bessel` single call                 |   6.2 ns/op  |  12.5 ns/op            |

[1] Rosetta bench for `wigner_d_matrix` not yet in baseline; arm64
number reflects the 3.5× algorithmic win from Varshalovich §4.4.1
symmetry exploitation in [`src/wigner_d.c`](src/wigner_d.c).

### vs e3nn

[`scripts/bench_vs_e3nn.py`](scripts/bench_vs_e3nn.py) times the same
NequIP-shape tensor product on libirrep and on `e3nn` 0.6.0's
`FullyConnectedTensorProduct` (Python 3.12, torch 2.11.0 CPU
single-thread, float64). Same descriptor, same inputs, no batching
(single-sample per call — the regime equivariant-NN inference engines
spend most of their time in):

| Side     | ns / call |
|----------|----------:|
| libirrep |     138.3 |
| e3nn     | 339,140.7 |

libirrep is **~2488× faster on a per-call basis**. Numerical agreement
on unweighted CG output: cosine similarity 1.0000000000, max abs
error 1.3e-8 (e3nn's internal CG tables are float32-precision;
libirrep's are double). Run the comparison yourself:

```
python3 -m venv /tmp/venv && /tmp/venv/bin/pip install e3nn
/tmp/venv/bin/python scripts/bench_vs_e3nn.py
```

Baseline JSON at
[`benchmarks/results/baseline/`](benchmarks/results/baseline/);
`bench_vs_e3nn` runs in
[`benchmarks/results/e3nn_compare/`](benchmarks/results/e3nn_compare/);
regressions gated by `scripts/perf_compare.sh` at 5%.

## Precision regime

Machine precision across every `j` regime the test suite exercises.
The algorithms were rewritten in Round-2/3 of the 1.3.0-alpha audit
cycle; stability is measured, not aspirational.

| Kernel        | Test | `j=20` | `j=50` | `j=80`  | `j=120` |
|---------------|------|-------:|-------:|--------:|--------:|
| Wigner-d      | unitarity err        | 8e-15 | 6e-14 | 1e-13 | (not tested) |
| CG / 3j       | sum-rule drift       |   0   | 2e-16 | 4e-16 |   0    |

Wigner-d via Edmonds Jacobi-polynomial form (NIST DLMF §18.9.1
recurrence); CG via Schulten–Gordon / Miller two-directional iteration
(Luscombe–Luban 1998). Algorithm choices cited in
[`docs/METHODS.md`](docs/METHODS.md).

## Capabilities

**Core kernels.** Complex / real / cartesian spherical harmonics up
to `l = 16`. Clebsch–Gordan and Wigner 3j / 6j / 9j / Racah W for
integer and half-integer spin. Wigner small-d and full-D matrices with
symmetry-exploited matrix builder. All pairwise conversions between
rotation representations. Markley-stable SO(3) log. Shoemake and
Karcher–Fréchet quaternion utilities. Lebedev (3, 5, 7) and
Gauss–Legendre tensor-product quadrature.

**Irrep algebra.** e3nn-compatible multiset type parsed from strings
like `"2x0e + 1x1o + 1x2e"`; doubled-integer companion
`irrep_multiset_2j_t` for half-integer labels. Parity, time-reversal
(integer and half-integer), `T² = ±1` Kramers detection.

**Tensor products.** Path-indexed UUU and UVW e3nn-style tensor
products, forward and backward, batched, with per-path scalar or
full-tensor weights. Sparse CG inner loop (≈22 % density → 3.8×
speedup vs dense contraction). Half-integer (`_2j`) variant for
spin-orbit and spinor coupling.

**NequIP layer.** `irrep_nequip_layer_t` composes edge SH + Bessel
radial basis + smooth cutoff + UVW tensor product + aggregation.
Forward, backward through hidden features and weights, and through
edge geometry (`_apply_forces`). Spec-string constructor accepts
e3nn-style layer descriptions.

**Point-group and space-group projection.** Six 2D / 3D point groups:
C₄ᵥ / D₆ / C₃ᵥ / D₃ (lattice point groups for 2D wallpaper) plus
T_d / O_h / O (cubic point groups for diamond / zincblende /
perovskite / FCC / BCC / chiral skyrmion magnets). Character tables
with pre-computed real-basis Wigner-D matrices (~80× speedup on
projection). 2D wallpaper groups p1 / p4mm / p6mm as site-permutation
tables. Configuration-space character-weighted projection `P_μ ψ(σ)`
including non-Γ Bloch-momentum sectors.

**3D Bravais lattices.** SC / BCC / FCC / Diamond / Pyrochlore in the
conventional cubic cell, sibling of the 2D `lattice.h`. Auto-discovered
NN / NNN bond directions per family — coordination numbers (6, 8, 12,
4, 6) verified against textbook values. FCC at 3³ lands at exactly
N = 108 sites — the SbNN/moonlab target size; pyrochlore at 1³ gives a
16-site frustrated cluster. Bond lists feed straight into
`irrep_heisenberg_new` for 3D ED, full-space, Γ-momentum, or
k-resolved (see `examples/pyrochlore16_heisenberg.c`,
`examples/lattice3d_kspace_ed.c`).

**Hamiltonians, RDM, entanglement.** Spin-½ Heisenberg apply-op
(`irrep_heisenberg_t`) plugs directly into `irrep_lanczos_eigvals`.
Partial trace; cyclic-Jacobi Hermitian eigendecomposition (no LAPACK
dependency); von Neumann / Rényi entropies; Kitaev–Preskill γ formula.

**Release engineering.** Stable ABI tracked by SHA-256 over the
exported-symbol set, baked into the binary (`irrep_abi_hash()`) and
the installed pkg-config variable. Portable across Mach-O / ELF / PE.
SPDX 2.3 SBOM per release. Runtime SIMD dispatch. CI matrix: macOS
arm64 / x86_64, Linux x86_64 (gcc + clang), Linux arm64, Windows MinGW.
Sanitizer builds, fuzz-smoke, ABI-drift gate.

## Physics application

The library's physics substrate is sized against the **Kagome
Heisenberg S = ½ ground-state-nature problem** — open since
Yan, Huse & White, *Science* **332**, 1173 (2011) — which asks whether
the spin gap `Δ_S(N→∞)` is `≈ 0.13 J` (gapped Z₂ spin liquid) or `0`
(gapless Dirac spin liquid). The answer needs symmetric
neural-quantum-state calculations at 108 sites with genuine p6mm
space-group projection.

### ED validation (12 / 18 / 24 sites)

End-to-end ED through the public C API (no LAPACK, no external
solver):

| N | E_0/N (J) | Δ_S (J) | runtime | cross-ref |
|---|-----------|---------|---------|-----------|
| 12 | −0.45374 | 0.3827 | 0.4 s | Elser 1989, Lecheminant 1997 |
| 18 | −0.44707 | 0.2835 | 27 s | Waldtmann 1998, Läuchli 2011 |
| 24 | −0.44833 | 0.2639 | 2.5 min | 2×4 torus published ≈ −0.441 |

A linear 1/N extrapolation of Δ_S gives `≈ 0.132 J`, consistent with
the YHW 2011 DMRG value — a sanity check, not a scaling study.
Reproducibility catalogued in
[`examples/EXPECTED_OUTPUT.md`](examples/EXPECTED_OUTPUT.md).

### 108-site substrate

[`examples/kagome_a1_projection.c`](examples/kagome_a1_projection.c)
builds the full 6×6-kagome p6mm group (432 permutation elements on
108 sites), enumerates the 432-image orbit of a trial classical-spin
configuration, and character-A₁-projects a trial amplitude in
~120 µs. The ED at this size is intractable (dim ≈ 3.25 × 10³²); the
substrate is consumed by downstream NQS / VMC code. The 1.3 primitives
that enable this — `lattice.h`, `space_group.h`, `config_project.h`,
`rdm.h`, `sym_group.h`, `spin_project.h`, `hamiltonian.h` — are all
shipped and tested.

### 3D extension

The 1.3-alpha cycle also adds a 3D substrate via `lattice3d.h` and the
cubic point groups T_d / O_h / O. The natural 3D analog of the kagome
problem is the **pyrochlore lattice** (corner-sharing tetrahedra,
3D analog of corner-sharing triangles, hosts spin-ice and quantum
spin-liquid candidates). `examples/pyrochlore16_heisenberg.c` runs ED
on the 16-site pyrochlore Heisenberg AFM:

| observable | value | interpretation |
|---|---|---|
| E₀ | −8.809 J | per site −0.551, per bond −0.184 |
| E₁ | −8.441 J | 4-fold degenerate (cubic-multiplet structure) |
| Δ_01 | 0.368 J | gap to first excited |
| frustration measure | 73 % of −12 J | independent-tetrahedron lower bound |

The same pipeline drives `lattice3d_sector_ed.c` (Γ-momentum reduction
on BCC 2³ — recovers the K_{8,8} closed form −20 J exactly with 39× dim
shrink) and `lattice3d_kspace_ed.c` (full BZ momentum scan — minimum
E₀(k) recovers the full-space ground state to 4×10⁻⁷). All without
3D-space-group infrastructure — pure translation symmetry on the 3D
lattice via `irrep_lattice3d_translate`.

Details in [`docs/PHYSICS_RESULTS.md`](docs/PHYSICS_RESULTS.md).

## More quickstart snippets

### NequIP layer from a spec string

```c
#include <irrep/nequip.h>

int main(void) {
    irrep_nequip_layer_t *layer = irrep_nequip_layer_from_spec(
        "2x0e + 1x1o -> 1x1o [sh=2, radial=8, r_cut=1.5, cutoff=polynomial(6)]");
    if (!layer) return 1;

    int nw = irrep_nequip_layer_num_weights(layer);
    double w[nw];                                   /* from your optimiser */
    /* irrep_nequip_layer_apply(layer, w, n_nodes, n_edges, ...) */
    /* irrep_nequip_layer_apply_backward(...) + _apply_forces(...)  */

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

double in[10], out[10];
/* ... fill in ... */
irrep_pg_project(c4v, /*mu=A1*/ 0, spec, in, out);

irrep_multiset_free(spec);
irrep_pg_table_free(c4v);
```

### Heisenberg ED through Lanczos

```c
#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>

irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
int N  = irrep_lattice_num_sites(L);        /* 12 */
int nb = irrep_lattice_num_bonds_nn(L);
int bi[216], bj[216];
irrep_lattice_fill_bonds_nn(L, bi, bj);

irrep_heisenberg_t *H = irrep_heisenberg_new(N, nb, bi, bj, 1.0);
long long dim = irrep_heisenberg_dim(H);
double _Complex seed[dim];                  /* fill in an Sz=0 pattern */
double eig[2];
irrep_lanczos_eigvals(irrep_heisenberg_apply, H, dim, 2, 80, seed, eig);

irrep_heisenberg_free(H);
irrep_lattice_free(L);
```

### ABI drift check at consumer compile time

```c
#include <irrep/version.h>
#include <string.h>

/* CFLAGS: -DIRREP_ABI_HASH_EXPECTED="$(pkg-config --variable=abi_hash libirrep)" */
int libirrep_abi_ok(void) {
    return strcmp(IRREP_ABI_HASH_EXPECTED, irrep_abi_hash()) == 0;
}
```

### Compile

```sh
# via pkg-config (after `make install`)
cc -o demo demo.c $(pkg-config --cflags --libs libirrep)

# or against an extracted release tarball
cc -o demo demo.c \
   -I/path/to/libirrep/include \
   -L/path/to/libirrep/lib \
   -Wl,-rpath,/path/to/libirrep/lib \
   -llibirrep -lm
```

## Deliberate non-goals

- **GPU kernels.** Compose against the stable C ABI. Dispatch-table
  ready but no Metal / CUDA / HIP kernels shipped.
- **Lie groups beyond SO(3) / SU(2) / O(3).** No SE(3), SU(3), higher-
  rank; they have their own natural libraries.
- **Autograd engine.** Every backward-capable primitive exposes forward
  + backward; composing them into a training loop is the caller's job.
- **Python bindings.** Consumers wrap via the stable C ABI. A fast
  `cffi` / `ctypes` shim is mechanically straightforward.

## Conventions (authoritative)

- **Angles** — radians throughout
- **Rotations** — active, right-handed
- **Euler** — ZYZ (Sakurai / physics convention)
- **Phase** — Condon–Shortley
- **Quaternions** — `{x, y, z, w}`, unit-norm, Shoemake with `w ≥ 0`
- **Complex** — `double _Complex` (`<complex.h>`)
- **Half-integer spin** — doubled-integer API, `_2j` suffix

Derivations and primary-source citations in
[`docs/PHYSICS_APPENDIX.md`](docs/PHYSICS_APPENDIX.md).

## Build and install

```
make             # lib + test + examples + check-headers + check-abi
make test        # run 42 test suites
make bench       # benchmark JSON under benchmarks/results/<utc>/
make asan        # rebuild and test with -fsanitize=address
make ubsan       # rebuild and test with -fsanitize=undefined
make coverage    # LLVM source-based coverage report
make fuzz-run    # non-libFuzzer driver, 1M iters per harness
make docs        # Doxygen HTML (if doxygen is installed)
make install     # install to $(PREFIX) (default /usr/local)
make release     # produce release/<VERSION>/ with per-triple tarball
```

CMake also supported:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Reproducible Linux environment via Docker (mirrors the CI runner):

```
docker build -t libirrep-dev .
docker run --rm -v "$(pwd):/w" -w /w libirrep-dev make all

# cross-arch from Apple Silicon to validate the x86_64 AVX2 path:
docker build --platform linux/amd64 -t libirrep-dev-x64 .
docker run --rm --platform linux/amd64 -v "$(pwd):/w" -w /w libirrep-dev-x64 make all
```

## Documentation

| Document | Purpose |
| -------- | ------- |
| [`docs/API.md`](docs/API.md) | Module-by-module API index |
| [`docs/DESIGN.md`](docs/DESIGN.md) | Architecture, module graph, ABI policy, threading |
| [`docs/METHODS.md`](docs/METHODS.md) | Algorithmic choices (Jacobi Wigner-d, Miller 3j, sparse CG TP, Lanczos, KP γ) |
| [`docs/PHYSICS_APPENDIX.md`](docs/PHYSICS_APPENDIX.md) | Mathematical foundations, convention derivations |
| [`docs/REFERENCES.md`](docs/REFERENCES.md) | Annotated bibliography (DOIs + edition specificity) |
| [`docs/MIGRATION_FROM_E3NN.md`](docs/MIGRATION_FROM_E3NN.md) | e3nn-to-libirrep mapping |
| [`docs/PHYSICS_RESULTS.md`](docs/PHYSICS_RESULTS.md) | Kagome ED numerical results |
| [`docs/tutorials/`](docs/tutorials/) | Eight walk-throughs (rotations → spherical harmonics → CG → Wigner-D → tensor products → equivariant NNs → kagome NQS substrate → materials-search pipeline) |
| [`examples/EXPECTED_OUTPUT.md`](examples/EXPECTED_OUTPUT.md) | Reproducibility catalogue |
| [`CHANGELOG.md`](CHANGELOG.md) | Per-release changes |

Per-symbol Doxygen blocks in the headers; `make docs` renders to HTML.

## Citation

```bibtex
@software{libirrep2026,
  author  = {tsotchke},
  title   = {libirrep: SO(3)/SU(2)/O(3) representation-theory primitives
             for equivariant neural networks and spin systems},
  year    = {2026},
  version = {1.3.0-alpha},
  url     = {https://github.com/tsotchke/libirrep},
  note    = {Pure C11; stable ABI; Jacobi-polynomial Wigner-D,
             Schulten–Gordon Miller-iterated 3j, sparse-CG tensor
             product. Also ships a physics substrate for
             symmetric-NQS work on the kagome Heisenberg S = 1/2
             ground-state-nature problem}
}
```

A GitHub-native `CITATION.cff` is also provided.

## License

MIT. See [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE). Bundled Lebedev
quadrature tables (Lebedev & Laikov 1999) are public domain;
attribution in `docs/REFERENCES.md`.
