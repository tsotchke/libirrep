# libirrep — Integration Report for Moonlab

**Source**: `libirrep` 1.3.0-alpha (commit-of-record ABI: `f58cb703…9469f05`)
**Audience**: Moonlab engineers (tensor-network + state-vector simulator at `~/Desktop/quantum_simulator`)
**Purpose**: specify the integration surface that lets Moonlab use libirrep's symmetry
infrastructure to target the spin-½ kagome-Heisenberg (KH) problem in joint work with SbNN.

This document is pragmatic, not aspirational. Every API named here exists and is tested.
Every gap called out is a gap — we flag what we will add, and what we need Moonlab to supply.

---

## 1. Executive summary

libirrep is a pure-C11 library that provides the **symmetry spine** for equivariant
numerical physics: SO(3)/SU(2) representation theory, 2D wallpaper-group projection, 
Clebsch–Gordan / Wigner-D machinery, and — new in 1.3 — the full (k, μ_k) sector
decomposition pipeline for 2D crystallographic lattices. It does **not** contain a
state-vector, tensor-network, or eigensolver framework; it is the layer that tells
those frameworks how to block-diagonalise their Hilbert spaces.

For KH in particular, libirrep supplies:

- the point-group and wallpaper-group tables on kagome (p6mm + lower-symmetry
  subgroups for finite clusters that break p6mm)
- the little-group at every Bloch momentum `k` (C_6v at Γ, C_3v at K, C_2v at M, …)
- named-irrep character rows for C_6v, C_3v, C_2v, C_4v (A_1, A_2, B_1, B_2, E, E_1, E_2)
- composite projectors `P_{k, μ_k}` that act on full computational-basis vectors
- RDM + entanglement-entropy diagnostics (von Neumann, Rényi-α,
  Kitaev-Preskill topological γ)
- matrix-free Lanczos (with and without full reorthogonalisation) and stable
  Heisenberg / J₁-J₂ / XY Hamiltonian apply-ops

What Moonlab brings: DMRG (2-site + subspace expansion), TDVP, KPM, MPS/MPO
infrastructure, state-vector ED to 32 qubits, kagome-compatible lattice-2D.

The joint KH attack requires two bridges that do not yet exist:

1. **libirrep → Moonlab MPO emission** of the composite projector `P_{k,μ_k}` —
   concretely: exporting the group-element / site-permutation / weight triples
   as a structure Moonlab's MPO builder can consume.
2. **libirrep orbit-representative bitstrings** — for Lanczos-in-sector on
   representatives rather than full-Hilbert vectors, so the torus ED that
   complements DMRG cylinder results can reach N = 27, 30, 36.

Both are scoped and tractable; §6 sketches them concretely.

---

## 2. What libirrep is (and isn't)

### Is

- **A pure-C11 library** with a stable C ABI, versioned SONAME, and a
  hash-pinned public symbol set. No C++, no runtime, no hidden state
  (one initialised-once SIMD dispatch table, acquired via C11 atomics).
- **Thread-safe on disjoint inputs**. Every public function is pure on its
  arguments except for the SIMD dispatch lazy-init. Lanczos, projection,
  and RDM can all be called concurrently from multiple threads with
  independent `(ctx, in, out)` tuples.
- **Bit-exact vectorised** on NEON (arm64) and AVX2+FMA (x86_64). Scalar
  reference is the same code path at `-ffp-contract=on`; tests assert
  bit-equality.
- **Numerically stable** where traditional closed forms fail:
  Wigner-d via Edmonds' Jacobi-polynomial DLMF §18.9.1 recurrence
  (unitarity ≤ 1e-12 through j = 80), CG via Miller-iterated
  Schulten–Gordon (machine precision to j ≳ 200).

### Is not

- **Not a tensor-network library.** No MPS, no MPO, no DMRG. It produces
  projectors; consumers apply them.
- **Not a state-vector simulator.** `irrep_lanczos_eigvals` takes a
  matrix-free apply callback over `double _Complex` vectors of length
  `2^N` — the caller owns the dimensionality, and at N ≳ 26 the caller
  must provide sparse apply-op machinery or live with dense scaling.
- **Not Python-bindable in-tree.** Consumers wrap via `cffi` / `ctypes`
  against the stable C ABI. Intentional non-goal.

---

## 3. Module map

31 public headers in `include/irrep/`, grouped by role. LOC counts for sizing.

### Representation-theory core
| Header | LOC | Purpose |
|---|---|---|
| `so3.h` | 136 | SO(3): quat/axis-angle/matrix conversions, SLERP, Fréchet mean |
| `su2.h` | 58 | SU(2): Pauli algebra, double cover of SO(3), exp/log |
| `types.h` | 108 | Irrep labels, status codes, library-wide limits |
| `parity.h` | 31 | O(3) parity on irreps |
| `time_reversal.h` | 52 | T operator, Kramers |

### Angular functions
| Header | LOC | Purpose |
|---|---|---|
| `spherical_harmonics.h` | 80 | `Y_l^m`: real, complex, cartesian forms; gradients |
| `solid_harmonics.h` | 45 | Real regular solid harmonics (homogeneous polys) |
| `wigner_d.h` | 83 | `d^j_{mm'}(β)`, full D, batched over β (NEON+AVX2) |
| `clebsch_gordan.h` | 70 | CG, 3j, dense tables; Miller-iterated for stability |
| `recoupling.h` | 40 | 6j, 9j, Racah-W |
| `quadrature.h` | 82 | Gauss-Legendre; Lebedev 3–41 via runtime registry |

### Crystallographic / spatial symmetry
| Header | LOC | Purpose |
|---|---|---|
| `lattice.h` | 158 | 2D Bravais (square, triangular, honeycomb, kagome), PBC, NN/NNN bonds |
| `lattice3d.h` | 184 | 3D Bravais (SC, BCC, FCC, Diamond, Pyrochlore), PBC, NN/NNN bonds |
| `space_group.h` | 153 | 2D wallpaper groups (p1, p2, p4, p4mm, p6, p6mm, p3m1, p31m; p4gm scaffolded) |
| `config_project.h` | 367 | Space-group projection; **little-group at k**; composite (k, μ_k) projectors |
| `point_group.h` | 123 | Discrete subgroups of O(3): C₄ᵥ / D₆ / C₃ᵥ / D₃ + cubic T_d / O_h / O |
| `dmi.h` | 148 | Bond-exchange-tensor symmetry analyzers — DMI vector + symmetric J^s |

### Multiset / equivariance
| Header | LOC | Purpose |
|---|---|---|
| `multiset.h` | 71 | Integer-spin irrep multisets (`"1x0e + 2x1o + 1x2e"`) |
| `multiset_2j.h` | 97 | Half-integer (doubled) irrep multisets |
| `tensor_product.h` | 253 | e3nn-style TP paths; forward/backward; dim-first batch variant (NEON+AVX2) |
| `equivariant_layers.h` | 63 | Linear mixing, RMS norm, gate activation |
| `radial.h` | 79 | Bessel/Gaussian RBFs, smooth cutoffs with analytic derivatives |
| `nequip.h` | 178 | NequIP-style E(3)-equivariant message-passing layer |

### Condensed-matter / MBQM
| Header | LOC | Purpose |
|---|---|---|
| `hamiltonian.h` | 115 | On-the-fly Heisenberg / J₁-J₂ / XY apply-ops |
| `rdm.h` | 244 | Partial trace, Hermitian eigen (Jacobi), VN/Rényi entropy, **Lanczos (±reorth)**, Kitaev-Preskill γ |
| `spin_project.h` | 90 | Total-J projection on N-site spin-½ systems |
| `sym_group.h` | 100 | S_N irreps, Young tableaux, hook-length, A/S projectors |

### ABI / infrastructure
| Header | LOC | Purpose |
|---|---|---|
| `export.h`, `irrep.h`, `version.h`, `simd.h` | <200 total | Symbol visibility, umbrella include, version/ABI hash |

---

## 4. ABI contract

### ABI hash

Every release bakes a SHA-256 of the sorted set of public `irrep_*` symbol
names (SIMD-suffixed kernels excluded — platform-specific, reached only
through internal dispatch pointers). Callers retrieve the compiled-in hash
at runtime:

```c
#include <irrep/version.h>
const char *hash = irrep_abi_hash();    /* baked at build time */
const char *ver  = irrep_version();     /* "1.3.0-alpha" */
```

**Current baseline** (1.3.0-alpha, post 3D-lattice + cubic-point-group + DMI-symmetry-analyzer additions): `113b96630e530184f37b5671955e402d299bf7607944b20dd36d0f878c6307b9`. Successive 1.3.0-alpha development baselines are documented in `release/BASELINE_ABI_HASH` and the [Unreleased] section of `CHANGELOG.md`; consumers tracking 1.3.0-alpha should re-pin on each refresh.

Moonlab should **pin the ABI hash** in a vendored `IRREP_EXPECTED_ABI_HASH`
constant and assert at program start. Any divergence is either an
accidental breaking change (our bug — file an issue) or an intentional
major-version bump (you must re-vet).

### Symbol visibility

- Built with `-fvisibility=hidden`; public API tagged `IRREP_API`.
- GCC/clang: `__attribute__((visibility("default")))`.
- Windows: `__declspec(dllexport)` when building, `dllimport` for consumers
  (auto-driven by `IRREP_BUILDING_DLL` macro).

### SONAME / dylib

| Platform | SONAME / install_name |
|---|---|
| Linux | `liblibirrep.so.1` |
| macOS | `liblibirrep.1.dylib`, `-install_name @rpath/liblibirrep.1.dylib` |
| Windows | `liblibirrep.dll` (+ `liblibirrep.dll.a` import lib, no SONAME versioning) |

### Versioning policy

- `1.3.x` minor bumps: additive ABI only. Existing symbols stable; new
  symbols may be added. Hash changes; no re-vet required beyond picking
  up the new baseline.
- `2.0.0`: first breaking change (no such change is planned near-term).

---

## 5. Build, consumption, supported triples

### Build

```bash
# Canonical (Makefile)
make all        # lib + test + examples + check-headers + check-abi
make test       # 37 suites, ~7s on M2 Ultra
make asan       # clean
make ubsan      # clean
make docs       # Doxygen HTML (0 warnings)
make release    # per-triple tarballs in release/<version>/

# CMake (mirrors Makefile; recommended for Moonlab's CMake tree)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Consume (Moonlab's vantage)

**Option A — pkg-config** (preferred on Linux/macOS):

```bash
pkg-config --cflags libirrep    # -I/path/to/include
pkg-config --libs   libirrep    # -L... -llibirrep -Wl,-rpath,...
```

**Option B — CMake** (preferred cross-platform):

```cmake
find_package(libirrep 1.3.0 REQUIRED)
target_link_libraries(moonlab_kagome PRIVATE libirrep::libirrep)
```

**Option C — vendored tarball**: per-triple prebuilt tarballs land on
GitHub Releases under `v1.3.0-alpha`. Moonlab's `vendor/libirrep/`
directory carries a pinned version file; CI fetches the matching tarball.

### Supported triples (CI matrix, all green)

- `macos-arm64-clang` (M-series native, Apple clang)
- `macos-x86_64-clang` (Intel, Apple clang, AVX2)
- `linux-x86_64-gcc`, `linux-x86_64-clang` (AVX2)
- `linux-arm64-gcc` (NEON)
- `windows-x86_64-mingw` (MSYS2 / MinGW64, AVX2; no NEON path needed)

---

## 6. The KH integration surface

This section is what Moonlab engineers spend most of their time on.

### 6.1 The sector-projection pipeline (exists)

Given a kagome cluster of size `Lx × Ly` sites, target ground state at
Bloch momentum `k = (kx/Lx) b1 + (ky/Ly) b2` inside little-group irrep
`μ_k` (e.g. C_6v A_1 at Γ, C_3v E at K, C_2v B_1 at M):

```c
/* 1. Build the lattice + space group. */
irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, Lx, Ly);
irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
int order_G = irrep_space_group_order(G);          /* = 12 · Lx · Ly */
int n_sites = irrep_space_group_num_sites(G);      /* = 3 · Lx · Ly  */

/* 2. The little group at k (automatically detects C_6v / C_3v / C_2v). */
irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, kx, ky);
int n_lg_point = irrep_sg_little_group_point_order(lg);

/* 3. The irrep μ_k. Use a builtin name, or supply characters explicitly. */
irrep_sg_little_group_irrep_t *mu_k =
    irrep_sg_little_group_irrep_named(lg, IRREP_LG_IRREP_B1);
/* Or: irrep_sg_little_group_irrep_new(lg, chi_of_slot, dim); */

/* 4. Project a caller-supplied "orbit vector" ψ(g), g ∈ G.
 *    Caller enumerates the orbit (or evaluates ψ lazily on each orbit
 *    member) and hands the order_G amplitudes to the library. */
double _Complex amplitude = irrep_sg_project_at_k(lg, mu_k, psi_of_g);

/* 5. Or: build a full orthonormal basis of the (k, μ_k) sector on the
 *    computational basis, for block-diagonal ED / DMRG warm-starts. */
int dim_sector = irrep_sg_adapted_basis_at_k(
    lg, mu_k, n_sites, /*local_dim=*/2, basis_out, /*n_max=*/basis_capacity);
```

**Character menu** (`irrep_lg_named_irrep_t` in `config_project.h`):

| Name | C_6v | C_3v | C_2v | C_4v |
|---|:-:|:-:|:-:|:-:|
| `IRREP_LG_IRREP_A1` | ✓ | ✓ | ✓ | ✓ |
| `IRREP_LG_IRREP_A2` | ✓ | ✓ | ✓ | ✓ |
| `IRREP_LG_IRREP_B1` | ✓ | — | ✓ | ✓ |
| `IRREP_LG_IRREP_B2` | ✓ | — | ✓ | ✓ |
| `IRREP_LG_IRREP_E`  | — | ✓ | — | — |
| `IRREP_LG_IRREP_E1` | ✓ | — | — | — |
| `IRREP_LG_IRREP_E2` | ✓ | — | — | — |
| `IRREP_LG_IRREP_E_C4V` | — | — | — | ✓ |

Invalid (name, little-group) pairs return `NULL` with a clean diagnostic.

### 6.2 Group-element site permutations (exists — the key handshake primitive)

The public API already exposes exactly what Moonlab's MPO builder needs:

```c
/* For each g ∈ [0, order_G), write the site-permutation π_g into perm[].
 * Site i maps to perm[i] under g. Perfect input for an MPO that acts
 * as a permutation of physical indices. */
void irrep_space_group_permutation(const irrep_space_group_t *G, int g,
                                   int *perm /*[num_sites]*/);

/* Forward apply on a length-num_sites real configuration. */
void irrep_space_group_apply_config(const irrep_space_group_t *G, int g,
                                    const double *in, double *out);
```

**Moonlab integration path**: the composite projector at `(k, μ_k)` is

```
P_{k,μ_k} = (d_{μ_k} / |G|) · Σ_{g ∈ G} conj(χ_{μ_k}(g)) · e^{-i k · t_g} · π_g
```

where `π_g` is the site permutation, `t_g` the space-group fractional
translation for element `g`, `χ_{μ_k}` the little-group character
zero-extended over G. Moonlab builds this as an MPO once per (k, μ_k);
applying it to an MPS is a standard bond-dimension cost `O(|G| · D^3 · N)`.

**Projector-weight export** (landed — `config_project.h`):

```c
/* Emit (character · Bloch phase · 1/|G| · d_μ) weights, one per g ∈ G,
 * for a given (lg, μ_k). Caller pairs with irrep_space_group_permutation()
 * to assemble the MPO. Elements outside the little group get exact 0. */
IRREP_API int irrep_sg_projector_weights(
    const irrep_sg_little_group_t       *lg,
    const irrep_sg_little_group_irrep_t *mu_k,
    double _Complex                     *weights_out /*[order_G]*/);
```

Bit-exact identity tested against `irrep_sg_project_at_k` on C_6v / C_3v /
C_4v at Γ and K on kagome and square lattices, including 2D irreps where
some little-group weights are legitimately zero on χ=0 classes.

### 6.3 Entanglement diagnostics (exists)

For the topological-γ test (Kitaev–Preskill 2006), which discriminates
gapped Z₂ (γ = log 2) from gapless U(1) Dirac (γ = 0):

```c
#include <irrep/rdm.h>

/* Partial trace: ρ_A = Tr_B |ψ⟩⟨ψ| for a bipartition A ⊔ B of sites. */
irrep_status_t irrep_partial_trace(
    int num_sites, int local_dim, const double _Complex *psi,
    const int *sites_A, int num_sites_A,
    double _Complex *rho_A_out /*[(local_dim^num_A)^2]*/);

/* Entropies */
double irrep_entropy_vonneumann(const double _Complex *rho, int n);
double irrep_entropy_renyi(const double _Complex *rho, int n, double alpha);

/* Kitaev–Preskill γ from three-region construction. */
double irrep_topological_entanglement_entropy(
    double S_A, double S_B, double S_C,
    double S_AB, double S_BC, double S_AC, double S_ABC);
```

For DMRG on cylinders, Moonlab should extract `|ψ⟩` on the relevant
sub-Hilbert space (up to ≤ 30 sites bipartite) and call these directly.
For larger bipartitions, ρ_A stays in MPO land — out of libirrep's scope.

### 6.4 Hamiltonian apply-ops (exists)

For reference implementations and cross-validation vs. Moonlab's own:

```c
/* Heisenberg H = J · Σ_{⟨i,j⟩} S_i · S_j on a caller-supplied bond list. */
irrep_heisenberg_t *H = irrep_heisenberg_new(N, num_bonds, bi, bj, J);
irrep_heisenberg_apply(psi_in, psi_out, /*opaque=*/H);   /* for Lanczos */

/* J₁-J₂ (frustration-sweep ready). */
irrep_heisenberg_t *H12 = irrep_heisenberg_j1j2_new(
    N, num_nn, nn_i, nn_j, J1, num_nnn, nnn_i, nnn_j, J2);
```

These are `O(|bonds| · 2^N)` per call — dense, one full-vector pass per
Lanczos iteration. They exist mainly as (a) a reference oracle against
which Moonlab's MPO-based H-apply can be cross-checked on N ≤ 24, and
(b) the backbone of our own N ≤ 24 ED suite.

### 6.5 Lanczos eigensolver (exists)

Two variants in `rdm.h`:

```c
/* Compact: 3-vector recurrence, no reorth. Fast, but ghost-eigenvalue
 * artefacts on near-degenerate spectra (kagome singlet tower!). */
irrep_status_t irrep_lanczos_eigvals(
    void (*apply_op)(const double _Complex *, double _Complex *, void *),
    void *ctx, long long dim, int max_iters, int n_eigvals,
    const double _Complex *seed /* or NULL for random */,
    double *eigvals_out /*[n_eigvals]*/);

/* With full Gram–Schmidt reorth: memory cost max_iters · dim, but
 * ghost-eigenvalue-free. Use this on any sector where the low spectrum
 * is nearly degenerate (Marshall-rule singlet tower, any kagome sector). */
irrep_status_t irrep_lanczos_eigvals_reorth(...);
```

Moonlab's state-vector ED path can call these directly for cross-checks;
our own use is inside examples such as `kagome12_k_resolved_ed.c`.

---

## 7. What's missing (and who closes each gap)

### 7.1 On libirrep's side

**7.1.1 Projector-weight export** (landed — `config_project.h`).
- API: `irrep_sg_projector_weights(lg, mu_k, weights_out[order_G])`.
- Unblocks: Moonlab MPO assembly for `P_{k, μ_k}`.

**7.1.2 Orbit-representative canonical form** (primitives landed —
`space_group.h`, `config_project.h`):

```c
/* Apply g to a bitstring config (bit i of input → bit perm[i] of output).
 * ≤ 64 sites. Reuses cached G->perm, no allocation. */
IRREP_API uint64_t irrep_space_group_apply_bits(
    const irrep_space_group_t *G, int g, uint64_t config_in);

/* Canonicalise: rep_out = lex-minimum bitstring in orbit(config_in);
 * g_idx_out is the smallest g mapping config_in → rep_out. */
IRREP_API void irrep_sg_canonicalise(
    const irrep_space_group_t *G, uint64_t config_in,
    uint64_t *rep_out, int *g_idx_out);

/* Orbit size = |G| / |stabiliser|. Needed for normalisation of
 * representative-basis coefficients. */
IRREP_API int irrep_sg_orbit_size(const irrep_space_group_t *G, uint64_t config);
```

6835 assertions passing; invariants tested include idempotence, orbit
invariance, Burnside sum rule (Σ orbit_size over reps = C(N, Sz)) on
kagome 12 Sz=0 and square 3×3 p1.

**Representative-table + rep→index lookup** (landed — `config_project.h`):

```c
typedef struct irrep_sg_rep_table irrep_sg_rep_table_t;

IRREP_API irrep_sg_rep_table_t *irrep_sg_rep_table_build(
    const irrep_space_group_t *G, int popcount);
IRREP_API void     irrep_sg_rep_table_free(irrep_sg_rep_table_t *T);
IRREP_API long long irrep_sg_rep_table_count(const irrep_sg_rep_table_t *T);
IRREP_API uint64_t irrep_sg_rep_table_get(const irrep_sg_rep_table_t *T, long long k);
IRREP_API int      irrep_sg_rep_table_orbit_size(const irrep_sg_rep_table_t *T, long long k);
IRREP_API long long irrep_sg_rep_table_index(const irrep_sg_rep_table_t *T, uint64_t rep);
```

Enumerates via Gosper's hack (iterates only C(N, p) popcount-p configs, not
all 2^N), stores sorted, binary-search rep→index. 2514 assertions covering
Burnside sum rule at every popcount 0..12, rep canonicity, sorted-order
invariant, rep→index round-trip on all 924 popcount-6 configs on kagome 2×2,
non-rep / wrong-popcount rejection, popcount-edge sectors. Open-addressed
hash can later replace the binary search for hot-loop speedup.

**Sparse H apply on rep basis, trivial sector** (landed — `hamiltonian.h`):

```c
IRREP_API void irrep_heisenberg_apply_in_sector(
    const irrep_heisenberg_t *H,
    const struct irrep_sg_rep_table *T,
    const double _Complex *psi_in, double _Complex *psi_out);
```

Matches the `irrep_lanczos_eigvals` apply-callback signature. Scope: the
totally-symmetric sector of whichever space group `T` was built on —
`(p1, Γ)` is pure-momentum Γ + Sz, `(p6mm, Γ, A_1)` is the full kagome
A_1 sector, etc. Complexity: `O(num_bonds · |reps|)` vs the dense path's
`O(num_bonds · 2^N)`; the reduction factor is roughly `|G|` (432 on
kagome 6×6 p6mm).

Formula (Sandvik / Weiße-Fehske):
  - Diagonal: `y_u += c_zz · (aligned? +1 : −1) · x_u` per bond
  - Off-diagonal on anti-aligned `(i, j)`:
    `y_v += c_pm · √(orbit_u / orbit_v) · x_u` where `v = canon(u ^ (1<<i) ^ (1<<j))`

1813 assertions: dense-vs-sparse matrix-element agreement to 1e-13 on 4-site
square p1 and full 30×30 kagome 12 Γ-A_1 Sz=0 block, Hermiticity to 1e-13,
ground state E_0 = −5.328392 J reproduced by sparse Lanczos to 1e-5.

**Cached sector-binding** (landed — `hamiltonian.h`):

```c
typedef struct irrep_sg_heisenberg_sector irrep_sg_heisenberg_sector_t;
IRREP_API irrep_sg_heisenberg_sector_t *irrep_sg_heisenberg_sector_build(
    const irrep_heisenberg_t *H, const struct irrep_sg_rep_table *T);
IRREP_API void irrep_sg_heisenberg_sector_free(irrep_sg_heisenberg_sector_t *S);
IRREP_API void irrep_sg_heisenberg_sector_apply(
    const double _Complex *psi_in, double _Complex *psi_out, void *S);
```

Canonicalises every `(rep, bond)` transition once, stores as CSR
`(diag[n_reps], rowptr[n_reps+1], col[nnz], coef[nnz])`, and the
matvec becomes pure memory-bound sparse-matrix × vector with zero
canonicalisation in the hot loop.

**Measured on Apple M2 Ultra**, p6mm kagome:
- N=12, Γ-A₁ Sz=0: 30 reps, cached build <1 ms, Lanczos (60 iters, 4 eigs) <1 ms
- N=27, Γ-A₁ Sz=-½: 186 616 reps, cached build 12.4 s, Lanczos 1.2 s
  (dense path infeasible at 2 GiB / vector)

**Stabiliser primitive** (landed — `config_project.h`):

```c
IRREP_API int irrep_sg_stabiliser(const irrep_space_group_t *G, uint64_t config,
                                  int *out_indices);
```

Emits the list of g fixing config, in ascending g order. Building block
for the composite-weight variant (next).

**Full (k, μ_k) sparse apply** (landed — `hamiltonian.h`):

```c
IRREP_API irrep_sg_heisenberg_sector_t *irrep_sg_heisenberg_sector_build_at_k(
    const irrep_heisenberg_t *H, const struct irrep_sg_rep_table *T,
    const struct irrep_sg_little_group *lg,
    const struct irrep_sg_little_group_irrep *mu_k);

IRREP_API long long irrep_sg_heisenberg_sector_dim(
    const irrep_sg_heisenberg_sector_t *S);
```

For each rep `u` computes `σ_u = Σ_{s ∈ Stab(u)} conj(w_s)` (the sector-
basis norm squared up to a factor); reps with `σ_u = 0` are annihilated
by the projector and filtered out. Off-diagonal CSR coefficient is
`c_b · |G| · conj(w_{g_idx}) · √(σ_v/σ_u)`, complex in general. At
(Γ, A₁) this reduces to the trivial-sector formula `c_b · √(orbit_u/orbit_v)`.

**Verified**: kagome 12 Γ-B_1 sparse Lanczos → E₀ = −5.44487522 J
(matches the absolute KH ground state, see `kagome12_k_resolved_ed`).
Γ-A_1 → −5.32839240 J (matches the trivial-sector path).

**Worked examples** (in `examples/`):
- `kagome27_sector_ed` — full 8-sector (k, μ_k) ED on N=27 kagome in ~120 s
  (dense path infeasible at 2 GiB/vector). Reproduces the Dirac-cone
  signature: E_0/N = −0.6132 at K-A_1 vs −0.4300 at Γ-A_1.
- `j1j2_square_sector_ed` — N=16 p4mm J₁-J₂ at J₂/J₁ ∈ {0, 0.5, 1.0}.
  12 sectors × 3 couplings in 0.6 s. Pure-J₁ GS matches Schulz 1989
  to 3e-4. DQCP regime shows Γ→X sector crossover.

**2D-irrep D_μ(g) matrix machinery** (landed — `config_project.h`):

```c
IRREP_API int irrep_sg_little_group_irrep_matrix(
    const irrep_sg_little_group_irrep_t *mu_k, int i, double _Complex *D_out);
```

Named 2D irreps (C_3v E, C_6v E_1, C_6v E_2, C_4v E) now carry full
per-element D_μ(g) matrices, populated automatically by
`irrep_sg_little_group_irrep_named`. Unitarity + character-consistency
verified across all named irreps in 100+ assertions. 1D irreps return
their character as a 1×1 matrix via the same accessor (uniform API).

Moonlab can consume these matrices directly for MPO-based 2D-sector
projectors: `P_{k, E} = Σ_g (d_μ/|G|) [D_μ(g)]_{00}* · π_g · e^{-ik·t_g}`
assembles as a standard weighted permutation-MPO sum.

**Dense reference sector builder** (landed — `hamiltonian.h`):

```c
IRREP_API int irrep_sg_heisenberg_sector_build_dense(
    const irrep_heisenberg_t *H, const irrep_sg_rep_table_t *T,
    const irrep_sg_little_group_t *lg,
    const irrep_sg_little_group_irrep_t *mu_k,
    double _Complex *H_out, int n_max);
```

First-principles construction: for each rep `u` in `T`, build
`|ũ⟩ = P_μ|u⟩ / √σ_u` as an explicit 2^N dense vector using the
projector weights; sandwich H through pairs to populate `H_out[j,i] =
⟨ũ_j|H|ũ_i⟩`. Returns actual dim; rejects N > 24 (memory cap).

**Validated for both 1D and 2D**:
- Kagome 12 Γ-A_1: E_0 = -5.328392 J (matches sparse 1D)
- Kagome 12 Γ-B_1: E_0 = -5.444875 J (matches sparse 1D, absolute KH GS)
- Kagome 12 Γ-E_1: E_0 = -3.299516 J (matches `kagome12_k_resolved_ed`)
- Kagome 12 Γ-E_2: E_0 = -4.177552 J (matches published value)
- Hermiticity error < 1e-15 across all tested sectors (1D and 2D).

**Usage guidance for Moonlab**:
- **`sector_build_at_k`** is the unified entry point for both 1D and 2D
  irreps. 1D uses the efficient sparse character formula directly
  (O(|reps|) memory). 2D internally routes through the dense reference
  for correctness, packs the result as CSR — same consumer apply interface.
- **Size limits**:
  - 1D sectors: feasible to the rep-table size limit (~N = 30 on kagome,
    workstation RAM).
  - 2D sectors: feasible to N ≤ 24 (dense builder memory cap).
- **N > 24 with 2D irreps**: currently infeasible. Genuine sparse 2D
  formula remains open; dense is the working oracle for debugging when
  we revisit it. For Moonlab's 2D-sector K-point work at N = 27 kagome,
  DMRG via the projector-MPO path (using `irrep_sg_projector_weights` +
  `irrep_space_group_permutation`) is the current route.

### 6.6 Persistence formats (for Moonlab integration)

Two on-disk binary formats ship as part of the sparse stack. Both are
little-endian, not cross-architecture-portable without care. Use for
research-workflow caching; do not rely on them as long-term archive
formats until explicit endian-normalisation lands.

**Rep table** (`irrep_sg_rep_table_save/_load`, `~1.9 MB` at N=27):
```
offset  type        field
0       8B          magic = "IRREP_RT"
8       uint32      version = 1
12      int32       popcount
16      int32       num_sites         (checked against caller's G at load)
20      int32       group_order       (checked against caller's G)
24      int64       count
32      count×u64   reps (sorted ascending)
...     count×i32   orbit_sizes (parallel to reps)
```

**Sector binding** (`irrep_sg_heisenberg_sector_save/_load`):
```
offset  type        field
0       8B          magic = "IRREP_SB"
8       uint32      version = 1
12      uint32      flavour   (0 = real/trivial sector, 1 = complex/(k, μ_k))
16      int64       n_reps
24      int64       nnz
32      n_reps×f64      diag
...     (n_reps+1)×i64  rowptr
...     nnz×i64         col
...     nnz×{f64|2f64}  coef (real if flavour=0, complex if 1)
```

Measured (Apple M2 Ultra): **N=27 cold-start 28.2 s → warm-cache 1.1 s**
(26× end-to-end speedup), bit-exact eigenvalues across runs. See
`examples/ed_cache_demo.c` for the recommended caller-side pattern.

**7.1.3 Reciprocal-lattice basis export** (already landed — `lattice.h:98`):
`irrep_lattice_reciprocal_vectors(L, b1, b2)` — caller reads `b_i · a_j
= 2π δ_ij` primitives directly.

### 7.2 On Moonlab's side

**7.2.1 Symmetry-projector MPO constructor**.
- Input: `(order_G, permutations[order_G][num_sites], weights[order_G])`
  — the triple that libirrep vends after 7.1.1 lands.
- Output: an MPO that acts as `P_{k, μ_k}` on the MPS physical indices.
- The standard "sum of permutation MPOs with complex weights" construction
  — bond dimension grows linearly with the number of distinct image-preimage
  pairs per site, worst-case `O(|G|)` but typically much smaller after
  Moonlab's existing MPO compression.

**7.2.2 DMRG driver that warm-starts from a libirrep sector basis**.
- Input: the `basis_out` buffer from `irrep_sg_adapted_basis_at_k` on a
  small (N ≤ 24) cluster — a dense `dim_sector × 2^N` matrix whose
  columns are orthonormal vectors in the sector.
- Output: compress to an MPS using Moonlab's existing SVD / variational
  compression, use as DMRG initial state.

**7.2.3 Cylinder-geometry analog of wallpaper groups**.
- Kagome YC cylinders (Y-cylinder geometry, standard in the DMRG
  literature) have reduced symmetry — the Y-periodic axis is intact,
  the X axis terminates. The Y-translation subgroup is what Moonlab
  would want to project onto. libirrep's 2D `space_group_build` assumes
  full PBC in both axes; for cylinder geometry Moonlab should use the
  *translation subgroup only* piece of the API. (`irrep_sg_bloch_basis`
  already isolates translations from point group — so this exists.)

### 7.3 On SbNN's side

SbNN already consumes libirrep's Γ-projector via the existing handshake
(inbox `agent-notes/inbox-spin`). Once 7.1.1 lands SbNN can call through
the same primitive to target non-Γ sectors.

---

## 8. Phased integration roadmap

Each phase produces a publishable deliverable even if the project halts there.

### Phase 0 — toolchain handshake (≤ 1 week after 7.1.1 lands)

- Moonlab vendors libirrep 1.3.x via one of the three consumption paths.
- Moonlab writes a 10-line smoke test calling `irrep_abi_hash()`,
  `irrep_space_group_build()`, `irrep_space_group_permutation()`,
  `irrep_sg_projector_weights()` on a 12-site kagome cluster. Pins
  the ABI hash in CI.
- libirrep publishes 1.3.0-alpha.1 (or 1.3.1) with 7.1.1 and 7.1.3.

**Success criterion**: Moonlab CI green against libirrep 1.3.x; smoke
test produces `A_1`-projector weights that sum to `|G|` (the canonical
normalisation check).

### Phase 1 — reproduce benchmarks (2–3 weeks)

- Moonlab DMRG on YC4 kagome cylinder, several lengths. Reproduce
  Yan–Huse–White 2011 E_0/N ≈ −0.4386 J ± 0.0001.
- libirrep ED at N = 12, 18, 24 (existing, in `kagome*_k_resolved_ed.c`).
  Moonlab state-vector ED reproduces the same numbers — that's the
  floor that validates Moonlab's lattice-2D + Heisenberg plumbing is
  byte-correct.
- SbNN NQS on N = 24 torus, (Γ, B_1) sector via libirrep projector.
  Converge to within 1e-3 J of the ED answer.

**Success criterion**: three-way agreement on N = 24 torus and
Moonlab-DMRG vs. YHW paper on YC4.

### Phase 2 — close the gap (1–2 months)

- libirrep 1.3.x lands 7.1.2 (orbit representatives + sparse
  Heisenberg apply). Torus ED reaches N = 27, 30, 36 with
  symmetry-resolved Lanczos.
- Moonlab DMRG on YC6, YC8, YC10 — full (k, μ_k) block structure from
  libirrep MPOs. Extract γ via Kitaev–Preskill construction at each
  width.
- SbNN NQS on N = 36, 48, 72 torus, all low-lying sectors.

**Success criterion**: three independent (width, torus-size) pairs
produce compatible γ extrapolation (log 2 or 0); scaling-collapse
plots demonstrate this.

### Phase 3 — verdict (1 month)

- Moonlab TDVP on the Phase-2 ground state. Compute `S(q, ω)` at K
  via KPM. Dirac-cone vs hard-threshold distinguishes the phases.
- Joint NQS + DMRG + torus-ED fit: spin gap vs 1/L, γ, `S(q, ω)` at K.
  Two of three observables agreeing ⇒ the answer. Publish.

---

## 9. Release coordination

- **libirrep cuts 1.3.0-alpha.N** (or 1.3.N) for each Moonlab-consumable
  change. ABI is additive within 1.3.x; Moonlab pins the hash.
- **Release cadence**: aligned with Phase-1 / 7.1.1. Expect the first
  post-alpha cut within 2 weeks of Moonlab deciding to integrate.
- **Coordination channel**: we can spin up `agent-notes/inbox-moonlab`
  parallel to the existing `inbox-spin` / `inbox-irrep` if Moonlab
  wants structured notes.

---

## 10. Contact points & references

- Repo: `github.com/tsotchke/libirrep`
- Current tag: `v1.3.0-alpha`
- CI: 5 triples, ASan + UBSan green on every commit, 7 fuzz targets
- Docs: `docs/{API,DESIGN,METHODS,PHYSICS_APPENDIX,PHYSICS_RESULTS,REFERENCES}.md`;
  Doxygen HTML bundled in release tarballs
- KH example in tree: `examples/kagome12_k_resolved_ed.c` — full (k, μ_k)
  block ED on 12-site kagome (18 sectors, E_0 = −5.44487522 J at (Γ, B_1))

For the deliverables in §7.1, either file issues against the libirrep
repo or open notes under `agent-notes/inbox-irrep/`. Expect same-week
turnaround on 7.1.1 and 7.1.3; 7.1.2 is the 1–2-week item.
