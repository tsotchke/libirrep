# libirrep QEC scoping — modern code families and supported targets

## Purpose

Identify the QEC code families libirrep should support to be useful for
current research (2023–2026), design an extensible abstract API around
them, and prioritise the implementation order. This is the gating
document for issue #299; the toric-code prototype (`toric_code.{h,c}`)
is the first instance.

## Code-family survey (current as of Q2 2026)

Modern QEC research has fragmented well beyond the original toric/surface
paradigm. The families that need first-class support, in priority order
of *near-term experimental relevance*:

### Tier 1 — workhorses

1. **Surface code** [Bravyi-Kitaev 1998, Fowler-Mariantoni-Martinis-Cleland
   2012]. Planar variant of the toric code with rough X / smooth Z
   boundaries. Encodes 1 logical qubit on a `d × d` patch with distance
   `d`. Workhorse of current Google and IBM superconducting devices.
   Already the single most-implemented QEC code in hardware experiments.

2. **Toric code** [Kitaev 2003] — already implemented in
   `toric_code.{h,c}`. Reference implementation; periodic boundary
   conditions; 2 logical qubits per torus.

3. **Color code** [Bombin–Martin-Delgado 2006]. Triangular lattice with
   3-coloring of faces; each face is a stabilizer of *all* qubits on its
   boundary. Transversal Clifford on the planar triangular code is the
   key advantage for fault-tolerant computation.

### Tier 2 — recent breakthroughs (2021–2024)

4. **Bivariate bicycle (BB) codes** [Bravyi-Cross-Gambetta-Maslov-Rall-
   Yoder, Nature 627 (2024) 778]. Parameterised by two polynomials
   `A(x,y), B(x,y) ∈ F₂[x,y]/(xˡ-1, yᵐ-1)`. Parity-check matrices
   `H_X = [A | B]`, `H_Z = [Bᵀ | Aᵀ]`. The `[[144, 12, 12]]` instance
   significantly outperforms surface code per logical qubit at IBM
   threshold (~0.7%) — the first clear "qLDPC beats surface" result on
   physical hardware.

5. **Honeycomb Floquet code** [Hastings-Haah, PRX Quantum 2 (2021)
   030328]. Honeycomb lattice with 3-coloring of edges; each round
   measures 2-body Paulis on edges of one color; after a 3-round cycle
   the instantaneous stabilizer group returns. Logical qubits encoded
   in the dynamics, not a static stabilizer group. Threshold ~0.2-0.3%
   under measurement-induced noise.

6. **CSS Floquet codes** [Davydova-Tantivasadakarn-Balasubramanian, PRX
   Quantum 4 (2023) 020341]. Generalises honeycomb Floquet to CSS
   structure with explicit logical operators.

### Tier 3 — qLDPC frameworks

7. **Hypergraph product codes** [Tillich-Zemor, IEEE Trans. Inf. Theory
   60 (2014)]. General CSS construction from a pair of classical
   parity-check matrices. Distance scaling `d ~ √n` rather than `√n`
   for surface code at fixed rate.

8. **Lifted product codes** [Panteleev-Kalachev, IEEE Trans. Inf. Theory
   68 (2022)]. Asymptotically good qLDPC: `[[n, kn, dn]]` with
   `k, d → constants > 0` as `n → ∞`. The first asymptotically good
   quantum codes.

### Tier 4 — research frontier (2024–2026)

9. **Single-shot qLDPC** [Quintavalle-Vasmer-Roffe-Campbell, PRX Quantum
   2 (2021) 020340; recent extensions]. Codes with single-shot syndrome
   extraction — eliminate the need for repeated measurements per round.

10. **Concatenated qLDPC** — recent work composing qLDPC outer codes
    with surface inner codes for fault-tolerant compilation.

11. **3D toric / X-cube / fracton codes** [Vijay-Haah-Fu 2016]. 3D
    stabilizer codes with restricted-mobility excitations; relevant for
    self-correcting quantum memories.

### Tier 5 — theoretical / specialised

12. **HaPPY holographic code** [Pastawski-Yoshida-Harlow-Preskill 2015].
    Tensor-network code on hyperbolic tilings; tied to AdS/CFT.

13. **Bacon-Shor / subsystem codes** [Bacon 2006]. Gauge-fixing variant;
    9-qubit example.

14. **Color-code subsystem variants** [Bombin 2010]. Includes triangular
    color code and gauge color code.

## Why this scoping matters

A naïve API exposes individual codes (`toric_init`, `surface_init`,
`color_init`, ...) — but qLDPC codes are defined by **arbitrary sparse
parity-check matrices** `H_X, H_Z`, and Floquet codes by **time-dependent
measurement schedules**. Hard-coding each family kills extensibility.

The right design is to have:

- A general **parity-check code** abstraction (`irrep_css_code_t`):
  a CSS code defined by two binary parity-check matrices `H_X` (rows are
  X-stabilizers) and `H_Z` (rows are Z-stabilizers), with the
  CSS-orthogonality constraint `H_X · H_Zᵀ = 0` as a runtime check.
  This subsumes toric, surface, color, hypergraph product, lifted product,
  bivariate bicycle.

- A **Floquet code** abstraction (`irrep_floquet_code_t`):
  a list of measurement rounds, each specifying a list of Pauli operators
  (their support + Pauli type) to measure. The "instantaneous stabilizer
  group" after each round is a derived object. This subsumes honeycomb
  Floquet, CSS Floquet.

- A **stabilizer-group** abstraction (`irrep_stabilizer_group_t`):
  list of Pauli operators (support + type) that pairwise commute, with
  utilities for: commutator parity, syndrome extraction at a given error,
  logical operator enumeration, distance computation.

## Implementation order

The plan is to build these in order of (complexity × utility):

| # | Module | What | Status | Est. LOC |
|---|---|---|---|---|
| 1 | `toric_code` | concrete reference, 2D periodic | **DONE** | 200 |
| 2 | `stabilizer_group` | abstract Pauli-string CSS / general | next | 400 |
| 3 | `surface_code` | planar with boundaries | next | 250 |
| 4 | `css_code` | general CSS from H_X, H_Z matrices | after #2 | 350 |
| 5 | `bivariate_bicycle` | BB[[n, k, d]] family | after #4 | 200 |
| 6 | `floquet_code` | measurement-schedule abstraction | medium-term | 500 |
| 7 | `honeycomb_floquet` | concrete instance of #6 | after #6 | 300 |
| 8 | `color_code` | triangular color code | medium-term | 300 |
| 9 | `hypergraph_product` | qLDPC from classical (H_a, H_b) | medium-term | 200 |
| 10 | `lifted_product` | asymptotically-good qLDPC | long-term | 400 |
| 11 | `subsystem_code` | gauge group + stabilizer group | long-term | 350 |
| 12 | `single_shot_lifter` | upgrade CSS to single-shot | ✅ shipped: `irrep_single_shot_lift` (F₂ left-nullspace via Gaussian elimination on [H \| I_m]) | 500 |

## API design highlights

The `stabilizer_group` abstraction will be the foundation. A Pauli
operator on n qubits is a triple `(x_support, z_support, sign)` where
`x_support, z_support ∈ {0,1}^n` are the X and Z parts (in symplectic
representation) and `sign ∈ {+1, -1, +i, -i}`.

Two Paulis `P, Q` commute iff their symplectic inner product is 0:
`P·Q = Q·P ⟺ x_P · z_Q + z_P · x_Q = 0 (mod 2)`.

This single primitive (`irrep_pauli_commute`) implements stabilizer
commutativity for ALL the CSS codes in the survey above. The toric
code's "shared-edge parity" is just a specialisation of this for the
specific (vertex × plaquette) Pauli supports.

## What this enables

Once Tier 1+2 codes (surface, color, BB) are implemented behind the
unified API:

- **Decoder benchmarking** across families on a common substrate
  (BP+OSD, MWPM, neural decoders all need the same syndrome
  extraction primitives)
- **Fault-tolerant compiler** experiments: lattice surgery between
  surface code and BB code, color code transversal Clifford
  composed with BB qLDPC outer code
- **Threshold studies** comparable across families (surface vs BB
  vs Floquet at fixed physical error rate)
- **Logical operator pushes** for runtime fault-tolerant gate synthesis

## Primary references

- Kitaev, *Ann. Phys.* 303 (2003) 2 [quant-ph/9707021] — toric code
- Bravyi-Kitaev (1998) — surface code
- Bombín-Martín-Delgado (2006) — color code
- Hastings-Haah, *PRX Quantum* 2 (2021) 030328 — honeycomb Floquet
- Davydova-Tantivasadakarn-Balasubramanian, *PRX Quantum* 4 (2023)
  020341 — CSS Floquet
- Bravyi-Cross-Gambetta-Maslov-Rall-Yoder, *Nature* 627 (2024) 778 —
  bivariate bicycle codes
- Tillich-Zemor, *IEEE Trans. Inf. Theory* 60 (2014) — hypergraph product
- Panteleev-Kalachev, *IEEE Trans. Inf. Theory* 68 (2022) — asymptotically
  good qLDPC
- Quintavalle-Vasmer-Roffe-Campbell, *PRX Quantum* 2 (2021) 020340 —
  single-shot qLDPC
- Pastawski-Yoshida-Harlow-Preskill, *JHEP* 06 (2015) 149 — HaPPY code
- Vijay-Haah-Fu, *Phys. Rev. B* 94 (2016) 235157 — X-cube / fracton
- Bombín, *Phys. Rev. A* 81 (2010) 032301 — gauge color code

## Adjacency to existing libirrep modules

- `lattice.h` (sites only) → `toric_code.h` (edges) → `stabilizer_group.h`
  (Pauli supports) is the lowest-level chain.
- `space_group.h` and `point_group.h` give the symmetry framework that
  Floquet codes' time-cycle and color codes' 3-coloring can hook into.
- `irrep.h` and `recoupling.h` are orthogonal to QEC but share the
  combinatorial-symmetry infrastructure that later qLDPC code-search
  algorithms will use.
