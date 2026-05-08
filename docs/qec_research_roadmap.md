# libirrep QEC research-grade development roadmap

## Status as of this session

- **17 production modules landed**, all with CSS-orthogonality and stabilizer-group commutativity verified, plus 7 codes with brute-force distance match (Steane, [[5,1,3]], surface d∈{2,3,4}, 2×2 toric, 3D toric L=2).
- **Sanitiser-clean**: ASan + UBSan green on all 62 test suites (30,363 wigner-d assertions + new QEC suites).
- **ICC-indexed**: 404 source files; the new QEC modules are now in the repo memory under `src/`, `include/irrep/`, `tests/`.
- **Adversarial production audit complete**, 7 high/medium-severity bugs found and fixed in this turn:
  1. `irrep_pauli_multiply` aliasing UB (`out == p`/`out == q`) — fixed by reading source words into locals first.
  2. `irrep_qec_distance_brute` `int total *= 3` overflow at weight ≥ 20 — switched to `uint64_t`.
  3. `lifted_product.h` H_Z formula out-of-sync with σ-antipode implementation — header corrected.
  4. `irrep_pauli_set` Y vs X·Z sign-bit convention divergence — header doc clarified.
  5. `irrep_single_shot_code_new` zombie-state on consumed CSS — full bookkeeping zeroed.
  6. `pm_offset` int overflow at large lifts — switched to `size_t`.
  7. `irrep_css_code_new` over-rejection of `m_X = 0` or `m_Z = 0` — degenerate sides now allowed.

## Deferred research-track items (this roadmap)

The remaining four research-grade items, all needing careful figure-decoding from the
literature:

| ID | Module | Source convention | Difficulty |
|----|--------|-------------------|------------|
| R1 | `[[19, 1, 5]]` hexagonal (6.6.6) triangular color code | Bombín-Martín-Delgado 2006; Landahl-Anderson-Rice 2011 Fig. 3 | Medium |
| R2 | `[[17, 1, 5]]` square-octagon (4.8.8) triangular color code | Pogorelov 2024 (Quantinuum); Krinner-Lacroix 2024 (IBM) | Medium |
| R3 | 3D Bombín gauge color code | Bombín 2010 PRA 81 032301 | High |
| R4 | Multi-tile HaPPY hyperbolic network | Pastawski-Yoshida-Harlow-Preskill 2015 JHEP 06 149 | High |

## Architectural framework

A single new module gets us most of the way there:

### `generic_color_code` — face-list-driven CSS

```c
typedef struct {
    int n_qubits;        /* Number of vertices (data qubits). */
    int n_faces;         /* Number of faces. */
    int *face_sizes;     /* face_sizes[f] = weight of face f. */
    int **face_qubits;   /* face_qubits[f][k] = k-th qubit of face f. */
    int *face_color;     /* face_color[f] ∈ {0=R, 1=G, 2=B}. Optional. */
} irrep_color_lattice_t;

irrep_status_t irrep_generic_color_build(const irrep_color_lattice_t *L,
                                         irrep_css_code_t *out);
```

This subsumes Steane, `[[19, 1, 5]]` hex, `[[17, 1, 5]]` 488, and any future 2D color code. Each face hosts both an X-stab and a Z-stab on its boundary qubits. CSS orthogonality follows from the fact that on a 3-colorable trivalent lattice every pair of faces meets in 0 or 2 qubits.

### Concrete instance: R1 `[[19, 1, 5]]` hex

19 qubits at vertices of a triangular patch of the honeycomb lattice with side L=2. Counts (general formula `|V_L| = 3L² + 3L + 1`, `|F_L| = 3L(L+1)/2`):

| L | distance d | n_qubits | n_faces | k |
|---|------------|----------|---------|---|
| 1 | 3 (Steane) | 7        | 3       | 1 |
| 2 | 5          | 19       | 9       | 1 |
| 3 | 7          | 37       | 18      | 1 |
| 4 | 9          | 61       | 30      | 1 |

Implementation: hard-coded face list for L=2, with a generator function for arbitrary L.

### Concrete instance: R2 `[[17, 1, 5]]` 488

17 qubits at vertices of a square-octagon (4.8.8) tiling patch. The interior octagon has weight 8; the 4 squares around it have weight 4; with truncated boundary squares to round out 8 total faces. Used by IBM and Quantinuum for distance-5 hardware demos in 2024.

Implementation: hard-coded face list per Pogorelov 2024 Fig. 2.

### R3 3D Bombín gauge color code

Subsystem code, not a static CSS. Lattice is the **rectified 3-cubic** tiling (also called the "3-colex"): each cell is a tetrahedron with 4-coloured vertices. Stabilizers come from cells (X-type and Z-type, both on the same 4-vertex support, weight 4). Gauge generators come from edges (weight-2). Smallest non-trivial instance: 3-torus with side L=1, ~15 qubits.

Implementation: extend `irrep_subsystem_code_t` with the cell + edge generator lists, and verify the stabilizer subgroup is in the centraliser of the gauge group. This is a **research-track** module — full lattice generation requires more work than fits in one session.

### R4 Multi-tile HaPPY network

Recursive `{5, 4}` hyperbolic tiling, each tile = `[[5, 1, 3]]` perfect tensor. Tiles in layer ℓ are contracted with their parent layer-(ℓ-1) tiles via shared boundary legs.

For depth `d` (number of tile layers from the central tile outward):
- Layer 0: 1 central tile (5 boundary legs + 1 bulk leg).
- Layer 1: 5 tiles attached to layer-0 boundary, each contributing new boundary legs.
- Layer ℓ: ~5·(3^(ℓ-1)) tiles (exponential growth, hyperbolic).

For `d = 1` (just the central tile): `[[5, 1, 3]]`. For `d = 2`: a 6-tile network with several boundary qubits. For `d = 3`: substantially more.

The HaPPY code is built as the *contraction* of the tensor network. After contraction, the
resulting object is a **stabilizer isometry** from bulk qubits to boundary qubits. We can
materialise it as an `irrep_stabilizer_group_t` on the boundary qubits.

Implementation: `irrep_happy_network_build(int depth, irrep_stabilizer_group_t *out)`. Internally enumerate tile placements, track boundary-leg ↔ tile-leg matching, and contract.

## Implementation order

This session targets **R1 + R2 + framework module** (the medium-effort items). R3 and R4 get *foundational primitives* delivered + explicit references + skeleton code; full-distance research instances follow later.

| Order | Item | Concrete deliverable | Complexity |
|-------|------|----------------------|------------|
| 1     | `generic_color_code` module | API + Steane re-derivation as sanity check | Low |
| 2     | R1: `[[19, 1, 5]]` hex via generic | n=19, 9 hexagons, distance verified | Medium |
| 3     | R2: `[[17, 1, 5]]` 488 via generic | n=17, 8 faces, distance verified | Medium |
| 4     | R3: 3D Bombin lattice primitive | Cell-vertex enumeration, stab + gauge generators on `L=1` 3-torus | Medium-high |
| 5     | R4: HaPPY depth-2 network | 6-tile contraction → stabilizer group on the boundary | Medium-high |

**Scope cap**: this turn ships items 1, 2, 3 (the explicit small-distance color codes), 4
foundation, 5 foundation. Full research-distance instances (e.g. `[[37, 1, 7]]` hex)
generated by formula + tested with brute-force distance up to weight 7.

## Verification plan

Each new code must pass:
- CSS-orthogonality (`irrep_css_code_verify` returns IRREP_OK).
- Materialised stabilizer-group full-pairwise commutativity.
- Brute-force distance match to expected (where computationally tractable: n ≤ 25).
- ASan + UBSan clean on all new tests.

Larger-distance instances are validated by **structural counts** (n_qubits, n_faces, k) matching the analytical formulas, plus CSS-orthogonality.

## References (full bibliography for this roadmap)

- **Bombín-Martín-Delgado**, *Topological quantum distillation*, PRL 97 (2006) 180501 — original 2D color code on the {6,3} hex lattice.
- **Landahl-Anderson-Rice**, *Fault-tolerant quantum computing with color codes*, arXiv:1108.5738 (2011) — `[[19,1,5]]` and `[[37,1,7]]` explicit layouts.
- **Bombín**, *Topological subsystem codes*, PRA 81 (2010) 032301 — 3D gauge color code.
- **Pogorelov et al. (Quantinuum)**, *Experimental fault-tolerant code-switching*, arXiv:2403.13732 (2024) — 488 `[[17,1,5]]` instance.
- **Pastawski-Yoshida-Harlow-Preskill**, *Holographic quantum error-correcting codes*, JHEP 06 (2015) 149 — HaPPY framework.
- **Hayden-Nezami-Qi-Thomas-Walter-Yang**, *Holographic duality from random tensor networks*, JHEP 11 (2016) 009 — HaPPY generalisations.
- **Kubica**, *The ABCs of the color code*, PhD thesis Caltech 2018 — comprehensive figures for all triangular hex distances.
