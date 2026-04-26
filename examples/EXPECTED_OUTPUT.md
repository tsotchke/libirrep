# Example outputs and reproducibility

Every non-toy example below takes a fixed RNG seed (inline in the source)
and produces deterministic output. This file documents the expected
values so a reader building from a fresh clone can verify the library
numerically matches this tree.

All runs measured on Apple M2 Ultra (aarch64, macOS 15). Wall-clock
times scale linearly with CPU single-thread perf; the numerical
outputs are platform-independent to the documented tolerance.

---

## `kagome12_ed`

12-site kagome Heisenberg S=½ ED via dense power iteration. Hilbert
space 2¹² = 4096.

**Seed:** uniform all-ones in the S_z = 0 subspace (dense power iteration
is not seed-sensitive on a non-degenerate ground state).

**Expected output (tolerance 1e-8 unless noted):**

```
E_0               = -5.44487522 J
E_0 / N           = -0.45374 J
sector            = B₁ (p6mm Γ-irrep)
⟨S_i · S_j⟩_NN    = -0.226870 J
S(Γ)              =  0.000000    (translation-invariant)
S(M_a)            =  0.395014
S(M_b)            =  0.395014
S(K)              =  0.395014    (C₃-symmetric across M / K)
spin gap Δ_S      =  0.382700 J  (power iteration in S_z = 1)
γ (bits)          = -0.476620    (finite-size dominated at 12 sites;
                                  not a physical γ estimate)
```

Cross-reference: Elser 1989, Lecheminant 1997.

## `kagome18_ed`

18-site kagome (2×3 torus) ED via power iteration. Hilbert space 2¹⁸.

**Seed:** `uint64_t rng = 0x123456789abcdefULL` (S_z = 0);
`uint64_t rng = 0xbeef1234beefULL` (S_z = 1 for triplet).

**Expected output (tolerance 1e-6):**

```
E_0 / nb          = -0.223533 J
S(1,0)            =  0.434779
S(0,1)            =  0.188174
S(1,1)            =  0.339669
S(0,2)            =  0.461956
S(1,2)            =  0.631402
E_triplet         = -7.76370191 J
spin gap Δ_S      = +0.28349302 J
```

Cross-reference: Waldtmann et al. 1998, Läuchli 2011.

## `kagome24_ed`

24-site kagome (2×4 torus) via sparse Lanczos. Hilbert space 2²⁴ ≈ 1.7 × 10⁷.

**Seed:** `uint64_t rng = 0x13579bdf13579ULL` (filled into all S_z = 0
basis states). Restricts Lanczos to the conservation-law-restricted
subspace of dimension C(24, 12) = 2 704 156.

**Expected output (tolerance 1e-4, Lanczos converges to fewer digits):**

```
E_0 / N           = -0.448330 J
spin gap Δ_S      =  0.2639 J
iterations        =  ~50 (max 80)
wall clock        =  ~2.5 min on M2 Ultra
peak memory       =  ~768 MB (3 Lanczos state vectors)
```

Cross-reference: Läuchli 2011 (comparable 2×L_y torus geometries).

## `kagome12_symmetry_ed`

12-site kagome, block ED across all six p6mm Γ-irreps.

**Deterministic:** builds symmetry-adapted basis per irrep via character
weighting. No RNG.

**Expected output:**

```
Γ-irrep multiplicities sum to = 593
  (= orbit-rep count under the full space group)
sector multiplicities by block:
  A₁ = 144,  A₂ = 76,   B₁ = 74,   B₂ = 74
  E₁ = 90,   E₂ = 135
ground-state sector            = B₁
ground-state energy            = -5.44487522 J
  (identical to the global power-iteration result)
```

## `kagome12_k_resolved_ed`

12-site kagome, block ED across every `(k, μ_k)` space-group sector —
Γ (C_6v, 6 irreps) plus the three M-points (C_2v, 4 irreps each).
Exercises the full libirrep 1.3 M15a little-group-projector stack:
`irrep_sg_little_group_build`, `_irrep_new`, `_adapted_basis_at_k`,
`element_matrix` introspection.

**Deterministic:** hand-coded C_6v / C_2v character tables,
conjugacy-classification keyed to the space-group builder's point-op
indexing. No RNG.

**Expected output (per-sector, J = 1):**

| k   | μ   | dim | E_min       | E_max       |
|-----|-----|----:|------------:|------------:|
| Γ   | A_1 | 144 |    -5.328392 |    +6.000000 |
| Γ   | A_2 |  76 |    -4.962435 |    +3.138647 |
| Γ   | B_1 |  74 |    **-5.44487522** |    +2.256488 |
| Γ   | B_2 |  74 |    -3.676094 |    +4.000000 |
| Γ   | E_1 |  90 |    -3.299516 |    +3.250000 |
| Γ   | E_2 | 135 |    -4.177552 |    +3.642602 |
| M_a | A_1 | 270 |    -4.902343 |    +4.224745 |
| M_a | A_2 | 210 |    -5.165293 |    +3.272791 |
| M_a | B_1 | 264 |    -4.972554 |    +5.000000 |
| M_a | B_2 | 264 |    -5.062207 |    +3.615360 |

(M_b and M_c are equivalent to M_a under C_6; same dimensions and
eigenvalues. Frobenius multiplicity sum = 593 + 3·1008 = 3617.)

Ground-state assignment: **(Γ, B_1)** at E_0 = −5.44487522 J, matching
the Γ-only projector and the dense power iteration.

## `kagome_a1_projection`

108-site kagome (6×6 × 3 sublattices) — the full
symmetric-neural-quantum-state target cluster. Constructs the lattice,
the p6mm space group (432 elements), the orbit enumeration of a trial
classical spin configuration, and the character-weighted A₁ projection.

**Seed:** classical-antiferromagnetic-like pattern, see source.

**Expected output:**

```
lattice           : kagome 6×6, 108 sites, 216 NN bonds, 216 NNN bonds
space group       : p6mm, order 432, point subgroup order 12
orbit             : 432 images of the 108-long configuration vector
A₁ projection time: a few ms
```

This example demonstrates the substrate is built; it does not
diagonalize 2^108. Full ED at this size is intractable (dim 3.25 × 10³²);
use of the substrate is downstream (VMC / symmetric-NQS).

## `spin_half_rotation`, `equivariant_mlp`, `torque_net_tp_paths`

Small single-check examples. Each prints one numeric assertion and a
PASS/FAIL. Expected output:

```
spin_half_rotation      : |↑⟩ → -|↑⟩ under 2π rotation    (Berry phase)
equivariant_mlp         : norm equivariance at 1e-12
torque_net_tp_paths     : 2 paths enumerated for l=1⊗l=1→l=1+2
```

## 3D Bravais lattice + cubic-point-group + DMI examples (1.3.0-alpha)

### `lattice3d_demo`

Geometry walkthrough for SC, BCC, FCC, Diamond, Pyrochlore. Verified
expected coordination numbers and NN distances:

```
SC          1 site/cell, 6 NN/site,  NN=1
BCC         2 sites/cell, 8 NN/site, NN=√3/2 ≈ 0.866
FCC         4 sites/cell, 12 NN/site, NN=√2/2 ≈ 0.707
Diamond     8 sites/cell, 4 NN/site, NN=√3/4 ≈ 0.433
Pyrochlore  16 sites/cell, 6 NN/site, NN=√2/4 ≈ 0.354
```

### `lattice3d_heisenberg_ed`

Three full-space Heisenberg AFM ED runs:

```
SC  2³ (N=8,  dim=256)   E_0 = -4.820089 J  (cube graph)
BCC 2³ (N=16, dim=65536) E_0 = -20.000000 J (K_{8,8} closed form, exact)
SC  4×2×2 (N=16, dim=65536) E_0 = -11.228 J
```

The BCC 2³ result hits the K_{8,8} closed form `−n(n+2)/4 = −8·10/4
= −20 J` to all printed digits — validates that the geometry / bond
list / Heisenberg apply / Lanczos pipeline are mutually consistent.

### `lattice3d_sector_ed`

Γ-momentum sector ED on BCC 2³. Sector dimension 1670 vs. full
65536 (39× reduction); recovers `E_0 = −20.000000 J` to all
printed digits.

### `lattice3d_kspace_ed`

k-resolved ED on SC 2³, scanning all 8 BZ corners; sum of dims = 70
= C(8,4); minimum E_0(k) = −4.820089, matches full-space ED to
4×10⁻⁷ (Lanczos precision).

### `pyrochlore16_heisenberg`

16-site pyrochlore Heisenberg AFM:

```
E_0 = -8.809084 J    (per site -0.551, per bond -0.184)
E_1 = E_2 = E_3 = E_4 = -8.440672 J   ← 4-fold cubic-multiplet
Δ_01 = 0.368412 J
```

### `pyrochlore16_j1j2`

J₁-J₂ phase sweep. Notable points:

```
J₂/J₁ = 0.0   E_0 = -8.809084   3-fold degenerate E_1=E_2=E_3
J₂/J₁ = 0.5   E_0 = -9.450      Δ_01 = 0.039 (near phase boundary)
J₂/J₁ = 1.0   E_0 = -12.000000  4-fold degenerate (hidden symmetry)
```

### `pyrochlore16_correlations`

Spin-spin correlation function on the GS:

```
⟨S_0 · S_NN⟩ = -0.183523   (all 6 NN equivalent to 13-digit precision)
⟨S_0 · S_r⟩  = +0.039015   (all 9 non-NN sites uniformly identical)
sum-rule:  48 × ⟨C_NN⟩ = -8.809084 = E_0  ✓ to 10⁻¹³
```

### `cubic_crystal_field`

Crystal-field decomposition of l-orbitals under O_h and T_d:

```
1x0e (s)   → A1g           1x1e   → T1g (axial vector under O_h)
1x1o (p)   → T1u (polar)   1x2e   → Eg + T2g (textbook eg / t_2g)
1x3o (l=3) → A2u + T1u + T2u    1x3e → A2g + T1g + T2g
```

### `dmi_pyrochlore_pattern`

Bond-exchange-tensor decomposition for pyrochlore NN bond
{(0,0,0), (¼,¼,0)}:

```
Site sym  DMI   D direction          J^s     |D · bond̂|
O_h       0-dim D = 0                3-dim   (no DMI)
O         1-dim D ∥ bond (1,1,0)/√2  3-dim   1.000000  ← Moriya rule e
T_d       1-dim D ⊥ bond             3-dim   0.000000
```

### `dmi_kagome_pattern`

Kagome NN bond under in-plane symmetries (relevant to RT magnets
Fe₃Sn₂ / Mn₃Sn / Co₃Sn₂S₂):

```
D_6  (chiral hex):   DMI 1-dim, D ∥ bond
C_3v (with σ_v):     DMI 1-dim, D ⊥ bond (in-plane, perp to σ_v)
D_3  (chiral, 6 el): DMI 1-dim, D ∥ bond
```

### `test_dmi`

All five Moriya-1960 rules independently verified, plus multi-symmetry
intersection cases, the symmetric-tensor analyzer, the magnetic-point-
group antiunitary sign rules, and the three-spin scalar chirality
selection rule. **65/65 assertions pass**.

### `kagome_triangle_chirality`

Scalar chirality verdict on the kagome triangle under candidate site
symmetries, mapping each verdict onto a real RT kagome magnet:

```
Identity-only:                χ ALLOWED
C₃ about centroid:            χ ALLOWED  (cyclic perm × proper)
σ_h IN triangle plane:        χ FORBIDDEN  (centrosymmetric kagome stacking)
T·σ_h "magnetic mirror":      χ ALLOWED  (Mn₃Sn-type 120° AFM route to TH)
C_3v (C₃ + 3 σ_v):            χ ALLOWED  (each σ_v transposes + improper = +1)
D_3h (C_3v + σ_h):            χ FORBIDDEN  (σ_h alone forbids)
```

### `pyrochlore_tetra_complete_catalog`

Full bilinear + trilinear exchange-tensor catalog for the pyrochlore
"up tetrahedron": 6 NN bonds × {DMI, J^s} × 3 site syms + 4 triangle
faces × χ × 3 site syms = **72 group-theoretic verdicts in one program**:

```
O_h:  6 × (DMI 0-dim, J^s 3-dim)   χ allowed × 4
O:    6 × (DMI 1-dim, J^s 3-dim)   χ allowed × 4
T_d:  6 × (DMI 1-dim, J^s 3-dim)   χ allowed × 4
```

Under O_h the 3-dim J^s is exactly the Curnoe-Ross-Kao parametrisation
(Ross 2011 PRX 1, 021002) for pyrochlore quantum-spin-ice. Under O / T_d
the DMI becomes 1-dim per bond. The catalog is the parameter scaffold
a downstream DFT or micromagnetic simulator consumes.

## Reproducibility bound

Across macOS arm64 (Apple clang) and Linux x86_64 / aarch64 (gcc + clang),
the library's numerical output agrees to machine precision (documented
per-kernel in the CHANGELOG and METHODS.md). Example-level output
agrees within the tolerances above; Lanczos convergence is the loosest
— the 24-site ground-state energy may vary by a few parts in 10⁶
across runs with different vector-starting noise magnitudes.

Any deviation from the expected outputs above on a standard CPU build
indicates a regression.
