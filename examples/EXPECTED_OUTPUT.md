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

## Reproducibility bound

Across macOS arm64 (Apple clang) and Linux x86_64 / aarch64 (gcc + clang),
the library's numerical output agrees to machine precision (documented
per-kernel in the CHANGELOG and METHODS.md). Example-level output
agrees within the tolerances above; Lanczos convergence is the loosest
— the 24-site ground-state energy may vary by a few parts in 10⁶
across runs with different vector-starting noise magnitudes.

Any deviation from the expected outputs above on a standard CPU build
indicates a regression.
