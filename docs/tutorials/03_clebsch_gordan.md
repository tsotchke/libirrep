# Tutorial 03 — Clebsch-Gordan coefficients

Clebsch-Gordan (CG) coefficients `⟨j₁ m₁; j₂ m₂ | J M⟩` are the unitary
change-of-basis matrix elements from the uncoupled product basis
`|j₁, m₁⟩ ⊗ |j₂, m₂⟩` to the coupled basis `|(j₁ j₂) J, M⟩`. They are the
central computational object in the coupling of angular momenta:

```
|(j₁ j₂) J, M⟩ = Σ_{m₁, m₂} ⟨j₁ m₁; j₂ m₂ | J M⟩ · |j₁, m₁⟩ ⊗ |j₂, m₂⟩.
```

Equivalently, CGs realise the Clebsch-Gordan decomposition of a tensor
product of SU(2) irreps:

```
j₁ ⊗ j₂ = |j₁ − j₂| ⊕ (|j₁ − j₂| + 1) ⊕ … ⊕ (j₁ + j₂).
```

Each allowed `J` on the right is present with multiplicity one — this is
the property that characterises the simply-reducible structure of SU(2).
The complete derivation and selection-rule analysis appears in
Varshalovich, Moskalev & Khersonskii §8 and Sakurai & Napolitano §3.8;
this tutorial assumes the definitions and walks through libirrep's
concrete API.

## 1. Selection rules

`⟨j₁ m₁; j₂ m₂ | J M⟩ = 0` unless **all** four conditions hold:

1. **Triangle:** `|j₁ − j₂| ≤ J ≤ j₁ + j₂`.
2. **z-component sum:** `m₁ + m₂ = M`.
3. **Ranges:** `|m₁| ≤ j₁`, `|m₂| ≤ j₂`, `|M| ≤ J`.
4. **Integer sum:** `2(j₁ + j₂ + J)` is an even integer (so that the
 Racah sum is over integer `k`).

Libirrep returns **exactly `0.0`** in every violation case. Callers can
sum CGs over index ranges freely without pre-filtering — forbidden terms
are already zero, so the sum is numerically correct.

## 2. The library surface

Two signatures for integer `j`:

```c
#include <irrep/clebsch_gordan.h>

/* ⟨1 0; 1 0 | 2 0⟩ = √(2/3) */
double cg = irrep_cg(1, 0, 1, 0, 2, 0);

/* Singlet projector ⟨1 1; 1 −1 | 0 0⟩ = 1/√3 */
double s = irrep_cg(1, 1, 1, -1, 0, 0);

/* Selection-rule violation (m-sum mismatch) — returns 0.0 exactly. */
double z = irrep_cg(1, 1, 1, 1, 2, 0);
```

Half-integer `j` via the **doubled-integer convention**: every
angular-momentum argument is passed as `2 · j`, so spin-½ enters as
`two_j = 1`:

```c
/* ⟨½ ½; ½ −½ | 1 0⟩ = 1/√2 */
double c = irrep_cg_2j(1, 1, 1, -1, 2, 0);
```

The `_2j` suffix is consistent across the library: `irrep_cg_2j`,
`irrep_wigner_d_small_2j`, `irrep_wigner_6j_2j`, etc.

## 3. Closed-form values for small j

These are worth committing to memory. libirrep's
`tests/reference_data/cg_reference.h` carries them as
`(sign, numerator, denominator)` rationals so the test compares
bit-against `sign · sqrt(num / den)` computed at runtime by libm.

**(½, ½) coupling** (spin-½ ⊗ spin-½ = singlet ⊕ triplet):

| ⟨j₁ m₁; j₂ m₂ \| J M⟩ | Value |
| --------------------------------- | -------: |
| ⟨½ ½; ½ ½ \| 1 1⟩ | 1 |
| ⟨½ ½; ½ −½ \| 1 0⟩ | 1/√2 |
| ⟨½ −½; ½ ½ \| 1 0⟩ | 1/√2 |
| ⟨½ −½; ½ −½ \| 1 −1⟩ | 1 |
| ⟨½ ½; ½ −½ \| 0 0⟩ | +1/√2 |
| ⟨½ −½; ½ ½ \| 0 0⟩ | −1/√2 |

The last two exhibit the antisymmetry of the spin singlet under particle
exchange.

**(1, 1) coupling** (vector ⊗ vector = scalar ⊕ pseudovector ⊕ traceless
symmetric tensor):

| ⟨1 m₁; 1 m₂ \| J M⟩ | J = 2 | J = 1 | J = 0 |
| ---------------------- | --------------- | -------------- | ------------- |
| m₁ = 1, m₂ = 1 | 1 | — | — |
| m₁ = 1, m₂ = 0 | 1/√2 | +1/√2 | — |
| m₁ = 1, m₂ = −1 | 1/√6 | +1/√2 | +1/√3 |
| m₁ = 0, m₂ = 0 | √(2/3) | 0 | −1/√3 |
| m₁ = −1, m₂ = 1 | 1/√6 | −1/√2 | +1/√3 |

Dashes indicate selection-rule zeros. The `J = 0` column is the
projector onto the scalar part of a vector tensor product — the formula
that libirrep's `irrep_tp_apply_uvw` uses, up to a `(1/√3)` prefactor,
to recover the cartesian dot product of two polar vectors (see
`examples/torque_net_tp_paths.c` and
[`tutorials/05_tensor_products.md`](05_tensor_products.md) §7).

## 4. Wigner 3j and the CG ↔ 3j relation

The Wigner 3j symbol is defined so that all three columns play
symmetric roles under permutation:

```
(j₁ j₂ j₃ ; m₁ m₂ m₃) =
 (−1)^{j₁ − j₂ − m₃} / √(2 j₃ + 1) · ⟨j₁ m₁; j₂ m₂ | j₃ (−m₃)⟩.
```

Column permutations produce at most a `(−1)^{j₁+j₂+j₃}` phase; the
explicit tables are in Varshalovich §8.5 (columns 1↔2, 1↔3) and §8.6
(cyclic). Libirrep provides both CG and 3j:

```c
double cg_val = irrep_cg(1, 0, 1, 0, 2, 0); /* √(2/3) */
double three_j = irrep_wigner_3j(1, 0, 1, 0, 2, 0); /* = √(2/3) · (-1)^{1-1+0} / √5 */
```

The `tests/test_clebsch_gordan.c` suite cross-checks the CG ↔ 3j
relation for `(j₁, j₂) ∈ {(½, ½), (1, ½), (1, 1), (3/2, 1), (2, 2),
(3, 2)}` to machine precision.

## 5. Orthogonality relations

Two orthogonality relations, both exact:

```
Σ_{m₁, m₂} ⟨j₁ m₁; j₂ m₂ | J M⟩ · ⟨j₁ m₁; j₂ m₂ | J' M'⟩ = δ_{J, J'} δ_{M, M'},
Σ_{J, M} ⟨j₁ m₁; j₂ m₂ | J M⟩ · ⟨j₁ m₁'; j₂ m₂' | J M⟩ = δ_{m₁, m₁'} δ_{m₂, m₂'}.
```

(Varshalovich §8.2). Libirrep's test suite verifies the first to
`1e-10` summing over `(m₁, m₂)` at fixed `(j₁, j₂, J)` for `j ≤ 20`.
Numerically:

```c
/* Σ_M ⟨j1 m1; j2 m2 | J M⟩² = 1 for any valid (j1, m1, j2, m2, J). */
double s = 0.0;
for (int M = -2; M <= 2; ++M) {
 double c = irrep_cg(1, 0, 1, 0, 2, M);
 s += c * c;
}
/* s ≈ 1.0 — forbidden M values contribute exact zero. */
```

## 6. Cached tables

For workloads that evaluate many CGs in the same `(j₁, j₂)` range —
tensor-product descriptor builds being the canonical example — the
`cg_table_t` cache amortises the per-call cost:

```c
cg_table_t *t = irrep_cg_table_build(/*j1_max=*/3, /*j2_max=*/3);
double v = irrep_cg_lookup(t, 2, 1, 1, -1, 2, 0); /* O(1) */
irrep_cg_table_free(t);
```

Table storage is `O(j₁_max · j₂_max · (j₁_max + j₂_max))` doubles, laid
out flat with explicit strides. Selection-rule violations and out-of-
range lookups both return `0.0`; the caller must free the table when
done.

## 7. Schulten–Gordon recurrence (implementation note)

The canonical closed-form expression for CG is Racah's single-sum
(1942, *Phys. Rev.* **62**, 438):

```
⟨j₁ m₁; j₂ m₂ | J M⟩ =
 δ_{M, m₁+m₂} · √(2J+1) · Δ(j₁, j₂, J) · √{N(j₁, m₁, j₂, m₂, J, M)}
 · Σ_k (−1)^k / [ k! · (j₁+j₂−J−k)! · (j₁−m₁−k)! · (j₂+m₂−k)!
 · (J−j₂+m₁+k)! · (J−j₁−m₂+k)! ].
```

Evaluated directly (or in log-gamma form), the alternating-sign sum
is cancellation-limited: individual terms can exceed the final result
by many orders of magnitude, and the subtraction leaves only the
error bits. In libirrep we measured this regime to produce NaN past
`j ≈ 60` near the triangle edge and `~10⁻⁹` sum-rule drift at `j = 50`.

The production implementation therefore computes CG indirectly via
the Wigner 3j symbol, using the Schulten–Gordon backward three-term
recurrence in `j₁` (Luscombe & Luban, *Phys. Rev. E* **57**, 7274,
1998):

```
j · E(j+1) · T(j+1) + F(j) · T(j) + (j+1) · E(j) · T(j−1) = 0
```

with `E(j)² = [j² − (j₂−j₃)²] · [(j₂+j₃+1)² − j²] · [j² − m₁²]` and
`F(j) = −(2j+1) · { [j₂(j₂+1) − j₃(j₃+1)] · m₁ + j(j+1) · (m₂ − m₃) }`.
The full `T(j)` series for `j ∈ [j_min, j_max]` is normalised by the
sum rule `Σ_j (2j+1) · T(j)² = 1` and sign-anchored at `j_max` by
`(−1)^{j₂−j₃−m₁}`. Clebsch-Gordan then derives from `T(j₁)` via
`⟨j₁ m₁; j₂ m₂ | J M⟩ = (−1)^{j₁−j₂+M} · √(2J+1) · T_{m₁, m₂, −M}(j₁)`.

Measured sum-rule precision at `(m₁, m₂) = (j/2, −j/2)`:

```
  j = 20 : 1e-16    j = 50 : 0
  j = 30 : 2e-16    j = 80 : 6e-4   (non-classical leakage, j > 80)
```

Backward-only recurrence picks up the minimal (physical) solution
cleanly through the classically allowed range; the `6 × 10⁻⁴` residual
at `j = 80` is subdominant-solution contamination in the deep
non-classical region. Miller two-directional iteration would recover
precision past `j ≈ 80`; tracked in `TODO.md`.

## 8. Dense blocks for kernel code

Tensor-product kernels want every `(m₁, m₂, M)` combination at a given
`(j₁, j₂, J)` triple flat in memory:

```c
int j1 = 1, j2 = 1, J = 2;
int size = (2*j1 + 1) * (2*j2 + 1) * (2*J + 1); /* = 45 */
double block[45];
irrep_cg_block(j1, j2, J, block);
/* index: [(m1 + j1) · (2j2+1) + (m2 + j2)] · (2J+1) + (M + J) */
```

Selection-rule violations sit in the block as **exact** zeros, so a
naïve sum over the whole block is numerically indistinguishable from a
pre-filtered sum.

## 9. References

- **Racah, G.** "Theory of Complex Spectra. II," *Phys. Rev.* **62**,
 438 (1942). DOI: [10.1103/PhysRev.62.438](https://doi.org/10.1103/PhysRev.62.438).
- **Sakurai, J. J. & Napolitano, J.** *Modern Quantum Mechanics* (3rd
 ed., Cambridge, 2020), §3.8 and Appendix A.
- **Varshalovich, D. A., Moskalev, A. N., & Khersonskii, V. K.**
 *Quantum Theory of Angular Momentum* (World Scientific, 1988), §8.
- **Edmonds, A. R.** *Angular Momentum in Quantum Mechanics* (Princeton,
 1957), §3–§4.

Full bibliography in [`REFERENCES.md`](../REFERENCES.md).
