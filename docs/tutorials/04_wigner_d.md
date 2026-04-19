# Tutorial 04 — Wigner-D matrices

The Wigner-D matrix `D^j(R)` is the matrix element of the rotation
operator on the spin-`j` irrep, expressed in the `|j, m⟩` basis. It is
the fundamental tool for rotating irrep-valued features: if a quantity
`ψ_m` transforms as `|j, m⟩`, then under rotation `R ∈ SO(3)` it
transforms by

```
ψ_m → Σ_{m'} D^j_{m', m}(R) · ψ_{m'}.
```

Libirrep computes `D^j` via the Sakurai direct sum (Sakurai & Napolitano
2020, §3.8) evaluated in log-gamma form for numerical stability past
`j = 50`, and provides scalar, matrix, half-integer, and multiset-level
forms.

## 1. Factorisation via Euler ZYZ

Every `R ∈ SO(3)` factors as `R = R_z(α) · R_y(β) · R_z(γ)` with
`α ∈ [0, 2π)`, `β ∈ [0, π]`, `γ ∈ [0, 2π)` (physics convention —
Sakurai §3.3, Varshalovich §1.4). The Wigner-D factors correspondingly:

```
D^j_{m', m}(α, β, γ) = e^{−i m' α} · d^j_{m', m}(β) · e^{−i m γ},
```

where `d^j(β)` — the **small Wigner-d matrix** — is real and depends
only on `β`. The library implements this factorisation directly.

## 2. Evaluating a single matrix element

```c
#include <irrep/wigner_d.h>

/* Real small-d: d^2_{1, 0}(0.7). */
double d = irrep_wigner_d_small(/*j=*/2, /*mp=*/1, /*m=*/0, /*beta=*/0.7);

/* Full complex D: D^2_{1, 0}(α=0.4, β=0.7, γ=1.1). */
double _Complex D = irrep_wigner_D(2, 1, 0, 0.4, 0.7, 1.1);
```

Half-integer `j` enters through the `_2j` suffix:

```c
/* d^{1/2}_{1/2, −1/2}(β) at β = 0.5. */
double d_half = irrep_wigner_d_small_2j(/*two_j=*/1, /*two_mp=*/1,
 /*two_m=*/ -1, /*beta=*/0.5);
```

## 3. Full matrices

For a fixed `j`, the `(2j+1) × (2j+1)` matrix is built by a single call:

```c
int j = 2, d = 2*j + 1; /* 5 */
double d_small[5 * 5];
double _Complex D_full[5 * 5];

irrep_wigner_d_matrix(j, d_small, /*beta=*/0.7);
irrep_wigner_D_matrix(j, D_full, /*α=*/0.4, /*β=*/0.7, /*γ=*/1.1);

/* Row index = m', column index = m, row-major layout. */
```

## 4. Block-diagonal D on an irrep multiset

The most-used API in equivariant networks: lift `D(R)` from a single
irrep to the full direct-sum feature space specified by an
`irrep_multiset_t`:

```c
#include <irrep/multiset.h>
#include <irrep/wigner_d.h>

irrep_multiset_t *h = irrep_multiset_parse("2x0e + 1x1o + 1x2e");
int n = h->total_dim; /* 2·1 + 1·3 + 1·5 = 10 */

double _Complex *D = malloc(n * n * sizeof(*D));
irrep_wigner_D_multiset(h, D, /*α=*/0.4, /*β=*/0.7, /*γ=*/1.1);
```

The result is block-diagonal: two `1×1` blocks (two copies of `l=0`),
one `3×3` block (`l=1`), one `5×5` block (`l=2`), each filled with the
corresponding `D^l(R)`. Multiplying a feature vector by this matrix
rotates every irrep component consistently — the foundational operation
behind every equivariance test in the library and every real-world
NequIP / MACE forward pass.

## 5. The Sakurai direct-sum formula

The closed-form expression for the small-d matrix is (Sakurai Eq. 3.8.33;
Varshalovich §4.3.4 Eq. 10):

```
d^j_{m', m}(β) =
 Σ_k (−1)^{k + m' − m} · √{(j+m)! (j−m)! (j+m')! (j−m')!}
 / [ k! · (j+m'−k)! · (j−m−k)! · (k+m−m')! ]
 · (cos β/2)^{2j + m − m' − 2k} · (sin β/2)^{2k + m' − m}.
```

The summation range is the intersection of the four non-negativity
constraints on the factorial arguments. For direct evaluation, the
`(2j)!`, `(j ± m)!`, and denominator factorials overflow at modest `j`;
the library refactors the expression to compute everything in log-gamma
space via libm's `lgamma`. This is the same numerical technique as the
Clebsch-Gordan implementation and remains stable past `j = 50` in
double precision.

## 6. Closed-form small-d values at low j

Useful for sanity checks and hand derivations. These are all tabulated
in Varshalovich §4.16 and reproduced to machine precision by
`irrep_wigner_d_small`:

**j = ½:**

```
d^{1/2}(β) = [ cos β/2 −sin β/2 ]
 [ sin β/2 cos β/2 ]
```

**j = 1:**

```
d^1(β) = [ cos²(β/2) −sin β / √2 sin²(β/2) ]
 [ sin β / √2 cos β −sin β / √2 ]
 [ sin²(β/2) sin β / √2 cos²(β/2) ]
```

The `d^j_{j, j}(β) = cos^{2j}(β/2)` identity that you can read off the
top-left entry generalises to every `j` and is a useful invariant to
check when debugging a custom rotation implementation.

## 7. Unitarity and composition

`D^j(R)` is unitary: `D^j(R) · D^j(R)† = I` up to `10⁻¹³` in double
precision. Composition satisfies `D^j(R_1 R_2) = D^j(R_1) · D^j(R_2)`.
Both identities are verified at random SO(3) pairs for `j ∈ {1, 2, 3,
4, 6}` in `tests/test_wigner_d.c`; the test hits every Wigner-D entry
at each of 1065 test pairs and asserts absolute error below `10⁻⁸`.

A quick hand check at `j = 1`:

```c
irrep_rot_matrix_t R1 = /* ... some rotation ... */;
irrep_rot_matrix_t R2 = /* ... another ... */;
irrep_rot_matrix_t R_composed = irrep_rot_compose(R1, R2);

double _Complex D1[9], D2[9], D_prod[9], D_composed[9];
irrep_euler_zyz_t e1 = irrep_euler_zyz_from_rot(R1);
irrep_euler_zyz_t e2 = irrep_euler_zyz_from_rot(R2);
irrep_euler_zyz_t e_composed = irrep_euler_zyz_from_rot(R_composed);

irrep_wigner_D_matrix(1, D1, e1.alpha, e1.beta, e1.gamma);
irrep_wigner_D_matrix(1, D2, e2.alpha, e2.beta, e2.gamma);
irrep_wigner_D_matrix(1, D_composed, e_composed.alpha, e_composed.beta, e_composed.gamma);

/* D_prod = D1 @ D2 (complex 3×3 matrix product) */
/* Assert max |D_prod - D_composed| < 1e-10 per entry */
```

## 8. Derivative `∂d/∂β`

Needed for gradient-based equivariant force evaluation (e.g., magnetic-
system LLG integrators, or NequIP's edge-geometry gradient). Libirrep
provides it via `irrep_wigner_d_small_dbeta`:

```c
double dd = irrep_wigner_d_small_dbeta(/*j=*/2, /*mp=*/1, /*m=*/0,
 /*beta=*/0.7);

/* Central-difference cross-check: */
double h = 1e-6;
double fd = (irrep_wigner_d_small(2, 1, 0, 0.7 + h)
 - irrep_wigner_d_small(2, 1, 0, 0.7 - h)) / (2 * h);
/* |dd - fd| < 1e-9 typically. */
```

The derivative is obtained analytically by term-wise differentiation of
the Sakurai sum — each factor `(cos β/2)^{p}` contributes
`−(p/2) (sin β/2) · (cos β/2)^{p−1}` and similarly for `(sin β/2)^{q}`.

## 9. Euler-angle vs. axis-angle inputs

If your rotation arrives as a quaternion or axis-angle, convert first:

```c
irrep_quaternion_t q = /* some unit quaternion */;
irrep_rot_matrix_t R = irrep_rot_from_quat(q);
irrep_euler_zyz_t e = irrep_euler_zyz_from_rot(R);

double _Complex D[9];
irrep_wigner_D_matrix(1, D, e.alpha, e.beta, e.gamma);
```

Near the gimbal-lock region (`|sin β| < 10⁻¹²` — see
[`PHYSICS_APPENDIX.md`](../PHYSICS_APPENDIX.md) §2.1), the Euler
extraction collapses `α + γ → α` and sets `γ = 0`. This preserves the
correct rotation (`D` is unchanged) but loses the ability to distinguish
`α` from `γ` if your pipeline depends on that — in which case, avoid
Euler entirely and work with quaternions through to `irrep_rot_apply`.

## 10. References

- **Sakurai, J. J. & Napolitano, J.** *Modern Quantum Mechanics* (3rd
 ed., Cambridge, 2020), §3.8 — the direct-sum formula (Eq. 3.8.33)
 and closed forms at small `j`.
- **Varshalovich, D. A., Moskalev, A. N., & Khersonskii, V. K.**
 *Quantum Theory of Angular Momentum* (World Scientific, 1988),
 §4.3.4 (Eq. 10) and §4.16 (tabulated closed forms).
- **Biedenharn, L. C. & Louck, J. D.** *Angular Momentum in Quantum
 Physics* (Addison-Wesley, 1981), §3 — representation-theoretic
 derivation of `D^j(R)`.
- **Wigner, E. P.** *Group Theory and its Application to the Quantum
 Mechanics of Atomic Spectra* (Academic Press, 1959).

Full bibliography in [`REFERENCES.md`](../REFERENCES.md).
