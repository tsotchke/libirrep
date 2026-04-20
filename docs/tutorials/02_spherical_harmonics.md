# Tutorial 02 — Spherical harmonics

Spherical harmonics `Y_l^m(θ, φ)` are the canonical orthonormal basis
on the 2-sphere. They arise as the angular part of separation of
variables for the Laplacian in spherical coordinates (Jackson §3.5),
as the matrix elements of the rotation operator in position space, and
— the view most relevant to libirrep — as a representation-theoretic
basis for the `(l)` irrep of SO(3) realised on `L²(S²)`.

This tutorial covers the three equivalent API surfaces libirrep
exposes, the conventions fixing signs and phases, and the numerical
considerations that matter at `l` above a few.

## 1. Three API surfaces

Libirrep provides SH in three equivalent forms:

| Function | Output | When to use |
| ------------------------------------- | ------- | ------------------------------------------ |
| `irrep_sph_harm(l, m, θ, φ)` | complex | Matching textbook formulas |
| `irrep_sph_harm_real(l, m, θ, φ)` | real | Real-basis work with polar inputs |
| `irrep_sph_harm_cart(l, out, r̂)` | real | GNN / MD hot path — cartesian unit vectors |

All three return the same mathematical object up to the complex↔real
basis change (§3). The cartesian form is the hot path for every
NequIP / MACE / Allegro message-passing layer: edge vectors arrive in
cartesian, and converting to `(θ, φ)` would cost a transcendental plus
introduce a pole singularity that the cartesian form avoids.

## 2. The Condon-Shortley phase convention

Libirrep uses the standard Condon-Shortley phase
(Sakurai & Napolitano 2020, Appendix A; Jackson 3rd ed., §3.6):

```
Y_l^m(θ, φ) = √{(2l+1)/(4π) · (l−m)!/(l+m)!} · P_l^m(cos θ) · e^{i m φ},
```

with the associated Legendre polynomial

```
P_l^m(x) = (−1)^m · (1 − x²)^{m/2} · (d^m / dx^m) P_l(x).
```

The `(−1)^m` factor lives in `P_l^m` — applied once, in
`irrep_legendre_assoc`. The Wigner-D construction does **not**
reintroduce it. An immediate consequence of this choice is the
relation `Y_l^{−m}(θ, φ) = (−1)^m · [Y_l^m(θ, φ)]*` without any
additional phase work.

The orthonormality relation

```
∫_{S²} Y_l^m(r̂) · [Y_{l'}^{m'}(r̂)]* dΩ = δ_{l, l'} δ_{m, m'}
```

holds on the unit sphere with `dΩ = sin θ dθ dφ`; libirrep verifies it
numerically via Lebedev-Laikov quadrature of order 35 (exact for
polynomials of total degree ≤ 35 on `S²`) to an absolute residual below
`10⁻¹⁰`.

## 3. Real spherical harmonics (e3nn convention)

Real-basis SH are built from the complex ones via a block-diagonal
unitary transform:

```
Y_{l, 0}^{real} = Y_l^0,
Y_{l, +|m|}^{real} = √2 · ℜ(Y_l^{|m|}) = (1/√2) (Y_l^{+|m|} + (−1)^m Y_l^{−|m|}),
Y_{l, −|m|}^{real} = √2 · ℑ(Y_l^{|m|}) = (i/√2) (Y_l^{+|m|} − (−1)^m Y_l^{−|m|}).
```

The phase is fixed so that **`Y_{1, +1}^{real} ∝ +x`** at the equator —
the e3nn / Wikipedia convention, also used by Geiger et al. (2022)
`e3nn` and by libirrep. Some quantum-chemistry packages (ORCA, NWChem)
flip the odd-m components' sign; the change-of-basis is then a
diagonal `±1` matrix that the consumer must apply explicitly.

The `irrep_sph_harm_complex_to_real` function returns the explicit
`(2l+1) × (2l+1)` complex change-of-basis matrix `U`:

```c
double _Complex U[(2*l+1) * (2*l+1)];
irrep_sph_harm_complex_to_real(l, U);
/* Y_real = U · Y_complex. U is sparse — only ±1, ±1/√2, ±i/√2 entries. */
```

## 4. The cartesian form

For a unit vector `r̂ ∈ S²`, `irrep_sph_harm_cart(l, out, r_hat)` writes
`2l + 1` real values corresponding to `m = −l, −l+1, …, +l`:

```c
#include <irrep/spherical_harmonics.h>

double r_hat[3] = {0.0, 0.0, 1.0}; /* north pole */
double y1[3]; /* 2l + 1 = 3 for l = 1 */
irrep_sph_harm_cart(1, y1, r_hat);
/*
 * y1[0] = Y_{1, -1}^{real}(r̂)
 * y1[1] = Y_{1, 0}^{real}(r̂)
 * y1[2] = Y_{1, +1}^{real}(r̂)
 *
 * At the pole only m = 0 is non-zero; the ±m harmonics carry a
 * (sin θ)^{|m|} prefactor that vanishes there.
 */
```

**Layout invariant:** within a single `l`, `out[m + l]` holds the
`m`-indexed entry, so `m = −l` is at offset 0 and `m = +l` is at
offset `2l`.

For the whole triangle `l = 0, 1, …, l_max`, call `_cart_all`:

```c
int l_max = 4;
int total = (l_max + 1) * (l_max + 1); /* = 25 */
double out[25];
irrep_sph_harm_cart_all(l_max, out, r_hat);

/* Layout:
 * out[0] = Y_{0, 0}
 * out[1..3] = Y_{1, -1..+1}
 * out[4..8] = Y_{2, -2..+2}
 * out[9..15] = Y_{3, -3..+3}
 * out[16..24] = Y_{4, -4..+4}
 */
```

For many unit vectors at once — the NequIP / MACE message-passing
regime — the `_batch` variant dispatches through the runtime SIMD
pointer (NEON on aarch64, scalar elsewhere):

```c
int l_max = 4;
int block = (l_max + 1) * (l_max + 1);
size_t N = 4096; /* many edges */
double *r_hats = malloc(N * 3 * sizeof(double));
double *out = malloc(N * block * sizeof(double));
/* ... fill r_hats with unit vectors ... */
irrep_sph_harm_cart_all_batch(l_max, N, r_hats, out);
```

Bit-exact against the per-edge scalar call across
`l_max ∈ {0, …, 6}` and every pole / near-pole edge case
(`tests/test_spherical_harmonics.c`). The 1.2 NEON kernel with
triangular-Legendre-grid caching delivers ~1.9× the pure-scalar loop
throughput.

## 5. Gradient

The force path in equivariant MD requires `∇_{r̂} Y_l^m`. Libirrep's
`irrep_sph_harm_cart_grad` writes `3 × (2l+1)` values with
axis-major layout:

```c
int l = 2;
int d = 2 * l + 1; /* = 5 */
double grad[3 * 5];
irrep_sph_harm_cart_grad(l, grad, r_hat);
/*
 * grad[axis * d + (m + l)] = ∂Y_{l, m}/∂r̂[axis] for axis ∈ {0:x, 1:y, 2:z}.
 */
```

A critical invariant: **the SH gradient is tangent to `S²`**. That is,

```
r̂ · ∇_{r̂} Y_l^m(r̂) = 0,
```

since `Y_l^m` is a function on the unit sphere and its radial derivative
is identically zero. Libirrep's `tests/test_spherical_harmonics.c`
verifies this to `10⁻⁸` for `l ≤ 3`. The tangency is what makes the
edge-geometry chain rule in `irrep_nequip_layer_apply_forces` collapse
to a single `(1/r) · ∂Y/∂r̂[axis]` term rather than the full
Jacobian-of-normalisation form.

A batched gradient variant `irrep_sph_harm_cart_all_grad_batch` produces
the gradient across the full `(l, m)` triangle for many edges at once,
laid out as `out[(edge · 3 + axis) · (l_max+1)² + (l² + m + l)]`. This
is the entry point used by `nequip_layer_apply_forces`.

## 6. Orthonormality verification (worked example)

The following reproduces the orthonormality check in
`tests/test_spherical_harmonics.c` using a Lebedev rule:

```c
#include <irrep/quadrature.h>
#include <irrep/spherical_harmonics.h>

int deg = 8; /* exact to l ≤ 8 on total degree */
int N = irrep_quadrature_sphere_size(deg);
double *nodes = malloc(N * 4 * sizeof(double)); /* [x, y, z, w] */
irrep_quadrature_sphere_fill(deg, nodes);

double integral = 0.0;
for (int i = 0; i < N; ++i) {
 double rhat[3] = { nodes[4*i + 0], nodes[4*i + 1], nodes[4*i + 2] };
 double w = nodes[4*i + 3];
 double buf[5]; /* 2·2+1 for l=2 */
 irrep_sph_harm_cart(2, buf, rhat);
 integral += w * buf[0] * buf[0]; /* ⟨Y_{2, -2} | Y_{2, -2}⟩ */
}
/* Note: irrep_quadrature_sphere_fill returns weights normalised so that
 * Σ w = 1 (not 4π). To compare with the textbook orthonormality
 * ∫ Y² dΩ = 1 multiply by 4π:
 *
 *     assert(fabs(integral * 4 * M_PI - 1.0) < 1e-10);
 *
 * The test suite in tests/test_spherical_harmonics.c uses a raw
 * Gauss–Legendre × uniform-φ product with un-normalised weights and
 * therefore compares directly to 1.0. */
```

The product rule's polynomial exactness degree (`deg`) determines the
smallest SH orthonormality that can be verified: a rule exact for
polynomials of total degree `2 · l_max` verifies orthonormality through
`l_max`. For pure Lebedev rules, libirrep currently ships orders
**3, 5, 7** (exact through `l = 1, 2, 3`); higher orders remain on the
roadmap and callers needing broader exactness should use
`irrep_quadrature_sphere_fill` (Gauss–Legendre × uniform-φ), which
gives arbitrary exactness at ≈ 2× the point count.

## 7. The addition theorem

A useful identity linking the SH at two different arguments:

```
P_l(cos γ) = (4π / (2l + 1)) · Σ_m Y_l^m(r̂₁) · [Y_l^m(r̂₂)]*,
```

where `cos γ = r̂₁ · r̂₂` (Jackson §3.6, Arfken §15.6). The real-basis
analogue uses `Σ_m Y_{l, m}^{real}(r̂₁) Y_{l, m}^{real}(r̂₂)` on the
right (no complex conjugate — real SH are self-conjugate). Libirrep's
test suite asserts the addition theorem against the Legendre polynomial
`P_l` at randomly sampled `(r̂₁, r̂₂)` pairs to absolute residual
below `5 × 10⁻¹⁴` for `l ≤ 8`.

## 8. Numerical considerations at the pole

At the pole (`r̂ → ẑ`), the spherical-coordinate prefactor `sin θ`
vanishes. In the cartesian API, the derived quantity `r_xy = √(x² + y²)`
is the input to the `cos(mφ), sin(mφ)` recurrence via `x/r_xy, y/r_xy`;
near the pole this is a `0/0`. Libirrep switches to
`(cp, sp) = (1, 0)` when `r_xy < 10⁻¹⁴`, which is exactly the correct
limit since every `Y_l^{m≠0}` carries a `(sin θ)^{|m|}` factor that
vanishes at the pole. The tests include exact-pole inputs
`(r̂ = ±ẑ)` and near-pole inputs on both sides of the threshold.

## 9. Associated Legendre via stable recurrence

For `l` above a few, the direct formula
`P_l^m = (−1)^m (1−x²)^{m/2} d^m P_l / dx^m` is numerically ill-posed:
the derivative of the unsigned Legendre polynomial is
catastrophically ill-conditioned at `|x| → 1`. Libirrep instead uses
the stable forward recurrence (Press et al. 2007, §6.7;
Limpanuparb & Milthorpe 2014):

```
P_l^m(x) = [ (2l − 1) · x · P_{l−1}^m(x) − (l + m − 1) · P_{l−2}^m(x) ] / (l − m),
```

evaluated from the closed-form initialisers

```
P_m^m(x) = (−1)^m · (2m − 1)!! · (1 − x²)^{m/2},
P_{m+1}^m(x) = x · (2m + 1) · P_m^m(x).
```

This forward recurrence is numerically stable for `x ∈ [−1, 1]` and
preserves full double precision up to `l = 50` in our testing. The
NEON SH batch kernel pre-computes the full triangular grid
`P_l^m(z)` in one sweep at O(l²) work (1.2 optimisation; replaces the
O(l³) per-call restarts).

## 10. Single precision

Single-precision variants `irrep_sph_harm_cart_f32` and `_all_f32` are
available. They currently wrap the `double` kernel with a cast — a
proper single-precision SIMD path is scheduled for 1.3. The precision
loss is the usual 7 decimal digits of binary32; sufficient for most
inference pipelines, insufficient for training where gradient noise
compounds.

## 11. References

- **Sakurai, J. J. & Napolitano, J.** *Modern Quantum Mechanics* (3rd
 ed., Cambridge, 2020), Appendix A — Condon-Shortley phase,
 tabulated `Y_l^m` at low `l`.
- **Jackson, J. D.** *Classical Electrodynamics* (3rd ed., Wiley,
 1998), §3.5–§3.6 — orthonormality, addition theorem.
- **Arfken, G. B., Weber, H. J., & Harris, F. E.** *Mathematical
 Methods for Physicists* (7th ed., Academic Press, 2013), §15.
- **Limpanuparb, T. & Milthorpe, J.** "Associated Legendre Polynomials
 and Spherical Harmonics Computation for Chemistry Applications,"
 arXiv:1410.1748 (2014).
- **Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B.
 P.** *Numerical Recipes* (3rd ed., Cambridge, 2007), §6.7.
- **Lebedev, V. I. & Laikov, D. N.** "A quadrature formula for the
 sphere of the 131st algebraic order of accuracy," *Doklady
 Mathematics* **59**(3), 477 (1999).
- **Geiger, M. et al.** "e3nn: Euclidean Neural Networks,"
 arXiv:2207.09453 (2022) — real-SH sign convention.

Full bibliography in [`REFERENCES.md`](../REFERENCES.md).
