# Tutorial 01 — Rotations

Every non-trivial operation in libirrep starts from an SO(3) rotation
supplied by the caller. This tutorial covers the four equivalent
representations the library accepts, the pairwise conversions among
them, the group operations, and the numerical gotchas that routinely
bite practitioners.

## 1. Four representations

| Type | Shape | When preferable |
| --------------------- | ------------------ | --------------------------------------------- |
| `irrep_axis_angle_t` | unit axis + angle | Physical torques, geometric intuition |
| `irrep_quaternion_t` | `{x, y, z, w}` | Composition, SLERP, numerical robustness |
| `irrep_rot_matrix_t` | row-major 3×3 | Matrix-vector apply, equivariance checks |
| `irrep_euler_zyz_t` | `(α, β, γ)` radians | Interop with quantum chemistry, Wigner-D input |

These are four coordinate charts on the same manifold: SO(3), the
compact Lie group of orientation-preserving orthogonal transformations
of `ℝ³`. As a manifold SO(3) is diffeomorphic to `ℝP³` (Hall 2015, §1.2),
connected with fundamental group `ℤ₂`. Every pair of representations has
a bijection; libirrep implements each pairwise conversion and tests them
to round-trip to `10⁻¹²`.

Concretely:

```c
#include <irrep/so3.h>
#include <math.h>

irrep_axis_angle_t aa = { .axis = {0, 0, 1}, .angle = M_PI / 4 };

irrep_quaternion_t q = irrep_quat_from_axis_angle(aa);
irrep_rot_matrix_t R = irrep_rot_from_quat(q);
irrep_euler_zyz_t e = irrep_euler_zyz_from_rot(R);
```

Each conversion is a pure function; the library provides every pair,
not just the four above.

## 2. Conventions

These are enforced throughout libirrep:

- **Active rotations.** A rotation `R` sends a vector `v` to `R · v`,
 leaving the coordinate axes fixed. Passive rotations — where instead
 the axes rotate and the vector's components change by `R⁻¹` — appear
 in older aerospace and crystallography literature; translate by
 inverting.
- **Right-handed basis.** Positive rotation about `ẑ` sends `x̂ → ŷ`.
- **Euler ordering: ZYZ.** The physics convention of Sakurai §3.3 and
 Varshalovich §1.4: `R(α, β, γ) = R_z(α) · R_y(β) · R_z(γ)`, with
 `α ∈ [0, 2π)`, `β ∈ [0, π]`, `γ ∈ [0, 2π)`. Read right to left — the
 inner factor `R_z(γ)` acts first on a vector being rotated.
- **Quaternion layout: scalar-last `{x, y, z, w}`.** The vector part
 occupies the first three `float64` lanes, which matches Eigen, glTF
 2.0, and most SIMD-ergonomic conventions. A unit quaternion
 `q = sin(θ/2) · n̂ + cos(θ/2) · ê_scalar` represents a rotation by
 angle `θ` about unit axis `n̂`. Both `q` and `−q` represent the same
 SO(3) rotation (double cover); Shoemake's uniform sampler
 canonicalises to `w ≥ 0`.

Formal derivations of each convention are in
[`PHYSICS_APPENDIX.md`](../PHYSICS_APPENDIX.md) §§2, 4; primary sources
in [`REFERENCES.md`](../REFERENCES.md).

## 3. Composition, inverse, and apply

```c
irrep_rot_matrix_t R1 = irrep_rot_from_axis_angle((irrep_axis_angle_t){{0,0,1}, 0.1});
irrep_rot_matrix_t R2 = irrep_rot_from_axis_angle((irrep_axis_angle_t){{1,0,0}, 0.2});
irrep_rot_matrix_t R = irrep_rot_compose(R1, R2); /* first R1, then R2 */

double v[3] = {1, 0, 0};
double w[3];
irrep_rot_apply(R, v, w); /* w = R · v */

/* Invariant: R · R⁻¹ = I to machine precision. */
irrep_rot_matrix_t I = irrep_rot_compose(R, irrep_rot_inverse(R));
```

The `irrep_rot_compose(a, b)` convention is "a then b" — read
left-to-right as a pipeline, matching the `R_z(α) · R_y(β) · R_z(γ)`
Euler expansion's matrix-multiplication order. Hamilton's quaternion
composition `q1 · q2` follows the same convention (`q1` applied first):

```c
irrep_quaternion_t q_composed = irrep_quat_compose(q1, q2);
```

## 4. Exponential and logarithm on SO(3)

The Lie-theoretic exponential `exp : so(3) → SO(3)` sends a
skew-symmetric generator `ω̂` to the rotation by angle `|ω|` about axis
`ω / |ω|`. Libirrep implements it via Rodrigues' rotation formula:

```
exp([ω]_×) = I + sin(|ω|) · [ω̂]_× + (1 − cos(|ω|)) · [ω̂]_×²
```

with a first-order Taylor fallback `I + [ω]_×` for `|ω| < 10⁻¹²`, since
`sin(θ)/θ` and `(1−cos θ)/θ²` both cancel catastrophically in that
regime.

```c
double omega[3] = {0.3, 0.1, -0.2};
irrep_rot_matrix_t R = irrep_rot_exp(omega);

double back[3];
irrep_rot_log(R, back);
/* back ≈ omega up to 1e-12 in each component. */
```

The logarithm `irrep_rot_log` is the local inverse of `exp`. The naïve
closed form

```
ω = (θ / (2 sin θ)) · axial(R − Rᵀ), θ = arccos((Tr R − 1) / 2)
```

breaks down at `θ = π` (rotations by a half-turn), where `R = Rᵀ` so
`axial(R − Rᵀ) = 0` — the axis has to be recovered from the symmetric
part `R + Rᵀ`. Libirrep routes through the quaternion representation
and uses Markley's branch-switching algorithm (Markley 2008) to select
the correct axis for every input, preserving accuracy to `10⁻¹²`
throughout.

## 5. Interpolation and distance

Two rotations are joined by a geodesic in SO(3) whose length — the
geodesic distance — is the angle of `R_a⁻¹ R_b`. The quaternion SLERP
(Shoemake 1985) parameterises this geodesic:

```c
irrep_quaternion_t q_start = irrep_quat_from_axis_angle((irrep_axis_angle_t){{0,0,1}, 0.0});
irrep_quaternion_t q_end = irrep_quat_from_axis_angle((irrep_axis_angle_t){{0,0,1}, 1.0});

irrep_quaternion_t q_mid = irrep_quat_slerp(q_start, q_end, 0.5);
/* q_mid is the midpoint of the geodesic — a rotation by 0.5 rad. */

double theta = irrep_rot_geodesic_distance(
 irrep_rot_from_quat(q_start), irrep_rot_from_quat(q_end));
/* theta = 1.0 (exactly). */
```

SLERP takes the short way around the sphere; libirrep's implementation
canonicalises the signs of the two quaternion endpoints so the
interpolation is always through the shorter geodesic.

## 6. Karcher-Fréchet mean

The weighted mean of `n` unit quaternions on `S³` — equivalently, on
SO(3) via the double cover — is the Riemannian centre of mass defined
by

```
μ = argmin_μ Σ_i w_i · d_SO(3)(μ, q_i)².
```

There is no closed form; Karcher's iteration (Karcher 1977; Buss &
Fillmore 2001) converges quadratically from a chordal initialisation:

```
μ_{k+1} = μ_k · exp( (1 / Σ w_i) · Σ_i w_i · log(μ_k⁻¹ · q_i) ).
```

Libirrep's `irrep_quat_frechet_mean` runs this up to 64 iterations with
a squared-step tolerance of `10⁻²⁴` (so the converged mean is accurate
to `10⁻¹²` rad, below the angular resolution of a unit quaternion in
double precision):

```c
irrep_quaternion_t qs[3] = { q_a, q_b, q_c };
double ws[3] = { 1.0, 1.0, 2.0 };
irrep_quaternion_t mu = irrep_quat_frechet_mean(qs, ws, 3);
```

The implementation flips each `q_i` to the same hemisphere as the
current `μ` before taking the log, so the mean is well-defined even
when the input quaternions come from sign-agnostic sources.

## 7. Uniform sampling on SO(3)

Shoemake's algorithm (Shoemake 1992) generates quaternions uniformly
on `S³`; modding by the `{q, −q}` equivalence gives a uniform
distribution on SO(3). Libirrep owns the RNG state explicitly so that
concurrent samplers can run without contention:

```c
uint64_t rng = 0x123456789abcdef0ULL;
for (int i = 0; i < 1000; ++i) {
 irrep_quaternion_t q = irrep_quat_random(&rng);
 /* q is uniformly distributed on SO(3); w ≥ 0 by convention. */
}
```

## 8. Numerical gotchas

These routinely produce wrong answers if ignored:

- **Gimbal lock in ZYZ.** `irrep_euler_zyz_from_rot` collapses `α` and
 `γ` into `α` when `|sin β| < 10⁻¹²`, leaving `γ = 0`. The recovered
 rotation is still exactly right (the Wigner-D factorisation is
 unchanged), but if your pipeline depends on distinguishing `α` from
 `γ` the information is gone. Avoid Euler for general code — use
 quaternions through to `irrep_rot_apply`.
- **Logarithm at `θ = π`.** The naïve `acos((Tr R − 1)/2)` formula
 loses precision as `Tr R → −1` — not because of a removable `sin`
 singularity, but because the axis extraction from the antisymmetric
 part vanishes. Libirrep's `irrep_rot_log` routes through the
 quaternion + Markley branch-switching; if you roll your own, make
 sure you do the same.
- **Near-zero `rot_exp`.** Rodrigues' full formula divides
 `(1 − cos θ)` and `sin θ` by `θ`, both cancelling at `θ → 0`.
 Libirrep's Taylor fallback at `|ω| < 10⁻¹²` is correct to
 `O(|ω|²) < 10⁻²⁴`, well below double-precision round-off on the
 linear term.
- **Non-unit quaternions.** The library's rotation API assumes
 unit-norm inputs. Non-unit quaternions will not produce rotations —
 `irrep_quat_normalize` projects to the unit sphere if you need it.
- **Scalar position.** Libirrep uses `{x, y, z, w}` with scalar `w`
 last. Some libraries (e.g. Boost, JPL astronomy code) use
 `{w, x, y, z}` — translate at the boundary.

## 9. References

- **Sakurai, J. J. & Napolitano, J.** *Modern Quantum Mechanics* (3rd
 ed., Cambridge, 2020), §3.3 — Euler ZYZ factorisation.
- **Hall, B. C.** *Lie Groups, Lie Algebras, and Representations* (2nd
 ed., Springer, 2015), §1.2 — SO(3) manifold structure and
 fundamental group.
- **Shuster, M. D.** "A Survey of Attitude Representations,"
 *J. Astronautical Sciences* **41**(4), 439–517 (1993).
- **Markley, F. L.** "Unit Quaternion from Rotation Matrix,"
 *J. Guidance, Control, and Dynamics* **31**(2), 440–442 (2008).
 DOI: [10.2514/1.31730](https://doi.org/10.2514/1.31730).
- **Shepperd, S. W.** "Quaternion from Rotation Matrix,"
 *J. Guidance and Control* **1**(3), 223–224 (1978).
- **Shoemake, K.** "Animating Rotation with Quaternion Curves,"
 *SIGGRAPH* 1985.
- **Shoemake, K.** "Uniform Random Rotations," in *Graphics Gems III*
 (1992), pp. 124–132.
- **Karcher, H.** "Riemannian Center of Mass and Mollifier Smoothing,"
 *Comm. Pure Appl. Math.* **30**(5), 509–541 (1977).
- **Buss, S. R. & Fillmore, J. P.** "Spherical Averages and Applications
 to Spherical Splines and Interpolation," *ACM Trans. Graph.* **20**(2),
 95–126 (2001).

Full bibliography in [`REFERENCES.md`](../REFERENCES.md).
