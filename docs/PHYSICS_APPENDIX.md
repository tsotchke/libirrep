# libirrep — Physics Conventions Appendix

The catalogue of conventions that the library adopts. Any deviation in a
user's own code will produce subtly wrong results; match exactly, or convert
explicitly at the boundary.

## Euler angles

- Ordering: **ZYZ** (physics / Sakurai).
- Ranges: α ∈ [0, 2π), β ∈ [0, π], γ ∈ [0, 2π).
- Gimbal-lock threshold: `|sin β| < 1e-12`. In that region the algorithm
  collapses α + γ into α and sets γ = 0.
- A rotation R acting actively on a vector v is, in ZYZ,
  `R = Rz(α) · Ry(β) · Rz(γ)` — read right to left.

## Phase

- **Condon-Shortley** factor `(-1)^m` is applied once, in the associated
  Legendre polynomial inside `irrep_sph_harm`. It is NOT applied a second time
  in the Wigner-D factor.
- Real spherical harmonics follow the e3nn sign convention:
  `Y^real_{l,0}   = Y^complex_{l,0}`,
  `Y^real_{l,+m}  = √2 · Re(Y^complex_{l,m})`,
  `Y^real_{l,-m}  = √2 · Im(Y^complex_{l,m})`.

## Rotations

- **Active** — a rotation sends vectors, not frames. `R v` rotates the
  vector; it does not re-coordinatize a fixed vector into a rotated frame.
- **Right-handed** — positive rotation about ẑ sends x̂ → ŷ.

## Quaternions

- Layout: `{x, y, z, w}` with the scalar `w` last.
- A unit quaternion `q = (x, y, z, w) = sin(θ/2) n̂ + cos(θ/2) ê` represents a
  rotation by angle θ about unit axis n̂.
- `q` and `-q` represent the same SO(3) rotation (double cover). Shoemake's
  sampler returns quaternions with `w ≥ 0` by convention.

## Half-integer spin

- The `_2j` suffix denotes doubled-integer arguments: all of `two_j`, `two_m`,
  `two_J`, `two_M` are twice the physical value, so spin-½ is `two_j = 1`.
- The API supports j up to `IRREP_TWO_J_MAX / 2 = 16`.

## Time reversal

- T is antiunitary: T (α |ψ⟩) = α* T |ψ⟩.
- For integer l: T acts as complex conjugation in the m-basis → identity
  on the real spherical-harmonic basis.
- For half-integer j: T = i σ_y K where K is complex conjugation. On a
  spin-½ state, T² = −1 (Kramers degeneracy).
- `irrep_time_reversal_square_sign()` returns +1 if every irrep in the
  multiset is integer-l, −1 if any is half-integer.

## Parity

- O(3) irreps are labelled `(l, p)` with `p ∈ {+1, −1}`.
- Under a tensor product, parity multiplies: `p_c = p_a · p_b`.
- `irrep_parity_filter_paths` prunes selection-rule-violating paths.

## Indexing

- `m` runs from `-l` to `+l` inclusive.
- Internal storage is `[m + l]`-indexed (so `m = -l` is at offset 0).
- Cartesian SH writes `(2l+1)` values per l; batched form writes
  `(l_max + 1)²` values total, laid out as `[l=0; l=1 (m=-1, 0, +1); l=2; …]`.
