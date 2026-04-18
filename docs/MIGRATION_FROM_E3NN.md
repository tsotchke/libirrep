# Migrating from Python e3nn to libirrep

A side-by-side cheat sheet for users fluent in Python `e3nn`. Conventions
match where possible; explicit differences are flagged.

## Irreps / multisets

| e3nn (Python) | libirrep (C) |
|---|---|
| `o3.Irreps("1x0e + 2x1o + 1x2e")` | `irrep_multiset_parse("1x0e + 2x1o + 1x2e")` |
| `irreps.dim` | `irrep_multiset_dim(m)` |
| `irreps.simplify()` | `irrep_multiset_simplify(m)` |
| `irreps + irreps2` (direct sum) | `irrep_multiset_direct_sum(m1, m2)` |

## Rotations

| e3nn | libirrep |
|---|---|
| `o3.rand_matrix()` | `irrep_quat_random` → `irrep_rot_from_quat` |
| `o3.angles_to_matrix(α, β, γ)` | `irrep_rot_from_euler_zyz({α, β, γ})` |
| `o3.matrix_to_angles(R)` | `irrep_euler_zyz_from_rot(R)` |

## Spherical harmonics

| e3nn | libirrep |
|---|---|
| `o3.spherical_harmonics(l, x, normalize=True)` | `irrep_sph_harm_cart_all(l_max, out, r_hat)` (input must be unit) |

## Wigner-D

| e3nn | libirrep |
|---|---|
| `o3.wigner_D(l, α, β, γ)` | `irrep_wigner_D_matrix(l, out, α, β, γ)` |

## Tensor products

| e3nn | libirrep |
|---|---|
| `o3.TensorProduct(irreps_in1, irreps_in2, irreps_out, instructions)` | `irrep_tp_build(a, b, c, selected_paths, num_selected_paths)` |
| `tp(a, b)` | `irrep_tp_apply(desc, a, b, c)` |
| `tp(a, b, weight)` | `irrep_tp_apply_weighted(desc, w, a, b, c)` |

## Known differences

- **Precision**: e3nn defaults to float32; libirrep defaults to float64
  (with `_f32` variants on SIMD hot paths).
- **Storage layout**: e3nn packs channels in the trailing dimension; libirrep
  uses channels-then-irrep-blocks by default (documented per function).
- **Normalisation of spherical harmonics**: libirrep uses physicist normalisation
  (orthonormal on S² under `∫ dΩ`). e3nn's `normalize=True` flag maps to this;
  `normalize=False` (which divides by √(4π) factors) is NOT directly supported
  — users should scale inputs.
- **Phase**: both use Condon-Shortley; this matches.
