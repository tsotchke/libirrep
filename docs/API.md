# libirrep — API Reference

This document is the human-readable entry point. Complete per-function
documentation is in the public headers under `include/irrep/`. `make docs`
generates Doxygen HTML under `build/docs/`.

## Modules

| Header | Module | Key types / functions |
|---|---|---|
| `<irrep/types.h>` | core types | `irrep_quaternion_t`, `irrep_rot_matrix_t`, `irrep_euler_zyz_t`, `irrep_axis_angle_t`, `irrep_label_t`, `irrep_multiset_t`, `irrep_status_t`, `irrep_strerror`, `irrep_last_error` |
| `<irrep/version.h>` | versioning | `IRREP_VERSION_*`, `irrep_version_*()`, `irrep_abi_hash()` |
| `<irrep/simd.h>` | CPU features | `irrep_cpu_features()`, `irrep_cpu_has_*()` |
| `<irrep/so3.h>` | SO(3) ops | rotation conversions, compose, inverse, exp, log, Shoemake random, SLERP, Fréchet mean |
| `<irrep/su2.h>` | SU(2) ops | Pauli matrices, SU(2) ↔ quaternion, SU(2) → SO(3) double cover |
| `<irrep/clebsch_gordan.h>` | CG coefficients | `irrep_cg`, `irrep_cg_2j`, `irrep_wigner_3j`, cached tables |
| `<irrep/spherical_harmonics.h>` | spherical harmonics | complex, real, cartesian, gradient, associated Legendre |
| `<irrep/wigner_d.h>` | Wigner-D | small-d, full D, matrix, block-diag on multiset |
| `<irrep/multiset.h>` | irrep-multiset algebra | parse "1x0e + 2x1o", simplify, direct sum |
| `<irrep/tensor_product.h>` | tensor products | e3nn-style path-indexed, weighted, batched, forward + backward |
| `<irrep/recoupling.h>` | 6j / 9j | `irrep_wigner_6j`, `irrep_wigner_9j`, `irrep_racah_w` |
| `<irrep/radial.h>` | radial basis | Bessel / Gaussian RBFs, cosine / polynomial cutoffs |
| `<irrep/quadrature.h>` | quadrature | Lebedev (S²), Gauss-Legendre |
| `<irrep/time_reversal.h>` | time reversal | integer / half-integer T operator, T² sign |
| `<irrep/parity.h>` | parity | parity, parity product, path filter |
| `<irrep/equivariant_layers.h>` | NN primitives | linear-on-irreps, RMS norm, gate |
| `<irrep/irrep.h>` | umbrella | pulls every module |

See `docs/tutorials/` for step-by-step walkthroughs.
