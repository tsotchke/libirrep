# libirrep

> SO(3) / SU(2) / O(3) irreducible-representation machinery in portable C11.

`libirrep` provides Clebsch-Gordan coefficients, real and complex spherical harmonics,
Wigner-D matrices, e3nn-style tensor products, and disciplined rotation math
(quaternions, Euler-ZYZ, axis-angle, Shoemake-uniform sampling) as a single
library with a stable ABI and per-triple release artifacts.

It is the mathematical backbone for E(3)-equivariant neural networks (NequIP,
MACE, Allegro, TENN, anything descended from `e3nn`), spin systems, and any C
code that needs a careful rotation library.

## Status

v1.0 in development. See `docs/DESIGN.md` for the spec, `CHANGELOG.md` for
release history, and `TODO.md` for remaining work.

## Scope

Provides, in pure C11:

- Real, complex, and cartesian spherical harmonics up to `l = 8` (aspirationally `l = 16`)
- Clebsch-Gordan coefficients for integer and half-integer spin, via the Racah formula
- Wigner-D matrices and small-d matrices, numerically stable at large `j` via Jacobi polynomials
- Wigner 3j / 6j / 9j symbols
- O(3) tensor products with `e3nn`-style path-indexed sparse decomposition, forward and backward
- SO(3) / SU(2) conversions, composition, exp/log (Rodrigues), SLERP, Fréchet mean, Shoemake uniform sampling
- Irrep-multiset algebra with `"1x0e + 2x1o + 1x2e"` string parsing
- Radial basis functions (Bessel, Gaussian) and smooth cutoffs, for NequIP/MACE-style networks
- Time-reversal and parity operators on irreps
- Minimal equivariant-layer primitives (linear-on-irreps, RMS norm, gate activation)
- Lebedev and Gauss-Legendre quadrature

Does not provide, in v1.0:

- GPU kernels (planned for v1.1; internal API is GPU-ready)
- SE(3), SE(2), SU(3), or other Lie groups
- Python bindings (by design; consumers wrap via the stable C ABI)

## Build

```
make            # build the library and run tests
make test       # run tests only
make bench      # run benchmarks, write JSON to benchmarks/results/
make examples   # build the three example programs
make release    # produce release/<VERSION>/ artifacts
```

Also supports CMake:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

## Quickstart

```c
#include <stdio.h>
#include <irrep/wigner_d.h>

int main(void) {
    /* Wigner-D matrix for j = 1 at (alpha, beta, gamma) = (0.1, 0.2, 0.3) */
    double _Complex D[9];
    irrep_wigner_D_matrix(1, D, 0.1, 0.2, 0.3);
    for (int mp = 0; mp < 3; ++mp) {
        for (int m = 0; m < 3; ++m) {
            double _Complex z = D[mp * 3 + m];
            printf("D[%+d][%+d] = % .6f %+ .6fi\n", mp - 1, m - 1,
                   (double)__real__ z, (double)__imag__ z);
        }
    }
    return 0;
}
```

Compile: `cc -o demo demo.c -I/path/to/libirrep/include -L/path/to/libirrep/lib -llibirrep -lm`

## Conventions

- **Euler angles**: ZYZ physics convention. `alpha ∈ [0, 2π)`, `beta ∈ [0, π]`, `gamma ∈ [0, 2π)`.
- **Phase**: Condon-Shortley `(-1)^m` throughout.
- **Rotations**: active, right-handed.
- **Quaternions**: `{x, y, z, w}` layout.
- **Precision**: `double` by default; `_f32` variants on SIMD hot paths.

See `docs/PHYSICS_APPENDIX.md` for the full conventions catalogue, and
`docs/MIGRATION_FROM_E3NN.md` for side-by-side mapping from Python `e3nn`.

## License

MIT. See `LICENSE`.
