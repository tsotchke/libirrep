# libirrep — Design

See also:

- `docs/PHYSICS_APPENDIX.md` — convention catalogue.
- `docs/REFERENCES.md` — full bibliography.

## Module graph

```
                     ┌────────────────┐
                     │    types.h     │  labels, multisets, rotation types, status
                     └──────┬─────────┘
            ┌───────────────┼────────────────────────────────┐
            ▼               ▼                                ▼
       ┌─────────┐   ┌─────────────┐                 ┌───────────────┐
       │  so3.h  │   │   simd.h    │                 │   multiset.h  │
       │  su2.h  │   │   version.h │                 │   parity.h    │
       └────┬────┘   └─────────────┘                 │ time_reversal │
            │                                        └───────┬───────┘
            ▼                                                ▼
   ┌──────────────────┐   ┌────────────────────┐   ┌──────────────────┐
   │  spherical_      │   │  clebsch_gordan.h  │   │  tensor_product  │
   │   harmonics.h    │   │  recoupling.h      │   │   .h             │
   └───────┬──────────┘   └─────────┬──────────┘   └─────────┬────────┘
           └──────────────┬─────────┘                        │
                          ▼                                  ▼
                   ┌──────────────┐                ┌────────────────────┐
                   │  wigner_d.h  │                │ equivariant_layers │
                   │  quadrature  │                │       .h           │
                   │  radial.h    │                └────────────────────┘
                   └──────────────┘
```

## Numerical conventions (authoritative)

- **Precision**: `double` default; `float` for `_f32` SIMD variants.
- **Euler**: ZYZ physics convention. `α ∈ [0, 2π)`, `β ∈ [0, π]`, `γ ∈ [0, 2π)`.
- **Phase**: Condon-Shortley `(-1)^m`, applied once (in associated Legendre).
- **Rotations**: active, right-handed.
- **Quaternion layout**: `{x, y, z, w}`, with `w ≥ 0` canonical sign.
- **Tests**: `1e-10` relative for double, `1e-5` for single.

## Threading

- All pure math functions are reentrant and thread-safe.
- Built tables (`cg_table_t`, `tp_descriptor_t`, `irrep_multiset_t`,
  `irrep_linear_t`) are safe to read concurrently once built; not safe to
  build/free concurrently with reads.
- `irrep_last_error()` is thread-local.
- SIMD feature table is populated once; readers see a fully-initialised value.

## ABI policy

- Public symbols have default visibility; everything else is hidden.
- Any change to exported symbol signatures or public struct layouts requires a
  MAJOR version bump.
- `irrep_abi_hash()` returns the SHA-256 of the canonicalised export surface,
  baked in at release time.

## GPU extensibility (v1.1 preview)

Hot paths are tagged `// HOT: GPU kernel target in v1.1`. The internal API
uses function-pointer dispatch so a Metal / CUDA / HIP / Vulkan backend can be
slotted in without ABI churn.

## Open items

Tracked in `TODO.md` and in the authoritative plan file.
