# Tutorials

Seven walkthroughs of libirrep's surface, ordered for a reader who has
already cloned, built, and passed `make test`:

| # | Tutorial                              | Core primitive                               |
|--:|---------------------------------------|----------------------------------------------|
| 1 | [Rotations](01_rotations.md)          | `irrep_quaternion_t`, SLERP, Fréchet mean   |
| 2 | [Spherical harmonics](02_spherical_harmonics.md) | `irrep_sph_harm_cart_all`         |
| 3 | [Clebsch-Gordan](03_clebsch_gordan.md) | `irrep_cg`, `irrep_wigner_3j`              |
| 4 | [Wigner-D](04_wigner_d.md)            | `irrep_wigner_d_matrix`, `irrep_wigner_D`    |
| 5 | [Tensor products](05_tensor_products.md) | `irrep_tp_build`, `irrep_tp_apply`        |
| 6 | [Equivariant NNs](06_equivariant_nn.md) | `irrep_linear_t`, `irrep_nequip_layer_t`   |
| 7 | [Kagome NQS substrate](07_kagome_nqs_substrate.md) | lattice + space-group + RDM pipeline |

## Running the code blocks

Every C code block in these tutorials has been compiled against
`build/lib/liblibirrep.a` and verified to produce the stated output.
Cold-read repro:

```
git clone https://github.com/tsotchke/libirrep.git && cd libirrep
make                # lib + tests + examples
cd docs/tutorials   # follow the numbered .md files in order
```

To verify a tutorial's snippet, copy the code block into a `.c` file,
compile against the built library, and compare output:

```
cc -std=c11 -I/path/to/libirrep/include snippet.c \
   /path/to/libirrep/build/lib/liblibirrep.a -lm -o /tmp/snip && /tmp/snip
```

## Numerical references

Exact numerical output for every shipped example is catalogued in
the top-level `examples/EXPECTED_OUTPUT.md`,
with the RNG seed, tolerance, and cross-reference citation per
example. If a build produces different numbers, that file is the
first place to look for what to expect.

## Conventions cheat-sheet

All tutorials assume the single library-wide convention:

- **Angles:** radians
- **Rotations:** active, right-handed
- **Euler:** ZYZ (Sakurai / physics convention)
- **Phase:** Condon-Shortley
- **Quaternions:** `{x, y, z, w}`, unit-norm, Shoemake-sampled with `w ≥ 0`
- **Complex:** `double _Complex` (`<complex.h>`)
- **Half-integer spin:** doubled-integer API, `_2j` suffix

Derivations in [`../PHYSICS_APPENDIX.md`](../PHYSICS_APPENDIX.md).
