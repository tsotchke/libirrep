# Tutorials

Eight numbered walkthroughs of libirrep's surface, ordered for a reader
who has already cloned, built, and passed `make test`:

| # | Tutorial                              | Core primitive                               |
|--:|---------------------------------------|----------------------------------------------|
| 1 | [Rotations](01_rotations.md)          | `irrep_quaternion_t`, SLERP, Fréchet mean   |
| 2 | [Spherical harmonics](02_spherical_harmonics.md) | `irrep_sph_harm_cart_all`         |
| 3 | [Clebsch-Gordan](03_clebsch_gordan.md) | `irrep_cg`, `irrep_wigner_3j`              |
| 4 | [Wigner-D](04_wigner_d.md)            | `irrep_wigner_d_matrix`, `irrep_wigner_D`    |
| 5 | [Tensor products](05_tensor_products.md) | `irrep_tp_build`, `irrep_tp_apply`        |
| 6 | [Equivariant NNs](06_equivariant_nn.md) | `irrep_linear_t`, `irrep_nequip_layer_t`   |
| 7 | [Kagome NQS substrate](07_kagome_nqs_substrate.md) | lattice + space-group + RDM pipeline |
| 8 | [Materials-search pipeline](08_materials_search.md) | `lattice3d.h` + cubic point groups + `dmi.h` (DMI / J^s symmetry analyzer) |

## 1.4 — quantum-error-correction substrate

The 1.4.0 cycle adds an 18-module quantum-error-correction layer (toric,
surface, color, BB qLDPC, hypergraph and lifted product, honeycomb and
CSS Floquet, 3D toric, X-cube fracton, HaPPY perfect tensor, single-shot
meta-checks, Bacon-Shor subsystem, Steane⊗Steane concatenation,
generic-color-code framework). These are not yet in the numbered
tutorial sequence; the entry points are:

| Topic | Header | Test (read as worked example) |
|---|---|---|
| Symplectic Pauli + commutativity | `<irrep/stabilizer_group.h>` | `tests/test_stabilizer_group.c` |
| CSS code from `(H_X, H_Z)` | `<irrep/css_code.h>` | `tests/test_css_code.c` |
| Pauli arithmetic + brute-force distance | `<irrep/qec_distance.h>` | `tests/test_qec_distance.c` |
| Surface code `[[d², 1, d]]` | `<irrep/surface_code.h>` | `tests/test_surface_code.c` |
| BB qLDPC `[[72, 12, 6]]` | `<irrep/bivariate_bicycle.h>` | `tests/test_bivariate_bicycle.c` |
| Honeycomb Floquet | `<irrep/honeycomb_floquet.h>` | `tests/test_honeycomb_floquet.c` |
| 3D toric | `<irrep/toric3d.h>` | `tests/test_toric3d.c` |
| X-cube fracton | `<irrep/xcube_code.h>` | `tests/test_xcube_code.c` |

See also `docs/qec_scoping.md` (5-tier code-family survey) and
`docs/qec_research_roadmap.md` (deferred research-track items).

## 1.3 additions — also covered by ready-to-run examples

The 1.3.0-alpha cycle landed three new modules (3D Bravais lattices,
cubic point groups, bond-exchange-tensor symmetry analyzer). These are
not yet in the numbered tutorial sequence but each ships with a
ready-to-run example that doubles as a tutorial:

| Module | Example | What it shows |
|---|---|---|
| `lattice3d.h` (5 Bravais families) | `examples/lattice3d_demo.c` | Geometry walkthrough — SC/BCC/FCC/Diamond/Pyrochlore |
| `lattice3d.h` + `hamiltonian.h` | `examples/lattice3d_heisenberg_ed.c` | 3D Heisenberg ED with K₈,₈ closed-form cross-check |
| `lattice3d.h` (Γ-momentum sector) | `examples/lattice3d_sector_ed.c` | Translation orbit canonicalisation |
| `lattice3d.h` (full BZ) | `examples/lattice3d_kspace_ed.c` | k-resolved tower-of-states diagnostic |
| Pyrochlore physics | `examples/pyrochlore16_heisenberg.c`, `_j1j2.c`, `_correlations.c` | 16-site frustrated 3D Heisenberg + sweep + correlations |
| Cubic point groups (T_d / O_h / O) | `examples/cubic_crystal_field.c` | Crystal-field d-orbital splitting eg + t_2g |
| `dmi.h` — DMI + symmetric exchange | `examples/dmi_pyrochlore_pattern.c`, `dmi_kagome_pattern.c` | Bond-exchange-tensor decomposition; B20 / kagome RT-magnet candidates |

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
- **Phase:** Condon–Shortley
- **Quaternions:** `{x, y, z, w}`, unit-norm, Shoemake-sampled with `w ≥ 0`
- **Complex:** `double _Complex` (`<complex.h>`)
- **Half-integer spin:** doubled-integer API, `_2j` suffix

Derivations in [`../PHYSICS_APPENDIX.md`](../PHYSICS_APPENDIX.md).
