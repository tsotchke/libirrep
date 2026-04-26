# Tutorial 8 — Materials-search pipeline: from space group to allowed exchange tensor

This tutorial walks through libirrep's materials-search workflow:
given a candidate crystal, derive the symmetry-allowed bond-exchange-
tensor structure (Heisenberg + DMI + Kitaev-Γ-type anisotropy) from
group theory alone, with no input from the materials literature. This
is the irreducible step that crystallographers do by hand from the
*International Tables for Crystallography* vol. A (Hahn 2005) and that
no DFT or micromagnetic package automates — libirrep handles it as a
single library call.

By the end of this tutorial you will:

1. Build a 3D crystal lattice with libirrep's geometry layer.
2. Pick a candidate site-symmetry point group.
3. Run the bond-exchange-tensor symmetry analyzer.
4. Interpret the result against textbook materials results (Bak-Jensen
   for B20 chiral magnets, Curnoe-Ross-Kao for pyrochlore).

We assume you've already worked through tutorial 07 (kagome NQS
substrate). The new modules in this tutorial are:

- `<irrep/lattice3d.h>` — 3D Bravais lattices (5 families)
- `<irrep/point_group.h>` — cubic point groups added in 1.3 (T_d, O_h, O)
- `<irrep/dmi.h>` — DMI vector + symmetric exchange tensor analyzers

## Why the bond-exchange-tensor decomposition matters

The full bilinear two-spin coupling on a bond is a 3×3 matrix `J_ij`:

```
H_ij = S_i · J_ij · S_j     (9 components in general)
```

Decomposed:

- **Symmetric part** `J^s_ij = (J_ij + J_ji)/2` (6 components):
  - 1 isotropic Heisenberg (trace, the `J` everyone knows)
  - 5 anisotropic (the Kitaev-Γ-type couplings that drive
    quantum-spin-liquid phases in α-RuCl₃, Tb₂Ti₂O₇, etc.)
- **Antisymmetric part** `J^a_ij`, encoded as the DMI vector
  `D_ij` via `J^a_{ab} = ε_{abc} D^c` (3 components):
  - This is the spin-orbit-coupled, parity-broken term that drives
    skyrmions in MnSi, Cu₂OSeO₃, and the W/Ta/CoFeB/MgO multilayers
    used in modern spintronic prototypes.

Crystal symmetry constrains which of these 9 components are
permitted. Moriya 1960 wrote down five rules for the antisymmetric
DMI vector; the analogous derivation for the symmetric tensor lives
in `docs/PHYSICS_APPENDIX.md` §13. libirrep returns the
symmetry-allowed `(D, J^s)` basis as orthonormal vectors / matrices —
the parameter scaffolding that DFT or micromagnetic codes need as
input.

## Worked example: kagome NN bond under D_6

The kagome lattice is the magnetic sublattice of room-temperature
metallic kagome ferromagnets (Fe₃Sn₂ at T_C ≈ 660 K, Co₃Sn₂S₂ at
T_C = 177 K) and antiferromagnets (Mn₃Sn / Mn₃Ge at T_N ≈ 380–420 K).
The in-plane site symmetry is D_6h in idealised symmetric stacking;
real materials break σ_h via bilayer stacking (Fe₃Sn₂) or by the
σ_h-perpendicular ordering vector (Co₃Sn₂S₂). The libirrep cubic-2D
group D_6 handles the chiral in-plane subgroup.

```c
/* 08_materials_search.c — full snippet, drop into any source dir */
#include <irrep/dmi.h>
#include <irrep/lattice.h>
#include <irrep/point_group.h>

#include <math.h>
#include <stdio.h>

int main(void) {
    /* Build a small kagome cluster (we only need positions). */
    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
    double xy0[2], xy1[2];
    irrep_lattice_site_position(L, 0, xy0);  /* sub 0 (corner of cell 0) */
    irrep_lattice_site_position(L, 1, xy1);  /* sub 1 (NN to sub 0) */

    /* Centre the bond at origin so the point group (which acts about
     * the origin) realises the bond's site stabiliser. Pad to 3D. */
    double mid[3] = {0.5*(xy0[0]+xy1[0]), 0.5*(xy0[1]+xy1[1]), 0};
    double r_a[3] = {xy0[0]-mid[0], xy0[1]-mid[1], 0};
    double r_b[3] = {xy1[0]-mid[0], xy1[1]-mid[1], 0};

    /* Run the analyzer under the chiral hexagonal point group D_6. */
    irrep_pg_table_t *D6 = irrep_pg_table_build(IRREP_PG_D6);
    double D_basis[9] = {0};
    int n_D = irrep_dmi_allowed_basis_from_pg(r_a, r_b, D6, 1e-6, D_basis);

    double Js_basis[54] = {0};
    int n_Js = irrep_exchange_symmetric_basis_from_pg(
        r_a, r_b, D6, 1e-6, Js_basis);

    printf("Kagome NN bond under D_6:\n");
    printf("  DMI: %d-dim", n_D);
    if (n_D > 0)
        printf("  D = (%+.4f, %+.4f, %+.4f)",
               D_basis[0], D_basis[1], D_basis[2]);
    printf("\n  J^s: %d-dim\n", n_Js);

    irrep_pg_table_free(D6);
    irrep_lattice_free(L);
    return 0;
}
```

Compile and run:

```
cc -std=c11 -I/path/to/libirrep/include 08_materials_search.c \
   /path/to/libirrep/build/lib/liblibirrep.a -lm -o /tmp/m && /tmp/m
```

Expected output:

```
Kagome NN bond under D_6:
  DMI: 1-dim  D = (+1.0000, +0.0000, +0.0000)
  J^s: 3-dim
```

The DMI vector is **along the bond direction** (since the libirrep D_6
puts the kagome NN along x). This is Moriya rule (e): a C₂' axis
along the bond pins D parallel to it. For Fe₃Sn₂ (where stacking
breaks σ_h, leaving D_6 as the effective in-plane group), this
reproduces the in-plane DMI that drives the observed RT skyrmion-
bubble phase.

The 3-dim symmetric tensor decomposes (in the analyzer's choice of
orthonormal basis) as Heisenberg + two anisotropic uniaxial /
Γ-type components — the parametrisation that DFT then has to
fix the magnitude on.

## What changes when the symmetry is more / less restrictive

Try the same bond under the chiral subgroup D_3 (= D_6 minus the
edge-midpoint C₂' axes):

```c
irrep_pg_table_t *D3 = irrep_pg_table_build(IRREP_PG_D3);
n_D  = irrep_dmi_allowed_basis_from_pg(r_a, r_b, D3, 1e-6, D_basis);
n_Js = irrep_exchange_symmetric_basis_from_pg(r_a, r_b, D3, 1e-6, Js_basis);
printf("Under D_3:  DMI %d-dim,  J^s %d-dim\n", n_D, n_Js);
```

Expected:

```
Under D_3:  DMI 1-dim,  J^s 4-dim
```

Same DMI direction, but J^s gains a component (4 instead of 3) —
because D_3 has fewer constraints than D_6. The extra freedom is a
Kitaev-Γ-type yz off-diagonal coupling that emerges when the in-plane
mirrors are absent.

For materials where the actual site symmetry is D_3 instead of D_6
(e.g., kagome compounds with intrinsic chirality from the substrate),
this 4th component is what you need to model. The analyzer tells you
exactly when to include it.

## The materials-design loop

The full loop, with libirrep's lane highlighted:

| step | tool | who runs it |
|------|------|-------------|
| 1. Pick a candidate space group from a crystal database | ICSD / COD / MAGNDATA | crystallographer |
| 2. Encode the lattice + site point group | libirrep `lattice.h` / `lattice3d.h` / `point_group.h` | **libirrep** |
| 3. Iterate over symmetry-distinct bond classes | libirrep | **libirrep** |
| 4. Derive symmetry-allowed (D, J^s) per bond class | **libirrep `dmi.h`** | **libirrep** |
| 5. Get magnitude on each allowed component | DFT (VASP / QE / Wannier-projected hopping) | DFT specialist |
| 6. Simulate skyrmion / domain-wall stability vs bias | Micromagnetic (mumax / OOMMF) | spintronics group |
| 7. Synthesise + measure | Experiment | materials physicist |

Step 4 is what was previously hand-derivation from the *International
Tables*. Automating it is the prerequisite for high-throughput
candidate-space-group enumeration rather than one-material-at-a-time
analysis. Together with the rest of the libirrep stack, you get
the full from-symmetry-to-Hamiltonian step in C, with no Python or
external dependencies.

## See also

- `examples/dmi_pyrochlore_pattern.c` — pyrochlore NN bond under
  O_h vs O vs T_d, demonstrates the inversion-driven D = 0 case
- `examples/dmi_kagome_pattern.c` — fuller kagome catalog (the
  source of this tutorial's snippet, with extended interpretation)
- `examples/cubic_crystal_field.c` — companion that decomposes
  l-orbitals under O_h / T_d (the Eg + T2g splitting that defines
  transition-metal complex spectra, derived from group theory)
- `docs/PHYSICS_APPENDIX.md` §13 — full mathematical derivation
  of Moriya's rules, the symmetric tensor projector, and the
  axial-vector argument that makes O_h and O agree on J^s
- `docs/PHYSICS_RESULTS.md` §7 — the validated reproductions of
  textbook patterns (Bak-Jensen, Curnoe-Ross-Kao) using only this
  pipeline
