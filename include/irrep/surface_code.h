/* SPDX-License-Identifier: MIT */
/** @file surface_code.h
 *  @brief Rotated planar surface code [[d², 1, d]].
 *
 *  The rotated surface code is the workhorse of contemporary
 *  superconducting fault-tolerant experiments (Google distance-7
 *  logical qubit, IBM, Quantinuum). For distance `d ≥ 2` it encodes
 *  exactly **one** logical qubit on a `d × d` square lattice of data
 *  qubits with two pairs of opposing boundaries: a *rough* (X-type)
 *  pair and a *smooth* (Z-type) pair.
 *
 *  ## Lattice convention
 *
 *  Data qubits live at integer coordinates `(r, c)` for `r, c ∈ [0, d)`,
 *  linearised as `q(r, c) = r·d + c`. Stabilizers fall into four
 *  geometric classes:
 *
 *  - **Bulk plaquettes** at `(a, b) ∈ [0, d-1)²`: weight-4 operator on
 *    `{q(a,b), q(a+1,b), q(a,b+1), q(a+1,b+1)}`. X-type if `a+b` is
 *    even, Z-type if odd. Total `(d-1)²` of them; X- and Z-counts
 *    differ by at most 1.
 *
 *  - **Top boundary X-stabs** (rough, r=0): for each `b ∈ [0, d-1)`
 *    with `b` odd, weight-2 X on `{q(0,b), q(0,b+1)}`.
 *
 *  - **Bottom boundary X-stabs** (rough, r=d-1): for each `b ∈ [0, d-1)`
 *    with `(d-1+b)` even, weight-2 X on `{q(d-1,b), q(d-1,b+1)}`.
 *
 *  - **Left boundary Z-stabs** (smooth, c=0): for each `a ∈ [0, d-1)`
 *    with `a` even, weight-2 Z on `{q(a,0), q(a+1,0)}`.
 *
 *  - **Right boundary Z-stabs** (smooth, c=d-1): for each `a ∈ [0, d-1)`
 *    with `(a+d-1)` odd, weight-2 Z on `{q(a,d-1), q(a+1,d-1)}`.
 *
 *  This gives `[[d², 1, d]]` for any `d ≥ 2`. The CSS orthogonality
 *  condition `H_X · H_Z^T = 0 (mod 2)` is verifiable by hand: every
 *  pair of stabilizers shares 0 or 2 data qubits.
 *
 *  ## Counts
 *
 *  | d  | n_qubits | n_X_stabs | n_Z_stabs | total stabs | k |
 *  |----|----------|-----------|-----------|-------------|---|
 *  | 2  | 4        | 1         | 2         | 3           | 1 |
 *  | 3  | 9        | 4         | 4         | 8           | 1 |
 *  | 5  | 25       | 12        | 12        | 24          | 1 |
 *  | 7  | 49       | 24        | 24        | 48          | 1 |
 *  | 9  | 81       | 40        | 40        | 80          | 1 |
 *
 *  ## Primary references
 *
 *  - Bravyi-Kitaev, *Quantum codes on a lattice with boundary*,
 *    arXiv:quant-ph/9811052 (1998).
 *  - Fowler-Mariantoni-Martinis-Cleland, *Surface codes*, Phys. Rev.
 *    A 86 (2012) 032324.
 *  - Litinski, *A game of surface codes*, Quantum 3 (2019) 128.
 *  - Acharya et al. (Google), *Suppressing quantum errors by scaling
 *    a surface code logical qubit*, Nature 614 (2023) 676 — distance-7.
 */
#ifndef IRREP_SURFACE_CODE_H
#define IRREP_SURFACE_CODE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/stabilizer_group.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Surface-code geometry parameters. */
typedef struct {
    int distance;     /**< Code distance d (≥ 2). */
    int n_qubits;     /**< Number of data qubits = d². */
    int n_X_stabs;    /**< Number of X-stabilizer generators. */
    int n_Z_stabs;    /**< Number of Z-stabilizer generators. */
} irrep_surface_params_t;

/** @brief Initialise surface-code parameters for distance d ≥ 2. */
IRREP_API irrep_status_t
irrep_surface_init(irrep_surface_params_t *out, int distance);

/** @brief Build the surface code as a CSS code.
 *
 *  Stabilizer ordering: bulk X (row-major over (a,b)), then top X,
 *  bottom X (so all of `H_X`); bulk Z, then left Z, right Z (all of
 *  `H_Z`). Caller must `irrep_css_code_free` when done. */
IRREP_API irrep_status_t
irrep_surface_build(const irrep_surface_params_t *p, irrep_css_code_t *out);

/** @brief Build the canonical logical X̄ operator for the rotated
 *  surface code: an X-string of length d running along column 0 from
 *  the top rough boundary to the bottom rough boundary.
 *
 *  L_X = ∏_{r=0..d-1} X on qubit q(r, 0).
 *
 *  Properties:
 *    - Weight d (matches the code distance).
 *    - Commutes with every stabilizer (each Z-stabilizer either has
 *      no support on column 0 or has exactly two qubits there).
 *    - Anti-commutes with `irrep_surface_logical_Z` on the shared
 *      corner qubit q(0, 0).
 *
 *  The Pauli `out` is allocated by this function; caller must
 *  `irrep_pauli_free` when done. */
IRREP_API irrep_status_t
irrep_surface_logical_X(const irrep_surface_params_t *p, irrep_pauli_t *out);

/** @brief Build the canonical logical Z̄ operator: a Z-string of
 *  length d running along row 0 from the left smooth boundary to the
 *  right smooth boundary.
 *
 *  L_Z = ∏_{c=0..d-1} Z on qubit q(0, c).
 *
 *  Properties: weight d, commutes with every stabilizer, anti-commutes
 *  with `irrep_surface_logical_X` on q(0, 0). */
IRREP_API irrep_status_t
irrep_surface_logical_Z(const irrep_surface_params_t *p, irrep_pauli_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SURFACE_CODE_H */
