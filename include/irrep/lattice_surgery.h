/* SPDX-License-Identifier: MIT */
/** @file lattice_surgery.h
 *  @brief Surface-code lattice-surgery primitives.
 *
 *  Lattice surgery is the standard procedure for performing fault-tolerant
 *  multi-qubit logical operations between adjacent surface-code patches.
 *  It works by merging two patches into one (turning on stabilizer
 *  generators along the shared boundary) and later splitting them back
 *  (turning those generators off). The merge implements a joint Pauli
 *  measurement: the type of parity measured is dictated by the type of
 *  boundary along which the two patches meet — *rough* (X-type) boundaries
 *  merge to measure `Z̄_A ⊗ Z̄_B`; *smooth* (Z-type) boundaries merge to
 *  measure `X̄_A ⊗ X̄_B`. Combined with single-patch operations these
 *  joint parities suffice to enact any Clifford circuit on the encoded
 *  qubits.
 *
 *  This module exposes the **smooth merge** of two `d × d` rotated
 *  surface-code patches A and B, glued along A's right smooth (Z-type)
 *  boundary and B's left smooth boundary. The merged geometry is a
 *  single `d × 2d` rectangular rotated surface code — a one-logical-
 *  qubit code on `2d²` data qubits — with stabilizer ordering that
 *  matches `<irrep/surface_code.h>` (bulk X then top, bottom X-boundary;
 *  bulk Z then left, right Z-boundary).
 *
 *  ## Qubit layout
 *
 *  Data qubits live at integer coordinates `(r, c)` for
 *  `r ∈ [0, d), c ∈ [0, 2d)`, linearised as `q(r, c) = r·(2d) + c`. The
 *  left half `c ∈ [0, d)` is patch A; the right half `c ∈ [d, 2d)` is
 *  patch B (B's local column `c' = c - d`).
 *
 *  ## Merged logical operators
 *
 *  The merged code's logicals follow `<irrep/surface_code.h>`:
 *
 *    `L̄_X(M) = ∏_{r=0..d-1} X_{q(r, 0)}`     (X-string on column 0;
 *                                              weight d; connects the
 *                                              top and bottom rough
 *                                              X-boundaries)
 *
 *    `L̄_Z(M) = ∏_{c=0..2d-1} Z_{q(0, c)}`    (Z-string on row 0;
 *                                              weight 2d; connects the
 *                                              left and right smooth
 *                                              Z-boundaries — now `2d`
 *                                              apart because of the merge)
 *
 *  Note: in the unmerged patch A, `L̄_X(A) = X` on column 0 has weight d;
 *  after merge this is `L̄_X(M)`. Similarly `L̄_X(B) = X` on column d of
 *  the merged frame is *also* a logical-X representative — but the two
 *  are stabilizer-equivalent (their product is a stabilizer of M).
 *
 *  ## Lattice-surgery semantics
 *
 *  Smooth merge is the joint-X parity measurement:
 *
 *    `P_joint = X̄_A · X̄_B
 *            = ∏_{r=0..d-1} X_{q(r, 0)}  ·  ∏_{r=0..d-1} X_{q(r, d)}`
 *
 *  Weight `2d`. This operator commutes with every stabilizer of the
 *  merged code AND with both merged logicals (`L̄_X(M)`, `L̄_Z(M)`),
 *  which forces `P_joint ∈ S_M` — i.e. the joint-X parity is a
 *  *stabilizer* of the merged code. Reading `P_joint`'s eigenvalue
 *  during the merge round is exactly the measurement of `X̄_A ⊗ X̄_B`
 *  on the pre-merge state.
 *
 *  (The joint-Z parity `Z̄_A · Z̄_B = ∏_{c=0..2d-1} Z_{q(0, c)} = L̄_Z(M)`
 *  is *not* a stabilizer — it equals the merged code's single logical
 *  Z. Smooth merge measures X⊗X; rough merge along a rough/rough
 *  interface would measure Z⊗Z and is the dual construction.)
 *
 *  ## Primary references
 *
 *  - Horsman-Fowler-Devitt-Van Meter, *Surface code quantum computing
 *    by lattice surgery*, New J. Phys. 14 (2012) 123011
 *    [arXiv:1111.4022].
 *  - Litinski, *A game of surface codes*, Quantum 3 (2019) 128
 *    [arXiv:1808.02892].
 *  - Fowler-Mariantoni-Martinis-Cleland, *Surface codes: Towards
 *    practical large-scale quantum computation*, Phys. Rev. A 86
 *    (2012) 032324 — operational picture of merge / split rounds.
 */
#ifndef IRREP_LATTICE_SURGERY_H
#define IRREP_LATTICE_SURGERY_H

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

/** @brief Geometry of a smooth-merge configuration: two `d × d` rotated
 *  surface-code patches glued along their smooth boundary to form a
 *  single `d × 2d` rectangular surface code. */
typedef struct {
    int d;            /**< Distance of each original patch (≥ 2). */
    int rows;         /**< Number of rows in the merged geometry = d. */
    int cols;         /**< Number of columns in the merged geometry = 2d. */
    int n_qubits;     /**< Total data qubits in the merged code = 2·d². */
    int n_X_stabs;    /**< Number of X-stabilizer generators in merged code. */
    int n_Z_stabs;    /**< Number of Z-stabilizer generators in merged code. */
} irrep_lattice_surgery_smooth_t;

/** @brief Initialise smooth-merge geometry for patches of distance d ≥ 2. */
IRREP_API irrep_status_t
irrep_lattice_surgery_smooth_init(irrep_lattice_surgery_smooth_t *out, int d);

/** @brief Build the merged CSS code = `d × 2d` rectangular surface code.
 *
 *  Stabilizer ordering matches `<irrep/surface_code.h>`: bulk X-plaquettes
 *  (row-major over `(a, b) ∈ [0, d-1) × [0, 2d-1)` with `a+b` even), then
 *  top X-boundary (r=0, b odd), then bottom X-boundary (r=d-1,
 *  (d-1+b) even); then bulk Z-plaquettes (`a+b` odd), left Z-boundary
 *  (c=0, a even), right Z-boundary (c=2d-1, (a+2d-1) odd). */
IRREP_API irrep_status_t
irrep_lattice_surgery_smooth_build(const irrep_lattice_surgery_smooth_t *p,
                                   irrep_css_code_t *out);

/** @brief Merged code's single logical X̄: X-string on column 0 of the
 *  merged frame, `L̄_X(M) = ∏_{r=0..d-1} X_{q(r, 0)}`. Weight d. */
IRREP_API irrep_status_t
irrep_lattice_surgery_smooth_logical_X(const irrep_lattice_surgery_smooth_t *p,
                                       irrep_pauli_t *out);

/** @brief Merged code's single logical Z̄: Z-string on row 0 spanning
 *  all `2d` columns, `L̄_Z(M) = ∏_{c=0..2d-1} Z_{q(0, c)}`. Weight 2d.
 *
 *  This equals `Z̄_A · Z̄_B` (the joint-Z parity of the pre-merge state),
 *  i.e. smooth merge collapses the two logical-Z operators into the
 *  single merged L̄_Z. */
IRREP_API irrep_status_t
irrep_lattice_surgery_smooth_logical_Z(const irrep_lattice_surgery_smooth_t *p,
                                       irrep_pauli_t *out);

/** @brief The joint-X parity operator `P_joint = X̄_A · X̄_B`.
 *
 *  Constructed as `X` on column 0 (qubits `q(r, 0)` for `r ∈ [0, d)`)
 *  AND column `d` (qubits `q(r, d)` for `r ∈ [0, d)`). Weight `2d`.
 *
 *  In the merged code this operator commutes with every stabilizer AND
 *  with both merged logicals; in a `[[2d², 1, d]]` code that forces
 *  `P_joint ∈ S_M`. Smooth merge therefore measures the joint-X parity
 *  of the two pre-merge logical qubits. */
IRREP_API irrep_status_t
irrep_lattice_surgery_smooth_joint_X_parity(
    const irrep_lattice_surgery_smooth_t *p, irrep_pauli_t *out);

/* ====================================================================
 * Rough merge — dual of smooth merge.
 *
 * Two `d × d` patches stacked vertically: A on top (rows `[0, d)`),
 * B on bottom (rows `[d, 2d)`). A's bottom rough (X-type) boundary
 * meets B's top rough boundary. Merged geometry: `2d × d` rectangular
 * surface code.
 *
 * Merged logicals:
 *   `L̄_X(M) = ∏_{r=0..2d-1} X_{q(r, 0)}`   (X-string on column 0
 *                                            spanning all 2d rows;
 *                                            weight 2d; equals
 *                                            `X̄_A · X̄_B`, which means
 *                                            the rough merge collapses
 *                                            both logical-X operators
 *                                            into the merged L̄_X).
 *   `L̄_Z(M) = ∏_{c=0..d-1} Z_{q(0, c)}`    (Z-string on row 0;
 *                                            weight d).
 *
 * Rough merge is the joint-Z parity measurement:
 *
 *   `P_joint = Z̄_A · Z̄_B
 *           = ∏_{c=0..d-1} Z_{q(0, c)} · ∏_{c=0..d-1} Z_{q(d, c)}`
 *
 * Weight `2d`. Commutes with every merged stabilizer AND with both
 * merged logicals — forcing `P_joint ∈ S_M`. Reading its eigenvalue
 * is the measurement of `Z̄_A ⊗ Z̄_B` on the pre-merge state.
 * ==================================================================== */

/** @brief Geometry of a rough-merge configuration. */
typedef struct {
    int d;            /**< Distance of each original patch (≥ 2). */
    int rows;         /**< Number of rows in the merged geometry = 2d. */
    int cols;         /**< Number of columns in the merged geometry = d. */
    int n_qubits;     /**< Total data qubits in the merged code = 2·d². */
    int n_X_stabs;    /**< Number of X-stabilizer generators in merged code. */
    int n_Z_stabs;    /**< Number of Z-stabilizer generators in merged code. */
} irrep_lattice_surgery_rough_t;

/** @brief Initialise rough-merge geometry for patches of distance d ≥ 2. */
IRREP_API irrep_status_t
irrep_lattice_surgery_rough_init(irrep_lattice_surgery_rough_t *out, int d);

/** @brief Build the merged CSS code = `2d × d` rectangular surface code. */
IRREP_API irrep_status_t
irrep_lattice_surgery_rough_build(const irrep_lattice_surgery_rough_t *p,
                                  irrep_css_code_t *out);

/** @brief Merged code's single logical X̄: X-string on column 0 spanning
 *  all `2d` rows. Weight `2d`. Equals `X̄_A · X̄_B`. */
IRREP_API irrep_status_t
irrep_lattice_surgery_rough_logical_X(const irrep_lattice_surgery_rough_t *p,
                                      irrep_pauli_t *out);

/** @brief Merged code's single logical Z̄: Z-string on row 0 spanning
 *  all `d` columns. Weight `d`. */
IRREP_API irrep_status_t
irrep_lattice_surgery_rough_logical_Z(const irrep_lattice_surgery_rough_t *p,
                                      irrep_pauli_t *out);

/** @brief Joint-Z parity `Z̄_A · Z̄_B = Z` on row 0 AND row `d` (each
 *  spanning all `d` columns). Weight `2d`. In the merged code this
 *  operator commutes with every stabilizer AND with both merged
 *  logicals, so it lies in `S_M`. Rough merge measures `Z̄_A ⊗ Z̄_B`. */
IRREP_API irrep_status_t
irrep_lattice_surgery_rough_joint_Z_parity(
    const irrep_lattice_surgery_rough_t *p, irrep_pauli_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_LATTICE_SURGERY_H */
