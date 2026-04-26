/* SPDX-License-Identifier: MIT */
/** @file dmi.h
 *  @brief Symmetry-allowed Dzyaloshinskii-Moriya vectors for a magnetic bond.
 *
 *  The Dzyaloshinskii-Moriya interaction (Dzyaloshinskii 1958, Moriya 1960)
 *  is the antisymmetric part of the spin-orbit-coupled exchange:
 *
 *      H_DMI = sum_{<i,j>}  D_ij · ( S_i × S_j )
 *
 *  with the antisymmetry constraint `D_ji = -D_ij`. DMI is responsible for
 *  spin-canting in centrosymmetric magnets, helimagnetism in chiral cubic
 *  B20 compounds (MnSi, FeGe), and skyrmion stabilisation in interfacial
 *  multilayers (Pt/Co/Ir, Ta/CoFeB/MgO). For each non-trivial bond class,
 *  the magnetic point group at the bond constrains `D_ij` to a specific
 *  sub-space — Moriya's five rules (PRB **120**, 91 (1960) §III):
 *
 *    1. Inversion centre at bond midpoint           → D = 0
 *    2. Mirror plane perpendicular to bond, through midpoint → D ⊥ bond
 *    3. Mirror plane containing the bond, through midpoint  → D ⊥ that plane
 *    4. C₂ axis through midpoint perpendicular to bond → D ⊥ that axis
 *    5. C_n axis along the bond (n ≥ 2)              → D ∥ that axis
 *
 *  These rules are the IUCr-international-tables data hand-applied by
 *  crystallographers when studying chiral magnets. This module automates
 *  them: given a bond and a list of point-group operations, return the
 *  symmetry-allowed `D`-vector subspace as an orthonormal basis (0, 1, 2,
 *  or 3 vectors).
 *
 *  Algorithm: build the 3×3 projector
 *
 *      P = (1/|S|) · sum_{g ∈ S} M_g(g)
 *
 *  with `M_g = +R_proper(g)` if `g` preserves the bond, `−R_proper(g)`
 *  if `g` reverses it. `R_proper(g)` is the orientation-preserving part
 *  of `g` (the axial-vector representation collapses out the `det(g)`
 *  factor). Diagonalise `P`; the +1 eigenspace is the allowed
 *  D-subspace.
 */
#ifndef IRREP_DMI_H
#define IRREP_DMI_H

#include <irrep/export.h>
#include <irrep/point_group.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief A single rigid symmetry operation acting on cartesian space.
 *
 *  The operation maps `r → det · R_proper · r`. `det = +1` is a proper
 *  rotation; `det = -1` is improper (reflection / S_n / inversion).
 *
 *  Set @p antiunitary to `1` for time-reversal-augmented operations
 *  (`T · g`), to `0` for plain spatial operations. For analyzers that
 *  ignore this field (the existing DMI / symmetric-exchange ones, when
 *  `antiunitary = 0`), the value is silently ignored. The
 *  symmetric-exchange tensor analyzer is invariant to the antiunitary
 *  flag — both pre- and post-`T` signs cancel for rank-2 axial-axial →
 *  polar tensors. The DMI analyzer flips the sign rule for antiunitary
 *  operations: see §14 of `docs/PHYSICS_APPENDIX.md`. */
typedef struct {
    double R_proper[9]; /**< 3×3 row-major proper-rotation matrix. */
    int    det;         /**< +1 (proper) or -1 (improper). */
    int    antiunitary; /**< 0 = unitary spatial; 1 = `T · g` (time-reversal augmented). */
} irrep_dmi_sym_op_t;

/** @brief Compute the symmetry-allowed DMI subspace for a bond.
 *
 *  Iterates over @p ops, applies each to the bond endpoints, classifies
 *  the operation as preserving (g(r_a) = r_a, g(r_b) = r_b) or
 *  reversing (g(r_a) = r_b, g(r_b) = r_a), and accumulates the projector.
 *  Operations that do not map `{r_a, r_b}` to itself (as an unordered
 *  set, within tolerance @p tol) are silently skipped — they are not in
 *  the bond's site-symmetry group.
 *
 *  At least the identity (in @p ops) should preserve any bond, so the
 *  worst-case output is `n_basis = 3` (no constraint). The strongest
 *  case is `n_basis = 0` (D ≡ 0).
 *
 *  @param r_a, r_b   Cartesian positions of the bond endpoints.
 *  @param ops        Array of @p n_ops candidate symmetry operations.
 *  @param n_ops      Length of @p ops.
 *  @param tol        Tolerance for matching positions (e.g. `1e-9`).
 *  @param basis_out  9-element caller buffer; on return, holds the
 *                    `n_basis` row vectors of the orthonormal allowed
 *                    basis (each row 3 doubles, contiguous).
 *  @return Dimension of the allowed D-subspace (0..3). */
IRREP_API int irrep_dmi_allowed_basis(const double r_a[3], const double r_b[3],
                                      const irrep_dmi_sym_op_t *ops, int n_ops, double tol,
                                      double basis_out[9]);

/** @brief Convenience wrapper: take the operations from a libirrep
 *         point-group table directly.
 *
 *  Identical to #irrep_dmi_allowed_basis except the operations are
 *  read from the table via #irrep_pg_element. The point-group origin
 *  is implicit: the elements rotate space about (0, 0, 0). The caller
 *  is responsible for centring the bond such that its site-symmetry
 *  is realised by elements fixing the origin. (Equivalently: pass a
 *  bond whose midpoint sits at the centre of the point-group axis
 *  system.) */
IRREP_API int irrep_dmi_allowed_basis_from_pg(const double r_a[3], const double r_b[3],
                                              const irrep_pg_table_t *pg, double tol,
                                              double basis_out[9]);

/** @brief Compute the symmetry-allowed *symmetric* part of the bond
 *         exchange tensor.
 *
 *  The full bilinear exchange Hamiltonian on a bond is
 *
 *  \f[ H_{ij} = S_i \cdot J_{ij} \cdot S_j \f]
 *
 *  with `J_ij` a real 3×3 matrix. Decompose `J = J^s + J^a` into
 *  symmetric and antisymmetric parts. The antisymmetric part is the
 *  DMI vector via `J^a_{ab} = ε_{abc} D^c` (handled by
 *  #irrep_dmi_allowed_basis). The symmetric part `J^s` carries the
 *  Heisenberg coupling (proportional to `I`), the dipolar /
 *  pseudodipolar anisotropy, and "Kitaev-Γ" type couplings — six
 *  components in general.
 *
 *  Under a symmetry operation `g` the symmetric tensor transforms as
 *  `J^s → R_g · J^s · R_g^T` for both bond-preserving AND bond-
 *  reversing operations (since `J^s = (J^s)^T`, the bond-reversal
 *  flip gives the same constraint as preservation). The allowed
 *  subspace is the joint kernel of the projectors
 *  `P_g = (R_g ⊗ R_g)(symmetrised)` averaged over the bond's site-
 *  symmetry.
 *
 *  @param r_a, r_b   Bond endpoints (cartesian).
 *  @param ops, n_ops List of candidate symmetry operations.
 *  @param tol        Position-matching tolerance.
 *  @param basis_out  Caller buffer of `6 × 9 = 54` doubles. On return,
 *                    the first `n × 9` entries hold `n` orthonormal
 *                    3×3 row-major symmetric matrices spanning the
 *                    allowed subspace under the Frobenius inner
 *                    product.
 *  @return Dimension of the allowed J^s subspace (0..6). */
IRREP_API int irrep_exchange_symmetric_basis(const double r_a[3], const double r_b[3],
                                             const irrep_dmi_sym_op_t *ops, int n_ops,
                                             double tol, double basis_out[54]);

/** @brief Convenience wrapper using a libirrep point-group table. */
IRREP_API int irrep_exchange_symmetric_basis_from_pg(const double r_a[3], const double r_b[3],
                                                     const irrep_pg_table_t *pg, double tol,
                                                     double basis_out[54]);

/** @brief Test whether the three-spin scalar chirality
 *         `χ_{ijk} = S_i · (S_j × S_k)` is symmetry-allowed on a triangle.
 *
 *  The scalar chirality is the next-order term beyond DMI in the
 *  spin-orbit-coupled hierarchy. It transforms as a **pseudoscalar**:
 *  invariant under proper rotations, sign-flipped under improper
 *  rotations and under time reversal. On a triangle `{r_a, r_b, r_c}`
 *  with site stabiliser `S`, `χ` is allowed iff every operation
 *  `g ∈ S` satisfies
 *
 *      `det(g) · σ_{perm}(g) · σ_T(g) = +1`,
 *
 *  where `σ_{perm} = ±1` is the sign of the induced permutation of
 *  `(r_a, r_b, r_c)` (cyclic = +1, anticyclic = −1), and
 *  `σ_T = -1 for antiunitary` (`T` flips a pseudoscalar) or `+1` for
 *  spatial-only.
 *
 *  Standard textbook consequences (drives topological-Hall responses
 *  in Mn₃Sn / Mn₃Ge / scalar-chiral kagome QSL candidates):
 *
 *    - Mirror **in** the triangle plane (σ_h, fixes all 3 sites,
 *      σ_perm = +1, det = −1) → product = −1 → χ = 0. The plane
 *      makes "above" and "below" indistinguishable.
 *    - Pure C_3 rotation about the triangle centroid (σ_perm = +1
 *      cyclic, det = +1) → product = +1 → trivial; χ allowed if no
 *      other element forbids it.
 *    - Pure spatial inversion of the triangle (no σ_h): permutes
 *      `(a, b, c) → ((-a, -b, -c))` and sends positions through the
 *      origin. For a triangle whose centroid is at the origin and
 *      that is invariant under inversion (e.g., a regular triangle
 *      paired with its inverted partner), the permutation involved
 *      depends on geometry. Caller responsibility to encode the
 *      operation correctly.
 *
 *  @param r_a, r_b, r_c  Cartesian positions of the triangle vertices.
 *  @param ops, n_ops     Candidate symmetry operations.
 *  @param tol            Position-matching tolerance.
 *  @return `1` if scalar chirality is symmetry-allowed,
 *          `0` if forbidden by some element of the stabiliser,
 *          `-1` on invalid input. */
IRREP_API int irrep_chirality_allowed(const double r_a[3], const double r_b[3],
                                      const double r_c[3], const irrep_dmi_sym_op_t *ops,
                                      int n_ops, double tol);

/** @brief Convenience wrapper for #irrep_chirality_allowed using a
 *         libirrep point-group table (all elements taken as unitary). */
IRREP_API int irrep_chirality_allowed_from_pg(const double r_a[3], const double r_b[3],
                                              const double r_c[3], const irrep_pg_table_t *pg,
                                              double tol);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_DMI_H */
