/* SPDX-License-Identifier: MIT */
/** @file toric_mtc.h
 *  @brief Z₂ × Z₂ toric-code modular tensor category.
 *
 *  The MTC underlying Kitaev's 2D toric code — the simplest non-trivial
 *  abelian MTC, with 4 simples {1, e, m, ψ} forming the abelian group
 *  Z₂ × Z₂ under fusion. Companion to `<irrep/crane_yetter_ising.h>`
 *  which provides the non-abelian Ising MTC.
 *
 *  ## Simples and fusion
 *
 *  Four simples, all with quantum dimension 1 (abelian):
 *    - `1` (vacuum)
 *    - `e` (electric charge / vertex defect)
 *    - `m` (magnetic flux / plaquette defect)
 *    - `ψ` = e·m (fermion / dyon)
 *
 *  Fusion rules (all multiplicities 1):
 *    `e × e = 1,  m × m = 1,  ψ × ψ = 1`
 *    `e × m = ψ,  e × ψ = m,  m × ψ = e`
 *
 *  All fusions are commutative; the simples form Z₂ × Z₂ under tensor.
 *
 *  ## Modular data
 *
 *  Global dimension: `D = √(1 + 1 + 1 + 1) = 2`. Central charge: `c = 0`
 *  (the Z₂×Z₂ toric code is a non-chiral / "doubled" topological order;
 *  the boundary CFT is trivial / has c = 0).
 *
 *  S-matrix (character table of Z₂ × Z₂, normalised so S is unitary):
 *
 *      S = (1/2) [[ 1,  1,  1,  1],
 *                 [ 1,  1, -1, -1],
 *                 [ 1, -1,  1, -1],
 *                 [ 1, -1, -1,  1]]
 *
 *  T-matrix (topological twists):
 *    `T_1 = 1, T_e = 1, T_m = 1, T_ψ = -1`
 *  (e and m are bosons; ψ is a fermion.)
 *
 *  R-symbols (braiding, F-symbols all = 1 since the MTC is abelian):
 *    Self-braids: R^{ee}_1 = R^{mm}_1 = 1, R^{ψψ}_1 = -1.
 *    Mutual braids: R^{em}_ψ = R^{me}_ψ = i (their product = -1, the
 *      mutual-monodromy phase characteristic of e ↔ m anyons).
 *    Derived: R^{eψ}_m = i, R^{ψe}_m = i, R^{mψ}_e = i, R^{ψm}_e = i.
 *
 *  ## Primary references
 *
 *  - Kitaev, *Fault-tolerant quantum computation by anyons*, Annals
 *    Phys. 303 (2003) 2 — the 2D toric code and its Z₂×Z₂ topological
 *    order.
 *  - Kitaev, *Anyons in an exactly solved model and beyond*, Annals
 *    Phys. 321 (2006) 2 — MTC framework, see Appendix B for Z₂×Z₂
 *    explicit modular data.
 */
#ifndef IRREP_TORIC_MTC_H
#define IRREP_TORIC_MTC_H

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Z₂ × Z₂ toric-code MTC simple objects. */
typedef enum {
    IRREP_TORIC_MTC_OBJ_1   = 0, /**< Vacuum / identity. */
    IRREP_TORIC_MTC_OBJ_E   = 1, /**< Electric charge (boson). */
    IRREP_TORIC_MTC_OBJ_M   = 2, /**< Magnetic flux (boson). */
    IRREP_TORIC_MTC_OBJ_PSI = 3, /**< Dyon ψ = e·m (fermion). */
} irrep_toric_mtc_object_t;

/** @brief Number of simples in the Z₂×Z₂ toric MTC. */
#define IRREP_TORIC_MTC_N_OBJECTS 4

/** @brief Quantum dimension d_a (= 1 for every simple of an abelian MTC). */
IRREP_API double
irrep_toric_mtc_quantum_dim(irrep_toric_mtc_object_t a);

/** @brief Global quantum dimension D = √(Σ d_a²) = 2. */
IRREP_API double
irrep_toric_mtc_global_dim(void);

/** @brief Central charge c = 0 (non-chiral / doubled topological order). */
IRREP_API double
irrep_toric_mtc_central_charge(void);

/** @brief Fusion coefficient N^{ab}_c ∈ {0, 1}. */
IRREP_API int
irrep_toric_mtc_fusion(irrep_toric_mtc_object_t a,
                      irrep_toric_mtc_object_t b,
                      irrep_toric_mtc_object_t c);

/** @brief Modular S-matrix entry. Unitary, normalised so S† S = I.
 *
 *      S = (1/2) [[1, 1, 1, 1], [1, 1, -1, -1],
 *                 [1, -1, 1, -1], [1, -1, -1, 1]]
 */
IRREP_API double _Complex
irrep_toric_mtc_S_matrix(irrep_toric_mtc_object_t a,
                        irrep_toric_mtc_object_t b);

/** @brief Modular T-matrix (topological twist) T_a = exp(2πi h_a):
 *
 *      h_1 = h_e = h_m = 0,  h_ψ = 1/2  →  T = diag(1, 1, 1, -1). */
IRREP_API double _Complex
irrep_toric_mtc_T_eigenvalue(irrep_toric_mtc_object_t a);

/** @brief F-symbol [F^{abc}_d]_{e,f} = 1 if allowed by fusion, 0 otherwise.
 *
 *  All non-zero F-symbols of an abelian MTC equal 1; there are no
 *  Hadamard-style mixers as in Ising. */
IRREP_API double _Complex
irrep_toric_mtc_F_symbol(irrep_toric_mtc_object_t a, irrep_toric_mtc_object_t b,
                         irrep_toric_mtc_object_t c, irrep_toric_mtc_object_t d,
                         irrep_toric_mtc_object_t e, irrep_toric_mtc_object_t f);

/** @brief R-symbol R^{ab}_c — the abelian braiding phase. */
IRREP_API double _Complex
irrep_toric_mtc_R_symbol(irrep_toric_mtc_object_t a,
                         irrep_toric_mtc_object_t b,
                         irrep_toric_mtc_object_t c);

/* ====================================================================
 * Consistency proofs (mirror the Ising-MTC layout).
 * ==================================================================== */

/** @brief Verify modular S² = I (self-dual MTC ⇒ charge conjugation = I).
 *  Returns max |S²_{ab} - δ_{ab}|. */
IRREP_API double
irrep_toric_mtc_S_squared_residual(void);

/** @brief Verify Verlinde formula `N^{ab}_c = Σ_x S_{ax} S_{bx} S*_{cx} / S_{0x}`
 *  across all 64 (a, b, c) triples. */
IRREP_API double
irrep_toric_mtc_verlinde_residual(void);

/** @brief Verify topological twist from R-symbol: θ_a = R^{aa}_1 for an
 *  abelian MTC (single-channel self-braid). Returns max |θ_a - T_a|. */
IRREP_API double
irrep_toric_mtc_twist_from_R_residual(void);

/* ====================================================================
 * Walker-Wang admissibility for Z₂ × Z₂.
 * ==================================================================== */

/** @brief Admissibility of a multiset of Z₂×Z₂ labels.
 *
 *  For an abelian MTC, the multiset {a_1, ..., a_n} fuses to vacuum
 *  iff their product (under Z₂×Z₂ multiplication) is the identity.
 *  Equivalently: viewing each simple as a pair (a_e, a_m) ∈ Z₂²,
 *  the multiset is admissible iff Σ a_e and Σ a_m are both even. */
IRREP_API int
irrep_toric_mtc_admissible(const irrep_toric_mtc_object_t *labels, int n);

/** @brief Walker-Wang ground-state dim on the 3-simplex (tetrahedron)
 *  for the Z₂×Z₂ toric-code MTC, vertex + face constrained.
 *
 *  Counts edge labelings of the 6-edge tetrahedron (4^6 = 4096
 *  configurations) for which every vertex's 3 incident edges AND
 *  every triangular face's 3 boundary edges fuse to vacuum.
 *
 *  Direct companion to `irrep_ising_walker_wang_simplex3_full_count`
 *  (= 16 for Ising); this gives the corresponding number for the
 *  abelian toric MTC. Different MTC, different count — showing that
 *  the WW ground-state dim depends on both the topology AND the
 *  underlying anyon model. */
IRREP_API long long
irrep_toric_mtc_walker_wang_simplex3_full_count(void);

/** @brief Walker-Wang ground-state dim on the unit cube
 *  (8 vertices 3-valent, 12 edges, 6 square faces) for Z₂×Z₂.
 *
 *  Enumerates 4¹² = 16,777,216 edge labelings under all 8 vertex
 *  + 6 face admissibility constraints. Runs in ~50 ms with
 *  early-out vertex pruning.
 *
 *  Cross-MTC fact: this differs from
 *  `irrep_ising_walker_wang_cube_full_count()` (= 120), proving
 *  that the WW state-sum depends on the underlying anyon model
 *  beyond just the global dimension D once the geometry is
 *  non-trivial enough. */
IRREP_API long long
irrep_toric_mtc_walker_wang_cube_full_count(void);

/* ====================================================================
 * Crane-Yetter 4-manifold invariant for the Z₂ × Z₂ MTC.
 * ==================================================================== */

/** @brief Crane-Yetter 4-manifold invariant for the Z₂×Z₂ MTC.
 *
 *  Closed-form for any MTC: `Z(M) = D^{-χ(M)} · exp(2πi c σ(M) / 8)`.
 *  For Z₂×Z₂: D = 2, c = 0, so the phase vanishes:
 *
 *      Z_Z2Z2(M) = 2^{-χ(M)}   (real, no σ-dependence)
 *
 *  This is simpler than the Ising case (where ζ = exp(iπ/8) gives the
 *  σ-twist phase). Distinguishes manifolds purely by Euler char χ. */
IRREP_API double _Complex
irrep_toric_mtc_invariant(int euler_char, int signature);

/** @brief Verify Z₂×Z₂ CY connected-sum multiplicativity at runtime.
 *
 *  Same axiom as the Ising version (`irrep_crane_yetter_connected_sum_residual`):
 *
 *      Z(M # N) = Z(M) · Z(N) / Z(S⁴)
 *
 *  For Z₂×Z₂, Z(S⁴) = 2^{-2} = 1/4, so Z(M#N) = 4·Z(M)·Z(N) — same
 *  factor as Ising (both have D = 2). Verified on the same 5-pair
 *  test set. Distinct from Ising's residual since Z₂×Z₂'s Z is real
 *  while Ising's is complex (ζ-phase).
 *
 *  @return Max absolute deviation over the test set; < 1e-12. */
IRREP_API double
irrep_toric_mtc_connected_sum_residual(void);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_TORIC_MTC_H */
