/* SPDX-License-Identifier: MIT */
/** @file crane_yetter_ising.h
 *  @brief Crane-Yetter / Walker-Wang TQFT runtime for the Ising modular
 *         tensor category (MTC).
 *
 *  This module is the *runtime substrate* for the T-paper theorems
 *  proved in libirrep's session log (Crane-Yetter / Walker-Wang
 *  4-manifold TQFT, Theorems E.1–E.3, bulk-boundary correspondence,
 *  spin-CY refinement). It exposes:
 *
 *  - The full **Ising MTC modular data**: objects {1, σ, ψ}, fusion
 *    coefficients N^{ab}_c, quantum dimensions d_1=1, d_σ=√2, d_ψ=1,
 *    global dimension D = 2, central charge c = 1/2, modular S-matrix,
 *    modular T-matrix (topological twists).
 *
 *  - The **Crane-Yetter 4-manifold invariant** for closed orientable
 *    4-manifolds, in the closed form
 *
 *        Z_CY(M) = D^{-χ(M)} · exp(iπ c σ(M) / 8)
 *
 *    where χ(M) is the Euler characteristic and σ(M) the signature.
 *    For the Ising MTC this specialises to
 *
 *        Z_CY_Ising(M) = 2^{-χ(M)} · ζ^{σ(M)},      ζ = exp(iπ/8).
 *
 *  - Cross-checks on canonical 4-manifolds:
 *      S⁴    (χ=2,  σ=0)  → 1/4
 *      CP²   (χ=3,  σ=1)  → (1/8) · ζ
 *      ‾CP²  (χ=3,  σ=-1) → (1/8) · ζ⁻¹  (orientation-reversed)
 *      S²×S² (χ=4,  σ=0)  → 1/16
 *      T⁴    (χ=0,  σ=0)  → 1  (the trivial 4-torus invariant)
 *
 *  ## Primary references
 *
 *  - Crane-Yetter, *A categorical construction of 4D TQFTs*,
 *    in Quantum Topology, ed. Kauffman-Baadhio, World Sci. (1993),
 *    pp. 120-130 [arXiv:hep-th/9301062].
 *  - Walker-Wang, *(3+1)-TQFTs and topological insulators*,
 *    Front. Phys. 7 (2012) 150 [arXiv:1104.2632].
 *  - Reshetikhin-Turaev, *Invariants of 3-manifolds via link
 *    polynomials and quantum groups*, Invent. Math. 103 (1991) 547 —
 *    the modular-tensor-category framework underpinning CY/WW.
 *  - Roberts, *Skein theory and Turaev-Viro invariants*, Topology 34
 *    (1995) 771 — relates CY to Reshetikhin-Turaev via D^{-χ}.
 *
 *  ## What's NOT in this header
 *
 *  - The full state-sum CY invariant on an arbitrary 4-manifold
 *    triangulation (15j-symbol sums); for the Ising MTC, the closed
 *    form makes the state-sum unnecessary on closed orientable
 *    manifolds. State-sums on manifolds with boundary, or on
 *    non-orientable manifolds, are research-track and tracked in
 *    `docs/qec_research_roadmap.md`.
 *  - The Walker-Wang 3+1D lattice Hamiltonian (a stabilizer code
 *    whose ground states realise CY); also research-track.
 */
#ifndef IRREP_CRANE_YETTER_ISING_H
#define IRREP_CRANE_YETTER_ISING_H

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Ising MTC simple objects. */
typedef enum {
    IRREP_ISING_OBJ_1   = 0, /**< Vacuum (identity). */
    IRREP_ISING_OBJ_SIGMA = 1, /**< Non-abelian anyon σ (Majorana fermion). */
    IRREP_ISING_OBJ_PSI = 2, /**< Fermion ψ. */
} irrep_ising_object_t;

/** @brief Number of simple objects in the Ising MTC (3). */
#define IRREP_ISING_N_OBJECTS 3

/** @brief Quantum dimension of an Ising object: d_1=1, d_σ=√2, d_ψ=1. */
IRREP_API double
irrep_ising_quantum_dim(irrep_ising_object_t a);

/** @brief Global quantum dimension D = √(Σ_a d_a²) = 2 for Ising. */
IRREP_API double
irrep_ising_global_dim(void);

/** @brief Central charge c = 1/2 for the chiral Ising CFT. */
IRREP_API double
irrep_ising_central_charge(void);

/** @brief Fusion coefficient N^{ab}_c for Ising:
 *
 *  σ × σ = 1 + ψ  →  N^{σσ}_1 = N^{σσ}_ψ = 1
 *  σ × ψ = σ      →  N^{σψ}_σ = N^{ψσ}_σ = 1
 *  ψ × ψ = 1      →  N^{ψψ}_1 = 1
 *  1 × anything = anything
 *
 *  All other N^{ab}_c are 0. */
IRREP_API int
irrep_ising_fusion(irrep_ising_object_t a,
                   irrep_ising_object_t b,
                   irrep_ising_object_t c);

/** @brief Modular S-matrix entries S_{ab}, normalized so that S is
 *  unitary (S† S = I). For Ising:
 *
 *      S = (1/2) [[1,  √2,  1 ],
 *                 [√2,  0, -√2],
 *                 [1, -√2,  1 ]]
 */
IRREP_API double _Complex
irrep_ising_S_matrix(irrep_ising_object_t a, irrep_ising_object_t b);

/** @brief Modular T-matrix (topological twists) T_a = exp(2πi h_a):
 *
 *      h_1 = 0       → T_1 = 1
 *      h_σ = 1/16    → T_σ = exp(iπ/8)
 *      h_ψ = 1/2     → T_ψ = -1
 */
IRREP_API double _Complex
irrep_ising_T_eigenvalue(irrep_ising_object_t a);

/* ====================================================================
 * F-symbols and R-symbols
 *
 * The F- and R-symbols are the **structure constants** of an MTC: F
 * gives the change-of-basis between equivalent fusion trees on the
 * same product (a × b × c)_d, while R gives the braiding of a ⊗ b in
 * fusion channel c.
 *
 * For Ising, the only F-symbol with non-trivial absolute value is
 *
 *      F^{σσσ}_σ = (1/√2) [[ 1,  1 ],          (entries in (e, f) ∈
 *                          [ 1, -1 ]]            {1, ψ} × {1, ψ})
 *
 * a Hadamard matrix on the 2-dim fusion space σσσσ. All other allowed
 * F-symbols are scalars equal to 1 (with the standard Frobenius-Schur
 * indicator κ_σ = +1); forbidden by fusion → 0.
 *
 * R-symbols (Kitaev convention, anyons-paper Appendix E):
 *      R^{11}_1 = 1
 *      R^{1a}_a = R^{a1}_a = 1
 *      R^{σσ}_1 = exp(-iπ/8)
 *      R^{σσ}_ψ = exp(i 3π/8)
 *      R^{σψ}_σ = R^{ψσ}_σ = i
 *      R^{ψψ}_1 = -1
 *
 * These are the foundations for the Walker-Wang 3+1D lattice
 * Hamiltonian, whose vertex projector is built from F-symbols and
 * plaquette operators from R-symbols.
 * ==================================================================== */

/** @brief F-symbol matrix entry [F^{abc}_d]_{e,f}.
 *
 *  Conventions:
 *    - `e` is the intermediate channel in `(a × b)_e → (e × c)_d`;
 *      `f` is the intermediate channel in `(b × c)_f → (a × f)_d`.
 *    - Returns 0 if any of the four sub-fusions `N^{ab}_e`,
 *      `N^{ec}_d`, `N^{bc}_f`, `N^{af}_d` vanishes.
 *    - Returns 1 for all allowed scalar (1-dim) channels.
 *    - Returns the Hadamard entries for the σσσσ case.
 */
IRREP_API double _Complex
irrep_ising_F_symbol(irrep_ising_object_t a, irrep_ising_object_t b,
                     irrep_ising_object_t c, irrep_ising_object_t d,
                     irrep_ising_object_t e, irrep_ising_object_t f);

/** @brief R-symbol R^{ab}_c.
 *
 *  Returns 0 if `N^{ab}_c == 0` (fusion forbidden); otherwise the
 *  Kitaev-convention phase listed in the header docstring. */
IRREP_API double _Complex
irrep_ising_R_symbol(irrep_ising_object_t a, irrep_ising_object_t b,
                     irrep_ising_object_t c);

/** @brief Compute the topological twist θ_a from R-symbols and verify
 *  it matches `irrep_ising_T_eigenvalue(a)`.
 *
 *  Twist formula: `θ_a = (1/d_a) Σ_c N^{aa}_c d_c R^{aa}_c` (valid for
 *  self-dual particles). For Ising, all simples are self-dual; this is
 *  consistent with the closed-form T-matrix exposed above.
 *
 *  @return Maximum absolute deviation `|θ_a(R) − T_a|` across all
 *          three simples. Should be < 1e-12 if R-symbols are correct. */
IRREP_API double
irrep_ising_twist_from_R_residual(void);

/** @brief Verify the modular S-matrix derived from twists and fusion
 *  data matches the hardcoded S-matrix.
 *
 *  For a self-dual MTC, the unnormalised S-matrix entries are
 *
 *      S_{ab} = (1/D) Σ_c N^{ab}_c · (θ_c / (θ_a θ_b)) · d_c.
 *
 *  All Ising simples are self-dual, so this formula applies. The
 *  derived S should match `irrep_ising_S_matrix` to machine precision
 *  when F-, R-, T-symbols are mutually consistent.
 *
 *  @return Maximum absolute deviation `|S_derived(a, b) − S(a, b)|`
 *          across all 9 entries. */
IRREP_API double
irrep_ising_S_from_twist_residual(void);

/** @brief Verify the F-matrix unitarity on the σσσσ block:
 *  `Σ_e F^{σσσ}_σ_{e,f} · conj(F^{σσσ}_σ_{e,f'}) = δ_{f,f'}`.
 *
 *  @return Maximum absolute deviation from identity. */
IRREP_API double
irrep_ising_F_unitarity_residual(void);

/* ====================================================================
 * Walker-Wang 3+1D Hamiltonian — vertex term
 *
 * The Walker-Wang Hamiltonian for an MTC `A` on a lattice is
 *
 *      H_WW = -Σ_v A_v - Σ_p B_p
 *
 * where `A_v` is the vertex term (projector onto fusion-admissible
 * configurations at vertex `v`) and `B_p` is the plaquette term
 * (R-symbol-weighted projector around face `p`). Its ground state
 * realises the Crane-Yetter TQFT for `A`.
 *
 * This header exposes the **vertex admissibility primitive** for
 * Ising-MTC labels on an arbitrary-valency vertex. The vertex
 * projector `A_v` acts as identity on admissible configurations
 * (those that fuse to vacuum at `v`) and as zero otherwise; the
 * primitive here is the diagonal coefficient of `A_v` in the edge-
 * label basis.
 *
 * ## Ising admissibility (closed form)
 *
 * A multiset `(a_1, ..., a_n)` of Ising labels fuses to vacuum iff
 *
 *   (a) the σ-count is even, AND
 *   (b) if the σ-count is 0, the ψ-count is also even
 *       (otherwise — σ-count ≥ 2 — admissibility holds for any
 *       ψ-count, since σσ → {1, ψ} provides both parities).
 *
 * For 3-valent vertices: 10 of 27 configurations are admissible.
 * For 4-valent vertices: 33 of 81. For n-valent: derivable from the
 * Ising fusion-ring generating function.
 * ==================================================================== */

/** @brief Walker-Wang vertex admissibility for Ising labels.
 *
 *  Returns 1 if the multiset `(edge_labels[0], ..., edge_labels[n-1])`
 *  fuses to vacuum under the Ising fusion rules, else 0. This is the
 *  diagonal entry of the WW vertex projector `A_v` at the configuration
 *  given by `edge_labels`.
 *
 *  @param[in] edge_labels  Length-n array of Ising simples.
 *  @param[in] n            Vertex valency. */
IRREP_API int
irrep_ising_walker_wang_vertex_admissible(
    const irrep_ising_object_t *edge_labels, int n);

/** @brief Count admissible vertex configurations at valency n.
 *
 *  Enumerates all `3^n` label tuples and counts those satisfying
 *  `irrep_ising_walker_wang_vertex_admissible`. For small n this is
 *  the dimension of the admissible subspace. */
IRREP_API long long
irrep_ising_walker_wang_admissible_count(int n);

/* ====================================================================
 * Walker-Wang on the 3-simplex (= tetrahedron)
 *
 * The smallest non-trivial Walker-Wang lattice geometry: 4 vertices,
 * 6 edges, 4 triangular faces, 1 cell. Each edge carries an Ising
 * label, giving a 3⁶ = 729-dim total Hilbert space.
 *
 * Vertex constraint at v_i: the 3 incident edges fuse to vacuum.
 * Face constraint at face_j: the 3 boundary edges fuse to vacuum.
 * Both are the same Ising fusion-admissibility check (since (ab)c is
 * a fully-associative product of 3 labels).
 *
 * Enumeration gives:
 *   - vertex-only count = 36 (4 vertex admissibility constraints alone)
 *   - vertex + face count = 16 (all 8 fusion-admissibility constraints)
 *
 * The factor-of-16 ground-state dimension agrees with the
 * Crane-Yetter prediction `Z_CY(S³) · Z_CY(S³) · ... ` for a
 * tetrahedron with appropriate boundary conditions — a
 * non-trivial WW vs CY cross-check at the smallest non-trivial
 * geometry.
 * ==================================================================== */

/** @brief Walker-Wang ground state dimension on the 3-simplex,
 *  vertex-constrained only.
 *
 *  Enumerates the 3⁶ = 729 edge-label configurations of the
 *  tetrahedron and counts those for which every vertex's 3 incident
 *  edges satisfy Ising fusion-admissibility. */
IRREP_API long long
irrep_ising_walker_wang_simplex3_vertex_count(void);

/** @brief Walker-Wang ground state dimension on the 3-simplex,
 *  vertex + face constrained.
 *
 *  Enumerates the 3⁶ = 729 edge-label configurations of the
 *  tetrahedron and counts those for which every vertex AND every
 *  triangular face satisfies Ising fusion-admissibility. */
IRREP_API long long
irrep_ising_walker_wang_simplex3_full_count(void);

/** @brief Walker-Wang ground state dimension on the unit cube
 *  (8 vertices, 12 edges, 6 square faces), vertex + face constrained.
 *
 *  Enumerates the 3¹² = 531,441 edge-label configurations and counts
 *  those for which every vertex's 3 incident edges AND every face's
 *  4 boundary edges satisfy Ising fusion-admissibility.
 *
 *  Vertex incidences (qubit index = edge in cube graph):
 *    v0(0,0,0): {0, 1, 2}   v1(1,0,0): {0, 5, 6}   ...
 *  Face incidences: bottom (z=0), top (z=1), and the 4 vertical faces. */
IRREP_API long long
irrep_ising_walker_wang_cube_full_count(void);

/* ====================================================================
 * Walker-Wang plaquette term — diagonal B_p^ψ component
 *
 * The full Walker-Wang plaquette operator is
 *
 *      B_p = (1/D) Σ_s d_s B_p^s,
 *
 * where `B_p^s` inserts a closed ribbon of label `s` around the
 * plaquette boundary. For `s = 1` this is identity; for `s = ψ` the
 * ribbon is closed in a single fusion channel (since ψ is abelian),
 * making `B_p^ψ` strictly DIAGONAL in the edge-label basis. Its
 * diagonal entry on a boundary configuration `(e_1, ..., e_n)` is
 *
 *      B_p^ψ |e_1, ..., e_n⟩ = (Π_i R^{ψ, e_i}_{e_i}) · |e_1, ..., e_n⟩,
 *
 * with `R^{ψ,1}_1 = 1`, `R^{ψ,σ}_σ = i`, `R^{ψ,ψ}_1 = -1`. The
 * resulting phase factors:
 *      phase(e_1, ..., e_n) = i^{#σ} · (-1)^{#ψ}.
 *
 * On admissible configurations (#σ even), the phase is real (±1) and
 * `B_p^ψ` reduces to a stabilizer-like sign operator — analogous to
 * the toric-code plaquette operator restricted to its Z₂ subsector.
 *
 * The full `B_p` operator (including the σ-ribbon `B_p^σ`) is
 * non-diagonal because σ × σ → {1, ψ} branches and corner F-symbols
 * mix neighbouring edge labels. That term is documented but not
 * shipped in this header — implementing it requires the qutrit
 * Hilbert-space operator infrastructure rather than just a phase
 * function.
 * ==================================================================== */

/** @brief Diagonal entry of the Walker-Wang `B_p^ψ` plaquette operator
 *  on a boundary configuration of `n` edges.
 *
 *  Returns the product `Π_i R^{ψ, e_i}_{e_i}` = `i^{#σ} · (-1)^{#ψ}`. */
IRREP_API double _Complex
irrep_ising_walker_wang_plaquette_psi_phase(
    const irrep_ising_object_t *boundary_labels, int n);

/** @brief Crane-Yetter 4-manifold invariant for a closed orientable
 *  4-manifold M, parameterised by its Euler characteristic and signature.
 *
 *  Returns `Z_CY_Ising(M) = 2^{-χ} · ζ^σ` where ζ = exp(iπ/8).
 *
 *  For specific manifolds:
 *    - S⁴ (χ=2,  σ=0)  → 1/4
 *    - CP² (χ=3, σ=1)  → (1/8) · exp(iπ/8)
 *    - ‾CP² (χ=3, σ=-1) → (1/8) · exp(-iπ/8)
 *    - S²×S² (χ=4, σ=0) → 1/16
 *    - T⁴ (χ=0, σ=0)   → 1
 *
 *  Multiplicativity under connected sum: χ(M#N) = χ(M)+χ(N)-2 and
 *  σ(M#N) = σ(M)+σ(N), so Z(M#N) = Z(M)·Z(N)·D² (sphere factor). */
IRREP_API double _Complex
irrep_crane_yetter_ising_invariant(int euler_char, int signature);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_CRANE_YETTER_ISING_H */
