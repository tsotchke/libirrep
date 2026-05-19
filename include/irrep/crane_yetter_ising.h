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

/** @brief Verify the F-matrix unitarity on the σσσσ block:
 *  `Σ_e F^{σσσ}_σ_{e,f} · conj(F^{σσσ}_σ_{e,f'}) = δ_{f,f'}`.
 *
 *  @return Maximum absolute deviation from identity. */
IRREP_API double
irrep_ising_F_unitarity_residual(void);

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
