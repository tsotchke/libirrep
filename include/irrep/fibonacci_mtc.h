/* SPDX-License-Identifier: MIT */
/** @file fibonacci_mtc.h
 *  @brief Fibonacci anyon modular tensor category.
 *
 *  The simplest non-abelian MTC with a universal braiding model:
 *  two simples {1, τ} and the single non-trivial fusion τ × τ = 1 + τ.
 *  Companion to `<irrep/crane_yetter_ising.h>` (Ising) and
 *  `<irrep/toric_mtc.h>` (Z₂×Z₂).
 *
 *  ## Simples and fusion
 *
 *  Two simples:
 *    - `1` (vacuum)
 *    - `τ` (Fibonacci anyon — non-abelian)
 *
 *  Single non-trivial fusion: `τ × τ = 1 + τ`. The fusion ring is
 *  Z[φ]/(φ² - φ - 1) where φ is the golden ratio, hence the name.
 *  Quantum dim: `d_τ = φ = (1 + √5) / 2 ≈ 1.618`.
 *  Global dim: `D² = 1 + φ² = 2 + φ`, so `D = √(2 + φ)`.
 *  Central charge: `c = 14/5` (chiral, this is the (G₂)₁ Yang-Lee
 *  minimal-model CFT central charge).
 *
 *  ## Modular data
 *
 *  S-matrix:
 *      S = (1/D) [[ 1,   φ  ],
 *                 [ φ,  -1  ]]
 *  with `D = √(2 + φ)`.
 *
 *  T-matrix (twists): `T_1 = 1`, `T_τ = exp(4πi/5)` (= θ_τ where
 *  `h_τ = 2/5`).
 *
 *  F-symbols: only `F^{τττ}_τ` is non-trivial, a 2×2 matrix in the
 *  `(1, τ)` intermediate-channel basis:
 *
 *      F^{τττ}_τ = [[ φ⁻¹,   φ⁻¹/²  ],
 *                   [ φ⁻¹/², -φ⁻¹   ]]
 *
 *  R-symbols (only non-trivial ones):
 *      R^{ττ}_1 = exp(-4πi/5)
 *      R^{ττ}_τ = exp(3πi/5)
 *
 *  ## Why this matters
 *
 *  Fibonacci anyons are the simplest model supporting UNIVERSAL
 *  topological quantum computation via braiding alone: every
 *  unitary gate on n logical qubits can be approximated to any
 *  desired precision using only braiding of Fibonacci anyons.
 *  Ising anyons, by contrast, give only a Clifford-restricted gate
 *  set without additional magic-state distillation.
 *
 *  ## Primary references
 *
 *  - Trebst-Troyer-Wang-Ludwig, *A short introduction to Fibonacci
 *    anyon models*, Prog. Theor. Phys. Suppl. 176 (2008) 384.
 *  - Kitaev, *Anyons in an exactly solved model and beyond*, Annals
 *    Phys. 321 (2006) 2 — MTC framework, Fibonacci as the canonical
 *    universal non-abelian example.
 *  - Nayak-Simon-Stern-Freedman-Das Sarma, *Non-Abelian anyons and
 *    topological quantum computation*, Rev. Mod. Phys. 80 (2008) 1083.
 */
#ifndef IRREP_FIBONACCI_MTC_H
#define IRREP_FIBONACCI_MTC_H

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Fibonacci MTC simple objects. */
typedef enum {
    IRREP_FIB_OBJ_1   = 0, /**< Vacuum. */
    IRREP_FIB_OBJ_TAU = 1, /**< Fibonacci anyon τ. */
} irrep_fib_object_t;

/** @brief Number of simples in the Fibonacci MTC (2). */
#define IRREP_FIB_N_OBJECTS 2

/** @brief Golden ratio φ = (1 + √5)/2. */
IRREP_API double irrep_fib_golden_ratio(void);

/** @brief Quantum dimension: d_1 = 1, d_τ = φ. */
IRREP_API double irrep_fib_quantum_dim(irrep_fib_object_t a);

/** @brief Global dim D = √(2 + φ). */
IRREP_API double irrep_fib_global_dim(void);

/** @brief Central charge c = 14/5 (chiral). */
IRREP_API double irrep_fib_central_charge(void);

/** @brief Fusion coefficient N^{ab}_c.
 *
 *      1 × 1 = 1,  1 × τ = τ,  τ × 1 = τ,  τ × τ = 1 + τ. */
IRREP_API int
irrep_fib_fusion(irrep_fib_object_t a, irrep_fib_object_t b,
                 irrep_fib_object_t c);

/** @brief Modular S-matrix entry S_{ab}.
 *
 *      S = (1/D) [[1, φ], [φ, -1]]. */
IRREP_API double _Complex
irrep_fib_S_matrix(irrep_fib_object_t a, irrep_fib_object_t b);

/** @brief Topological twist T_a = exp(2πi h_a):
 *      T_1 = 1,  T_τ = exp(4πi/5)  (h_τ = 2/5). */
IRREP_API double _Complex
irrep_fib_T_eigenvalue(irrep_fib_object_t a);

/** @brief F-symbol [F^{abc}_d]_{e,f}.
 *
 *  Non-trivial only for `a = b = c = d = τ` where it is the 2×2
 *  Fibonacci F-matrix; otherwise 1 or 0. */
IRREP_API double _Complex
irrep_fib_F_symbol(irrep_fib_object_t a, irrep_fib_object_t b,
                   irrep_fib_object_t c, irrep_fib_object_t d,
                   irrep_fib_object_t e, irrep_fib_object_t f);

/** @brief R-symbol R^{ab}_c.
 *
 *      R^{ττ}_1 = exp(-4πi/5),  R^{ττ}_τ = exp(3πi/5),
 *      all other allowed R-symbols = 1. */
IRREP_API double _Complex
irrep_fib_R_symbol(irrep_fib_object_t a, irrep_fib_object_t b,
                   irrep_fib_object_t c);

/* ====================================================================
 * Consistency proofs.
 * ==================================================================== */

/** @brief Verify F-matrix unitarity on the τττ Hadamard-like 2×2 block. */
IRREP_API double
irrep_fib_F_unitarity_residual(void);

/** @brief Verify modular S² = C (= I for self-dual Fibonacci). */
IRREP_API double
irrep_fib_S_squared_residual(void);

/** @brief Verify Verlinde formula N^{ab}_c = sum_x S_ax S_bx conj(S_cx) / S_0x
 *  across all 8 (a, b, c) triples. */
IRREP_API double
irrep_fib_verlinde_residual(void);

/** @brief Verify topological twist from R: θ_a = (1/d_a) Σ_c N^{aa}_c
 *  d_c R^{aa}_c reproduces the T-matrix. */
IRREP_API double
irrep_fib_twist_from_R_residual(void);

/* ====================================================================
 * Walker-Wang fusion-admissibility primitive for Fibonacci.
 *
 * Multiset {a_1, ..., a_n} of Fibonacci labels fuses to vacuum iff
 * τ^(#τ in multiset) contains 1 with non-zero multiplicity.
 *
 * Closed form (F_k = Fibonacci numbers, F_0 = 0, F_1 = 1, F_2 = 1, ...):
 *   τ^n = F_{n-1} · 1 + F_n · τ
 *
 * So coefficient of 1 in τ^n is F_{n-1}, which is > 0 iff n - 1 ≥ 1,
 * i.e., n ≠ 1. Hence:
 *
 *   Admissibility: multiset fuses to vacuum  ⇔  #τ ≠ 1.
 *
 * (Slightly counterintuitive: #τ = 0 trivially admissible; #τ ≥ 2 also
 * admissible since fusion multiplicities give non-zero coefficient of 1;
 * only #τ = 1 is non-admissible.)
 * ==================================================================== */

/** @brief Fibonacci multiset admissibility. Returns 1 if `#τ ≠ 1`
 *  (= can fuse to vacuum), else 0. */
IRREP_API int
irrep_fib_admissible(const irrep_fib_object_t *labels, int n);

/** @brief Walker-Wang admissibility count on the 3-simplex
 *  (tetrahedron). 4 vertices (3-valent) + 4 triangular faces. */
IRREP_API long long
irrep_fib_walker_wang_simplex3_full_count(void);

/** @brief Walker-Wang admissibility count on the unit cube
 *  (8 vertices 3-valent, 12 edges, 6 square faces). */
IRREP_API long long
irrep_fib_walker_wang_cube_full_count(void);

/** @brief Walker-Wang admissibility count on the regular octahedron
 *  (6 vertices 4-valent, 12 edges, 8 triangular faces).
 *  Dual to the cube; should give the same count. */
IRREP_API long long
irrep_fib_walker_wang_octahedron_full_count(void);

/** @brief Walker-Wang admissibility count on the triangular bipyramid. */
IRREP_API long long
irrep_fib_walker_wang_tri_bipyramid_full_count(void);

/** @brief Walker-Wang admissibility count on the triangular prism
 *  (dual to the bipyramid). */
IRREP_API long long
irrep_fib_walker_wang_tri_prism_full_count(void);

/* ====================================================================
 * Fibonacci anyonic braiding — universal quantum computation via braids.
 *
 * For 4 Fibonacci anyons fused pairwise to vacuum, the fusion space has
 * dimension F_3 = 2 (the τ-channel fusion has multiplicity 2 from
 * τ × τ × τ × τ → 1). This gives a "logical qubit" with two basis
 * states determined by the intermediate fusion channels (1 or τ
 * for the first pair).
 *
 * In this basis, the braid generators σ_1, σ_2, σ_3 (braiding anyons
 * 1↔2, 2↔3, 3↔4) act as 2×2 unitaries:
 *
 *   σ_1 = R⊕R₂ (braid in-place, no F-move): diagonal in the basis,
 *     entries (R^{ττ}_1, R^{ττ}_τ).
 *
 *   σ_2 = F · σ_1 · F⁻¹ (change basis, braid, change back): involves
 *     the τττ F-matrix conjugation.
 *
 *   σ_3 = R again, diagonal in the basis with entries
 *     (R^{ττ}_1, R^{ττ}_τ).
 *
 * Repeated applications of σ_1, σ_2, σ_3, σ_1⁻¹, ... approximate any
 * U(2) gate to any desired precision — this is the universal-braiding
 * property of Fibonacci anyons.
 * ==================================================================== */

/** @brief Apply braid generator σ_1 (= σ_3 by symmetry) to the 2-dim
 *  fusion space of 4 Fibonacci anyons fusing to vacuum.
 *
 *  σ_1 is diagonal in the (1, τ) intermediate-channel basis:
 *    σ_1 |1⟩ = R^{ττ}_1 |1⟩ = exp(-4πi/5) |1⟩
 *    σ_1 |τ⟩ = R^{ττ}_τ |τ⟩ = exp(3πi/5)  |τ⟩
 *
 *  @param[in,out] state  Length-2 complex amplitude vector (basis
 *                        order: |1⟩, |τ⟩). Modified in place.
 *  @param[in]     inverse  If non-zero, apply σ_1⁻¹ instead. */
IRREP_API void
irrep_fib_braid_sigma_1(double _Complex state[2], int inverse);

/** @brief Apply braid generator σ_2 to the 2-dim fusion space.
 *
 *  σ_2 = F · diag(R^{ττ}_1, R^{ττ}_τ) · F⁻¹ — braid in the rotated
 *  fusion-tree basis. Performs a full 2×2 unitary on the state vector. */
IRREP_API void
irrep_fib_braid_sigma_2(double _Complex state[2], int inverse);

/** @brief Verify that σ_1 and σ_2 are both 2×2 unitaries.
 *
 *  Constructs the matrix forms by applying each generator to the two
 *  basis vectors and checks `U · U† = I` to machine precision.
 *
 *  @return Max |U·U† - I| across both generators. */
IRREP_API double
irrep_fib_braid_unitarity_residual(void);

/** @brief Verify the Yang-Baxter relation: σ_1·σ_2·σ_1 = σ_2·σ_1·σ_2.
 *
 *  This is the defining relation of the braid group B_3 and must hold
 *  for any valid braid representation.
 *
 *  @return Max |LHS - RHS| over the four 2×2 matrix entries. */
IRREP_API double
irrep_fib_braid_yang_baxter_residual(void);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_FIBONACCI_MTC_H */
