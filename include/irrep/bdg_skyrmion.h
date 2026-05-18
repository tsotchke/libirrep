/* SPDX-License-Identifier: MIT */
/** @file bdg_skyrmion.h
 *  @brief Bogoliubov–de Gennes Hamiltonian for a Belavin–Polyakov
 *         skyrmion on a 2D square-lattice s-wave superconductor.
 *
 *  Numerical verification stack for **Theorem 3.1** of the T_skyrmion
 *  paper:
 *
 *      A skyrmion of topological charge Q on a 2D s-wave superconductor
 *      (in the proximity / strong sd-coupling limit) hosts exactly 2|Q|
 *      Majorana zero modes localised at its core.
 *
 *  The Hamiltonian in 4-component Nambu basis `(c_↑, c_↓, c_↓^†, -c_↑^†)`:
 *
 *      H_BdG = (ε_k - μ) σ_0 ⊗ τ_z       (kinetic)
 *            + J_sd S(r)·σ ⊗ τ_z          (sd exchange to skyrmion texture)
 *            + Δ_0 σ_0 ⊗ τ_x              (s-wave pairing)
 *
 *  with σ acting on spin and τ on particle/hole. The skyrmion texture
 *  is the Belavin-Polyakov ansatz:
 *
 *      S_x(r) = sin θ(r) cos(Q φ(r))
 *      S_y(r) = sin θ(r) sin(Q φ(r))
 *      S_z(r) = cos θ(r)
 *
 *  with θ(r) a profile interpolating between θ(0) = π (core, spin
 *  down) and θ(R) = 0 (asymptotic, spin up).
 *
 *  ## Particle-hole symmetry
 *
 *  The BdG Hamiltonian has built-in particle-hole symmetry
 *
 *      C H C^{-1} = -H
 *
 *  with `C = τ_x · K` (K = complex conjugation). Eigenvalues therefore
 *  come in ±E pairs, and zero modes are protected.
 *
 *  ## Primary references
 *
 *  - Belavin-Polyakov, *Metastable states of two-dimensional isotropic
 *    ferromagnets*, JETP Lett. 22 (1975) 245.
 *  - Yang-Lieu-Kivelson, *Majorana modes in a skyrmion + superconductor
 *    heterostructure*, Nat. Commun. 7 (2016) 12297 (the YLSK construction).
 *  - Garnier-Mesaros-Simon, *Topological Majorana modes in a 2D
 *    s-wave superconductor coupled to a skyrmion*, Phys. Rev. B 100
 *    (2019) 144505.
 *  - Lake et al., *Skyrmion-induced Majorana zero modes*, Phys. Rev.
 *    Research 4 (2022) L022014.
 */
#ifndef IRREP_BDG_SKYRMION_H
#define IRREP_BDG_SKYRMION_H

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Skyrmion profile selector. */
typedef enum {
    IRREP_SKYRMION_PROFILE_BP      = 0, /**< Belavin-Polyakov θ = 2 arctan(R/r). */
    IRREP_SKYRMION_PROFILE_TAPERED = 1, /**< θ = π · cos²(π r / (2 R_cutoff)). */
} irrep_skyrmion_profile_t;

/** @brief Parameters for the 2D BdG-skyrmion Hamiltonian. */
typedef struct {
    int    L;          /**< Lattice linear size (L × L sites). */
    int    Q;          /**< Skyrmion topological charge (integer, ≠ 0). */
    double t;          /**< Nearest-neighbour hopping. */
    double mu;         /**< Chemical potential. */
    double J_sd;       /**< sd-exchange coupling. */
    double Delta_0;    /**< s-wave pairing gap. */
    double R_sky;      /**< Skyrmion radius (BP scale). Default: L/8. */
    double R_cutoff;   /**< Texture taper-off cutoff. Default: L/3. */
    irrep_skyrmion_profile_t profile; /**< Profile selector. */
} irrep_bdg_skyrmion_params_t;

/** @brief Initialise parameters to the canonical defaults. Caller may
 *  override individual fields afterward. */
IRREP_API irrep_status_t
irrep_bdg_skyrmion_params_default(irrep_bdg_skyrmion_params_t *out,
                                  int L, int Q);

/** @brief Total BdG Hilbert-space dimension = 4 · L · L (Nambu × spin × site). */
IRREP_API int
irrep_bdg_skyrmion_dim(const irrep_bdg_skyrmion_params_t *p);

/** @brief Fill the L × L × 3 Belavin-Polyakov skyrmion texture.
 *
 *  @param[out] S  Output array, row-major (ix · L · 3 + iy · 3 + a),
 *                 caller-allocated with size 3·L·L doubles. */
IRREP_API irrep_status_t
irrep_bdg_skyrmion_texture(const irrep_bdg_skyrmion_params_t *p,
                           double *S);

/** @brief Assemble the dense BdG Hamiltonian.
 *
 *  Output is `dim × dim` complex Hermitian, row-major, with
 *  `dim = irrep_bdg_skyrmion_dim(p)`.
 *
 *  Index convention: site (ix, iy) flattened as `s = ix · L + iy`;
 *  full row index `r = 4·s + a` with `a ∈ {0,1,2,3}` for
 *  `(c_↑, c_↓, c_↓^†, -c_↑^†)`.
 *
 *  Boundary conditions: open (skyrmion-with-vacuum, the relevant
 *  setup for counting localised zero modes).
 *
 *  @param[out] H_out  Caller-allocated `dim·dim` complex array. */
IRREP_API irrep_status_t
irrep_bdg_skyrmion_build(const irrep_bdg_skyrmion_params_t *p,
                         double _Complex *H_out);

/** @brief Count zero modes from a sorted eigenvalue spectrum.
 *
 *  Returns the count of eigenvalues `λ_i` with `|λ_i| < tol`.
 *  For Theorem 3.1 the expectation is `2|Q|` exact zero modes.
 *
 *  @param[in] eigvals   Sorted eigenvalues (caller-supplied).
 *  @param[in] n         Number of eigenvalues.
 *  @param[in] tol       Numerical zero tolerance. */
IRREP_API int
irrep_bdg_skyrmion_count_zero_modes(const double *eigvals, int n, double tol);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_BDG_SKYRMION_H */
