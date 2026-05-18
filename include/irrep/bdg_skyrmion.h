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

/* ====================================================================
 * Multi-skyrmion lattice texture.
 *
 * A spatial arrangement of multiple skyrmions on the same `L × L` SC
 * lattice. Each skyrmion has its own center, charge, and radius; the
 * total magnetisation at site `r` is the sum of individual Belavin-
 * Polyakov spin vectors, normalised back to unit length:
 *
 *      m_tot(r)  = ∑_i m_BP^{(Q_i, R_i)}(r - r_i)
 *      m̂_tot(r) = m_tot(r) / |m_tot(r)|.
 *
 * The composite BdG Hamiltonian then enters via the same s-wave +
 * J_sd coupling on m̂_tot. By topological additivity (each isolated
 * skyrmion contributing 2|Q_i| Majorana zero modes per Theorem 3.1),
 * a lattice of N skyrmions hosts up to `2 · ∑_i |Q_i|` Majorana zero
 * modes when the centers are sufficiently separated for the inter-
 * skyrmion overlap to vanish.
 *
 * The simple "sum then normalise" rule does NOT preserve the
 * topological charge of m̂_tot exactly — for closely-spaced skyrmions
 * the integral charge can dip below ∑ |Q_i| due to interference of
 * the BP tail fields. In practice, set inter-center distance ≥ 4 R_sky
 * to keep the MZM count = 2 ∑ |Q_i| up to finite-size effects.
 * ==================================================================== */

/** @brief One skyrmion center in a lattice. */
typedef struct {
    int    x0, y0;   /**< Integer-lattice center position. */
    int    Q;        /**< Topological charge (≠ 0). */
    double R_sky;    /**< BP radius for this skyrmion. */
    double R_cutoff; /**< Taper-off cutoff. */
    irrep_skyrmion_profile_t profile;
} irrep_skyrmion_center_t;

/** @brief Multi-skyrmion lattice parameters. */
typedef struct {
    int    L;                                 /**< Linear lattice size. */
    int    n_skyrmions;                       /**< Number of centers. */
    const irrep_skyrmion_center_t *centers;   /**< Length `n_skyrmions`. */
    double t;                                 /**< NN hopping. */
    double mu;                                /**< Chemical potential. */
    double J_sd;                              /**< sd-coupling. */
    double Delta_0;                           /**< Pairing gap. */
} irrep_bdg_skyrmion_lattice_t;

/** @brief Fill the `L × L × 3` magnetisation texture from the
 *  superposition of all skyrmion centers. */
IRREP_API irrep_status_t
irrep_bdg_skyrmion_lattice_texture(const irrep_bdg_skyrmion_lattice_t *p,
                                   double *S);

/** @brief Build the BdG Hamiltonian on the multi-skyrmion texture.
 *  Dim = `4·L²`; same Nambu basis and BC convention as
 *  `irrep_bdg_skyrmion_build`. */
IRREP_API irrep_status_t
irrep_bdg_skyrmion_lattice_build(const irrep_bdg_skyrmion_lattice_t *p,
                                 double _Complex *H_out);

/* ====================================================================
 * Sparse (matrix-free) BdG diagonalisation via Lanczos.
 *
 * The BdG Hamiltonian on the L × L lattice has ≈ 16 non-zero entries
 * per row (4 on-site Nambu coupling + 4·2 hopping with both particle
 * and hole blocks) out of `dim = 4 L²`, so storage scales as `O(L²)`
 * rather than `O(L⁴)`. Dense Jacobi at L = 20 needs ~10 GB of memory
 * and several minutes per diagonalisation; the sparse-apply + Lanczos
 * path delivers the K-lowest-|E| eigenvalues in seconds at L = 32+.
 *
 * MZMs are interior eigenvalues near E = 0 — not extremal — so we
 * Lanczos on `-H²` rather than `H` directly. The eigenvalues of `-H²`
 * are `-|E_i|²`, and Lanczos converges on the algebraically-largest
 * (= |E|-smallest) end of the spectrum efficiently. After
 * convergence we take `|E_i| = √(-λ_i)` and report the K smallest.
 * ==================================================================== */

/** @brief Compute the K smallest distinct |E| BdG-spectrum eigenvalues
 *  by matrix-free Lanczos on `H²`.
 *
 *  Internally:
 *    - Precomputes the magnetisation texture once.
 *    - Runs `irrep_lanczos_eigvals_reorth` with a sparse-apply
 *      callback that applies `H² · x = H · (H · x)` site by site.
 *    - Returns the K algebraically-smallest eigenvalues of `H²`
 *      (= smallest `|E|²`); reports `|E_i| = sqrt(λ_i)` in ascending
 *      order.
 *
 *  ## PH-pair collapse semantic
 *
 *  BdG eigenvalues come in PH-mandated ±E pairs, so each distinct
 *  `|E|² > 0` has multiplicity ≥ 2 in the spectrum of `H²`. Lanczos
 *  on `H²` typically converges to ONE Ritz value per `|E|²` cluster at
 *  moderate iteration counts (e.g. `max_iters ≤ 0.5 · dim`); only at
 *  near-full Krylov depth does it resolve the multiplicity. The
 *  returned array therefore contains DISTINCT `|E|` values — one
 *  representative per ±E pair. Multiply each returned value's BdG
 *  multiplicity by 2 when comparing to a dense spectrum.
 *
 *  ## Iteration budget
 *
 *  Recommended `max_iters` ≈ `8 · k_wanted` for well-separated low
 *  modes; the call typically converges to 1e-10 in 50-150 iterations
 *  for `k_wanted` ≤ 8 and L ≤ 32.
 *
 *  @param[in]  p              BdG-skyrmion parameters.
 *  @param[in]  k_wanted       Number of distinct lowest |E| values.
 *  @param[in]  max_iters      Maximum Lanczos iterations (≥ k_wanted).
 *  @param[out] abs_eigvals_out  Caller-allocated length-k_wanted, sorted
 *                              ascending. */
IRREP_API irrep_status_t
irrep_bdg_skyrmion_lanczos_lowest_abs_eigvals(
    const irrep_bdg_skyrmion_params_t *p,
    int k_wanted,
    int max_iters,
    double *abs_eigvals_out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_BDG_SKYRMION_H */
