/* SPDX-License-Identifier: MIT */
/** @file config_project.h
 *  @brief Configuration-space projection onto a space-group irrep.
 *
 *  The projection operator on functions of a configuration `σ` (a spin
 *  assignment, an occupation pattern, a permutation; whatever the caller's
 *  Hilbert-space basis looks like) is
 *
 *  \f[ \mathcal P_\mu \psi(\sigma) \;=\; \frac{d_\mu}{|G|}\;\sum_{g\in G}\;
 *      \chi_\mu^*(g) \; \psi(g\cdot\sigma) \f]
 *
 *  where `G` is the space group, `μ` is a target irrep of `G` with dimension
 *  `d_μ` and character `χ_μ(g)`, and `g·σ` is the action of `g` on `σ` (a
 *  site permutation — see @ref irrep/space_group.h).
 *
 *  Downstream usage pattern (symmetric neural-quantum-state ansatz):
 *
 *    1. caller samples σ via MCMC;
 *    2. caller enumerates the orbit `{g·σ : g ∈ G}` with
 *       @ref irrep_sg_enumerate_orbit — `order · num_sites` doubles;
 *    3. caller evaluates the wavefunction amplitude `ψ(g·σ)` at each orbit
 *       member — this is the expensive step (one NN forward per element);
 *    4. caller reduces with @ref irrep_sg_project_amplitude or the
 *       convenience shortcut @ref irrep_sg_project_A1 (totally-symmetric).
 *
 *  The library supplies the book-keeping and the character-weighted reduce
 *  but not the wavefunction — that is specific to the ansatz family and
 *  lives downstream. For the 6×6 kagome A₁ projection the orbit has 432
 *  images, and the reduce is a single 432-term dot product.
 *
 *  Distinct from @ref irrep/point_group.h which exposes the analogous
 *  projection on *irrep-space feature vectors* (acting on the internal
 *  hidden state of an equivariant layer). This one acts on *configurations*
 *  (the spin degrees of freedom on the lattice), which is what NQS
 *  symmetrisation requires.
 */
#ifndef IRREP_CONFIG_PROJECT_H
#define IRREP_CONFIG_PROJECT_H

#include <stddef.h>
#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>
#include <irrep/space_group.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Row of a space-group character table for one target irrep. */
typedef struct irrep_sg_irrep irrep_sg_irrep_t;

/** @brief Build a handle for an arbitrary irrep given its character row.
 *
 *  @param G           space group whose elements the characters index
 *  @param characters  array of length `irrep_space_group_order(G)` with
 *                     `χ_μ(g)` at index `g`
 *  @param irrep_dim   dimension `d_μ` of the target irrep
 *  @return            new handle, or `NULL` on OOM / bad input. */
IRREP_API irrep_sg_irrep_t *irrep_sg_irrep_new(const irrep_space_group_t *G,
                                               const double _Complex *characters, int irrep_dim);

/** @brief Release a handle returned by @ref irrep_sg_irrep_new or
 *         @ref irrep_sg_trivial. */
IRREP_API void irrep_sg_irrep_free(irrep_sg_irrep_t *mu);

/** @brief Convenience builder for the totally-symmetric (trivial / A₁) irrep:
 *         `χ(g) = 1` for every `g`, `d_μ = 1`. */
IRREP_API irrep_sg_irrep_t *irrep_sg_trivial(const irrep_space_group_t *G);

/** @brief Convenience builder for the sign-representation (A₂ for p4mm and
 *         p6mm; reduces to the trivial irrep for p1). Characters are `+1` on
 *         all translations and on the pure rotations, `−1` on reflections.
 *
 *  This irrep is the natural companion to @ref irrep_sg_trivial for NQS
 *  ansätze whose wavefunction changes sign under a mirror — e.g. chiral
 *  spin-liquid candidates on the kagome lattice. */
IRREP_API irrep_sg_irrep_t *irrep_sg_sign_rep(const irrep_space_group_t *G);

/** @brief Reduce a list of amplitudes `ψ_g = ψ(g·σ)` into the projected
 *         amplitude `P_μ ψ(σ)`.
 *
 *  @param mu         irrep handle
 *  @param psi_of_g   length-`order(G)` amplitudes, one per group element
 *  @return           projected amplitude on success; `NaN + NaN·i` on
 *                    `NULL` input (zero is a legitimate result, so we use
 *                    IEEE-754 NaN to distinguish error from answer; caller
 *                    may check with `isnan(creal(·))`). */
IRREP_API double _Complex irrep_sg_project_amplitude(const irrep_sg_irrep_t *mu,
                                                     const double _Complex  *psi_of_g);

/** @brief Shortcut: totally-symmetric projection,
 *         `(1/|G|) Σ_g ψ(g·σ)`.
 *         Returns `NaN + NaN·i` on `NULL` input. */
IRREP_API double _Complex irrep_sg_project_A1(const irrep_space_group_t *G,
                                              const double _Complex     *psi_of_g);

/** @brief Enumerate the orbit of a real-valued configuration `σ` under `G`.
 *
 *  Writes `out_orbit[g · num_sites + s] = σ[g⁻¹ · s]` for every group
 *  element `g`: the pullback of `σ` by `g`, i.e. `(g·σ)(s) = σ(g⁻¹·s)`.
 *
 *  The output buffer must have capacity `order(G) · num_sites`. */
IRREP_API void irrep_sg_enumerate_orbit(const irrep_space_group_t *G, const double *sigma,
                                        double *out_orbit);

/** @brief Build an orthonormal basis for the `μ`-irrep sector of a
 *         `local_dim^num_sites`-dimensional Hilbert space on which `G` acts
 *         by site permutation.
 *
 *  For each computational basis state `|s⟩` (`s ∈ [0, local_dim^num_sites)`),
 *  the function computes `|v⟩ = P_μ|s⟩` via the character-weighted sum
 *  and Gram-Schmidt-orthogonalises against previously-accepted basis
 *  vectors. Non-zero residuals (norm > 1e-9) are normalised and appended.
 *
 *  The resulting basis spans the `μ`-isotypic component of the Hilbert
 *  space. Its dimension is `m_μ · d_μ` where `m_μ` is the multiplicity of
 *  `μ` in the permutation representation.
 *
 *  This is the core primitive for symmetry-adapted block exact
 *  diagonalisation: once the basis is built, a Hamiltonian commuting with
 *  `G` has matrix elements only within each block, so diagonalising each
 *  block separately costs O(b³) for block size `b` instead of `O(D³)` for
 *  total Hilbert-space dimension `D`. On 18- and 24-site kagome clusters
 *  where naive ED is intractable, the per-sector blocks stay small enough
 *  for `irrep_hermitian_eigvals` to handle.
 *
 *  ### Parameters
 *  @param G           Space group whose elements permute sites.
 *  @param mu          Target irrep handle (from @ref irrep_sg_irrep_new,
 *                     @ref irrep_sg_trivial, or @ref irrep_sg_sign_rep).
 *  @param num_sites   Number of sites in the configuration space.
 *                     Must equal `irrep_space_group_num_sites(G)`.
 *  @param local_dim   Per-site Hilbert-space dimension (typically 2 for
 *                     spin-½).
 *  @param basis_out   Output buffer of length `n_max · local_dim^num_sites`.
 *                     Each row is one basis vector, written contiguously
 *                     row-major.
 *  @param n_max       Capacity of @p basis_out in rows. An upper bound on
 *                     the sector dimension is `dim V · d_μ / |G_point|`
 *                     for 1D irreps and `2 · dim V · d_μ / |G_point|` for
 *                     2D irreps; passing `dim V` is always safe.
 *
 *  @return            Actual basis size (≤ `n_max`); `-1` on error
 *                     (NULL pointer, overflow, OOM). */
IRREP_API int irrep_sg_adapted_basis(const irrep_space_group_t *G, const irrep_sg_irrep_t *mu,
                                     int num_sites, int local_dim, double _Complex *basis_out,
                                     int n_max);

/** @brief Project an amplitude onto Bloch momentum `k = (kx/Lx) b1 + (ky/Ly) b2`.
 *
 *  Sums the translation-subgroup amplitudes with the Bloch phase:
 *
 *  \f[ P_k\,\psi(\sigma) \;=\; \frac{1}{N_T}\sum_{t\in T} e^{-i k\cdot t}\,\psi(T_t\cdot\sigma) \f]
 *
 *  Only the identity point element at each translation is read (`g = tidx ·
 *  point_order + 0` in the group element layout). This is the *abelian*
 *  momentum projector — it ignores any extra point-group little-group
 *  structure at high-symmetry k, so at Γ/M/K the result still needs a
 *  subsequent point-group projection for full isotypic decomposition.
 *
 *  @param G         space group (translation subgroup used; point subgroup ignored)
 *  @param kx, ky    Bloch indices; any integer is accepted and canonicalised
 *                   to `[0, Lx) × [0, Ly)` modulo the cluster
 *  @param psi_of_g  length-`irrep_space_group_order(G)` amplitudes as produced
 *                   by @ref irrep_sg_enumerate_orbit followed by a caller
 *                   wavefunction evaluation
 *  @return          projected amplitude on success; `NaN + NaN·i` on NULL
 *                   input or if the space group has no lattice handle. */
IRREP_API double _Complex irrep_sg_bloch_amplitude(const irrep_space_group_t *G, int kx, int ky,
                                                   const double _Complex *psi_of_g);

/** @brief Build an orthonormal basis for the Bloch-momentum-`k` sector of
 *         `local_dim^num_sites`-dimensional Hilbert space.
 *
 *  For each computational basis state `|s⟩`, computes
 *  `|v⟩ = (1/N_T) Σ_t e^{-i k·t} T_t|s⟩` and Gram-Schmidt-orthogonalises
 *  against previously accepted vectors. Non-zero residuals are normalised
 *  and appended.
 *
 *  The sector dimension sums to the full Hilbert-space dimension across the
 *  `N_T = Lx · Ly` momentum sectors, so the average sector is `D/N_T`. On
 *  6×6 kagome (108 sites, `D = 2^108`) this is impractical; on 18-site
 *  kagome in the `S_z = 0` subspace the per-k block is ~900 states.
 *
 *  @param G          space group (translation subgroup used)
 *  @param kx, ky     integer Bloch indices; canonicalised mod `(Lx, Ly)`
 *  @param num_sites  number of sites (must equal `irrep_space_group_num_sites(G)`)
 *  @param local_dim  per-site Hilbert-space dimension (2 for spin-½)
 *  @param basis_out  output buffer of length `n_max · local_dim^num_sites`
 *  @param n_max      row capacity of @p basis_out; pass the full Hilbert-space
 *                    dimension for safety
 *  @return           actual basis size, or `-1` on error */
IRREP_API int irrep_sg_bloch_basis(const irrep_space_group_t *G, int kx, int ky, int num_sites,
                                   int local_dim, double _Complex *basis_out, int n_max);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_CONFIG_PROJECT_H */
