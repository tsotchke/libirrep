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

/* ----- Little-group machinery at Bloch momentum k ------------------- */

/** @brief Opaque handle for the little group at momentum `k = (kx/Lx)b1 +
 *         (ky/Ly)b2`: the subset of space-group elements `g` for which
 *         `R(g) · k ≡ k (mod reciprocal lattice)`, i.e. the point-group
 *         stabiliser of `k` semidirect-producted with the full translation
 *         subgroup.
 *
 *  On p6mm kagome at Γ the little point group is the full `C_6v` (order 12);
 *  at the three M-points it is `C_2v` (order 4); at the two K-points
 *  (only present on 3-divisible clusters like 6×6) it is `C_3v` (order 6).
 *
 *  Used to build `(k, μ_k)`-indexed space-group irreps: non-trivial
 *  little-group irreps at M and K expose the sector structure needed to
 *  distinguish gapped Z₂ and gapless Dirac spin-liquid candidates on the
 *  kagome Heisenberg model. */
typedef struct irrep_sg_little_group irrep_sg_little_group_t;

/** @brief Build the little group of @p G at Bloch momentum `(kx, ky)`.
 *
 *  @param G         parent space group (lattice must be attached; p4mm / p6mm
 *                   supported; p1 returns a handle with point_order = 1).
 *  @param kx, ky    Bloch indices; canonicalised mod `(Lx, Ly)` like
 *                   @ref irrep_sg_bloch_amplitude.
 *  @return          new handle, or `NULL` on OOM / bad input. Reason is
 *                   available via @ref irrep_last_error. */
IRREP_API irrep_sg_little_group_t *irrep_sg_little_group_build(const irrep_space_group_t *G, int kx,
                                                               int ky);

/** @brief Release a little-group handle. */
IRREP_API void irrep_sg_little_group_free(irrep_sg_little_group_t *lg);

/** @brief Number of point-group elements in the little point group. For Γ
 *         this equals `irrep_space_group_point_order(G)`; for M / K it is
 *         strictly smaller. */
IRREP_API int irrep_sg_little_group_point_order(const irrep_sg_little_group_t *lg);

/** @brief Total order of the little group, `point_order · (Lx · Ly)`. */
IRREP_API int irrep_sg_little_group_order(const irrep_sg_little_group_t *lg);

/** @brief Write the parent space-group point-op indices that form the little
 *         point group. Sorted ascending; length = @ref
 *         irrep_sg_little_group_point_order. */
IRREP_API void irrep_sg_little_group_point_ops(const irrep_sg_little_group_t *lg, int *out_indices);

/** @brief Recover the `(kx, ky)` this little group was built at. */
IRREP_API void irrep_sg_little_group_k(const irrep_sg_little_group_t *lg, int *out_kx,
                                       int *out_ky);

/** @brief Borrow the parent space group. */
IRREP_API const irrep_space_group_t *irrep_sg_little_group_parent(const irrep_sg_little_group_t *lg);

/** @brief Write the 2×2 integer rotation matrix of element @p i in
 *         lattice-cell basis. Column @p j is the image of `a_{j+1}` — for
 *         a proper rotation `det = +1`, for a mirror `det = −1`. Useful
 *         for assembling character rows on the caller side: on C_{n}v
 *         little groups the tuple `(order, det)` together with the
 *         reflection axis separates every conjugacy class. */
IRREP_API void irrep_sg_little_group_element_matrix(const irrep_sg_little_group_t *lg, int i,
                                                    int out_M[2][2]);

/** @brief Opaque handle for an irrep of the little POINT group (C_2v at M,
 *         C_3v at K, etc.). The Bloch phase `e^{-i k·t}` is carried
 *         separately by @ref irrep_sg_project_at_k; this object records
 *         only the point-group characters. */
typedef struct irrep_sg_little_group_irrep irrep_sg_little_group_irrep_t;

/** @brief Build a little-group-irrep handle from a character row.
 *
 *  @param lg          little-group parent
 *  @param characters  length-`point_order(lg)` array of characters,
 *                     one per element in the same order as
 *                     @ref irrep_sg_little_group_point_ops returns
 *  @param dim         dimension `d_μ` of the irrep (1 for C_2v, C_3v
 *                     one-dimensional irreps; 2 for C_3v's E irrep; etc.)
 *  @return            new handle or `NULL` on OOM / bad input. */
IRREP_API irrep_sg_little_group_irrep_t *
irrep_sg_little_group_irrep_new(const irrep_sg_little_group_t *lg, const double _Complex *characters,
                                int dim);

/** @brief Release an irrep handle. */
IRREP_API void irrep_sg_little_group_irrep_free(irrep_sg_little_group_irrep_t *mu_k);

/** @brief Dimension `d_μ` of the irrep. */
IRREP_API int irrep_sg_little_group_irrep_dim(const irrep_sg_little_group_irrep_t *mu_k);

/** @brief Named irreps of the little point groups that appear on p1 /
 *         p4mm / p6mm. Used by @ref irrep_sg_little_group_irrep_named to
 *         construct a handle without the caller assembling a character
 *         row by hand. */
typedef enum {
    IRREP_LG_IRREP_A1,  /**< Totally symmetric; valid for every little group. */
    IRREP_LG_IRREP_A2,  /**< Sign on reflections; C_2v / C_3v / C_4v / C_6v. */
    IRREP_LG_IRREP_B1,  /**< C_2v / C_4v / C_6v (not C_3v). */
    IRREP_LG_IRREP_B2,  /**< C_2v / C_4v / C_6v (not C_3v). */
    IRREP_LG_IRREP_E,   /**< 2D irrep on C_3v (at K on kagome). */
    IRREP_LG_IRREP_E1,  /**< First 2D irrep on C_6v (at Γ on kagome). */
    IRREP_LG_IRREP_E2,  /**< Second 2D irrep on C_6v. */
    IRREP_LG_IRREP_E_C4V /**< 2D irrep on C_4v (Γ / M on p4mm). */
} irrep_lg_named_irrep_t;

/** @brief Build a named irrep of the little point group at `lg`'s base
 *         point, auto-classifying each element by its conjugacy class
 *         within the abstract little group.
 *
 *  Supported so far: Γ on p6mm (C_6v) and K on p6mm (C_3v) — the two
 *  little groups the gapped-vs-gapless kagome-Heisenberg protocol
 *  consumes. Other little groups (C_2v at M on p6mm, C_4v at Γ/M on
 *  p4mm, C_s, C_1, …) return `NULL` with an informative error in
 *  @ref irrep_last_error pending follow-up support; the lower-level
 *  @ref irrep_sg_little_group_irrep_new remains available for every
 *  case.
 *
 *  @param lg    little group built via @ref irrep_sg_little_group_build.
 *  @param name  named irrep to construct.
 *  @return      handle or `NULL` on mismatch (e.g. `E_1` on a `C_3v`
 *               little group). Caller frees with
 *               @ref irrep_sg_little_group_irrep_free. */
IRREP_API irrep_sg_little_group_irrep_t *
irrep_sg_little_group_irrep_named(const irrep_sg_little_group_t *lg, irrep_lg_named_irrep_t name);

/** @brief Composite Bloch + little-group projection:
 *
 *  \f[ P_{k,\mu_k}\,\psi(\sigma) \;=\; \frac{d_{\mu_k}}{|G_k|}\;
 *      \sum_{t \in T}\;\sum_{p \in P_k}\; e^{-i k \cdot t}\;\chi_{\mu_k}^*(p)\;
 *      \psi\!\big((\tau_t p)\cdot\sigma\big) \f]
 *
 *  @param lg        little-group handle built at `(kx, ky)`
 *  @param mu_k      irrep of the little point group `P_k`
 *  @param psi_of_g  length-`order(G)` orbit amplitudes (same convention as
 *                   @ref irrep_sg_project_amplitude — element index
 *                   `g = tidx * point_order + p`)
 *  @return          projected amplitude, or `NaN + NaN·i` on NULL input. */
IRREP_API double _Complex irrep_sg_project_at_k(const irrep_sg_little_group_t *lg,
                                                const irrep_sg_little_group_irrep_t *mu_k,
                                                const double _Complex               *psi_of_g);

/** @brief Build an orthonormal basis for the `(k, μ_k)`-irrep sector of a
 *         `local_dim^num_sites`-dimensional configuration Hilbert space.
 *
 *  Applies @ref irrep_sg_project_at_k to each computational basis state
 *  and Gram-Schmidt-orthogonalises against previously accepted vectors.
 *  Non-zero residuals (norm > 1e-9) are normalised and appended.
 *
 *  Block-diagonal ED recipe: build the basis at each `(k, μ_k)`, sandwich
 *  the Hamiltonian against it, diagonalise each dense block. Sector
 *  dimensions sum to the full Hilbert-space dimension across all `(k, μ_k)`
 *  — a finer decomposition than the pure @ref irrep_sg_bloch_basis (which
 *  only resolves translations).
 *
 *  @param lg          little group at `(kx, ky)` (carries the parent
 *                     space group internally)
 *  @param mu_k        target irrep of the little point group
 *  @param num_sites   configuration length (must equal
 *                     `irrep_space_group_num_sites(parent)`)
 *  @param local_dim   per-site Hilbert-space dimension (2 for spin-½)
 *  @param basis_out   output buffer of length `n_max · local_dim^num_sites`
 *  @param n_max       row capacity of @p basis_out; pass the full
 *                     Hilbert-space dimension for safety
 *  @return            actual basis size, or `-1` on error (NULL input,
 *                     overflow, OOM). */
IRREP_API int irrep_sg_adapted_basis_at_k(const irrep_sg_little_group_t       *lg,
                                          const irrep_sg_little_group_irrep_t *mu_k, int num_sites,
                                          int local_dim, double _Complex *basis_out, int n_max);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_CONFIG_PROJECT_H */
