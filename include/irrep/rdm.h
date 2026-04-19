/* SPDX-License-Identifier: MIT */
/** @file rdm.h
 *  @brief Reduced density matrix, entanglement-entropy, and Kitaev-Preskill
 *         topological-entanglement-entropy primitives.
 *
 *
 *  ### Scope
 *
 *  Primitives for the Kagome-Heisenberg ground-state-nature programme:
 *
 *    - @ref irrep_partial_trace — compute `ρ_A = Tr_B |ψ⟩⟨ψ|` from a full
 *      state vector on `num_sites` sites of local dimension `local_dim`.
 *      Generic in `local_dim`; limited to `num_sites ≤ 30` by the
 *      `int64_t` basis-state index.
 *    - @ref irrep_hermitian_eigvals — diagonalise a Hermitian matrix via
 *      cyclic-Jacobi sweeps; returns eigenvalues only.
 *    - @ref irrep_entropy_vonneumann_spectrum and
 *      @ref irrep_entropy_renyi_spectrum — take an eigenvalue list and
 *      return `S_VN = − Σ λ ln λ` or `S_α = (1/(1−α)) ln Σ λ^α`.
 *    - @ref irrep_entropy_vonneumann / @ref irrep_entropy_renyi — fused
 *      convenience: diagonalise `ρ_A`, then return the entropy.
 *    - @ref irrep_topological_entanglement_entropy — assemble `γ` from the
 *      seven entropies of the Kitaev-Preskill tripartition:
 *      `γ = S_A + S_B + S_C − S_AB − S_BC − S_AC + S_ABC`.
 *
 *  Entropies use natural log (nats). Multiply by `1/ln 2` for bits.
 *
 *  ### Kagome diagnosis
 *
 *  For a gapped Z₂ spin liquid on kagome, `γ = ln 2 ≈ 0.693 15`. For a
 *  trivially ordered state, `γ = 0`. For a gapless Dirac liquid, `γ`
 *  flows logarithmically with region size and is not well defined as a
 *  constant. Measuring `γ` on sufficiently large regions and extrapolating
 *  is the primary observable.
 *
 *  ### Scaling
 *
 *  The full-state API is only tractable at exact-diagonalisation scale
 *  (say `num_sites ≤ 24` on spin-½). For the 108-site 6×6 target,
 *  the partial trace is assembled from MCMC or MPS contractions outside
 *  libirrep; downstream callers pass the resulting `ρ_A` directly to
 *  @ref irrep_entropy_vonneumann.
 */
#ifndef IRREP_RDM_H
#define IRREP_RDM_H

#include <stddef.h>
#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Compute the reduced density matrix `ρ_A = Tr_B |ψ⟩⟨ψ|` for a pure
 *         state `|ψ⟩` on `num_sites` sites of local dimension `local_dim`.
 *
 *  Basis convention: the composite index `i ∈ [0, local_dim^{num_sites})`
 *  encodes site digits in little-endian order,
 *  `i = Σ_s digit_s · local_dim^s`.
 *
 *  @param num_sites    total site count (1 ≤ num_sites ≤ 30)
 *  @param local_dim    local Hilbert-space dimension per site (≥ 2)
 *  @param psi          length `local_dim^num_sites` amplitudes, unit norm
 *                      preferred; normalisation is *not* enforced
 *  @param sites_A      sorted ascending list of A-subsystem site indices
 *  @param nA           length of @p sites_A (0 ≤ nA ≤ num_sites)
 *  @param rho_A        output `dA × dA` matrix, row-major,
 *                      `dA = local_dim^nA`
 *  @return             #IRREP_OK, or #IRREP_ERR_INVALID_ARG. */
IRREP_API irrep_status_t
irrep_partial_trace(int num_sites, int local_dim,
                    const double _Complex *psi,
                    const int *sites_A, int nA,
                    double _Complex *rho_A);

/** @brief Cyclic-Jacobi Hermitian eigendecomposition. Destroys @p A (writes
 *         a diagonal matrix on exit). @p eigvals must have room for `n`
 *         entries. Converged when the off-diagonal Frobenius norm drops
 *         below `n · 1e-14 · tr(A)` or after 50 · n sweeps, whichever
 *         comes first. Returns #IRREP_OK or #IRREP_ERR_INVALID_ARG. */
IRREP_API irrep_status_t
irrep_hermitian_eigvals(int n, double _Complex *A, double *eigvals);

/** @brief von Neumann entropy `− Σ λ_i ln λ_i` from an eigenvalue list.
 *         Treats `λ ≤ 1e-15` as zero (since `0 ln 0 = 0`). */
IRREP_API double
irrep_entropy_vonneumann_spectrum(const double *eigvals, int n);

/** @brief Rényi entropy of order `α ≠ 1`,
 *         `S_α = (1/(1−α)) ln Σ λ_i^α`. For `α = 1` falls back to
 *         @ref irrep_entropy_vonneumann_spectrum. */
IRREP_API double
irrep_entropy_renyi_spectrum(const double *eigvals, int n, double alpha);

/** @brief Fused: diagonalise @p rho (copy internally so @p rho survives),
 *         then take von Neumann entropy. */
IRREP_API double
irrep_entropy_vonneumann(const double _Complex *rho, int n);

/** @brief Fused Rényi counterpart to @ref irrep_entropy_vonneumann. */
IRREP_API double
irrep_entropy_renyi(const double _Complex *rho, int n, double alpha);

/** @brief Sparse Hermitian eigenvalue extraction via 3-term-recurrence
 *         Lanczos. Ground-state energy and a few low-lying excited states
 *         from a matrix-free Hermitian operator, without materialising
 *         the full matrix.
 *
 *  Algorithm: Lanczos tridiagonalisation with no reorthogonalisation (3
 *  state vectors held at a time, α / β scalars stored). The resulting
 *  tridiagonal is diagonalised via @ref irrep_hermitian_eigvals. For large
 *  `max_iters` ghost eigenvalues can appear as duplicates of converged
 *  Ritz values — filter by multiplicity if needed.
 *
 *  Memory: `3 · dim` complex doubles; no O(k · dim) Krylov basis.
 *
 *  ### Parameters
 *  @param apply_op    Callback that computes `y = H · x` on vectors of
 *                     length @p dim. Caller threads state through @p ctx.
 *  @param ctx         Opaque context passed to @p apply_op.
 *  @param dim         Hilbert-space dimension (vector length).
 *  @param k_wanted    Number of lowest eigenvalues to return.
 *  @param max_iters   Lanczos iteration count (typical 50–200).
 *                     For accurate extremal eigenvalues `max_iters` should
 *                     exceed `2·k_wanted` with margin; larger values refine
 *                     at the cost of potential ghost eigenvalues (always
 *                     appearing as near-duplicate extremal values).
 *  @param seed        Length-`dim` initial vector; if `NULL`, the builder
 *                     uses a deterministic pseudo-random seed. The seed
 *                     should have non-zero overlap with the eigenvectors
 *                     of interest — pathological seeds (e.g. uniformly
 *                     spread in a momentum sector that doesn't contain
 *                     the ground state) silently miss the ground state.
 *  @param eigvals_out Length-`k_wanted` output, sorted ascending. */
IRREP_API irrep_status_t
irrep_lanczos_eigvals(
    void (*apply_op)(const double _Complex *x,
                     double _Complex *y,
                     void *ctx),
    void *ctx,
    long long dim,
    int k_wanted,
    int max_iters,
    const double _Complex *seed,
    double *eigvals_out);

/** @brief Kitaev-Preskill topological entanglement entropy `γ`:
 *
 *      γ = S_A + S_B + S_C − S_AB − S_BC − S_AC + S_ABC
 *
 *  The caller supplies the seven entropies of the tripartition. A value of
 *  `γ = ln 2` signals Z₂ topological order; `γ = 0` signals no topological
 *  order. */
IRREP_API double
irrep_topological_entanglement_entropy(double SA, double SB, double SC,
                                       double SAB, double SBC, double SAC,
                                       double SABC);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_RDM_H */
