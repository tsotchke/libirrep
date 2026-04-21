/* SPDX-License-Identifier: MIT */
/** @file sym_group.h
 *  @brief Symmetric group `S_N`, Young tableaux, hook-length dimension formula,
 *         and antisymmetric / totally-symmetric projectors on tensor-factored
 *         wavefunctions.
 *
 *  Required for the 2D-Hubbard / fermion-sign-problem target: the
 *  antisymmetric projector implements the Slater-determinant structure of a
 *  fermion wavefunction, and the hook-length formula gives the dimension of
 *  each `S_N` irrep as scheduled in the multi-flavour extension.
 *
 *  ### Scope
 *  `S_N` for `N ≤ 10` is fully tabulated; full-permutation enumeration
 *  costs `N! · N` words (3.6 M entries at N = 10). The hook-length formula
 *  is exact for any `N`. The antisymmetric / symmetric projectors operate
 *  on a tensor state `ψ[i_0, i_1, …, i_{N-1}]` of total size
 *  `local_dim^N`; practical limits are `N ≤ 10, local_dim ≤ 8`.
 *
 *  Higher-N use cases (sparse Slater-determinant bases in Hubbard) should
 *  compose this module's primitives with external sparse storage; libirrep
 *  exposes the group-theoretic machinery rather than the representation
 *  of the quantum state.
 *
 *  Cited: James & Kerber, *The Representation Theory of the Symmetric
 *  Group* (Addison-Wesley 1981); Fulton, *Young Tableaux* (Cambridge 1997).
 */
#ifndef IRREP_SYM_GROUP_H
#define IRREP_SYM_GROUP_H

#include <stddef.h>
#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief `n!` for `n ∈ [0, 20]`. Returns `-1` for out-of-range input to
 *         signal overflow; the largest representable factorial in int64
 *         is 20! = 2 432 902 008 176 640 000. */
IRREP_API long long irrep_factorial(int n);

/** @brief Sign of a permutation: `(-1)` raised to the number of inversions.
 *         Returns `+1` or `−1`; invalid input (duplicates, out-of-range
 *         entries) returns `0`. */
IRREP_API int irrep_permutation_sign(const int *perm, int n);

/** @brief Write every permutation of `[0, n)` into a buffer of shape
 *         `[n factorial][n]`, in lexicographic order (Heap's algorithm
 *         gives a different order; we use a canonical lex sweep so the
 *         identity permutation is at index 0). Requires `n ≤ 10`. Returns
 *         #IRREP_OK or an error. */
IRREP_API irrep_status_t irrep_permutations_all(int n, int *out_perms);

/** @brief Hook-length dimension formula for the irrep of `S_N` labelled
 *         by integer partition `λ`: `dim(λ) = (N factorial) divided by the
 *         product of hook lengths over cells of the Young diagram`.
 *
 *  @param partition  weakly-decreasing positive parts
 *  @param n_parts    length of @p partition
 *  @return           irrep dimension, or `-1` on invalid partition
 *                    (non-decreasing, non-positive, or arithmetic overflow
 *                    for `N > 20`). */
IRREP_API long long irrep_young_dim(const int *partition, int n_parts);

/** @brief Apply the antisymmetric (fermion, sign-representation) projector
 *
 *  \f[ \mathcal A \psi(i_0,\dots,i_{N-1}) \;=\; \frac{1}{N!}
 *      \sum_\sigma \mathrm{sign}(\sigma)\,\psi(i_{\sigma(0)},\dots,i_{\sigma(N-1)}) \f]
 *
 *  to an array-encoded tensor state `ψ[i]` with `i = Σ_k i_k · d^k` and
 *  `d = local_dim`. Total state dimension `local_dim^N`.
 *
 *  The output is zero on the "Pauli-forbidden" basis states where any two
 *  indices coincide. For `N ≤ 10`, the projector is applied by full
 *  permutation enumeration.
 *
 *  @param N          number of tensor factors (1 ≤ N ≤ 10)
 *  @param local_dim  per-factor dimension (≥ 2)
 *  @param psi_in     input state, length `local_dim^N`
 *  @param psi_out    output state, length `local_dim^N` (may alias)
 *  @return           #IRREP_OK, #IRREP_ERR_INVALID_ARG, or #IRREP_ERR_OUT_OF_MEMORY. */
IRREP_API irrep_status_t irrep_sym_group_antisymmetrize(int N, int local_dim,
                                                        const double _Complex *psi_in,
                                                        double _Complex       *psi_out);

/** @brief Totally-symmetric (boson, trivial-representation) projector.
 *         Same conventions as @ref irrep_sym_group_antisymmetrize but
 *         sums over permutations without sign. */
IRREP_API irrep_status_t irrep_sym_group_symmetrize(int N, int local_dim,
                                                    const double _Complex *psi_in,
                                                    double _Complex       *psi_out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SYM_GROUP_H */
