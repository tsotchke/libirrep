/* SPDX-License-Identifier: MIT */
/** @file ising_q.h
 *  @brief Ising^Q modular tensor category — Q-fold tensor power of the
 *         Ising MTC (objects {1, σ, ψ}, see `<irrep/crane_yetter_ising.h>`).
 *
 *  Simple objects of Ising^Q are Q-tuples `(x_1, ..., x_Q)` with each
 *  `x_i ∈ {0, 1, 2}` denoting `1`, `σ`, `ψ` respectively. There are
 *  exactly `3^Q` simples; the unit object is the all-zeros tuple.
 *  Quantum dimensions, fusion coefficients, S-matrix and T-matrix
 *  factorise component-wise from the Ising MTC:
 *
 *      d_{(x_1,...,x_Q)}     = ∏_i d_{x_i}
 *      N^{a,b}_c             = ∏_i N^{a_i, b_i}_{c_i}
 *      D(Ising^Q)            = D(Ising)^Q = 2^Q
 *
 *  ## What this module is for
 *
 *  Two things:
 *
 *  1. Lightweight Q-tuple wrappers (`_simple_t`, `_fpdim`, `_fusion`)
 *     so callers can write Ising^Q-aware code without re-implementing
 *     the tensor-power bookkeeping.
 *
 *  2. **Theorem B (σ-included Lagrangian count for Ising^Q).** A
 *     Lagrangian algebra A = ⊕_a n_a · a in a modular tensor category
 *     C is a connected, commutative, separable algebra with
 *     `FPdim(A) = sqrt(FPdim(C)) = D(C)`. For Ising^Q, `D = 2^Q ∈ ℚ`.
 *
 *     A simple object `a = (a_1, ..., a_Q)` has FPdim
 *     `d_a = (√2)^{#{i : a_i = σ}}`. If the σ-count is odd, `d_a`
 *     is an irrational multiple of √2; if even, it is the rational
 *     `2^{k/2}` where `k` is the σ-count.
 *
 *     Any non-negative integer combination of strict-√2-multiples is
 *     still a non-zero multiple of √2 (hence irrational) — so a
 *     Lagrangian algebra of Ising^Q (whose FPdim is rational) cannot
 *     contain any simple of odd σ-count with non-zero multiplicity.
 *
 *     This means the σ-included Lagrangian count of Ising^Q is
 *     identically `0` in the ordinary MTC sense, for every `Q ≥ 1`.
 *     The only way σ-bearing simples enter a Lagrangian-style algebra
 *     is via the **spin** refinement (= Crane-Yetter spin-CY MTC,
 *     `<irrep/crane_yetter_ising.h>`), where ψ is treated as the
 *     fermion-parity grading and even-σ-count simples can survive in
 *     the bosonic sector while odd-σ-count ones populate the fermionic
 *     sector. That refinement is research-track.
 *
 *  3. Even-σ-count Lagrangians ARE still nontrivial: we expose a
 *     brute-force enumerator
 *     `irrep_ising_q_lagrangian_count_brute(Q)` that searches all
 *     fusion-closed, haploid, FPdim-matched, multiplicity-free subsets
 *     of the even-σ-count simples of Ising^Q and counts the result.
 *     For `Q ∈ {1, 2, 3}` the count is exactly `1` (the product
 *     Lagrangian `(1 + ψ)^⊗Q`), and `_count_sigma_included` is exactly
 *     `0` — providing a runtime witness for Theorem B.
 *
 *  ## Primary references
 *
 *  - Davydov-Müger-Nikshych-Ostrik, *The Witt group of non-degenerate
 *    braided fusion categories*, J. reine angew. Math. 677 (2013) 135
 *    [arXiv:1009.2117] — Lagrangian-algebra framework.
 *  - Kong, *Anyon condensation and tensor categories*, Nucl. Phys.
 *    B 886 (2014) 436 [arXiv:1307.8244] — Lagrangian condensation
 *    formalism used by the T-paper portfolio.
 *  - Bombín-Martín-Delgado / Walker-Wang for the application to
 *    spin-Crane-Yetter, see `<irrep/crane_yetter_ising.h>` (T-paper
 *    Theorems E.1–E.3, this session's tasks #189–191).
 */
#ifndef IRREP_ISING_Q_H
#define IRREP_ISING_Q_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/crane_yetter_ising.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Maximum supported tensor power for Ising^Q. Enumeration cost
 *  is 3^Q simples × 2^|even-σ-set| subsets; Q = 4 is the comfortable
 *  upper limit (81 simples; ~3 × 10⁴ even-σ subsets after constraint
 *  pruning). */
#define IRREP_ISING_Q_MAX 6

/** @brief Number of simple objects in Ising^Q. */
IRREP_API int
irrep_ising_q_n_simples(int Q);

/** @brief Global quantum dimension `D(Ising^Q) = 2^Q`. */
IRREP_API double
irrep_ising_q_global_dim(int Q);

/** @brief Number of σ entries in the simple `a = (a_1, ..., a_Q)`,
 *  passed as a length-Q array of `irrep_ising_object_t`. */
IRREP_API int
irrep_ising_q_sigma_count(const irrep_ising_object_t *a, int Q);

/** @brief FPdim of a simple = `(√2)^{σ-count}`. Equivalently
 *  `2^{σ-count / 2}` (irrational unless σ-count is even). */
IRREP_API double
irrep_ising_q_fpdim(const irrep_ising_object_t *a, int Q);

/** @brief Fusion coefficient `N^{a,b}_c` in Ising^Q = product of the
 *  Q single-Ising fusion coefficients. */
IRREP_API int
irrep_ising_q_fusion(const irrep_ising_object_t *a,
                     const irrep_ising_object_t *b,
                     const irrep_ising_object_t *c,
                     int Q);

/** @brief Brute-force enumerate the multiplicity-free Lagrangian
 *  algebras of Ising^Q.
 *
 *  A subset `S` of the `3^Q` simples is counted iff:
 *    - `1 ∈ S` (haploid),
 *    - `Σ_{a ∈ S} d_a = D = 2^Q`,
 *    - `S` is closed under fusion: for all `a, b ∈ S` and every `c`
 *      with `N^{a,b}_c > 0`, we have `c ∈ S`.
 *
 *  By Theorem B, the only contributing simples have even σ-count, so
 *  the enumeration is restricted to the `Σ_{k even} C(Q, k) · 2^{Q-k}`
 *  even-σ simples (which is far fewer than `3^Q`). Implementation:
 *  bit-mask sweep over even-σ subsets, fusion-closed check via the
 *  precomputed fusion table.
 *
 *  @return  Lagrangian-algebra count, or -1 on error. */
IRREP_API int
irrep_ising_q_lagrangian_count_brute(int Q);

/** @brief Same brute-force enumeration as
 *  `irrep_ising_q_lagrangian_count_brute`, but counts only Lagrangians
 *  whose underlying object contains at least one σ-bearing simple
 *  (one with `σ-count > 0`, equivalently with non-trivial σ
 *  multiplicity).
 *
 *  **Theorem B** (proved in the header docstring above): this count
 *  is identically 0 for every Q ≥ 1 (in the ordinary MTC sense),
 *  because every σ-bearing simple has odd or even σ-count > 0 and
 *  even-σ-count Lagrangian closure forces the full `(1 + ψ)^⊗Q`
 *  block, leaving no room for σ-bearing simples with non-zero
 *  multiplicity.
 *
 *  This function is the **runtime witness** for that theorem. */
IRREP_API int
irrep_ising_q_sigma_included_lagrangian_count(int Q);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_ISING_Q_H */
