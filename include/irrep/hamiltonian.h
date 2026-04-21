/* SPDX-License-Identifier: MIT */
/** @file hamiltonian.h
 *  @brief On-the-fly Hamiltonian apply-operators for ED / Lanczos.
 *
 *  Every example under `examples/` that does exact diagonalisation
 *  previously re-implemented the same Heisenberg `apply_op` callback
 *  by hand. This header promotes the canonical forms to library
 *  primitives, so downstream consumers plug them directly into
 *  `irrep_lanczos_eigvals` without writing their own bit-twiddling.
 *
 *  Convention for spin-½ computational basis: integer `s ∈ [0, 2^N)`
 *  represents the `N`-bit Ising configuration with bit `i = (s >> i) & 1`
 *  = {1 → up, 0 → down}. `S_z_total = (popcount(s) − N/2)/2` is exactly
 *  conserved by every model here; callers who want fixed-magnetisation
 *  ED should seed the Lanczos starting vector on a single `S_z` sector
 *  (the apply functions preserve it automatically).
 *
 *  All apply functions have the signature required by the sparse
 *  eigensolver `irrep_lanczos_eigvals`: `(const double _Complex *psi,
 *  double _Complex *out, void *opaque)`. Pass the returned opaque
 *  handle as the `ctx` argument to Lanczos.
 */
#ifndef IRREP_HAMILTONIAN_H
#define IRREP_HAMILTONIAN_H

#include <complex.h>
#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Opaque handle for a spin-½ Heisenberg Hamiltonian
 *  `H = J · Σ_{⟨i,j⟩ ∈ bonds} S_i · S_j`.
 *
 *  Built from an explicit NN-bond list (as produced by
 *  @ref irrep_lattice_fill_bonds_nn) plus a coupling constant. */
typedef struct irrep_heisenberg irrep_heisenberg_t;

/** @brief Construct a Heisenberg operator on `num_sites` spin-½ sites.
 *
 *  The bond list `{(bi[b], bj[b]) : b ∈ [0, num_bonds)}` is copied into
 *  the handle — the caller is free to free its arrays after this call.
 *
 *  @param num_sites  Number of sites (≤ 30 for typical 2^N Hilbert spaces).
 *  @param num_bonds  Length of the bond arrays.
 *  @param bi, bj     NN bond endpoints; `bi[b] < bj[b]` is not required
 *                    but duplicates are not deduplicated.
 *  @param J          Coupling constant (positive = antiferromagnetic).
 *  @return Owned handle, or NULL on OOM / invalid input. */
IRREP_API irrep_heisenberg_t *
irrep_heisenberg_new(int num_sites, int num_bonds,
                     const int *bi, const int *bj,
                     double J);

/** @brief Release a handle returned by @ref irrep_heisenberg_new. */
IRREP_API void irrep_heisenberg_free(irrep_heisenberg_t *H);

/** @brief Apply callback with the signature required by
 *         @ref irrep_lanczos_eigvals.
 *
 *  Computes `out = H · psi` on length-`2^num_sites` complex vectors.
 *  For each bond `(i, j)`:
 *    - `S_i^z · S_j^z |s⟩ = ±¼ |s⟩`   (sign = +1 if aligned, −1 if anti)
 *    - `½(S_i^+ S_j^- + S_i^- S_j^+) |s⟩ = ½ |s ⊕ 2^i ⊕ 2^j⟩`
 *      when the two spins are anti-aligned, zero otherwise.
 *
 *  Time complexity: `O(num_bonds · 2^num_sites)` per call. Memory:
 *  writes into `out`; reads `psi`. Thread-safe on disjoint `(psi, out)`
 *  pairs with a shared `opaque`. */
IRREP_API void
irrep_heisenberg_apply(const double _Complex *psi,
                       double _Complex       *out,
                       void                  *opaque);

/** @brief Number of sites the handle was built with. */
IRREP_API int irrep_heisenberg_num_sites(const irrep_heisenberg_t *H);

/** @brief Hilbert-space dimension (`2^num_sites`). Convenience for
 *         sizing Lanczos buffers. */
IRREP_API long long irrep_heisenberg_dim(const irrep_heisenberg_t *H);

/** @brief Build a J₁-J₂ Heisenberg on two bond sets (NN + NNN):
 *
 *    H = J₁ · Σ_{⟨i,j⟩ ∈ nn}  S_i · S_j
 *      + J₂ · Σ_{⟨i,j⟩ ∈ nnn} S_i · S_j
 *
 *  Produced as an `irrep_heisenberg_t` so the apply callback is shared
 *  with the pure-NN case; internally stored as a concatenated bond list
 *  with per-bond effective coupling. Callers fill the nn and nnn bond
 *  arrays from @ref irrep_lattice_fill_bonds_nn and
 *  @ref irrep_lattice_fill_bonds_nnn.
 *
 *  The J₁-J₂ square Heisenberg is one of the standard frustrated-
 *  magnetism benchmarks; J₂/J₁ ≈ 0.5 sits near the proposed spin-liquid
 *  window and is the secondary target the 1.3 substrate supports. */
IRREP_API irrep_heisenberg_t *
irrep_heisenberg_j1j2_new(int num_sites,
                          int num_bonds_nn,  const int *nn_i,  const int *nn_j,  double J1,
                          int num_bonds_nnn, const int *nnn_i, const int *nnn_j, double J2);

/** @brief Build a spin-½ XY chain / lattice:
 *
 *    H = J · Σ_{⟨i,j⟩}  (S_i^x S_j^x + S_i^y S_j^y)  =  (J/2) · Σ  (S_i^+ S_j^− + h.c.)
 *
 *  No S^z-S^z term. S_z_total is still conserved, but the eigenvalue
 *  spectrum is continuous-like for large N (Bethe-ansatz soluble in 1D).
 *  Same apply callback as Heisenberg. */
IRREP_API irrep_heisenberg_t *
irrep_xy_new(int num_sites, int num_bonds,
             const int *bi, const int *bj,
             double J);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_HAMILTONIAN_H */
