/* SPDX-License-Identifier: MIT */
/** @file dmi_hamiltonian.h
 *  @brief Spin-½ Dzyaloshinskii–Moriya Hamiltonian apply operator.
 *
 *  Companion to @ref hamiltonian.h. Where `irrep_heisenberg_t` handles
 *  the symmetric `S_i · S_j` term, this module handles the antisymmetric
 *  DMI term:
 *
 *      H_DMI = sum_{<i,j>}  D_ij · ( S_i × S_j )
 *
 *  Combined with the symmetry analyzer in @ref dmi.h, this closes the
 *  loop from group-theoretic prediction to actual exact diagonalisation:
 *
 *      1. analyzer:   magnetic-point-group → allowed D-vector subspace
 *      2. parameterisation: pick a magnitude on each allowed component
 *         (DFT input, or assumed for screening)
 *      3. apply:      build `irrep_dmi_hamiltonian_t` + couple with
 *                     `irrep_heisenberg_apply` to get the full
 *                     Heisenberg-DMI Hamiltonian on a finite cluster
 *      4. Lanczos:    pass a combined-apply callback to
 *                     `irrep_lanczos_eigvals_reorth` for the GS energy
 *                     and ordering vector.
 *
 *  Spin-½ basis convention: `bit i of s = 1` ↔ spin-up at site i.
 *  S_z = (bit_i − ½), S_+ = |↑⟩⟨↓|, S_− = |↓⟩⟨↑|.
 *
 *  DMI matrix element on a bond `(a, b)` for an arbitrary 3-vector
 *  D = (D_x, D_y, D_z):
 *
 *      D · (S_a × S_b)
 *        = D_x (S_a^y S_b^z − S_a^z S_b^y)
 *        + D_y (S_a^z S_b^x − S_a^x S_b^z)
 *        + D_z (S_a^x S_b^y − S_a^y S_b^x)
 *
 *  In raising/lowering basis (S_+ = S_x + i S_y, S_− = S_x − i S_y):
 *      S_x = (S_+ + S_−)/2,    S_y = (S_+ − S_−)/(2i)
 *
 *  Substituting and collecting terms:
 *
 *      D · (S_a × S_b) =
 *           ½ (D_x − i D_y) [ S_a^z S_b^+ − S_a^+ S_b^z ]
 *         + ½ (D_x + i D_y) [ S_a^- S_b^z − S_a^z S_b^- ]
 *         + (i/2) D_z       [ S_a^+ S_b^- − S_a^- S_b^+ ]
 *
 *  Each term flips at most two bits and is straightforward to evaluate
 *  on the spin-½ computational basis. */
#ifndef IRREP_DMI_HAMILTONIAN_H
#define IRREP_DMI_HAMILTONIAN_H

#include <complex.h>
#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Opaque DMI-Hamiltonian handle. Built from a per-bond D-vector
 *         list; applied via the Lanczos-compatible @ref irrep_dmi_apply
 *         callback. */
typedef struct irrep_dmi_hamiltonian irrep_dmi_hamiltonian_t;

/** @brief Construct a DMI Hamiltonian on `num_sites` spin-½ sites.
 *
 *  The bond list `{(bi[b], bj[b]) : b ∈ [0, num_bonds)}` is copied;
 *  each bond carries a 3-vector D = (D_x[b], D_y[b], D_z[b]). The
 *  caller is free to free the input arrays after this call.
 *
 *  Choosing the orientation of D matters: `D_ij = -D_ji` is the
 *  standard antisymmetric convention. The caller is responsible for
 *  feeding the bond list with consistent orientation (typically
 *  `bi[b] < bj[b]` and a single D-vector per undirected bond).
 *
 *  @return Owned handle, or NULL on OOM / invalid input. */
IRREP_API irrep_dmi_hamiltonian_t *irrep_dmi_hamiltonian_new(int num_sites, int num_bonds,
                                                              const int *bi, const int *bj,
                                                              const double *D_x, const double *D_y,
                                                              const double *D_z);

/** @brief Release a handle returned by @ref irrep_dmi_hamiltonian_new. */
IRREP_API void irrep_dmi_hamiltonian_free(irrep_dmi_hamiltonian_t *H);

/** @brief Apply callback compatible with @ref irrep_lanczos_eigvals.
 *
 *  Computes `out = H_DMI · psi` on length-`2^num_sites` complex
 *  vectors. The result is generally complex (DMI is **not** purely
 *  real even on a real-coupling-coefficient input). Combined-apply
 *  patterns (e.g., Heisenberg + DMI) sum the contributions from each
 *  apply callback into a shared `out` buffer.
 *
 *  Time complexity: `O(num_bonds · 2^num_sites)` per call. Memory:
 *  writes `out`, reads `psi`. Thread-safe across disjoint
 *  `(psi, out)` pairs sharing the same `opaque`. */
IRREP_API void irrep_dmi_apply(const double _Complex *psi, double _Complex *out, void *opaque);

/** @brief Number of sites the handle was built with. */
IRREP_API int irrep_dmi_hamiltonian_num_sites(const irrep_dmi_hamiltonian_t *H);

/** @brief Hilbert-space dimension (`2^num_sites`). Convenience for
 *         sizing Lanczos buffers. */
IRREP_API long long irrep_dmi_hamiltonian_dim(const irrep_dmi_hamiltonian_t *H);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_DMI_HAMILTONIAN_H */
