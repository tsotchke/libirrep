/* SPDX-License-Identifier: MIT */
/** @file magnon.h
 *  @brief Linearized spin-wave theory: magnon dispersion ω(k) and Berry
 *         curvature for an arbitrary Heisenberg + DMI + anisotropy
 *         Hamiltonian on a libirrep lattice.
 *
 *  Closes the loop from libirrep's *algebraic* layer (symmetry analyzer,
 *  bilinear bond Hamiltonian) to *measurable* physics: magnon dispersion,
 *  Berry curvature, Chern numbers, topological edge modes. These are the
 *  observables that inelastic neutron scattering, thermal Hall transport,
 *  and Brillouin-light-scattering measurements actually probe.
 *
 *  ## Theory
 *
 *  For a magnetically ordered ground state with collinear spins along an
 *  ordering axis, linearised spin-wave theory expands `S_i = S - a_i^† a_i`
 *  (Holstein-Primakoff) about the classical ground state and keeps the
 *  bilinear bosonic Hamiltonian
 *
 *      H_LSW = sum_k psi_k^† H(k) psi_k
 *
 *  where `psi_k = (a_{1,k}, ..., a_{n,k})` runs over the magnetic
 *  sublattices. `H(k)` is a Hermitian n×n matrix per k whose eigenvalues
 *  are the magnon energies ω_b(k). The eigenvectors |u_b(k)⟩ carry
 *  Berry connection
 *
 *      A_a^b(k) = i ⟨u_b(k)| ∂_a |u_b(k)⟩
 *
 *  whose curl is the Berry curvature
 *
 *      Ω_b(k) = ∂_kx A_y^b - ∂_ky A_x^b
 *
 *  Integrated over the Brillouin zone divided by 2π gives the band
 *  Chern number — the magnon analog of the QHE invariant. Non-zero
 *  Chern numbers on magnon bands give rise to a *thermal* Hall effect
 *  that has been measured in Lu₂V₂O₇, Fe₃Sn₂, and others.
 *
 *  ## Conventions
 *
 *  - Spins of magnitude `S` (caller-supplied; spin-½ is `S = 0.5`).
 *  - Classical ground state: collinear ferromagnet (the simplest case;
 *    AFM and helimagnet ground states require a unit-cell expansion
 *    handled by the caller — pass an enlarged "magnetic unit cell"
 *    bond list).
 *  - Heisenberg + DMI + uniaxial anisotropy K_z (S_i^z)^2:
 *
 *      H = sum_<ij> J_ij S_i · S_j + sum_<ij> D_ij · (S_i × S_j)
 *          - K_z sum_i (S_i^z)^2
 *
 *    K_z > 0 favours easy-axis along z (FM stabilisation).
 *
 *  - Lattice + bond list comes from the caller; this module is
 *    geometry-agnostic, only consumes (n_sublattices, bond list with
 *    per-bond couplings, primitive vectors).
 */
#ifndef IRREP_MAGNON_H
#define IRREP_MAGNON_H

#include <complex.h>
#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Per-bond coupling spec for a magnetic ground state. */
typedef struct {
    int    bi;        /**< source sublattice (within the magnetic unit cell) */
    int    bj;        /**< destination sublattice */
    int    delta_x;   /**< inter-cell displacement along a₁ */
    int    delta_y;   /**< inter-cell displacement along a₂ */
    double J;         /**< Heisenberg coupling magnitude */
    double D[3];      /**< DMI vector D_ij */
} irrep_magnon_bond_t;

/** @brief Linearised-spin-wave handle. Built from the lattice geometry
 *         + bond list with per-bond (J, D); applied to compute ω(k) and
 *         Berry curvature on demand. Currently 2D (caller passes
 *         primitive vectors a₁, a₂). 3D is a straightforward extension
 *         once the caller supplies a₃; deferred for v1.4. */
typedef struct irrep_magnon_lsw irrep_magnon_lsw_t;

/** @brief Construct an LSW handle.
 *
 *  @param n_sub          number of magnetic sublattices in the unit cell
 *  @param S              total spin per site (e.g., 0.5 for spin-½, 1 for S=1)
 *  @param a1             2D primitive vector a₁ (cartesian)
 *  @param a2             2D primitive vector a₂
 *  @param n_bonds        number of bonds in the bond list
 *  @param bonds          array of irrep_magnon_bond_t describing the
 *                        Heisenberg + DMI coupling on each bond. Bonds
 *                        are NOT pre-canonicalised; pass each undirected
 *                        bond once, in either orientation (the LSW
 *                        construction is symmetric under bond reversal
 *                        for J, antisymmetric for D).
 *  @param Kz             uniaxial anisotropy strength (per site, positive
 *                        favours easy-axis along z)
 *  @return owned handle, or NULL on OOM / invalid input. */
IRREP_API irrep_magnon_lsw_t *irrep_magnon_lsw_new(int n_sub, double S, const double a1[2],
                                                   const double a2[2], int n_bonds,
                                                   const irrep_magnon_bond_t *bonds, double Kz);

/** @brief Release an LSW handle. */
IRREP_API void irrep_magnon_lsw_free(irrep_magnon_lsw_t *L);

/** @brief Compute magnon dispersion ω_b(k) at a single k point.
 *
 *  Diagonalises the n_sub × n_sub LSW Hamiltonian H(k); writes the
 *  n_sub eigenvalues sorted ascending into omega_out. Eigenvectors are
 *  written into u_out (n_sub × n_sub complex matrix, row-major; row b
 *  is the eigenvector at band b).
 *
 *  Result is real (ω) and complex (u) per band. Negative or zero
 *  eigenvalues indicate the assumed FM ground state is *not* the true
 *  classical ground state — the caller should pick a different ordering
 *  ansatz.
 *
 *  @param L          LSW handle
 *  @param kx, ky     momentum (cartesian, dimensions of inverse length)
 *  @param omega_out  caller buffer of size n_sub doubles
 *  @param u_out      caller buffer of size n_sub * n_sub complex doubles */
IRREP_API irrep_status_t irrep_magnon_dispersion(const irrep_magnon_lsw_t *L, double kx,
                                                  double ky, double *omega_out,
                                                  double _Complex *u_out);

/** @brief Compute the Berry curvature Ω_b(k) of every band at a single
 *         k point via the gauge-invariant 4-point formula
 *
 *      F_b(k) = arg ⟨u_b(k₁)|u_b(k₂)⟩ ⟨u_b(k₂)|u_b(k₃)⟩
 *                  ⟨u_b(k₃)|u_b(k₄)⟩ ⟨u_b(k₄)|u_b(k₁)⟩
 *
 *  on a small plaquette at k of side `delta_k`. The curvature is
 *  `Ω_b ≈ F_b / delta_k²`. This is gauge-invariant (no fixing-of-phase
 *  needed) and converges quadratically in delta_k for smooth bands.
 *
 *  @param L           LSW handle
 *  @param kx, ky      central momentum
 *  @param delta_k     plaquette half-side (typical 1e-3 .. 1e-2)
 *  @param berry_out   caller buffer of size n_sub doubles */
IRREP_API irrep_status_t irrep_magnon_berry(const irrep_magnon_lsw_t *L, double kx, double ky,
                                             double delta_k, double *berry_out);

/** @brief Integrate Berry curvature over the Brillouin zone to get the
 *         band Chern number. Uses a uniform Nx × Ny grid in reduced
 *         momentum coordinates; integrates the gauge-invariant plaquette
 *         flux over each grid cell and sums. The result for a topologically
 *         non-trivial band converges to an integer.
 *
 *  @param L            LSW handle
 *  @param Nx, Ny       integration grid (typical 50–200)
 *  @param chern_out    caller buffer of size n_sub doubles; on return
 *                      holds Chern number per band (should be near-integer
 *                      for gapped bands). */
IRREP_API irrep_status_t irrep_magnon_chern(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                             double *chern_out);

/** @brief Read out the number of magnon bands (= number of magnetic
 *         sublattices). */
IRREP_API int irrep_magnon_lsw_num_bands(const irrep_magnon_lsw_t *L);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_MAGNON_H */
