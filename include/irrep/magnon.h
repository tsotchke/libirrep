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

/** @brief Compute the magnon dispersion ω_b(k) on top of a *general*
 *         collinear ground state in which sublattice α has its spin
 *         along σ_α · ẑ where σ_α ∈ {+1, -1}. AFM, ferrimagnetic,
 *         and any collinear order with ↑/↓ sublattice mixture is
 *         supported.
 *
 *  This is the *Bogoliubov-Colpa* generalisation of LSW to non-FM
 *  ground states. For each unique bond (i, j, t) with sublattice
 *  signs σ_i and σ_j, the bilinear-in-bosons Hamiltonian acquires
 *  three different structures:
 *
 *    - σ_i σ_j = +1 (parallel sublattices): same as the FM case —
 *      hopping a_i^† a_j with coefficient +S·J·exp(i k·t) and
 *      anomalous-DMI hopping −i S D_z exp(i k·t).
 *    - σ_i σ_j = -1 (antiparallel sublattices): the off-diagonal
 *      Heisenberg piece becomes a particle-creation pairing
 *      a_i a_j with coefficient +S·J·exp(i k·t), and the diagonal
 *      "self-energy" flips sign relative to FM (for Néel
 *      stability). The DMI z-component becomes ineffective at
 *      bilinear order on antiparallel-sublattice bonds.
 *
 *  We assemble a 2n × 2n Hermitian *Bogoliubov-de-Gennes* matrix
 *  M(k) on the basis (a_{1,k}, ..., a_{n,k}, a_{1,-k}^†, ...,
 *  a_{n,-k}^†) and diagonalise it via the Colpa (1978) prescription:
 *  Cholesky M = K^† K, then diagonalise K η K^† (Hermitian, with η =
 *  diag(I_n, −I_n)) — the positive eigenvalues are the n magnon
 *  energies ω_b(k). Returns IRREP_ERR_INVALID_ARG if M(k) is not
 *  positive-definite (the assumed ground state is unstable for the
 *  given bond list).
 *
 *  @param L                LSW handle
 *  @param sublattice_signs array of length n_sub with values +1 or
 *                          -1, indicating the spin direction of each
 *                          sublattice in the classical ground state
 *  @param kx, ky           momentum (cartesian)
 *  @param omega_out        caller buffer of size n_sub doubles —
 *                          magnon energies sorted ascending */
IRREP_API irrep_status_t irrep_magnon_dispersion_general(const irrep_magnon_lsw_t *L,
                                                          const int *sublattice_signs, double kx,
                                                          double ky, double *omega_out);

/** @brief Compute the magnon dispersion on a strip geometry: open
 *         boundary conditions along a₁ (Lx unit cells stacked) and
 *         periodic boundary conditions along a₂ (momentum k_y).
 *
 *  This is the *bulk-boundary correspondence* probe. For a 2D model
 *  with non-zero band Chern numbers, the strip-projected dispersion
 *  exhibits gapless chiral edge modes crossing every bulk gap. The
 *  net chirality of the edge modes at one edge equals the sum of
 *  Chern numbers of all bands below the gap (Hatsugai 1993).
 *
 *  The strip Hamiltonian H_strip(k_y) is an (Lx · n_sub) × (Lx ·
 *  n_sub) Hermitian matrix obtained by:
 *
 *    - keeping intra-cell bonds (delta_x = 0) on every cell;
 *    - connecting cells `i` and `i+1` via bonds with delta_x = ±1;
 *    - applying Bloch phase exp(i k_y · t_y) for the y-projection of
 *      every bond.
 *
 *  Bonds with |delta_x| > 1 are not currently supported (most NN /
 *  NNN bond lists fit within ±1).
 *
 *  @param L                LSW handle
 *  @param Lx               number of unit cells along a₁ (open BC).
 *                          Must be ≥ 2.
 *  @param ky               momentum along a₂ (cartesian)
 *  @param omega_out        caller buffer of size Lx · n_sub doubles —
 *                          eigenvalues sorted ascending
 *  @param edge_weight_out  caller buffer of size Lx · n_sub doubles
 *                          (may be NULL). On return, holds Σ_{x in
 *                          left half} |ψ_x|² for each mode — values
 *                          near 1 indicate left-edge localisation,
 *                          near 0 right-edge, near 0.5 bulk. */
IRREP_API irrep_status_t irrep_magnon_strip_dispersion(const irrep_magnon_lsw_t *L, int Lx,
                                                        double ky, double *omega_out,
                                                        double *edge_weight_out);

/** @brief Magnon thermal Hall conductivity κ_xy(T) in the
 *         Matsumoto-Murakami formulation (PRL 106, 197202, 2011 /
 *         PRB 89, 054420, 2014).
 *
 *  For magnons with Berry curvature Ω_b(k), the thermal Hall response is
 *
 *      κ_xy(T) = -(k_B² T / (ℏ V_uc · (2π)²))
 *                · Σ_b ∫_BZ d²k · c₂(n_B(ω_b/T)) · Ω_b(k)
 *
 *  where n_B(x) = 1/(e^x − 1) is the Bose-Einstein distribution and
 *
 *      c₂(g) = (1+g) [ln((1+g)/g)]² − (ln g)² − 2 Li₂(−g)
 *
 *  is a special function with the limits c₂(g) ≈ 2g for g → 0 and
 *  c₂(g) → π²/3 for g → ∞. The high-T limit makes κ_xy proportional
 *  to Σ_b C_b which vanishes for any closed band system; the low-T
 *  limit makes κ_xy → 0 exponentially because no magnons are
 *  populated. The peak temperature is set by the lowest-band gap or
 *  Dirac mass — typically O(D·sin(2π/3)) for kagome FM with DMI.
 *
 *  Returned value is in dimensionless natural units of (k_B² / ℏ)·T·A_uc⁻¹
 *  where T is also in natural units (k_B = ℏ = 1). For a real material
 *  with |J| = J_meV meV and unit-cell area A_uc_Å² in Å², multiply by
 *  J_meV² · 1.6e-22 / A_uc_Å² to get SI units of W/(K·m).
 *
 *  @param L            LSW handle
 *  @param T            temperature (in same energy units as ω)
 *  @param Nx, Ny       BZ integration grid (typical 50–200)
 *  @return κ_xy / (k_B² / ℏ) in natural units, or NaN on error. */
IRREP_API double irrep_magnon_thermal_hall_kxy(const irrep_magnon_lsw_t *L, double T, int Nx,
                                                int Ny);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_MAGNON_H */
