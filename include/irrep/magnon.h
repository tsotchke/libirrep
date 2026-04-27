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

/** @brief Per-bond coupling spec for a magnetic ground state.
 *
 *  The `delta_z` field defaults to 0 (zero-init) and is unused by the
 *  existing 2D LSW functions, which only consume `delta_x` and
 *  `delta_y`. It carries the inter-cell displacement along a₃ for
 *  the 3D extension (`irrep_magnon_dispersion_3d`). */
typedef struct {
    int    bi;        /**< source sublattice (within the magnetic unit cell) */
    int    bj;        /**< destination sublattice */
    int    delta_x;   /**< inter-cell displacement along a₁ */
    int    delta_y;   /**< inter-cell displacement along a₂ */
    int    delta_z;   /**< inter-cell displacement along a₃ (3D only; 0 for 2D) */
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

/** @brief Magnon density of states D(ω) via BZ histogram.
 *
 *  Samples ω_b(k) on an Nx × Ny grid in reduced momentum coordinates,
 *  bins the eigenvalues into n_bins between ω_min and ω_max, and
 *  normalises so that
 *
 *      ∫ D(ω) dω = n_sub                     (per unit cell)
 *
 *  i.e., the integrated DOS equals the number of bands per unit cell.
 *  Useful for comparison with magnetic specific heat (C_V(T) = ∫ ω
 *  D(ω) ∂_T n_B(ω/T) dω) and magnon-band van-Hove singularities. The
 *  DOS of a 2D FM with a quadratic Goldstone has a step edge at ω
 *  → 0⁺; saddle-point van-Hove peaks appear at zone-boundary
 *  energies; ω_max is the dispersion peak.
 *
 *  @param L           LSW handle
 *  @param Nx, Ny      BZ integration grid (typical 50–200)
 *  @param omega_min   lower edge of histogram window
 *  @param omega_max   upper edge of histogram window
 *  @param n_bins      number of bins; bin width = (ω_max−ω_min)/n_bins
 *  @param dos_out     caller buffer of size n_bins; D(ω) = states per
 *                     unit cell per (energy unit). Bin i is centred at
 *                     ω_min + (i + ½) · bin_width. */
IRREP_API irrep_status_t irrep_magnon_dos(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                           double omega_min, double omega_max, int n_bins,
                                           double *dos_out);

/** @brief Compute the band-resolved transverse dynamic spin structure
 *         factor S⊥_b(q) of a collinear-FM-along-z LSW magnet.
 *
 *  For a one-magnon excitation at momentum q on band b, the transverse
 *  spin operator S^- = √(2S)·a^† gives the matrix element
 *
 *      ⟨b, q | S^-_q | 0⟩ = √(2S) · Σ_α u_b(q)_α
 *
 *  in the Wannier-centred Bloch convention used by the rest of this
 *  module (sublattice positions absorbed into the cell-origin
 *  convention; caller can post-multiply the band-resolved result by
 *  intra-cell phase factors |Σ_α e^{i q·r_α} u_b(q)_α|² /
 *  |Σ_α u_b(q)_α|² for full INS prediction). The band-resolved
 *  spectral weight is
 *
 *      S⊥_b(q) = 2S · |Σ_α u_b(q)_α|²
 *
 *  This is the *intensity* an inelastic-neutron experiment measures
 *  along the magnon-band line ω = ω_b(q). Bands with `S⊥_b(q) = 0`
 *  are *dark* — invisible to transverse-channel neutron scattering at
 *  that q.
 *
 *  For a single-sublattice FM, |u_b|² = 1 by normalisation, so S⊥ =
 *  2S identically (q-independent intensity, all bands fully bright).
 *
 *  For the kagome FM at Γ, the lowest band is the uniform-mode
 *  Goldstone u_1 = (1, 1, 1)/√3 with |Σ u| = √3, giving S⊥ = 2S·3 =
 *  6S. The upper bands are orthogonal to the uniform mode and have
 *  S⊥ = 0 — dark at Γ. As q moves through the BZ, spectral weight
 *  redistributes between bands; the *integrated* sum-rule is
 *  Σ_b S⊥_b(q) = 2S·n_sub independent of q.
 *
 *  @param L            LSW handle
 *  @param qx, qy       momentum (cartesian)
 *  @param omega_out    caller buffer of size n_sub doubles — band
 *                      energies sorted ascending
 *  @param S_perp_out   caller buffer of size n_sub doubles — band-
 *                      resolved transverse spectral weight S⊥_b(q) */
IRREP_API irrep_status_t irrep_magnon_structure_factor(const irrep_magnon_lsw_t *L, double qx,
                                                        double qy, double *omega_out,
                                                        double *S_perp_out);

/** @brief Compute band Chern numbers on a 2D BZ slice at fixed k_z
 *         of a 3D LSW model. Useful for layered van der Waals magnets
 *         (CrI₃ stacks, kagome layer × c-axis), B20 chiral magnets in
 *         a FM-approximated ground state, or any 3D system whose
 *         topology localises on horizontal cuts of the 3D BZ.
 *
 *  Same FHS-plaquette algorithm as `irrep_magnon_chern`, but using
 *  `_dispersion_3d` to sample the (k_x, k_y) plane at fixed k_z. For
 *  a layered kagome topological FM, the Chern signature (-1, 0, +1)
 *  at k_z = 0 typically persists across the k_z BZ until band-gap
 *  closures at zone-boundary planes redistribute the integers.
 *
 *  @param L           LSW handle
 *  @param a3          third primitive vector
 *  @param kz          fixed k_z value of the slice
 *  @param Nx, Ny      integration grid (typical 50–200)
 *  @param chern_out   caller buffer of size n_sub doubles */
IRREP_API irrep_status_t irrep_magnon_chern_3d_slice_kz(const irrep_magnon_lsw_t *L,
                                                         const double a3[3], double kz, int Nx,
                                                         int Ny, double *chern_out);

/** @brief Compute the magnon dispersion ω_b(k) on a *3D* lattice.
 *
 *  The 2D `irrep_magnon_dispersion` only consumes a₁, a₂ and the
 *  bond's (delta_x, delta_y). For 3D systems — simple cubic FM,
 *  honeycomb-stack van der Waals magnets like CrI₃, B20 helimagnets
 *  approximated by their FM ground state — pass the third primitive
 *  vector a₃ here, and use `delta_z` on each bond to specify the
 *  inter-cell displacement along a₃. The bond translation vector
 *  becomes t = delta_x·a₁ + delta_y·a₂ + delta_z·a₃ in cartesian
 *  3D space.
 *
 *  All other parameters mirror `irrep_magnon_dispersion`. The
 *  returned ω_out and u_out have the same layout (n_sub
 *  eigenvalues / n_sub × n_sub eigenvectors).
 *
 *  @param L          LSW handle
 *  @param a3         third primitive vector (cartesian 3D)
 *  @param kx, ky, kz momentum (cartesian 3D)
 *  @param omega_out  caller buffer of size n_sub doubles
 *  @param u_out      caller buffer of size n_sub × n_sub complex doubles */
IRREP_API irrep_status_t irrep_magnon_dispersion_3d(const irrep_magnon_lsw_t *L,
                                                     const double a3[3], double kx, double ky,
                                                     double kz, double *omega_out,
                                                     double _Complex *u_out);

/** @brief Compute the *Abelian* Wilson loop spectrum θ_b(k_x) at fixed
 *         k_x. For each band b, the Wilson loop along k_y is
 *
 *      W_b(k_x) = ∏_{i=0}^{N_y-1} ⟨u_b(k_x, k_y_i) | u_b(k_x, k_y_{i+1})⟩
 *
 *  where the k_y values are taken at N_y equally-spaced points around
 *  the second BZ direction (b₂). The Wilson-loop phase is θ_b(k_x) =
 *  arg(W_b(k_x)) ∈ (−π, π].
 *
 *  This is the *Wannier-center flow* of the band b: as k_x sweeps
 *  across its BZ (-π/|b₁|, +π/|b₁|], θ_b(k_x) traces out a curve. The
 *  winding number of θ_b(k_x) — i.e., the number of times the curve
 *  crosses the branch cut at θ = ±π modulo 2π — equals the Chern
 *  number C_b of the band. This is the Soluyanov-Vanderbilt 2011
 *  formulation of band topology — a sharper probe than just the
 *  integer C_b because the *shape* of the curve distinguishes
 *  fragile vs stable topology and gives access to higher-order
 *  topological invariants.
 *
 *  For the canonical kagome FM with (-1, 0, +1) Chern signature:
 *  band 0 (lower) has θ(k_x) winding once; band 1 (mid) is flat;
 *  band 2 (upper) winds the opposite direction.
 *
 *  @param L          LSW handle
 *  @param kx         momentum along a₁ (cartesian)
 *  @param Ny         k_y discretisation (typical 32–128)
 *  @param theta_out  caller buffer of size n_sub doubles; on return
 *                    holds θ_b(k_x) ∈ (−π, π] for each band */
IRREP_API irrep_status_t irrep_magnon_wilson_spectrum(const irrep_magnon_lsw_t *L, double kx,
                                                       int Ny, double *theta_out);

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
