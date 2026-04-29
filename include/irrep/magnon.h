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

/** @brief Factory callback for Chern parameter sweeps. Constructs a
 *         fresh LSW handle from a single scalar parameter value. The
 *         caller of `irrep_magnon_chern_sweep` will free the returned
 *         handle via `irrep_magnon_lsw_free` after each Chern call.
 *
 *  @param param   the parameter value
 *  @param ud      user-data pointer forwarded from the sweep call
 *  @return        a fresh LSW handle, or NULL on construction failure */
typedef irrep_magnon_lsw_t *(*irrep_magnon_lsw_factory_fn)(double param, void *ud);

/** @brief Sweep band Chern numbers across a 1D parameter axis. For
 *         each value `params[i]`, calls `factory(params[i], ud)` to
 *         build an LSW handle, runs `irrep_magnon_chern` on a
 *         Nx × Ny BZ grid, and stores the per-band Chern numbers in
 *         `chern_out`. Frees each LSW handle after use.
 *
 *  Output `chern_out` is row-major n_params × n_sub: entry
 *  `chern_out[i * n_sub + b]` is the Chern number of band b at
 *  parameter value `params[i]`.
 *
 *  Use case: parameter scans of a topological phase diagram (vary D,
 *  J', K_z, applied field, etc.). Combine with
 *  `irrep_magnon_chern_detect_boundaries` to find Chern-changing
 *  transitions automatically.
 *
 *  @param factory     factory callback (see typedef above)
 *  @param ud          user data forwarded to factory
 *  @param params      array of n_params parameter values
 *  @param n_params    number of parameter samples (≥ 1)
 *  @param n_sub       number of bands per LSW handle (must match
 *                     factory output)
 *  @param Nx, Ny      BZ grid for the Chern integration
 *  @param chern_out   caller buffer of size n_params · n_sub doubles */
IRREP_API irrep_status_t irrep_magnon_chern_sweep(irrep_magnon_lsw_factory_fn factory, void *ud,
                                                    const double *params, int n_params, int n_sub,
                                                    int Nx, int Ny, double *chern_out);

/** @brief Detect Chern-integer-changing boundaries in a 1D parameter
 *         sweep. Given a Chern table `chern_in` (n_params × n_sub),
 *         returns the parameter-axis indices where the rounded integer
 *         Chern signature changes between adjacent samples — i.e.
 *         where `round(C[i, *])` and `round(C[i+1, *])` differ.
 *
 *  Each boundary index i marks the start of an interval `[params[i],
 *  params[i+1]]` containing a Chern-changing gap closing. Refine with
 *  a denser sweep on that sub-interval if a sharper boundary location
 *  is needed (the resolution of this detector is the parameter grid
 *  spacing).
 *
 *  Output `boundary_indices_out[k]` is the lower-index `i` of the
 *  k-th detected boundary. At most `max_boundaries` are written; if
 *  more boundaries exist than the buffer can hold, the function
 *  returns `IRREP_ERR_INVALID_ARG`.
 *
 *  @param chern_in              n_params × n_sub Chern table (row-major)
 *  @param n_params, n_sub       table dimensions
 *  @param boundary_indices_out  caller buffer for boundary indices
 *  @param max_boundaries        size of `boundary_indices_out`
 *  @param n_boundaries_out      number of boundaries found */
IRREP_API irrep_status_t irrep_magnon_chern_detect_boundaries(const double *chern_in, int n_params,
                                                                int n_sub,
                                                                int *boundary_indices_out,
                                                                int max_boundaries,
                                                                int *n_boundaries_out);

/** @brief Read out the number of magnon bands (= number of magnetic
 *         sublattices). */
IRREP_API int irrep_magnon_lsw_num_bands(const irrep_magnon_lsw_t *L);

/** @brief Inelastic-neutron Q-ω intensity map I(q, ω) for a path of
 *         momenta — the direct experimental observable that neutron
 *         papers publish.
 *
 *  For each path momentum q_i, computes the band-resolved transverse
 *  structure factor and sums Lorentzian-broadened delta functions:
 *
 *      I(q_i, ω_j) = Σ_b S_⊥_b(q_i) · L(ω_j − ω_b(q_i); η)
 *
 *  with L(x; η) = (1/π) · η / (x² + η²) the unit-area Lorentzian and
 *  η the user-supplied energy resolution (FWHM / 2). The result is a
 *  Q × ω heatmap directly comparable to the canonical inelastic-
 *  neutron-scattering plots.
 *
 *  Caller supplies:
 *    - a momentum path through the BZ as an array of (qx, qy) pairs;
 *    - an energy axis: omega_min, omega_max, n_omega bins;
 *    - the broadening η (typical 0.02 - 0.05 in J-units).
 *
 *  Output: row-major intensity[i_q · n_omega + j_omega], in
 *  natural units (S_⊥ × Lorentzian).
 *
 *  @param L           LSW handle
 *  @param qpath       2D array [n_q][2] of (qx, qy) pairs
 *  @param n_q         number of path points
 *  @param omega_min   lower edge of energy axis
 *  @param omega_max   upper edge of energy axis
 *  @param n_omega     number of energy bins
 *  @param eta         Lorentzian half-width
 *  @param intensity_out caller buffer of size n_q · n_omega doubles */
IRREP_API irrep_status_t irrep_magnon_neutron_qomega_map(const irrep_magnon_lsw_t *L,
                                                          const double (*qpath)[2], int n_q,
                                                          double omega_min, double omega_max,
                                                          int n_omega, double eta,
                                                          double *intensity_out);

/** @brief Locate the *softest mode* — the k-point + band where the
 *         dispersion attains its global minimum across a sampled BZ.
 *
 *  Useful for identifying:
 *
 *    - **Ordering wavevectors**: an AFM-canting / helimagnet
 *      instability soft-modes at a finite k* before condensing into
 *      a new ordered phase. The softest-mode k tells you what
 *      ground-state ansatz to switch to.
 *
 *    - **Dirac points**: when the gap between bands b and b+1
 *      closes at some k*, the softest mode of the upper band is
 *      degenerate with the highest of the lower — a Dirac touching.
 *
 *    - **Flat bands**: if the softest-mode value is reached on a
 *      large region of the BZ rather than a single k-point, the
 *      band is flat at that energy.
 *
 *  The `exclude_below` knob skips numerical-noise Goldstone modes:
 *  set to a small ε ~ 10⁻⁶ to find the *next* softest above a true
 *  Goldstone; set to -∞ for the global minimum.
 *
 *  Returns (kx, ky) at floating-point precision of the underlying
 *  grid (Nx × Ny resolution).
 *
 *  @param L              LSW handle
 *  @param Nx, Ny         BZ grid (typical 50-200)
 *  @param exclude_below  modes with ω ≤ this are skipped
 *  @param kx_out, ky_out caller buffers for the soft-mode momentum
 *  @param omega_out      caller buffer for the soft-mode energy
 *  @param band_out       caller buffer for the band index (0..n_sub-1)
 *                        of the softest mode */
IRREP_API irrep_status_t irrep_magnon_softest_mode(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                                    double exclude_below, double *kx_out,
                                                    double *ky_out, double *omega_out,
                                                    int *band_out);

/** @brief Scan the dispersion across the BZ and return the global
 *         minimum and maximum magnon energies.
 *
 *  The two outputs are the *most-quoted* dispersion summary
 *  statistics:
 *
 *    - ω_min: the *spin gap* (= 0 for unbroken-symmetry Goldstone
 *      modes; > 0 for systems with anisotropy K_z, DMI gap, or
 *      AFM-canting gap). Distinguishes gapped from gapless ground
 *      states in a single number.
 *
 *    - ω_max: the *band-top* energy. Bandwidth = ω_max − ω_min is
 *      the dispersion span and sets the energy scale for thermal
 *      population (T ~ bandwidth → all magnon bands are thermally
 *      populated).
 *
 *  The user can exclude numerical-noise Goldstone modes by setting
 *  `exclude_below`: any sampled mode with ω ≤ exclude_below is
 *  skipped. For a clean LSW with a true gap, set
 *  exclude_below = 0; for systems with a Goldstone mode that you
 *  want to skip when reporting the *next* gap, set
 *  exclude_below = some small ε (typical 10⁻⁶).
 *
 *  @param L              LSW handle
 *  @param Nx, Ny         BZ grid (typical 50-200)
 *  @param exclude_below  modes with ω ≤ this are skipped
 *  @param omega_min_out  on return, lowest-band energy in the
 *                        sampled BZ above exclude_below
 *  @param omega_max_out  on return, highest-band energy in the
 *                        sampled BZ */
IRREP_API irrep_status_t irrep_magnon_band_extrema(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                                    double exclude_below, double *omega_min_out,
                                                    double *omega_max_out);

/** @brief Compute the band-resolved magnon-dispersion Hessian
 *         H_ij(k) = ∂²ω_b/∂k_i ∂k_j at a single k point.
 *
 *  At a band extremum (Γ for collinear FM), the inverse Hessian is
 *  the *effective mass tensor*:
 *
 *      m*_b_ij(k_min) = ℏ² · (H_b⁻¹)_ij
 *
 *  For an isotropic FM, the Hessian at Γ is diag(2D, 2D) where D is
 *  the *spin stiffness* — the prefactor in the Bloch T^{d/2} law:
 *
 *      M(T) − M(0) = −ζ(d/2) / Γ(d/2) · (T / 4πD)^{d/2}
 *
 *  reported in inelastic-neutron papers as the slope ω/k² in the
 *  low-k Goldstone limit. Anisotropic FMs have non-degenerate
 *  H_xx ≠ H_yy and possibly non-zero H_xy.
 *
 *  At saddle points (van-Hove singularities), H has mixed-sign
 *  eigenvalues. At band maxima, both eigenvalues are negative.
 *
 *  5-point central-difference estimator with step `h` (typical
 *  10⁻³ - 10⁻²). Caveats around band crossings same as for
 *  `_group_velocity`.
 *
 *  @param L         LSW handle
 *  @param kx, ky    central momentum
 *  @param h         finite-difference step
 *  @param hxx_out   ∂²ω/∂k_x² per band (size n_sub)
 *  @param hyy_out   ∂²ω/∂k_y² per band (size n_sub)
 *  @param hxy_out   ∂²ω/∂k_x ∂k_y per band (size n_sub) */
IRREP_API irrep_status_t irrep_magnon_hessian(const irrep_magnon_lsw_t *L, double kx, double ky,
                                               double h, double *hxx_out, double *hyy_out,
                                               double *hxy_out);

/** @brief AFM-aware two-magnon S⁽²⁾(q, ω). Same as
 *         `_two_magnon_qomega` but uses the Bogoliubov-aware
 *         structure factor `_structure_factor_general` instead of
 *         the FM-only Hermitian form.
 *
 *  Required for AFM materials where the Bogoliubov pairing is
 *  essential to the 2-magnon spectral weight (e.g., cuprate INS
 *  continuum above the band-top includes AFM 2-magnon
 *  contributions).
 *
 *  Same parameters and output format as `_two_magnon_qomega` plus
 *  the `sublattice_signs` array. */
IRREP_API irrep_status_t irrep_magnon_two_magnon_qomega_general(const irrep_magnon_lsw_t *L,
                                                                  const int *sublattice_signs,
                                                                  const double (*qpath)[2],
                                                                  int n_q, int Nx, int Ny,
                                                                  double omega_min,
                                                                  double omega_max, int n_omega,
                                                                  double eta,
                                                                  double *intensity_out);

/** @brief Two-magnon dynamical structure factor S⁽²⁾(q, ω) — the
 *         q-resolved 2-magnon continuum, a real beyond-1-magnon
 *         INS observable.
 *
 *      S⁽²⁾(q, ω) = (1/N) Σ_{k, b₁, b₂}
 *                   |M_{b₁}(k)|² · |M_{b₂}(q−k)|²
 *                   · L(ω − ω_{b₁}(k) − ω_{b₂}(q−k); η)
 *
 *  with M_b(k) = Σ_α u_b(k)_α (the uniform-mode form factor used in
 *  the 1-magnon `_structure_factor`).
 *
 *  Goes beyond 1-magnon LSW: this is a TWO-PARTICLE observable
 *  computed at the tree level. Combined with 1-magnon
 *  `_neutron_qomega_map`, gives the full INS prediction for
 *  comparison with high-energy continuum data above the 1-magnon
 *  band-top.
 *
 *  For cuprates: above the 1-magnon ω_max ≈ 2J·S, INS data shows
 *  a continuum extending to 4J·S — this 2-magnon S⁽²⁾(q, ω) is
 *  what fits that continuum.
 *
 *  Output is row-major n_q × n_omega; same format as
 *  `_neutron_qomega_map`.
 *
 *  @param L              LSW handle
 *  @param qpath          n_q × 2 array of (qx, qy) momenta
 *  @param n_q            number of q-points
 *  @param Nx, Ny         BZ integration grid (cost is quadratic)
 *  @param omega_min      lower edge of energy axis
 *  @param omega_max      upper edge of energy axis
 *  @param n_omega        energy bins
 *  @param eta            Lorentzian half-width
 *  @param intensity_out  caller buffer of size n_q × n_omega doubles. */
IRREP_API irrep_status_t irrep_magnon_two_magnon_qomega(const irrep_magnon_lsw_t *L,
                                                          const double (*qpath)[2], int n_q,
                                                          int Nx, int Ny, double omega_min,
                                                          double omega_max, int n_omega,
                                                          double eta, double *intensity_out);

/** @brief AFM-aware one-magnon S⁽¹⁾(q, ω) — Bogoliubov-Colpa
 *         structure factor in the (u, v) basis. Same parameters as
 *         `_one_magnon_qomega` plus `sublattice_signs` (σ_α ∈
 *         {+1, -1}).
 *
 *  Required for AFM materials where the Bogoliubov pairing controls
 *  the spectral weight (cuprate La₂CuO₄, kagome 120° Néel, ferri).
 *
 *  Sum rule: same — ∫ S⁽¹⁾(q, ω) dω = 2S · n_sub.
 *
 *  For FM-track ground states (single sublattice, or all signs +1),
 *  use `_neutron_qomega_map`. */
IRREP_API irrep_status_t irrep_magnon_one_magnon_qomega_general(const irrep_magnon_lsw_t *L,
                                                                  const int *sublattice_signs,
                                                                  const double (*qpath)[2],
                                                                  int n_q, double omega_min,
                                                                  double omega_max, int n_omega,
                                                                  double eta,
                                                                  double *intensity_out);

/** @brief Total dynamical structure factor S(q, ω) = S⁽¹⁾ + S⁽²⁾
 *         on a q-path — the full LSW prediction for inelastic neutron
 *         scattering intensity in (q, ω) space.
 *
 *      S(q, ω) = Σ_b S_⊥_b(q) · (η/π) / ((ω − ω_b(q))² + η²)         [1-magnon]
 *              + (1/N_BZ) Σ_{k, b₁, b₂} S_⊥_{b₁}(k) · S_⊥_{b₂}(q−k)
 *                            · (η/π) / ((ω − ω_{b₁}(k) − ω_{b₂}(q−k))² + η²)   [2-magnon]
 *
 *  Equivalent to calling `_one_magnon_qomega` and `_two_magnon_qomega`
 *  separately and summing — provided as a single call for the common
 *  case where the full (1 + 2)-magnon prediction is wanted.
 *
 *  For AFM/ferri ground states use `_dynamical_structure_factor_general`.
 *
 *  @param L              LSW handle
 *  @param qpath          n_q × 2 array of (qx, qy) momenta
 *  @param n_q            number of q-points
 *  @param Nx, Ny         BZ grid for 2-magnon convolution
 *  @param omega_min      lower edge of energy axis
 *  @param omega_max      upper edge of energy axis
 *  @param n_omega        energy bins
 *  @param eta            Lorentzian half-width
 *  @param intensity_out  caller buffer of size n_q × n_omega doubles. */
IRREP_API irrep_status_t irrep_magnon_dynamical_structure_factor(const irrep_magnon_lsw_t *L,
                                                                   const double (*qpath)[2],
                                                                   int n_q, int Nx, int Ny,
                                                                   double omega_min,
                                                                   double omega_max, int n_omega,
                                                                   double eta,
                                                                   double *intensity_out);

/** @brief Finite-T 1-magnon dynamical structure factor S^(1)(q, ω, T)
 *         — Stokes channel with the (1 + n_B(ω, T)) Bose enhancement
 *         that real inelastic-neutron experiments measure.
 *
 *      S^(1)(q, ω; T) = Σ_b (1 + n_B(ω_b(q), T)) · S_⊥_b(q)
 *                          · (η/π) / ((ω − ω_b(q))² + η²)
 *
 *  At T → 0, n_B → 0 and this reduces to `_neutron_qomega_map`.
 *  At higher T, soft modes (ω_b(q) ≪ k_B T) are strongly amplified —
 *  a real effect in finite-T INS that affects scaling fits and
 *  sub-bandwidth weight estimates.
 *
 *  This handles the Stokes (ω > 0, magnon creation) channel; the
 *  anti-Stokes (magnon-annihilation) side is suppressed by n_B / (1
 *  + n_B) = exp(-ω/T) and is not returned here.
 *
 *  For AFM ground states use `_dynamical_structure_factor_T_general`.
 *
 *  @param L              LSW handle
 *  @param qpath          n_q × 2 momenta
 *  @param n_q            number of q-points
 *  @param omega_min      lower edge of energy axis (ω > 0 for Stokes)
 *  @param omega_max      upper edge of energy axis
 *  @param n_omega        number of energy bins
 *  @param eta            Lorentzian half-width
 *  @param T              temperature in J units (k_B = 1)
 *  @param intensity_out  caller buffer of size n_q × n_omega doubles. */
IRREP_API irrep_status_t irrep_magnon_dynamical_structure_factor_T(
    const irrep_magnon_lsw_t *L, const double (*qpath)[2], int n_q, double omega_min,
    double omega_max, int n_omega, double eta, double T, double *intensity_out);

/** @brief AFM-aware finite-T S^(1)(q, ω, T). Same as
 *         `_dynamical_structure_factor_T` but uses Bogoliubov-Colpa
 *         structure factors. */
IRREP_API irrep_status_t irrep_magnon_dynamical_structure_factor_T_general(
    const irrep_magnon_lsw_t *L, const int *sublattice_signs, const double (*qpath)[2], int n_q,
    double omega_min, double omega_max, int n_omega, double eta, double T,
    double *intensity_out);

/** @brief Anti-Stokes channel of the finite-T 1-magnon spectral
 *         function — magnon-annihilation contribution at ω < 0.
 *
 *      S^(1)_anti(q, ω; T) = Σ_b n_B(ω_b(q), T) · S_⊥_b(q)
 *                              · (η/π) / ((ω + ω_b(q))² + η²)
 *
 *  At T = 0 returns identically 0 (no thermal magnons to annihilate).
 *  Together with `_dynamical_structure_factor_T` (Stokes), gives the
 *  full-window spectral function consistent with detailed balance:
 *
 *      S(q, ω, T) / S(q, -ω, T) = exp(ω / T)
 *
 *  For TR-symmetric ground states the anti-Stokes peak is at
 *  ω = -ω_b(q); for DMI-broken systems this would shift to
 *  -ω_b(-q), which the FM-track function does NOT account for —
 *  for chiral magnets use the AFM-aware variant which works in the
 *  Bogoliubov basis.
 *
 *  Caller is expected to choose ω_min < 0 < ω_max to see both
 *  channels in a single histogram window. */
IRREP_API irrep_status_t irrep_magnon_dynamical_structure_factor_T_anti_stokes(
    const irrep_magnon_lsw_t *L, const double (*qpath)[2], int n_q, double omega_min,
    double omega_max, int n_omega, double eta, double T, double *intensity_out);

/** @brief AFM-aware anti-Stokes channel — same as
 *         `_dynamical_structure_factor_T_anti_stokes` but with
 *         Bogoliubov structure factors. */
IRREP_API irrep_status_t irrep_magnon_dynamical_structure_factor_T_anti_stokes_general(
    const irrep_magnon_lsw_t *L, const int *sublattice_signs, const double (*qpath)[2], int n_q,
    double omega_min, double omega_max, int n_omega, double eta, double T,
    double *intensity_out);

/** @brief AFM-aware total S(q, ω) — same as `_dynamical_structure_factor`
 *         but uses Bogoliubov-aware structure factors for both 1- and
 *         2-magnon contributions. Required for AFM/ferri materials. */
IRREP_API irrep_status_t irrep_magnon_dynamical_structure_factor_general(
    const irrep_magnon_lsw_t *L, const int *sublattice_signs, const double (*qpath)[2], int n_q,
    int Nx, int Ny, double omega_min, double omega_max, int n_omega, double eta,
    double *intensity_out);

/** @brief Fit Heisenberg couplings J's to observed band-resolved
 *         dispersion data via gradient descent with backtracking.
 *
 *  Given a bond list with initial-guess J values and a set of
 *  measured points (q_i, ω_i, b_i), iteratively adjusts each bond's J
 *  to minimise χ² = Σ_i (ω_b_i(q_i; J) − ω_i)². The DMI fields D[3]
 *  are NOT fit — pass them as known constants in the bond list.
 *
 *  Algorithm: gradient descent with backtracking line search. Each
 *  iteration computes the gradient ∂χ²/∂J_k by finite-differencing
 *  with ε = 1e-4 (rebuilds the LSW handle inside). Step size adapts
 *  via (×1.1 on accept, ×0.5 on reject) backtracking.
 *
 *  Each bond's J is fit INDEPENDENTLY. To share a J across symmetry-
 *  related bonds, callers should manually average post-hoc, or set
 *  identical initial guesses and let the symmetry of the data drive
 *  the fit toward identical optima (typically does within ~1e-4 in
 *  noise-free synthetic data).
 *
 *  @param n_sub        number of magnetic sublattices
 *  @param S            spin per site
 *  @param a1, a2       primitive vectors
 *  @param bonds        bond list — J's are starting guesses, overwritten with fit
 *  @param n_bonds      number of bonds
 *  @param Kz           uniaxial anisotropy (held fixed)
 *  @param q_obs        n_obs × 2 momenta of observations
 *  @param omega_obs    n_obs observed dispersion energies
 *  @param band_obs     n_obs band indices (0 .. n_sub-1)
 *  @param n_obs        number of observations
 *  @param max_iter     maximum gradient-descent iterations
 *  @param tol          convergence tolerance on |Δχ²|
 *  @param chi2_out     final χ² (caller may pass NULL)
 *  @return IRREP_OK on success, IRREP_ERR_* otherwise. */
IRREP_API irrep_status_t irrep_magnon_fit_J(int n_sub, double S, const double a1[2],
                                              const double a2[2],
                                              irrep_magnon_bond_t *bonds, int n_bonds,
                                              double Kz, const double (*q_obs)[2],
                                              const double *omega_obs, const int *band_obs,
                                              int n_obs, int max_iter, double tol,
                                              double *chi2_out);

/** @brief Fit Heisenberg J + DMI D[3] for every bond simultaneously.
 *         Same gradient-descent / backtracking algorithm as `_fit_J`,
 *         but with 4× the parameters: each bond contributes one J
 *         and three DMI components (D[0], D[1], D[2]) to the fit.
 *
 *  Use this when the dispersion data has enough q-coverage to
 *  constrain the antisymmetric exchange: typically requires INS
 *  data along a path that crosses the DMI-induced gap or band
 *  splitting (e.g., near zone-boundary K-points on a kagome FM).
 *
 *  Parameters and convergence behaviour identical to `_fit_J`.
 *
 *  @param n_sub        number of magnetic sublattices
 *  @param S            spin per site
 *  @param a1, a2       primitive vectors
 *  @param bonds        bond list — J and D's are starting guesses, overwritten with fit
 *  @param n_bonds      number of bonds
 *  @param Kz           uniaxial anisotropy (held fixed)
 *  @param q_obs        n_obs × 2 momenta of observations
 *  @param omega_obs    n_obs observed energies
 *  @param band_obs     n_obs band indices
 *  @param n_obs        number of observations
 *  @param max_iter     maximum iterations
 *  @param tol          convergence tolerance on |Δχ²|
 *  @param chi2_out     final χ² (caller may pass NULL)
 *  @return IRREP_OK on success. */
IRREP_API irrep_status_t irrep_magnon_fit_J_and_DMI(
    int n_sub, double S, const double a1[2], const double a2[2],
    irrep_magnon_bond_t *bonds, int n_bonds, double Kz, const double (*q_obs)[2],
    const double *omega_obs, const int *band_obs, int n_obs, int max_iter, double tol,
    double *chi2_out);

/** @brief Three-magnon dynamical structure factor S^(3)(q, ω) — the
 *         q-resolved 3-magnon continuum, providing spectral weight
 *         for high-energy INS features above the 2-magnon band-top.
 *
 *      S^(3)(q, ω) = (1/N²) Σ_{k1,k2,b1,b2,b3}
 *                       S_⊥_b1(k1) · S_⊥_b2(k2) · S_⊥_b3(q − k1 − k2)
 *                       · (η/π) / ((ω − ω_{b1}(k1) − ω_{b2}(k2) − ω_{b3}(q − k1 − k2))² + η²)
 *
 *  Output is row-major n_q × n_omega; same format as the 2-magnon
 *  family. Cost per q-point is O(Nx² · Ny² · n_sub³ · n_omega) —
 *  significantly more expensive than `_two_magnon_qomega`. Use modest
 *  Nx, Ny (≤ 16) for n_q ≥ 4.
 *
 *  For AFM ground states use `_three_magnon_qomega_general`.
 *
 *  @param L              LSW handle
 *  @param qpath          n_q × 2 momenta
 *  @param n_q            number of q-points
 *  @param Nx, Ny         BZ grid (typical 8-16; quartic in cost)
 *  @param omega_min      lower edge of energy axis
 *  @param omega_max      upper edge (typical 3·ω_max of 1-magnon)
 *  @param n_omega        number of energy bins
 *  @param eta            Lorentzian half-width
 *  @param intensity_out  caller buffer of size n_q × n_omega doubles. */
IRREP_API irrep_status_t irrep_magnon_three_magnon_qomega(const irrep_magnon_lsw_t *L,
                                                            const double (*qpath)[2], int n_q,
                                                            int Nx, int Ny, double omega_min,
                                                            double omega_max, int n_omega,
                                                            double eta,
                                                            double *intensity_out);

/** @brief AFM-aware q-resolved three-magnon S^(3)(q, ω). Same as
 *         `_three_magnon_qomega` but uses the Bogoliubov-aware
 *         structure factor `_structure_factor_general` for the SF
 *         weighting. Required for AFM/ferri ground states where the
 *         pairing controls the spectral weight. */
IRREP_API irrep_status_t irrep_magnon_three_magnon_qomega_general(
    const irrep_magnon_lsw_t *L, const int *sublattice_signs, const double (*qpath)[2],
    int n_q, int Nx, int Ny, double omega_min, double omega_max, int n_omega, double eta,
    double *intensity_out);

/** @brief Three-magnon density of states D^(3)(ω) — spectral support
 *         for the 3-magnon continuum.
 *
 *      D^(3)(ω) = (1/N²) Σ_{k1,k2,b1,b2,b3} δ(ω − ω_{b1}(k1) − ω_{b2}(k2) − ω_{b3}(k3))
 *
 *  with momentum conservation k1 + k2 + k3 = 0 (k3 is determined).
 *  Like the 2-magnon DOS but one more convolution. The 3-magnon
 *  continuum extends to 3·ω_max and provides the spectral support
 *  for high-energy INS continuum features above the 2-magnon
 *  band-top ≈ 2·ω_max.
 *
 *  Sum rule:   ∫ D^(3)(ω) dω = n_sub³
 *
 *  Cost is O(Nx² · Ny² · n_sub³) — quadratic in the BZ grid (one
 *  order more than the 2-magnon DOS), so use modest Nx, Ny (≤ 32)
 *  unless you have time to spare.
 *
 *  @param L              LSW handle
 *  @param Nx, Ny         BZ grid (typical 16-32; quartic in cost)
 *  @param omega_min      lower edge of histogram window
 *  @param omega_max      upper edge (typical 3·ω_max of 1-magnon)
 *  @param n_bins         number of energy bins
 *  @param dos_out        caller buffer of n_bins doubles. */
IRREP_API irrep_status_t irrep_magnon_three_magnon_dos(const irrep_magnon_lsw_t *L, int Nx,
                                                         int Ny, double omega_min,
                                                         double omega_max, int n_bins,
                                                         double *dos_out);

/** @brief User-provided cubic vertex matrix-element-squared
 *         |V_3(k, k1, k-k1; b, b1, b2)|², in units of J².
 *
 *  Called by `irrep_magnon_born_decay_rate` once per (k, k1, b, b1,
 *  b2) tuple in the phase-space sum. Allows callers to plug in any
 *  model-specific vertex (Mourigal-Chernyshev for non-collinear, the
 *  generic Heisenberg cubic for canted FMs, etc.) without the
 *  library committing to a specific Holstein-Primakoff expansion.
 *
 *  Return |V_3|² (NOT V_3 itself) — must be non-negative. The
 *  library multiplies by π and convolves with a Lorentzian over the
 *  energy delta. */
typedef double (*irrep_magnon_cubic_vertex_fn)(int b, int b1, int b2, double kx, double ky,
                                                 double k1x, double k1y, double k2x,
                                                 double k2y, void *user_data);

/** @brief Heisenberg + DMI cubic-vertex Born decay rate for
 *         non-collinear LSW models — full Bogoliubov-mixed Born
 *         self-energy.
 *
 *  Implements the Born self-energy by enumerating all γ†γ†γ
 *  contributions to the matrix element from the four cubic operator
 *  types (a_α·a_β†·a_β, a_α†·a_β†·a_β, a_α†·a_α·a_β, a_α†·a_α·a_β†)
 *  in the local-frame Holstein-Primakoff expansion of the bond
 *  Hamiltonian. Each Type yields 6 γ†γ†γ sub-cases per bond
 *  orientation, giving 48 sub-cases per bond when both orientations
 *  are summed.
 *
 *  Bond cubic coefficient combines Heisenberg + DMI contributions:
 *
 *      M_eff[μν] = J · (R_α^T R_β)_{μν} + D · (R_α[:,μ] × R_β[:,ν])
 *
 *  i.e., the rotation-projected Heisenberg M-matrix plus the
 *  DMI cross-product tensor. Both contribute linearly to the cubic
 *  vertex.
 *
 *  Each sub-case is an ordered product (creator-1, creator-2,
 *  annihilator) of (u, v) Bogoliubov amplitudes from the
 *  noncollinear-LSW BdG diagonalization, with appropriate momentum
 *  conservation and phase factor e^{i k·t}.
 *
 *  Goldstone regularization: the BdG Cholesky uses eps_psd = η to
 *  bound Bogoliubov amplitudes near gap-closing modes. Without this,
 *  the BZ integral on Goldstone-having models (kagome 120° Néel,
 *  triangular AFM) is contaminated by 1/√ω divergences.
 *
 *  TR-symmetric BZ grid: the BZ sum uses a centered grid q ∈ [-1/2,
 *  1/2)·b1 + [-1/2, 1/2)·b2 (no half-shift) so that Γ_b(k) = Γ_b(-k)
 *  to discretization precision rather than half-shift asymmetry.
 *
 *  Validations passing in tests/test_magnon.c:
 *    - Γ ≡ 0 on collinear ground states with U(1) symmetry preserved
 *      (Heisenberg-only, OR Heisenberg + out-of-plane DMI on
 *      collinear-ẑ FM)
 *    - Γ > 0 on kagome 120° Néel (genuine non-collinear cubic vertex)
 *    - Γ > 0 on collinear-ẑ FM with in-plane DMI (U(1) broken)
 *    - Γ_b(k) = Γ_b(-k) (TR symmetry, machine precision)
 *    - Γ invariant under bond reversal ⟨α,β,t,J,D⟩ ↔ ⟨β,α,-t,J,-D⟩
 *      (1e-9 precision)
 *    - Γ_b(k) ≥ 0 per band (unitarity)
 *
 *  ON THE `gap_floor` PARAMETER (READ THIS):
 *
 *  The Born self-energy is logarithmically IR-divergent on models
 *  with gapless modes (Goldstones). |V_3|² scales as |u|⁶ near
 *  gap-closing modes where the Bogoliubov amplitudes diverge as
 *  1/√ω; the on-shell δ-function constraint reduces but does not
 *  cure the singularity, so the BZ integral has a log-divergent
 *  contribution from gap-closing regions on truly-gapless models.
 *
 *  `gap_floor` is the IR cutoff: it adds eps_psd = gap_floor to the
 *  BdG mass matrix in the Cholesky regularization, which bounds the
 *  Bogoliubov amplitudes at √(1/gap_floor). The user picks
 *  gap_floor based on physical context:
 *    - For NATURALLY GAPPED models (single-ion anisotropy, applied
 *      field, finite-T self-consistent gap): gap_floor can be small
 *      (1e-6 .. 1e-10) — the BdG amplitudes are bounded by the
 *      physical gap, not the regularization.
 *    - For GOLDSTONE-HAVING models (kagome 120° Néel, triangular
 *      AFM, helical magnets without applied field): gap_floor IS
 *      the regularization, and Γ depends on it. Common physical
 *      interpretations: dissipation width, finite-T soft-mode gap,
 *      or self-consistent-Born-approximation gap. The integral
 *      converges to a well-defined value at any FIXED gap_floor as
 *      Nx, Ny → ∞ (verified by `test_heisenberg_decay_rate_grid_convergence`).
 *      The gap_floor → 0 limit is genuinely divergent and reflects
 *      the underlying physics, not an algorithmic flaw.
 *
 *  At fixed gap_floor, Γ:
 *    - Converges to a finite value as Nx, Ny → ∞ (grid convergence)
 *    - Has well-defined symmetries (TR, bond-reversal invariance,
 *      U(1) protection on collinear → Γ = 0)
 *    - Is physically interpretable: the regularization is exposed
 *      explicitly so the user knows what they're computing
 *
 *  For matching published linewidth maps (e.g., Mourigal-Chernyshev
 *  PRL 2013 for kagome 120° Néel), pick gap_floor consistent with
 *  the paper's IR treatment (typically a self-consistent or
 *  experimental gap value).
 *
 *  Implements the Born self-energy
 *
 *      Γ_b(k) = π Σ_{q, b1, b2} |V_3(k → q, k-q; b → b1, b2)|²
 *                 · L_η(ω_b(k) - ω_{b1}(q) - ω_{b2}(k-q))
 *
 *  where V_3 is the Heisenberg cubic vertex from local-frame
 *  Holstein-Primakoff expansion. Per bond ⟨α(R), β(R+t)⟩ with
 *  rotation-induced relative orientation M_αβ = R_α^T R_β, the
 *  cubic Hamiltonian's 1→2 piece is
 *
 *      H_3^{1→2}_⟨αβ⟩ = -J √(S/2) · {
 *          (M[XZ] + i M[YZ]) · a_α†(R) a_β†(R+t) a_β(R+t)
 *        + (M[ZX] + i M[ZY]) · a_α†(R) a_α(R) a_β†(R+t) }
 *
 *  After Fourier transform and summing over the bond list (both
 *  orientations), the band-basis matrix element is computed by
 *  sandwiching with the Bogoliubov u-amplitudes from
 *  `_dispersion_noncollinear_full`.
 *
 *  IMPORTANT — Bogoliubov approximation: this implementation uses
 *  ONLY the leading u·u·u* term from the (u, v) Bogoliubov expansion
 *  of the bare-boson cubic vertex. The exact result includes v-mixed
 *  contributions (8 total terms per bond contribution) that are
 *  parametrically suppressed for "FM-like" non-collinear states (low
 *  pairing amplitudes) but become O(1) corrections for strongly-
 *  paired AFMs. For quantitative agreement on canonical references
 *  like the kagome 120° Néel (Mourigal-Chernyshev PRL 2013), the
 *  v contributions need to be added — this is documented and left
 *  for follow-on work.
 *
 *  By construction Γ → 0 on COLLINEAR ground states (M = identity
 *  → no cubic vertex by U(1)), and Γ > 0 on genuine non-collinear
 *  ground states.
 *
 *  @param L           non-collinear LSW handle
 *  @param n_vectors   3·n_sub doubles — local ẑ direction per sublattice
 *  @param kpath       n_k × 2 momenta
 *  @param n_k         number of k-points
 *  @param Nx, Ny      BZ grid for the phase-space sum
 *  @param eta         Lorentzian half-width on the on-shell δ-function
 *                     (energy-axis broadening; independent from gap_floor)
 *  @param gap_floor   IR cutoff added to the BdG mass matrix (eps_psd in
 *                     the Cholesky regularization). Controls the minimum
 *                     gap enforced on the dispersion. For naturally
 *                     gapped models can be set very small (1e-10 to
 *                     1e-6); for Goldstone-having models the caller
 *                     picks based on physics — anisotropy gap, applied
 *                     field gap, dissipation width, finite-temperature
 *                     soft-mode gap, etc. The Born self-energy is
 *                     logarithmically IR-divergent on truly-gapless
 *                     models, so this parameter is intentionally
 *                     required (no default) and the regularization-
 *                     dependence of Γ is the caller's responsibility
 *                     to interpret.
 *  @param gamma_out   n_k × n_sub doubles — Γ_b(k) in J units */
IRREP_API irrep_status_t irrep_magnon_heisenberg_decay_rate(const irrep_magnon_lsw_t *L,
                                                              const double *n_vectors,
                                                              const double (*kpath)[2],
                                                              int n_k, int Nx, int Ny,
                                                              double eta, double gap_floor,
                                                              double *gamma_out);

/** @brief Born-approximation magnon decay rate Γ_b(k) from a user-
 *         provided cubic vertex, via Fermi golden rule:
 *
 *      Γ_b(k) = π Σ_{q, b1, b2} |V_3(k, q, k-q; b, b1, b2)|²
 *                 · (η/π) / ((ω_b(k) - ω_{b1}(q) - ω_{b2}(k-q))² + η²)
 *
 *  This is the matrix-element-aware companion to
 *  `_kinematic_damping`: where the kinematic version returns the
 *  phase-space upper bound (|V_3|² ≡ 1), this returns the actual
 *  linewidth for whatever vertex the caller supplies.
 *
 *  Callers building |V_3|² for a specific model:
 *    - Collinear FM Heisenberg: vertex vanishes by U(1) → return 0.
 *    - Collinear AFM Heisenberg: vertex vanishes by U(1) → return 0.
 *    - Canted FM (field-tilted): vertex ∝ sin(canting)·J · γ-function.
 *    - Kagome 120° Néel: see Mourigal-Chernyshev PRL 2013 for the
 *      explicit form including DMI corrections.
 *
 *  Output units are J (matches the dispersion energy scale).
 *
 *  @param L              LSW handle (used for ω_b(k) and BZ geometry)
 *  @param kpath          n_k × 2 momenta
 *  @param n_k            number of k-points
 *  @param Nx, Ny         BZ grid for the phase-space sum
 *  @param eta            Lorentzian half-width on the energy delta
 *  @param vertex_fn      user callback returning |V_3|² (units J²)
 *  @param vertex_data    opaque pointer passed to vertex_fn
 *  @param gamma_out      n_k × n_sub doubles. */
IRREP_API irrep_status_t irrep_magnon_born_decay_rate(const irrep_magnon_lsw_t *L,
                                                       const double (*kpath)[2], int n_k,
                                                       int Nx, int Ny, double eta,
                                                       irrep_magnon_cubic_vertex_fn vertex_fn,
                                                       void *vertex_data, double *gamma_out);

/** @brief Kinematic upper bound on magnon decay rate Γ_kin(k, b) —
 *         the 2-magnon phase-space estimate of the linewidth.
 *
 *      Γ_kin_b(k) = π · D⁽²⁾(k, ω_b(k))
 *
 *  where D⁽²⁾(k, ω) is the q-resolved 2-magnon density of states
 *  evaluated at the 1-magnon dispersion energy. This is the phase-
 *  space available for 1→2 magnon decay via cubic anharmonicity (the
 *  matrix-element-stripped Born approximation). It vanishes when the
 *  decay channel is kinematically forbidden (k below the 2-magnon
 *  band-bottom) and is finite when allowed.
 *
 *  Caller multiplies by their estimate of |V_3(k)|² (the cubic-vertex
 *  squared amplitude — material- and ground-state-dependent) for a
 *  full Born linewidth. For order-of-magnitude estimates the kinematic
 *  factor alone is what determines whether decays are *possible*.
 *
 *  Concrete example: in a collinear FM Heisenberg model the cubic
 *  vertex vanishes by U(1) symmetry → Γ_real = 0 even when Γ_kin > 0.
 *  In a non-collinear AFM (kagome 120° Néel) the cubic vertex is
 *  O(sin(canting)) and the full damping is Γ_kin · (matrix element).
 *
 *  Output units: same as ω (matches the dispersion energy scale).
 *
 *  @param L           LSW handle
 *  @param kpath       n_k × 2 momenta to evaluate
 *  @param n_k         number of k-points
 *  @param Nx, Ny      BZ grid for the 2-magnon convolution
 *  @param eta         Lorentzian half-width on the energy axis
 *  @param gamma_out   caller buffer of size n_k × n_sub doubles. */
IRREP_API irrep_status_t irrep_magnon_kinematic_damping(const irrep_magnon_lsw_t *L,
                                                         const double (*kpath)[2], int n_k,
                                                         int Nx, int Ny, double eta,
                                                         double *gamma_out);

/** @brief Two-magnon density of states D⁽²⁾(ω) — spectral support
 *         for the 2-magnon continuum.
 *
 *      D⁽²⁾(ω) = (1/N²) Σ_{k₁, k₂, b₁, b₂} δ(ω − ω_{b₁}(k₁) − ω_{b₂}(k₂))
 *
 *  The 2-magnon DOS gives the energy range where 2-magnon scattering
 *  has spectral support — a real observable in INS spectra as the
 *  "continuum shoulder" above the sharp 1-magnon peak. For example,
 *  in cuprate INS data above the 2J·S band-top, intensity persists
 *  due to the 2-magnon continuum extending up to 4J·S = 2·ω_max.
 *
 *  Goes beyond simple 1-magnon LSW: this is a TWO-PARTICLE observable
 *  that requires summing over pairs (k₁, k₂) of magnons. Even at the
 *  tree level (no 1/S corrections), 2-magnon DOS reveals features
 *  invisible to 1-magnon DOS:
 *    - Total bandwidth doubled (2·ω_max)
 *    - Peak shifted from middle of 1-magnon DOS
 *    - Convolution suppresses van-Hove singularities
 *
 *  Sum rule:   ∫ D⁽²⁾(ω) dω = n_sub²
 *
 *  Implementation: histogram the energy sums ω(k₁) + ω(k₂) over all
 *  Nx² × Ny² pairs of BZ points.
 *
 *  @param L         LSW handle
 *  @param Nx, Ny    BZ grid (typical 32-64; quadratic in cost)
 *  @param omega_min lower bound of energy histogram
 *  @param omega_max upper bound (typical 2·ω_max of 1-magnon)
 *  @param n_bins    number of energy bins
 *  @param dos_out   caller buffer of n_bins doubles. */
IRREP_API irrep_status_t irrep_magnon_two_magnon_dos(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                                      double omega_min, double omega_max,
                                                      int n_bins, double *dos_out);

/** @brief Hartree-Fock-renormalised magnon dispersion at finite T —
 *         the leading 1/S correction beyond LSW.
 *
 *  LSW (linear spin-wave theory) keeps only bilinear-in-bosons terms
 *  in the Holstein-Primakoff expansion. The next term is QUARTIC,
 *  encoding magnon-magnon scattering. Hartree-Fock decoupling of
 *  these quartic terms gives a finite-T renormalised dispersion:
 *
 *      ω_HF(k, T) = ω_LSW(k) · Z(T)
 *
 *  where Z(T) = 1 − ⟨n⟩(T)/S is the standard Bloch-Dyson factor,
 *  with ⟨n⟩(T) = (1/N_BZ) Σ_b ∫ n_BE(ω_LSW/T) d²k the average
 *  thermal magnon population per cell.
 *
 *  Limits:
 *    - T → 0:           Z → 1 (no correction; recover LSW)
 *    - T → T_c (= O(JS)): Z → 0 (bands collapse, LSW breakdown signal)
 *
 *  At T < ~T_c/3, Z is within a few percent of 1 and LSW is
 *  quantitatively reliable. Near T_c, Z becomes O(1) and full
 *  beyond-LSW (Schwinger-boson MF, classical MC) is needed.
 *
 *  This function returns the scalar Z(T); the caller multiplies any
 *  LSW dispersion to get the HF-renormalised one. Useful as both a
 *  beyond-LSW correction AND a *diagnostic* for LSW validity:
 *  Z >> 0.7 means LSW is fine; Z < 0.5 means breakdown.
 *
 *  References:
 *    - Bloch (1930): T^{3/2} magnetisation law
 *    - Dyson (1956): leading 1/S spin-wave theory beyond LSW
 *    - Holstein-Primakoff (1940): the bosonisation
 *
 *  @param L      LSW handle
 *  @param T      temperature (same units as ω)
 *  @param Nx, Ny BZ grid (typical 32-128)
 *  @return       Z(T) ∈ [0, 1] (or NaN on error). */
IRREP_API double irrep_magnon_hartree_renormalisation(const irrep_magnon_lsw_t *L, double T,
                                                       int Nx, int Ny);

/** @brief Compute the magnon group velocity v_g(k) = ∇_k ω_b(k) for
 *         each band at a single k point via central finite difference.
 *
 *  The group velocity tells you how fast a magnon wave-packet
 *  propagates and is the natural input for spin-current transport
 *  calculations:
 *
 *      j_s = (1/V) Σ_b ∫ d²k/(2π)² · v_g_b(k) · n_B(ω_b/T)
 *
 *  is the magnon spin-current density. Regions of the BZ where
 *  v_g is large dominate spin transport; regions with v_g → 0 (band
 *  bottoms / tops, flat bands) are spectator modes. Linearly
 *  dispersing magnons (Dirac points, AFM Goldstone) have constant
 *  |v_g| ~ c_s (spin-wave velocity).
 *
 *  Central-difference estimator:
 *
 *      v_x_b(k) = [ω_b(k_x + h, k_y) − ω_b(k_x − h, k_y)] / (2h)
 *      v_y_b(k) = [ω_b(k_x, k_y + h) − ω_b(k_x, k_y − h)] / (2h)
 *
 *  Step `h` (typical 1e-3 - 1e-2): smaller h gives sharper local
 *  resolution but more numerical noise; larger h smooths over band
 *  curvature. Near band crossings, the eigenvector ordering can
 *  flip and the central difference picks up a discontinuity — the
 *  result there is unreliable; the user should be aware of band-
 *  crossing locations from the dispersion plot.
 *
 *  @param L         LSW handle
 *  @param kx, ky    central momentum
 *  @param h         finite-difference step (typical 1e-3 - 1e-2)
 *  @param vx_out    caller buffer of size n_sub doubles
 *  @param vy_out    caller buffer of size n_sub doubles */
IRREP_API irrep_status_t irrep_magnon_group_velocity(const irrep_magnon_lsw_t *L, double kx,
                                                      double ky, double h, double *vx_out,
                                                      double *vy_out);

/** @brief Magnon contribution to the per-cell internal energy U(T)
 *         in natural units (k_B = ℏ = 1).
 *
 *      U(T) = Σ_b ∫_BZ d²k/(2π)² · ω_b(k) · n_B(ω_b/T)
 *
 *  This is the *correction* to the classical-ground-state energy
 *  E_cl from thermally-populated magnons; total internal energy of
 *  a magnetic system is E_cl + U(T).
 *
 *  Limits:
 *
 *    - **T → 0**: U → 0 like Σ_b ∫ ω·e^{−ω/T} (vanishes
 *      exponentially for gapped systems; in 3D FM with quadratic
 *      Goldstone, U ∝ T^{5/2} via the Bloch DOS).
 *
 *    - **T ≫ bandwidth**: equipartition gives U → n_sub · T (each
 *      band per cell contributes 1·k_B·T in classical limit). This
 *      is the same regime where C_V → n_sub k_B; LSW remains
 *      reliable while T < bandwidth.
 *
 *  Thermodynamic consistency:
 *
 *      U(T) = F(T) + T·S_th(T)        with S_th = −∂F/∂T
 *      C_V(T) = ∂U/∂T                  (Maxwell relation)
 *
 *  Direct application: the magnon contribution to the Joule-Thomson
 *  cooling coefficient and to magnetic-cooling device specifications
 *  (energy storage capacity per unit cell vs operating temperature).
 *
 *  Same Goldstone-mode (ω < 10⁻¹⁰) and BE-overflow (x > 700) regul-
 *  arisations as the rest of the thermodynamic API.
 *
 *  @param L      LSW handle
 *  @param T      temperature (same units as ω)
 *  @param Nx, Ny BZ grid (typical 50-200)
 *  @return       U(T) per unit cell, or NaN on error. */
IRREP_API double irrep_magnon_internal_energy(const irrep_magnon_lsw_t *L, double T, int Nx,
                                               int Ny);

/** @brief Magnon contribution to the Helmholtz free energy F(T) per
 *         unit cell, in natural units (k_B = ℏ = 1). For a Bose gas
 *         of non-interacting magnons,
 *
 *      F(T) = T · Σ_b ∫_BZ d²k/(2π)² · ln(1 − e^{−ω_b/T})
 *
 *  (negative for any T > 0). This is the *correction* to the
 *  classical-ground-state energy from thermally-populated magnons;
 *  the total Helmholtz free energy of a magnetic system is E_cl +
 *  F(T) (caller-supplied E_cl, since LSW alone does not fix the
 *  classical ground state).
 *
 *  Limits:
 *
 *    - **T → 0**: F → 0 like −T·Σ e^{−ω_b/T} (vanishes exponentially
 *      for gapped systems; in 3D FM with quadratic Goldstone, F ∝
 *      −T^{5/2} via the Bloch-law DOS).
 *
 *    - **T ≫ bandwidth**: F ∝ T · ln(T) — unphysical regime, like
 *      χ(T) and M(T). HP bosonisation breaks down well before this.
 *
 *  Connections to other thermodynamic functions:
 *
 *      U(T) = F(T) + T·S_th(T)        (internal energy)
 *      S_th(T) = −∂F/∂T               (entropy)
 *      C_V(T) = T · ∂²F/∂T² · (−1)    (= ∂U/∂T; computed independently
 *                                        by `_specific_heat`).
 *
 *  Free-energy *differences* between competing magnetic phases (FM
 *  vs Néel vs canted, etc.) drop out the classical-ground-state E_cl
 *  if both phases have the same total spin per cell, and F(T) alone
 *  is enough to pick the energetically favoured phase at finite T.
 *
 *  @param L      LSW handle
 *  @param T      temperature (same units as ω)
 *  @param Nx, Ny BZ grid (typical 50-200)
 *  @return       F(T) per unit cell, or NaN on error. */
IRREP_API double irrep_magnon_free_energy(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny);

/** @brief Magnon contribution to the longitudinal magnetic
 *         susceptibility χ(T) per unit cell, in natural units.
 *
 *  For a Bose gas of magnons in a Zeeman-coupled FM, the
 *  fluctuation-dissipation relation gives
 *
 *      χ(T) = (1/T) · Σ_b ∫_BZ d²k/(2π)² · n_B(ω_b/T) · (1 + n_B(ω_b/T))
 *
 *  in units of (g·μ_B)²/k_B. Limits:
 *
 *    - **T → 0, gapped systems** (uniaxial K_z > 0 or DMI gap):
 *      χ → 0 exponentially (no thermal magnons).
 *    - **T → 0, gapless systems** (no anisotropy): χ diverges as
 *      Σ_k 1/(βω_k)² ∝ T^(d/z − 2) for z = 2 quadratic Goldstone
 *      (FM): 2D log-divergent (Mermin-Wagner echo); 3D: T^{−1/2}
 *      finite divergence.
 *    - **T ≫ bandwidth**: LSW formula gives χ ∝ T (Bose enhancement
 *      unbounded). This is *unphysical*: HP bosonisation breaks down
 *      well before T ~ bandwidth, and the result should not be
 *      trusted at high T. Use χ only at T ≲ bandwidth.
 *
 *  Within the LSW regime, χ(T) rises steeply once T crosses the
 *  spin-wave gap and continues to grow. The susceptibility peak
 *  identifying T_c (Curie or Néel) requires beyond-LSW physics
 *  (Schwinger-boson MF, RPA, classical Monte-Carlo).
 *
 *  Same Goldstone (ω < 10⁻¹⁰) and BE-overflow (x > 700) regularis-
 *  ations as `_specific_heat`.
 *
 *  @param L      LSW handle
 *  @param T      temperature (same units as ω)
 *  @param Nx, Ny BZ grid (typical 50-200)
 *  @return       χ(T) per unit cell, or NaN on error. */
IRREP_API double irrep_magnon_susceptibility(const irrep_magnon_lsw_t *L, double T, int Nx,
                                              int Ny);

/** @brief Finite-temperature sublattice-averaged magnetisation M(T)
 *         per spin, in units of (g·μ_B). For a collinear-FM ground
 *         state with spin S per site, the LSW result at T > 0 is
 *
 *      M(T) = S − (1/n_sub) · Σ_b ∫_BZ d²k/(2π)² n_B(ω_b(k)/T)
 *
 *  where the BZ integral counts the thermally-excited magnon
 *  population per cell. The leading low-T behaviour is the famous
 *  *Bloch law*:
 *
 *      M(T) − M(0) ∝ −T^{d/2}                     (FM, quadratic Goldstone)
 *
 *  with d the spatial dimension. d = 3 (cubic FM): the canonical
 *  T^{3/2} signature. d = 2: log-divergent at any T > 0 (Mermin-
 *  Wagner — no long-range FM order in 2D Heisenberg without
 *  anisotropy); the function still returns a numerical value for
 *  *gapped* 2D FM models (with K_z > 0 or DMI gap), where the
 *  integral converges.
 *
 *  LSW assumes small spin fluctuations; the result is reliable when
 *  S − M(T) ≪ S. At higher temperatures the user should switch to
 *  Schwinger-boson or Monte-Carlo treatments.
 *
 *  Goldstone-mode regularisation (ω < 10⁻¹⁰) and BE overflow guard
 *  (x = ω/T > 700) are identical to `_specific_heat`.
 *
 *  @param L      LSW handle
 *  @param T      temperature (same units as ω)
 *  @param Nx, Ny BZ grid (typical 50–200)
 *  @return       M(T) per spin (in (g·μ_B) units), or NaN on error. */
IRREP_API double irrep_magnon_magnetization(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny);

/** @brief Magnon contribution to the magnetic specific heat C_V(T)
 *         per unit cell, in units of k_B.
 *
 *  For a Bose gas of non-interacting magnons,
 *
 *      C_V(T) = Σ_b ∫_BZ d²k/(2π)² · (ω_b/T)² · n_B(ω_b/T)·(1 + n_B(ω_b/T))
 *
 *  with n_B(x) = 1/(eˣ − 1). Limits:
 *
 *    - T → 0: C_V vanishes exponentially (BE-suppressed). For
 *      gapless modes, the leading power follows the dispersion
 *      shape: T (2D quadratic FM Goldstone) ∝ T (∫ const · T·dT),
 *      T^{3/2} (3D quadratic FM, Bloch law), T^d (d-dim linear AFM
 *      Goldstone).
 *
 *    - T → ∞: C_V → n_sub · k_B (equipartition: each band per cell
 *      contributes 1 k_B of heat capacity).
 *
 *  Sampling on Nx × Ny grid; the Goldstone-mode contribution at
 *  k = 0 is regularised by skipping ω < 10⁻¹⁰ to avoid the BE
 *  singularity (the missing measure-zero contribution is well
 *  below grid noise).
 *
 *  @param L      LSW handle
 *  @param T      temperature (same units as ω)
 *  @param Nx, Ny BZ grid (typical 50–200)
 *  @return       C_V(T) per unit cell (dimensionless in k_B = ℏ = 1
 *                units), or NaN on error. */
IRREP_API double irrep_magnon_specific_heat(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny);

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

/** @brief Powder-averaged 1-magnon spectrum S_powder(ω) — the BZ-
 *         integrated, structure-factor-weighted DOS.
 *
 *      S_powder(ω) = (1/N_BZ) Σ_{k, b} S_⊥_b(k) · δ(ω - ω_b(k))
 *
 *  This is the angular-averaged inelastic-neutron intensity measured
 *  on POWDER samples — distinct from the plain magnon DOS in that
 *  bands with low spectral weight (dark on transverse INS) are
 *  suppressed. For a kagome FM, the upper bands are dark at Γ but
 *  carry weight elsewhere in the BZ; the powder spectrum picks that
 *  up correctly.
 *
 *  Sum rule:  ∫ S_powder(ω) dω = 2S · n_sub.
 *
 *  Use `_powder_spectrum_general` for AFM/ferri ground states; the
 *  Bogoliubov-Colpa structure factor is required for correct
 *  spectral weight there.
 *
 *  @param L         LSW handle
 *  @param Nx, Ny    BZ integration grid
 *  @param omega_min, omega_max  energy histogram window
 *  @param n_bins    number of energy bins
 *  @param S_out     caller buffer of n_bins doubles. */
IRREP_API irrep_status_t irrep_magnon_powder_spectrum(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                                       double omega_min, double omega_max,
                                                       int n_bins, double *S_out);

/** @brief AFM-aware powder-averaged 1-magnon spectrum — same as
 *         `_powder_spectrum` but with sublattice signs and the
 *         Bogoliubov structure factor. */
IRREP_API irrep_status_t irrep_magnon_powder_spectrum_general(const irrep_magnon_lsw_t *L,
                                                               const int *sublattice_signs,
                                                               int Nx, int Ny, double omega_min,
                                                               double omega_max, int n_bins,
                                                               double *S_out);

/** @brief Band-resolved transverse structure factor for a *general*
 *         (collinear FM/AFM/ferri) ground state — uses the Bogoliubov-
 *         Colpa machinery instead of the FM-only Hermitian path.
 *
 *  For AFM/ferri ground states, the transverse INS observable involves
 *  both u (particle) and v (hole) Bogoliubov amplitudes:
 *
 *      S_⊥_b(q) = 2S · |Σ_α (u_α^b(q) + v_α^b(q))|²
 *
 *  where (u, v) come from the Bogoliubov-Colpa transformation of the
 *  generalised BdG matrix M(k). For collinear FM (all σ_α = +1),
 *  v = 0 and this reduces to `_structure_factor`.
 *
 *  Sum rule:  Σ_b S_⊥_b(q) = 2S · n_sub  for every q (just like the
 *  FM-track function).
 *
 *  Useful for INS spectrum prediction on AFM materials where the
 *  Bogoliubov pairing is essential to the spectral weight (e.g.,
 *  cuprate La₂CuO₄ structure factor below the band-top).
 *
 *  @param L                LSW handle
 *  @param sublattice_signs σ_α ∈ {+1, -1} per sublattice
 *  @param qx, qy           momentum (cartesian)
 *  @param omega_out        n_sub doubles — band energies sorted asc
 *  @param S_perp_out       n_sub doubles — band-resolved S_⊥_b(q) */
IRREP_API irrep_status_t irrep_magnon_structure_factor_general(const irrep_magnon_lsw_t *L,
                                                                 const int *sublattice_signs,
                                                                 double qx, double qy,
                                                                 double *omega_out,
                                                                 double *S_perp_out);

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

/** @brief Compute the band-resolved transverse structure factor with
 *         intra-cell INS form factor — proper q-dependent phases for
 *         each sublattice's position within the unit cell.
 *
 *      S_⊥_b(q) = 2S · |Σ_α e^{i q·r_α} u_α^b(q)|²
 *
 *  This is the *full* INS-measurable spectral weight on a band: the
 *  existing `_structure_factor` returns the Wannier-centred form
 *  |Σ u|² (sublattice positions absorbed into the unit-cell origin),
 *  while this function returns what neutron experiments actually see
 *  at finite q — the intra-cell phase factor e^{iq·r_α} can swap
 *  bright vs. dark bands at zone-boundary momenta.
 *
 *  At q = Γ (q = 0) the phase factor is unity for all α and this
 *  reduces to `_structure_factor`. At finite q the difference is
 *  band- and lattice-specific (kagome: pronounced; bipartite square:
 *  absorbed into the existing magnetic-BZ folding).
 *
 *  Sum rule preserved: integrated over the BZ, Σ_b ⟨S_⊥_b(q)⟩ =
 *  2S · n_sub (independent of positions).
 *
 *  @param L              LSW handle
 *  @param positions      n_sub × 2 array of cartesian sublattice
 *                        positions r_α inside the unit cell
 *  @param qx, qy         momentum (cartesian)
 *  @param omega_out      n_sub doubles — band energies sorted ascending
 *  @param S_perp_out     n_sub doubles — band-resolved form-factor SF */
IRREP_API irrep_status_t irrep_magnon_structure_factor_with_form_factor(
    const irrep_magnon_lsw_t *L, const double (*positions)[2], double qx, double qy,
    double *omega_out, double *S_perp_out);

/** @brief Polarisation-resolved 1-magnon structure factors —
 *         S^xx_b(q), S^yy_b(q), S^zz_b(q) per band, separately.
 *
 *  Polarised INS experiments measure these channels directly via
 *  spin-flip / non-spin-flip discrimination. For a collinear FM
 *  along ẑ at T=0:
 *
 *      S^xx_b(q) = S^yy_b(q) = ½ · 2S · |Σ_α u_α^b(q)|² = ½ S_⊥_b
 *      S^zz_b(q) = 0   (longitudinal — no 1-magnon weight)
 *
 *  For non-collinear ground states (specified by sublattice direction
 *  vectors `n_hat`) the rotation matrices R_α mix the channels:
 *
 *      S^aa_b(q) = (S/2) · |Σ_α c_α^a · u_α^b(q)|²
 *      with c_α^a = R_α[a, X] - i R_α[a, Y]
 *
 *  All three channels can carry 1-magnon weight in the non-collinear
 *  case — including the longitudinal S^zz_b(q) which vanishes only
 *  in the collinear limit.
 *
 *  Pass `n_hat = NULL` for the collinear-FM-along-ẑ default.
 *
 *  @param L           LSW handle
 *  @param n_hat       n_sub × 3 ground-state directions, or NULL
 *  @param qx, qy      momentum
 *  @param omega_out   n_sub doubles — band energies (asc)
 *  @param Sxx_out     n_sub doubles — S^xx per band
 *  @param Syy_out     n_sub doubles — S^yy per band
 *  @param Szz_out     n_sub doubles — S^zz per band */
IRREP_API irrep_status_t irrep_magnon_polarized_structure_factor(const irrep_magnon_lsw_t *L,
                                                                   const double (*n_hat)[3],
                                                                   double qx, double qy,
                                                                   double *omega_out,
                                                                   double *Sxx_out,
                                                                   double *Syy_out,
                                                                   double *Szz_out);

/** @brief AFM-aware version of `_structure_factor_with_form_factor`.
 *         Uses the Bogoliubov-Colpa amplitudes (u, v) and applies the
 *         intra-cell form factor exp(i q·r_α) before summing. */
IRREP_API irrep_status_t irrep_magnon_structure_factor_general_with_form_factor(
    const irrep_magnon_lsw_t *L, const int *sublattice_signs, const double (*positions)[2],
    double qx, double qy, double *omega_out, double *S_perp_out);

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

/** @brief Non-collinear LSW magnon dispersion for an arbitrary
 *         classical-ground-state direction per sublattice.
 *
 *  Generalises `_dispersion_general` (which only handles σ_α ∈
 *  {+1, -1} collinear) to ANY collinear, non-collinear, or canted
 *  ground state — required for skyrmions, helimagnets, 120° Néel
 *  on triangular/kagome AFMs, and chiral magnets like MnSi /
 *  Cu₂OSeO₃ / FeGe in their helical phases.
 *
 *  Each sublattice α carries a unit-vector ground state n_α (3D
 *  real, |n_α| = 1). HP bosonisation is performed in a LOCAL FRAME
 *  per sublattice with the +ẑ axis rotated to n_α; bond Hamiltonian
 *  terms are expressed in terms of M_ij = R_i^T R_j (the relative
 *  rotation between local frames) and the DMI matrix N_ij computed
 *  similarly. The bilinear-in-bosons Hamiltonian then has FOUR
 *  channels per bond: hopping (a^†_i a_j + h.c.), pairing (a_i a_j +
 *  a^†_i a^†_j), DMI hopping (anti-Hermitian off-diagonal), and
 *  anti-symmetric pairing — all derived from the entries of c_{de}
 *  = J·M_{de} + DMI contribution.
 *
 *  Linear-in-bosons terms (from S̃^x S̃^z and S̃^y S̃^z mixed
 *  channels) are assumed to cancel by the equilibrium condition of
 *  the supplied ground state. Caller responsibility: provide a
 *  physically valid n_α set (the function does not check stability).
 *
 *  Reduces to:
 *    - `_dispersion`         when all n_α = +ẑ (FM)
 *    - `_dispersion_general` (σ_α=±1) when n_α = ±ẑ
 *
 *  @param L                LSW handle
 *  @param n_vectors        n_sub × 3 array of unit ground-state
 *                          vectors (row-major: n_vectors[3*α + i] is
 *                          the i-th component of n_α)
 *  @param kx, ky           momentum (cartesian)
 *  @param omega_out        caller buffer of size n_sub doubles —
 *                          magnon energies sorted ascending */
IRREP_API irrep_status_t irrep_magnon_dispersion_noncollinear(const irrep_magnon_lsw_t *L,
                                                                const double *n_vectors, double kx,
                                                                double ky, double *omega_out);

/** @brief Non-collinear LSW dispersion AND Bogoliubov amplitudes.
 *         Same as `_dispersion_noncollinear` but additionally returns
 *         the (u, v) BdG eigenvectors per band — required infrastructure
 *         for cubic-vertex calculations and other beyond-LSW work that
 *         needs the full Bogoliubov transformation explicitly.
 *
 *  Output `bogoliubov_uv` is row-major n_sub × 2·n_sub complex.
 *  Row b is the BdG eigenvector for band b: the first n_sub entries
 *  are the u-amplitudes (one per sublattice), the next n_sub are
 *  the v-amplitudes. Magnon-creation operator γ_b†(k) is then
 *
 *      γ_b†(k) = Σ_α [u_α^b(k) · a_α†(k) + v_α^b(k) · a_α(-k)]
 *
 *  Bands are sorted ascending by ω.
 *
 *  @param L              non-collinear LSW handle (built via the
 *                        `_dispersion_noncollinear` route)
 *  @param n_vectors      flat 3·n_sub doubles — local ẑ direction
 *                        per sublattice (cartesian)
 *  @param kx, ky         momentum
 *  @param omega_out      n_sub doubles — band energies (asc)
 *  @param bogoliubov_uv  n_sub × 2·n_sub complex — BdG eigenvectors */
IRREP_API irrep_status_t irrep_magnon_dispersion_noncollinear_full(
    const irrep_magnon_lsw_t *L, const double *n_vectors, double kx, double ky,
    double *omega_out, double _Complex *bogoliubov_uv);

/** @brief 3D non-collinear LSW dispersion. Same as
 *         `irrep_magnon_dispersion_noncollinear` but consumes a third
 *         primitive vector a₃ and per-bond `delta_z`.
 *
 *  Closes the capability matrix:
 *    - 2D FM:     `_dispersion`
 *    - 2D AFM:    `_dispersion_general`
 *    - 2D non-collinear: `_dispersion_noncollinear`
 *    - 3D FM:     `_dispersion_3d`
 *    - 3D AFM:    `_dispersion_general_3d`
 *    - 3D non-collinear: `_dispersion_noncollinear_3d` (this function)
 *
 *  Real applications: B20 cubic chiral magnets MnSi/Cu₂OSeO₃/FeGe in
 *  3D-helical-spiral phase, 3D pyrochlore AFM with non-collinear all-
 *  in / all-out ordering, layered helimagnet stacks.
 *
 *  @param L                LSW handle
 *  @param n_vectors        n_sub × 3 array of unit ground-state vectors
 *  @param a3               third primitive vector (cartesian 3D)
 *  @param kx, ky, kz       momentum (cartesian 3D)
 *  @param omega_out        caller buffer of size n_sub doubles —
 *                          magnon energies sorted ascending */
IRREP_API irrep_status_t irrep_magnon_dispersion_noncollinear_3d(const irrep_magnon_lsw_t *L,
                                                                   const double *n_vectors,
                                                                   const double a3[3], double kx,
                                                                   double ky, double kz,
                                                                   double *omega_out);

/** @brief 3D Bogoliubov-Colpa AFM dispersion. Same as
 *         `irrep_magnon_dispersion_general` but consumes a third
 *         primitive vector a₃ and per-bond `delta_z` to handle 3D
 *         lattices: 3D-cubic AFM (NaCl-like Néel), 3D-AFM stacks of
 *         2D layers, pyrochlore AFM with collinear ground state, etc.
 *
 *  Closes the (FM, AFM) × (2D, 3D) capability matrix:
 *    - 2D FM:  irrep_magnon_dispersion          (Hermitian)
 *    - 2D AFM: irrep_magnon_dispersion_general  (Bogoliubov-Colpa)
 *    - 3D FM:  irrep_magnon_dispersion_3d       (Hermitian)
 *    - 3D AFM: irrep_magnon_dispersion_general_3d (this function)
 *
 *  @param L                LSW handle
 *  @param sublattice_signs σ_α ∈ {+1, -1} per sublattice
 *  @param a3               third primitive vector (cartesian 3D)
 *  @param kx, ky, kz       momentum (cartesian 3D)
 *  @param omega_out        caller buffer of size n_sub doubles —
 *                          magnon energies sorted ascending */
IRREP_API irrep_status_t irrep_magnon_dispersion_general_3d(const irrep_magnon_lsw_t *L,
                                                             const int *sublattice_signs,
                                                             const double a3[3], double kx,
                                                             double ky, double kz,
                                                             double *omega_out);

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

/** @brief Quantum zero-point reduction of the sublattice magnetisation
 *         on a generic collinear (FM/AFM/ferrimagnetic) ground state.
 *
 *  For an antiferromagnet, the classical Néel state is *not* an
 *  eigenstate of the Heisenberg Hamiltonian; quantum fluctuations
 *  reduce the staggered magnetisation per site at T = 0 by the
 *  Anderson (1952) correction:
 *
 *      ⟨n_α⟩_GS = Σ_b ∫_BZ d²k/(2π)² · |v_α^b(k)|²
 *
 *  where v_α^b are the "anomalous" components of the Bogoliubov-
 *  Colpa transformation T = K⁻¹ U. The actual sublattice magneti-
 *  sation at T = 0 is M_α(T=0) = S − ⟨n_α⟩_GS, falsifiable by
 *  μSR Knight-shift / Mössbauer-spectroscopy / NMR data.
 *
 *  Textbook result: spin-½ square AFM gives ⟨n_α⟩_GS ≈ 0.1966
 *  (Anderson 1952), so M(T=0) = 0.5 − 0.197 = 0.303 — observed in
 *  La₂CuO₄ to within a few-percent ring-exchange correction.
 *
 *  Pure ferromagnets (all `sublattice_signs[α] = +1`) have no
 *  anomalous pairing and ⟨n_α⟩_GS = 0 identically — verified as
 *  a smoke-test in the test suite.
 *
 *  Internally re-runs the Cholesky-Colpa diagonalisation from
 *  `_dispersion_general` and computes the v-coefficients via
 *  back-substitution K · T = ψ_+.
 *
 *  @param L                LSW handle
 *  @param sublattice_signs σ_α ∈ {+1, -1} per sublattice
 *  @param Nx, Ny           BZ grid (typical 32-128)
 *  @param delta_m_out      caller buffer of size n_sub doubles —
 *                          per-sublattice ⟨n_α⟩_GS quantum reduction. */
IRREP_API irrep_status_t irrep_magnon_afm_zero_point(const irrep_magnon_lsw_t *L,
                                                      const int *sublattice_signs, int Nx, int Ny,
                                                      double *delta_m_out);

/** @brief Compute the 2D topological charge (skyrmion number) of a
 *         spin-vector field on a square lattice via the Berg-Lüscher
 *         signed-solid-angle method.
 *
 *  The continuum topological charge of a 2D unit-vector field n̂(r)
 *  is
 *
 *      Q = (1/4π) ∫ d²r · n̂ · (∂_x n̂ × ∂_y n̂)
 *
 *  classifying the field by its π₂(S²) homotopy class (Q ∈ ℤ for
 *  smooth boundary-condition-respecting textures). On a discrete
 *  lattice the Berg-Lüscher (1981) discretisation evaluates the
 *  integral as a sum of signed solid angles, one per triangular
 *  half-plaquette:
 *
 *      Ω(a, b, c) = 2 · atan2( a · (b × c),
 *                              1 + a·b + b·c + c·a )
 *
 *      Q = (1/4π) Σ_plaquettes [Ω(n_00, n_10, n_11) + Ω(n_00, n_11, n_01)]
 *
 *  This formula gives ℤ-valued results for smoothly-varying spin
 *  fields up to lattice-discretisation errors on rapidly-varying
 *  textures.
 *
 *  Use cases:
 *    - Skyrmion-number measurement on real-space spin textures
 *    - Verification of topological-charge conservation in dynamical
 *      simulations
 *    - Identification of monopole / skyrmion / hopfion features in
 *      lattice-spin output from MC, MD, or LSW relaxation
 *    - Numerical tests of homotopy classification (e.g., π₂(S²)
 *      degree-counting on interpolating Hamiltonians)
 *
 *  PERIODIC BOUNDARY CONDITIONS are assumed: the spin field is
 *  treated as wrapping around in both directions.
 *
 *  @param n_field    flat 3·Nx·Ny doubles — unit spin vectors per
 *                    site in row-major order (site index = iy·Nx + ix);
 *                    each (n_x, n_y, n_z) triple should be unit-length
 *                    (caller's responsibility)
 *  @param Nx, Ny     lattice dimensions
 *  @param charge_out caller-allocated double — receives Q
 *  @return IRREP_OK on success */
IRREP_API irrep_status_t irrep_magnon_topological_charge_2d(const double *n_field, int Nx,
                                                              int Ny, double *charge_out);

/** @brief Compute the discrete Whitehead-integral Hopf charge H of a
 *         3D unit-spin field on a periodic cubic lattice:
 *
 *             H = ∫ A · F d³r,
 *
 *         where F^α(r) = (1/4π) ε^αβγ ∂_β m̂ · (∂_γ m̂ × m̂) is the
 *         (1/4π)-normalised pullback of the S² area 2-form and A is
 *         the vector potential satisfying ∇×A = F in Coulomb gauge
 *         ∇·A = 0. A is solved by Jacobi iteration on the discrete
 *         Poisson equation ∇²A = −∇×F using **4th-order central
 *         differences** on the periodic grid (∂_α m̂(r) ≈
 *         (-m̂(r+2ê_α) + 8m̂(r+ê_α) - 8m̂(r-ê_α) + m̂(r-2ê_α))/12,
 *         O(h⁴) error).
 *
 *  Convention: F has 1/(4π) baked in, so the Bott-Tu formula
 *  H = (1/(4π)²) ∫ A_FN · F_FN reduces to ∫ A · F directly here
 *  (no extra prefactor).
 *
 *  Validated structural properties:
 *    - **Zero on trivial textures**: uniform m̂(r) = ẑ, smooth small
 *      perturbations of ẑ, and untwisted skyrmion tubes (no z-
 *      dependence) all return H = 0 to machine precision.
 *    - **Linear in twist count**: a skyrmion tube with n full twists
 *      along the z-axis returns H proportional to n.
 *    - **Antisymmetric under twist sign-flip**: H(−n) = −H(n).
 *
 *  Calibration status: validated against the analytic charge-1
 *  Hopfion via stereographic R³ → S³ → S²: m̂_x = (8xz + 4y(r²−1))/D²,
 *  m̂_y = (8yz − 4x(r²−1))/D², m̂_z = (−r⁴ + 6r² − 8z² − 1)/D² with
 *  D = 1 + r². This ansatz has m̂(0) = m̂(∞) = −ẑ uniformly in any
 *  direction, so the texture closes cleanly under PBC; the integer
 *  H = +1 (linking number of the unit-circle m̂ = +ẑ preimage and
 *  the z-axis m̂ = −ẑ preimage). With 4th-order central differences
 *  the numerical result is H = +0.9936 on 64³ with R = 6 — within
 *  0.6% of integer; H = +0.987 at R = 8. Worst-case across R ∈
 *  [4, 12] is ~8% (boundary truncation at R = 12). For textures
 *  that do not compactify cleanly on T³ (e.g., twisted skyrmion
 *  tubes with direction-dependent asymptotic behaviour), H is the
 *  well-defined discrete Whitehead integral but is not required to
 *  be a clean integer.
 *
 *  Use cases:
 *    - Relative Hopf-charge comparison between texture configurations
 *    - Linking-detection in skyrmion-tube linking simulations
 *    - 3D analog of the 2D `_topological_charge_2d` for π₃(S²)
 *
 *  Convergence: the Poisson solve uses Jacobi relaxation with up to
 *  `max_iter` iterations and a relative ℓ²-residual tolerance of
 *  `tol`. Suggested defaults: `tol = 1e-9`, `max_iter = 8·N²`
 *  where N = max(Nx, Ny, Nz). Jacobi convergence rate scales as
 *  (1 − π²/N²) per iteration, so larger lattices need quadratically
 *  more iterations.
 *
 *  @param n_field    flat 3·Nx·Ny·Nz doubles, row-major in
 *                    (iz, iy, ix), with the inner 3 = (m_x, m_y, m_z)
 *                    and |m̂| = 1 (caller responsibility)
 *  @param Nx, Ny, Nz lattice extents (each ≥ 4)
 *  @param tol        Poisson solver convergence tolerance (rel ℓ²)
 *  @param max_iter   Poisson solver iteration cap
 *  @param hopf_out   caller-allocated double — receives H
 *  @return IRREP_OK on success */
IRREP_API irrep_status_t irrep_magnon_hopf_charge_3d(const double *n_field, int Nx, int Ny, int Nz,
                                                       double tol, int max_iter,
                                                       double *hopf_out);

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

/** @brief Magnon spin-Nernst coefficient α^s_xy(T) — the spin-current
 *         analog of the thermal Hall conductivity. For magnons of
 *         spin σ_b along the magnetisation axis (σ_b = ±ℏ),
 *
 *      α^s_xy(T) = -(k_B / ℏ V_uc · (2π)²)
 *                  · Σ_b ∫_BZ d²k · c_1(n_B(ω_b/T)) · Ω_b(k)
 *
 *  where c_1(g) = (1+g)·ln(1+g) − g·ln(g) is the Matsumoto-Murakami
 *  spin-current weight function (distinct from the c_2 thermal-Hall
 *  weight). c_1 limits: c_1(g → 0) → −g ln g → 0 (low T, BE-
 *  suppressed), c_1(g → ∞) → ln g + 1 (high T, log growth).
 *
 *  α^s_xy is the response coefficient to a transverse temperature
 *  gradient: ⟨j_s_x⟩ = −α^s_xy · ∂_y T. It captures the magnon-
 *  Berry-curvature spin transport that, alongside κ_xy, is the
 *  experimental signature of topological magnons. The two
 *  coefficients differ by which weight function dresses Σ_b ∫ Ω_b.
 *
 *  Returned in natural units (k_B = ℏ = 1, A_uc in unit-cell
 *  geometry units).
 *
 *  @param L      LSW handle
 *  @param T      temperature (same units as ω)
 *  @param Nx, Ny BZ grid (typical 50–200)
 *  @return       α^s_xy(T) in natural units, or NaN on error. */
IRREP_API double irrep_magnon_spin_nernst(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny);

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
