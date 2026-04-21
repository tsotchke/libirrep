/* SPDX-License-Identifier: MIT */
/** @file spin_project.h
 *  @brief Total-J projection on spin-½ wavefunctions over N sites.
 *
 *  Projects an arbitrary wavefunction |ψ⟩ onto a fixed total-J subspace
 *  via the Wigner-D character-weighted SU(2) integral
 *
 *  \f[ \mathcal{P}_J = \frac{2J+1}{8\pi^2} \int\!\!\mathrm{d}\Omega\;
 *      \chi_J^*(\Omega)\; \hat R(\Omega) \f]
 *
 *  where `R(Ω) = [D^{½}(α,β,γ)]^{⊗N}` is the spin-½ tensor-product rotation
 *  on all sites, and `χ_J(Ω) = Tr D^J(α,β,γ)` is the character of the
 *  target irrep.
 *
 *  The integral is discretised as a tensor-product quadrature: uniform
 *  Gaussian-like sampling in `α ∈ [0, 2π)` and `γ ∈ [0, 2π)`, and
 *  Gauss-Legendre (with transformation) in `cos β ∈ [-1, 1]`. See
 *  Varshalovich §4.2 for the integration measure.
 *
 *  ### Heisenberg ground state
 *  On an even number of spin-½ sites with an antiferromagnetic Heisenberg
 *  Hamiltonian, the ground state has `J = 0` (Marshall sign rule; Lieb &
 *  Mattis 1962). The caller sets `two_J_target = 0` to restrict the
 *  variational ansatz to the singlet sector.
 *
 *  ### Cost and scaling
 *  Each rotation costs O(N · 2^N) via in-place single-qubit updates (the
 *  N-fold tensor product structure of D^{½⊗N}). One projection with an
 *  `n_α × n_β × n_γ` grid costs `n_α · n_β · n_γ · O(N · 2^N)`.
 *  For N=14 and a 16·8·16 grid, that's ~6M rotations: about 0.5 s on a
 *  modern desktop. Not practical beyond ED-scale; production NQS work
 *  substitutes MCMC estimation of the character-weighted sum.
 *
 *  ### Convergence
 *  The integrand is a polynomial of degree `2·J_max` in `cos β, cos α,
 *  cos γ` where `J_max = N/2`. Gauss-Legendre in `cos β` is exact at
 *  order `⌈J_max + 1⌉`; uniform sampling in `α, γ` is exact at that
 *  order too (the integrand is a trigonometric polynomial). The caller
 *  picks a grid large enough for their J_max.
 */
#ifndef IRREP_SPIN_PROJECT_H
#define IRREP_SPIN_PROJECT_H

#include <stddef.h>
#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Apply the spin-½ tensor-product rotation
 *         `[D^{½}(α,β,γ)]^{⊗N}` to @p psi_in (length `2^N`),
 *         writing the result to @p psi_out. In-place allowed
 *         (`psi_in == psi_out`).
 *
 *         The rotation follows the physics ZYZ Euler-angle convention used
 *         throughout libirrep. */
IRREP_API void
irrep_spin_half_apply_rotation(int N,
                               double alpha, double beta, double gamma,
                               const double _Complex *psi_in,
                               double _Complex *psi_out);

/** @brief Project @p psi_in onto the total-`J = two_J_target/2` subspace
 *         on `N` spin-½ sites. Uses a tensor-product quadrature of size
 *         `n_alpha × n_beta × n_gamma`.
 *
 *         Output is accumulated into @p psi_out (zeroed internally);
 *         its norm-squared gives the weight of the projected component
 *         of `|ψ⟩` in the target sector.
 *
 *  @param two_J_target doubled target J (0 for singlet, 2 for triplet, …)
 *  @param N            number of spin-½ sites (2 ≤ N ≤ 24)
 *  @param n_alpha      `α` grid size (≥ ⌈J_max + 1⌉)
 *  @param n_beta       `β` grid size (≥ ⌈J_max + 1⌉)
 *  @param n_gamma      `γ` grid size (≥ ⌈J_max + 1⌉)
 *  @param psi_in       input wavefunction, length `2^N`
 *  @param psi_out      output wavefunction, length `2^N`
 *  @return             #IRREP_OK, or an error code. */
IRREP_API irrep_status_t
irrep_spin_project_spin_half(int two_J_target, int N,
                             int n_alpha, int n_beta, int n_gamma,
                             const double _Complex *psi_in,
                             double _Complex *psi_out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SPIN_PROJECT_H */
