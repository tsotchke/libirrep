/* SPDX-License-Identifier: MIT */
/* Total-J projection on spin-½ wavefunctions over N sites.
 *
 * The spin-½ rotation D^{½}(α,β,γ) is the 2×2 matrix
 *
 *   [  e^{-i(α+γ)/2} cos(β/2)   -e^{-i(α-γ)/2} sin(β/2) ]
 *   [  e^{+i(α-γ)/2} sin(β/2)    e^{+i(α+γ)/2} cos(β/2) ]
 *
 * (Sakurai 3.3.21 / Varshalovich 4.7.1). On N sites the rotation operator
 * is the tensor power D^{½⊗N}; applying it factors as a cache-friendly
 * sweep: for each site k ∈ [0, N), for each pair of basis states differing
 * only in qubit k, form the 2×2 matrix-vector product in place.
 *
 * For the character χ_J(α,β,γ) = Tr D^J(α,β,γ) we use the closed form
 *
 *   χ_J(α,β,γ) = sin((2J+1) ω/2) / sin(ω/2),
 *
 * where `ω` is the rotation angle encoded by (α,β,γ) via
 * cos(ω/2) = cos(β/2) cos((α+γ)/2)    (Rose 4.64; Varshalovich 4.7.12).
 *
 * Integration measure: ∫ dα dβ dγ sin(β) / (8π²). For Gauss-Legendre in
 * cos β ∈ [-1,1] we use the standard quadrature; for α and γ, the periodic
 * integrand with its degree ≤ 2·J_max trigonometric content integrates
 * exactly under uniform sampling once n ≥ 2J_max + 1.
 */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/spin_project.h>
#include <irrep/quadrature.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

extern void irrep_set_error_(const char *fmt, ...);

/* -------------------------------------------------------------------------- *
 * Single-qubit rotation applied to all N sites sequentially.                 *
 * For each site k, rewrite ψ as pairs of amplitudes differing only in qubit k*
 * and apply the 2×2 matrix (u00, u01, u10, u11).                             *
 * -------------------------------------------------------------------------- */

static void apply_qubit_(int N, int k, const double _Complex U[4],
                         double _Complex *psi) {
    long long stride = 1LL << k;
    long long dim    = 1LL << N;
    for (long long base = 0; base < dim; base += stride << 1) {
        for (long long off = 0; off < stride; ++off) {
            long long i0 = base + off;
            long long i1 = i0 + stride;
            double _Complex a = psi[i0];
            double _Complex b = psi[i1];
            psi[i0] = U[0] * a + U[1] * b;
            psi[i1] = U[2] * a + U[3] * b;
        }
    }
}

void irrep_spin_half_apply_rotation(int N,
                                    double alpha, double beta, double gamma,
                                    const double _Complex *psi_in,
                                    double _Complex *psi_out) {
    if (N < 1 || !psi_in || !psi_out) return;

    /* The 2×2 D^{½} matrix (Sakurai 3.3.21). */
    double cb = cos(0.5 * beta);
    double sb = sin(0.5 * beta);
    double _Complex e_m_apg = cexp(-0.5 * I * (alpha + gamma));
    double _Complex e_m_amg = cexp(-0.5 * I * (alpha - gamma));
    double _Complex e_p_apg = cexp( 0.5 * I * (alpha + gamma));
    double _Complex e_p_amg = cexp( 0.5 * I * (alpha - gamma));
    double _Complex U[4] = {
        e_m_apg * cb,       -e_m_amg * sb,
        e_p_amg * sb,        e_p_apg * cb
    };

    if (psi_in != psi_out) {
        memcpy(psi_out, psi_in, (size_t)(1LL << N) * sizeof(double _Complex));
    }
    for (int k = 0; k < N; ++k) {
        apply_qubit_(N, k, U, psi_out);
    }
}

/* -------------------------------------------------------------------------- *
 * Character χ_J(α, β, γ).                                                    *
 * ω is the total rotation angle; cos(ω/2) = cos(β/2) cos((α+γ)/2).           *
 * -------------------------------------------------------------------------- */

static double character_J_(int two_J, double alpha, double beta, double gamma) {
    double half_omega_cos = cos(0.5 * beta) * cos(0.5 * (alpha + gamma));
    /* clamp to [-1, 1] to absorb rounding */
    if (half_omega_cos >  1.0) half_omega_cos =  1.0;
    if (half_omega_cos < -1.0) half_omega_cos = -1.0;
    double half_omega = acos(half_omega_cos);

    /* χ_J = sin((2J+1) ω/2) / sin(ω/2) */
    double denom = sin(half_omega);
    if (fabs(denom) < 1e-14) {
        /* ω → 0: χ_J → 2J+1 (trace of the identity D^J) */
        return (double)(two_J + 1);
    }
    double numer = sin((double)(two_J + 1) * half_omega);
    return numer / denom;
}

/* -------------------------------------------------------------------------- *
 * Total-J projection                                                         *
 * -------------------------------------------------------------------------- */

irrep_status_t
irrep_spin_project_spin_half(int two_J_target, int N,
                             int n_alpha, int n_beta, int n_gamma,
                             const double _Complex *psi_in,
                             double _Complex *psi_out) {
    if (N < 1 || N > 24 || two_J_target < 0 ||
        n_alpha < 1 || n_beta < 1 || n_gamma < 1 ||
        !psi_in || !psi_out)
        return IRREP_ERR_INVALID_ARG;

    long long dim = 1LL << N;

    /* β grid via Gauss-Legendre in cos β ∈ [-1, 1] */
    double *cosb = malloc((size_t)n_beta * sizeof(double));
    double *wb   = malloc((size_t)n_beta * sizeof(double));
    if (!cosb || !wb) { free(cosb); free(wb); return IRREP_ERR_OUT_OF_MEMORY; }
    if (!irrep_gauss_legendre(n_beta, cosb, wb)) {
        free(cosb); free(wb); return IRREP_ERR_NOT_IMPLEMENTED;
    }

    /* Working buffers */
    double _Complex *psi_rot = malloc((size_t)dim * sizeof(double _Complex));
    if (!psi_rot) {
        free(cosb); free(wb); return IRREP_ERR_OUT_OF_MEMORY;
    }

    memset(psi_out, 0, (size_t)dim * sizeof(double _Complex));

    /* α, γ in [0, 2π) with uniform-step rectangle rule (periodic → exact
     * for trig polynomials of degree < n). */
    double d_alpha = 2.0 * M_PI / (double)n_alpha;
    double d_gamma = 2.0 * M_PI / (double)n_gamma;

    /* Measure:  (2J+1) / (8π²) · dα · dγ · d(cos β) · w  (Gauss)
     * where α and γ grids carry each an implicit (2π / n) weight. */
    double prefactor = (double)(two_J_target + 1) / (8.0 * M_PI * M_PI)
                     * d_alpha * d_gamma;

    for (int ia = 0; ia < n_alpha; ++ia) {
        double alpha = ia * d_alpha;
        for (int ib = 0; ib < n_beta; ++ib) {
            double beta = acos(cosb[ib]);
            for (int ig = 0; ig < n_gamma; ++ig) {
                double gamma = ig * d_gamma;

                irrep_spin_half_apply_rotation(N, alpha, beta, gamma,
                                               psi_in, psi_rot);

                double chi = character_J_(two_J_target, alpha, beta, gamma);
                double weight = prefactor * wb[ib] * chi;

                for (long long k = 0; k < dim; ++k) {
                    psi_out[k] += weight * psi_rot[k];
                }
            }
        }
    }

    free(psi_rot);
    free(cosb); free(wb);
    return IRREP_OK;
}
