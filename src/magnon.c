/* SPDX-License-Identifier: MIT */
/* Linearized spin-wave theory and Berry curvature.
 *
 * Holstein-Primakoff bosonisation of a Heisenberg + DMI + anisotropy
 * Hamiltonian on top of a collinear FM ordering along z. The LSW
 * Hamiltonian H(k) is an n_sub × n_sub Hermitian matrix per k. Magnon
 * dispersion ω_b(k) = eigenvalues of H(k); Berry curvature is computed
 * via the gauge-invariant Wilson-loop method on a 4-point plaquette in
 * k-space (Fukui-Hatsugai-Suzuki 2005). */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/magnon.h>
#include <irrep/rdm.h> /* irrep_hermitian_eigvals */

extern void irrep_set_error_(const char *fmt, ...);

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct irrep_magnon_lsw {
    int                  n_sub;
    double               S;
    double               a1[2], a2[2];
    /* reciprocal vectors, computed at build time */
    double               b1[2], b2[2];
    int                  n_bonds;
    irrep_magnon_bond_t *bonds;
    double               Kz;
};

/* Reciprocal of a 2D primitive basis. */
static void reciprocal_2d_(const double a1[2], const double a2[2], double b1[2], double b2[2]) {
    double det = a1[0] * a2[1] - a1[1] * a2[0];
    double s = (2.0 * M_PI) / det;
    b1[0] = s * a2[1];
    b1[1] = -s * a2[0];
    b2[0] = -s * a1[1];
    b2[1] = s * a1[0];
}

irrep_magnon_lsw_t *irrep_magnon_lsw_new(int n_sub, double S, const double a1[2],
                                         const double a2[2], int n_bonds,
                                         const irrep_magnon_bond_t *bonds, double Kz) {
    if (n_sub <= 0 || n_sub > 64 || S <= 0 || n_bonds < 0) {
        irrep_set_error_("irrep_magnon_lsw_new: invalid arguments");
        return NULL;
    }
    irrep_magnon_lsw_t *L = calloc(1, sizeof(*L));
    if (!L)
        return NULL;
    L->n_sub = n_sub;
    L->S = S;
    L->a1[0] = a1[0]; L->a1[1] = a1[1];
    L->a2[0] = a2[0]; L->a2[1] = a2[1];
    reciprocal_2d_(a1, a2, L->b1, L->b2);
    L->Kz = Kz;
    L->n_bonds = n_bonds;
    if (n_bonds > 0) {
        L->bonds = malloc((size_t)n_bonds * sizeof(*L->bonds));
        if (!L->bonds) {
            free(L);
            return NULL;
        }
        memcpy(L->bonds, bonds, (size_t)n_bonds * sizeof(*L->bonds));
    }
    return L;
}

void irrep_magnon_lsw_free(irrep_magnon_lsw_t *L) {
    if (!L)
        return;
    free(L->bonds);
    free(L);
}

int irrep_magnon_lsw_num_bands(const irrep_magnon_lsw_t *L) {
    return L ? L->n_sub : 0;
}

/* Build the LSW Hamiltonian H(k) — a Hermitian n_sub × n_sub matrix.
 *
 * For a collinear FM along z, Holstein-Primakoff gives:
 *   S_i^z = S - a_i^† a_i,   S_i^+ = √(2S) a_i,   S_i^- = √(2S) a_i^†.
 *
 * Convention: H = J Σ_<ij> S_i · S_j + Σ_<ij> D · (S_i × S_j) - K_z Σ (S_i^z)².
 * For Heisenberg FM stability we need J < 0; the LSW dispersion is positive
 * iff the assumed FM ground state is energetically favoured.
 *
 * Bilinear HP expansion:
 *
 *   S_i · S_j = S² - S(n_i + n_j) + S(a_i^† a_j + a_j^† a_i) + (n_i n_j: drop)
 *
 *   (S_i × S_j)_z = (1/2i)(S_i^- S_j^+ - S_i^+ S_j^-)
 *                ≈ -i S (a_i^† a_j - a_j^† a_i)
 *
 *   -K_z (S_i^z)²  ≈ -K_z S² + 2 K_z S n_i      [diagonal: +2 K_z S]
 *
 * So per unique bond (i, j, t) with t = r_j - r_i:
 *   H[i,i] += -S·J    H[j,j] += -S·J          (from -S J (n_i + n_j))
 *   H[i,j] += +S·J·e^{ik·t} - i S D_z e^{ik·t}
 *   H[j,i] += conj(...)                        (Hermiticity; D_ji = -D_ij
 *                                                 already handled by conj)
 *
 * Both orientations are added explicitly so the matrix is Hermitian after
 * a single sweep of the bond list. */
static void build_H_(const irrep_magnon_lsw_t *L, double kx, double ky,
                     double _Complex *H_out) {
    int n = L->n_sub;
    memset(H_out, 0, (size_t)n * n * sizeof(double _Complex));
    double S = L->S;
    /* Diagonal: uniaxial anisotropy. */
    for (int a = 0; a < n; ++a)
        H_out[a * n + a] += 2.0 * L->Kz * S;

    for (int b = 0; b < L->n_bonds; ++b) {
        const irrep_magnon_bond_t *bd = &L->bonds[b];
        int a_idx = bd->bi, b_idx = bd->bj;
        double tx = bd->delta_x * L->a1[0] + bd->delta_y * L->a2[0];
        double ty = bd->delta_x * L->a1[1] + bd->delta_y * L->a2[1];
        double phase_arg = kx * tx + ky * ty;
        double _Complex eikt = cos(phase_arg) + I * sin(phase_arg);

        /* Heisenberg J: hopping coefficient is +S·J on the off-diagonal
         * (sign matches the +J(a_i^†a_j + a_j^†a_i) bilinear). */
        double _Complex Jhop = +S * bd->J * eikt;
        /* DMI z-component: bilinear -i S D_z (a_i^†a_j - a_j^†a_i) gives
         * -i S D_z to H[i,j] and +i S D_z to H[j,i]. */
        double _Complex Dhop = -I * S * bd->D[2] * eikt;

        H_out[a_idx * n + b_idx] += Jhop + Dhop;
        H_out[b_idx * n + a_idx] += conj(Jhop) + conj(Dhop);

        /* Diagonal "self-energy" -S·J(n_i + n_j) per unique bond. */
        H_out[a_idx * n + a_idx] += -S * bd->J;
        H_out[b_idx * n + b_idx] += -S * bd->J;
    }
}

/* Diagonalise a Hermitian n×n matrix using libirrep's existing
 * cyclic-Jacobi solver. The solver returns only eigenvalues; for
 * eigenvectors we re-implement a Hermitian Jacobi sweep here. ~80
 * lines but kept local to magnon.c since `irrep_hermitian_eigvals`
 * doesn't expose eigenvectors. */
static void hermitian_eig_(int n, const double _Complex *A_in, double *evals,
                           double _Complex *evecs) {
    /* Allocate a working copy. */
    double _Complex *A = malloc((size_t)n * n * sizeof *A);
    double _Complex *V = malloc((size_t)n * n * sizeof *V);
    memcpy(A, A_in, (size_t)n * n * sizeof *A);
    for (int i = 0; i < n * n; ++i)
        V[i] = 0;
    for (int i = 0; i < n; ++i)
        V[i * n + i] = 1.0;

    /* Standard Jacobi-Hermitian sweep (Wilkinson 1965). */
    for (int sweep = 0; sweep < 200; ++sweep) {
        double max_off = 0;
        int    p = 0, q = 1;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                double aij = cabs(A[i * n + j]);
                if (aij > max_off) {
                    max_off = aij;
                    p = i; q = j;
                }
            }
        if (max_off < 1e-15)
            break;

        /* Hermitian 2×2 sub-block: app, aqq real (diag); apq complex. */
        double          app = creal(A[p * n + p]);
        double          aqq = creal(A[q * n + q]);
        double _Complex apq = A[p * n + q];
        double          alpha = cabs(apq);
        if (alpha < 1e-300)
            continue;
        double _Complex phase = apq / alpha; /* unit-modulus phase */

        double theta = 0.5 * atan2(2.0 * alpha, app - aqq);
        double c = cos(theta), s = sin(theta);

        /* Update the rows / columns p, q of A and V. */
        for (int k = 0; k < n; ++k) {
            double _Complex Akp = A[k * n + p];
            double _Complex Akq = A[k * n + q];
            A[k * n + p] = c * Akp + s * conj(phase) * Akq;
            A[k * n + q] = -s * phase * Akp + c * Akq;
        }
        for (int k = 0; k < n; ++k) {
            double _Complex Apk = A[p * n + k];
            double _Complex Aqk = A[q * n + k];
            A[p * n + k] = c * Apk + s * phase * Aqk;
            A[q * n + k] = -s * conj(phase) * Apk + c * Aqk;
        }
        for (int k = 0; k < n; ++k) {
            double _Complex Vkp = V[k * n + p];
            double _Complex Vkq = V[k * n + q];
            V[k * n + p] = c * Vkp + s * conj(phase) * Vkq;
            V[k * n + q] = -s * phase * Vkp + c * Vkq;
        }
    }

    /* Eigenvalues are now on the diagonal. Sort ascending. */
    int *idx = malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; ++i)
        idx[i] = i;
    for (int i = 0; i < n - 1; ++i)
        for (int j = i + 1; j < n; ++j)
            if (creal(A[idx[j] * n + idx[j]]) < creal(A[idx[i] * n + idx[i]])) {
                int t = idx[i]; idx[i] = idx[j]; idx[j] = t;
            }
    for (int i = 0; i < n; ++i) {
        evals[i] = creal(A[idx[i] * n + idx[i]]);
        for (int k = 0; k < n; ++k)
            evecs[i * n + k] = V[k * n + idx[i]];
    }
    free(idx);
    free(A);
    free(V);
}

irrep_status_t irrep_magnon_dispersion(const irrep_magnon_lsw_t *L, double kx, double ky,
                                        double *omega_out, double _Complex *u_out) {
    if (!L || !omega_out || !u_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    double _Complex *H = malloc((size_t)n * n * sizeof *H);
    if (!H)
        return IRREP_ERR_OUT_OF_MEMORY;
    build_H_(L, kx, ky, H);
    hermitian_eig_(n, H, omega_out, u_out);
    free(H);
    return IRREP_OK;
}

/* Inner product of two band eigenvectors (length n_sub). */
static double _Complex inner_(const double _Complex *u, const double _Complex *v, int n) {
    double _Complex s = 0;
    for (int i = 0; i < n; ++i)
        s += conj(u[i]) * v[i];
    return s;
}

irrep_status_t irrep_magnon_berry(const irrep_magnon_lsw_t *L, double kx, double ky,
                                  double delta_k, double *berry_out) {
    if (!L || !berry_out || delta_k <= 0)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    /* 4 corners of the plaquette around (kx, ky). */
    double _Complex *u = malloc((size_t)4 * n * n * sizeof *u);
    double          *w = malloc((size_t)4 * n * sizeof *w);
    if (!u || !w) { free(u); free(w); return IRREP_ERR_OUT_OF_MEMORY; }

    double dk = delta_k;
    double kxs[4] = {kx - dk, kx + dk, kx + dk, kx - dk};
    double kys[4] = {ky - dk, ky - dk, ky + dk, ky + dk};
    for (int p = 0; p < 4; ++p) {
        irrep_magnon_dispersion(L, kxs[p], kys[p], w + p * n, u + p * n * n);
    }

    /* Phase-fix: align each pair (u_p, u_{p+1}) by a global phase to
     * preserve gauge invariance through the loop. We use the
     * Fukui-Hatsugai-Suzuki link variable formulation. */
    for (int b = 0; b < n; ++b) {
        double _Complex z = 1.0;
        for (int p = 0; p < 4; ++p) {
            int q = (p + 1) % 4;
            double _Complex ip = inner_(u + p * n * n + b * n, u + q * n * n + b * n, n);
            double mod = cabs(ip);
            if (mod < 1e-300) {
                /* Band crossing — Berry curvature ill-defined; signal NaN. */
                berry_out[b] = NAN;
                z = 0;
                break;
            }
            z *= ip / mod;
        }
        if (z != 0) {
            /* Loop integral = arg(z); curvature = arg / area. */
            berry_out[b] = atan2(cimag(z), creal(z)) / (4.0 * dk * dk);
        }
    }
    free(u);
    free(w);
    return IRREP_OK;
}

/* Dilogarithm Li_2(-x) for x ≥ 0. Taylor series for x < 1; for x > 1,
 * use the Landen transformation Li₂(-x) + Li₂(-1/x) = -π²/6 - ½(ln x)².
 * Convergence ~ 18 terms for x = 0.99; one recursion deep for x > 1. */
static double dilog_neg_(double x) {
    if (x < 1e-15)
        return 0.0;
    if (x < 1.0) {
        double sum = 0.0;
        double term = -x;
        for (int n = 1; n < 200; ++n) {
            sum += term / ((double)n * (double)n);
            term *= -x;
            if (fabs(term) < 1e-18 * fabs(sum))
                break;
        }
        return sum;
    }
    /* x > 1: invert. */
    double l = log(x);
    return -M_PI * M_PI / 6.0 - 0.5 * l * l - dilog_neg_(1.0 / x);
}

/* Matsumoto-Murakami c₁ function for the spin-Nernst coefficient:
 *   c₁(g) = (1+g) ln(1+g) − g ln g
 * Limits: c₁(g→0) → -g ln g → 0 (BE-suppressed), c₁(g→∞) → ln g + 1
 * (logarithmic growth). Smooth and monotonic. */
static double c1_mm_(double g) {
    if (g < 1e-300)
        return 0.0;
    /* For g < 1e-15, -g ln g is what dominates; (1+g) ln(1+g) ≈ g - g²/2. */
    if (g < 1e-15)
        return -g * log(g);
    return (1.0 + g) * log1p(g) - g * log(g);
}

/* Matsumoto-Murakami c₂ function:
 *   c₂(g) = (1+g)·[ln((1+g)/g)]² − (ln g)² − 2 Li₂(−g)
 * Limits: c₂(g→0) → 2g (linear), c₂(g→∞) → π²/3. Smooth and monotonic. */
static double c2_mm_(double g) {
    if (g < 1e-15)
        return 2.0 * g; /* leading order: avoids ln(0) */
    if (g > 1.0e15)
        return M_PI * M_PI / 3.0; /* asymptotic limit */
    double l_inv = log((1.0 + g) / g);
    double l = log(g);
    return (1.0 + g) * l_inv * l_inv - l * l - 2.0 * dilog_neg_(g);
}

irrep_status_t irrep_magnon_chern(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                   double *chern_out) {
    if (!L || !chern_out || Nx <= 0 || Ny <= 0)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    /* Sample u_b(k) on the Nx × Ny grid in reduced coordinates
     * f_x ∈ [0, 1) and f_y ∈ [0, 1), with k = f_x · b1 + f_y · b2.
     * For each grid plaquette, compute the FHS link variable and
     * accumulate. */
    double _Complex *u_grid = malloc((size_t)Nx * Ny * n * n * sizeof *u_grid);
    double          *w_grid = malloc((size_t)Nx * Ny * n * sizeof *w_grid);
    if (!u_grid || !w_grid) {
        free(u_grid);
        free(w_grid);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            int    p = (iy * Nx + ix);
            irrep_magnon_dispersion(L, kx, ky, w_grid + p * n, u_grid + p * n * n);
        }

    for (int b = 0; b < n; ++b) {
        double total = 0;
        for (int iy = 0; iy < Ny; ++iy) {
            for (int ix = 0; ix < Nx; ++ix) {
                int p00 = (iy * Nx + ix);
                int p10 = (iy * Nx + (ix + 1) % Nx);
                int p11 = (((iy + 1) % Ny) * Nx + (ix + 1) % Nx);
                int p01 = (((iy + 1) % Ny) * Nx + ix);
                double _Complex z = 1.0;
                int             skipped = 0;
                int             pairs[4][2] = {{p00, p10}, {p10, p11}, {p11, p01}, {p01, p00}};
                for (int pp = 0; pp < 4; ++pp) {
                    double _Complex ip = inner_(u_grid + pairs[pp][0] * n * n + b * n,
                                                u_grid + pairs[pp][1] * n * n + b * n, n);
                    double          mod = cabs(ip);
                    if (mod < 1e-300) {
                        skipped = 1;
                        break;
                    }
                    z *= ip / mod;
                }
                if (!skipped)
                    total += atan2(cimag(z), creal(z));
            }
        }
        chern_out[b] = total / (2.0 * M_PI);
    }
    free(u_grid);
    free(w_grid);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_neutron_qomega_map(const irrep_magnon_lsw_t *L,
                                                const double (*qpath)[2], int n_q,
                                                double omega_min, double omega_max, int n_omega,
                                                double eta, double *intensity_out) {
    if (!L || !qpath || n_q <= 0 || n_omega <= 0 || eta <= 0 || omega_max <= omega_min ||
        !intensity_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(intensity_out, 0, (size_t)n_q * (size_t)n_omega * sizeof *intensity_out);

    double *omega = malloc((size_t)n * sizeof *omega);
    double *S_perp = malloc((size_t)n * sizeof *S_perp);
    if (!omega || !S_perp) {
        free(omega);
        free(S_perp);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    double dw = (omega_max - omega_min) / n_omega;

    for (int iq = 0; iq < n_q; ++iq) {
        irrep_magnon_structure_factor(L, qpath[iq][0], qpath[iq][1], omega, S_perp);
        for (int jw = 0; jw < n_omega; ++jw) {
            double w = omega_min + (jw + 0.5) * dw;
            double total = 0;
            for (int b = 0; b < n; ++b) {
                double dx = w - omega[b];
                /* Unit-area Lorentzian */
                total += S_perp[b] * (eta / M_PI) / (dx * dx + eta * eta);
            }
            intensity_out[iq * n_omega + jw] = total;
        }
    }
    free(omega);
    free(S_perp);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_group_velocity(const irrep_magnon_lsw_t *L, double kx, double ky,
                                            double h, double *vx_out, double *vy_out) {
    if (!L || !vx_out || !vy_out || h <= 0)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    double *omega_pp = malloc((size_t)n * sizeof *omega_pp);
    double *omega_mm = malloc((size_t)n * sizeof *omega_mm);
    double _Complex *u_dummy = malloc((size_t)n * n * sizeof *u_dummy);
    if (!omega_pp || !omega_mm || !u_dummy) {
        free(omega_pp);
        free(omega_mm);
        free(u_dummy);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    /* x-direction */
    irrep_magnon_dispersion(L, kx + h, ky, omega_pp, u_dummy);
    irrep_magnon_dispersion(L, kx - h, ky, omega_mm, u_dummy);
    for (int b = 0; b < n; ++b)
        vx_out[b] = (omega_pp[b] - omega_mm[b]) / (2.0 * h);

    /* y-direction */
    irrep_magnon_dispersion(L, kx, ky + h, omega_pp, u_dummy);
    irrep_magnon_dispersion(L, kx, ky - h, omega_mm, u_dummy);
    for (int b = 0; b < n; ++b)
        vy_out[b] = (omega_pp[b] - omega_mm[b]) / (2.0 * h);

    free(omega_pp);
    free(omega_mm);
    free(u_dummy);
    return IRREP_OK;
}

double irrep_magnon_free_energy(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny) {
    if (!L || T <= 0 || Nx <= 0 || Ny <= 0) {
        irrep_set_error_("irrep_magnon_free_energy: invalid arguments");
        return NAN;
    }
    int     n = L->n_sub;
    double *omega = malloc((size_t)n * sizeof *omega);
    double _Complex *u = malloc((size_t)n * n * sizeof *u);
    if (!omega || !u) {
        free(omega);
        free(u);
        return NAN;
    }
    double accum = 0;
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            irrep_magnon_dispersion(L, kx, ky, omega, u);
            for (int b = 0; b < n; ++b) {
                double w = omega[b];
                if (w < 1e-10)
                    continue;
                double x = w / T;
                /* For very large x: ln(1 − e^{−x}) ≈ −e^{−x}, contribution
                 * vanishes exponentially. */
                if (x > 700.0)
                    continue;
                /* For very small x: ln(1 − e^{−x}) → ln(x − x²/2) ≈ ln x.
                 * Use log1p(-exp(-x)) for numerical stability. */
                accum += log1p(-exp(-x));
            }
        }
    free(omega);
    free(u);
    /* F(T) = T · (1/N_BZ) · Σ_b Σ_k ln(1 - exp(-ω_b/T)) per cell. */
    return T * accum / ((double)Nx * (double)Ny);
}

double irrep_magnon_susceptibility(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny) {
    if (!L || T <= 0 || Nx <= 0 || Ny <= 0) {
        irrep_set_error_("irrep_magnon_susceptibility: invalid arguments");
        return NAN;
    }
    int     n = L->n_sub;
    double *omega = malloc((size_t)n * sizeof *omega);
    double _Complex *u = malloc((size_t)n * n * sizeof *u);
    if (!omega || !u) {
        free(omega);
        free(u);
        return NAN;
    }
    double accum = 0;
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            irrep_magnon_dispersion(L, kx, ky, omega, u);
            for (int b = 0; b < n; ++b) {
                double w = omega[b];
                if (w < 1e-10)
                    continue;
                double x = w / T;
                if (x > 700.0)
                    continue;
                double n_BE = 1.0 / (exp(x) - 1.0);
                accum += n_BE * (1.0 + n_BE);
            }
        }
    free(omega);
    free(u);
    /* χ(T) = (1/T) · (1/N_BZ) · Σ_b Σ_k n_B(ω_b/T)(1+n_B(ω_b/T))
     * BZ sum / N_BZ approximates ∫ d²k/(2π)² over the BZ. */
    return accum / (T * (double)Nx * (double)Ny);
}

double irrep_magnon_magnetization(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny) {
    if (!L || T <= 0 || Nx <= 0 || Ny <= 0) {
        irrep_set_error_("irrep_magnon_magnetization: invalid arguments");
        return NAN;
    }
    int     n = L->n_sub;
    double *omega = malloc((size_t)n * sizeof *omega);
    double _Complex *u = malloc((size_t)n * n * sizeof *u);
    if (!omega || !u) {
        free(omega);
        free(u);
        return NAN;
    }
    double thermal_pop = 0;
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            irrep_magnon_dispersion(L, kx, ky, omega, u);
            for (int b = 0; b < n; ++b) {
                double w = omega[b];
                if (w < 1e-10)
                    continue;
                double x = w / T;
                if (x > 700.0)
                    continue;
                thermal_pop += 1.0 / (exp(x) - 1.0);
            }
        }
    free(omega);
    free(u);
    /* M(T) = S − (1/n_sub) · (1/N_BZ) · Σ_b Σ_k n_B(ω_b(k)/T)
     * The BZ sum / N_BZ approximates ∫ d²k/(2π)² over the BZ. */
    return L->S - thermal_pop / ((double)Nx * (double)Ny * (double)n);
}

double irrep_magnon_specific_heat(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny) {
    if (!L || T <= 0 || Nx <= 0 || Ny <= 0) {
        irrep_set_error_("irrep_magnon_specific_heat: invalid arguments");
        return NAN;
    }
    int     n = L->n_sub;
    double *omega = malloc((size_t)n * sizeof *omega);
    double _Complex *u = malloc((size_t)n * n * sizeof *u);
    if (!omega || !u) {
        free(omega);
        free(u);
        return NAN;
    }
    double cv = 0;
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            irrep_magnon_dispersion(L, kx, ky, omega, u);
            for (int b = 0; b < n; ++b) {
                double w = omega[b];
                if (w < 1e-10)
                    continue; /* BE singularity; measure-zero, skip */
                double x = w / T;
                /* For very large x (low T regime), exp(x) overflows but
                 * the contribution is exponentially suppressed; treat as 0. */
                if (x > 700.0)
                    continue;
                double ex = exp(x);
                double n_BE = 1.0 / (ex - 1.0);
                cv += x * x * n_BE * (1.0 + n_BE);
            }
        }
    free(omega);
    free(u);
    /* Average over BZ grid → integral ∫ d²k/(2π)² in our discretisation.
     * The grid has Nx · Ny points; each represents a fraction 1/(Nx·Ny)
     * of the BZ. So sum/Nx/Ny is the per-cell C_V. */
    return cv / ((double)Nx * (double)Ny);
}

irrep_status_t irrep_magnon_dos(const irrep_magnon_lsw_t *L, int Nx, int Ny, double omega_min,
                                 double omega_max, int n_bins, double *dos_out) {
    if (!L || Nx <= 0 || Ny <= 0 || n_bins <= 0 || omega_max <= omega_min || !dos_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(dos_out, 0, (size_t)n_bins * sizeof *dos_out);

    double          bin_w = (omega_max - omega_min) / n_bins;
    double         *omega = malloc((size_t)n * sizeof *omega);
    double _Complex *u = malloc((size_t)n * n * sizeof *u);
    if (!omega || !u) {
        free(omega);
        free(u);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    long long total_evals = 0;
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            irrep_magnon_dispersion(L, kx, ky, omega, u);
            for (int b = 0; b < n; ++b) {
                int idx = (int)((omega[b] - omega_min) / bin_w);
                if (idx >= 0 && idx < n_bins)
                    dos_out[idx] += 1.0;
                ++total_evals;
            }
        }
    /* Normalise so ∫ D(ω) dω = n_sub.
     * total_evals = Nx · Ny · n_sub. Bin sum currently = total_evals.
     * Want sum · bin_w = n_sub, so divide by (Nx · Ny · bin_w). */
    double norm = 1.0 / ((double)Nx * (double)Ny * bin_w);
    for (int i = 0; i < n_bins; ++i)
        dos_out[i] *= norm;
    (void)total_evals;
    free(omega);
    free(u);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_structure_factor(const irrep_magnon_lsw_t *L, double qx, double qy,
                                              double *omega_out, double *S_perp_out) {
    if (!L || !omega_out || !S_perp_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    double _Complex *u = malloc((size_t)n * n * sizeof *u);
    if (!u)
        return IRREP_ERR_OUT_OF_MEMORY;
    irrep_status_t st = irrep_magnon_dispersion(L, qx, qy, omega_out, u);
    if (st != IRREP_OK) {
        free(u);
        return st;
    }
    /* For each band b, S_perp_b = 2S · |Σ_α u_b(q)_α|² . Eigenvectors
     * are stored as rows of u: u[b*n + α]. */
    for (int b = 0; b < n; ++b) {
        double _Complex sum = 0;
        for (int a = 0; a < n; ++a)
            sum += u[b * n + a];
        S_perp_out[b] = 2.0 * L->S * (creal(sum) * creal(sum) + cimag(sum) * cimag(sum));
    }
    free(u);
    return IRREP_OK;
}

/* Build H(k) for the 3D extension: t = delta_x·a₁ + delta_y·a₂ +
 * delta_z·a₃. Otherwise identical to build_H_(). */
static void build_H_3d_(const irrep_magnon_lsw_t *L, const double a3[3], double kx, double ky,
                         double kz, double _Complex *H_out) {
    int n = L->n_sub;
    memset(H_out, 0, (size_t)n * n * sizeof(double _Complex));
    double S = L->S;
    for (int a = 0; a < n; ++a)
        H_out[a * n + a] += 2.0 * L->Kz * S;

    for (int b = 0; b < L->n_bonds; ++b) {
        const irrep_magnon_bond_t *bd = &L->bonds[b];
        int    a_idx = bd->bi, b_idx = bd->bj;
        /* a₁, a₂ have only x, y components (2D structure). a₃ is fully
         * 3D. */
        double tx = bd->delta_x * L->a1[0] + bd->delta_y * L->a2[0] + bd->delta_z * a3[0];
        double ty = bd->delta_x * L->a1[1] + bd->delta_y * L->a2[1] + bd->delta_z * a3[1];
        double tz = bd->delta_z * a3[2];
        double phase_arg = kx * tx + ky * ty + kz * tz;
        double _Complex eikt = cos(phase_arg) + I * sin(phase_arg);

        double _Complex Jhop = +S * bd->J * eikt;
        double _Complex Dhop = -I * S * bd->D[2] * eikt;

        H_out[a_idx * n + b_idx] += Jhop + Dhop;
        H_out[b_idx * n + a_idx] += conj(Jhop) + conj(Dhop);

        H_out[a_idx * n + a_idx] += -S * bd->J;
        H_out[b_idx * n + b_idx] += -S * bd->J;
    }
}

irrep_status_t irrep_magnon_chern_3d_slice_kz(const irrep_magnon_lsw_t *L, const double a3[3],
                                               double kz, int Nx, int Ny, double *chern_out) {
    if (!L || !a3 || !chern_out || Nx <= 0 || Ny <= 0)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    double _Complex *u_grid = malloc((size_t)Nx * Ny * n * n * sizeof *u_grid);
    double          *w_grid = malloc((size_t)Nx * Ny * n * sizeof *w_grid);
    if (!u_grid || !w_grid) {
        free(u_grid);
        free(w_grid);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    /* Sample (k_x, k_y) on a uniform grid in reduced (b1, b2)
     * coordinates, with k_z held fixed. */
    double _Complex *H = malloc((size_t)n * n * sizeof *H);
    if (!H) {
        free(u_grid);
        free(w_grid);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            int    p = iy * Nx + ix;
            build_H_3d_(L, a3, kx, ky, kz, H);
            hermitian_eig_(n, H, w_grid + p * n, u_grid + p * n * n);
        }
    free(H);

    for (int b = 0; b < n; ++b) {
        double total = 0;
        for (int iy = 0; iy < Ny; ++iy) {
            for (int ix = 0; ix < Nx; ++ix) {
                int p00 = iy * Nx + ix;
                int p10 = iy * Nx + (ix + 1) % Nx;
                int p11 = ((iy + 1) % Ny) * Nx + (ix + 1) % Nx;
                int p01 = ((iy + 1) % Ny) * Nx + ix;
                int pairs[4][2] = {{p00, p10}, {p10, p11}, {p11, p01}, {p01, p00}};
                double _Complex z = 1.0;
                int             skipped = 0;
                for (int pp = 0; pp < 4; ++pp) {
                    double _Complex ip = inner_(u_grid + pairs[pp][0] * n * n + b * n,
                                                u_grid + pairs[pp][1] * n * n + b * n, n);
                    double          mod = cabs(ip);
                    if (mod < 1e-300) {
                        skipped = 1;
                        break;
                    }
                    z *= ip / mod;
                }
                if (!skipped)
                    total += atan2(cimag(z), creal(z));
            }
        }
        chern_out[b] = total / (2.0 * M_PI);
    }
    free(u_grid);
    free(w_grid);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_dispersion_3d(const irrep_magnon_lsw_t *L, const double a3[3],
                                           double kx, double ky, double kz, double *omega_out,
                                           double _Complex *u_out) {
    if (!L || !a3 || !omega_out || !u_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    double _Complex *H = malloc((size_t)n * n * sizeof *H);
    if (!H)
        return IRREP_ERR_OUT_OF_MEMORY;
    build_H_3d_(L, a3, kx, ky, kz, H);
    hermitian_eig_(n, H, omega_out, u_out);
    free(H);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_wilson_spectrum(const irrep_magnon_lsw_t *L, double kx, int Ny,
                                             double *theta_out) {
    if (!L || !theta_out || Ny < 4)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    /* Sample u_b(k_x, k_y_i) at Ny equally-spaced points along the b₂
     * BZ direction (cartesian: k_y_i = (i/Ny) · b₂_y, plus the
     * k_x-shift along b₂ if a₁ and a₂ are non-orthogonal). For
     * Wannier-center flow we want a closed loop in cartesian k-space
     * along the b₂ direction at fixed (k_x · b₁ / |b₁|²). */
    double _Complex *u_grid = malloc((size_t)Ny * n * n * sizeof *u_grid);
    double          *w_grid = malloc((size_t)Ny * n * sizeof *w_grid);
    if (!u_grid || !w_grid) {
        free(u_grid);
        free(w_grid);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    /* The Wilson loop is along the b₂ direction. For each band b,
     * compute u_b at points k = k_fixed + (i/Ny)·b₂ (closed loop since
     * adding b₂ wraps the BZ). The fixed k component is the projection
     * of (kx, 0) onto b₁: f₁ = (kx · b₁_x + 0 · b₁_y) / |b₁|² (we just
     * use kx along x for simplicity). */
    for (int i = 0; i < Ny; ++i) {
        double fy = (double)i / Ny;
        double k1 = kx + fy * L->b2[0];
        double k2 = 0.0 + fy * L->b2[1];
        irrep_magnon_dispersion(L, k1, k2, w_grid + i * n, u_grid + i * n * n);
    }

    for (int b = 0; b < n; ++b) {
        double _Complex W = 1.0;
        for (int i = 0; i < Ny; ++i) {
            int    j = (i + 1) % Ny;
            double _Complex ip = 0;
            for (int k = 0; k < n; ++k)
                ip += conj(u_grid[i * n * n + b * n + k]) * u_grid[j * n * n + b * n + k];
            double mod = cabs(ip);
            if (mod < 1e-300) {
                /* Band crossing — Wilson loop phase ill-defined. */
                theta_out[b] = NAN;
                W = 0;
                break;
            }
            W *= ip / mod; /* normalise for unitarity */
        }
        if (W != 0.0)
            theta_out[b] = atan2(cimag(W), creal(W));
    }
    free(u_grid);
    free(w_grid);
    return IRREP_OK;
}

/* Cholesky factorisation M = K^† K for n×n positive-definite Hermitian M.
 * K is upper triangular. Returns 0 on success, -1 if a non-positive
 * pivot is found (M not positive-definite — caller's assumed ground
 * state is unstable). */
static int cholesky_(int n, const double _Complex *M, double _Complex *K) {
    memset(K, 0, (size_t)n * n * sizeof *K);
    for (int j = 0; j < n; ++j) {
        double _Complex s = M[j * n + j];
        for (int k = 0; k < j; ++k)
            s -= conj(K[k * n + j]) * K[k * n + j];
        double d = creal(s);
        if (d <= 1e-12)
            return -1;
        K[j * n + j] = sqrt(d);
        for (int i = j + 1; i < n; ++i) {
            double _Complex sum = M[j * n + i];
            for (int k = 0; k < j; ++k)
                sum -= conj(K[k * n + j]) * K[k * n + i];
            K[j * n + i] = sum / K[j * n + j];
        }
    }
    return 0;
}

irrep_status_t irrep_magnon_dispersion_general(const irrep_magnon_lsw_t *L,
                                                const int *sublattice_signs, double kx, double ky,
                                                double *omega_out) {
    if (!L || !sublattice_signs || !omega_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    int N = 2 * n;
    /* Validate signs. */
    for (int a = 0; a < n; ++a)
        if (sublattice_signs[a] != +1 && sublattice_signs[a] != -1) {
            irrep_set_error_("irrep_magnon_dispersion_general: sublattice_signs[%d] must be ±1", a);
            return IRREP_ERR_INVALID_ARG;
        }

    double _Complex *M = calloc((size_t)N * N, sizeof *M);
    if (!M)
        return IRREP_ERR_OUT_OF_MEMORY;
    double S = L->S;

    /* Anisotropy on every site of particle and hole blocks. */
    for (int a = 0; a < n; ++a) {
        M[a * N + a] += 2.0 * L->Kz * S;
        M[(a + n) * N + (a + n)] += 2.0 * L->Kz * S;
    }

    for (int b = 0; b < L->n_bonds; ++b) {
        const irrep_magnon_bond_t *bd = &L->bonds[b];
        int i = bd->bi, j = bd->bj;
        int sigma = sublattice_signs[i] * sublattice_signs[j];
        double tx = bd->delta_x * L->a1[0] + bd->delta_y * L->a2[0];
        double ty = bd->delta_x * L->a1[1] + bd->delta_y * L->a2[1];
        double phase = kx * tx + ky * ty;
        double _Complex eikt = cos(phase) + I * sin(phase);
        double _Complex eikt_neg = cos(phase) - I * sin(phase);

        if (sigma > 0) {
            /* Parallel sublattices: same as FM construction. */
            M[i * N + i] += -S * bd->J;
            M[j * N + j] += -S * bd->J;
            M[(i + n) * N + (i + n)] += -S * bd->J;
            M[(j + n) * N + (j + n)] += -S * bd->J;

            double _Complex hop_part = (S * bd->J - I * S * bd->D[2]) * eikt;
            M[i * N + j] += hop_part;
            M[j * N + i] += conj(hop_part);

            /* Hole block A(-k)^T: J part same, D_z sign flipped. */
            double _Complex hop_hole = (S * bd->J + I * S * bd->D[2]) * eikt;
            M[(i + n) * N + (j + n)] += hop_hole;
            M[(j + n) * N + (i + n)] += conj(hop_hole);
        } else {
            /* Antiparallel sublattices: HP gives anomalous pairing
             * a_i a_j + a_i^† a_j^†, plus +S·J diagonal "self-energy". */
            M[i * N + i] += S * bd->J;
            M[j * N + j] += S * bd->J;
            M[(i + n) * N + (i + n)] += S * bd->J;
            M[(j + n) * N + (j + n)] += S * bd->J;

            /* Anomalous block: B[i,j](k) and B[j,i](k) = B[i,j](-k). */
            double _Complex pair_fwd = (S * bd->J - I * S * bd->D[2]) * eikt;
            double _Complex pair_bwd = (S * bd->J - I * S * bd->D[2]) * eikt_neg;
            M[i * N + (j + n)] += pair_fwd;
            M[j * N + (i + n)] += pair_bwd;
            M[(j + n) * N + i] += conj(pair_fwd);
            M[(i + n) * N + j] += conj(pair_bwd);
        }
    }

    /* Colpa diagonalisation: M = K^† K (Cholesky), then diagonalise
     * W = K η K^† where η = diag(I, -I). Eigenvalues come in ±ω pairs;
     * the n positive ones are the magnon energies. */
    double _Complex *K = malloc((size_t)N * N * sizeof *K);
    double _Complex *W = malloc((size_t)N * N * sizeof *W);
    if (!K || !W) {
        free(M);
        free(K);
        free(W);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    /* At gap-closing points (Goldstone modes), M is positive-semi-
     * definite, not strictly positive-definite, and Cholesky fails.
     * Add a tiny regularisation to lift the zero eigenvalues into a
     * numerically tractable positive range. The introduced systematic
     * error in ω is O(√ε) ~ 1e-5, negligible compared to the LSW
     * approximation itself. */
    static const double EPS_PSD = 1e-10;
    for (int a = 0; a < N; ++a)
        M[a * N + a] += EPS_PSD;

    if (cholesky_(N, M, K) != 0) {
        free(M);
        free(K);
        free(W);
        irrep_set_error_(
            "irrep_magnon_dispersion_general: M(k) is not positive-definite even with "
            "regularisation — assumed ground state is unstable for the given bond list "
            "at k=(%g, %g)",
            kx, ky);
        return IRREP_ERR_INVALID_ARG;
    }

    /* W[i,j] = Σ_l σ_l K[i,l] conj(K[j,l]) where σ_l = +1 (l<n), -1 (l≥n). */
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double _Complex s = 0;
            for (int l = 0; l < N; ++l) {
                double sigma_l = (l < n) ? +1.0 : -1.0;
                s += sigma_l * K[i * N + l] * conj(K[j * N + l]);
            }
            W[i * N + j] = s;
        }

    /* Diagonalise W. Discard eigenvectors; keep eigenvalues. */
    double          *eigs_2n = malloc((size_t)N * sizeof *eigs_2n);
    double _Complex *evecs_2n = malloc((size_t)N * N * sizeof *evecs_2n);
    if (!eigs_2n || !evecs_2n) {
        free(M);
        free(K);
        free(W);
        free(eigs_2n);
        free(evecs_2n);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    hermitian_eig_(N, W, eigs_2n, evecs_2n);

    /* Pick positive eigenvalues, sort ascending. */
    int count = 0;
    for (int i = 0; i < N; ++i)
        if (eigs_2n[i] > 0)
            omega_out[count++] = eigs_2n[i];
    /* Should have exactly n positive eigenvalues. If not, the matrix
     * had degenerate ±ω pairs at zero — pad with zeros. */
    while (count < n)
        omega_out[count++] = 0.0;
    /* hermitian_eig_ already returns ascending; positive ones are the
     * upper half, also ascending. */

    free(M);
    free(K);
    free(W);
    free(eigs_2n);
    free(evecs_2n);
    return IRREP_OK;
}

/* Build the strip Hamiltonian H_strip(k_y) of size (Lx·n_sub)². Each
 * undirected bond contributes both forward and backward orientations.
 * Bonds with delta_x = 0 stay within a cell; |delta_x|=1 connect to the
 * next/previous cell (skip boundary wraps to enforce open BC); |delta_x|>1
 * is rejected upstream. */
static int build_H_strip_(const irrep_magnon_lsw_t *L, int Lx, double ky,
                          double _Complex *H_out) {
    int n = L->n_sub;
    int N = Lx * n;
    memset(H_out, 0, (size_t)N * N * sizeof(double _Complex));
    double S = L->S;

    /* Diagonal anisotropy on every site of every cell. */
    for (int cx = 0; cx < Lx; ++cx)
        for (int a = 0; a < n; ++a) {
            int idx = cx * n + a;
            H_out[idx * N + idx] += 2.0 * L->Kz * S;
        }

    for (int b = 0; b < L->n_bonds; ++b) {
        const irrep_magnon_bond_t *bd = &L->bonds[b];
        int dx = bd->delta_x;
        if (dx < -1 || dx > 1)
            return -1;
        /* y-component of bond translation. */
        double ty = bd->delta_x * L->a1[1] + bd->delta_y * L->a2[1];
        double phase_arg = ky * ty;
        double _Complex eikt = cos(phase_arg) + I * sin(phase_arg);

        double _Complex Jhop = +S * bd->J * eikt;
        double _Complex Dhop = -I * S * bd->D[2] * eikt;
        double _Complex hop = Jhop + Dhop;

        /* For each strip cell cx, the bond connects (cx, bi) to
         * (cx + dx, bj). If cx + dx is outside [0, Lx), drop the bond
         * (open BC). */
        for (int cx = 0; cx < Lx; ++cx) {
            int cx2 = cx + dx;
            if (cx2 < 0 || cx2 >= Lx)
                continue;
            int idx_a = cx * n + bd->bi;
            int idx_b = cx2 * n + bd->bj;
            H_out[idx_a * N + idx_b] += hop;
            H_out[idx_b * N + idx_a] += conj(hop);
            /* Diagonal "self-energy" -S·J(n_a + n_b) per unique bond.
             * Only counted on bonds that actually connect — equivalently,
             * once per undirected bond endpoint that lies inside the
             * strip. Since we already gate cx + dx ∈ [0, Lx), both
             * endpoints exist when this branch fires. */
            H_out[idx_a * N + idx_a] += -S * bd->J;
            H_out[idx_b * N + idx_b] += -S * bd->J;
        }
    }
    return 0;
}

irrep_status_t irrep_magnon_strip_dispersion(const irrep_magnon_lsw_t *L, int Lx, double ky,
                                              double *omega_out, double *edge_weight_out) {
    if (!L || Lx < 2 || !omega_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    int N = Lx * n;
    double _Complex *H = malloc((size_t)N * N * sizeof *H);
    double _Complex *V = malloc((size_t)N * N * sizeof *V);
    if (!H || !V) {
        free(H);
        free(V);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    if (build_H_strip_(L, Lx, ky, H) < 0) {
        free(H);
        free(V);
        irrep_set_error_("irrep_magnon_strip_dispersion: bonds with |delta_x| > 1 not supported");
        return IRREP_ERR_INVALID_ARG;
    }
    hermitian_eig_(N, H, omega_out, V);

    if (edge_weight_out) {
        /* For each band b, weight = Σ_{x in left half} |ψ_x|². Modes
         * localised on the left edge get values near 1; right-edge near
         * 0; bulk modes near 0.5. The "left half" is the first Lx/2
         * unit cells. */
        int half = Lx / 2;
        for (int b = 0; b < N; ++b) {
            double w = 0;
            for (int cx = 0; cx < half; ++cx)
                for (int a = 0; a < n; ++a) {
                    int idx = cx * n + a;
                    /* hermitian_eig_ writes evec b into row b,
                     * components 0..N-1 in the V buffer. */
                    double _Complex c = V[b * N + idx];
                    w += creal(c) * creal(c) + cimag(c) * cimag(c);
                }
            edge_weight_out[b] = w;
        }
    }
    free(H);
    free(V);
    return IRREP_OK;
}

double irrep_magnon_spin_nernst(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny) {
    if (!L || T <= 0 || Nx <= 0 || Ny <= 0) {
        irrep_set_error_("irrep_magnon_spin_nernst: invalid arguments");
        return NAN;
    }
    int n = L->n_sub;
    double _Complex *u_grid = malloc((size_t)Nx * Ny * n * n * sizeof *u_grid);
    double          *w_grid = malloc((size_t)Nx * Ny * n * sizeof *w_grid);
    if (!u_grid || !w_grid) {
        free(u_grid);
        free(w_grid);
        return NAN;
    }
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            int    p = iy * Nx + ix;
            irrep_magnon_dispersion(L, kx, ky, w_grid + p * n, u_grid + p * n * n);
        }

    double A_uc = fabs(L->a1[0] * L->a2[1] - L->a1[1] * L->a2[0]);
    double sum = 0;
    for (int b = 0; b < n; ++b) {
        for (int iy = 0; iy < Ny; ++iy) {
            for (int ix = 0; ix < Nx; ++ix) {
                int p00 = iy * Nx + ix;
                int p10 = iy * Nx + (ix + 1) % Nx;
                int p11 = ((iy + 1) % Ny) * Nx + (ix + 1) % Nx;
                int p01 = ((iy + 1) % Ny) * Nx + ix;
                int pairs[4][2] = {{p00, p10}, {p10, p11}, {p11, p01}, {p01, p00}};
                double _Complex z = 1.0;
                int             skipped = 0;
                for (int pp = 0; pp < 4; ++pp) {
                    double _Complex ip = inner_(u_grid + pairs[pp][0] * n * n + b * n,
                                                u_grid + pairs[pp][1] * n * n + b * n, n);
                    double          mod = cabs(ip);
                    if (mod < 1e-300) {
                        skipped = 1;
                        break;
                    }
                    z *= ip / mod;
                }
                if (skipped)
                    continue;
                double F_p = atan2(cimag(z), creal(z));
                double w_b = w_grid[p00 * n + b];
                double w_eff = w_b > 1e-12 ? w_b : 1e-12;
                double g = 1.0 / (exp(w_eff / T) - 1.0);
                sum += c1_mm_(g) * F_p;
            }
        }
    }
    free(u_grid);
    free(w_grid);
    /* Sign convention: α^s_xy = -(1/(A_uc · (2π)²)) · Σ_b Σ_p c₁(g) · F_b. */
    return -sum / (A_uc * 4.0 * M_PI * M_PI);
}

double irrep_magnon_thermal_hall_kxy(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny) {
    if (!L || T <= 0 || Nx <= 0 || Ny <= 0) {
        irrep_set_error_("irrep_magnon_thermal_hall_kxy: invalid arguments");
        return NAN;
    }
    int n = L->n_sub;
    double _Complex *u_grid = malloc((size_t)Nx * Ny * n * n * sizeof *u_grid);
    double          *w_grid = malloc((size_t)Nx * Ny * n * sizeof *w_grid);
    if (!u_grid || !w_grid) {
        free(u_grid);
        free(w_grid);
        return NAN;
    }
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            int    p = iy * Nx + ix;
            irrep_magnon_dispersion(L, kx, ky, w_grid + p * n, u_grid + p * n * n);
        }

    /* 2D unit cell area = |a1 × a2|. */
    double A_uc = fabs(L->a1[0] * L->a2[1] - L->a1[1] * L->a2[0]);

    double sum = 0;
    for (int b = 0; b < n; ++b) {
        for (int iy = 0; iy < Ny; ++iy) {
            for (int ix = 0; ix < Nx; ++ix) {
                int p00 = iy * Nx + ix;
                int p10 = iy * Nx + (ix + 1) % Nx;
                int p11 = ((iy + 1) % Ny) * Nx + (ix + 1) % Nx;
                int p01 = ((iy + 1) % Ny) * Nx + ix;
                int pairs[4][2] = {{p00, p10}, {p10, p11}, {p11, p01}, {p01, p00}};
                double _Complex z = 1.0;
                int             skipped = 0;
                for (int pp = 0; pp < 4; ++pp) {
                    double _Complex ip = inner_(u_grid + pairs[pp][0] * n * n + b * n,
                                                u_grid + pairs[pp][1] * n * n + b * n, n);
                    double          mod = cabs(ip);
                    if (mod < 1e-300) {
                        skipped = 1;
                        break;
                    }
                    z *= ip / mod;
                }
                if (skipped)
                    continue;
                double F_p = atan2(cimag(z), creal(z));
                double w_b = w_grid[p00 * n + b];
                /* Below the Goldstone gap (negative or zero ω due to FP
                 * noise) the BE distribution diverges; clamp ω to a
                 * small positive ε so the prefactor stays finite. */
                double w_eff = w_b > 1e-12 ? w_b : 1e-12;
                double g = 1.0 / (exp(w_eff / T) - 1.0);
                sum += c2_mm_(g) * F_p;
            }
        }
    }

    free(u_grid);
    free(w_grid);
    /* κ_xy = -(T / (A_uc · (2π)²)) · Σ_b Σ_p c₂(g_b(k_p)) · F_b(p) */
    return -T * sum / (A_uc * 4.0 * M_PI * M_PI);
}
