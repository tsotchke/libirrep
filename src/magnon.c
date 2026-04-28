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

irrep_status_t irrep_magnon_softest_mode(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                          double exclude_below, double *kx_out, double *ky_out,
                                          double *omega_out, int *band_out) {
    if (!L || Nx <= 0 || Ny <= 0 || !kx_out || !ky_out || !omega_out || !band_out)
        return IRREP_ERR_INVALID_ARG;
    int     n = L->n_sub;
    double *omega = malloc((size_t)n * sizeof *omega);
    double _Complex *u = malloc((size_t)n * n * sizeof *u);
    if (!omega || !u) {
        free(omega);
        free(u);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    double w_best = 1e300;
    double kx_best = 0, ky_best = 0;
    int    b_best = 0;
    int    found = 0;
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            irrep_magnon_dispersion(L, kx, ky, omega, u);
            for (int b = 0; b < n; ++b) {
                double w = omega[b];
                if (w > exclude_below && w < w_best) {
                    w_best = w;
                    kx_best = kx;
                    ky_best = ky;
                    b_best = b;
                    found = 1;
                }
            }
        }
    free(omega);
    free(u);
    if (!found) {
        irrep_set_error_("irrep_magnon_softest_mode: no modes above exclude_below");
        return IRREP_ERR_INVALID_ARG;
    }
    *kx_out = kx_best;
    *ky_out = ky_best;
    *omega_out = w_best;
    *band_out = b_best;
    return IRREP_OK;
}

irrep_status_t irrep_magnon_band_extrema(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                          double exclude_below, double *omega_min_out,
                                          double *omega_max_out) {
    if (!L || Nx <= 0 || Ny <= 0 || !omega_min_out || !omega_max_out)
        return IRREP_ERR_INVALID_ARG;
    int     n = L->n_sub;
    double *omega = malloc((size_t)n * sizeof *omega);
    double _Complex *u = malloc((size_t)n * n * sizeof *u);
    if (!omega || !u) {
        free(omega);
        free(u);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    double w_min = 1e300;
    double w_max = -1e300;
    int    found_any_above_threshold = 0;
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            irrep_magnon_dispersion(L, kx, ky, omega, u);
            for (int b = 0; b < n; ++b) {
                double w = omega[b];
                if (w > exclude_below) {
                    if (w < w_min)
                        w_min = w;
                    if (w > w_max)
                        w_max = w;
                    found_any_above_threshold = 1;
                }
            }
        }
    free(omega);
    free(u);
    if (!found_any_above_threshold) {
        irrep_set_error_("irrep_magnon_band_extrema: no modes above exclude_below threshold");
        return IRREP_ERR_INVALID_ARG;
    }
    *omega_min_out = w_min;
    *omega_max_out = w_max;
    return IRREP_OK;
}

irrep_status_t irrep_magnon_hessian(const irrep_magnon_lsw_t *L, double kx, double ky, double h,
                                     double *hxx_out, double *hyy_out, double *hxy_out) {
    if (!L || !hxx_out || !hyy_out || !hxy_out || h <= 0)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    double *w_c = malloc((size_t)n * sizeof *w_c);   /* ω at (kx, ky) */
    double *w_xp = malloc((size_t)n * sizeof *w_xp); /* (kx+h, ky)    */
    double *w_xm = malloc((size_t)n * sizeof *w_xm);
    double *w_yp = malloc((size_t)n * sizeof *w_yp);
    double *w_ym = malloc((size_t)n * sizeof *w_ym);
    double *w_pp = malloc((size_t)n * sizeof *w_pp); /* (kx+h, ky+h)  */
    double *w_pm = malloc((size_t)n * sizeof *w_pm); /* (kx+h, ky-h)  */
    double *w_mp = malloc((size_t)n * sizeof *w_mp); /* (kx-h, ky+h)  */
    double *w_mm = malloc((size_t)n * sizeof *w_mm); /* (kx-h, ky-h)  */
    double _Complex *u_dummy = malloc((size_t)n * n * sizeof *u_dummy);
    if (!w_c || !w_xp || !w_xm || !w_yp || !w_ym || !w_pp || !w_pm || !w_mp || !w_mm ||
        !u_dummy) {
        free(w_c); free(w_xp); free(w_xm); free(w_yp); free(w_ym);
        free(w_pp); free(w_pm); free(w_mp); free(w_mm); free(u_dummy);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    irrep_magnon_dispersion(L, kx, ky, w_c, u_dummy);
    irrep_magnon_dispersion(L, kx + h, ky, w_xp, u_dummy);
    irrep_magnon_dispersion(L, kx - h, ky, w_xm, u_dummy);
    irrep_magnon_dispersion(L, kx, ky + h, w_yp, u_dummy);
    irrep_magnon_dispersion(L, kx, ky - h, w_ym, u_dummy);
    irrep_magnon_dispersion(L, kx + h, ky + h, w_pp, u_dummy);
    irrep_magnon_dispersion(L, kx + h, ky - h, w_pm, u_dummy);
    irrep_magnon_dispersion(L, kx - h, ky + h, w_mp, u_dummy);
    irrep_magnon_dispersion(L, kx - h, ky - h, w_mm, u_dummy);

    double h2 = h * h;
    double four_h2 = 4.0 * h2;
    for (int b = 0; b < n; ++b) {
        hxx_out[b] = (w_xp[b] - 2.0 * w_c[b] + w_xm[b]) / h2;
        hyy_out[b] = (w_yp[b] - 2.0 * w_c[b] + w_ym[b]) / h2;
        hxy_out[b] = (w_pp[b] - w_pm[b] - w_mp[b] + w_mm[b]) / four_h2;
    }

    free(w_c); free(w_xp); free(w_xm); free(w_yp); free(w_ym);
    free(w_pp); free(w_pm); free(w_mp); free(w_mm); free(u_dummy);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_two_magnon_qomega_general(const irrep_magnon_lsw_t *L,
                                                        const int *sublattice_signs,
                                                        const double (*qpath)[2], int n_q,
                                                        int Nx, int Ny, double omega_min,
                                                        double omega_max, int n_omega,
                                                        double eta, double *intensity_out) {
    if (!L || !sublattice_signs || !qpath || n_q <= 0 || Nx <= 0 || Ny <= 0 || n_omega <= 0 ||
        eta <= 0 || omega_max <= omega_min || !intensity_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(intensity_out, 0, (size_t)n_q * (size_t)n_omega * sizeof *intensity_out);
    double dw = (omega_max - omega_min) / n_omega;
    double inv_NBZ = 1.0 / ((double)Nx * (double)Ny);

    int N_grid = Nx * Ny;
    double *omega_grid = malloc((size_t)N_grid * n * sizeof *omega_grid);
    double *S_perp_grid = malloc((size_t)N_grid * n * sizeof *S_perp_grid);
    if (!omega_grid || !S_perp_grid) {
        free(omega_grid);
        free(S_perp_grid);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    /* Use AFM-aware Bogoliubov structure factor on the BZ grid. */
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = ((double)ix + 0.5) / Nx; /* half-shift to avoid Goldstone */
            double fy = ((double)iy + 0.5) / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            int    p = iy * Nx + ix;
            irrep_magnon_structure_factor_general(L, sublattice_signs, kx, ky,
                                                    omega_grid + p * n, S_perp_grid + p * n);
        }

    for (int iq = 0; iq < n_q; ++iq) {
        double qx = qpath[iq][0];
        double qy = qpath[iq][1];
        for (int p1 = 0; p1 < N_grid; ++p1) {
            int    ix1 = p1 % Nx;
            int    iy1 = p1 / Nx;
            double fx1 = ((double)ix1 + 0.5) / Nx;
            double fy1 = ((double)iy1 + 0.5) / Ny;
            double kx1 = fx1 * L->b1[0] + fy1 * L->b2[0];
            double ky1 = fx1 * L->b1[1] + fy1 * L->b2[1];
            double k2x = qx - kx1;
            double k2y = qy - ky1;
            double det = L->b1[0] * L->b2[1] - L->b1[1] * L->b2[0];
            double f2x = (k2x * L->b2[1] - k2y * L->b2[0]) / det;
            double f2y = (-k2x * L->b1[1] + k2y * L->b1[0]) / det;
            f2x = f2x - floor(f2x);
            f2y = f2y - floor(f2y);
            int ix2 = (int)(f2x * Nx + 0.5) % Nx;
            int iy2 = (int)(f2y * Ny + 0.5) % Ny;
            int p2 = iy2 * Nx + ix2;

            for (int b1 = 0; b1 < n; ++b1)
                for (int b2 = 0; b2 < n; ++b2) {
                    double w_total = omega_grid[p1 * n + b1] + omega_grid[p2 * n + b2];
                    double weight = S_perp_grid[p1 * n + b1] * S_perp_grid[p2 * n + b2];
                    for (int jw = 0; jw < n_omega; ++jw) {
                        double w = omega_min + (jw + 0.5) * dw;
                        double dx = w - w_total;
                        intensity_out[iq * n_omega + jw] +=
                            inv_NBZ * weight * (eta / M_PI) / (dx * dx + eta * eta);
                    }
                }
        }
    }

    free(omega_grid);
    free(S_perp_grid);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_one_magnon_qomega_general(const irrep_magnon_lsw_t *L,
                                                        const int *sublattice_signs,
                                                        const double (*qpath)[2], int n_q,
                                                        double omega_min, double omega_max,
                                                        int n_omega, double eta,
                                                        double *intensity_out) {
    if (!L || !sublattice_signs || !qpath || n_q <= 0 || n_omega <= 0 || eta <= 0 ||
        omega_max <= omega_min || !intensity_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(intensity_out, 0, (size_t)n_q * (size_t)n_omega * sizeof *intensity_out);
    double dw = (omega_max - omega_min) / n_omega;
    double *omega = malloc((size_t)n * sizeof *omega);
    double *Sb    = malloc((size_t)n * sizeof *Sb);
    if (!omega || !Sb) {
        free(omega);
        free(Sb);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iq = 0; iq < n_q; ++iq) {
        irrep_magnon_structure_factor_general(L, sublattice_signs, qpath[iq][0], qpath[iq][1],
                                                omega, Sb);
        for (int b = 0; b < n; ++b) {
            for (int jw = 0; jw < n_omega; ++jw) {
                double w  = omega_min + (jw + 0.5) * dw;
                double dx = w - omega[b];
                intensity_out[iq * n_omega + jw] +=
                    Sb[b] * (eta / M_PI) / (dx * dx + eta * eta);
            }
        }
    }
    free(omega);
    free(Sb);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_dynamical_structure_factor_T(const irrep_magnon_lsw_t *L,
                                                           const double (*qpath)[2], int n_q,
                                                           double omega_min, double omega_max,
                                                           int n_omega, double eta, double T,
                                                           double *intensity_out) {
    if (!L || !qpath || n_q <= 0 || n_omega <= 0 || eta <= 0 || omega_max <= omega_min ||
        T < 0 || !intensity_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(intensity_out, 0, (size_t)n_q * (size_t)n_omega * sizeof *intensity_out);
    double dw = (omega_max - omega_min) / n_omega;
    double *omega = malloc((size_t)n * sizeof *omega);
    double *Sb    = malloc((size_t)n * sizeof *Sb);
    if (!omega || !Sb) {
        free(omega);
        free(Sb);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iq = 0; iq < n_q; ++iq) {
        irrep_magnon_structure_factor(L, qpath[iq][0], qpath[iq][1], omega, Sb);
        for (int b = 0; b < n; ++b) {
            /* (1 + n_B) Stokes factor. Clamp tiny ω to avoid 0/0; at
             * T=0 n_B=0 → factor 1 (recovers the T=0 path). */
            double bose1 = 1.0;
            if (T > 0 && omega[b] > 1e-12) {
                double x = omega[b] / T;
                bose1 = (x > 50) ? 1.0 : 1.0 / (1.0 - exp(-x));
            }
            double w_b = bose1 * Sb[b];
            for (int jw = 0; jw < n_omega; ++jw) {
                double w  = omega_min + (jw + 0.5) * dw;
                double dx = w - omega[b];
                intensity_out[iq * n_omega + jw] +=
                    w_b * (eta / M_PI) / (dx * dx + eta * eta);
            }
        }
    }
    free(omega);
    free(Sb);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_dynamical_structure_factor_T_general(
    const irrep_magnon_lsw_t *L, const int *sublattice_signs, const double (*qpath)[2], int n_q,
    double omega_min, double omega_max, int n_omega, double eta, double T,
    double *intensity_out) {
    if (!L || !sublattice_signs || !qpath || n_q <= 0 || n_omega <= 0 || eta <= 0 ||
        omega_max <= omega_min || T < 0 || !intensity_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(intensity_out, 0, (size_t)n_q * (size_t)n_omega * sizeof *intensity_out);
    double dw = (omega_max - omega_min) / n_omega;
    double *omega = malloc((size_t)n * sizeof *omega);
    double *Sb    = malloc((size_t)n * sizeof *Sb);
    if (!omega || !Sb) {
        free(omega);
        free(Sb);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iq = 0; iq < n_q; ++iq) {
        irrep_magnon_structure_factor_general(L, sublattice_signs, qpath[iq][0], qpath[iq][1],
                                                omega, Sb);
        for (int b = 0; b < n; ++b) {
            double bose1 = 1.0;
            if (T > 0 && omega[b] > 1e-12) {
                double x = omega[b] / T;
                bose1 = (x > 50) ? 1.0 : 1.0 / (1.0 - exp(-x));
            }
            double w_b = bose1 * Sb[b];
            for (int jw = 0; jw < n_omega; ++jw) {
                double w  = omega_min + (jw + 0.5) * dw;
                double dx = w - omega[b];
                intensity_out[iq * n_omega + jw] +=
                    w_b * (eta / M_PI) / (dx * dx + eta * eta);
            }
        }
    }
    free(omega);
    free(Sb);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_dynamical_structure_factor_T_anti_stokes(
    const irrep_magnon_lsw_t *L, const double (*qpath)[2], int n_q, double omega_min,
    double omega_max, int n_omega, double eta, double T, double *intensity_out) {
    if (!L || !qpath || n_q <= 0 || n_omega <= 0 || eta <= 0 || omega_max <= omega_min ||
        T < 0 || !intensity_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(intensity_out, 0, (size_t)n_q * (size_t)n_omega * sizeof *intensity_out);
    if (T <= 0) return IRREP_OK; /* no thermal magnons → anti-Stokes vanishes */
    double dw = (omega_max - omega_min) / n_omega;
    double *omega = malloc((size_t)n * sizeof *omega);
    double *Sb    = malloc((size_t)n * sizeof *Sb);
    if (!omega || !Sb) {
        free(omega);
        free(Sb);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iq = 0; iq < n_q; ++iq) {
        irrep_magnon_structure_factor(L, qpath[iq][0], qpath[iq][1], omega, Sb);
        for (int b = 0; b < n; ++b) {
            if (omega[b] < 1e-12) continue; /* Goldstone — n_B undefined; skip */
            double x = omega[b] / T;
            double bose = (x > 50) ? 0.0 : 1.0 / (exp(x) - 1.0);
            double w_b = bose * Sb[b];
            for (int jw = 0; jw < n_omega; ++jw) {
                double w  = omega_min + (jw + 0.5) * dw;
                double dx = w + omega[b]; /* anti-Stokes peak at -ω_b */
                intensity_out[iq * n_omega + jw] +=
                    w_b * (eta / M_PI) / (dx * dx + eta * eta);
            }
        }
    }
    free(omega);
    free(Sb);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_dynamical_structure_factor_T_anti_stokes_general(
    const irrep_magnon_lsw_t *L, const int *sublattice_signs, const double (*qpath)[2],
    int n_q, double omega_min, double omega_max, int n_omega, double eta, double T,
    double *intensity_out) {
    if (!L || !sublattice_signs || !qpath || n_q <= 0 || n_omega <= 0 || eta <= 0 ||
        omega_max <= omega_min || T < 0 || !intensity_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(intensity_out, 0, (size_t)n_q * (size_t)n_omega * sizeof *intensity_out);
    if (T <= 0) return IRREP_OK;
    double dw = (omega_max - omega_min) / n_omega;
    double *omega = malloc((size_t)n * sizeof *omega);
    double *Sb    = malloc((size_t)n * sizeof *Sb);
    if (!omega || !Sb) {
        free(omega);
        free(Sb);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iq = 0; iq < n_q; ++iq) {
        irrep_magnon_structure_factor_general(L, sublattice_signs, qpath[iq][0], qpath[iq][1],
                                                omega, Sb);
        for (int b = 0; b < n; ++b) {
            if (omega[b] < 1e-12) continue;
            double x = omega[b] / T;
            double bose = (x > 50) ? 0.0 : 1.0 / (exp(x) - 1.0);
            double w_b = bose * Sb[b];
            for (int jw = 0; jw < n_omega; ++jw) {
                double w  = omega_min + (jw + 0.5) * dw;
                double dx = w + omega[b];
                intensity_out[iq * n_omega + jw] +=
                    w_b * (eta / M_PI) / (dx * dx + eta * eta);
            }
        }
    }
    free(omega);
    free(Sb);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_dynamical_structure_factor(const irrep_magnon_lsw_t *L,
                                                         const double (*qpath)[2], int n_q,
                                                         int Nx, int Ny, double omega_min,
                                                         double omega_max, int n_omega,
                                                         double eta, double *intensity_out) {
    if (!intensity_out) return IRREP_ERR_INVALID_ARG;
    irrep_status_t st = irrep_magnon_neutron_qomega_map(L, qpath, n_q, omega_min, omega_max,
                                                          n_omega, eta, intensity_out);
    if (st != IRREP_OK) return st;
    double *I2 = malloc((size_t)n_q * (size_t)n_omega * sizeof *I2);
    if (!I2) return IRREP_ERR_OUT_OF_MEMORY;
    st = irrep_magnon_two_magnon_qomega(L, qpath, n_q, Nx, Ny, omega_min, omega_max, n_omega,
                                          eta, I2);
    if (st != IRREP_OK) {
        free(I2);
        return st;
    }
    for (size_t i = 0; i < (size_t)n_q * (size_t)n_omega; ++i) intensity_out[i] += I2[i];
    free(I2);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_dynamical_structure_factor_general(
    const irrep_magnon_lsw_t *L, const int *sublattice_signs, const double (*qpath)[2],
    int n_q, int Nx, int Ny, double omega_min, double omega_max, int n_omega, double eta,
    double *intensity_out) {
    if (!intensity_out) return IRREP_ERR_INVALID_ARG;
    irrep_status_t st = irrep_magnon_one_magnon_qomega_general(
        L, sublattice_signs, qpath, n_q, omega_min, omega_max, n_omega, eta, intensity_out);
    if (st != IRREP_OK) return st;
    double *I2 = malloc((size_t)n_q * (size_t)n_omega * sizeof *I2);
    if (!I2) return IRREP_ERR_OUT_OF_MEMORY;
    st = irrep_magnon_two_magnon_qomega_general(L, sublattice_signs, qpath, n_q, Nx, Ny,
                                                  omega_min, omega_max, n_omega, eta, I2);
    if (st != IRREP_OK) {
        free(I2);
        return st;
    }
    for (size_t i = 0; i < (size_t)n_q * (size_t)n_omega; ++i) intensity_out[i] += I2[i];
    free(I2);
    return IRREP_OK;
}

/* Internal helper: build LSW handle, evaluate χ² over the observation
 * set, free handle. Returns +∞ if construction fails. */
static double fit_chi2_(int n_sub, double S, const double a1[2], const double a2[2],
                         irrep_magnon_bond_t *bonds, int n_bonds, double Kz,
                         const double (*q_obs)[2], const double *omega_obs,
                         const int *band_obs, int n_obs, double *omega_buf,
                         double _Complex *u_buf) {
    irrep_magnon_lsw_t *Lh = irrep_magnon_lsw_new(n_sub, S, a1, a2, n_bonds, bonds, Kz);
    if (!Lh) return INFINITY;
    double chi2 = 0;
    for (int i = 0; i < n_obs; ++i) {
        if (band_obs[i] < 0 || band_obs[i] >= n_sub) continue;
        irrep_magnon_dispersion(Lh, q_obs[i][0], q_obs[i][1], omega_buf, u_buf);
        double r = omega_buf[band_obs[i]] - omega_obs[i];
        chi2 += r * r;
    }
    irrep_magnon_lsw_free(Lh);
    return chi2;
}

irrep_status_t irrep_magnon_fit_J(int n_sub, double S, const double a1[2], const double a2[2],
                                    irrep_magnon_bond_t *bonds, int n_bonds, double Kz,
                                    const double (*q_obs)[2], const double *omega_obs,
                                    const int *band_obs, int n_obs, int max_iter, double tol,
                                    double *chi2_out) {
    if (!bonds || !q_obs || !omega_obs || !band_obs || n_bonds <= 0 || n_obs <= 0 ||
        max_iter <= 0 || tol < 0)
        return IRREP_ERR_INVALID_ARG;

    double *omega = malloc((size_t)n_sub * sizeof *omega);
    double _Complex *u = malloc((size_t)n_sub * n_sub * sizeof *u);
    double *grad = malloc((size_t)n_bonds * sizeof *grad);
    double *J_save = malloc((size_t)n_bonds * sizeof *J_save);
    if (!omega || !u || !grad || !J_save) {
        free(omega);
        free(u);
        free(grad);
        free(J_save);
        return IRREP_ERR_OUT_OF_MEMORY;
    }

    double chi2_curr = fit_chi2_(n_sub, S, a1, a2, bonds, n_bonds, Kz, q_obs, omega_obs,
                                  band_obs, n_obs, omega, u);
    if (!isfinite(chi2_curr)) {
        free(omega);
        free(u);
        free(grad);
        free(J_save);
        return IRREP_ERR_INVALID_ARG;
    }
    double eps   = 1e-4;
    double alpha = 0.01;

    for (int iter = 0; iter < max_iter; ++iter) {
        /* Save current J's. */
        for (int k = 0; k < n_bonds; ++k) J_save[k] = bonds[k].J;

        /* Forward-difference gradient. */
        for (int k = 0; k < n_bonds; ++k) {
            bonds[k].J = J_save[k] + eps;
            double chi2_plus = fit_chi2_(n_sub, S, a1, a2, bonds, n_bonds, Kz, q_obs,
                                          omega_obs, band_obs, n_obs, omega, u);
            bonds[k].J = J_save[k];
            grad[k]    = (chi2_plus - chi2_curr) / eps;
        }

        /* Backtracking line search. */
        double chi2_new = chi2_curr;
        int    accepted = 0;
        for (int bt = 0; bt < 30; ++bt) {
            for (int k = 0; k < n_bonds; ++k) bonds[k].J = J_save[k] - alpha * grad[k];
            chi2_new = fit_chi2_(n_sub, S, a1, a2, bonds, n_bonds, Kz, q_obs, omega_obs,
                                  band_obs, n_obs, omega, u);
            if (isfinite(chi2_new) && chi2_new < chi2_curr) {
                alpha *= 1.1;
                accepted = 1;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            /* Stuck — restore and stop. */
            for (int k = 0; k < n_bonds; ++k) bonds[k].J = J_save[k];
            break;
        }
        if (fabs(chi2_curr - chi2_new) < tol) {
            chi2_curr = chi2_new;
            break;
        }
        chi2_curr = chi2_new;
    }

    if (chi2_out) *chi2_out = chi2_curr;
    free(omega);
    free(u);
    free(grad);
    free(J_save);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_kinematic_damping(const irrep_magnon_lsw_t *L,
                                                const double (*kpath)[2], int n_k, int Nx,
                                                int Ny, double eta, double *gamma_out) {
    if (!L || !kpath || n_k <= 0 || Nx <= 0 || Ny <= 0 || eta <= 0 || !gamma_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(gamma_out, 0, (size_t)n_k * (size_t)n * sizeof *gamma_out);

    int N_grid = Nx * Ny;
    double *omega_grid = malloc((size_t)N_grid * n * sizeof *omega_grid);
    double *omega_k    = malloc((size_t)n * sizeof *omega_k);
    double _Complex *u_dummy = malloc((size_t)n * n * sizeof *u_dummy);
    if (!omega_grid || !omega_k || !u_dummy) {
        free(omega_grid);
        free(omega_k);
        free(u_dummy);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    /* Pre-compute 1-magnon dispersion on the BZ grid (half-shift to
     * stay clear of Goldstone). */
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = ((double)ix + 0.5) / Nx;
            double fy = ((double)iy + 0.5) / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            irrep_magnon_dispersion(L, kx, ky, omega_grid + (iy * Nx + ix) * n, u_dummy);
        }
    double inv_NBZ = 1.0 / ((double)Nx * (double)Ny);
    double det = L->b1[0] * L->b2[1] - L->b1[1] * L->b2[0];

    for (int ik = 0; ik < n_k; ++ik) {
        irrep_magnon_dispersion(L, kpath[ik][0], kpath[ik][1], omega_k, u_dummy);
        for (int b = 0; b < n; ++b) {
            double w_target = omega_k[b];
            double accum    = 0.0;
            for (int p1 = 0; p1 < N_grid; ++p1) {
                int    ix1 = p1 % Nx;
                int    iy1 = p1 / Nx;
                double fx1 = ((double)ix1 + 0.5) / Nx;
                double fy1 = ((double)iy1 + 0.5) / Ny;
                double kx1 = fx1 * L->b1[0] + fy1 * L->b2[0];
                double ky1 = fx1 * L->b1[1] + fy1 * L->b2[1];
                /* k₂ = k - k₁; nearest-neighbour grid lookup. */
                double k2x = kpath[ik][0] - kx1;
                double k2y = kpath[ik][1] - ky1;
                double f2x = (k2x * L->b2[1] - k2y * L->b2[0]) / det;
                double f2y = (-k2x * L->b1[1] + k2y * L->b1[0]) / det;
                f2x        = f2x - floor(f2x);
                f2y        = f2y - floor(f2y);
                int ix2    = (int)(f2x * Nx + 0.5) % Nx;
                int iy2    = (int)(f2y * Ny + 0.5) % Ny;
                int p2     = iy2 * Nx + ix2;
                for (int b1 = 0; b1 < n; ++b1)
                    for (int b2 = 0; b2 < n; ++b2) {
                        double dx = w_target - omega_grid[p1 * n + b1] -
                                    omega_grid[p2 * n + b2];
                        accum += inv_NBZ * (eta / M_PI) / (dx * dx + eta * eta);
                    }
            }
            /* Γ_kin = π · D^(2)(k, ω_b(k)). */
            gamma_out[ik * n + b] = M_PI * accum;
        }
    }
    free(omega_grid);
    free(omega_k);
    free(u_dummy);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_two_magnon_qomega(const irrep_magnon_lsw_t *L,
                                                const double (*qpath)[2], int n_q, int Nx,
                                                int Ny, double omega_min, double omega_max,
                                                int n_omega, double eta, double *intensity_out) {
    if (!L || !qpath || n_q <= 0 || Nx <= 0 || Ny <= 0 || n_omega <= 0 || eta <= 0 ||
        omega_max <= omega_min || !intensity_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(intensity_out, 0, (size_t)n_q * (size_t)n_omega * sizeof *intensity_out);
    double dw = (omega_max - omega_min) / n_omega;
    double inv_NBZ = 1.0 / ((double)Nx * (double)Ny);

    /* Pre-compute 1-magnon dispersion + structure factor on BZ grid. */
    int N_grid = Nx * Ny;
    double *omega_grid = malloc((size_t)N_grid * n * sizeof *omega_grid);
    double *S_perp_grid = malloc((size_t)N_grid * n * sizeof *S_perp_grid);
    if (!omega_grid || !S_perp_grid) {
        free(omega_grid);
        free(S_perp_grid);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            int    p = iy * Nx + ix;
            irrep_magnon_structure_factor(L, kx, ky, omega_grid + p * n, S_perp_grid + p * n);
        }

    /* For each q on the path, sum over k, b₁, b₂. The (q-k) momentum
     * needs interpolation from the grid; we use nearest-neighbour. */
    for (int iq = 0; iq < n_q; ++iq) {
        double qx = qpath[iq][0];
        double qy = qpath[iq][1];
        for (int p1 = 0; p1 < N_grid; ++p1) {
            int    ix1 = p1 % Nx;
            int    iy1 = p1 / Nx;
            double fx1 = (double)ix1 / Nx;
            double fy1 = (double)iy1 / Ny;
            double kx1 = fx1 * L->b1[0] + fy1 * L->b2[0];
            double ky1 = fx1 * L->b1[1] + fy1 * L->b2[1];
            /* Compute (q - k₁) and find nearest grid point. */
            double k2x = qx - kx1;
            double k2y = qy - ky1;
            /* Project onto reciprocal-vector basis to get fractional grid coords.
             * Use the formula (b1, b2)·(fx, fy) = (kx, ky), invert: */
            double det = L->b1[0] * L->b2[1] - L->b1[1] * L->b2[0];
            double f2x = (k2x * L->b2[1] - k2y * L->b2[0]) / det;
            double f2y = (-k2x * L->b1[1] + k2y * L->b1[0]) / det;
            /* Wrap to [0, 1) */
            f2x = f2x - floor(f2x);
            f2y = f2y - floor(f2y);
            int ix2 = (int)(f2x * Nx + 0.5) % Nx;
            int iy2 = (int)(f2y * Ny + 0.5) % Ny;
            int p2 = iy2 * Nx + ix2;

            for (int b1 = 0; b1 < n; ++b1)
                for (int b2 = 0; b2 < n; ++b2) {
                    double w_total =
                        omega_grid[p1 * n + b1] + omega_grid[p2 * n + b2];
                    double weight = S_perp_grid[p1 * n + b1] * S_perp_grid[p2 * n + b2];
                    /* Convolve with Lorentzian over the energy axis. */
                    for (int jw = 0; jw < n_omega; ++jw) {
                        double w = omega_min + (jw + 0.5) * dw;
                        double dx = w - w_total;
                        intensity_out[iq * n_omega + jw] +=
                            inv_NBZ * weight * (eta / M_PI) / (dx * dx + eta * eta);
                    }
                }
        }
    }

    free(omega_grid);
    free(S_perp_grid);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_two_magnon_dos(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                            double omega_min, double omega_max, int n_bins,
                                            double *dos_out) {
    if (!L || Nx <= 0 || Ny <= 0 || n_bins <= 0 || omega_max <= omega_min || !dos_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(dos_out, 0, (size_t)n_bins * sizeof *dos_out);
    double bin_w = (omega_max - omega_min) / n_bins;

    /* Pre-compute 1-magnon dispersion on the BZ grid. Avoids redundant
     * diagonalisation inside the (k₁, k₂) double loop. */
    int N_grid = Nx * Ny;
    double *omega_grid = malloc((size_t)N_grid * n * sizeof *omega_grid);
    double _Complex *u_grid = malloc((size_t)N_grid * n * n * sizeof *u_grid);
    if (!omega_grid || !u_grid) {
        free(omega_grid);
        free(u_grid);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = (double)ix / Nx;
            double fy = (double)iy / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            int    p = iy * Nx + ix;
            irrep_magnon_dispersion(L, kx, ky, omega_grid + p * n, u_grid + p * n * n);
        }

    /* Convolve: for each pair (p1, p2) and each (b1, b2), histogram
     * ω(k₁, b₁) + ω(k₂, b₂) into the 2-magnon DOS. */
    for (int p1 = 0; p1 < N_grid; ++p1)
        for (int p2 = 0; p2 < N_grid; ++p2)
            for (int b1 = 0; b1 < n; ++b1)
                for (int b2 = 0; b2 < n; ++b2) {
                    double w_sum = omega_grid[p1 * n + b1] + omega_grid[p2 * n + b2];
                    int    idx = (int)((w_sum - omega_min) / bin_w);
                    if (idx >= 0 && idx < n_bins)
                        dos_out[idx] += 1.0;
                }

    /* Normalise so ∫ D⁽²⁾(ω) dω = n_sub². Bin sum is N_grid² · n²;
     * after dividing by N_grid² and bin_w, integral = n²·bin_w·n_bins
     * (if all bins occupied) — actually we need ∫ D dω = n_sub², so:
     * D bin = (count) / (N_grid² · bin_w) gives ∫ = (total count) / (N_grid² · bin_w) · bin_w
     *       = total count / N_grid² = n². */
    double norm = 1.0 / ((double)N_grid * (double)N_grid * bin_w);
    for (int i = 0; i < n_bins; ++i)
        dos_out[i] *= norm;

    free(omega_grid);
    free(u_grid);
    return IRREP_OK;
}

double irrep_magnon_hartree_renormalisation(const irrep_magnon_lsw_t *L, double T, int Nx,
                                              int Ny) {
    if (!L || T <= 0 || Nx <= 0 || Ny <= 0) {
        irrep_set_error_("irrep_magnon_hartree_renormalisation: invalid arguments");
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
    /* ⟨n⟩(T) = (1/N_BZ) · Σ_b Σ_k n_BE(ω_LSW(k, b)/T)
     * (per cell — sum over bands AND momentum, divided by N_BZ only). */
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
                if (x > 700)
                    continue;
                accum += 1.0 / (exp(x) - 1.0);
            }
        }
    free(omega);
    free(u);
    /* Per cell, normalised to N_BZ; <n>/S is the dimensionless
     * "thermal magnon density relative to the maximum (S)". */
    double n_avg = accum / ((double)Nx * (double)Ny);
    double Z = 1.0 - n_avg / L->S;
    /* Clamp to [0, 1]: at extreme temperatures the formula gives
     * Z < 0 which is unphysical (LSW has broken down completely). */
    if (Z < 0) Z = 0;
    if (Z > 1) Z = 1;
    return Z;
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

double irrep_magnon_internal_energy(const irrep_magnon_lsw_t *L, double T, int Nx, int Ny) {
    if (!L || T <= 0 || Nx <= 0 || Ny <= 0) {
        irrep_set_error_("irrep_magnon_internal_energy: invalid arguments");
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
                accum += w * n_BE;
            }
        }
    free(omega);
    free(u);
    /* U(T) = (1/N_BZ) · Σ_b Σ_k ω_b · n_B(ω_b/T) per cell. */
    return accum / ((double)Nx * (double)Ny);
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

irrep_status_t irrep_magnon_powder_spectrum(const irrep_magnon_lsw_t *L, int Nx, int Ny,
                                              double omega_min, double omega_max, int n_bins,
                                              double *S_out) {
    if (!L || Nx <= 0 || Ny <= 0 || n_bins <= 0 || omega_max <= omega_min || !S_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(S_out, 0, (size_t)n_bins * sizeof *S_out);
    double  bin_w = (omega_max - omega_min) / n_bins;
    double *omega = malloc((size_t)n * sizeof *omega);
    double *Sb    = malloc((size_t)n * sizeof *Sb);
    if (!omega || !Sb) {
        free(omega);
        free(Sb);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = ((double)ix + 0.5) / Nx;
            double fy = ((double)iy + 0.5) / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            irrep_magnon_structure_factor(L, kx, ky, omega, Sb);
            for (int b = 0; b < n; ++b) {
                int idx = (int)((omega[b] - omega_min) / bin_w);
                if (idx >= 0 && idx < n_bins)
                    S_out[idx] += Sb[b];
            }
        }
    double norm = 1.0 / ((double)Nx * (double)Ny * bin_w);
    for (int i = 0; i < n_bins; ++i)
        S_out[i] *= norm;
    free(omega);
    free(Sb);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_powder_spectrum_general(const irrep_magnon_lsw_t *L,
                                                      const int *sublattice_signs, int Nx, int Ny,
                                                      double omega_min, double omega_max,
                                                      int n_bins, double *S_out) {
    if (!L || !sublattice_signs || Nx <= 0 || Ny <= 0 || n_bins <= 0 ||
        omega_max <= omega_min || !S_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    memset(S_out, 0, (size_t)n_bins * sizeof *S_out);
    double  bin_w = (omega_max - omega_min) / n_bins;
    double *omega = malloc((size_t)n * sizeof *omega);
    double *Sb    = malloc((size_t)n * sizeof *Sb);
    if (!omega || !Sb) {
        free(omega);
        free(Sb);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = ((double)ix + 0.5) / Nx; /* avoid Goldstone at Γ */
            double fy = ((double)iy + 0.5) / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            irrep_magnon_structure_factor_general(L, sublattice_signs, kx, ky, omega, Sb);
            for (int b = 0; b < n; ++b) {
                int idx = (int)((omega[b] - omega_min) / bin_w);
                if (idx >= 0 && idx < n_bins)
                    S_out[idx] += Sb[b];
            }
        }
    double norm = 1.0 / ((double)Nx * (double)Ny * bin_w);
    for (int i = 0; i < n_bins; ++i)
        S_out[i] *= norm;
    free(omega);
    free(Sb);
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

irrep_status_t irrep_magnon_structure_factor_with_form_factor(const irrep_magnon_lsw_t *L,
                                                                const double (*positions)[2],
                                                                double qx, double qy,
                                                                double *omega_out,
                                                                double *S_perp_out) {
    if (!L || !positions || !omega_out || !S_perp_out) return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    double _Complex *u = malloc((size_t)n * n * sizeof *u);
    if (!u) return IRREP_ERR_OUT_OF_MEMORY;
    irrep_status_t st = irrep_magnon_dispersion(L, qx, qy, omega_out, u);
    if (st != IRREP_OK) {
        free(u);
        return st;
    }
    /* Pre-compute intra-cell phase factors. */
    double _Complex *phases = malloc((size_t)n * sizeof *phases);
    if (!phases) {
        free(u);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int a = 0; a < n; ++a) {
        double arg = qx * positions[a][0] + qy * positions[a][1];
        phases[a]  = cos(arg) + I * sin(arg);
    }
    for (int b = 0; b < n; ++b) {
        double _Complex sum = 0;
        for (int a = 0; a < n; ++a)
            sum += phases[a] * u[b * n + a];
        S_perp_out[b] = 2.0 * L->S * (creal(sum) * creal(sum) + cimag(sum) * cimag(sum));
    }
    free(u);
    free(phases);
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

/* Build the BdG matrix M(k) for a generic collinear ground state.
 * Same structure as irrep_magnon_dispersion_general embeds inline.
 * Caller-supplied output buffer M of size N*N where N = 2*n_sub. */
static void build_M_bdg_(const irrep_magnon_lsw_t *L, const int *sublattice_signs, double kx,
                          double ky, double _Complex *M) {
    int n = L->n_sub;
    int N = 2 * n;
    memset(M, 0, (size_t)N * N * sizeof(double _Complex));
    double S = L->S;

    for (int a = 0; a < n; ++a) {
        M[a * N + a] += 2.0 * L->Kz * S;
        M[(a + n) * N + (a + n)] += 2.0 * L->Kz * S;
    }

    for (int b = 0; b < L->n_bonds; ++b) {
        const irrep_magnon_bond_t *bd = &L->bonds[b];
        int    i = bd->bi, j = bd->bj;
        int    sigma = sublattice_signs[i] * sublattice_signs[j];
        double tx = bd->delta_x * L->a1[0] + bd->delta_y * L->a2[0];
        double ty = bd->delta_x * L->a1[1] + bd->delta_y * L->a2[1];
        double phase = kx * tx + ky * ty;
        double _Complex eikt = cos(phase) + I * sin(phase);
        double _Complex eikt_neg = cos(phase) - I * sin(phase);

        if (sigma > 0) {
            M[i * N + i] += -S * bd->J;
            M[j * N + j] += -S * bd->J;
            M[(i + n) * N + (i + n)] += -S * bd->J;
            M[(j + n) * N + (j + n)] += -S * bd->J;
            double _Complex hop_part = (S * bd->J - I * S * bd->D[2]) * eikt;
            M[i * N + j] += hop_part;
            M[j * N + i] += conj(hop_part);
            double _Complex hop_hole = (S * bd->J + I * S * bd->D[2]) * eikt;
            M[(i + n) * N + (j + n)] += hop_hole;
            M[(j + n) * N + (i + n)] += conj(hop_hole);
        } else {
            M[i * N + i] += S * bd->J;
            M[j * N + j] += S * bd->J;
            M[(i + n) * N + (i + n)] += S * bd->J;
            M[(j + n) * N + (j + n)] += S * bd->J;
            double _Complex pair_fwd = (S * bd->J - I * S * bd->D[2]) * eikt;
            double _Complex pair_bwd = (S * bd->J - I * S * bd->D[2]) * eikt_neg;
            M[i * N + (j + n)] += pair_fwd;
            M[j * N + (i + n)] += pair_bwd;
            M[(j + n) * N + i] += conj(pair_fwd);
            M[(i + n) * N + j] += conj(pair_bwd);
        }
    }
}

/* Compute the 3×3 rotation matrix R that takes +ẑ to a given unit
 * vector n. Column convention: R[3*row + col] (row-major). The
 * rotation is the canonical "tilt then twist" form. For n = ẑ, R = I.
 * For n = -ẑ, R = diag(1, -1, -1) (180° about x-axis). */
static void rotate_zhat_to_n_(const double n[3], double R[9]) {
    double nx = n[0], ny = n[1], nz = n[2];
    double r2_xy = nx * nx + ny * ny;
    if (r2_xy < 1e-24) {
        /* n is along +ẑ or -ẑ: use a simple R. */
        if (nz > 0) {
            R[0] = 1; R[1] = 0; R[2] = 0;
            R[3] = 0; R[4] = 1; R[5] = 0;
            R[6] = 0; R[7] = 0; R[8] = 1;
        } else {
            R[0] = 1; R[1] = 0; R[2] = 0;
            R[3] = 0; R[4] = -1; R[5] = 0;
            R[6] = 0; R[7] = 0; R[8] = -1;
        }
        return;
    }
    double r_xy = sqrt(r2_xy);
    double cos_phi = nx / r_xy;
    double sin_phi = ny / r_xy;
    double cos_theta = nz;
    double sin_theta = r_xy;
    /* Standard rotation: first rotate about y by θ (tilting ẑ toward n's
     * azimuth), then about z by φ. The columns of R are e_x', e_y', e_z'
     * (the local frame's basis vectors expressed in global coords): */
    R[0] = cos_phi * cos_theta;  R[1] = -sin_phi;       R[2] = cos_phi * sin_theta;
    R[3] = sin_phi * cos_theta;  R[4] =  cos_phi;       R[5] = sin_phi * sin_theta;
    R[6] =          -sin_theta;  R[7] =  0;             R[8] =          cos_theta;
}

/* Compute c[3][3] = J·M + DMI for a bond between sublattices with
 * rotation matrices R_i, R_j, Heisenberg J, and DMI vector D[3].
 * c_{de} = J·(R_i^T R_j)_{de} + Σ_klm D^k ε^{klm} R_i^{ld} R_j^{me} */
static void bond_c_matrix_(const double R_i[9], const double R_j[9], double J,
                            const double D[3], double c[9]) {
    /* M = R_i^T R_j. Note R[3*row + col]. Then M_de = Σ_k R_i[k][d] R_j[k][e]. */
    for (int d = 0; d < 3; ++d)
        for (int e = 0; e < 3; ++e) {
            double s = 0;
            for (int k = 0; k < 3; ++k)
                s += R_i[3 * k + d] * R_j[3 * k + e];
            c[3 * d + e] = J * s;
        }

    /* DMI: c_{de} += Σ_klm D^k ε^{klm} R_i^{ld} R_j^{me}.
     * ε^{klm} non-zero for (012, 120, 201) = +1 and (021, 102, 210) = -1. */
    static const int eps_idx[6][3] = {
        {0, 1, 2}, {1, 2, 0}, {2, 0, 1},  /* +1 */
        {0, 2, 1}, {1, 0, 2}, {2, 1, 0},  /* -1 */
    };
    static const double eps_sign[6] = {+1, +1, +1, -1, -1, -1};
    for (int d = 0; d < 3; ++d)
        for (int e = 0; e < 3; ++e) {
            double s = 0;
            for (int idx = 0; idx < 6; ++idx) {
                int k = eps_idx[idx][0];
                int l = eps_idx[idx][1];
                int m = eps_idx[idx][2];
                s += eps_sign[idx] * D[k] * R_i[3 * l + d] * R_j[3 * m + e];
            }
            c[3 * d + e] += s;
        }
}

irrep_status_t irrep_magnon_dispersion_noncollinear(const irrep_magnon_lsw_t *L,
                                                     const double *n_vectors, double kx, double ky,
                                                     double *omega_out) {
    if (!L || !n_vectors || !omega_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    int N = 2 * n;

    /* Compute rotation matrices for each sublattice. */
    double *R_arr = malloc((size_t)n * 9 * sizeof *R_arr);
    if (!R_arr)
        return IRREP_ERR_OUT_OF_MEMORY;
    for (int a = 0; a < n; ++a)
        rotate_zhat_to_n_(n_vectors + 3 * a, R_arr + 9 * a);

    /* Build BdG matrix M(k) of size 2n × 2n. */
    double _Complex *M = calloc((size_t)N * N, sizeof *M);
    if (!M) {
        free(R_arr);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    double S = L->S;

    /* Anisotropy on every site (in particle and hole blocks, both diagonals). */
    for (int a = 0; a < n; ++a) {
        M[a * N + a] += 2.0 * L->Kz * S;
        M[(a + n) * N + (a + n)] += 2.0 * L->Kz * S;
    }

    for (int b = 0; b < L->n_bonds; ++b) {
        const irrep_magnon_bond_t *bd = &L->bonds[b];
        int    i = bd->bi, j = bd->bj;
        double tx = bd->delta_x * L->a1[0] + bd->delta_y * L->a2[0];
        double ty = bd->delta_x * L->a1[1] + bd->delta_y * L->a2[1];
        double phase = kx * tx + ky * ty;
        double _Complex eikt = cos(phase) + I * sin(phase);

        /* Compute the 3×3 c matrix for this bond. */
        double c[9];
        bond_c_matrix_(R_arr + 9 * i, R_arr + 9 * j, bd->J, bd->D, c);
        double c_xx = c[3 * 0 + 0], c_xy = c[3 * 0 + 1], c_yx = c[3 * 1 + 0],
               c_yy = c[3 * 1 + 1], c_zz = c[3 * 2 + 2];

        /* Decompose c into 4 bilinear channels:
         * α_A: hopping a^†_i a_j coefficient
         * α_B: pairing a_i a_j coefficient
         */
        double _Complex alpha_A = 0.5 * S * (c_xx + c_yy) - 0.5 * I * S * (c_xy - c_yx);
        double _Complex alpha_B = 0.5 * S * (c_xx - c_yy) - 0.5 * I * S * (c_xy + c_yx);

        /* Diagonal "self-energy" from c_zz·S̃^z_i S̃^z_j: -c_zz·S per
         * endpoint, in BOTH particle and hole blocks. */
        M[i * N + i] += -c_zz * S;
        M[j * N + j] += -c_zz * S;
        M[(i + n) * N + (i + n)] += -c_zz * S;
        M[(j + n) * N + (j + n)] += -c_zz * S;

        /* Hopping (particle block): A_{i,j} += α_A·e^{ikt}, A_{j,i} += conj(α_A)·e^{-ikt}. */
        M[i * N + j] += alpha_A * eikt;
        M[j * N + i] += conj(alpha_A * eikt);

        /* Hopping (hole block): A^*_{i,j} += conj(α_A)·e^{ikt}. The hole block carries
         * the conjugate hopping (this is the BdG redundancy). */
        M[(i + n) * N + (j + n)] += conj(alpha_A) * eikt;
        M[(j + n) * N + (i + n)] += alpha_A * conj(eikt);

        /* Pairing (anomalous block): B_{i,j} += α_B·e^{ikt}, with the
         * canonical-commutation symmetric structure B_{j,i} = α_B·e^{-ikt}. */
        M[i * N + (j + n)] += alpha_B * eikt;
        M[j * N + (i + n)] += alpha_B * conj(eikt);
        M[(j + n) * N + i] += conj(alpha_B * eikt);
        M[(i + n) * N + j] += conj(alpha_B * conj(eikt));
    }

    /* Cholesky + Colpa diagonalisation (same machinery as
     * `_dispersion_general`). Regularise PSD matrix at gap-closing
     * Goldstone points. */
    double _Complex *K = malloc((size_t)N * N * sizeof *K);
    double _Complex *W = malloc((size_t)N * N * sizeof *W);
    double          *eigs = malloc((size_t)N * sizeof *eigs);
    double _Complex *evecs = malloc((size_t)N * N * sizeof *evecs);
    if (!K || !W || !eigs || !evecs) {
        free(R_arr); free(M); free(K); free(W); free(eigs); free(evecs);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    static const double EPS_PSD = 1e-10;
    for (int a = 0; a < N; ++a)
        M[a * N + a] += EPS_PSD;
    if (cholesky_(N, M, K) != 0) {
        free(R_arr); free(M); free(K); free(W); free(eigs); free(evecs);
        irrep_set_error_("irrep_magnon_dispersion_noncollinear: M(k) not positive-definite "
                          "at k=(%g, %g) — supplied n_α set may not be a stable ground state",
                          kx, ky);
        return IRREP_ERR_INVALID_ARG;
    }
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double _Complex s = 0;
            for (int l = 0; l < N; ++l) {
                double sg = (l < n) ? +1.0 : -1.0;
                s += sg * K[i * N + l] * conj(K[j * N + l]);
            }
            W[i * N + j] = s;
        }
    hermitian_eig_(N, W, eigs, evecs);
    int count = 0;
    for (int i = 0; i < N; ++i)
        if (eigs[i] > 0)
            omega_out[count++] = eigs[i];
    while (count < n)
        omega_out[count++] = 0.0;

    free(R_arr); free(M); free(K); free(W); free(eigs); free(evecs);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_dispersion_noncollinear_3d(const irrep_magnon_lsw_t *L,
                                                         const double *n_vectors,
                                                         const double a3[3], double kx, double ky,
                                                         double kz, double *omega_out) {
    if (!L || !n_vectors || !a3 || !omega_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    int N = 2 * n;

    double *R_arr = malloc((size_t)n * 9 * sizeof *R_arr);
    if (!R_arr)
        return IRREP_ERR_OUT_OF_MEMORY;
    for (int a = 0; a < n; ++a)
        rotate_zhat_to_n_(n_vectors + 3 * a, R_arr + 9 * a);

    double _Complex *M = calloc((size_t)N * N, sizeof *M);
    if (!M) {
        free(R_arr);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    double S = L->S;

    for (int a = 0; a < n; ++a) {
        M[a * N + a] += 2.0 * L->Kz * S;
        M[(a + n) * N + (a + n)] += 2.0 * L->Kz * S;
    }

    for (int b = 0; b < L->n_bonds; ++b) {
        const irrep_magnon_bond_t *bd = &L->bonds[b];
        int    i = bd->bi, j = bd->bj;
        /* 3D bond translation: t = δ_x·a₁ + δ_y·a₂ + δ_z·a₃ */
        double tx = bd->delta_x * L->a1[0] + bd->delta_y * L->a2[0] + bd->delta_z * a3[0];
        double ty = bd->delta_x * L->a1[1] + bd->delta_y * L->a2[1] + bd->delta_z * a3[1];
        double tz = bd->delta_z * a3[2];
        double phase = kx * tx + ky * ty + kz * tz;
        double _Complex eikt = cos(phase) + I * sin(phase);

        double c[9];
        bond_c_matrix_(R_arr + 9 * i, R_arr + 9 * j, bd->J, bd->D, c);
        double c_xx = c[0 * 3 + 0], c_xy = c[0 * 3 + 1], c_yx = c[1 * 3 + 0],
               c_yy = c[1 * 3 + 1], c_zz = c[2 * 3 + 2];

        double _Complex alpha_A = 0.5 * S * (c_xx + c_yy) - 0.5 * I * S * (c_xy - c_yx);
        double _Complex alpha_B = 0.5 * S * (c_xx - c_yy) - 0.5 * I * S * (c_xy + c_yx);

        M[i * N + i] += -c_zz * S;
        M[j * N + j] += -c_zz * S;
        M[(i + n) * N + (i + n)] += -c_zz * S;
        M[(j + n) * N + (j + n)] += -c_zz * S;

        M[i * N + j] += alpha_A * eikt;
        M[j * N + i] += conj(alpha_A * eikt);
        M[(i + n) * N + (j + n)] += conj(alpha_A) * eikt;
        M[(j + n) * N + (i + n)] += alpha_A * conj(eikt);

        M[i * N + (j + n)] += alpha_B * eikt;
        M[j * N + (i + n)] += alpha_B * conj(eikt);
        M[(j + n) * N + i] += conj(alpha_B * eikt);
        M[(i + n) * N + j] += conj(alpha_B * conj(eikt));
    }

    double _Complex *K = malloc((size_t)N * N * sizeof *K);
    double _Complex *W = malloc((size_t)N * N * sizeof *W);
    double          *eigs = malloc((size_t)N * sizeof *eigs);
    double _Complex *evecs = malloc((size_t)N * N * sizeof *evecs);
    if (!K || !W || !eigs || !evecs) {
        free(R_arr); free(M); free(K); free(W); free(eigs); free(evecs);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    static const double EPS_PSD = 1e-10;
    for (int a = 0; a < N; ++a)
        M[a * N + a] += EPS_PSD;
    if (cholesky_(N, M, K) != 0) {
        free(R_arr); free(M); free(K); free(W); free(eigs); free(evecs);
        irrep_set_error_("irrep_magnon_dispersion_noncollinear_3d: M(k) not positive-definite "
                          "at k=(%g, %g, %g)", kx, ky, kz);
        return IRREP_ERR_INVALID_ARG;
    }
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double _Complex s = 0;
            for (int l = 0; l < N; ++l) {
                double sg = (l < n) ? +1.0 : -1.0;
                s += sg * K[i * N + l] * conj(K[j * N + l]);
            }
            W[i * N + j] = s;
        }
    hermitian_eig_(N, W, eigs, evecs);
    int count = 0;
    for (int i = 0; i < N; ++i)
        if (eigs[i] > 0)
            omega_out[count++] = eigs[i];
    while (count < n)
        omega_out[count++] = 0.0;

    free(R_arr); free(M); free(K); free(W); free(eigs); free(evecs);
    return IRREP_OK;
}

/* 3D BdG matrix builder. Same as build_M_bdg_ but uses 3D primitive
 * vector a₃ and per-bond delta_z. */
static void build_M_bdg_3d_(const irrep_magnon_lsw_t *L, const int *sublattice_signs,
                             const double a3[3], double kx, double ky, double kz,
                             double _Complex *M) {
    int n = L->n_sub;
    int N = 2 * n;
    memset(M, 0, (size_t)N * N * sizeof(double _Complex));
    double S = L->S;

    for (int a = 0; a < n; ++a) {
        M[a * N + a] += 2.0 * L->Kz * S;
        M[(a + n) * N + (a + n)] += 2.0 * L->Kz * S;
    }

    for (int b = 0; b < L->n_bonds; ++b) {
        const irrep_magnon_bond_t *bd = &L->bonds[b];
        int    i = bd->bi, j = bd->bj;
        int    sigma = sublattice_signs[i] * sublattice_signs[j];
        double tx = bd->delta_x * L->a1[0] + bd->delta_y * L->a2[0] + bd->delta_z * a3[0];
        double ty = bd->delta_x * L->a1[1] + bd->delta_y * L->a2[1] + bd->delta_z * a3[1];
        double tz = bd->delta_z * a3[2];
        double phase = kx * tx + ky * ty + kz * tz;
        double _Complex eikt = cos(phase) + I * sin(phase);
        double _Complex eikt_neg = cos(phase) - I * sin(phase);

        if (sigma > 0) {
            M[i * N + i] += -S * bd->J;
            M[j * N + j] += -S * bd->J;
            M[(i + n) * N + (i + n)] += -S * bd->J;
            M[(j + n) * N + (j + n)] += -S * bd->J;
            double _Complex hop_part = (S * bd->J - I * S * bd->D[2]) * eikt;
            M[i * N + j] += hop_part;
            M[j * N + i] += conj(hop_part);
            double _Complex hop_hole = (S * bd->J + I * S * bd->D[2]) * eikt;
            M[(i + n) * N + (j + n)] += hop_hole;
            M[(j + n) * N + (i + n)] += conj(hop_hole);
        } else {
            M[i * N + i] += S * bd->J;
            M[j * N + j] += S * bd->J;
            M[(i + n) * N + (i + n)] += S * bd->J;
            M[(j + n) * N + (j + n)] += S * bd->J;
            double _Complex pair_fwd = (S * bd->J - I * S * bd->D[2]) * eikt;
            double _Complex pair_bwd = (S * bd->J - I * S * bd->D[2]) * eikt_neg;
            M[i * N + (j + n)] += pair_fwd;
            M[j * N + (i + n)] += pair_bwd;
            M[(j + n) * N + i] += conj(pair_fwd);
            M[(i + n) * N + j] += conj(pair_bwd);
        }
    }
}

irrep_status_t irrep_magnon_dispersion_general_3d(const irrep_magnon_lsw_t *L,
                                                   const int *sublattice_signs,
                                                   const double a3[3], double kx, double ky,
                                                   double kz, double *omega_out) {
    if (!L || !sublattice_signs || !a3 || !omega_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    int N = 2 * n;
    for (int a = 0; a < n; ++a)
        if (sublattice_signs[a] != +1 && sublattice_signs[a] != -1)
            return IRREP_ERR_INVALID_ARG;

    double _Complex *M = malloc((size_t)N * N * sizeof *M);
    double _Complex *K = malloc((size_t)N * N * sizeof *K);
    double _Complex *W = malloc((size_t)N * N * sizeof *W);
    double          *eigs = malloc((size_t)N * sizeof *eigs);
    double _Complex *evecs = malloc((size_t)N * N * sizeof *evecs);
    if (!M || !K || !W || !eigs || !evecs) {
        free(M); free(K); free(W); free(eigs); free(evecs);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    static const double EPS_PSD = 1e-10;

    build_M_bdg_3d_(L, sublattice_signs, a3, kx, ky, kz, M);
    for (int a = 0; a < N; ++a)
        M[a * N + a] += EPS_PSD;
    if (cholesky_(N, M, K) != 0) {
        free(M); free(K); free(W); free(eigs); free(evecs);
        irrep_set_error_("irrep_magnon_dispersion_general_3d: M(k) is not positive-definite "
                          "even with regularisation at k=(%g, %g, %g)", kx, ky, kz);
        return IRREP_ERR_INVALID_ARG;
    }
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double _Complex s = 0;
            for (int l = 0; l < N; ++l) {
                double sg = (l < n) ? +1.0 : -1.0;
                s += sg * K[i * N + l] * conj(K[j * N + l]);
            }
            W[i * N + j] = s;
        }
    hermitian_eig_(N, W, eigs, evecs);
    int count = 0;
    for (int i = 0; i < N; ++i)
        if (eigs[i] > 0)
            omega_out[count++] = eigs[i];
    while (count < n)
        omega_out[count++] = 0.0;

    free(M); free(K); free(W); free(eigs); free(evecs);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_afm_zero_point(const irrep_magnon_lsw_t *L,
                                            const int *sublattice_signs, int Nx, int Ny,
                                            double *delta_m_out) {
    if (!L || !sublattice_signs || Nx <= 0 || Ny <= 0 || !delta_m_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    int N = 2 * n;
    for (int a = 0; a < n; ++a)
        if (sublattice_signs[a] != +1 && sublattice_signs[a] != -1)
            return IRREP_ERR_INVALID_ARG;
    for (int a = 0; a < n; ++a)
        delta_m_out[a] = 0.0;

    double _Complex *M = malloc((size_t)N * N * sizeof *M);
    double _Complex *K = malloc((size_t)N * N * sizeof *K);
    double _Complex *W = malloc((size_t)N * N * sizeof *W);
    double          *eigs = malloc((size_t)N * sizeof *eigs);
    double _Complex *evecs = malloc((size_t)N * N * sizeof *evecs);
    double _Complex *v_band = malloc((size_t)N * sizeof *v_band);
    if (!M || !K || !W || !eigs || !evecs || !v_band) {
        free(M); free(K); free(W); free(eigs); free(evecs); free(v_band);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    static const double EPS_PSD = 1e-10;

    /* Half-shifted grid: f = (i + 0.5)/N. Avoids exact Goldstone-mode
     * sampling at k = 0 where M is singular and the |v|² integrand
     * diverges (integrably) in 2D AFM systems. */
    for (int iy = 0; iy < Ny; ++iy)
        for (int ix = 0; ix < Nx; ++ix) {
            double fx = ((double)ix + 0.5) / Nx;
            double fy = ((double)iy + 0.5) / Ny;
            double kx = fx * L->b1[0] + fy * L->b2[0];
            double ky = fx * L->b1[1] + fy * L->b2[1];
            build_M_bdg_(L, sublattice_signs, kx, ky, M);
            for (int a = 0; a < N; ++a)
                M[a * N + a] += EPS_PSD;
            if (cholesky_(N, M, K) != 0)
                continue; /* skip k-points where M is non-PD */
            /* W = K η K^† */
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    double _Complex s = 0;
                    for (int l = 0; l < N; ++l) {
                        double sg = (l < n) ? +1.0 : -1.0;
                        s += sg * K[i * N + l] * conj(K[j * N + l]);
                    }
                    W[i * N + j] = s;
                }
            hermitian_eig_(N, W, eigs, evecs);
            /* Holstein-Primakoff convention: c_α = a_α at every site
             * (magnon = +1 boson regardless of σ_α — the parallel/anti-
             * parallel split affects ONLY the bilinear hopping/pairing
             * structure, not the boson↔magnon mapping itself).
             *
             * With Ψ = (a_1, ..., a_n, a_1^†, ..., a_n^†)^T = T β
             * where β = (β_1, ..., β_n, β_1^†, ..., β_n^†)^T, the
             * coefficient of β^†_b in a_α is T[α, n + b].
             *
             * In the Colpa formulation, eigenvectors with *negative*
             * eigenvalues correspond to β^† columns (Holstein-Primakoff
             * particle-hole conjugation). hermitian_eig_ sorts ascending,
             * so negative eigenvalues are at indices 0..n-1.
             *
             * The unit-normalised W-eigenvectors give T = K⁻¹ψ with
             * T^† η T = sign(λ_b)/|λ_b| per column. To make T paraunitary
             * (T^† η T = ±1), multiply each column by √|λ_b| = √ω_b. So
             * the physical |v|² weight per band is |λ_b| · |K⁻¹ψ_b[α]|². */
            for (int b = 0; b < n; ++b) {
                for (int i = N - 1; i >= 0; --i) {
                    double _Complex sum = evecs[b * N + i];
                    for (int kk = i + 1; kk < N; ++kk)
                        sum -= K[i * N + kk] * v_band[kk];
                    v_band[i] = sum / K[i * N + i];
                }
                double abs_lambda = fabs(eigs[b]);
                for (int a = 0; a < n; ++a) {
                    double re = creal(v_band[a]);
                    double im = cimag(v_band[a]);
                    delta_m_out[a] += abs_lambda * (re * re + im * im);
                }
            }
        }

    /* Average over BZ. */
    double inv_NBZ = 1.0 / ((double)Nx * (double)Ny);
    for (int a = 0; a < n; ++a)
        delta_m_out[a] *= inv_NBZ;

    free(M); free(K); free(W); free(eigs); free(evecs); free(v_band);
    return IRREP_OK;
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

irrep_status_t irrep_magnon_structure_factor_general(const irrep_magnon_lsw_t *L,
                                                       const int *sublattice_signs, double qx,
                                                       double qy, double *omega_out,
                                                       double *S_perp_out) {
    if (!L || !sublattice_signs || !omega_out || !S_perp_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    int N = 2 * n;
    for (int a = 0; a < n; ++a)
        if (sublattice_signs[a] != +1 && sublattice_signs[a] != -1)
            return IRREP_ERR_INVALID_ARG;

    double _Complex *M = malloc((size_t)N * N * sizeof *M);
    double _Complex *K = malloc((size_t)N * N * sizeof *K);
    double _Complex *W = malloc((size_t)N * N * sizeof *W);
    double          *eigs = malloc((size_t)N * sizeof *eigs);
    double _Complex *evecs = malloc((size_t)N * N * sizeof *evecs);
    double _Complex *T_col = malloc((size_t)N * sizeof *T_col);
    if (!M || !K || !W || !eigs || !evecs || !T_col) {
        free(M); free(K); free(W); free(eigs); free(evecs); free(T_col);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    static const double EPS_PSD = 1e-10;

    build_M_bdg_(L, sublattice_signs, qx, qy, M);
    for (int a = 0; a < N; ++a)
        M[a * N + a] += EPS_PSD;
    if (cholesky_(N, M, K) != 0) {
        free(M); free(K); free(W); free(eigs); free(evecs); free(T_col);
        irrep_set_error_("structure_factor_general: M(k) not positive-definite");
        return IRREP_ERR_INVALID_ARG;
    }
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double _Complex s = 0;
            for (int l = 0; l < N; ++l) {
                double sg = (l < n) ? +1.0 : -1.0;
                s += sg * K[i * N + l] * conj(K[j * N + l]);
            }
            W[i * N + j] = s;
        }
    hermitian_eig_(N, W, eigs, evecs);

    /* Positive eigenvalues are at indices n..2n-1 (ascending sort).
     * For each such band b, compute T_b = K^{-1} ψ_b via back-substitution.
     * The Bogoliubov amplitudes:
     *   For σ_α = +1: u^b_α = T_b[α], v^b_α = T_b[α+n]
     *   For σ_α = -1: u^b_α = T_b[α+n], v^b_α = T_b[α]
     * Transverse INS: S_⊥_b(q) = 2S · |Σ_α (u + v)|² with paraunitary
     * normalisation by |λ_b|. */
    int count = 0;
    for (int idx = 0; idx < N; ++idx) {
        if (eigs[idx] <= 0) continue;
        for (int i = N - 1; i >= 0; --i) {
            double _Complex sum = evecs[idx * N + i];
            for (int kk = i + 1; kk < N; ++kk)
                sum -= K[i * N + kk] * T_col[kk];
            T_col[i] = sum / K[i * N + i];
        }
        double _Complex amp = 0;
        for (int a = 0; a < n; ++a) {
            int idx_u = (sublattice_signs[a] == +1) ? a : (a + n);
            int idx_v = (sublattice_signs[a] == +1) ? (a + n) : a;
            amp += T_col[idx_u] + T_col[idx_v];
        }
        double mag2 = creal(amp) * creal(amp) + cimag(amp) * cimag(amp);
        double abs_lambda = fabs(eigs[idx]);
        S_perp_out[count] = 2.0 * L->S * abs_lambda * mag2;
        omega_out[count] = eigs[idx];
        ++count;
    }
    while (count < n) {
        omega_out[count] = 0;
        S_perp_out[count] = 0;
        ++count;
    }

    free(M); free(K); free(W); free(eigs); free(evecs); free(T_col);
    return IRREP_OK;
}

irrep_status_t irrep_magnon_structure_factor_general_with_form_factor(
    const irrep_magnon_lsw_t *L, const int *sublattice_signs, const double (*positions)[2],
    double qx, double qy, double *omega_out, double *S_perp_out) {
    if (!L || !sublattice_signs || !positions || !omega_out || !S_perp_out)
        return IRREP_ERR_INVALID_ARG;
    int n = L->n_sub;
    int N = 2 * n;
    for (int a = 0; a < n; ++a)
        if (sublattice_signs[a] != +1 && sublattice_signs[a] != -1)
            return IRREP_ERR_INVALID_ARG;

    double _Complex *M = malloc((size_t)N * N * sizeof *M);
    double _Complex *K = malloc((size_t)N * N * sizeof *K);
    double _Complex *W = malloc((size_t)N * N * sizeof *W);
    double          *eigs = malloc((size_t)N * sizeof *eigs);
    double _Complex *evecs = malloc((size_t)N * N * sizeof *evecs);
    double _Complex *T_col = malloc((size_t)N * sizeof *T_col);
    double _Complex *phases = malloc((size_t)n * sizeof *phases);
    if (!M || !K || !W || !eigs || !evecs || !T_col || !phases) {
        free(M); free(K); free(W); free(eigs); free(evecs); free(T_col); free(phases);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    static const double EPS_PSD = 1e-10;

    build_M_bdg_(L, sublattice_signs, qx, qy, M);
    for (int a = 0; a < N; ++a)
        M[a * N + a] += EPS_PSD;
    if (cholesky_(N, M, K) != 0) {
        free(M); free(K); free(W); free(eigs); free(evecs); free(T_col); free(phases);
        irrep_set_error_("structure_factor_general_with_form_factor: M(k) not positive-definite");
        return IRREP_ERR_INVALID_ARG;
    }
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double _Complex s = 0;
            for (int l = 0; l < N; ++l) {
                double sg = (l < n) ? +1.0 : -1.0;
                s += sg * K[i * N + l] * conj(K[j * N + l]);
            }
            W[i * N + j] = s;
        }
    hermitian_eig_(N, W, eigs, evecs);

    for (int a = 0; a < n; ++a) {
        double arg = qx * positions[a][0] + qy * positions[a][1];
        phases[a]  = cos(arg) + I * sin(arg);
    }

    int count = 0;
    for (int idx = 0; idx < N; ++idx) {
        if (eigs[idx] <= 0) continue;
        for (int i = N - 1; i >= 0; --i) {
            double _Complex sum = evecs[idx * N + i];
            for (int kk = i + 1; kk < N; ++kk)
                sum -= K[i * N + kk] * T_col[kk];
            T_col[i] = sum / K[i * N + i];
        }
        double _Complex amp = 0;
        for (int a = 0; a < n; ++a) {
            int idx_u = (sublattice_signs[a] == +1) ? a : (a + n);
            int idx_v = (sublattice_signs[a] == +1) ? (a + n) : a;
            amp += phases[a] * (T_col[idx_u] + T_col[idx_v]);
        }
        double mag2       = creal(amp) * creal(amp) + cimag(amp) * cimag(amp);
        double abs_lambda = fabs(eigs[idx]);
        S_perp_out[count] = 2.0 * L->S * abs_lambda * mag2;
        omega_out[count]  = eigs[idx];
        ++count;
    }
    while (count < n) {
        omega_out[count]  = 0;
        S_perp_out[count] = 0;
        ++count;
    }

    free(M); free(K); free(W); free(eigs); free(evecs); free(T_col); free(phases);
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
