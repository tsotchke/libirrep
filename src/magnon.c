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
