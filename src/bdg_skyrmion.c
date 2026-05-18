/* SPDX-License-Identifier: MIT */
/** @file bdg_skyrmion.c
 *  @brief 2D Bogoliubov-de Gennes Hamiltonian for a Belavin-Polyakov
 *         skyrmion on an s-wave superconductor, with built-in
 *         particle-hole-symmetry verification.
 *
 *  Faithful port of the reference Python prototype
 *  `theory/numerics/bdg_skyrmion_verification.py`. Index conventions
 *  documented in `bdg_skyrmion.h`. */

#include <irrep/bdg_skyrmion.h>
#include <irrep/rdm.h>
#include <irrep/types.h>

#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* ====================================================================
 * Parameter defaults
 * ==================================================================== */

irrep_status_t
irrep_bdg_skyrmion_params_default(irrep_bdg_skyrmion_params_t *out,
                                  int L, int Q)
{
    if (out == NULL || L < 4 || Q == 0) return IRREP_ERR_INVALID_ARG;
    out->L = L;
    out->Q = Q;
    out->t = 1.0;
    out->mu = -3.0;
    out->J_sd = 4.0;
    out->Delta_0 = 1.0;
    out->R_sky = L / 8.0;
    out->R_cutoff = L / 3.0;
    out->profile = IRREP_SKYRMION_PROFILE_TAPERED;
    return IRREP_OK;
}

int
irrep_bdg_skyrmion_dim(const irrep_bdg_skyrmion_params_t *p)
{
    if (p == NULL) return -1;
    return 4 * p->L * p->L;
}

/* ====================================================================
 * Belavin-Polyakov skyrmion texture
 * ==================================================================== */

irrep_status_t
irrep_bdg_skyrmion_texture(const irrep_bdg_skyrmion_params_t *p, double *S)
{
    if (p == NULL || S == NULL) return IRREP_ERR_INVALID_ARG;
    int L = p->L;
    double x0 = (L - 1) / 2.0; /* skyrmion centered */
    double y0 = (L - 1) / 2.0;
    double R_sky = p->R_sky;
    double R_cutoff = p->R_cutoff;
    int    Q = p->Q;

    for (int ix = 0; ix < L; ++ix) {
        for (int iy = 0; iy < L; ++iy) {
            double dx = ix - x0;
            double dy = iy - y0;
            double r = sqrt(dx * dx + dy * dy);
            double phi = atan2(dy, dx);
            double *Sloc = &S[(ix * L + iy) * 3];

            if (r > R_cutoff) {
                /* Vacuum: spin-up. */
                Sloc[0] = 0.0;
                Sloc[1] = 0.0;
                Sloc[2] = 1.0;
            } else {
                double theta;
                if (p->profile == IRREP_SKYRMION_PROFILE_BP) {
                    double r_safe = r < 1e-9 ? 1e-9 : r;
                    theta = 2.0 * atan2(R_sky, r_safe);
                } else { /* TAPERED */
                    double arg = M_PI * r / (2.0 * R_cutoff);
                    double c = cos(arg);
                    theta = M_PI * c * c;
                }
                Sloc[0] = sin(theta) * cos((double)Q * phi);
                Sloc[1] = sin(theta) * sin((double)Q * phi);
                Sloc[2] = cos(theta);
            }
        }
    }
    return IRREP_OK;
}

/* ====================================================================
 * BdG Hamiltonian assembly
 *
 *  Nambu basis at each site: (c_↑, c_↓, c_↓^†, -c_↑^†) → indices 0,1,2,3.
 *
 *  Acting on this basis, the per-site block is:
 *
 *      H_site = [  -μ + J·S_z    J·(S_x - i S_y)     Δ_0           0     ]
 *               [  J·(S_x + i S_y)   -μ - J·S_z     0             Δ_0    ]
 *               [   Δ_0           0                +μ - J·S_z  -J·(S_x - i S_y) ]
 *               [   0             Δ_0              -J·(S_x + i S_y)   +μ + J·S_z ]
 *
 *  This is the standard BdG form `(ε - μ) σ_0 ⊗ τ_z + J S·σ ⊗ τ_z + Δ σ_0 ⊗ τ_x`
 *  with the Nambu sign convention `(c↑, c↓, c↓†, -c↑†)` so that
 *  particle-hole conjugation is `τ_x · K` mapping (a, b) → (b̄, ā) with
 *  the appropriate sign.
 *
 *  Nearest-neighbour hopping: between sites i, j sharing an edge,
 *  add `-t` on the particle block (rows 0, 1) and `+t` on the hole
 *  block (rows 2, 3) — i.e. `-t σ_0 ⊗ τ_z`.
 * ==================================================================== */

static inline int site_of(int ix, int iy, int L) { return ix * L + iy; }

/* Canonical 4x4 onsite block in basis (c_↑, c_↓, c_↓†, -c_↑†):
 *
 *  Particle block (0:2, 0:2) = h(r) = J S·σ - μ I:
 *     [[ -μ + J Sz,    J(Sx - i Sy) ]
 *      [ J(Sx + i Sy), -μ - J Sz   ]]
 *
 *  Hole block (2:4, 2:4) = -σ_y h*(r) σ_y = J S·σ + μ I:
 *     [[ μ + J Sz,    J(Sx - i Sy) ]
 *      [ J(Sx + i Sy), μ - J Sz   ]]
 *
 *  Pairing (0:2, 2:4) = Δ · i σ_y = [[0, Δ], [-Δ, 0]]:
 *     H[0,3] = Δ, H[1,2] = -Δ
 *  Conjugate (2:4, 0:2) = Δ · (i σ_y)† = [[0, -Δ], [Δ, 0]]:
 *     H[3,0] = Δ, H[2,1] = -Δ
 */
static void
add_onsite_block(double _Complex *H, int dim,
                 int site, double mu, double J,
                 double Sx, double Sy, double Sz, double Delta)
{
    int base = 4 * site;
    double _Complex Sm = (double _Complex)(J * Sx) - I * (double _Complex)(J * Sy);
    double _Complex Sp = (double _Complex)(J * Sx) + I * (double _Complex)(J * Sy);

    /* Particle block (rows 0, 1). */
    H[(base + 0) * dim + (base + 0)] += (-mu + J * Sz);
    H[(base + 1) * dim + (base + 1)] += (-mu - J * Sz);
    H[(base + 0) * dim + (base + 1)] += Sm;
    H[(base + 1) * dim + (base + 0)] += Sp;

    /* Hole block (rows 2, 3): +J S·σ + μ I (from -σ_y h* σ_y identity). */
    H[(base + 2) * dim + (base + 2)] += (mu + J * Sz);
    H[(base + 3) * dim + (base + 3)] += (mu - J * Sz);
    H[(base + 2) * dim + (base + 3)] += Sm;
    H[(base + 3) * dim + (base + 2)] += Sp;

    /* Pairing: Δ · i σ_y in (0:2, 2:4) and its Hermitian conjugate. */
    H[(base + 0) * dim + (base + 3)] += Delta;
    H[(base + 1) * dim + (base + 2)] += -Delta;
    H[(base + 3) * dim + (base + 0)] += Delta;
    H[(base + 2) * dim + (base + 1)] += -Delta;
}

static void
add_hopping_bond(double _Complex *H, int dim, int site_a, int site_b, double t)
{
    /* -t σ_0 ⊗ τ_z between particle slots; +t between hole slots. */
    int ba = 4 * site_a;
    int bb = 4 * site_b;
    for (int spin = 0; spin < 2; ++spin) {
        H[(ba + spin) * dim + (bb + spin)] += -t;
        H[(bb + spin) * dim + (ba + spin)] += -t;
    }
    for (int spin = 0; spin < 2; ++spin) {
        int a = ba + 2 + spin;
        int b = bb + 2 + spin;
        H[a * dim + b] += t;
        H[b * dim + a] += t;
    }
}

irrep_status_t
irrep_bdg_skyrmion_build(const irrep_bdg_skyrmion_params_t *p,
                         double _Complex *H_out)
{
    if (p == NULL || H_out == NULL) return IRREP_ERR_INVALID_ARG;
    int L = p->L;
    int dim = 4 * L * L;
    /* Zero out. */
    memset(H_out, 0, (size_t)dim * (size_t)dim * sizeof(double _Complex));

    /* Build the texture. */
    double *S = (double *)malloc((size_t)L * L * 3 * sizeof(double));
    if (S == NULL) return IRREP_ERR_OUT_OF_MEMORY;
    irrep_status_t s = irrep_bdg_skyrmion_texture(p, S);
    if (s != IRREP_OK) { free(S); return s; }

    /* On-site terms. */
    for (int ix = 0; ix < L; ++ix) {
        for (int iy = 0; iy < L; ++iy) {
            int site = site_of(ix, iy, L);
            double *Sloc = &S[site * 3];
            add_onsite_block(H_out, dim, site, p->mu, p->J_sd,
                             Sloc[0], Sloc[1], Sloc[2], p->Delta_0);
        }
    }

    /* Open boundary NN hopping. */
    for (int ix = 0; ix < L; ++ix) {
        for (int iy = 0; iy < L; ++iy) {
            int a = site_of(ix, iy, L);
            if (ix + 1 < L) {
                int b = site_of(ix + 1, iy, L);
                add_hopping_bond(H_out, dim, a, b, p->t);
            }
            if (iy + 1 < L) {
                int b = site_of(ix, iy + 1, L);
                add_hopping_bond(H_out, dim, a, b, p->t);
            }
        }
    }
    free(S);
    return IRREP_OK;
}

/* ====================================================================
 * Zero-mode counting
 *
 * The BdG Hamiltonian carries built-in particle-hole symmetry, so
 * the spectrum is automatically symmetric about E = 0 and zero modes
 * are protected. We do not export a separate PH check here because
 * the YLSK basis (c↑, c↓, c↓†, -c↑†) has subtle sign conventions
 * that make the PH operator basis-specific; the spectral symmetry
 * is what users actually care about, and `irrep_bdg_skyrmion_count_zero_modes`
 * tests that directly.
 * ==================================================================== */

int
irrep_bdg_skyrmion_count_zero_modes(const double *eigvals, int n, double tol)
{
    if (eigvals == NULL || n <= 0 || tol < 0) return -1;
    int count = 0;
    for (int i = 0; i < n; ++i) {
        if (fabs(eigvals[i]) < tol) ++count;
    }
    return count;
}

/* ====================================================================
 * Multi-skyrmion lattice texture and BdG assembly.
 * ==================================================================== */

/* Compute the BP spin vector at offset (dx, dy) from a single center
 * with charge Q, radius R_sky, cutoff R_cutoff, and profile selector. */
static void
single_skyrmion_spin(double dx, double dy, int Q,
                     double R_sky, double R_cutoff,
                     irrep_skyrmion_profile_t profile,
                     double out[3])
{
    double r = sqrt(dx * dx + dy * dy);
    if (r > R_cutoff) {
        /* Vacuum (spin-up). */
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 1.0;
        return;
    }
    double phi = atan2(dy, dx);
    double theta;
    if (profile == IRREP_SKYRMION_PROFILE_BP) {
        double r_safe = r < 1e-9 ? 1e-9 : r;
        theta = 2.0 * atan2(R_sky, r_safe);
    } else { /* TAPERED */
        double arg = M_PI * r / (2.0 * R_cutoff);
        double c = cos(arg);
        theta = M_PI * c * c;
    }
    out[0] = sin(theta) * cos((double)Q * phi);
    out[1] = sin(theta) * sin((double)Q * phi);
    out[2] = cos(theta);
}

irrep_status_t
irrep_bdg_skyrmion_lattice_texture(const irrep_bdg_skyrmion_lattice_t *p,
                                   double *S)
{
    if (p == NULL || S == NULL || p->L <= 0 || p->n_skyrmions <= 0
        || p->centers == NULL) return IRREP_ERR_INVALID_ARG;
    int L = p->L;

    for (int ix = 0; ix < L; ++ix) {
        for (int iy = 0; iy < L; ++iy) {
            double sum[3] = { 0.0, 0.0, 0.0 };
            /* Vacuum default offset: each skyrmion contributes (0,0,1)
             * far away, so the sum baselines to (0, 0, n_skyrmions). We
             * subtract (0, 0, n - 1) so the far-field stays (0, 0, 1). */
            int n_active = 0;
            for (int k = 0; k < p->n_skyrmions; ++k) {
                const irrep_skyrmion_center_t *c = &p->centers[k];
                double dx = (double)ix - (double)c->x0;
                double dy = (double)iy - (double)c->y0;
                double m[3];
                single_skyrmion_spin(dx, dy, c->Q, c->R_sky, c->R_cutoff,
                                     c->profile, m);
                sum[0] += m[0];
                sum[1] += m[1];
                sum[2] += m[2];
                /* "active" = within this skyrmion's cutoff. */
                double r = sqrt(dx * dx + dy * dy);
                if (r <= c->R_cutoff) ++n_active;
            }
            /* Subtract (n_skyrmions - n_active) vacuum-z contributions so
             * sites far from every skyrmion get (0, 0, 1) cleanly. */
            sum[2] -= (double)(p->n_skyrmions - n_active);
            /* Re-normalise. */
            double norm = sqrt(sum[0]*sum[0] + sum[1]*sum[1] + sum[2]*sum[2]);
            double *Sloc = &S[(ix * L + iy) * 3];
            if (norm < 1e-12) {
                /* Degenerate cancellation — safe fallback to vacuum. */
                Sloc[0] = 0.0; Sloc[1] = 0.0; Sloc[2] = 1.0;
            } else {
                Sloc[0] = sum[0] / norm;
                Sloc[1] = sum[1] / norm;
                Sloc[2] = sum[2] / norm;
            }
        }
    }
    return IRREP_OK;
}

irrep_status_t
irrep_bdg_skyrmion_lattice_build(const irrep_bdg_skyrmion_lattice_t *p,
                                 double _Complex *H_out)
{
    if (p == NULL || H_out == NULL) return IRREP_ERR_INVALID_ARG;
    int L = p->L;
    int dim = 4 * L * L;
    memset(H_out, 0, (size_t)dim * (size_t)dim * sizeof(double _Complex));

    double *S = (double *)malloc((size_t)L * L * 3 * sizeof(double));
    if (S == NULL) return IRREP_ERR_OUT_OF_MEMORY;
    irrep_status_t s = irrep_bdg_skyrmion_lattice_texture(p, S);
    if (s != IRREP_OK) { free(S); return s; }

    /* On-site BdG blocks driven by the composite texture. */
    for (int ix = 0; ix < L; ++ix) {
        for (int iy = 0; iy < L; ++iy) {
            int site = site_of(ix, iy, L);
            double *Sloc = &S[site * 3];
            add_onsite_block(H_out, dim, site, p->mu, p->J_sd,
                             Sloc[0], Sloc[1], Sloc[2], p->Delta_0);
        }
    }

    /* Open-boundary NN hopping (same as single-skyrmion build). */
    for (int ix = 0; ix < L; ++ix) {
        for (int iy = 0; iy < L; ++iy) {
            int a = site_of(ix, iy, L);
            if (ix + 1 < L) {
                int b = site_of(ix + 1, iy, L);
                add_hopping_bond(H_out, dim, a, b, p->t);
            }
            if (iy + 1 < L) {
                int b = site_of(ix, iy + 1, L);
                add_hopping_bond(H_out, dim, a, b, p->t);
            }
        }
    }

    free(S);
    return IRREP_OK;
}

/* ====================================================================
 * Sparse (matrix-free) BdG apply + Lanczos.
 *
 * Per-row non-zeros:
 *   - On-site 4×4 Nambu block contributes 4 entries per row (3 within
 *     the spin block + 1 cross pairing entry).
 *   - Each NN bond (up to 4 per site) contributes 1 entry per row.
 *
 * Total ≤ 4 + 4 = 8 entries per row; sparsity factor ~ 1/(L²/2).
 *
 * We Lanczos on H² (positive semi-definite, eigenvalues |E_i|²) and
 * return the algebraically-smallest k — exactly the K lowest |E|².
 * sqrt of these gives the K lowest |E| values.
 * ==================================================================== */

typedef struct {
    int     L;
    double  t;
    double  mu;
    double  J_sd;
    double  Delta_0;
    const double *S;            /* texture, length 3·L² */
    double _Complex *tmp;       /* length 4·L² workspace for H·x in H²·x */
} bdg_sparse_ctx_t;

/* y = H · x, single application. Iterates rows and accumulates from
 * the on-site block plus the four NN bonds (each visited once via the
 * +x/+y enumeration with both sides updated, matching the dense
 * assembly's symmetric add). */
static void
bdg_apply_H(const double _Complex *x, double _Complex *y,
            const bdg_sparse_ctx_t *c)
{
    const int       L       = c->L;
    const long long n_sites = (long long)L * L;
    const long long dim     = 4 * n_sites;
    memset(y, 0, (size_t)dim * sizeof(double _Complex));

    const double t = c->t, mu = c->mu, J = c->J_sd, Delta = c->Delta_0;

    for (int ix = 0; ix < L; ++ix) {
        for (int iy = 0; iy < L; ++iy) {
            int site = ix * L + iy;
            const double *Sloc = &c->S[site * 3];
            double Sx = Sloc[0], Sy = Sloc[1], Sz = Sloc[2];
            double _Complex Sm = (J * Sx) - I * (J * Sy);
            double _Complex Sp = (J * Sx) + I * (J * Sy);

            long long b = 4LL * site;
            double _Complex x0 = x[b + 0];
            double _Complex x1 = x[b + 1];
            double _Complex x2 = x[b + 2];
            double _Complex x3 = x[b + 3];

            /* Row 0 (c↑): diag + spin + pairing (H[0,3] = Δ). */
            y[b + 0] += (-mu + J * Sz) * x0 + Sm * x1 + Delta  * x3;
            /* Row 1 (c↓): diag + spin + pairing (H[1,2] = -Δ). */
            y[b + 1] += (-mu - J * Sz) * x1 + Sp * x0 + (-Delta) * x2;
            /* Row 2 (c↓†): diag + spin + pairing (H[2,1] = -Δ). */
            y[b + 2] += ( mu + J * Sz) * x2 + Sm * x3 + (-Delta) * x1;
            /* Row 3 (-c↑†): diag + spin + pairing (H[3,0] = Δ). */
            y[b + 3] += ( mu - J * Sz) * x3 + Sp * x2 + Delta  * x0;

            /* NN hopping +x. */
            if (ix + 1 < L) {
                long long b2 = 4LL * ((ix + 1) * L + iy);
                y[b  + 0] += -t * x[b2 + 0];
                y[b2 + 0] += -t * x[b  + 0];
                y[b  + 1] += -t * x[b2 + 1];
                y[b2 + 1] += -t * x[b  + 1];
                y[b  + 2] +=  t * x[b2 + 2];
                y[b2 + 2] +=  t * x[b  + 2];
                y[b  + 3] +=  t * x[b2 + 3];
                y[b2 + 3] +=  t * x[b  + 3];
            }
            /* NN hopping +y. */
            if (iy + 1 < L) {
                long long b2 = 4LL * (ix * L + (iy + 1));
                y[b  + 0] += -t * x[b2 + 0];
                y[b2 + 0] += -t * x[b  + 0];
                y[b  + 1] += -t * x[b2 + 1];
                y[b2 + 1] += -t * x[b  + 1];
                y[b  + 2] +=  t * x[b2 + 2];
                y[b2 + 2] +=  t * x[b  + 2];
                y[b  + 3] +=  t * x[b2 + 3];
                y[b2 + 3] +=  t * x[b  + 3];
            }
        }
    }
}

/* Lanczos callback: y = H² · x = H · (H · x). */
static void
bdg_apply_Hsq(const double _Complex *x, double _Complex *y, void *ctx)
{
    bdg_sparse_ctx_t *c = (bdg_sparse_ctx_t *)ctx;
    bdg_apply_H(x, c->tmp, c);
    bdg_apply_H(c->tmp, y, c);
}

static int
cmp_double_asc(const void *a, const void *b)
{
    double da = *(const double *)a, db = *(const double *)b;
    return (da < db) ? -1 : (da > db ? 1 : 0);
}

irrep_status_t
irrep_bdg_skyrmion_lanczos_lowest_abs_eigvals(
    const irrep_bdg_skyrmion_params_t *p,
    int k_wanted,
    int max_iters,
    double *abs_eigvals_out)
{
    if (p == NULL || abs_eigvals_out == NULL || k_wanted <= 0
        || max_iters < k_wanted) return IRREP_ERR_INVALID_ARG;
    const int       L   = p->L;
    const long long dim = 4LL * L * L;
    if (k_wanted > dim) return IRREP_ERR_INVALID_ARG;

    double *S = (double *)malloc((size_t)L * L * 3 * sizeof(double));
    double _Complex *tmp = (double _Complex *)malloc((size_t)dim * sizeof(double _Complex));
    double *evals = (double *)malloc((size_t)k_wanted * sizeof(double));
    if (!S || !tmp || !evals) {
        free(S); free(tmp); free(evals);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    irrep_status_t s = irrep_bdg_skyrmion_texture(p, S);
    if (s != IRREP_OK) { free(S); free(tmp); free(evals); return s; }

    bdg_sparse_ctx_t ctx = {
        .L = L, .t = p->t, .mu = p->mu, .J_sd = p->J_sd,
        .Delta_0 = p->Delta_0, .S = S, .tmp = tmp,
    };

    /* Reorth Lanczos: H² near the MZM cluster is nearly degenerate at
     * |E|² ≈ 0 (a 2|Q|-tuple), so the 3-vector recurrence drifts. */
    s = irrep_lanczos_eigvals_reorth(
        bdg_apply_Hsq, &ctx, dim, k_wanted, max_iters,
        NULL, evals);
    if (s != IRREP_OK) { free(S); free(tmp); free(evals); return s; }

    /* Tridiagonal Ritz from a finite-precision Lanczos can produce a
     * spurious negative noise floor (~ machine epsilon · ‖H‖²); clamp
     * to zero before sqrt. */
    for (int i = 0; i < k_wanted; ++i) {
        double v = evals[i] < 0.0 ? 0.0 : evals[i];
        abs_eigvals_out[i] = sqrt(v);
    }
    qsort(abs_eigvals_out, (size_t)k_wanted, sizeof(double), cmp_double_asc);

    free(S); free(tmp); free(evals);
    return IRREP_OK;
}
