/* SPDX-License-Identifier: MIT */
/* Symmetry-allowed Dzyaloshinskii-Moriya vector analysis.
 *
 * For each candidate symmetry operation, classify whether it preserves
 * or reverses the bond, build the contribution to the projector
 *   P = (1/|S|) Σ_g M_g    where M_g = ±R_proper(g)
 * (axial-vector representation: det collapses out), and extract the
 * allowed D subspace by diagonalising the resulting 3×3 matrix. */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/dmi.h>

extern void irrep_set_error_(const char *fmt, ...);

/* Apply (det · R) to v, writing the result into out. */
static void apply_op_(const irrep_dmi_sym_op_t *op, const double v[3], double out[3]) {
    double r[3];
    r[0] = op->R_proper[0] * v[0] + op->R_proper[1] * v[1] + op->R_proper[2] * v[2];
    r[1] = op->R_proper[3] * v[0] + op->R_proper[4] * v[1] + op->R_proper[5] * v[2];
    r[2] = op->R_proper[6] * v[0] + op->R_proper[7] * v[1] + op->R_proper[8] * v[2];
    if (op->det == -1) {
        out[0] = -r[0];
        out[1] = -r[1];
        out[2] = -r[2];
    } else {
        out[0] = r[0];
        out[1] = r[1];
        out[2] = r[2];
    }
}

static int vec_close_(const double a[3], const double b[3], double tol) {
    double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return (dx * dx + dy * dy + dz * dz) < tol * tol;
}

/* Symmetric 3×3 eigendecomposition by Jacobi rotations.
 * On return, `evals[i]` is sorted descending and `evecs` holds the
 * eigenvectors as ROWS (row-major). Tolerance ~ machine eps. */
static void sym3_eig_(const double M_in[9], double evals[3], double evecs[9]) {
    double A[9];
    memcpy(A, M_in, 9 * sizeof(double));
    /* Symmetrise just in case of round-off asymmetry. */
    A[1] = A[3] = 0.5 * (A[1] + A[3]);
    A[2] = A[6] = 0.5 * (A[2] + A[6]);
    A[5] = A[7] = 0.5 * (A[5] + A[7]);

    double V[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; /* identity */

    for (int sweep = 0; sweep < 50; ++sweep) {
        /* Find largest off-diagonal in absolute value. */
        int    p = 0, q = 1;
        double maxv = fabs(A[0 * 3 + 1]);
        if (fabs(A[0 * 3 + 2]) > maxv) {
            maxv = fabs(A[0 * 3 + 2]);
            p = 0;
            q = 2;
        }
        if (fabs(A[1 * 3 + 2]) > maxv) {
            maxv = fabs(A[1 * 3 + 2]);
            p = 1;
            q = 2;
        }
        if (maxv < 1e-15)
            break;

        double app = A[p * 3 + p], aqq = A[q * 3 + q], apq = A[p * 3 + q];
        double theta = 0.5 * atan2(2.0 * apq, app - aqq);
        double c = cos(theta), s = sin(theta);

        /* Rotate A and V. */
        double new_app = c * c * app + 2 * s * c * apq + s * s * aqq;
        double new_aqq = s * s * app - 2 * s * c * apq + c * c * aqq;
        A[p * 3 + p] = new_app;
        A[q * 3 + q] = new_aqq;
        A[p * 3 + q] = A[q * 3 + p] = 0.0;

        for (int k = 0; k < 3; ++k) {
            if (k == p || k == q)
                continue;
            double akp = A[k * 3 + p], akq = A[k * 3 + q];
            A[k * 3 + p] = A[p * 3 + k] = c * akp + s * akq;
            A[k * 3 + q] = A[q * 3 + k] = -s * akp + c * akq;
        }
        for (int k = 0; k < 3; ++k) {
            double vkp = V[k * 3 + p], vkq = V[k * 3 + q];
            V[k * 3 + p] = c * vkp + s * vkq;
            V[k * 3 + q] = -s * vkp + c * vkq;
        }
    }

    /* Sort descending by eigenvalue. V columns are eigenvectors → transpose to rows. */
    int idx[3] = {0, 1, 2};
    double diag[3] = {A[0], A[4], A[8]};
    /* Selection sort. */
    for (int i = 0; i < 2; ++i) {
        int best = i;
        for (int j = i + 1; j < 3; ++j)
            if (diag[idx[j]] > diag[idx[best]])
                best = j;
        if (best != i) {
            int t = idx[i];
            idx[i] = idx[best];
            idx[best] = t;
        }
    }
    for (int i = 0; i < 3; ++i) {
        evals[i] = diag[idx[i]];
        int col = idx[i];
        evecs[i * 3 + 0] = V[0 * 3 + col];
        evecs[i * 3 + 1] = V[1 * 3 + col];
        evecs[i * 3 + 2] = V[2 * 3 + col];
    }
}

int irrep_dmi_allowed_basis(const double r_a[3], const double r_b[3],
                            const irrep_dmi_sym_op_t *ops, int n_ops, double tol,
                            double basis_out[9]) {
    if (!r_a || !r_b || !ops || n_ops <= 0 || !basis_out) {
        irrep_set_error_("irrep_dmi_allowed_basis: invalid arguments");
        return -1;
    }

    /* Accumulate sum of M_g across symmetry-preserving operations.
     *
     * Sign rule for the DMI projector (axial 3-vector representation):
     *   unitary, preserving:    +R       (D = R·D constraint)
     *   unitary, reversing:     -R       (D = -R·D — Moriya rule)
     *   antiunitary, preserving: -R      (T flips axial vector — opposite of unitary)
     *   antiunitary, reversing:  +R      (T flip + bond-reversal flip cancel)
     *
     * See `docs/PHYSICS_APPENDIX.md` §14.3 for the derivation. */
    double M_sum[9] = {0};
    int    n_kept = 0;
    for (int k = 0; k < n_ops; ++k) {
        double g_a[3], g_b[3];
        apply_op_(&ops[k], r_a, g_a);
        apply_op_(&ops[k], r_b, g_b);
        int preserves = vec_close_(g_a, r_a, tol) && vec_close_(g_b, r_b, tol);
        int reverses = vec_close_(g_a, r_b, tol) && vec_close_(g_b, r_a, tol);
        if (!preserves && !reverses)
            continue;
        double sign = preserves ? +1.0 : -1.0;
        if (ops[k].antiunitary)
            sign = -sign;
        for (int i = 0; i < 9; ++i)
            M_sum[i] += sign * ops[k].R_proper[i];
        ++n_kept;
    }
    if (n_kept == 0) {
        irrep_set_error_("irrep_dmi_allowed_basis: no operation preserves the bond");
        return -1;
    }
    double inv_n = 1.0 / (double)n_kept;
    for (int i = 0; i < 9; ++i)
        M_sum[i] *= inv_n;

    /* P = M_sum is now a (real) projector onto the allowed subspace.
     * Diagonalise; eigenvectors with eigenvalue ~1 form the basis. */
    double evals[3], evecs[9];
    sym3_eig_(M_sum, evals, evecs);

    int n_basis = 0;
    for (int i = 0; i < 3; ++i) {
        if (fabs(evals[i] - 1.0) < 1e-9) {
            basis_out[n_basis * 3 + 0] = evecs[i * 3 + 0];
            basis_out[n_basis * 3 + 1] = evecs[i * 3 + 1];
            basis_out[n_basis * 3 + 2] = evecs[i * 3 + 2];
            ++n_basis;
        }
    }
    /* Zero-fill the unused rows for caller convenience. */
    for (int i = n_basis; i < 3; ++i) {
        basis_out[i * 3 + 0] = 0.0;
        basis_out[i * 3 + 1] = 0.0;
        basis_out[i * 3 + 2] = 0.0;
    }
    return n_basis;
}

int irrep_dmi_allowed_basis_from_pg(const double r_a[3], const double r_b[3],
                                    const irrep_pg_table_t *pg, double tol,
                                    double basis_out[9]) {
    if (!pg) {
        irrep_set_error_("irrep_dmi_allowed_basis_from_pg: pg is NULL");
        return -1;
    }
    int order = irrep_pg_order(pg);
    if (order <= 0) {
        irrep_set_error_("irrep_dmi_allowed_basis_from_pg: empty point group");
        return -1;
    }
    irrep_dmi_sym_op_t *ops = calloc((size_t)order, sizeof(*ops));
    if (!ops) {
        irrep_set_error_("irrep_dmi_allowed_basis_from_pg: out of memory");
        return -1;
    }
    for (int i = 0; i < order; ++i) {
        irrep_pg_element(pg, i, ops[i].R_proper, &ops[i].det);
        ops[i].antiunitary = 0;
    }

    int rc = irrep_dmi_allowed_basis(r_a, r_b, ops, order, tol, basis_out);
    free(ops);
    return rc;
}

/* -------------------------------------------------------------------------- *
 * Symmetric exchange tensor analyzer                                         *
 *                                                                            *
 * The symmetric basis on R^{3×3}_sym (six elements, orthonormal under        *
 * Frobenius inner product `<A, B> = trace(A^T B)`):                          *
 *                                                                            *
 *   e_0 = diag(1, 0, 0)                  (xx)                                *
 *   e_1 = diag(0, 1, 0)                  (yy)                                *
 *   e_2 = diag(0, 0, 1)                  (zz)                                *
 *   e_3 = (E_xy + E_yx) / √2             (xy/yx symmetric)                   *
 *   e_4 = (E_xz + E_zx) / √2             (xz/zx symmetric)                   *
 *   e_5 = (E_yz + E_zy) / √2             (yz/zy symmetric)                   *
 *                                                                            *
 * Under an operation R: e_α → R · e_α · R^T. Express in the basis as a       *
 * 6×6 matrix `M_g[β, α] = trace(e_β · R · e_α · R^T)`. Sum over the bond's   *
 * site-symmetry, normalise, diagonalise — eigenvalue-1 eigenvectors are      *
 * the symmetry-allowed J^s tensors. */

#define N_SYM_BASIS 6

static const double SYM_BASIS[N_SYM_BASIS][9] = {
    /* e0 */ {1, 0, 0, 0, 0, 0, 0, 0, 0},
    /* e1 */ {0, 0, 0, 0, 1, 0, 0, 0, 0},
    /* e2 */ {0, 0, 0, 0, 0, 0, 0, 0, 1},
    /* e3 */ {0, 0.70710678118654752440, 0, 0.70710678118654752440, 0, 0, 0, 0, 0},
    /* e4 */ {0, 0, 0.70710678118654752440, 0, 0, 0, 0.70710678118654752440, 0, 0},
    /* e5 */ {0, 0, 0, 0, 0, 0.70710678118654752440, 0, 0.70710678118654752440, 0}};

/* M = R · A · R^T for 3×3 matrices, all row-major flat. */
static void rart_3x3_(const double R[9], const double A[9], double out[9]) {
    double RA[9];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0;
            for (int k = 0; k < 3; ++k)
                s += R[i * 3 + k] * A[k * 3 + j];
            RA[i * 3 + j] = s;
        }
    /* RA · R^T */
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0;
            for (int k = 0; k < 3; ++k)
                s += RA[i * 3 + k] * R[j * 3 + k]; /* R^T[k,j] = R[j,k] */
            out[i * 3 + j] = s;
        }
}

static double frob_dot_(const double A[9], const double B[9]) {
    double s = 0;
    for (int i = 0; i < 9; ++i)
        s += A[i] * B[i];
    return s;
}

/* Generic NxN symmetric Jacobi diagonalisation. evals sorted descending. */
static void sym_eig_n_(int n, const double *M_in, double *evals, double *evecs) {
    double *A = malloc((size_t)n * n * sizeof(double));
    double *V = malloc((size_t)n * n * sizeof(double));
    memcpy(A, M_in, (size_t)n * n * sizeof(double));
    for (int i = 0; i < n * n; ++i)
        V[i] = 0;
    for (int i = 0; i < n; ++i)
        V[i * n + i] = 1.0;

    for (int sweep = 0; sweep < 80; ++sweep) {
        int    p = 0, q = 1;
        double maxv = 0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (fabs(A[i * n + j]) > maxv) {
                    maxv = fabs(A[i * n + j]);
                    p = i;
                    q = j;
                }
        if (maxv < 1e-15)
            break;

        double app = A[p * n + p], aqq = A[q * n + q], apq = A[p * n + q];
        double theta = 0.5 * atan2(2.0 * apq, app - aqq);
        double c = cos(theta), s = sin(theta);

        A[p * n + p] = c * c * app + 2 * s * c * apq + s * s * aqq;
        A[q * n + q] = s * s * app - 2 * s * c * apq + c * c * aqq;
        A[p * n + q] = A[q * n + p] = 0.0;

        for (int k = 0; k < n; ++k) {
            if (k == p || k == q)
                continue;
            double akp = A[k * n + p], akq = A[k * n + q];
            A[k * n + p] = A[p * n + k] = c * akp + s * akq;
            A[k * n + q] = A[q * n + k] = -s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
            double vkp = V[k * n + p], vkq = V[k * n + q];
            V[k * n + p] = c * vkp + s * vkq;
            V[k * n + q] = -s * vkp + c * vkq;
        }
    }

    int    *idx = malloc((size_t)n * sizeof(int));
    double *diag = malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; ++i) {
        idx[i] = i;
        diag[i] = A[i * n + i];
    }
    /* Selection sort descending. */
    for (int i = 0; i < n - 1; ++i) {
        int best = i;
        for (int j = i + 1; j < n; ++j)
            if (diag[idx[j]] > diag[idx[best]])
                best = j;
        if (best != i) {
            int t = idx[i];
            idx[i] = idx[best];
            idx[best] = t;
        }
    }
    for (int i = 0; i < n; ++i) {
        evals[i] = diag[idx[i]];
        for (int k = 0; k < n; ++k)
            evecs[i * n + k] = V[k * n + idx[i]]; /* row i = eigenvector i */
    }
    free(idx);
    free(diag);
    free(A);
    free(V);
}

int irrep_exchange_symmetric_basis(const double r_a[3], const double r_b[3],
                                   const irrep_dmi_sym_op_t *ops, int n_ops, double tol,
                                   double basis_out[54]) {
    if (!r_a || !r_b || !ops || n_ops <= 0 || !basis_out) {
        irrep_set_error_("irrep_exchange_symmetric_basis: invalid arguments");
        return -1;
    }

    /* Build 6×6 projector M_sum = Σ_g M_g over preserving + reversing ops. */
    double M_sum[N_SYM_BASIS * N_SYM_BASIS] = {0};
    int    n_kept = 0;
    for (int k = 0; k < n_ops; ++k) {
        double g_a[3], g_b[3];
        apply_op_(&ops[k], r_a, g_a);
        apply_op_(&ops[k], r_b, g_b);
        int preserves = vec_close_(g_a, r_a, tol) && vec_close_(g_b, r_b, tol);
        int reverses = vec_close_(g_a, r_b, tol) && vec_close_(g_b, r_a, tol);
        if (!preserves && !reverses)
            continue;
        /* Note: det collapses for axial-rank-2 (= axial × axial ↦ proper rank-2),
         * leaving R_proper · J · R_proper^T regardless of sign — same as for the
         * cross product in the DMI case. Both preserving and reversing give the
         * same constraint here because J^s is symmetric. */
        const double *R = ops[k].R_proper;
        for (int alpha = 0; alpha < N_SYM_BASIS; ++alpha) {
            double R_e_RT[9];
            rart_3x3_(R, SYM_BASIS[alpha], R_e_RT);
            for (int beta = 0; beta < N_SYM_BASIS; ++beta)
                M_sum[beta * N_SYM_BASIS + alpha] += frob_dot_(SYM_BASIS[beta], R_e_RT);
        }
        ++n_kept;
    }
    if (n_kept == 0) {
        irrep_set_error_("irrep_exchange_symmetric_basis: no operation preserves the bond");
        return -1;
    }
    double inv_n = 1.0 / (double)n_kept;
    for (int i = 0; i < N_SYM_BASIS * N_SYM_BASIS; ++i)
        M_sum[i] *= inv_n;

    /* Diagonalise; +1 eigenvectors are the allowed J^s coordinates in SYM_BASIS. */
    double evals[N_SYM_BASIS];
    double evecs[N_SYM_BASIS * N_SYM_BASIS];
    sym_eig_n_(N_SYM_BASIS, M_sum, evals, evecs);

    int n_basis = 0;
    for (int i = 0; i < N_SYM_BASIS; ++i) {
        if (fabs(evals[i] - 1.0) < 1e-9) {
            /* Reconstruct the 3×3 matrix from the 6-component eigenvector. */
            double J[9] = {0};
            for (int alpha = 0; alpha < N_SYM_BASIS; ++alpha) {
                double c = evecs[i * N_SYM_BASIS + alpha];
                for (int k = 0; k < 9; ++k)
                    J[k] += c * SYM_BASIS[alpha][k];
            }
            memcpy(basis_out + n_basis * 9, J, 9 * sizeof(double));
            ++n_basis;
        }
    }
    /* Zero-fill the unused matrices for caller convenience. */
    for (int i = n_basis; i < N_SYM_BASIS; ++i)
        memset(basis_out + i * 9, 0, 9 * sizeof(double));
    return n_basis;
}

/* -------------------------------------------------------------------------- *
 * Three-spin scalar chirality: χ_ijk = S_i · (S_j × S_k)                     *
 *                                                                            *
 * χ is a pseudoscalar: invariant under proper rotations, sign-flipped under  *
 * improper rotations AND under time reversal AND under odd permutations of   *
 * (i, j, k). For triangle invariance under each stabiliser element g:        *
 *                                                                            *
 *   sign(g) := det_spatial(g) · σ_perm(g) · σ_T(g)                           *
 *           where σ_perm = +1 if g cycles (a, b, c) → (b, c, a) etc.         *
 *                       = -1 if g transposes two and fixes one               *
 *                 σ_T   = -1 if antiunitary, +1 if unitary                   *
 *                                                                            *
 * χ is allowed iff sign(g) = +1 for every g preserving the triple as a set.  *
 * Any g with sign(g) = -1 forces χ = 0.                                      *
 * -------------------------------------------------------------------------- */

static int perm_sign_3_(int p0, int p1, int p2) {
    /* Sign of the permutation (p0, p1, p2) of (0, 1, 2). */
    int n_inv = 0;
    if (p0 > p1)
        ++n_inv;
    if (p0 > p2)
        ++n_inv;
    if (p1 > p2)
        ++n_inv;
    return (n_inv % 2 == 0) ? +1 : -1;
}

int irrep_chirality_allowed(const double r_a[3], const double r_b[3], const double r_c[3],
                            const irrep_dmi_sym_op_t *ops, int n_ops, double tol) {
    if (!r_a || !r_b || !r_c || !ops || n_ops <= 0)
        return -1;

    const double *r[3] = {r_a, r_b, r_c};
    int           n_kept = 0;
    for (int k = 0; k < n_ops; ++k) {
        double g_r[3][3];
        for (int i = 0; i < 3; ++i)
            apply_op_(&ops[k], r[i], g_r[i]);

        /* For each transformed g(r_i), find which input r_j it equals. */
        int perm[3] = {-1, -1, -1};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (vec_close_(g_r[i], r[j], tol)) {
                    perm[i] = j;
                    break;
                }
            }
            if (perm[i] < 0) {
                /* g(r_i) is not in the triple — operation not in stabiliser. */
                perm[0] = -1;
                break;
            }
        }
        if (perm[0] < 0)
            continue;
        /* Validate that perm is a true permutation (each index appears once). */
        int seen[3] = {0, 0, 0};
        int ok = 1;
        for (int i = 0; i < 3; ++i) {
            if (seen[perm[i]]) {
                ok = 0;
                break;
            }
            seen[perm[i]] = 1;
        }
        if (!ok)
            continue;

        ++n_kept;

        int sigma_perm = perm_sign_3_(perm[0], perm[1], perm[2]);
        int sigma_T = ops[k].antiunitary ? -1 : +1;
        int total_sign = ops[k].det * sigma_perm * sigma_T;
        if (total_sign != +1)
            return 0; /* Forbidden by this stabiliser element. */
    }
    if (n_kept == 0) {
        irrep_set_error_("irrep_chirality_allowed: no operation preserves the triple");
        return -1;
    }
    return 1; /* Every preserving op satisfies the constraint. */
}

int irrep_chirality_allowed_from_pg(const double r_a[3], const double r_b[3],
                                    const double r_c[3], const irrep_pg_table_t *pg,
                                    double tol) {
    if (!pg)
        return -1;
    int order = irrep_pg_order(pg);
    irrep_dmi_sym_op_t *ops = calloc((size_t)order, sizeof(*ops));
    if (!ops)
        return -1;
    for (int i = 0; i < order; ++i) {
        irrep_pg_element(pg, i, ops[i].R_proper, &ops[i].det);
        ops[i].antiunitary = 0;
    }
    int rc = irrep_chirality_allowed(r_a, r_b, r_c, ops, order, tol);
    free(ops);
    return rc;
}

int irrep_exchange_symmetric_basis_from_pg(const double r_a[3], const double r_b[3],
                                           const irrep_pg_table_t *pg, double tol,
                                           double basis_out[54]) {
    if (!pg) {
        irrep_set_error_("irrep_exchange_symmetric_basis_from_pg: pg is NULL");
        return -1;
    }
    int order = irrep_pg_order(pg);
    irrep_dmi_sym_op_t *ops = calloc((size_t)order, sizeof(*ops));
    if (!ops) {
        irrep_set_error_("irrep_exchange_symmetric_basis_from_pg: out of memory");
        return -1;
    }
    for (int i = 0; i < order; ++i) {
        irrep_pg_element(pg, i, ops[i].R_proper, &ops[i].det);
        ops[i].antiunitary = 0;
    }

    int rc = irrep_exchange_symmetric_basis(r_a, r_b, ops, order, tol, basis_out);
    free(ops);
    return rc;
}
