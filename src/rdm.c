/* SPDX-License-Identifier: MIT */
/* Reduced density matrix, entanglement entropies, and Kitaev-Preskill γ.
 *
 * Algorithms:
 *   - Partial trace: accumulate ρ_A[α, α'] = Σ_β ψ[α,β] conj(ψ[α',β]) by
 *     grouping full-Hilbert basis states by their B-digit tuple. For each
 *     β bucket we take the rank-1 outer product v_β v_β† and sum.
 *   - Hermitian eigendecomposition: cyclic-Jacobi sweeps with complex
 *     rotations. Eigenvalues accumulate on the diagonal as off-diagonal
 *     entries are rotated to zero. Terminates when the Frobenius norm of
 *     the off-diagonal part falls below a relative tolerance.
 *   - Entropies: natural-log convention; `0 ln 0 = 0` via a small-λ guard.
 *   - Kitaev-Preskill: literal sum of seven entropies.
 *
 * Cited:
 *   - Nielsen & Chuang, "Quantum Computation and Quantum Information"
 *     (Cambridge 2010), §11.3 for the partial-trace formula.
 *   - Kitaev & Preskill, Phys. Rev. Lett. 96, 110404 (2006) for the γ
 *     subtraction formula.
 *   - Golub & Van Loan, "Matrix Computations" (4e Johns Hopkins 2013),
 *     §8.5 for cyclic-Jacobi convergence bounds.
 */

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/rdm.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

extern void irrep_set_error_(const char *fmt, ...);

/* -------------------------------------------------------------------------- *
 * Partial trace                                                              *
 * -------------------------------------------------------------------------- */

static long long ipow_(int base, int exp) {
    long long r = 1;
    for (int k = 0; k < exp; ++k) r *= (long long)base;
    return r;
}

irrep_status_t irrep_partial_trace(int num_sites, int local_dim,
                                   const double _Complex *psi,
                                   const int *sites_A, int nA,
                                   double _Complex *rho_A) {
    if (num_sites < 1 || num_sites > 30 || local_dim < 2 || !psi || !rho_A)
        return IRREP_ERR_INVALID_ARG;
    if (nA < 0 || nA > num_sites || (nA > 0 && !sites_A))
        return IRREP_ERR_INVALID_ARG;

    /* Classify sites: flag[s] = 1 if s ∈ A, 0 if s ∈ B. A-indices are in
     * the order supplied by the caller; B-indices are taken in ascending
     * order of site index. */
    int in_A[32] = {0};
    for (int k = 0; k < nA; ++k) {
        int s = sites_A[k];
        if (s < 0 || s >= num_sites || in_A[s]) return IRREP_ERR_INVALID_ARG;
        in_A[s] = 1;
    }
    int nB = num_sites - nA;
    int sites_B[32];
    {
        int bi = 0;
        for (int s = 0; s < num_sites; ++s) if (!in_A[s]) sites_B[bi++] = s;
    }

    long long dA = ipow_(local_dim, nA);
    long long dB = ipow_(local_dim, nB);

    /* Precompute site weights w_s = local_dim^s so that i = Σ digit_s · w_s. */
    long long weight[32];
    for (int s = 0; s < num_sites; ++s) weight[s] = ipow_(local_dim, s);

    /* Reset output */
    memset(rho_A, 0, (size_t)(dA * dA) * sizeof(double _Complex));

    /* Scratch vector v_β of length dA. Reused each β. */
    double _Complex *v = malloc((size_t)dA * sizeof(double _Complex));
    if (!v) return IRREP_ERR_OUT_OF_MEMORY;

    for (long long beta = 0; beta < dB; ++beta) {
        /* Build the B-digit tuple from beta */
        long long b_remaining = beta;
        int b_digits[32] = {0};
        for (int k = 0; k < nB; ++k) {
            b_digits[k] = (int)(b_remaining % local_dim);
            b_remaining /= local_dim;
        }
        /* i_base is the full-basis index with all A-digits zero and B-digits
         * set. */
        long long i_base = 0;
        for (int k = 0; k < nB; ++k) i_base += (long long)b_digits[k] * weight[sites_B[k]];

        /* Fill v[alpha] = ψ[i_base + Σ alpha_digit_k · weight[sites_A[k]]]. */
        for (long long alpha = 0; alpha < dA; ++alpha) {
            long long a_remaining = alpha;
            long long offset = 0;
            for (int k = 0; k < nA; ++k) {
                int d = (int)(a_remaining % local_dim);
                a_remaining /= local_dim;
                offset += (long long)d * weight[sites_A[k]];
            }
            v[alpha] = psi[i_base + offset];
        }

        /* Accumulate ρ_A += v v† */
        for (long long a = 0; a < dA; ++a) {
            double _Complex va = v[a];
            double _Complex *row = rho_A + a * dA;
            for (long long c = 0; c < dA; ++c) {
                row[c] += va * conj(v[c]);
            }
        }
    }

    free(v);
    return IRREP_OK;
}

/* -------------------------------------------------------------------------- *
 * Cyclic-Jacobi Hermitian eigendecomposition.                                *
 *                                                                            *
 * Each rotation zeros one off-diagonal entry A[p,q] in two conceptual steps: *
 *                                                                            *
 *   1. Phase reduction. Let c = A[p,q] = |c| e^{iφ}. Apply V = diag(I, e^{-iφ}*
 *      at index q, I) as V†AV. This leaves A[q,q] and A[p,p] unchanged and   *
 *      scales row q by e^{iφ} (except (q,q)) and column q by e^{-iφ} (except *
 *      (q,q)); afterwards A[p,q] = A[q,p] = |c| (real).                      *
 *                                                                            *
 *   2. Real Givens. With the (p,q) block now real symmetric [[a, r],[r, b]], *
 *      the Givens rotation G = [[cos ψ, sin ψ],[−sin ψ, cos ψ]] with         *
 *      tan 2ψ = −2r / (a − b) gives G^T [[a,r],[r,b]] G diagonal. Applied as *
 *      G^T A G on rows/cols p, q leaves the rest of A untouched structurally *
 *      while zeroing A[p,q] exactly.                                         *
 *                                                                            *
 * Cyclic sweeps continue until the Frobenius norm of the off-diagonal part   *
 * drops below `tol · max(tr A, 1)`. Golub & Van Loan §8.5 quadratic          *
 * convergence near the fixed point gives ~10 sweeps for ~10⁻¹⁴ relative.     *
 * -------------------------------------------------------------------------- */

irrep_status_t irrep_hermitian_eigvals(int n, double _Complex *A, double *eigvals) {
    if (n < 1 || !A || !eigvals) return IRREP_ERR_INVALID_ARG;

    if (n == 1) { eigvals[0] = creal(A[0]); return IRREP_OK; }

    /* Symmetrise numerically (in case input carries Hermiticity noise). */
    for (int i = 0; i < n; ++i) {
        A[(size_t)i * n + i] = creal(A[(size_t)i * n + i]);
        for (int j = i + 1; j < n; ++j) {
            double _Complex avg = 0.5 * (A[(size_t)i*n + j] + conj(A[(size_t)j*n + i]));
            A[(size_t)i*n + j] = avg;
            A[(size_t)j*n + i] = conj(avg);
        }
    }

    const int    max_sweeps = 60 * n;
    const double tol        = 1e-14;

    double norm = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = creal(A[(size_t)i * n + i]);
        norm += d * d;
    }
    norm = sqrt(norm);
    double tol_abs = (norm > 1e-300 ? norm : 1.0) * tol;

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        double off = 0.0;
        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double _Complex apq = A[(size_t)p * n + q];
                off += 2.0 * (creal(apq)*creal(apq) + cimag(apq)*cimag(apq));
            }
        }
        if (sqrt(off) <= tol_abs) break;

        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double _Complex apq = A[(size_t)p * n + q];
                double r = cabs(apq);
                if (r < 1e-300) continue;

                /* ------- Step 1: phase reduction ------- */
                double _Complex e_pos = apq / r;       /* e^{iφ} */
                double _Complex e_neg = conj(e_pos);   /* e^{-iφ} */
                /* Row q scaled by e_pos (except the diagonal (q,q)),
                 * column q scaled by e_neg (except (q,q)). */
                for (int k = 0; k < n; ++k) {
                    if (k == q) continue;
                    A[(size_t)q * n + k] *= e_pos;
                    A[(size_t)k * n + q] *= e_neg;
                }
                /* After phase, A[p,q] = A[q,p] = r (real). */

                /* ------- Step 2: real Givens rotation ------- */
                double a = creal(A[(size_t)p * n + p]);
                double b = creal(A[(size_t)q * n + q]);
                double psi = (fabs(a - b) < 1e-300) ? (-M_PI / 4.0)
                                                    : 0.5 * atan2(-2.0 * r, a - b);
                double c_p = cos(psi), s_p = sin(psi);

                /* Row-update: A[p,k] ← c·A[p,k] − s·A[q,k]
                 *             A[q,k] ← s·A[p,k] + c·A[q,k] */
                for (int k = 0; k < n; ++k) {
                    double _Complex rp = A[(size_t)p * n + k];
                    double _Complex rq = A[(size_t)q * n + k];
                    A[(size_t)p * n + k] = c_p * rp - s_p * rq;
                    A[(size_t)q * n + k] = s_p * rp + c_p * rq;
                }
                /* Column-update: A[k,p] ← c·A[k,p] − s·A[k,q]
                 *                A[k,q] ← s·A[k,p] + c·A[k,q] */
                for (int k = 0; k < n; ++k) {
                    double _Complex cp = A[(size_t)k * n + p];
                    double _Complex cq = A[(size_t)k * n + q];
                    A[(size_t)k * n + p] = c_p * cp - s_p * cq;
                    A[(size_t)k * n + q] = s_p * cp + c_p * cq;
                }

                /* Clamp the now-zero (p,q) and (q,p) entries. */
                A[(size_t)p * n + q] = 0.0 + 0.0 * I;
                A[(size_t)q * n + p] = 0.0 + 0.0 * I;
            }
        }
    }

    for (int i = 0; i < n; ++i) eigvals[i] = creal(A[(size_t)i * n + i]);

    /* Sort descending (density-matrix convention). */
    for (int i = 0; i < n - 1; ++i) {
        int maxk = i;
        for (int k = i + 1; k < n; ++k) if (eigvals[k] > eigvals[maxk]) maxk = k;
        if (maxk != i) { double t = eigvals[i]; eigvals[i] = eigvals[maxk]; eigvals[maxk] = t; }
    }
    return IRREP_OK;
}

/* -------------------------------------------------------------------------- *
 * Entropies                                                                  *
 * -------------------------------------------------------------------------- */

double irrep_entropy_vonneumann_spectrum(const double *eigvals, int n) {
    if (!eigvals || n <= 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        double l = eigvals[i];
        if (l > 1e-15) s -= l * log(l);
    }
    return s;
}

double irrep_entropy_renyi_spectrum(const double *eigvals, int n, double alpha) {
    if (!eigvals || n <= 0) return 0.0;
    if (fabs(alpha - 1.0) < 1e-12) return irrep_entropy_vonneumann_spectrum(eigvals, n);
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        double l = eigvals[i];
        if (l > 1e-15) s += pow(l, alpha);
    }
    if (s <= 0.0) return 0.0;
    return log(s) / (1.0 - alpha);
}

static double entropy_fused_(const double _Complex *rho, int n, double alpha) {
    if (!rho || n <= 0) return 0.0;
    double _Complex *M = malloc((size_t)n * n * sizeof(double _Complex));
    double *ev = malloc((size_t)n * sizeof(double));
    if (!M || !ev) { free(M); free(ev); return NAN; }
    memcpy(M, rho, (size_t)n * n * sizeof(double _Complex));
    if (irrep_hermitian_eigvals(n, M, ev) != IRREP_OK) {
        free(M); free(ev); return NAN;
    }
    double s = (fabs(alpha - 1.0) < 1e-12)
             ? irrep_entropy_vonneumann_spectrum(ev, n)
             : irrep_entropy_renyi_spectrum    (ev, n, alpha);
    free(M); free(ev);
    return s;
}

double irrep_entropy_vonneumann(const double _Complex *rho, int n) {
    return entropy_fused_(rho, n, 1.0);
}

double irrep_entropy_renyi(const double _Complex *rho, int n, double alpha) {
    return entropy_fused_(rho, n, alpha);
}

/* -------------------------------------------------------------------------- *
 * Kitaev-Preskill topological entanglement entropy                           *
 * -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- *
 * Sparse Hermitian Lanczos: 3-term recurrence, no reorthogonalisation.       *
 *                                                                            *
 * Builds an orthogonal sequence v_0, v_1, … v_{k−1} via                      *
 *     β_{i+1} v_{i+1} = H v_i − α_i v_i − β_i v_{i−1}                         *
 *     α_i = ⟨v_i | H | v_i⟩,    β_{i+1} = ‖w‖                                 *
 * where w is the unnormalised RHS. Keeping only 3 state vectors at a time    *
 * keeps the algorithm memory-lean. The α, β scalars form a real symmetric   *
 * tridiagonal whose eigenvalues (Ritz values) converge to the extremal      *
 * spectrum of H, with the ground-state eigenvalue typically converging in   *
 * ~50 iterations for well-separated spectra.                                 *
 * -------------------------------------------------------------------------- */

irrep_status_t
irrep_lanczos_eigvals(
    void (*apply_op)(const double _Complex *x,
                     double _Complex *y,
                     void *ctx),
    void *ctx,
    long long dim,
    int k_wanted,
    int max_iters,
    const double _Complex *seed,
    double *eigvals_out) {
    if (!apply_op || dim <= 0 || k_wanted < 1 ||
        max_iters < 2 * k_wanted || !eigvals_out) return IRREP_ERR_INVALID_ARG;

    double _Complex *v_prev = calloc((size_t)dim, sizeof(double _Complex));
    double _Complex *v_curr = malloc((size_t)dim * sizeof(double _Complex));
    double _Complex *v_next = malloc((size_t)dim * sizeof(double _Complex));
    double          *alpha  = malloc((size_t)max_iters * sizeof(double));
    /* beta has one extra slot because the recurrence writes beta[j+1] at the
     * end of iteration j; the value at beta[max_iters] is never consumed but
     * must be within allocation. */
    double          *beta   = malloc((size_t)(max_iters + 1) * sizeof(double));
    if (!v_prev || !v_curr || !v_next || !alpha || !beta) {
        free(v_prev); free(v_curr); free(v_next); free(alpha); free(beta);
        return IRREP_ERR_OUT_OF_MEMORY;
    }

    /* Seed */
    if (seed) {
        memcpy(v_curr, seed, (size_t)dim * sizeof(double _Complex));
    } else {
        uint64_t rng = 0xdeadbeefcafeULL;
        for (long long i = 0; i < dim; ++i) {
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            double re = (double)(rng >> 32) / (double)0xFFFFFFFFULL - 0.5;
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            double im = (double)(rng >> 32) / (double)0xFFFFFFFFULL - 0.5;
            v_curr[i] = re + im * I;
        }
    }
    double norm = 0.0;
    for (long long i = 0; i < dim; ++i)
        norm += creal(v_curr[i]) * creal(v_curr[i]) + cimag(v_curr[i]) * cimag(v_curr[i]);
    norm = sqrt(norm);
    if (norm < 1e-300) {
        free(v_prev); free(v_curr); free(v_next); free(alpha); free(beta);
        return IRREP_ERR_PRECONDITION;
    }
    for (long long i = 0; i < dim; ++i) v_curr[i] /= norm;

    int n_iters = 0;
    beta[0] = 0.0;

    for (int j = 0; j < max_iters; ++j) {
        /* w = H · v_curr */
        apply_op(v_curr, v_next, ctx);

        /* α_j = Re ⟨v_curr | w⟩ (α is real because H is Hermitian). */
        double a = 0.0;
        for (long long i = 0; i < dim; ++i)
            a += creal(conj(v_curr[i]) * v_next[i]);
        alpha[j] = a;

        /* w ← w − α_j · v_curr − β_j · v_prev */
        for (long long i = 0; i < dim; ++i)
            v_next[i] -= a * v_curr[i] + beta[j] * v_prev[i];

        /* β_{j+1} = ‖w‖ */
        double b2 = 0.0;
        for (long long i = 0; i < dim; ++i)
            b2 += creal(v_next[i]) * creal(v_next[i]) + cimag(v_next[i]) * cimag(v_next[i]);
        double b = sqrt(b2);
        ++n_iters;

        if (b < 1e-14) break;   /* happy breakdown — invariant subspace found */
        beta[j + 1] = b;

        /* v_prev ← v_curr, v_curr ← v_next / β */
        double _Complex *tmp = v_prev;
        v_prev = v_curr;
        v_curr = v_next;
        v_next = tmp;
        for (long long i = 0; i < dim; ++i) v_curr[i] /= b;
    }

    /* Diagonalise the (n_iters × n_iters) real symmetric tridiagonal via
     * our existing Jacobi path, reusing a Hermitian complex wrapper. */
    double _Complex *T = calloc((size_t)n_iters * n_iters, sizeof(double _Complex));
    if (!T) {
        free(v_prev); free(v_curr); free(v_next); free(alpha); free(beta);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int i = 0; i < n_iters; ++i) {
        T[(size_t)i * n_iters + i] = alpha[i];
        if (i + 1 < n_iters) {
            T[(size_t)i * n_iters + (i + 1)] = beta[i + 1];
            T[(size_t)(i + 1) * n_iters + i] = beta[i + 1];
        }
    }
    double *ritz = malloc((size_t)n_iters * sizeof(double));
    if (!ritz) {
        free(T); free(v_prev); free(v_curr); free(v_next); free(alpha); free(beta);
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    irrep_hermitian_eigvals(n_iters, T, ritz);

    /* ritz is sorted descending; we want ascending for "lowest k". */
    if (k_wanted > n_iters) k_wanted = n_iters;
    for (int k = 0; k < k_wanted; ++k) {
        eigvals_out[k] = ritz[n_iters - 1 - k];
    }

    free(ritz);
    free(T);
    free(v_prev); free(v_curr); free(v_next); free(alpha); free(beta);
    return IRREP_OK;
}

double irrep_topological_entanglement_entropy(double SA, double SB, double SC,
                                              double SAB, double SBC, double SAC,
                                              double SABC) {
    return SA + SB + SC - SAB - SBC - SAC + SABC;
}
