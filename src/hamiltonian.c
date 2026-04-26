/* SPDX-License-Identifier: MIT */
/* Spin-½ Hamiltonian apply operators.
 *
 * Internal representation: a flat bond list with *per-bond* diagonal
 * (S^z S^z) and off-diagonal (½·(S^+ S^- + h.c.)) couplings. Every
 * exchange-type model in this module is a specialisation:
 *
 *   Heisenberg:  coeff_zz[b] = J/4,   coeff_pm[b] = J/2
 *   XY:          coeff_zz[b] = 0,     coeff_pm[b] = J/2
 *   J₁–J₂:       coeff_zz / coeff_pm vary between NN and NNN segments
 *
 * Same apply loop for all three — data-driven, no branching on model
 * type inside the hot path. This is the canonical payoff of promoting
 * the ad-hoc apply_H_heisenberg to a library primitive: extending the
 * model surface does not cost a new apply, just a new constructor.
 *
 * Hot loop (irrep_heisenberg_apply):
 *   for each bond b, endpoints (i, j):
 *     for each basis state s:
 *       zz_sign = (bit_i ^ bit_j) ? −1 : +1
 *       out[s]             += coeff_zz[b] · zz_sign · psi[s]          (diag)
 *       if bit_i != bit_j and coeff_pm[b] != 0:
 *         out[s ⊕ (1<<i) ⊕ (1<<j)] += coeff_pm[b] · psi[s]             (flip)
 */

#include <complex.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/hamiltonian.h>

extern void irrep_set_error_(const char *fmt, ...);

struct irrep_heisenberg {
    int     num_sites;
    int     num_bonds;
    int    *bi;
    int    *bj;
    double *coeff_zz; /* per-bond S^z_i S^z_j coefficient  */
    double *coeff_pm; /* per-bond ½(S^+_i S^−_j + h.c.)    */
};

/* Shared allocator. Copies the bond arrays + per-bond coefficient
 * arrays. The caller owns the input arrays. */
static irrep_heisenberg_t *build_(int num_sites, int num_bonds, const int *bi, const int *bj,
                                  const double *coeff_zz, const double *coeff_pm) {
    if (num_sites < 1 || num_sites > 62 || num_bonds < 0 || !bi || !bj) {
        irrep_set_error_("irrep_heisenberg build: invalid arguments");
        return NULL;
    }
    for (int b = 0; b < num_bonds; ++b) {
        if (bi[b] < 0 || bi[b] >= num_sites || bj[b] < 0 || bj[b] >= num_sites || bi[b] == bj[b]) {
            irrep_set_error_("irrep_heisenberg build: bond %d out of range", b);
            return NULL;
        }
    }
    irrep_heisenberg_t *H = calloc(1, sizeof *H);
    if (!H) {
        irrep_set_error_("irrep_heisenberg build: OOM");
        return NULL;
    }
    H->num_sites = num_sites;
    H->num_bonds = num_bonds;
    if (num_bonds > 0) {
        H->bi = malloc((size_t)num_bonds * sizeof(int));
        H->bj = malloc((size_t)num_bonds * sizeof(int));
        H->coeff_zz = malloc((size_t)num_bonds * sizeof(double));
        H->coeff_pm = malloc((size_t)num_bonds * sizeof(double));
        if (!H->bi || !H->bj || !H->coeff_zz || !H->coeff_pm) {
            free(H->bi);
            free(H->bj);
            free(H->coeff_zz);
            free(H->coeff_pm);
            free(H);
            irrep_set_error_("irrep_heisenberg build: OOM (bond arrays)");
            return NULL;
        }
        memcpy(H->bi, bi, (size_t)num_bonds * sizeof(int));
        memcpy(H->bj, bj, (size_t)num_bonds * sizeof(int));
        memcpy(H->coeff_zz, coeff_zz, (size_t)num_bonds * sizeof(double));
        memcpy(H->coeff_pm, coeff_pm, (size_t)num_bonds * sizeof(double));
    }
    return H;
}

irrep_heisenberg_t *irrep_heisenberg_new(int num_sites, int num_bonds, const int *bi, const int *bj,
                                         double J) {
    if (num_bonds < 0)
        num_bonds = 0;
    double *czz = malloc((size_t)(num_bonds ? num_bonds : 1) * sizeof(double));
    double *cpm = malloc((size_t)(num_bonds ? num_bonds : 1) * sizeof(double));
    if (!czz || !cpm) {
        free(czz);
        free(cpm);
        return NULL;
    }
    for (int b = 0; b < num_bonds; ++b) {
        czz[b] = 0.25 * J;
        cpm[b] = 0.5 * J;
    }
    irrep_heisenberg_t *H = build_(num_sites, num_bonds, bi, bj, czz, cpm);
    free(czz);
    free(cpm);
    return H;
}

irrep_heisenberg_t *irrep_heisenberg_j1j2_new(int num_sites, int num_bonds_nn, const int *nn_i,
                                              const int *nn_j, double J1, int num_bonds_nnn,
                                              const int *nnn_i, const int *nnn_j, double J2) {
    if (num_bonds_nn < 0 || num_bonds_nnn < 0) {
        irrep_set_error_("irrep_heisenberg_j1j2_new: negative bond count");
        return NULL;
    }
    int     total = num_bonds_nn + num_bonds_nnn;
    int    *bi = malloc((size_t)(total ? total : 1) * sizeof(int));
    int    *bj = malloc((size_t)(total ? total : 1) * sizeof(int));
    double *czz = malloc((size_t)(total ? total : 1) * sizeof(double));
    double *cpm = malloc((size_t)(total ? total : 1) * sizeof(double));
    if (!bi || !bj || !czz || !cpm) {
        free(bi);
        free(bj);
        free(czz);
        free(cpm);
        irrep_set_error_("irrep_heisenberg_j1j2_new: OOM");
        return NULL;
    }
    for (int b = 0; b < num_bonds_nn; ++b) {
        bi[b] = nn_i[b];
        bj[b] = nn_j[b];
        czz[b] = 0.25 * J1;
        cpm[b] = 0.5 * J1;
    }
    for (int b = 0; b < num_bonds_nnn; ++b) {
        bi[num_bonds_nn + b] = nnn_i[b];
        bj[num_bonds_nn + b] = nnn_j[b];
        czz[num_bonds_nn + b] = 0.25 * J2;
        cpm[num_bonds_nn + b] = 0.5 * J2;
    }
    irrep_heisenberg_t *H = build_(num_sites, total, bi, bj, czz, cpm);
    free(bi);
    free(bj);
    free(czz);
    free(cpm);
    return H;
}

irrep_heisenberg_t *irrep_xy_new(int num_sites, int num_bonds, const int *bi, const int *bj,
                                 double J) {
    if (num_bonds < 0)
        num_bonds = 0;
    double *czz = malloc((size_t)(num_bonds ? num_bonds : 1) * sizeof(double));
    double *cpm = malloc((size_t)(num_bonds ? num_bonds : 1) * sizeof(double));
    if (!czz || !cpm) {
        free(czz);
        free(cpm);
        return NULL;
    }
    for (int b = 0; b < num_bonds; ++b) {
        czz[b] = 0.0; /* XY has no S^z S^z term */
        cpm[b] = 0.5 * J;
    }
    irrep_heisenberg_t *H = build_(num_sites, num_bonds, bi, bj, czz, cpm);
    free(czz);
    free(cpm);
    return H;
}

void irrep_heisenberg_free(irrep_heisenberg_t *H) {
    if (!H)
        return;
    free(H->bi);
    free(H->bj);
    free(H->coeff_zz);
    free(H->coeff_pm);
    free(H);
}

int irrep_heisenberg_num_sites(const irrep_heisenberg_t *H) {
    return H ? H->num_sites : 0;
}

long long irrep_heisenberg_dim(const irrep_heisenberg_t *H) {
    return H ? (1LL << H->num_sites) : 0;
}

#include <irrep/config_project.h>
#include <irrep/space_group.h>
#include <math.h>

void irrep_heisenberg_apply_in_sector(const irrep_heisenberg_t *H,
                                      const struct irrep_sg_rep_table *T,
                                      const double _Complex *psi_in,
                                      double _Complex *psi_out) {
    if (!H || !T || !psi_in || !psi_out)
        return;
    long long n_reps = irrep_sg_rep_table_count(T);
    if (n_reps <= 0)
        return;
    const irrep_space_group_t *G = irrep_sg_rep_table_space_group(T);
    if (!G)
        return;
    memset(psi_out, 0, (size_t)n_reps * sizeof(double _Complex));

    for (long long k = 0; k < n_reps; ++k) {
        uint64_t u       = irrep_sg_rep_table_get(T, k);
        int      orbit_u = irrep_sg_rep_table_orbit_size(T, k);
        double   sqrt_u  = sqrt((double)orbit_u);
        double _Complex xu = psi_in[k];
        if (xu == 0.0)
            continue;

        double _Complex diag_acc = 0.0 + 0.0 * I;
        for (int b = 0; b < H->num_bonds; ++b) {
            int i = H->bi[b], j = H->bj[b];
            int zi = (int)((u >> i) & 1ULL);
            int zj = (int)((u >> j) & 1ULL);
            /* Diagonal S^z_i S^z_j: H_d is translation-invariant so the
             * sector-basis diagonal equals the bit-pattern expectation on
             * the representative u itself. */
            diag_acc += ((zi ^ zj) ? -H->coeff_zz[b] : H->coeff_zz[b]) * xu;

            /* Off-diagonal S^+ S^- + h.c. on an anti-aligned pair:
             *   ⟨ṽ|H|ũ⟩ = c_pm · K(u→v) · √(orbit_u / orbit_v)
             * per bond contribution — source orbit size in numerator. */
            if ((zi ^ zj) && H->coeff_pm[b] != 0.0) {
                uint64_t uflip = u ^ ((uint64_t)1 << i) ^ ((uint64_t)1 << j);
                uint64_t v; int g_idx;
                irrep_sg_canonicalise(G, uflip, &v, &g_idx);
                long long kv = irrep_sg_rep_table_index(T, v);
                if (kv >= 0) {
                    int    orbit_v = irrep_sg_rep_table_orbit_size(T, kv);
                    double factor  = sqrt_u / sqrt((double)orbit_v);
                    psi_out[kv] += H->coeff_pm[b] * factor * xu;
                }
            }
        }
        psi_out[k] += diag_acc;
    }
}

struct irrep_sg_heisenberg_sector {
    long long n_reps;             /* sector dimension (≤ rep_table count) */
    double   *diag;               /* [n_reps] per-rep total diagonal coefficient */
    long long *rowptr;            /* [n_reps + 1] CSR row pointers */
    long long *col;               /* [nnz] target sector-basis indices */
    double   *coef_real;          /* [nnz] real coefs (trivial sector); NULL otherwise */
    double _Complex *coef_cplx;   /* [nnz] complex coefs ((k, μ_k) sector); NULL otherwise */
};

irrep_sg_heisenberg_sector_t *
irrep_sg_heisenberg_sector_build(const irrep_heisenberg_t *H,
                                 const struct irrep_sg_rep_table *T) {
    if (!H || !T)
        return NULL;
    const irrep_space_group_t *G = irrep_sg_rep_table_space_group(T);
    if (!G)
        return NULL;
    long long n_reps = irrep_sg_rep_table_count(T);
    if (n_reps <= 0)
        return NULL;

    irrep_sg_heisenberg_sector_t *S = calloc(1, sizeof(*S));
    if (!S)
        return NULL;
    S->n_reps = n_reps;
    S->diag   = calloc((size_t)n_reps, sizeof(double));
    S->rowptr = malloc((size_t)(n_reps + 1) * sizeof(long long));
    int *counts = calloc((size_t)n_reps, sizeof(int));
    if (!S->diag || !S->rowptr || !counts) {
        free(counts);
        irrep_sg_heisenberg_sector_free(S);
        return NULL;
    }

    /* Two-pass construction: (1) count + diag, (2) fill CSR. Both passes
     * are parallel over k — each thread writes to a distinct k-indexed
     * slot with no cross-thread ordering. OpenMP is opt-in via the
     * USE_OPENMP=1 build flag; without it this reduces to sequential
     * execution (pragmas become no-ops). */
#pragma omp parallel for schedule(dynamic, 64)
    for (long long k = 0; k < n_reps; ++k) {
        uint64_t u = irrep_sg_rep_table_get(T, k);
        double   diag_acc = 0.0;
        int      cnt      = 0;

        for (int b = 0; b < H->num_bonds; ++b) {
            int i = H->bi[b], j = H->bj[b];
            int zi = (int)((u >> i) & 1ULL);
            int zj = (int)((u >> j) & 1ULL);
            diag_acc += ((zi ^ zj) ? -H->coeff_zz[b] : H->coeff_zz[b]);

            if ((zi ^ zj) && H->coeff_pm[b] != 0.0) {
                uint64_t uflip = u ^ ((uint64_t)1 << i) ^ ((uint64_t)1 << j);
                uint64_t v; int g_idx;
                irrep_sg_canonicalise(G, uflip, &v, &g_idx);
                long long kv = irrep_sg_rep_table_index(T, v);
                if (kv >= 0)
                    ++cnt;
            }
        }
        S->diag[k] = diag_acc;
        counts[k]  = cnt;
    }

    /* Prefix-sum → rowptr; total nnz at the tail. */
    S->rowptr[0] = 0;
    for (long long k = 0; k < n_reps; ++k)
        S->rowptr[k + 1] = S->rowptr[k] + counts[k];
    long long nnz = S->rowptr[n_reps];

    S->col       = malloc((size_t)(nnz ? nnz : 1) * sizeof(long long));
    S->coef_real = malloc((size_t)(nnz ? nnz : 1) * sizeof(double));
    if ((!S->col || !S->coef_real) && nnz > 0) {
        free(counts);
        irrep_sg_heisenberg_sector_free(S);
        return NULL;
    }

    /* Fill pass: each thread writes to rows [rowptr[k], rowptr[k+1]). */
#pragma omp parallel for schedule(dynamic, 64)
    for (long long k = 0; k < n_reps; ++k) {
        uint64_t u       = irrep_sg_rep_table_get(T, k);
        int      orbit_u = irrep_sg_rep_table_orbit_size(T, k);
        double   sqrt_u  = sqrt((double)orbit_u);
        long long pos    = S->rowptr[k];

        for (int b = 0; b < H->num_bonds; ++b) {
            int i = H->bi[b], j = H->bj[b];
            int zi = (int)((u >> i) & 1ULL);
            int zj = (int)((u >> j) & 1ULL);
            if ((zi ^ zj) && H->coeff_pm[b] != 0.0) {
                uint64_t uflip = u ^ ((uint64_t)1 << i) ^ ((uint64_t)1 << j);
                uint64_t v; int g_idx;
                irrep_sg_canonicalise(G, uflip, &v, &g_idx);
                long long kv = irrep_sg_rep_table_index(T, v);
                if (kv >= 0) {
                    int    orbit_v = irrep_sg_rep_table_orbit_size(T, kv);
                    double factor  = sqrt_u / sqrt((double)orbit_v);
                    S->col[pos]       = kv;
                    S->coef_real[pos] = H->coeff_pm[b] * factor;
                    ++pos;
                }
            }
        }
    }

    free(counts);
    return S;
}

irrep_sg_heisenberg_sector_t *
irrep_sg_heisenberg_sector_build_at_k(const irrep_heisenberg_t *H,
                                      const struct irrep_sg_rep_table *T,
                                      const struct irrep_sg_little_group *lg,
                                      const struct irrep_sg_little_group_irrep *mu_k) {
    if (!H || !T || !lg || !mu_k)
        return NULL;
    int d_mu = irrep_sg_little_group_irrep_dim(mu_k);
    if (d_mu != 1 && d_mu != 2)
        return NULL;

    /* For 2D irreps, build the dense reference matrix (which is known-
     * correct from first principles) and pack it into CSR. Internal-
     * complexity trades memory for formula robustness: the sparse
     * character-projector formula for 2D was subtly non-Hermitian in
     * my first-pass derivation, and going through the dense builder
     * sidesteps that entirely. N ≤ 24 is the dense-builder cap. */
    if (d_mu == 2) {
        const irrep_space_group_t *Gcheck = irrep_sg_rep_table_space_group(T);
        if (!Gcheck) return NULL;
        int num_sites = irrep_space_group_num_sites(Gcheck);
        if (num_sites > 24) return NULL;

        long long n_reps_all = irrep_sg_rep_table_count(T);
        double _Complex *Hdense = malloc((size_t)n_reps_all * n_reps_all * sizeof(double _Complex));
        if (!Hdense) return NULL;

        int dim_2d = irrep_sg_heisenberg_sector_build_dense(H, T, lg, mu_k,
                                                             Hdense, (int)n_reps_all);
        if (dim_2d < 0) {
            free(Hdense);
            return NULL;
        }

        /* Pack dense into CSR. Every non-zero entry in the top-left
         * dim_2d × dim_2d block becomes an entry. */
        irrep_sg_heisenberg_sector_t *S = calloc(1, sizeof(*S));
        if (!S) { free(Hdense); return NULL; }
        S->n_reps = dim_2d;
        S->diag   = calloc((size_t)(dim_2d ? dim_2d : 1), sizeof(double));
        S->rowptr = malloc((size_t)(dim_2d + 1) * sizeof(long long));
        if (!S->diag || !S->rowptr) {
            free(Hdense); irrep_sg_heisenberg_sector_free(S);
            return NULL;
        }
        /* Count non-zeros (off-diagonal). */
        long long nnz = 0;
        for (int r = 0; r < dim_2d; ++r) {
            for (int c = 0; c < dim_2d; ++c) {
                if (r == c) continue;
                if (cabs(Hdense[(size_t)r * n_reps_all + c]) > 1e-15)
                    ++nnz;
            }
        }
        S->col       = malloc((size_t)(nnz ? nnz : 1) * sizeof(long long));
        S->coef_cplx = malloc((size_t)(nnz ? nnz : 1) * sizeof(double _Complex));
        if (!S->col || !S->coef_cplx) {
            free(Hdense); irrep_sg_heisenberg_sector_free(S);
            return NULL;
        }
        /* Fill: iterate row (source) k, collect columns (targets). Note
         * the dense matrix is H[row j, col i] = ⟨ũ_j|H|ũ_i⟩; for CSR
         * with source k in outer loop we want entries at (target, k) =
         * M[target, k] = Hdense[target, k]. Column-k of dense matrix. */
        long long pos = 0;
        S->rowptr[0] = 0;
        for (int k_src = 0; k_src < dim_2d; ++k_src) {
            S->diag[k_src] = creal(Hdense[(size_t)k_src * n_reps_all + k_src]);
            for (int k_tgt = 0; k_tgt < dim_2d; ++k_tgt) {
                if (k_tgt == k_src) continue;
                double _Complex m = Hdense[(size_t)k_tgt * n_reps_all + k_src];
                if (cabs(m) > 1e-15) {
                    S->col[pos]       = k_tgt;
                    S->coef_cplx[pos] = m;
                    ++pos;
                }
            }
            S->rowptr[k_src + 1] = pos;
        }
        free(Hdense);
        return S;
    }

    const irrep_space_group_t *G = irrep_sg_rep_table_space_group(T);
    if (!G)
        return NULL;
    int group_order = irrep_space_group_order(G);
    long long n_reps_all = irrep_sg_rep_table_count(T);
    if (n_reps_all <= 0)
        return NULL;

    /* Vend the per-g projector weights w_g once. */
    double _Complex *w = malloc((size_t)group_order * sizeof(double _Complex));
    if (!w)
        return NULL;
    if (irrep_sg_projector_weights(lg, mu_k, w) != 0) {
        free(w);
        return NULL;
    }

    /* Per-rep σ_u = Σ_{s ∈ Stab(u)} conj(w_s). Filter reps with σ_u = 0.
     * Parallel over k — each rep's stabiliser walk is independent. */
    double _Complex *sigma = malloc((size_t)n_reps_all * sizeof(double _Complex));
    long long       *fwd   = malloc((size_t)n_reps_all * sizeof(long long));
    if (!sigma || !fwd) {
        free(w); free(sigma); free(fwd);
        return NULL;
    }
    /* Per-thread stabiliser scratch. A thread whose malloc fails would
     * silently call irrep_sg_stabiliser with out_indices=NULL (which
     * returns 0, producing σ=0 and wrongly filtering the rep out).
     * Detect and propagate the failure so we don't ship corrupt sectors. */
    volatile int alloc_failed = 0;
#pragma omp parallel
    {
        int *stab_local = malloc((size_t)group_order * sizeof(int));
        if (!stab_local)
            alloc_failed = 1;
#pragma omp for schedule(dynamic, 64)
        for (long long k = 0; k < n_reps_all; ++k) {
            if (!stab_local) continue; /* failed-thread no-op; caller aborts after barrier */
            uint64_t u = irrep_sg_rep_table_get(T, k);
            int ns = irrep_sg_stabiliser(G, u, stab_local);
            double _Complex s_sum = 0.0;
            for (int j = 0; j < ns; ++j)
                s_sum += conj(w[stab_local[j]]);
            sigma[k] = s_sum;
        }
        free(stab_local);
    }
    if (alloc_failed) {
        free(w); free(sigma); free(fwd);
        return NULL;
    }

    /* Sequential: build fwd[] (cannot parallelise — serialised filtered index). */
    long long dim = 0;
    for (long long k = 0; k < n_reps_all; ++k) {
        if (cabs(sigma[k]) > 1e-10)
            fwd[k] = dim++;
        else
            fwd[k] = -1;
    }

    irrep_sg_heisenberg_sector_t *S = calloc(1, sizeof(*S));
    if (!S) {
        free(w); free(sigma); free(fwd);
        return NULL;
    }
    S->n_reps = dim;
    S->diag   = calloc((size_t)(dim ? dim : 1), sizeof(double));
    S->rowptr = malloc((size_t)(dim + 1) * sizeof(long long));
    int *counts = calloc((size_t)(dim ? dim : 1), sizeof(int));
    if (!S->diag || !S->rowptr || !counts) {
        free(counts); free(w); free(sigma); free(fwd);
        irrep_sg_heisenberg_sector_free(S);
        return NULL;
    }

    /* Pass 1: count + diag (parallel over filtered reps via fwd[]). */
#pragma omp parallel for schedule(dynamic, 64)
    for (long long k_all = 0; k_all < n_reps_all; ++k_all) {
        long long k_filt = fwd[k_all];
        if (k_filt < 0) continue;

        uint64_t u = irrep_sg_rep_table_get(T, k_all);
        double diag_acc = 0.0;
        int    cnt      = 0;

        for (int b = 0; b < H->num_bonds; ++b) {
            int i = H->bi[b], j = H->bj[b];
            int zi = (int)((u >> i) & 1ULL);
            int zj = (int)((u >> j) & 1ULL);
            diag_acc += ((zi ^ zj) ? -H->coeff_zz[b] : H->coeff_zz[b]);

            if ((zi ^ zj) && H->coeff_pm[b] != 0.0) {
                uint64_t uflip = u ^ ((uint64_t)1 << i) ^ ((uint64_t)1 << j);
                uint64_t v; int g_idx;
                irrep_sg_canonicalise(G, uflip, &v, &g_idx);
                long long kv_orig = irrep_sg_rep_table_index(T, v);
                if (kv_orig >= 0 && fwd[kv_orig] >= 0)
                    ++cnt;
            }
        }
        S->diag[k_filt] = diag_acc;
        counts[k_filt]  = cnt;
    }

    /* Prefix sum. */
    S->rowptr[0] = 0;
    for (long long k = 0; k < dim; ++k)
        S->rowptr[k + 1] = S->rowptr[k] + counts[k];
    long long nnz = S->rowptr[dim];

    S->col       = malloc((size_t)(nnz ? nnz : 1) * sizeof(long long));
    S->coef_cplx = malloc((size_t)(nnz ? nnz : 1) * sizeof(double _Complex));
    if ((!S->col || !S->coef_cplx) && nnz > 0) {
        free(counts); free(w); free(sigma); free(fwd);
        irrep_sg_heisenberg_sector_free(S);
        return NULL;
    }

    /* Pass 2: fill CSR (parallel over filtered reps). */
#pragma omp parallel for schedule(dynamic, 64)
    for (long long k_all = 0; k_all < n_reps_all; ++k_all) {
        long long k_filt = fwd[k_all];
        if (k_filt < 0) continue;

        uint64_t u = irrep_sg_rep_table_get(T, k_all);
        double _Complex sigma_u = sigma[k_all];
        double _Complex sqrt_sigma_u = csqrt(sigma_u);
        long long pos = S->rowptr[k_filt];

        for (int b = 0; b < H->num_bonds; ++b) {
            int i = H->bi[b], j = H->bj[b];
            int zi = (int)((u >> i) & 1ULL);
            int zj = (int)((u >> j) & 1ULL);
            if ((zi ^ zj) && H->coeff_pm[b] != 0.0) {
                uint64_t uflip = u ^ ((uint64_t)1 << i) ^ ((uint64_t)1 << j);
                uint64_t v; int g_idx;
                irrep_sg_canonicalise(G, uflip, &v, &g_idx);
                long long kv_orig = irrep_sg_rep_table_index(T, v);
                if (kv_orig < 0) continue;
                long long kv_filt = fwd[kv_orig];
                if (kv_filt < 0) continue;

                double _Complex sigma_v = sigma[kv_orig];
                double _Complex coef = H->coeff_pm[b] * (double)group_order
                                       * conj(w[g_idx]) * csqrt(sigma_v) / sqrt_sigma_u;
                S->col[pos]       = kv_filt;
                S->coef_cplx[pos] = coef;
                ++pos;
            }
        }
    }

    free(counts); free(w); free(sigma); free(fwd);
    return S;
}

long long irrep_sg_heisenberg_sector_dim(const irrep_sg_heisenberg_sector_t *S) {
    return S ? S->n_reps : 0;
}

int irrep_sg_heisenberg_sector_build_dense(const irrep_heisenberg_t *H,
                                           const struct irrep_sg_rep_table *T,
                                           const struct irrep_sg_little_group *lg,
                                           const struct irrep_sg_little_group_irrep *mu_k,
                                           double _Complex *H_out, int n_max) {
    if (!H || !T || !lg || !mu_k || !H_out || n_max <= 0)
        return -1;
    const irrep_space_group_t *G = irrep_sg_rep_table_space_group(T);
    if (!G)
        return -1;
    int num_sites = irrep_space_group_num_sites(G);
    if (num_sites <= 0 || num_sites > 24)
        return -1; /* 2^24 × 16 bytes ≈ 256 MB per dense vector — practical cap */
    long long D = 1LL << num_sites;
    int group_order = irrep_space_group_order(G);
    long long n_reps = irrep_sg_rep_table_count(T);
    if (n_reps <= 0)
        return -1;

    /* Projector weights w_g for composite (k, μ_k). */
    double _Complex *w = malloc((size_t)group_order * sizeof(double _Complex));
    if (!w)
        return -1;
    if (irrep_sg_projector_weights(lg, mu_k, w) != 0) {
        free(w);
        return -1;
    }

    /* Allocate a single 2^N dense vector scratch + per-rep |ũ⟩ storage.
     * For each kept rep we store the full dense basis vector so the
     * sandwich step can compute ⟨ũ_j|H|ũ_i⟩ without recomputing. */
    double _Complex *basis    = calloc((size_t)n_reps * D, sizeof(double _Complex));
    double _Complex *Hb       = malloc((size_t)D * sizeof(double _Complex));
    int             *keep     = calloc((size_t)n_reps, sizeof(int));
    if (!basis || !Hb || !keep) {
        free(w); free(basis); free(Hb); free(keep);
        return -1;
    }

    /* Pass 1: build |ũ⟩ for each rep. Filter out orbits where σ_u = 0. */
    int dim = 0;
    for (long long k = 0; k < n_reps; ++k) {
        uint64_t u = irrep_sg_rep_table_get(T, k);
        double _Complex *vec = basis + k * D;
        /* vec[x] = Σ_{g: g·u = x} w_g */
        for (int g = 0; g < group_order; ++g) {
            uint64_t gx = irrep_space_group_apply_bits(G, g, u);
            vec[(long long)gx] += w[g];
        }
        /* Norm² of raw projected vector. */
        double norm2 = 0.0;
        for (long long s = 0; s < D; ++s)
            norm2 += creal(vec[s] * conj(vec[s]));
        if (norm2 < 1e-20) {
            /* Null projection — orbit annihilated by μ_k. */
            memset(vec, 0, (size_t)D * sizeof(double _Complex));
            continue;
        }
        double inv_norm = 1.0 / sqrt(norm2);
        for (long long s = 0; s < D; ++s)
            vec[s] *= inv_norm;
        keep[k] = 1;
        ++dim;
    }

    if (dim > n_max) {
        free(w); free(basis); free(Hb); free(keep);
        return -1;
    }

    /* Compact kept basis into contiguous row-block [0, dim). */
    long long write = 0;
    long long *kept_orig = malloc((size_t)(dim ? dim : 1) * sizeof(long long));
    if (!kept_orig) {
        free(w); free(basis); free(Hb); free(keep);
        return -1;
    }
    for (long long k = 0; k < n_reps; ++k) {
        if (!keep[k]) continue;
        if (write != k) {
            memmove(basis + write * D, basis + k * D, (size_t)D * sizeof(double _Complex));
        }
        kept_orig[write] = k;
        ++write;
    }

    /* Pass 2: sandwich. H_out[j, i] = ⟨ũ_j | H | ũ_i⟩.
     * Compute column-by-column: y = H · |ũ_i⟩, then H_out[j, i] = ⟨ũ_j | y⟩. */
    for (int i = 0; i < dim; ++i) {
        const double _Complex *vec_i = basis + (long long)i * D;
        irrep_heisenberg_apply(vec_i, Hb, (void *)H);
        for (int j = 0; j < dim; ++j) {
            const double _Complex *vec_j = basis + (long long)j * D;
            double _Complex acc = 0.0;
            for (long long s = 0; s < D; ++s)
                acc += conj(vec_j[s]) * Hb[s];
            H_out[(size_t)j * n_max + i] = acc;
        }
    }

    free(kept_orig);
    free(w); free(basis); free(Hb); free(keep);
    return dim;
}

/* On-disk sector-binding format:
 *   magic      : 8 bytes "IRREP_SB"
 *   version    : uint32 (currently 1)
 *   flavour    : uint32 (0 = real/trivial, 1 = complex/(k,μ_k))
 *   n_reps     : int64
 *   nnz        : int64
 *   diag       : n_reps × double
 *   rowptr     : (n_reps+1) × int64
 *   col        : nnz × int64
 *   coef_real  : nnz × double        (flavour == 0)
 *   coef_cplx  : nnz × double _Complex (flavour == 1) */

static const char SECTOR_BIN_MAGIC[8]   = {'I','R','R','E','P','_','S','B'};
static const uint32_t SECTOR_BIN_VERSION = 1;

int irrep_sg_heisenberg_sector_save(const irrep_sg_heisenberg_sector_t *S,
                                    const char *path) {
    if (!S || !path)
        return -1;
    FILE *f = fopen(path, "wb");
    if (!f)
        return -1;

    uint32_t flav = S->coef_cplx ? 1u : 0u;
    long long n   = S->n_reps;
    long long nnz = (n > 0) ? S->rowptr[n] : 0;

    if (fwrite(SECTOR_BIN_MAGIC, 1, 8, f) != 8) goto fail;
    if (fwrite(&SECTOR_BIN_VERSION, sizeof(uint32_t), 1, f) != 1) goto fail;
    if (fwrite(&flav, sizeof(uint32_t), 1, f) != 1) goto fail;
    if (fwrite(&n,    sizeof(long long), 1, f) != 1) goto fail;
    if (fwrite(&nnz,  sizeof(long long), 1, f) != 1) goto fail;

    if (n > 0) {
        if ((long long)fwrite(S->diag,   sizeof(double), (size_t)n, f) != n) goto fail;
        if ((long long)fwrite(S->rowptr, sizeof(long long), (size_t)(n + 1), f) != n + 1) goto fail;
        if (nnz > 0) {
            if ((long long)fwrite(S->col, sizeof(long long), (size_t)nnz, f) != nnz) goto fail;
            if (flav == 0) {
                if ((long long)fwrite(S->coef_real, sizeof(double), (size_t)nnz, f) != nnz) goto fail;
            } else {
                if ((long long)fwrite(S->coef_cplx, sizeof(double _Complex), (size_t)nnz, f) != nnz) goto fail;
            }
        }
    }
    fclose(f);
    return 0;
fail:
    fclose(f);
    return -1;
}

irrep_sg_heisenberg_sector_t *
irrep_sg_heisenberg_sector_load(const char *path, long long expected_dim) {
    if (!path)
        return NULL;
    FILE *f = fopen(path, "rb");
    if (!f)
        return NULL;

    char magic[8];
    uint32_t ver, flav;
    long long n, nnz;

    if (fread(magic, 1, 8, f) != 8) goto fail;
    if (memcmp(magic, SECTOR_BIN_MAGIC, 8) != 0) goto fail;
    if (fread(&ver,  sizeof(uint32_t), 1, f) != 1) goto fail;
    if (ver != SECTOR_BIN_VERSION) goto fail;
    if (fread(&flav, sizeof(uint32_t), 1, f) != 1) goto fail;
    if (flav > 1) goto fail;
    if (fread(&n,   sizeof(long long), 1, f) != 1) goto fail;
    if (fread(&nnz, sizeof(long long), 1, f) != 1) goto fail;
    if (n < 0 || nnz < 0) goto fail;
    if (expected_dim >= 0 && expected_dim != n) goto fail;

    irrep_sg_heisenberg_sector_t *S = calloc(1, sizeof(*S));
    if (!S) goto fail;
    S->n_reps = n;
    if (n > 0) {
        S->diag   = malloc((size_t)n * sizeof(double));
        S->rowptr = malloc((size_t)(n + 1) * sizeof(long long));
        if (!S->diag || !S->rowptr) {
            irrep_sg_heisenberg_sector_free(S);
            goto fail;
        }
        if ((long long)fread(S->diag, sizeof(double), (size_t)n, f) != n) {
            irrep_sg_heisenberg_sector_free(S); goto fail;
        }
        if ((long long)fread(S->rowptr, sizeof(long long), (size_t)(n + 1), f) != n + 1) {
            irrep_sg_heisenberg_sector_free(S); goto fail;
        }
        /* Paranoia: the rowptr tail must equal the declared nnz. */
        if (S->rowptr[n] != nnz) {
            irrep_sg_heisenberg_sector_free(S); goto fail;
        }
        if (nnz > 0) {
            S->col = malloc((size_t)nnz * sizeof(long long));
            if (!S->col) { irrep_sg_heisenberg_sector_free(S); goto fail; }
            if ((long long)fread(S->col, sizeof(long long), (size_t)nnz, f) != nnz) {
                irrep_sg_heisenberg_sector_free(S); goto fail;
            }
            if (flav == 0) {
                S->coef_real = malloc((size_t)nnz * sizeof(double));
                if (!S->coef_real) { irrep_sg_heisenberg_sector_free(S); goto fail; }
                if ((long long)fread(S->coef_real, sizeof(double), (size_t)nnz, f) != nnz) {
                    irrep_sg_heisenberg_sector_free(S); goto fail;
                }
            } else {
                S->coef_cplx = malloc((size_t)nnz * sizeof(double _Complex));
                if (!S->coef_cplx) { irrep_sg_heisenberg_sector_free(S); goto fail; }
                if ((long long)fread(S->coef_cplx, sizeof(double _Complex), (size_t)nnz, f) != nnz) {
                    irrep_sg_heisenberg_sector_free(S); goto fail;
                }
            }
        }
    }
    fclose(f);
    return S;
fail:
    fclose(f);
    return NULL;
}

void irrep_sg_heisenberg_sector_free(irrep_sg_heisenberg_sector_t *S) {
    if (!S)
        return;
    free(S->diag);
    free(S->rowptr);
    free(S->col);
    free(S->coef_real);
    free(S->coef_cplx);
    free(S);
}

void irrep_sg_heisenberg_sector_apply(const double _Complex *psi_in,
                                      double _Complex *psi_out, void *opaque) {
    if (!psi_in || !psi_out || !opaque)
        return;
    const irrep_sg_heisenberg_sector_t *S = (const irrep_sg_heisenberg_sector_t *)opaque;
    long long n_reps = S->n_reps;
    memset(psi_out, 0, (size_t)n_reps * sizeof(double _Complex));

    /* Diagonal: y[k] = diag[k] · x[k] */
    for (long long k = 0; k < n_reps; ++k)
        psi_out[k] = S->diag[k] * psi_in[k];

    /* Off-diagonal CSR matvec. Dispatch on which coef array is populated. */
    if (S->coef_real) {
        for (long long k = 0; k < n_reps; ++k) {
            double _Complex xk = psi_in[k];
            if (xk == 0.0)
                continue;
            long long lo = S->rowptr[k], hi = S->rowptr[k + 1];
            for (long long j = lo; j < hi; ++j)
                psi_out[S->col[j]] += S->coef_real[j] * xk;
        }
    } else if (S->coef_cplx) {
        for (long long k = 0; k < n_reps; ++k) {
            double _Complex xk = psi_in[k];
            if (xk == 0.0)
                continue;
            long long lo = S->rowptr[k], hi = S->rowptr[k + 1];
            for (long long j = lo; j < hi; ++j)
                psi_out[S->col[j]] += S->coef_cplx[j] * xk;
        }
    }
}

void irrep_heisenberg_apply(const double _Complex *psi, double _Complex *out, void *opaque) {
    if (!psi || !out || !opaque)
        return;
    const irrep_heisenberg_t *H = (const irrep_heisenberg_t *)opaque;
    const long long           dim = 1LL << H->num_sites;

    memset(out, 0, (size_t)dim * sizeof(double _Complex));
    for (int b = 0; b < H->num_bonds; ++b) {
        const int       i = H->bi[b], j = H->bj[b];
        const long long mi = 1LL << i, mj = 1LL << j;
        const double    czz = H->coeff_zz[b];
        const double    cpm = H->coeff_pm[b];
        for (long long s = 0; s < dim; ++s) {
            const int    zi = (int)((s >> i) & 1);
            const int    zj = (int)((s >> j) & 1);
            const double zz_sign = (zi ^ zj) ? -1.0 : +1.0;
            out[s] += (czz * zz_sign) * psi[s];
            if (zi != zj && cpm != 0.0) {
                const long long s_flip = s ^ mi ^ mj;
                out[s_flip] += cpm * psi[s];
            }
        }
    }
}
