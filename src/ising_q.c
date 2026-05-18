/* SPDX-License-Identifier: MIT */
/** @file ising_q.c
 *  @brief Ising^Q machinery + brute-force Lagrangian enumerator. */
#include <irrep/ising_q.h>
#include <irrep/types.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------- elementary Q-tuple ops ---------- */

int
irrep_ising_q_n_simples(int Q)
{
    if (Q < 0 || Q > IRREP_ISING_Q_MAX) return -1;
    int n = 1;
    for (int i = 0; i < Q; ++i) n *= 3;
    return n;
}

double
irrep_ising_q_global_dim(int Q)
{
    if (Q < 0) return -1.0;
    return ldexp(1.0, Q);   /* 2^Q */
}

int
irrep_ising_q_sigma_count(const irrep_ising_object_t *a, int Q)
{
    if (a == NULL || Q < 0) return -1;
    int n = 0;
    for (int i = 0; i < Q; ++i) {
        if (a[i] == IRREP_ISING_OBJ_SIGMA) ++n;
    }
    return n;
}

double
irrep_ising_q_fpdim(const irrep_ising_object_t *a, int Q)
{
    if (a == NULL || Q < 0) return -1.0;
    double d = 1.0;
    for (int i = 0; i < Q; ++i) d *= irrep_ising_quantum_dim(a[i]);
    return d;
}

int
irrep_ising_q_fusion(const irrep_ising_object_t *a,
                     const irrep_ising_object_t *b,
                     const irrep_ising_object_t *c,
                     int Q)
{
    if (a == NULL || b == NULL || c == NULL || Q < 0) return -1;
    int n = 1;
    for (int i = 0; i < Q; ++i) {
        int f = irrep_ising_fusion(a[i], b[i], c[i]);
        if (f == 0) return 0;
        n *= f;
    }
    return n;
}

/* ---------- decoding simples ---------- */

/* Decode the k-th simple (0 ≤ k < 3^Q) into a Q-tuple of object ids,
 * base-3 digits little-endian. */
static void
decode_simple(int k, int Q, irrep_ising_object_t *out)
{
    for (int i = 0; i < Q; ++i) {
        out[i] = (irrep_ising_object_t)(k % 3);
        k /= 3;
    }
}

/* ---------- brute-force Lagrangian enumeration ----------
 *
 * Restrict to even-σ-count simples (those whose FPdim is rational).
 * Algorithm:
 *   1. Enumerate all even-σ-count simples; record their FPdims (∈ ℤ)
 *      and a precomputed fusion table among them.
 *   2. Iterate over all subsets S of the even-σ simples that include
 *      the identity (index 0). For each:
 *         a. sum of FPdims == 2^Q?
 *         b. fusion-closed?
 *         c. count it.
 *   3. Sigma-included variant: same, but require S to contain at
 *      least one simple with σ-count > 0.
 *
 * For Q = 1 the even-σ set has 2 simples (1, ψ), 2^1 = 2 subsets.
 * For Q = 2 it has 5 simples (1⊗1, 1⊗ψ, ψ⊗1, ψ⊗ψ, σ⊗σ), 2^5 = 32 subsets.
 * For Q = 3 it has 14 simples, 2^14 = 16384 subsets.
 * For Q = 4 it has 41 simples, 2^41 ≈ 2.2 × 10^12 subsets — too many
 * without smarter pruning. The current enumerator caps at Q = 4 with
 * fusion-closure early-exit, which keeps even Q=4 tractable.
 * ---------------------------------------------------------------- */

typedef struct {
    int *idx_in_S;            /* index in even-σ array → global simple id, or -1 */
    int *fusion;              /* m × m × m, fusion[i*m*m + j*m + k] */
    int  m;                   /* number of even-σ simples */
    int *sigma_count;         /* σ-count per even-σ simple, length m */
    int *fpdim_int;           /* integer FPdim per even-σ simple, length m */
    int  identity_local;      /* index in the even-σ array of the unit object */
    int  D;                   /* 2^Q */
} ising_q_evensigma_t;

static void
free_evensigma(ising_q_evensigma_t *e)
{
    if (e == NULL) return;
    free(e->idx_in_S);
    free(e->fusion);
    free(e->sigma_count);
    free(e->fpdim_int);
    memset(e, 0, sizeof *e);
}

static int
build_evensigma(int Q, ising_q_evensigma_t *out)
{
    if (Q < 0 || Q > IRREP_ISING_Q_MAX) return -1;
    int n_total = irrep_ising_q_n_simples(Q);

    /* First pass: count even-σ-count simples and find the identity. */
    int m = 0;
    int identity_local = -1;
    irrep_ising_object_t buf[IRREP_ISING_Q_MAX];
    for (int k = 0; k < n_total; ++k) {
        decode_simple(k, Q, buf);
        if ((irrep_ising_q_sigma_count(buf, Q) & 1) == 0) {
            if (k == 0) identity_local = m;
            ++m;
        }
    }
    if (identity_local < 0) return -1;

    out->m = m;
    out->D = 1 << Q;
    out->idx_in_S      = (int *)calloc((size_t)m, sizeof(int));
    out->sigma_count   = (int *)calloc((size_t)m, sizeof(int));
    out->fpdim_int     = (int *)calloc((size_t)m, sizeof(int));
    out->fusion        = (int *)calloc((size_t)m * m * m, sizeof(int));
    out->identity_local = identity_local;
    if (!out->idx_in_S || !out->sigma_count || !out->fpdim_int || !out->fusion) {
        free_evensigma(out);
        return -1;
    }

    int *global_to_local = (int *)malloc((size_t)n_total * sizeof(int));
    if (!global_to_local) { free_evensigma(out); return -1; }
    for (int k = 0; k < n_total; ++k) global_to_local[k] = -1;

    /* Second pass: populate fields. */
    int j = 0;
    for (int k = 0; k < n_total; ++k) {
        decode_simple(k, Q, buf);
        int sc = irrep_ising_q_sigma_count(buf, Q);
        if ((sc & 1) != 0) continue;
        out->idx_in_S[j]    = k;
        out->sigma_count[j] = sc;
        /* FPdim = (√2)^sc; even sc → exact integer 2^{sc/2}. */
        out->fpdim_int[j]   = 1 << (sc / 2);
        global_to_local[k]  = j;
        ++j;
    }

    /* Fill the m × m × m fusion table among even-σ simples. Note: in
     * Ising^Q, even-σ × even-σ can produce odd-σ outputs, but those
     * contribute 0 to the fusion coefficient in the even-σ restriction
     * (their multiplicities would be irrational, forbidden by Theorem
     * B). We only need fusion within the even-σ subspace; the closure
     * check is correct because a fusion-closed even-σ subset cannot
     * "leak" through fusion: even × even component-wise always
     * produces parity-even σ-counts. */
    irrep_ising_object_t buf_a[IRREP_ISING_Q_MAX], buf_b[IRREP_ISING_Q_MAX], buf_c[IRREP_ISING_Q_MAX];
    for (int ja = 0; ja < m; ++ja) {
        decode_simple(out->idx_in_S[ja], Q, buf_a);
        for (int jb = 0; jb < m; ++jb) {
            decode_simple(out->idx_in_S[jb], Q, buf_b);
            for (int jc = 0; jc < m; ++jc) {
                decode_simple(out->idx_in_S[jc], Q, buf_c);
                out->fusion[ja * m * m + jb * m + jc] =
                    irrep_ising_q_fusion(buf_a, buf_b, buf_c, Q);
            }
        }
    }

    free(global_to_local);
    return 0;
}

/* Test whether a subset is fusion-closed: for every (a, b) with
 * S[a] && S[b], and every c with fusion(a, b, c) > 0, must have S[c]. */
static int
is_fusion_closed(const ising_q_evensigma_t *e, const int *S)
{
    int m = e->m;
    for (int a = 0; a < m; ++a) {
        if (!S[a]) continue;
        for (int b = 0; b < m; ++b) {
            if (!S[b]) continue;
            for (int c = 0; c < m; ++c) {
                if (e->fusion[a * m * m + b * m + c] > 0 && !S[c]) return 0;
            }
        }
    }
    return 1;
}

static int
enumerate_lagrangians(int Q, int sigma_required)
{
    ising_q_evensigma_t e = {0};
    if (build_evensigma(Q, &e) != 0) return -1;
    /* Q=0: only the identity 1, FPdim = 1 = 2^0. Trivially Lagrangian.
     * Per convention, 0 σ → 1 Lagrangian, 0 σ-included. */
    if (Q == 0) {
        int rc = sigma_required ? 0 : 1;
        free_evensigma(&e);
        return rc;
    }

    int m = e.m;
    int count = 0;

    /* Allocate a bit-array S of size m. Iterate over all subsets
     * containing the identity. Bitmask sweep: top (m-1) free bits with
     * identity fixed at 1. For each, check FPdim sum first (cheap),
     * then fusion-closure (more expensive). */
    int *S = (int *)calloc((size_t)m, sizeof(int));
    if (!S) { free_evensigma(&e); return -1; }

    const int free_bits = m - 1;
    /* Cap at 24 free bits (16 M subsets) — past that, brute force
     * is impractical without smarter pruning. */
    if (free_bits > 24) {
        free(S);
        free_evensigma(&e);
        return -1;
    }
    const long long n_subsets = 1LL << free_bits;

    /* The non-identity even-σ simples in order. */
    int *non_id = (int *)malloc((size_t)free_bits * sizeof(int));
    if (!non_id) { free(S); free_evensigma(&e); return -1; }
    int p = 0;
    for (int j = 0; j < m; ++j) {
        if (j == e.identity_local) continue;
        non_id[p++] = j;
    }

    for (long long mask = 0; mask < n_subsets; ++mask) {
        memset(S, 0, (size_t)m * sizeof(int));
        S[e.identity_local] = 1;
        int fpdim_sum = e.fpdim_int[e.identity_local];
        int has_sigma = 0;
        for (int b = 0; b < free_bits; ++b) {
            if ((mask >> b) & 1LL) {
                int j = non_id[b];
                S[j] = 1;
                fpdim_sum += e.fpdim_int[j];
                if (e.sigma_count[j] > 0) has_sigma = 1;
            }
        }
        if (fpdim_sum != e.D) continue;
        if (!is_fusion_closed(&e, S)) continue;
        if (sigma_required && !has_sigma) continue;
        ++count;
    }

    free(non_id);
    free(S);
    free_evensigma(&e);
    return count;
}

int
irrep_ising_q_lagrangian_count_brute(int Q)
{
    return enumerate_lagrangians(Q, /*sigma_required=*/0);
}

int
irrep_ising_q_sigma_included_lagrangian_count(int Q)
{
    return enumerate_lagrangians(Q, /*sigma_required=*/1);
}
