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


/* Backtracking-with-incremental-closure enumerator.
 *
 * State: a partial subset S (with current FPdim sum and σ-included
 * flag) and a sweep cursor `idx` into the array of non-identity even-σ
 * simples. At each cursor position we choose either to skip the simple
 * or to include it (along with all forced-by-fusion simples not yet
 * in S).
 *
 * Pruning:
 *   - FPdim budget: if the current sum > D, abort the branch.
 *   - Forced-inclusion budget: when including a simple, compute the
 *     transitive closure under fusion with the current S; abort if the
 *     closure exceeds D.
 *
 * Worst case is still O(2^m), but the FPdim and closure pruning
 * eliminate the vast majority of subsets in practice. Q=4 (41 simples,
 * 2^40 raw) runs in seconds with this scheme; Q=5 (122 simples, 2^121
 * raw) remains research-track. */

/* Compute the fusion-closure of S under the fusion table: extend S
 * by adding every simple c with N^{ab}_c > 0 for some a, b ∈ S.
 * Returns the new total FPdim sum, or -1 if the closure overflows
 * the target D or tries to add a simple marked in `forbidden`.
 *
 * The `forbidden` mask is the canonicalisation device: a simple gets
 * marked forbidden when the recursion's SKIP branch at that idx
 * deliberately excludes it. If closure later tries to add a forbidden
 * simple, the current branch is rejected — that subset belongs to a
 * different recursion path. This guarantees each Lagrangian is counted
 * exactly once. */
static int
fusion_close(const ising_q_evensigma_t *e, int *S, int fpdim_sum, int D,
             const int *forbidden)
{
    int m = e->m;
    int changed = 1;
    while (changed) {
        changed = 0;
        for (int a = 0; a < m; ++a) {
            if (!S[a]) continue;
            for (int b = 0; b < m; ++b) {
                if (!S[b]) continue;
                for (int c = 0; c < m; ++c) {
                    if (e->fusion[a * m * m + b * m + c] > 0 && !S[c]) {
                        if (forbidden && forbidden[c]) return -1;
                        S[c] = 1;
                        fpdim_sum += e->fpdim_int[c];
                        changed = 1;
                        if (fpdim_sum > D) return -1;
                    }
                }
            }
        }
    }
    return fpdim_sum;
}

/* Snapshot S → buffer for backtracking. */
static void
snapshot_S(const int *S, int *buf, int m)
{ memcpy(buf, S, (size_t)m * sizeof(int)); }

/* Recursive enumerator. Modifies `*count` in place.
 *
 * The `forbidden` array tracks simples that were explicitly SKIPped
 * earlier in the path. Fusion-closure refuses to add forbidden
 * simples — this is the canonicalisation step that ensures each
 * Lagrangian is reached by exactly one path. */
static void
enum_recurse(const ising_q_evensigma_t *e, int *S, int *forbidden,
             int fpdim_sum, int idx, int sigma_required,
             int has_sigma, int *count)
{
    int m = e->m;
    int D = e->D;
    if (idx >= m) {
        if (fpdim_sum != D) return;
        if (sigma_required && !has_sigma) return;
        ++(*count);
        return;
    }
    if (idx == e->identity_local) {
        /* Identity already in S; advance. */
        enum_recurse(e, S, forbidden, fpdim_sum, idx + 1, sigma_required,
                     has_sigma, count);
        return;
    }

    if (S[idx]) {
        /* Already in S (forced by an earlier closure). No SKIP/INCLUDE
         * choice here — just advance. */
        enum_recurse(e, S, forbidden, fpdim_sum, idx + 1, sigma_required,
                     has_sigma, count);
        return;
    }

    /* Branch 1: SKIP simple `idx`. Mark it forbidden so closure can't
     * later add it back via fusion of larger-idx simples. */
    forbidden[idx] = 1;
    enum_recurse(e, S, forbidden, fpdim_sum, idx + 1, sigma_required,
                 has_sigma, count);
    forbidden[idx] = 0;  /* backtrack */

    /* Branch 2: INCLUDE simple `idx`. */
    if (fpdim_sum + e->fpdim_int[idx] > D) return;  /* FPdim prune */

    int *snapshot = (int *)malloc((size_t)m * sizeof(int));
    if (!snapshot) return;
    snapshot_S(S, snapshot, m);

    S[idx] = 1;
    int new_fpdim = fpdim_sum + e->fpdim_int[idx];
    int new_has_sigma = has_sigma || (e->sigma_count[idx] > 0);
    new_fpdim = fusion_close(e, S, new_fpdim, D, forbidden);
    if (new_fpdim >= 0) {
        /* Update has_sigma to account for closure-forced σ-bearing
         * inclusions. */
        for (int j = 0; j < m; ++j) {
            if (S[j] && !snapshot[j] && e->sigma_count[j] > 0) {
                new_has_sigma = 1;
            }
        }
        enum_recurse(e, S, forbidden, new_fpdim, idx + 1, sigma_required,
                     new_has_sigma, count);
    }

    /* Backtrack: restore S. */
    snapshot_S(snapshot, S, m);
    free(snapshot);
}

static int
enumerate_lagrangians(int Q, int sigma_required)
{
    ising_q_evensigma_t e = {0};
    if (build_evensigma(Q, &e) != 0) return -1;
    if (Q == 0) {
        int rc = sigma_required ? 0 : 1;
        free_evensigma(&e);
        return rc;
    }

    int m = e.m;
    int *S = (int *)calloc((size_t)m, sizeof(int));
    int *forbidden = (int *)calloc((size_t)m, sizeof(int));
    if (!S || !forbidden) {
        free(S); free(forbidden); free_evensigma(&e);
        return -1;
    }
    S[e.identity_local] = 1;
    int fpdim_sum = e.fpdim_int[e.identity_local];

    /* Safety cap. Backtracking with canonical-forbidden pruning
     * handles much larger m than brute-force bit-mask sweep, but
     * Q=5 (m=122) and beyond remain research-track. */
    if (m > 60) {
        free(S); free(forbidden); free_evensigma(&e);
        return -1;
    }

    int count = 0;
    enum_recurse(&e, S, forbidden, fpdim_sum, /*idx=*/0, sigma_required,
                 /*has_sigma=*/0, &count);

    free(S);
    free(forbidden);
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
