/* SPDX-License-Identifier: MIT */
/* S_N primitives: factorial, permutation sign, full-permutation enumeration,
 * Young-tableaux hook-length formula, and (anti)symmetric projectors on
 * tensor-factored wavefunctions.
 *
 * References:
 *   James & Kerber, "The Representation Theory of the Symmetric Group",
 *     Addison-Wesley 1981, §2 for the hook-length formula.
 *   Fulton, "Young Tableaux", Cambridge 1997, §4 for the partition /
 *     tableau combinatorics.
 */

#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/sym_group.h>

extern void irrep_set_error_(const char *fmt, ...);

/* -------------------------------------------------------------------------- *
 * Factorial                                                                  *
 * -------------------------------------------------------------------------- */

long long irrep_factorial(int n) {
    if (n < 0 || n > 20)
        return -1;
    long long r = 1;
    for (int k = 2; k <= n; ++k)
        r *= k;
    return r;
}

/* -------------------------------------------------------------------------- *
 * Permutation sign                                                           *
 * -------------------------------------------------------------------------- */

int irrep_permutation_sign(const int *perm, int n) {
    if (!perm || n < 1)
        return 0;
    /* Validate: each of [0, n) exactly once. */
    int seen[32] = {0};
    if (n > 32)
        return 0;
    for (int i = 0; i < n; ++i) {
        if (perm[i] < 0 || perm[i] >= n || seen[perm[i]])
            return 0;
        seen[perm[i]] = 1;
    }
    /* Count inversions via simple O(n²). */
    int inv = 0;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (perm[i] > perm[j])
                ++inv;
        }
    }
    return (inv & 1) ? -1 : +1;
}

/* -------------------------------------------------------------------------- *
 * Lexicographic permutation enumeration                                      *
 * -------------------------------------------------------------------------- */

static int next_lex_perm_(int *a, int n) {
    int i = n - 2;
    while (i >= 0 && a[i] >= a[i + 1])
        --i;
    if (i < 0)
        return 0;
    int j = n - 1;
    while (a[j] <= a[i])
        --j;
    {
        int t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
    /* reverse a[i+1..n-1] */
    int lo = i + 1, hi = n - 1;
    while (lo < hi) {
        int t = a[lo];
        a[lo] = a[hi];
        a[hi] = t;
        ++lo;
        --hi;
    }
    return 1;
}

irrep_status_t irrep_permutations_all(int n, int *out_perms) {
    if (n < 1 || n > 10 || !out_perms)
        return IRREP_ERR_INVALID_ARG;
    long long fact = irrep_factorial(n);
    int       a[16];
    for (int i = 0; i < n; ++i)
        a[i] = i;
    long long r = 0;
    do {
        memcpy(out_perms + r * n, a, (size_t)n * sizeof(int));
        ++r;
    } while (next_lex_perm_(a, n));
    if (r != fact)
        return IRREP_ERR_PRECONDITION;
    return IRREP_OK;
}

/* -------------------------------------------------------------------------- *
 * Hook-length dimension formula                                              *
 * -------------------------------------------------------------------------- */

long long irrep_young_dim(const int *partition, int n_parts) {
    if (!partition || n_parts < 1 || n_parts > 32)
        return -1;
    int N = 0;
    for (int i = 0; i < n_parts; ++i) {
        if (partition[i] <= 0)
            return -1;
        if (i > 0 && partition[i] > partition[i - 1])
            return -1;
        N += partition[i];
    }
    if (N > 20)
        return -1; /* factorial overflow */

    /* Compute hook length h(i, j) for each cell (i, j) in the Young diagram:
     *   h(i, j) = (# cells to the right in row i) + (# cells below in col j) + 1
     */
    long long hook_prod = 1;
    for (int i = 0; i < n_parts; ++i) {
        int row_len = partition[i];
        for (int j = 0; j < row_len; ++j) {
            int arm = row_len - j - 1;
            int leg = 0;
            for (int ii = i + 1; ii < n_parts; ++ii) {
                if (partition[ii] > j)
                    ++leg;
                else
                    break;
            }
            long long h = arm + leg + 1;
            hook_prod *= h;
            if (hook_prod <= 0)
                return -1; /* overflow */
        }
    }

    long long Nfact = irrep_factorial(N);
    if (Nfact < 0 || hook_prod == 0)
        return -1;
    return Nfact / hook_prod;
}

/* -------------------------------------------------------------------------- *
 * (Anti)symmetric projectors on tensor-factored states.                      *
 *                                                                            *
 * For a state ψ[i_0, i_1, …, i_{N-1}] with i = Σ i_k · d^k, the permutation  *
 * σ acts as                                                                  *
 *   (σ · ψ)(i_0, …, i_{N-1}) = ψ(i_{σ(0)}, …, i_{σ(N-1)}).                   *
 * The projector is (1/N!) Σ_σ [sign(σ)] · (σ · ψ).                           *
 * -------------------------------------------------------------------------- */

static long long ipow_ll_(int base, int exp) {
    long long r = 1;
    for (int k = 0; k < exp; ++k)
        r *= (long long)base;
    return r;
}

/* Permute the digits of `i` by σ: new_digits[k] = old_digits[σ(k)]. */
static long long permute_index_(long long i, int N, int d, const int *sigma,
                                const long long *weight) {
    int digits[16];
    for (int k = 0; k < N; ++k) {
        digits[k] = (int)(i % d);
        i /= d;
    }
    long long out = 0;
    for (int k = 0; k < N; ++k) {
        out += (long long)digits[sigma[k]] * weight[k];
    }
    return out;
}

static irrep_status_t apply_projector_(int N, int d, int use_sign, const double _Complex *psi_in,
                                       double _Complex *psi_out) {
    if (N < 1 || N > 10 || d < 2 || !psi_in || !psi_out)
        return IRREP_ERR_INVALID_ARG;

    long long fact = irrep_factorial(N);
    long long dim = ipow_ll_(d, N);

    int      *perms = malloc((size_t)fact * (size_t)N * sizeof(int));
    if (!perms)
        return IRREP_ERR_OUT_OF_MEMORY;
    if (irrep_permutations_all(N, perms) != IRREP_OK) {
        free(perms);
        return IRREP_ERR_PRECONDITION;
    }

    long long weight[16];
    for (int k = 0; k < N; ++k)
        weight[k] = ipow_ll_(d, k);

    /* Alias-safe: copy if needed. */
    const double _Complex *src = psi_in;
    double _Complex       *psi_in_copy = NULL;
    if (psi_in == psi_out) {
        psi_in_copy = malloc((size_t)dim * sizeof(double _Complex));
        if (!psi_in_copy) {
            free(perms);
            return IRREP_ERR_OUT_OF_MEMORY;
        }
        memcpy(psi_in_copy, psi_in, (size_t)dim * sizeof(double _Complex));
        src = psi_in_copy;
    }

    memset(psi_out, 0, (size_t)dim * sizeof(double _Complex));

    for (long long s = 0; s < fact; ++s) {
        const int *sigma = perms + s * N;
        int        sign = use_sign ? irrep_permutation_sign(sigma, N) : 1;
        double _Complex coef = ((double)sign) / (double)fact;
        for (long long i = 0; i < dim; ++i) {
            long long j = permute_index_(i, N, d, sigma, weight);
            psi_out[i] += coef * src[j];
        }
    }

    free(perms);
    free(psi_in_copy);
    return IRREP_OK;
}

irrep_status_t irrep_sym_group_antisymmetrize(int N, int local_dim, const double _Complex *psi_in,
                                              double _Complex *psi_out) {
    return apply_projector_(N, local_dim, /*use_sign=*/1, psi_in, psi_out);
}

irrep_status_t irrep_sym_group_symmetrize(int N, int local_dim, const double _Complex *psi_in,
                                          double _Complex *psi_out) {
    return apply_projector_(N, local_dim, /*use_sign=*/0, psi_in, psi_out);
}
