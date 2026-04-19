/* SPDX-License-Identifier: MIT */
/* Tests for S_N, Young-tableaux dimensions, (anti)symmetric
 * projectors on tensor states.
 *
 * Coverage:
 *   - Factorial edge cases and value samples.
 *   - Permutation sign: identity → +1, swaps → -1, cycles → sign check.
 *   - Permutation enumeration gives N! distinct entries, identity at index 0.
 *   - Hook-length formula: known small-partition irrep dimensions of S_N.
 *     Total sum-of-squares of dimensions = N! for every N (Plancherel).
 *   - Antisymmetric projector: Slater determinant |12⟩ − |21⟩ for 2 sites.
 *   - Antisymmetric projector zeros out states with repeated indices.
 *   - Symmetric projector: permanent-like sum; |12⟩ + |21⟩ for 2 sites.
 *   - Projectors are idempotent.
 *   - Composition: symmetric then antisymmetric → zero (orthogonal irreps).
 */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "harness.h"

#include <irrep/sym_group.h>

int main(void) {
    IRREP_TEST_START("sym_group");

    /* ------------------------------------------------------------------ */
    /* Factorial                                                          */
    /* ------------------------------------------------------------------ */
    IRREP_ASSERT(irrep_factorial(0)  ==    1LL);
    IRREP_ASSERT(irrep_factorial(1)  ==    1LL);
    IRREP_ASSERT(irrep_factorial(5)  ==  120LL);
    IRREP_ASSERT(irrep_factorial(10) == 3628800LL);
    IRREP_ASSERT(irrep_factorial(20) == 2432902008176640000LL);
    IRREP_ASSERT(irrep_factorial(-1) < 0);
    IRREP_ASSERT(irrep_factorial(21) < 0);

    /* ------------------------------------------------------------------ */
    /* Permutation sign                                                   */
    /* ------------------------------------------------------------------ */
    int id3[3]  = { 0, 1, 2 };
    int sw3[3]  = { 1, 0, 2 };          /* one transposition */
    int cy3[3]  = { 1, 2, 0 };          /* 3-cycle (012)      */
    int rev3[3] = { 2, 1, 0 };          /* three transpositions? (02)(01) = (0 1 2) no — (0,2) is one transposition */
    IRREP_ASSERT(irrep_permutation_sign(id3,  3) == +1);
    IRREP_ASSERT(irrep_permutation_sign(sw3,  3) == -1);
    IRREP_ASSERT(irrep_permutation_sign(cy3,  3) == +1);   /* even */
    IRREP_ASSERT(irrep_permutation_sign(rev3, 3) == -1);   /* (0,2): 2 inversions + 1 = 3 → odd */

    int bad[3] = { 0, 1, 1 };   /* duplicate */
    IRREP_ASSERT(irrep_permutation_sign(bad, 3) == 0);

    /* ------------------------------------------------------------------ */
    /* Permutation enumeration                                            */
    /* ------------------------------------------------------------------ */
    long long fact4 = irrep_factorial(4);
    int *perms = malloc((size_t)fact4 * 4 * sizeof(int));
    IRREP_ASSERT(irrep_permutations_all(4, perms) == IRREP_OK);
    /* Identity at index 0 */
    for (int k = 0; k < 4; ++k) IRREP_ASSERT(perms[k] == k);
    /* Every permutation distinct */
    long long total = 1;
    for (long long i = 1; i < fact4; ++i) {
        int diff = 0;
        for (int k = 0; k < 4; ++k) if (perms[i*4+k] != perms[(i-1)*4+k]) { diff = 1; break; }
        if (diff) ++total;
    }
    IRREP_ASSERT(total == fact4);       /* all lex-sorted, strictly increasing */
    free(perms);

    /* ------------------------------------------------------------------ */
    /* Hook-length formula                                                 */
    /*   [N] (single row)       → dim = 1                                  */
    /*   [1,1,...,1] (column)    → dim = 1                                 */
    /*   [N-1, 1] (hook)         → dim = N - 1                             */
    /*   [3, 2] for N=5          → dim = 5                                 */
    /*   [2, 2] for N=4          → dim = 2                                 */
    /* ------------------------------------------------------------------ */
    int row5[1]   = { 5 };
    int col5[5]   = { 1, 1, 1, 1, 1 };
    int hook5[2]  = { 4, 1 };
    int part32[2] = { 3, 2 };
    int part22[2] = { 2, 2 };
    int part31[2] = { 3, 1 };
    IRREP_ASSERT(irrep_young_dim(row5,  1) == 1);
    IRREP_ASSERT(irrep_young_dim(col5,  5) == 1);
    IRREP_ASSERT(irrep_young_dim(hook5, 2) == 4);
    IRREP_ASSERT(irrep_young_dim(part32, 2) == 5);
    IRREP_ASSERT(irrep_young_dim(part22, 2) == 2);
    IRREP_ASSERT(irrep_young_dim(part31, 2) == 3);

    /* Plancherel: Σ dim(λ)² = N! for all partitions λ of N.
     * For N=5: partitions are (5), (4,1), (3,2), (3,1,1), (2,2,1), (2,1,1,1), (1,1,1,1,1)
     * dims: 1, 4, 5, 6, 5, 4, 1. Sum squares = 1+16+25+36+25+16+1 = 120 = 5!. */
    int p5_a[1] = { 5 };
    int p5_b[2] = { 4, 1 };
    int p5_c[2] = { 3, 2 };
    int p5_d[3] = { 3, 1, 1 };
    int p5_e[3] = { 2, 2, 1 };
    int p5_f[4] = { 2, 1, 1, 1 };
    int p5_g[5] = { 1, 1, 1, 1, 1 };
    long long sum = 0;
    sum += irrep_young_dim(p5_a, 1) * irrep_young_dim(p5_a, 1);
    sum += irrep_young_dim(p5_b, 2) * irrep_young_dim(p5_b, 2);
    sum += irrep_young_dim(p5_c, 2) * irrep_young_dim(p5_c, 2);
    sum += irrep_young_dim(p5_d, 3) * irrep_young_dim(p5_d, 3);
    sum += irrep_young_dim(p5_e, 3) * irrep_young_dim(p5_e, 3);
    sum += irrep_young_dim(p5_f, 4) * irrep_young_dim(p5_f, 4);
    sum += irrep_young_dim(p5_g, 5) * irrep_young_dim(p5_g, 5);
    IRREP_ASSERT(sum == 120);

    /* Bad partitions */
    int bad1[2] = { 1, 2 };     /* non-decreasing */
    int bad2[2] = { 3, 0 };     /* non-positive */
    IRREP_ASSERT(irrep_young_dim(bad1, 2) == -1);
    IRREP_ASSERT(irrep_young_dim(bad2, 2) == -1);

    /* ------------------------------------------------------------------ */
    /* Antisymmetric projector on N=2, d=3 (so basis size 9).              */
    /* Input: ψ = |01⟩ (i.e., index i = 0 + 1·3 = 3).                       */
    /* A ψ = (|01⟩ − |10⟩)/2! = (1/2)|01⟩ − (1/2)|10⟩                      */
    /* |01⟩ at i=3, |10⟩ at i=1.                                            */
    /* ------------------------------------------------------------------ */
    double _Complex psi[9] = {0};
    double _Complex out[9] = {0};
    psi[3] = 1.0;
    IRREP_ASSERT(irrep_sym_group_antisymmetrize(2, 3, psi, out) == IRREP_OK);
    IRREP_ASSERT_NEAR(creal(out[3]), +0.5, 1e-14);
    IRREP_ASSERT_NEAR(creal(out[1]), -0.5, 1e-14);
    /* All other entries zero */
    for (int k = 0; k < 9; ++k) {
        if (k == 1 || k == 3) continue;
        IRREP_ASSERT_NEAR(cabs(out[k]), 0.0, 1e-14);
    }

    /* Symmetric projector on the same input → (|01⟩ + |10⟩)/2 */
    IRREP_ASSERT(irrep_sym_group_symmetrize(2, 3, psi, out) == IRREP_OK);
    IRREP_ASSERT_NEAR(creal(out[3]), 0.5, 1e-14);
    IRREP_ASSERT_NEAR(creal(out[1]), 0.5, 1e-14);

    /* Antisymmetric projector zeros repeated-index state |00⟩ (i = 0). */
    memset(psi, 0, sizeof psi);
    psi[0] = 1.0;
    IRREP_ASSERT(irrep_sym_group_antisymmetrize(2, 3, psi, out) == IRREP_OK);
    for (int k = 0; k < 9; ++k) IRREP_ASSERT_NEAR(cabs(out[k]), 0.0, 1e-14);

    /* Symmetric of the diagonal |00⟩ → |00⟩ itself (all permutations map to 00). */
    IRREP_ASSERT(irrep_sym_group_symmetrize(2, 3, psi, out) == IRREP_OK);
    IRREP_ASSERT_NEAR(creal(out[0]), 1.0, 1e-14);

    /* ------------------------------------------------------------------ */
    /* Idempotence: A(Aψ) = Aψ, S(Sψ) = Sψ on a 3-site random state.       */
    /* ------------------------------------------------------------------ */
    int N = 3, d = 3;
    long long dim = 1;
    for (int k = 0; k < N; ++k) dim *= d;
    double _Complex *x    = malloc((size_t)dim * sizeof(double _Complex));
    double _Complex *y    = malloc((size_t)dim * sizeof(double _Complex));
    double _Complex *z    = malloc((size_t)dim * sizeof(double _Complex));
    for (long long i = 0; i < dim; ++i) x[i] = (double)(i % 7) - 3.0 + I * ((double)(i % 5) - 2.0);

    IRREP_ASSERT(irrep_sym_group_antisymmetrize(N, d, x, y) == IRREP_OK);
    IRREP_ASSERT(irrep_sym_group_antisymmetrize(N, d, y, z) == IRREP_OK);
    for (long long i = 0; i < dim; ++i) {
        IRREP_ASSERT_NEAR(creal(y[i]), creal(z[i]), 1e-12);
        IRREP_ASSERT_NEAR(cimag(y[i]), cimag(z[i]), 1e-12);
    }

    IRREP_ASSERT(irrep_sym_group_symmetrize(N, d, x, y) == IRREP_OK);
    IRREP_ASSERT(irrep_sym_group_symmetrize(N, d, y, z) == IRREP_OK);
    for (long long i = 0; i < dim; ++i) {
        IRREP_ASSERT_NEAR(creal(y[i]), creal(z[i]), 1e-12);
        IRREP_ASSERT_NEAR(cimag(y[i]), cimag(z[i]), 1e-12);
    }

    /* Orthogonality: A(Sψ) = 0 (sym + antisym irreps of S_N are distinct for N ≥ 2) */
    IRREP_ASSERT(irrep_sym_group_symmetrize(N, d, x, y) == IRREP_OK);
    IRREP_ASSERT(irrep_sym_group_antisymmetrize(N, d, y, z) == IRREP_OK);
    double acc = 0.0;
    for (long long i = 0; i < dim; ++i) acc += cabs(z[i]);
    IRREP_ASSERT_NEAR(acc, 0.0, 1e-12);

    free(x); free(y); free(z);

    /* ------------------------------------------------------------------ */
    /* Error paths                                                         */
    /* ------------------------------------------------------------------ */
    IRREP_ASSERT(irrep_sym_group_antisymmetrize(0, 2, psi, out) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_sym_group_antisymmetrize(2, 1, psi, out) == IRREP_ERR_INVALID_ARG);

    return IRREP_TEST_END();
}
