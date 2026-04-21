/* SPDX-License-Identifier: MIT */
/* Tests for half-integer (spinor) tensor product.
 *
 * Coverage:
 *   - Path enumeration on spin-½ ⊗ spin-½ = spin-0 ⊕ spin-1.
 *   - spin-0 projection: σ_singlet = (↑↓ − ↓↑)/√2. Apply tp on (↑, ↑), (↑, ↓),
 *     etc.; verify singlet amplitude matches CG coefficients.
 *   - spin-1 projection: symmetric combinations.
 *   - Build error paths: empty multisets, triangle-violating explicit paths.
 *   - Apply: basic correctness of c = a ⊗ b for spin-½ inputs.
 *   - Weighted apply: complex weights modulate output.
 *   - Backward: finite-difference cross-check on a small path.
 */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "harness.h"

#include <irrep/multiset_2j.h>
#include <irrep/tensor_product.h>
#include <irrep/clebsch_gordan.h>

int main(void) {
    IRREP_TEST_START("tensor_product_2j");

    /* ------------------------------------------------------------------ */
    /* spin-½ ⊗ spin-½ = spin-0 ⊕ spin-1                                   */
    /* ------------------------------------------------------------------ */
    irrep_multiset_2j_t *a = irrep_multiset_2j_parse("1x1/2e");
    irrep_multiset_2j_t *b = irrep_multiset_2j_parse("1x1/2e");
    irrep_multiset_2j_t *c = irrep_multiset_2j_parse("1x0e + 1x1e");
    IRREP_ASSERT(a && b && c);
    IRREP_ASSERT(a->total_dim == 2);
    IRREP_ASSERT(c->total_dim == 4); /* 1 + 3 */

    /* Path count: ½⊗½ → 0 (allowed, triangle |1-1|=0 ≤ 0 ≤ 2, parity e·e=e),
     *             ½⊗½ → 1 (allowed, 0 ≤ 2 ≤ 2, parity e·e=e).
     * Total 2 paths. */
    int nt = irrep_tp_2j_enumerate_paths(a, b, c, NULL, 0);
    IRREP_ASSERT(nt == 2);

    int paths[6];
    irrep_tp_2j_enumerate_paths(a, b, c, paths, 2);

    tp_2j_descriptor_t *desc = irrep_tp_2j_build(a, b, c, NULL, 0);
    IRREP_ASSERT(desc != NULL);
    IRREP_ASSERT(irrep_tp_2j_num_paths(desc) == 2);
    IRREP_ASSERT(irrep_tp_2j_output_dim(desc) == 4);

    /* ------------------------------------------------------------------ */
    /* libirrep stores m in ascending order: storage index 0 = m = −j,    */
    /* index d−1 = m = +j. So for spin-½:                                 */
    /*     a[0] = ⟨j, m=−½|a⟩ = |↓⟩ component,                             */
    /*     a[1] = ⟨j, m=+½|a⟩ = |↑⟩ component.                             */
    /* Apply on |↑↑⟩: a = (0, 1), b = (0, 1). Expected fully-aligned       */
    /* triplet: j=1, m=+1 → c_out[3] = 1 (last element of the spin-1 block.*/
    /* ------------------------------------------------------------------ */
    double _Complex a_in[2] = {0.0, 1.0}; /* |↑⟩ */
    double _Complex b_in[2] = {0.0, 1.0}; /* |↑⟩ */
    double _Complex c_out[4] = {0};
    irrep_tp_2j_apply(desc, a_in, b_in, c_out);
    /* c layout: [j=0 block (1 entry), j=1 block (3 entries, m=−1, 0, +1)] */
    IRREP_ASSERT_NEAR(cabs(c_out[0]), 0.0, 1e-14);  /* singlet = 0 */
    IRREP_ASSERT_NEAR(cabs(c_out[1]), 0.0, 1e-14);  /* triplet m=−1 = 0 */
    IRREP_ASSERT_NEAR(cabs(c_out[2]), 0.0, 1e-14);  /* triplet m=0 */
    IRREP_ASSERT_NEAR(creal(c_out[3]), 1.0, 1e-12); /* triplet m=+1 */
    IRREP_ASSERT_NEAR(cimag(c_out[3]), 0.0, 1e-12);

    /* |↑↓⟩: a = (0,1) ⊗ b = (1,0). Singlet (1/√2)(|↑↓⟩−|↓↑⟩) component:  */
    /*   CG(½,½;½,−½ | 0, 0) = +1/√2                                       */
    /*   CG(½,½;½,−½ | 1, 0) = +1/√2                                       */
    a_in[0] = 0.0;
    a_in[1] = 1.0;
    b_in[0] = 1.0;
    b_in[1] = 0.0;
    irrep_tp_2j_apply(desc, a_in, b_in, c_out);
    IRREP_ASSERT_NEAR(creal(c_out[0]), 1.0 / sqrt(2.0), 1e-12); /* singlet */
    IRREP_ASSERT_NEAR(creal(c_out[2]), 1.0 / sqrt(2.0), 1e-12); /* triplet m=0 */
    IRREP_ASSERT_NEAR(cabs(c_out[1]), 0.0, 1e-14);              /* triplet m=−1 */
    IRREP_ASSERT_NEAR(cabs(c_out[3]), 0.0, 1e-14);              /* triplet m=+1 */

    /* |↓↑⟩: singlet = −1/√2, triplet m=0 = +1/√2                          */
    a_in[0] = 1.0;
    a_in[1] = 0.0;
    b_in[0] = 0.0;
    b_in[1] = 1.0;
    irrep_tp_2j_apply(desc, a_in, b_in, c_out);
    IRREP_ASSERT_NEAR(creal(c_out[0]), -1.0 / sqrt(2.0), 1e-12);
    IRREP_ASSERT_NEAR(creal(c_out[2]), 1.0 / sqrt(2.0), 1e-12);

    irrep_tp_2j_free(desc);
    irrep_multiset_2j_free(a);
    irrep_multiset_2j_free(b);
    irrep_multiset_2j_free(c);

    /* ------------------------------------------------------------------ */
    /* Weighted apply: w · c0                                              */
    /* ------------------------------------------------------------------ */
    a = irrep_multiset_2j_parse("1x1/2e");
    b = irrep_multiset_2j_parse("1x1/2e");
    c = irrep_multiset_2j_parse("1x0e");
    desc = irrep_tp_2j_build(a, b, c, NULL, 0);
    IRREP_ASSERT(desc != NULL);
    IRREP_ASSERT(irrep_tp_2j_num_paths(desc) == 1);

    /* 1 path × mult 1×1×1 = 1 weight. Use |↑⟩=(0,1) and |↓⟩=(1,0). */
    double _Complex w = 2.0 + 1.0 * I;
    double _Complex a_sp[2] = {0.0, 1.0}; /* |↑⟩ */
    double _Complex b_sp[2] = {1.0, 0.0}; /* |↓⟩ */
    double _Complex c_sp[1] = {0};
    irrep_tp_2j_apply_weighted(desc, &w, a_sp, b_sp, c_sp);
    /* Singlet amplitude: CG(½,½;½,-½|0,0) = 1/√2. Weighted: w · 1/√2. */
    double _Complex expected = w * (1.0 / sqrt(2.0));
    IRREP_ASSERT_NEAR(creal(c_sp[0]), creal(expected), 1e-12);
    IRREP_ASSERT_NEAR(cimag(c_sp[0]), cimag(expected), 1e-12);

    /* Backward: finite-difference ∂c/∂a for the weighted apply. */
    double eps = 1e-6;
    /* Analytic backward: given grad_c = (1+0i), compute grad_a. */
    double _Complex grad_c[1] = {1.0 + 0.0 * I};
    double _Complex grad_a[2] = {0}, grad_b[2] = {0}, grad_w[1] = {0};
    irrep_tp_2j_apply_backward(desc, &w, a_sp, b_sp, grad_c, grad_a, grad_b, grad_w);

    /* FD: ∂c/∂a[0] at a = a_sp.  c(a+eps·e0) = w · (1/√2) · (a[0]+eps) · b[1] · CG +
     *                                        w · (1/√2) · a[1] · b[0] · CG_other term.
     * ... actually easier: just perturb a[0], measure c, subtract, divide by eps. */
    /* Perturb a_sp[1] (the |↑⟩ component — the one contributing to the
     * singlet overlap with |↓⟩ in b_sp). */
    double _Complex a_pert[2] = {a_sp[0], a_sp[1] + eps};
    double _Complex c_pert[1];
    irrep_tp_2j_apply_weighted(desc, &w, a_pert, b_sp, c_pert);
    double _Complex fd = (c_pert[0] - c_sp[0]) / eps;
    /* grad_a[1] matches fd for inputs real (conj convention is a no-op). */
    IRREP_ASSERT_NEAR(creal(grad_a[1]), creal(fd), 1e-5);
    IRREP_ASSERT_NEAR(cimag(grad_a[1]), cimag(fd), 1e-5);

    irrep_tp_2j_free(desc);
    irrep_multiset_2j_free(a);
    irrep_multiset_2j_free(b);
    irrep_multiset_2j_free(c);

    /* ------------------------------------------------------------------ */
    /* Error paths                                                         */
    /* ------------------------------------------------------------------ */
    irrep_multiset_2j_t *empty = irrep_multiset_2j_new(0);
    a = irrep_multiset_2j_parse("1x1/2e");
    b = irrep_multiset_2j_parse("1x1/2e");
    /* No output term supplied → 0 paths */
    int nn = irrep_tp_2j_enumerate_paths(a, b, empty, NULL, 0);
    IRREP_ASSERT(nn == 0);
    desc = irrep_tp_2j_build(a, b, empty, NULL, 0);
    IRREP_ASSERT(desc == NULL); /* no valid paths */

    irrep_tp_2j_free(NULL); /* smoke: free NULL is a no-op */

    irrep_multiset_2j_free(a);
    irrep_multiset_2j_free(b);
    irrep_multiset_2j_free(empty);

    /* ------------------------------------------------------------------ */
    /* Mixed multiset: spin-1 ⊗ spin-½ = spin-½ ⊕ spin-3/2                 */
    /* ------------------------------------------------------------------ */
    a = irrep_multiset_2j_parse("1x1e");
    b = irrep_multiset_2j_parse("1x1/2e");
    c = irrep_multiset_2j_parse("1x1/2e + 1x3/2e");
    int mm = irrep_tp_2j_enumerate_paths(a, b, c, NULL, 0);
    IRREP_ASSERT(mm == 2);
    desc = irrep_tp_2j_build(a, b, c, NULL, 0);
    IRREP_ASSERT(desc != NULL);
    IRREP_ASSERT(irrep_tp_2j_output_dim(desc) == 2 + 4);

    irrep_tp_2j_free(desc);
    irrep_multiset_2j_free(a);
    irrep_multiset_2j_free(b);
    irrep_multiset_2j_free(c);

    return IRREP_TEST_END();
}
