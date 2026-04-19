/* SPDX-License-Identifier: MIT */
/* Exercises the error / boundary branches that pass-fail tests skip:
 * NULL inputs, out-of-range l / j / m, bad multiset strings, parse
 * failures, builder failures. Each branch should return the documented
 * failure signal (NULL, 0.0, or IRREP_ERR_*) without crashing. */

#include "harness.h"
#include <irrep/clebsch_gordan.h>
#include <irrep/multiset.h>
#include <irrep/nequip.h>
#include <irrep/radial.h>
#include <irrep/recoupling.h>
#include <irrep/solid_harmonics.h>
#include <irrep/spherical_harmonics.h>
#include <irrep/tensor_product.h>
#include <irrep/time_reversal.h>
#include <irrep/types.h>
#include <irrep/wigner_d.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
    IRREP_TEST_START("error_paths");

    /* ---- irrep_cg: invalid args return 0.0, not NaN ---- */
    IRREP_ASSERT(irrep_cg(-1,  0,  0, 0, 0, 0) == 0.0);  /* j1 < 0 */
    IRREP_ASSERT(irrep_cg( 0,  5,  0, 0, 0, 0) == 0.0);  /* |m1| > j1 */
    IRREP_ASSERT(irrep_cg( 1,  0,  1, 0, 5, 0) == 0.0);  /* J outside triangle */
    IRREP_ASSERT(irrep_cg( 1,  1,  1, 1, 2, 0) == 0.0);  /* m-sum mismatch */

    /* ---- irrep_wigner_d_small: invalid args ---- */
    IRREP_ASSERT(irrep_wigner_d_small(-1, 0, 0, 0.0) == 0.0);
    IRREP_ASSERT(irrep_wigner_d_small(1, 5, 0, 0.0) == 0.0);

    /* ---- irrep_wigner_6j: selection rule violations ---- */
    IRREP_ASSERT(irrep_wigner_6j(-1, 0, 0, 0, 0, 0) == 0.0);

    /* ---- spherical harmonics: bad l / m returns 0 ---- */
    IRREP_ASSERT(irrep_sph_harm_real(-1, 0, 0.1, 0.2) == 0.0);
    IRREP_ASSERT(irrep_sph_harm_real(2,  5, 0.1, 0.2) == 0.0);
    IRREP_ASSERT(irrep_legendre_assoc(-1, 0, 0.3) == 0.0);
    IRREP_ASSERT(irrep_legendre_assoc(1, 5, 0.3) == 0.0);

    /* ---- solid harmonics: l out of range is a no-op (not UB) ---- */
    {
        double r[3] = { 0.5, 0.4, 0.3 };
        double buf[64];
        for (int i = 0; i < 64; ++i) buf[i] = -12345.0;        /* sentinel */
        irrep_solid_harm_cart(-1, buf, r);                     /* must leave untouched */
        IRREP_ASSERT(buf[0] == -12345.0);
        irrep_solid_harm_cart(IRREP_SOLID_L_MAX + 1, buf, r);  /* out-of-range */
        IRREP_ASSERT(buf[0] == -12345.0);
        irrep_solid_harm_cart_grad(-1, buf, r);
        IRREP_ASSERT(buf[0] == -12345.0);
        irrep_solid_harm_cart_grad(IRREP_SOLID_L_MAX + 1, buf, r);
        IRREP_ASSERT(buf[0] == -12345.0);
        irrep_solid_harm_cart_all(-1, buf, r);
        IRREP_ASSERT(buf[0] == -12345.0);
        irrep_solid_harm_cart_all(IRREP_SOLID_L_MAX + 1, buf, r);
        IRREP_ASSERT(buf[0] == -12345.0);
    }

    /* ---- Bessel / cutoffs: zero / negative / out-of-range r ---- */
    IRREP_ASSERT(irrep_rbf_bessel(1, -1.0, 2.0) == 0.0);
    IRREP_ASSERT(irrep_rbf_bessel(1,  3.0, 2.0) == 0.0);
    IRREP_ASSERT(irrep_cutoff_cosine(-0.1, 2.0) == 0.0);
    IRREP_ASSERT(irrep_cutoff_cosine(3.0, 2.0)  == 0.0);
    IRREP_ASSERT(irrep_cutoff_cosine(1.0, 0.0)  == 0.0);  /* r_cut ≤ 0 */
    IRREP_ASSERT(irrep_cutoff_polynomial(1.0, 2.0, 0) == 0.0);   /* p < 1 */
    IRREP_ASSERT(irrep_cutoff_polynomial_d(-0.1, 2.0, 6) == 0.0);

    /* ---- multiset parser: garbage strings return NULL ----
     * Empty input is documented as a valid empty multiset (matches the e3nn
     * convention for the trivial space), so it gets its own separate check. */
    {
        irrep_multiset_t *empty = irrep_multiset_parse("");
        IRREP_ASSERT(empty != NULL);
        IRREP_ASSERT(empty->num_terms == 0);
        IRREP_ASSERT(empty->total_dim == 0);
        irrep_multiset_free(empty);

        IRREP_ASSERT(irrep_multiset_parse("foo")             == NULL);
        IRREP_ASSERT(irrep_multiset_parse("-1x0e")           == NULL);
        IRREP_ASSERT(irrep_multiset_parse("1x0x")            == NULL);  /* bad parity */
        IRREP_ASSERT(irrep_multiset_parse("1x")              == NULL);  /* missing l */
        IRREP_ASSERT(irrep_multiset_parse("1x0e + ")         == NULL);  /* dangling + */

        /* Valid round-trips stay valid. */
        irrep_multiset_t *m = irrep_multiset_parse("2x0e + 1x1o");
        IRREP_ASSERT(m != NULL);
        IRREP_ASSERT(m->total_dim == 2 * 1 + 1 * 3);
        irrep_multiset_free(m);
    }

    /* ---- tensor_product build: mismatched / empty inputs ---- */
    {
        IRREP_ASSERT(irrep_tp_build(NULL, NULL, NULL, NULL, 0) == NULL);
        IRREP_ASSERT(irrep_tp_build_uvw(NULL, NULL, NULL, NULL, 0) == NULL);

        irrep_multiset_t *a = irrep_multiset_parse("1x0e");
        irrep_multiset_t *b = irrep_multiset_parse("1x0e");
        irrep_multiset_t *c = irrep_multiset_parse("1x0e");
        /* Bad path count (negative) */
        int bogus_path[3] = { 5, 5, 5 };  /* indices out of range */
        tp_descriptor_t *d = irrep_tp_build(a, b, c, bogus_path, 1);
        IRREP_ASSERT(d == NULL);
        irrep_multiset_free(a);
        irrep_multiset_free(b);
        irrep_multiset_free(c);
    }

    /* ---- NequIP build: bad params fail cleanly ---- */
    {
        irrep_multiset_t *h = irrep_multiset_parse("1x0e + 1x1o");
        IRREP_ASSERT(h != NULL);
        /* r_cut = 0 */
        IRREP_ASSERT(irrep_nequip_layer_build(h, 2, 4, 0.0,
            IRREP_NEQUIP_CUTOFF_POLYNOMIAL, 6, h) == NULL);
        /* n_radial = 0 */
        IRREP_ASSERT(irrep_nequip_layer_build(h, 2, 0, 3.0,
            IRREP_NEQUIP_CUTOFF_POLYNOMIAL, 6, h) == NULL);
        /* l_sh_max negative */
        IRREP_ASSERT(irrep_nequip_layer_build(h, -1, 4, 3.0,
            IRREP_NEQUIP_CUTOFF_POLYNOMIAL, 6, h) == NULL);
        irrep_multiset_free(h);
    }

    /* ---- NULL output buffer must not crash ---- */
    irrep_sph_harm_cart(1, NULL, (double[3]){1, 0, 0});           /* silent return */
    irrep_solid_harm_cart(1, NULL, (double[3]){1, 0, 0});

    return IRREP_TEST_END();
}
