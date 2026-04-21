/* SPDX-License-Identifier: MIT */
/* Tests the `irrep_nequip_layer_from_spec` parser per coord doc PR#5 draft.
 *
 * Structure:
 *   - 5 valid cases: parse → apply → compare bit-exact to `_build` output.
 *   - 5 malformed cases: parse returns NULL, irrep_last_error() is non-empty.
 *
 * The fixture is a small deterministic graph (4 nodes, 8 edges) with
 * splitmix-seeded h_in features. Shared across every valid-case check so
 * both layer instances see identical inputs. */

#include "harness.h"
#include <irrep/multiset.h>
#include <irrep/nequip.h>
#include <irrep/types.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- splitmix64 fixture generator ---------------------------------------- */
static uint64_t sm64_(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static double sm64_unit_(uint64_t *s) {
    return (double)(sm64_(s) >> 11) / (double)(1ULL << 53);
}

#define FIX_N_NODES 4
#define FIX_N_EDGES 8
static int    fix_src[FIX_N_EDGES] = {0, 1, 1, 2, 2, 3, 3, 0};
static int    fix_dst[FIX_N_EDGES] = {1, 0, 2, 1, 3, 2, 0, 3};
static double fix_edge_vec[FIX_N_EDGES * 3];

static void   fill_fixture_(void) {
    uint64_t s = 0xA5A5A5A5A5A5A5A5ULL;
    for (int i = 0; i < FIX_N_EDGES * 3; ++i) {
        fix_edge_vec[i] = 2.0 * sm64_unit_(&s) - 1.0;
    }
    /* Rescale to unit-ish magnitudes; keep r ∈ (0, 1). */
    for (int e = 0; e < FIX_N_EDGES; ++e) {
        double *v = fix_edge_vec + 3 * e;
        double  n = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (n == 0.0) {
            v[0] = 0.1;
            v[1] = 0.2;
            v[2] = 0.3;
            n = sqrt(0.14);
        }
        double target = 0.3 + 0.4 * sm64_unit_(&s); /* ∈ [0.3, 0.7] */
        double k = target / n;
        v[0] *= k;
        v[1] *= k;
        v[2] *= k;
    }
}

/* --- helper: apply both layers on the fixture and assert bit-equality ---- */
static int layers_agree_(irrep_nequip_layer_t *A, irrep_nequip_layer_t *B, int h_in_dim,
                         int h_out_dim) {
    int nwA = irrep_nequip_layer_num_weights(A);
    int nwB = irrep_nequip_layer_num_weights(B);
    if (nwA != nwB)
        return 0;

    double  *w = calloc((size_t)nwA, sizeof(double));
    double  *h_in = calloc((size_t)FIX_N_NODES * h_in_dim, sizeof(double));
    double  *h_outA = calloc((size_t)FIX_N_NODES * h_out_dim, sizeof(double));
    double  *h_outB = calloc((size_t)FIX_N_NODES * h_out_dim, sizeof(double));

    uint64_t s = 0x12345678DEADBEEFULL;
    for (int i = 0; i < nwA; ++i)
        w[i] = 2.0 * sm64_unit_(&s) - 1.0;
    for (int i = 0; i < FIX_N_NODES * h_in_dim; ++i)
        h_in[i] = 2.0 * sm64_unit_(&s) - 1.0;

    irrep_nequip_layer_apply(A, w, FIX_N_NODES, FIX_N_EDGES, fix_src, fix_dst, fix_edge_vec, h_in,
                             h_outA);
    irrep_nequip_layer_apply(B, w, FIX_N_NODES, FIX_N_EDGES, fix_src, fix_dst, fix_edge_vec, h_in,
                             h_outB);

    int ok = 1;
    for (int i = 0; i < FIX_N_NODES * h_out_dim; ++i) {
        if (h_outA[i] != h_outB[i]) {
            ok = 0;
            break;
        }
    }
    free(w);
    free(h_in);
    free(h_outA);
    free(h_outB);
    return ok;
}

/* Build the reference layer with the verbose API, for round-trip comparison. */
static irrep_nequip_layer_t *build_verbose_(const char *hi_spec, const char *ho_spec, int sh,
                                            int n_radial, double r_cut,
                                            irrep_nequip_cutoff_t cutoff_kind, int cutoff_p) {

    irrep_multiset_t     *hi = irrep_multiset_parse(hi_spec);
    irrep_multiset_t     *ho = irrep_multiset_parse(ho_spec);
    irrep_nequip_layer_t *layer =
        irrep_nequip_layer_build(hi, sh, n_radial, r_cut, cutoff_kind, cutoff_p, ho);
    irrep_multiset_free(hi);
    irrep_multiset_free(ho);
    return layer;
}

/* Valid-case test driver. */
static void check_valid_(const char *spec, const char *hi_spec, const char *ho_spec, int sh,
                         int n_radial, double r_cut, irrep_nequip_cutoff_t cutoff_kind,
                         int cutoff_p, int hi_dim, int ho_dim) {
    irrep_nequip_layer_t *A = irrep_nequip_layer_from_spec(spec);
    IRREP_ASSERT(A != NULL);
    irrep_nequip_layer_t *B =
        build_verbose_(hi_spec, ho_spec, sh, n_radial, r_cut, cutoff_kind, cutoff_p);
    IRREP_ASSERT(B != NULL);
    IRREP_ASSERT(layers_agree_(A, B, hi_dim, ho_dim));
    irrep_nequip_layer_free(A);
    irrep_nequip_layer_free(B);
}

/* Malformed-case driver — expect NULL + non-empty last-error. */
static void check_malformed_(const char *spec) {
    irrep_nequip_layer_t *L = irrep_nequip_layer_from_spec(spec);
    IRREP_ASSERT(L == NULL);
    const char *err = irrep_last_error();
    IRREP_ASSERT(err != NULL);
    IRREP_ASSERT(err[0] != '\0');
}

int main(void) {
    IRREP_TEST_START("nequip_from_spec");
    fill_fixture_();

    /* ---- 1–5 valid cases ---- */

    /* (1) Default everything. */
    check_valid_("1x0e + 1x1o -> 1x1o", "1x0e + 1x1o", "1x1o",
                 /*sh=*/2, /*radial=*/8, /*r_cut=*/1.0, IRREP_NEQUIP_CUTOFF_POLYNOMIAL, /*p=*/6,
                 /*hi_dim=*/1 + 3, /*ho_dim=*/3);

    /* (2) Override sh. */
    check_valid_("1x0e + 1x1o -> 1x1o [sh=3]", "1x0e + 1x1o", "1x1o", 3, 8, 1.0,
                 IRREP_NEQUIP_CUTOFF_POLYNOMIAL, 6, 4, 3);

    /* (3) Override radial. */
    check_valid_("1x0e + 1x1o -> 1x1o [radial=16]", "1x0e + 1x1o", "1x1o", 2, 16, 1.0,
                 IRREP_NEQUIP_CUTOFF_POLYNOMIAL, 6, 4, 3);

    /* (4) Override r_cut. */
    check_valid_("1x0e + 1x1o -> 1x1o [r_cut=1.5]", "1x0e + 1x1o", "1x1o", 2, 8, 1.5,
                 IRREP_NEQUIP_CUTOFF_POLYNOMIAL, 6, 4, 3);

    /* (5) Override everything (cosine cutoff, large shape). */
    check_valid_("4x0e + 2x1o + 1x2e -> 2x0e + 1x1o "
                 "[sh=3, radial=16, r_cut=1.5, cutoff=cosine]",
                 "4x0e + 2x1o + 1x2e", "2x0e + 1x1o", 3, 16, 1.5, IRREP_NEQUIP_CUTOFF_COSINE, 0,
                 /*hi_dim=*/4 + 2 * 3 + 5, /*ho_dim=*/2 + 3);

    /* ---- 6–10 malformed cases ---- */
    check_malformed_("1x0e + 1x1o 1x1o");                    /* missing arrow */
    check_malformed_("1x0e + 1x1o -> 1x1o +");               /* trailing operator */
    check_malformed_("1q0e + 1x1o -> 1x1o");                 /* invalid irrep */
    check_malformed_("1x0e + 1x1o -> 1x1o [cutoff=linear]"); /* unknown cutoff */
    check_malformed_("1x0e + 1x1o -> 1x1o [zoom=2]");        /* unknown option */

    /* ---- Additional parser boundary coverage (post-audit) ---- */
    check_malformed_("1x0e + 1x1o -> 1x1o [cutoff=polynomial(0)]");  /* poly order < 1 */
    check_malformed_("1x0e + 1x1o -> 1x1o [cutoff=polynomial(-3)]"); /* negative poly order */
    check_malformed_("1x0e + 1x1o -> 1x1o [sh=-1]");                 /* negative sh */
    check_malformed_("1x0e + 1x1o -> 1x1o [radial=0]");              /* radial < 1 */
    check_malformed_("1x0e + 1x1o -> 1x1o [r_cut=-1.0]");            /* negative r_cut */
    check_malformed_("1x0e + 1x1o -> 1x1o [r_cut=0.0]");             /* zero r_cut */
    check_malformed_("1x0e + 1x1o -> 1x1o [sh=999999999999]");       /* overflows INT_MAX */
    check_malformed_("   -> 1x1o");                                  /* empty hidden_in */
    check_malformed_("1x0e + 1x1o ->");                              /* empty hidden_out */

    /* ---- Empty options block is accepted (documented) ---- */
    {
        irrep_nequip_layer_t *A = irrep_nequip_layer_from_spec("1x0e + 1x1o -> 1x1o []");
        IRREP_ASSERT(A != NULL);
        irrep_nequip_layer_free(A);
    }

    /* ---- Duplicate options: last value wins (documented behaviour) ----
     * No diagnostic; caller gets the most-recent assignment for each key.
     * Lock this behaviour in so future parser changes don't silently alter it. */
    {
        irrep_nequip_layer_t *A = irrep_nequip_layer_from_spec("1x0e + 1x1o -> 1x1o [sh=2, sh=4]");
        IRREP_ASSERT(A != NULL);
        /* Can't introspect sh directly from the opaque layer, but if sh=4 took
         * effect the weight count will match the sh=4 reference. */
        irrep_nequip_layer_t *B = build_verbose_("1x0e + 1x1o", "1x1o", /*sh=*/4, /*radial=*/8,
                                                 /*r_cut=*/1.0, IRREP_NEQUIP_CUTOFF_POLYNOMIAL, 6);
        IRREP_ASSERT(B != NULL);
        IRREP_ASSERT(irrep_nequip_layer_num_weights(A) == irrep_nequip_layer_num_weights(B));
        irrep_nequip_layer_free(A);
        irrep_nequip_layer_free(B);
    }

    /* ---- Bonus: whitespace-stripped equivalent round-trip ---- */
    {
        irrep_nequip_layer_t *A =
            irrep_nequip_layer_from_spec("8x0e+4x1o+2x2e+1x3o->4x0e+2x1o+1x2e[sh=4,radial=16]");
        IRREP_ASSERT(A != NULL);
        irrep_nequip_layer_t *B = build_verbose_("8x0e + 4x1o + 2x2e + 1x3o", "4x0e + 2x1o + 1x2e",
                                                 4, 16, 1.0, IRREP_NEQUIP_CUTOFF_POLYNOMIAL, 6);
        IRREP_ASSERT(B != NULL);
        IRREP_ASSERT(layers_agree_(A, B,
                                   /*hi_dim=*/8 + 4 * 3 + 2 * 5 + 1 * 7,
                                   /*ho_dim=*/4 + 2 * 3 + 1 * 5));
        irrep_nequip_layer_free(A);
        irrep_nequip_layer_free(B);
    }

    /* ---- Bonus: NULL input doesn't crash ---- */
    {
        irrep_nequip_layer_t *L = irrep_nequip_layer_from_spec(NULL);
        IRREP_ASSERT(L == NULL);
    }

    return IRREP_TEST_END();
}
