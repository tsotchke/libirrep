/* SPDX-License-Identifier: MIT */
/* Worked example from `docs/tutorials/05_tensor_products.md` §7 — computing
 * the five cartesian basis terms of a spin-torque-net's edge interaction
 * via `irrep_tp_apply_uvw` paths, and comparing to the hand-rolled form
 * with the documented prefactors.
 *
 * The conventions:
 *   (1o ⊗ 1o → 0e) path produces  (1/√3) · (a · b)
 *   (1o ⊗ 1o → 1e) path produces  (1/√2) · (a × b)   [in the (y,z,x) layout]
 *
 * Parity is multiplicative: odd × odd = even, so the cross product of two
 * polar vectors (both 1o) is correctly an axial pseudovector (1e). The
 * libirrep TP path filter enforces this — attempting `1o ⊗ 1o → 1o` would
 * silently have zero paths.
 *
 * This example builds the TP descriptor, applies it to two l=1 inputs,
 * rescales, and asserts bit-equivalence with the direct cross / dot. */

#include <irrep/multiset.h>
#include <irrep/tensor_product.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* libirrep's real-basis l=1 layout is (y, z, x), so the cartesian vector
 * (x, y, z) enters the TP as (y, z, x) and the l=1 output (y, z, x) must
 * be permuted back to cartesian on the way out. */
static void cart_to_real_l1_(const double cart[3], double real[3]) {
    real[0] = cart[1];  /* y */
    real[1] = cart[2];  /* z */
    real[2] = cart[0];  /* x */
}
static void real_l1_to_cart_(const double real[3], double cart[3]) {
    cart[0] = real[2];  /* x */
    cart[1] = real[0];  /* y */
    cart[2] = real[1];  /* z */
}

static void cross3_(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

static double dot3_(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

int main(void) {
    /* Two fixed l=1 inputs in cartesian. */
    const double a_cart[3] = {  0.7, -0.3,  0.5 };
    const double b_cart[3] = { -0.2,  0.4,  0.8 };

    /* -------- Build TP descriptor: 1o ⊗ 1o → 0e ⊕ 1o -------- */
    irrep_multiset_t *ms_1o  = irrep_multiset_parse("1x1o");
    irrep_multiset_t *ms_out = irrep_multiset_parse("1x0e + 1x1e");  /* parity: odd·odd = even */
    if (!ms_1o || !ms_out) return 1;

    /* Auto-enumerate paths (two: 1⊗1→0 and 1⊗1→1). */
    int n_paths = irrep_tp_enumerate_paths(ms_1o, ms_1o, ms_out, NULL, 0);
    int *paths  = malloc((size_t)n_paths * 3 * sizeof(int));
    irrep_tp_enumerate_paths(ms_1o, ms_1o, ms_out, paths, n_paths);

    tp_descriptor_t *tp = irrep_tp_build_uvw(ms_1o, ms_1o, ms_out, paths, n_paths);
    if (!tp) { fprintf(stderr, "tp_build_uvw failed\n"); return 1; }

    int nw = irrep_tp_num_weights_uvw(tp);
    double *w = calloc((size_t)nw, sizeof(double));
    for (int i = 0; i < nw; ++i) w[i] = 1.0;         /* identity weights */

    /* Convert cartesian → real-SH basis. */
    double a_real[3], b_real[3];
    cart_to_real_l1_(a_cart, a_real);
    cart_to_real_l1_(b_cart, b_real);

    /* TP output: 1 scalar (l=0) + 3 values (l=1). */
    double c[4] = { 0 };
    irrep_tp_apply_uvw(tp, w, a_real, b_real, c);

    /* Scale to recover the bare cartesian operations.
     *
     * Libirrep's real-basis `(1,1,1)` UVW path, combined with the i-phase
     * convention and the (y, z, x) real-SH layout, comes out with a
     * √2-scaled cross product of flipped overall sign vs the naïve
     * cartesian `a × b`. So the recovery factor is `−√2`, applied to each
     * of the three real-basis components. The remaining (y, z, x) →
     * (x, y, z) permutation then puts the result back into cartesian. */
    const double sqrt3 = sqrt(3.0);
    const double sqrt2 = sqrt(2.0);
    double tp_dot    = sqrt3  * c[0];         /* should equal a · b  */
    double tp_cross_real[3] = { sqrt2 * c[1], sqrt2 * c[2], sqrt2 * c[3] };
    double tp_cross[3];
    real_l1_to_cart_(tp_cross_real, tp_cross); /* should equal a × b (up to sign — see below) */

    /* Reference. */
    double ref_dot = dot3_(a_cart, b_cart);
    double ref_cross[3];
    cross3_(a_cart, b_cart, ref_cross);

    /* -------- Report -------- */
    printf("a        = (% .4f, % .4f, % .4f)\n", a_cart[0], a_cart[1], a_cart[2]);
    printf("b        = (% .4f, % .4f, % .4f)\n\n", b_cart[0], b_cart[1], b_cart[2]);

    printf("a · b (cartesian)            = % .10f\n", ref_dot);
    printf("√3 · tp[0e path]             = % .10f   (residual %.1e)\n",
           tp_dot, fabs(ref_dot - tp_dot));

    printf("\na × b (cartesian)            = (% .6f, % .6f, % .6f)\n",
           ref_cross[0], ref_cross[1], ref_cross[2]);
    printf("√2 · tp[1e path] (cartesian) = (% .6f, % .6f, % .6f)\n",
           tp_cross[0], tp_cross[1], tp_cross[2]);
    double max_cross_err = 0.0;
    for (int i = 0; i < 3; ++i) {
        double e = fabs(ref_cross[i] - tp_cross[i]);
        if (e > max_cross_err) max_cross_err = e;
    }
    printf("max residual                  = %.1e\n", max_cross_err);

    int ok = (fabs(ref_dot - tp_dot) < 1e-12) && (max_cross_err < 1e-12);
    printf("\n%s\n", ok ? "OK  — TP paths reproduce cartesian basis to 1e-12."
                        : "FAIL — residuals exceed threshold.");

    free(w); free(paths);
    irrep_tp_free(tp);
    irrep_multiset_free(ms_1o);
    irrep_multiset_free(ms_out);
    return ok ? 0 : 1;
}
