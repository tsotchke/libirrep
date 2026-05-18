/* SPDX-License-Identifier: MIT */
/* Tests for `<irrep/se3.h>`: SE(3) group, se(3) Lie algebra, E(3). */
#include "harness.h"
#include <irrep/se3.h>

#include <math.h>
#include <stdio.h>

#define TOL_TIGHT 1e-12
#define TOL_LOOSE 1e-10

/* ---------- helpers ---------- */

static double diff_R(const double a[9], const double b[9]) {
    double d = 0;
    for (int i = 0; i < 9; ++i) d = fmax(d, fabs(a[i] - b[i]));
    return d;
}
static double diff_t(const double a[3], const double b[3]) {
    double d = 0;
    for (int i = 0; i < 3; ++i) d = fmax(d, fabs(a[i] - b[i]));
    return d;
}
static double diff_se3(const irrep_se3_t *a, const irrep_se3_t *b) {
    return fmax(diff_R(a->R, b->R), diff_t(a->t, b->t));
}
static double diff_twist(const irrep_se3_twist_t *a, const irrep_se3_twist_t *b) {
    double d = 0;
    for (int i = 0; i < 3; ++i) {
        d = fmax(d, fabs(a->omega[i] - b->omega[i]));
        d = fmax(d, fabs(a->v[i]     - b->v[i]));
    }
    return d;
}

/* ---------- group axioms ---------- */

static int test_identity_and_inverse(void) {
    IRREP_TEST_START("se3_identity_and_inverse");
    irrep_se3_t id = irrep_se3_identity();
    IRREP_ASSERT_NEAR(id.R[0], 1.0, TOL_TIGHT);
    IRREP_ASSERT_NEAR(id.R[4], 1.0, TOL_TIGHT);
    IRREP_ASSERT_NEAR(id.R[8], 1.0, TOL_TIGHT);
    IRREP_ASSERT_NEAR(id.t[0], 0.0, TOL_TIGHT);

    /* g · g^{-1} = id for a generic element. */
    double axis[3] = { 0.3, -0.5, 0.8 };
    double t[3] = { 1.5, -2.0, 0.7 };
    irrep_se3_t g = irrep_se3_from_axis_angle_t(axis, 0.6, t);
    irrep_se3_t inv = irrep_se3_inverse(&g);
    irrep_se3_t prod = irrep_se3_compose(&g, &inv);
    IRREP_ASSERT_NEAR(diff_se3(&prod, &id), 0.0, TOL_TIGHT);

    /* g^{-1} · g = id. */
    irrep_se3_t prod2 = irrep_se3_compose(&inv, &g);
    IRREP_ASSERT_NEAR(diff_se3(&prod2, &id), 0.0, TOL_TIGHT);
    return IRREP_TEST_END();
}

static int test_composition_associativity(void) {
    IRREP_TEST_START("se3_composition_associativity");
    /* Three random-ish elements: check (a ∘ b) ∘ c = a ∘ (b ∘ c). */
    double ax_a[3] = { 1, 0, 0 }, ta[3] = { 0.5, 0, 0 };
    double ax_b[3] = { 0, 1, 0 }, tb[3] = { 0, 0.7, 0 };
    double ax_c[3] = { 0, 0, 1 }, tc[3] = { 0, 0, 1.0 };
    irrep_se3_t a = irrep_se3_from_axis_angle_t(ax_a, 0.3, ta);
    irrep_se3_t b = irrep_se3_from_axis_angle_t(ax_b, 0.5, tb);
    irrep_se3_t c = irrep_se3_from_axis_angle_t(ax_c, 0.7, tc);

    irrep_se3_t ab    = irrep_se3_compose(&a, &b);
    irrep_se3_t ab_c  = irrep_se3_compose(&ab, &c);
    irrep_se3_t bc    = irrep_se3_compose(&b, &c);
    irrep_se3_t a_bc  = irrep_se3_compose(&a, &bc);

    IRREP_ASSERT_NEAR(diff_se3(&ab_c, &a_bc), 0.0, TOL_TIGHT);
    return IRREP_TEST_END();
}

static int test_action_consistency(void) {
    IRREP_TEST_START("se3_action_consistency");
    /* (a ∘ b) applied to p = a applied to (b applied to p). */
    double ax_a[3] = { 0.6, -0.4, 0.7 }, ta[3] = { 1, 2, 3 };
    double ax_b[3] = { -0.2, 0.9, 0.1 }, tb[3] = { -1, 0.5, 2 };
    irrep_se3_t a = irrep_se3_from_axis_angle_t(ax_a, 0.4, ta);
    irrep_se3_t b = irrep_se3_from_axis_angle_t(ax_b, 0.9, tb);

    double p[3] = { 0.7, -1.5, 2.2 };
    double bp[3]; irrep_se3_apply(&b, p, bp);
    double abp_chained[3]; irrep_se3_apply(&a, bp, abp_chained);

    irrep_se3_t ab = irrep_se3_compose(&a, &b);
    double abp_composed[3]; irrep_se3_apply(&ab, p, abp_composed);

    IRREP_ASSERT_NEAR(diff_t(abp_chained, abp_composed), 0.0, TOL_TIGHT);
    return IRREP_TEST_END();
}

/* ---------- Lie algebra ---------- */

static int test_exp_log_roundtrip(void) {
    IRREP_TEST_START("se3_exp_log_roundtrip");
    /* log(exp(ξ)) = ξ for generic twists. */
    irrep_se3_twist_t xi = {
        .omega = { 0.4, -0.2, 0.5 },
        .v     = { 1.2,  0.3, -0.6 },
    };
    irrep_se3_t g = irrep_se3_exp(&xi);
    irrep_se3_twist_t xi2 = irrep_se3_log(&g);
    IRREP_ASSERT_NEAR(diff_twist(&xi, &xi2), 0.0, TOL_LOOSE);
    return IRREP_TEST_END();
}

static int test_exp_zero_is_identity(void) {
    IRREP_TEST_START("se3_exp_zero_is_identity");
    irrep_se3_twist_t xi = { .omega = {0,0,0}, .v = {0,0,0} };
    irrep_se3_t g = irrep_se3_exp(&xi);
    irrep_se3_t id = irrep_se3_identity();
    IRREP_ASSERT_NEAR(diff_se3(&g, &id), 0.0, TOL_TIGHT);
    return IRREP_TEST_END();
}

static int test_exp_pure_rotation(void) {
    IRREP_TEST_START("se3_exp_pure_rotation");
    /* exp((ω, 0)) = (R(ω), 0): pure rotation, zero translation. */
    irrep_se3_twist_t xi = { .omega = { 0, 0, M_PI / 2 }, .v = {0,0,0} };
    irrep_se3_t g = irrep_se3_exp(&xi);
    /* R is a +π/2 rotation about z. */
    IRREP_ASSERT_NEAR(g.R[0],  0.0, TOL_TIGHT);
    IRREP_ASSERT_NEAR(g.R[1], -1.0, TOL_TIGHT);
    IRREP_ASSERT_NEAR(g.R[3],  1.0, TOL_TIGHT);
    IRREP_ASSERT_NEAR(g.R[4],  0.0, TOL_TIGHT);
    IRREP_ASSERT_NEAR(g.R[8],  1.0, TOL_TIGHT);
    IRREP_ASSERT_NEAR(g.t[0],  0.0, TOL_TIGHT);
    IRREP_ASSERT_NEAR(g.t[1],  0.0, TOL_TIGHT);
    IRREP_ASSERT_NEAR(g.t[2],  0.0, TOL_TIGHT);
    return IRREP_TEST_END();
}

static int test_exp_pure_translation(void) {
    IRREP_TEST_START("se3_exp_pure_translation");
    /* exp((0, v)) = (I, v): V(0) = I, so translation is straight. */
    irrep_se3_twist_t xi = { .omega = {0,0,0}, .v = { 1.5, -0.7, 0.3 } };
    irrep_se3_t g = irrep_se3_exp(&xi);
    irrep_se3_t id = irrep_se3_identity();
    /* R = I. */
    IRREP_ASSERT_NEAR(diff_R(g.R, id.R), 0.0, TOL_TIGHT);
    /* t = v. */
    IRREP_ASSERT_NEAR(g.t[0], 1.5, TOL_TIGHT);
    IRREP_ASSERT_NEAR(g.t[1], -0.7, TOL_TIGHT);
    IRREP_ASSERT_NEAR(g.t[2], 0.3, TOL_TIGHT);
    return IRREP_TEST_END();
}

static int test_exp_small_omega_stability(void) {
    IRREP_TEST_START("se3_exp_small_omega_stability");
    /* Near-zero |ω|: the closed-form trig has 0/0; verify the Taylor
     * fallback delivers machine precision at θ ~ 1e-8. */
    irrep_se3_twist_t xi = { .omega = { 1e-10, 1e-10, 1e-10 }, .v = { 1, 2, 3 } };
    irrep_se3_t g = irrep_se3_exp(&xi);
    irrep_se3_twist_t back = irrep_se3_log(&g);
    IRREP_ASSERT_NEAR(diff_twist(&xi, &back), 0.0, TOL_LOOSE);
    return IRREP_TEST_END();
}

/* ---------- adjoint ---------- */

static int test_adjoint_satisfies_conjugation(void) {
    IRREP_TEST_START("se3_adjoint_conjugation");
    /* For ξ ∈ se(3) and g ∈ SE(3):  g · exp(ξ) · g^{-1} = exp(Ad_g · ξ). */
    double ax[3] = { 0.4, -0.3, 0.8 };
    double t[3] = { 1.0, 2.0, -0.5 };
    irrep_se3_t g = irrep_se3_from_axis_angle_t(ax, 0.6, t);

    irrep_se3_twist_t xi = { .omega = { 0.2, 0.1, -0.3 }, .v = { 0.5, -0.4, 0.1 } };

    /* LHS: g · exp(ξ) · g^{-1}. */
    irrep_se3_t exp_xi = irrep_se3_exp(&xi);
    irrep_se3_t g_inv = irrep_se3_inverse(&g);
    irrep_se3_t tmp = irrep_se3_compose(&g, &exp_xi);
    irrep_se3_t lhs = irrep_se3_compose(&tmp, &g_inv);

    /* RHS: exp(Ad_g · ξ). */
    double Ad[36];
    irrep_se3_adjoint(&g, Ad);
    double xi_vec[6] = { xi.omega[0], xi.omega[1], xi.omega[2],
                          xi.v[0],     xi.v[1],     xi.v[2] };
    double Ad_xi[6] = {0};
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            Ad_xi[i] += Ad[i*6 + j] * xi_vec[j];
    irrep_se3_twist_t Ad_xi_twist = {
        .omega = { Ad_xi[0], Ad_xi[1], Ad_xi[2] },
        .v     = { Ad_xi[3], Ad_xi[4], Ad_xi[5] },
    };
    irrep_se3_t rhs = irrep_se3_exp(&Ad_xi_twist);

    IRREP_ASSERT_NEAR(diff_se3(&lhs, &rhs), 0.0, TOL_LOOSE);
    return IRREP_TEST_END();
}

/* ---------- hat / vee ---------- */

static int test_hat_vee_roundtrip(void) {
    IRREP_TEST_START("se3_hat_vee_roundtrip");
    irrep_se3_twist_t xi = { .omega = { 0.5, -0.7, 1.1 }, .v = { 0.3, 0.2, -0.9 } };
    double mat[16];
    irrep_se3_hat(&xi, mat);
    irrep_se3_twist_t xi2;
    irrep_se3_vee(mat, &xi2);
    IRREP_ASSERT_NEAR(diff_twist(&xi, &xi2), 0.0, TOL_TIGHT);
    return IRREP_TEST_END();
}

/* ---------- renormalise ---------- */

static int test_renormalise_fixes_drift(void) {
    IRREP_TEST_START("se3_renormalise_fixes_drift");
    /* Compose 1000 small rotations; without renormalise, R drifts off SO(3).
     * After renormalise the result should still be orthonormal. */
    irrep_se3_t g = irrep_se3_identity();
    double ax[3] = { 1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3) };
    double t0[3] = {0,0,0};
    irrep_se3_t step = irrep_se3_from_axis_angle_t(ax, 0.001, t0);
    for (int i = 0; i < 1000; ++i) g = irrep_se3_compose(&g, &step);
    irrep_se3_renormalise(&g);
    /* R · R^T = I after renormalise. */
    double RRT[9] = {0};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                RRT[i*3+j] += g.R[i*3+k] * g.R[j*3+k];
    double id[9] = {1,0,0, 0,1,0, 0,0,1};
    IRREP_ASSERT_NEAR(diff_R(RRT, id), 0.0, TOL_LOOSE);
    return IRREP_TEST_END();
}

/* ---------- E(3) ---------- */

static int test_e3_parity_compose(void) {
    IRREP_TEST_START("e3_parity_compose");
    irrep_e3_t id = irrep_e3_identity();
    irrep_e3_t inv = irrep_e3_inversion();
    IRREP_ASSERT(id.parity == +1);
    IRREP_ASSERT(inv.parity == -1);

    irrep_e3_t inv2 = irrep_e3_compose(&inv, &inv);
    IRREP_ASSERT(inv2.parity == +1);
    /* Inversion squared = identity. */
    IRREP_ASSERT_NEAR(diff_R(inv2.se3.R, id.se3.R), 0.0, TOL_TIGHT);
    return IRREP_TEST_END();
}

static int test_e3_inversion_acts_as_negation(void) {
    IRREP_TEST_START("e3_inversion_acts_as_negation");
    irrep_e3_t inv = irrep_e3_inversion();
    double p[3] = { 1.5, -0.7, 2.0 };
    double out[3];
    irrep_e3_apply(&inv, p, out);
    IRREP_ASSERT_NEAR(out[0], -1.5, TOL_TIGHT);
    IRREP_ASSERT_NEAR(out[1],  0.7, TOL_TIGHT);
    IRREP_ASSERT_NEAR(out[2], -2.0, TOL_TIGHT);
    return IRREP_TEST_END();
}

int main(void) {
    int rc = 0;
    rc |= test_identity_and_inverse();
    rc |= test_composition_associativity();
    rc |= test_action_consistency();
    rc |= test_exp_log_roundtrip();
    rc |= test_exp_zero_is_identity();
    rc |= test_exp_pure_rotation();
    rc |= test_exp_pure_translation();
    rc |= test_exp_small_omega_stability();
    rc |= test_adjoint_satisfies_conjugation();
    rc |= test_hat_vee_roundtrip();
    rc |= test_renormalise_fixes_drift();
    rc |= test_e3_parity_compose();
    rc |= test_e3_inversion_acts_as_negation();
    return rc;
}
