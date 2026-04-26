/* SPDX-License-Identifier: MIT */
/* Verify Moriya's five DMI symmetry rules (Moriya 1960 §III) one by one.
 * Each test constructs a minimal symmetry that exercises a single rule
 * and checks that the analyzer returns the textbook D-subspace. */

#include <irrep/dmi.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int total = 0, failed = 0;

#define ASSERT(cond, msg)                                                                          \
    do {                                                                                           \
        ++total;                                                                                   \
        if (!(cond)) {                                                                             \
            fprintf(stderr, "  FAIL  %s:%d  %s\n", __FILE__, __LINE__, msg);                       \
            ++failed;                                                                              \
        }                                                                                          \
    } while (0)

#define ASSERT_NEAR(a, b, tol, msg)                                                                \
    do {                                                                                           \
        ++total;                                                                                   \
        if (fabs((a) - (b)) > (tol)) {                                                             \
            fprintf(stderr, "  FAIL  %s:%d  %s  got %.6g vs %.6g\n", __FILE__, __LINE__, msg, (a), \
                    (b));                                                                          \
            ++failed;                                                                              \
        }                                                                                          \
    } while (0)

static void identity_op(irrep_dmi_sym_op_t *op) {
    memset(op, 0, sizeof *op);
    op->R_proper[0] = 1;
    op->R_proper[4] = 1;
    op->R_proper[8] = 1;
    op->det = +1;
}

/* Rotation about cartesian axis n (unit vector) by angle θ — proper rotation matrix. */
static void axis_rot(irrep_dmi_sym_op_t *op, double nx, double ny, double nz, double theta,
                     int det) {
    memset(op, 0, sizeof *op); /* zero antiunitary + any future fields */
    double n_norm = sqrt(nx * nx + ny * ny + nz * nz);
    nx /= n_norm;
    ny /= n_norm;
    nz /= n_norm;
    double c = cos(theta), s = sin(theta), C = 1 - c;
    double *R = op->R_proper;
    R[0] = c + nx * nx * C;
    R[1] = nx * ny * C - nz * s;
    R[2] = nx * nz * C + ny * s;
    R[3] = ny * nx * C + nz * s;
    R[4] = c + ny * ny * C;
    R[5] = ny * nz * C - nx * s;
    R[6] = nz * nx * C - ny * s;
    R[7] = nz * ny * C + nx * s;
    R[8] = c + nz * nz * C;
    op->det = det;
}

/* Inversion: R_proper = I, det = -1. */
static void inversion(irrep_dmi_sym_op_t *op) {
    identity_op(op);
    op->det = -1;
}

/* Mirror plane perpendicular to axis n. In (R_proper, det) form this is
 * R_proper = R(n, π), det = -1. */
static void mirror_perp_to(irrep_dmi_sym_op_t *op, double nx, double ny, double nz) {
    axis_rot(op, nx, ny, nz, M_PI, -1);
}

/* C₂ rotation about axis n. R_proper = R(n, π), det = +1. */
static void c2_axis(irrep_dmi_sym_op_t *op, double nx, double ny, double nz) {
    axis_rot(op, nx, ny, nz, M_PI, +1);
}

/* C_n rotation about z by 2π/n, det = +1. */
static void cn_z(irrep_dmi_sym_op_t *op, int n) {
    axis_rot(op, 0, 0, 1, 2.0 * M_PI / n, +1);
}

int main(void) {
    fprintf(stderr, "test_dmi: Moriya's five rules\n");

    /* Place a bond from (-1, 0, 0) to (+1, 0, 0). The bond direction is x;
     * midpoint at the origin. We add the identity to every test so the
     * projector is non-trivial. */
    double r_a[3] = {-1.0, 0.0, 0.0};
    double r_b[3] = {+1.0, 0.0, 0.0};
    double basis[9];

    /* Rule (a): inversion at midpoint → D = 0. */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        inversion(&ops[1]);
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 2, 1e-9, basis);
        ASSERT(n == 0, "(a) inversion at midpoint: D = 0");
    }

    /* Rule (b): mirror plane perpendicular to bond → D ⊥ bond.
     *           Mirror perp to x-axis swaps r_a ↔ r_b (reversing). */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        mirror_perp_to(&ops[1], 1, 0, 0);
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 2, 1e-9, basis);
        ASSERT(n == 2, "(b) mirror perp to bond: D 2-dim (perp plane)");
        /* Both basis vectors should have x-component = 0. */
        ASSERT_NEAR(basis[0], 0.0, 1e-9, "(b) basis[0].x = 0");
        ASSERT_NEAR(basis[3], 0.0, 1e-9, "(b) basis[1].x = 0");
    }

    /* Rule (c): mirror plane CONTAINING the bond → D ⊥ that plane.
     *           A plane through the x-axis; choose the xz-plane (normal y).
     *           Mirror perp to y leaves x and z fixed → preserves bond. */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        mirror_perp_to(&ops[1], 0, 1, 0); /* mirror plane = xz-plane */
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 2, 1e-9, basis);
        ASSERT(n == 1, "(c) mirror containing bond: D 1-dim (perp to plane)");
        /* Basis vector should be along y. */
        ASSERT_NEAR(fabs(basis[1]), 1.0, 1e-9, "(c) basis along y");
        ASSERT_NEAR(basis[0], 0.0, 1e-9, "(c) basis.x = 0");
        ASSERT_NEAR(basis[2], 0.0, 1e-9, "(c) basis.z = 0");
    }

    /* Rule (d): C₂ axis perpendicular to bond, through midpoint → D ⊥ axis.
     *           C₂ along y reverses the x-component → reverses bond. */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        c2_axis(&ops[1], 0, 1, 0);
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 2, 1e-9, basis);
        ASSERT(n == 2, "(d) C₂ perp to bond: D 2-dim (perp to axis)");
        ASSERT_NEAR(basis[1], 0.0, 1e-9, "(d) basis[0].y = 0");
        ASSERT_NEAR(basis[4], 0.0, 1e-9, "(d) basis[1].y = 0");
    }

    /* Rule (e): C₂ axis along bond → D ∥ bond.
     *           C₂ along x preserves x-component, flips y, z → preserves
     *           both endpoints. */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        c2_axis(&ops[1], 1, 0, 0);
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 2, 1e-9, basis);
        ASSERT(n == 1, "(e) C₂ along bond: D 1-dim (along axis)");
        ASSERT_NEAR(fabs(basis[0]), 1.0, 1e-9, "(e) basis along x");
        ASSERT_NEAR(basis[1], 0.0, 1e-9, "(e) basis.y = 0");
        ASSERT_NEAR(basis[2], 0.0, 1e-9, "(e) basis.z = 0");
    }

    /* Rule (e) extended: C_n along bond for n ≥ 2 → D ∥ bond.
     *           Use a 4-fold rotation axis along the bond. (Simplification:
     *           use C₄ along x — operationally same as a 90° rotation.) */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        axis_rot(&ops[1], 1, 0, 0, M_PI / 2.0, +1); /* C₄ along x */
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 2, 1e-9, basis);
        ASSERT(n == 1, "(e+) C₄ along bond: D 1-dim (along axis)");
        ASSERT_NEAR(fabs(basis[0]), 1.0, 1e-9, "(e+) basis along x");
    }

    /* Trivial case: only identity → D unconstrained, dim = 3. */
    {
        irrep_dmi_sym_op_t ops[1];
        identity_op(&ops[0]);
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 1, 1e-9, basis);
        ASSERT(n == 3, "trivial: D 3-dim (unconstrained)");
    }

    /* Combined: C₂ along bond + mirror containing bond.
     *           C₂(x) preserves bond and forces D ∥ x.
     *           Mirror perp to y preserves bond and forces D ⊥ xz-plane (∥ y).
     *           Intersection: D ∥ x AND D ∥ y — IMPOSSIBLE except D = 0.
     *           This is the textbook situation in centrosymmetric crystals
     *           where multiple bond stabilisers leave only D = 0. */
    {
        irrep_dmi_sym_op_t ops[3];
        identity_op(&ops[0]);
        c2_axis(&ops[1], 1, 0, 0);
        mirror_perp_to(&ops[2], 0, 1, 0);
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 3, 1e-9, basis);
        ASSERT(n == 0, "(combined C2 + mirror): conflicting constraints → D = 0");
    }

    /* Two-mirror combination: σ ⊥ bond  +  σ containing bond.
     *           σ ⊥ bond (perp to x): forces D ⊥ bond (D in yz-plane).
     *           σ in xz (perp to y): forces D ⊥ xz (D ∥ y).
     *           Intersection: D ∥ y. */
    {
        irrep_dmi_sym_op_t ops[3];
        identity_op(&ops[0]);
        mirror_perp_to(&ops[1], 1, 0, 0);
        mirror_perp_to(&ops[2], 0, 1, 0);
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 3, 1e-9, basis);
        ASSERT(n == 1, "(σ⊥bond + σ‖bond): D 1-dim");
        ASSERT_NEAR(fabs(basis[1]), 1.0, 1e-9, "two-mirror: basis along y");
    }

    /* ---- Symmetric exchange tensor analyzer ---- */

    /* Trivial-symmetry case: J^s unconstrained → 6-dim. */
    {
        irrep_dmi_sym_op_t ops[1];
        identity_op(&ops[0]);
        double sym_basis[54] = {0};
        int    n = irrep_exchange_symmetric_basis(r_a, r_b, ops, 1, 1e-9, sym_basis);
        ASSERT(n == 6, "trivial-sym J^s 6-dim (unconstrained)");
    }

    /* C₂ along bond: forces J^s to commute with R(x, π) = diag(1, -1, -1).
     * The +1 eigenspace of M_g action on symmetric J^s is 4-dim:
     *   diag(1,0,0), diag(0,1,0), diag(0,0,1), (yz+zy)/√2.
     * The xy and xz off-diagonals are forbidden. */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        c2_axis(&ops[1], 1, 0, 0);
        double sym_basis[54] = {0};
        int    n = irrep_exchange_symmetric_basis(r_a, r_b, ops, 2, 1e-9, sym_basis);
        ASSERT(n == 4, "C2 along bond: J^s 4-dim (xy and xz off-diag forbidden)");
    }

    /* Mirror perp to bond AND containing bond: forces J^s diagonal in (x, y, z).
     * Three mirrors (perp to x, perp to y, perp to z) → only the three diagonal
     * components survive → J^s 3-dim. */
    {
        irrep_dmi_sym_op_t ops[4];
        identity_op(&ops[0]);
        mirror_perp_to(&ops[1], 1, 0, 0);
        mirror_perp_to(&ops[2], 0, 1, 0);
        mirror_perp_to(&ops[3], 0, 0, 1);
        double sym_basis[54] = {0};
        int    n = irrep_exchange_symmetric_basis(r_a, r_b, ops, 4, 1e-9, sym_basis);
        ASSERT(n == 3, "three-mirror: J^s 3-dim (diagonal only)");
        /* All three basis matrices should be diagonal: off-diag entries 0. */
        for (int b_i = 0; b_i < 3; ++b_i) {
            const double *J = sym_basis + b_i * 9;
            ASSERT_NEAR(J[1], 0.0, 1e-9, "diagonal-only basis: J[xy] = 0");
            ASSERT_NEAR(J[2], 0.0, 1e-9, "diagonal-only basis: J[xz] = 0");
            ASSERT_NEAR(J[5], 0.0, 1e-9, "diagonal-only basis: J[yz] = 0");
        }
    }

    /* Tetragonal site sym (C₄ along bond + perpendicular C₂s + mirrors):
     * forces J^s to diagonal + J_yy = J_zz, leaving 2 free components
     * (the isotropic trace and the uniaxial xx vs yy=zz anisotropy). */
    {
        irrep_dmi_sym_op_t ops[7];
        identity_op(&ops[0]);
        axis_rot(&ops[1], 1, 0, 0, M_PI / 2.0, +1);  /* C₄ along bond */
        axis_rot(&ops[2], 1, 0, 0, -M_PI / 2.0, +1); /* C₄³ along bond */
        c2_axis(&ops[3], 1, 0, 0);                   /* C₂ along bond */
        c2_axis(&ops[4], 0, 1, 0);                   /* C₂ perp (reversing) */
        c2_axis(&ops[5], 0, 0, 1);                   /* C₂ perp (reversing) */
        mirror_perp_to(&ops[6], 0, 1, 0);            /* mirror containing bond */
        double sym_basis[54] = {0};
        int    n = irrep_exchange_symmetric_basis(r_a, r_b, ops, 7, 1e-9, sym_basis);
        ASSERT(n == 2, "tetragonal-bond J^s 2-dim (isotropic + uniaxial)");
        /* Both basis matrices: diagonal with J_yy = J_zz. */
        for (int b_i = 0; b_i < n; ++b_i) {
            const double *J = sym_basis + b_i * 9;
            ASSERT_NEAR(J[1], 0.0, 1e-9, "tetragonal: J[xy] = 0");
            ASSERT_NEAR(J[2], 0.0, 1e-9, "tetragonal: J[xz] = 0");
            ASSERT_NEAR(J[5], 0.0, 1e-9, "tetragonal: J[yz] = 0");
            ASSERT_NEAR(J[4], J[8], 1e-9, "tetragonal: J_yy = J_zz");
        }
    }

    /* Body-diagonal bond + 3-fold cubic axis: a bond from (-d,-d,-d)/√3 to
     * (+d,+d,+d)/√3 along (1,1,1). The C₃ rotation about (1,1,1) cycles
     * (x,y,z) → (y,z,x), which preserves the body-diagonal bond. Combined
     * with C₂ ALONG the bond (= R((1,1,1)/√3, π)), this enforces
     * J_xx = J_yy = J_zz and J_xy = J_xz = J_yz. The allowed J^s subspace
     * is 2-dim: an isotropic component (Heisenberg) plus a single off-
     * diagonal "all-equal" component (a particular Kitaev-Γ scalar
     * identifying axes). */
    {
        double r_a3[3] = {-1, -1, -1};
        double r_b3[3] = {+1, +1, +1};
        irrep_dmi_sym_op_t ops[3];
        identity_op(&ops[0]);
        axis_rot(&ops[1], 1, 1, 1, 2.0 * M_PI / 3.0, +1);  /* C₃ along (1,1,1) */
        axis_rot(&ops[2], 1, 1, 1, -2.0 * M_PI / 3.0, +1); /* C₃² along (1,1,1) */
        double sym_basis[54] = {0};
        int    n = irrep_exchange_symmetric_basis(r_a3, r_b3, ops, 3, 1e-9, sym_basis);
        ASSERT(n == 2, "C₃ body-diagonal: J^s 2-dim (iso + Γ)");
        /* Each basis matrix should have J_xx = J_yy = J_zz and J_xy = J_xz = J_yz. */
        for (int b_i = 0; b_i < n; ++b_i) {
            const double *J = sym_basis + b_i * 9;
            ASSERT_NEAR(J[0], J[4], 1e-9, "C₃ body-diag: J_xx = J_yy");
            ASSERT_NEAR(J[4], J[8], 1e-9, "C₃ body-diag: J_yy = J_zz");
            ASSERT_NEAR(J[1], J[2], 1e-9, "C₃ body-diag: J_xy = J_xz");
            ASSERT_NEAR(J[2], J[5], 1e-9, "C₃ body-diag: J_xz = J_yz");
        }
    }

    /* Note: a "fully cubic" bond stabiliser (3 perpendicular C₂ + cubic-
     * permutation symmetry) reduces J^s to 1-dim (Heisenberg only). This
     * happens for bonds at the centre of a body-diagonal axis. We don't
     * construct that here. */

    /* ---- Antiunitary (time-reversal-augmented) operations ---- */

    /* Pure time reversal alone (T·E): preserving the bond, antiunitary.
     * Sign rule: D = -R·D = -I·D = -D → D = 0.
     * This is the textbook result that pure T forbids DMI on a bond
     * with no other symmetry; physically, DMI requires SOC + broken
     * inversion AND not T-symmetric per-site magnetic structure. */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        identity_op(&ops[1]);
        ops[1].antiunitary = 1;        /* T · E */
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 2, 1e-9, basis);
        ASSERT(n == 0, "T·E pure time-reversal: D = 0 (sign opposite to E)");
    }

    /* T combined with inversion at midpoint (T·i): reversing, antiunitary.
     * Sign rule: +R = +I, so M = +I. D = +D — trivially allowed.
     * Cf. ordinary inversion at midpoint (rule a, M = -I) which forces
     * D = 0. The T flip exactly compensates the bond-reversal flip. */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        inversion(&ops[1]);
        ops[1].antiunitary = 1;        /* T · i — reversing, antiunitary */
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 2, 1e-9, basis);
        ASSERT(n == 3, "T·i: D unconstrained (T cancels inversion's bond reversal)");
    }

    /* Symmetric exchange tensor: invariant to antiunitary flag (both
     * pre- and post-T signs cancel for rank-2 axial-axial → polar). */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        identity_op(&ops[1]);
        ops[1].antiunitary = 1;
        double sym_basis[54] = {0};
        int    n = irrep_exchange_symmetric_basis(r_a, r_b, ops, 2, 1e-9, sym_basis);
        ASSERT(n == 6, "T·E on J^s: 6-dim, no constraint (T-invariant for rank-2)");
    }

    /* T·C₂ along bond: combined with C₂ along bond (preserving, unitary)
     * which forces D ∥ bond. T·C₂ adds the constraint D = -R·D where
     * R = R_x(π). Eigenvalues of -R_x(π) are -1, +1, +1 along (x, y, z).
     * So T·C₂ wants D ⊥ bond. Combined with C₂ wanting D ∥ bond:
     * intersection D = 0. */
    {
        irrep_dmi_sym_op_t ops[3];
        identity_op(&ops[0]);
        c2_axis(&ops[1], 1, 0, 0);     /* C₂ along bond, preserving, unitary */
        c2_axis(&ops[2], 1, 0, 0);
        ops[2].antiunitary = 1;        /* T·C₂ along bond, preserving, antiunitary */
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 3, 1e-9, basis);
        ASSERT(n == 0, "C₂ ∥ bond + T·C₂ ∥ bond: contradictory → D = 0");
    }

    /* Single antiunitary C₂ along bond (no plain C₂): D = -R·D, eigenvalues
     * +1, -1, -1 in (-R_x(π)) along (y, z, x) wait — let me recompute.
     * R_x(π) = diag(1, -1, -1). -R_x(π) = diag(-1, +1, +1).
     * +1 eigenspace = perp to bond. So D ⊥ bond. */
    {
        irrep_dmi_sym_op_t ops[2];
        identity_op(&ops[0]);
        c2_axis(&ops[1], 1, 0, 0);
        ops[1].antiunitary = 1;
        int n = irrep_dmi_allowed_basis(r_a, r_b, ops, 2, 1e-9, basis);
        ASSERT(n == 2, "T·C₂ ∥ bond alone: D ⊥ bond (sign-flipped vs. plain C₂)");
        ASSERT_NEAR(basis[0], 0.0, 1e-9, "T·C₂: basis[0].x = 0");
        ASSERT_NEAR(basis[3], 0.0, 1e-9, "T·C₂: basis[1].x = 0");
    }

    /* ---- Three-spin scalar chirality ---- */

    /* Equilateral triangle in xy-plane centred at origin. Vertices at
     * (1, 0, 0), (-½, +√3/2, 0), (-½, -√3/2, 0). */
    {
        double t_a[3] = {1.0, 0.0, 0.0};
        double t_b[3] = {-0.5, 0.5 * sqrt(3.0), 0.0};
        double t_c[3] = {-0.5, -0.5 * sqrt(3.0), 0.0};

        /* Just identity → trivially allowed (perm=identity, det=+1, sign=+1). */
        {
            irrep_dmi_sym_op_t ops[1];
            identity_op(&ops[0]);
            int rc = irrep_chirality_allowed(t_a, t_b, t_c, ops, 1, 1e-9);
            ASSERT(rc == 1, "chirality: identity-only stabiliser → allowed");
        }

        /* C₃ about z-axis (cycles a→b→c→a). Cyclic perm (+1) × det(+1) = +1.
         * χ remains allowed. */
        {
            irrep_dmi_sym_op_t ops[2];
            identity_op(&ops[0]);
            cn_z(&ops[1], 3);
            int rc = irrep_chirality_allowed(t_a, t_b, t_c, ops, 2, 1e-9);
            ASSERT(rc == 1, "chirality: C₃ along z (cyclic) → still allowed");
        }

        /* σ_h (mirror perpendicular to z, i.e. mirror in the triangle plane).
         * Fixes all three vertices individually (perm = identity, +1).
         * Improper, det = -1. Total sign = +1 · -1 = -1 → forbidden. */
        {
            irrep_dmi_sym_op_t ops[2];
            identity_op(&ops[0]);
            mirror_perp_to(&ops[1], 0, 0, 1);
            int rc = irrep_chirality_allowed(t_a, t_b, t_c, ops, 2, 1e-9);
            ASSERT(rc == 0, "chirality: mirror IN triangle plane → forbidden");
        }

        /* Mirror perpendicular to x-axis (σ_v): swaps b ↔ c, fixes a.
         * Permutation (a→a, b→c, c→b) = transposition, σ_perm = -1.
         * Improper, det = -1. Total sign = -1 · -1 = +1 → allowed. */
        {
            irrep_dmi_sym_op_t ops[2];
            identity_op(&ops[0]);
            mirror_perp_to(&ops[1], 1, 0, 0);
            int rc = irrep_chirality_allowed(t_a, t_b, t_c, ops, 2, 1e-9);
            ASSERT(rc == 1, "chirality: σ_v through vertex a (transp + det=-1) → allowed");
        }

        /* Combination C₃ + σ_v generates C_3v. The σ_v is a transposition
         * (= odd perm) with det = -1, total +1 (allowed by σ_v alone).
         * Combined with C₃ (cyclic + proper, total +1), still allowed. */
        {
            irrep_dmi_sym_op_t ops[3];
            identity_op(&ops[0]);
            cn_z(&ops[1], 3);
            mirror_perp_to(&ops[2], 1, 0, 0);
            int rc = irrep_chirality_allowed(t_a, t_b, t_c, ops, 3, 1e-9);
            ASSERT(rc == 1, "chirality: C_3v (C₃ + σ_v through vertex) → allowed");
        }

        /* Pure time reversal alone: T·E. T flips χ as a pseudoscalar.
         * Operation E preserves the triple, perm = identity (+1).
         * det = +1, antiunitary → σ_T = -1. Total = +1 · +1 · -1 = -1 → forbidden.
         * Physically: an "instantaneously time-reversed" magnetic structure
         * cannot have a non-zero scalar chirality if it must coincide with
         * the original — the symmetry forces χ_after = -χ_before but
         * χ_after = χ_before requires χ = 0. */
        {
            irrep_dmi_sym_op_t ops[2];
            identity_op(&ops[0]);
            identity_op(&ops[1]);
            ops[1].antiunitary = 1;
            int rc = irrep_chirality_allowed(t_a, t_b, t_c, ops, 2, 1e-9);
            ASSERT(rc == 0, "chirality: T·E pure time-reversal → forbidden");
        }

        /* T·σ_h: σ_h was previously a forbidder (sign -1). T flips it again,
         * giving total +1 → allowed. This is the "magnetic mirror" structure
         * that admits scalar chirality even with σ_h present. Relevant to
         * Mn_3Sn-type non-collinear AFM where the spatial mirror is broken
         * by the magnetic structure but T·σ_h survives as a magnetic-point-
         * group element. */
        {
            irrep_dmi_sym_op_t ops[2];
            identity_op(&ops[0]);
            mirror_perp_to(&ops[1], 0, 0, 1);
            ops[1].antiunitary = 1;
            int rc = irrep_chirality_allowed(t_a, t_b, t_c, ops, 2, 1e-9);
            ASSERT(rc == 1, "chirality: T·σ_h (magnetic mirror) → allowed");
        }
    }

    /* unused warning suppression */
    (void)cn_z;

    fprintf(stderr, "  %d / %d assertions passed\n", total - failed, total);
    return failed == 0 ? 0 : 1;
}
