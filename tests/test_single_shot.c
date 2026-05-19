/* SPDX-License-Identifier: MIT */
/* Tests for the single-shot QEC framework, validated against the
 * 3D toric code (Quintavalle-Vasmer-Roffe-Campbell 2021).
 *
 * Verifies:
 *  - Meta-check matrix shape: M_X has m_X_meta rows × n_X cols, etc.
 *  - 3D toric meta-checks built explicitly:
 *      M_X = single row of all-1s (the global vertex-redundancy:
 *            product of all vertex stars = identity).
 *      M_Z = one row per cube, with 1s on the 6 face-stab indices
 *            (the cube-redundancy: product of 6 face stabs = identity).
 *  - Verify property: M_X · H_X = 0 and M_Z · H_Z = 0 over F₂.
 *  - Meta-syndrome of zero data error = zero (sanity).
 *  - Meta-syndrome decoupling: a fabricated measurement-error pattern
 *    produces a non-trivial meta-syndrome, while a fabricated data-error
 *    pattern (its true syndrome) produces a zero meta-syndrome.
 */
#include "harness.h"
#include <irrep/css_code.h>
#include <irrep/single_shot.h>
#include <irrep/toric3d.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* For the 3D toric on Lx×Ly×Lz, populate M_X and M_Z with the
 * vertex-redundancy and cube-redundancy meta-checks. */
static int build_meta_for_toric3d(irrep_single_shot_code_t *ss,
                                  int Lx, int Ly, int Lz) {
    int N = Lx * Ly * Lz;
    /* M_X: 1 row, all-1s of length n_X = N. */
    for (int i = 0; i < N; ++i) {
        irrep_parity_matrix_set(&ss->M_X, 0, i);
    }
    /* M_Z: N rows (one per cube), with 1s on the 6 face-stab indices.
     * Face-stab ordering in our toric3d build: xy-block [0, N), yz-block
     * [N, 2N), xz-block [2N, 3N). Within each block, row-major (z, y, x).
     *
     * Cube (cx, cy, cz) faces:
     *   bottom-xy = xy(cx, cy, cz)         → index linear_idx(cx, cy, cz)
     *   top-xy    = xy(cx, cy, cz+1)       → linear_idx(cx, cy, cz+1) (PBC)
     *   left-yz   = yz(cx, cy, cz)         → N + linear_idx(cx, cy, cz)
     *   right-yz  = yz(cx+1, cy, cz)       → N + linear_idx(cx+1, cy, cz)
     *   front-xz  = xz(cx, cy, cz)         → 2N + linear_idx(cx, cy, cz)
     *   back-xz   = xz(cx, cy+1, cz)       → 2N + linear_idx(cx, cy+1, cz)
     */
    for (int cz = 0; cz < Lz; ++cz) {
        for (int cy = 0; cy < Ly; ++cy) {
            for (int cx = 0; cx < Lx; ++cx) {
                int xp1 = (cx + 1) % Lx;
                int yp1 = (cy + 1) % Ly;
                int zp1 = (cz + 1) % Lz;
                int row = cz * Lx * Ly + cy * Lx + cx;

                int b_xy  = cz * Lx * Ly + cy * Lx + cx;
                int t_xy  = zp1 * Lx * Ly + cy * Lx + cx;
                int l_yz  = N + (cz * Lx * Ly + cy * Lx + cx);
                int r_yz  = N + (cz * Lx * Ly + cy * Lx + xp1);
                int f_xz  = 2 * N + (cz * Lx * Ly + cy * Lx + cx);
                int bk_xz = 2 * N + (cz * Lx * Ly + yp1 * Lx + cx);

                irrep_parity_matrix_set(&ss->M_Z, row, b_xy);
                irrep_parity_matrix_set(&ss->M_Z, row, t_xy);
                irrep_parity_matrix_set(&ss->M_Z, row, l_yz);
                irrep_parity_matrix_set(&ss->M_Z, row, r_yz);
                irrep_parity_matrix_set(&ss->M_Z, row, f_xz);
                irrep_parity_matrix_set(&ss->M_Z, row, bk_xz);
            }
        }
    }
    return 0;
}

static int test_3d_toric_meta_verify(int Lx, int Ly, int Lz) {
    irrep_toric3d_params_t p;
    irrep_toric3d_init(&p, Lx, Ly, Lz);
    irrep_css_code_t css;
    if (irrep_toric3d_build(&p, &css) != IRREP_OK) return 1;
    int N = Lx * Ly * Lz;
    irrep_single_shot_code_t ss;
    if (irrep_single_shot_code_new(&ss, &css, /*M_X meta rows=*/1, /*M_Z meta rows=*/N) != IRREP_OK) {
        irrep_css_code_free(&css);
        return 1;
    }
    int rc = 0;
    if (build_meta_for_toric3d(&ss, Lx, Ly, Lz)) rc = 1;
    if (irrep_single_shot_verify_meta(&ss) != IRREP_OK) rc = 1;
    irrep_single_shot_code_free(&ss);
    return rc;
}

static int test_meta_syndrome_zero(int Lx, int Ly, int Lz) {
    irrep_toric3d_params_t p;
    irrep_toric3d_init(&p, Lx, Ly, Lz);
    irrep_css_code_t css;
    if (irrep_toric3d_build(&p, &css) != IRREP_OK) return 1;
    int N = Lx * Ly * Lz;
    irrep_single_shot_code_t ss;
    irrep_single_shot_code_new(&ss, &css, 1, N);
    build_meta_for_toric3d(&ss, Lx, Ly, Lz);

    /* All-zero data error → all-zero syndrome → all-zero meta-syndrome. */
    int n_syndrome_words = (3 * N + 63) / 64;
    uint64_t *syndrome = (uint64_t *)calloc(n_syndrome_words, sizeof(uint64_t));
    int n_meta_words = (N + 63) / 64;
    uint64_t *meta = (uint64_t *)calloc(n_meta_words, sizeof(uint64_t));
    irrep_single_shot_meta_syndrome_Z(&ss, syndrome, meta);
    int rc = 0;
    for (int w = 0; w < n_meta_words; ++w) if (meta[w] != 0) rc = 1;
    free(syndrome);
    free(meta);
    irrep_single_shot_code_free(&ss);
    return rc;
}

/* Decoupling: a single faulty syndrome bit (a measurement error) produces
 * a non-trivial meta-syndrome, while a true syndrome from a single-qubit
 * data error produces a zero meta-syndrome. */
static int test_meta_syndrome_decoupling(int Lx, int Ly, int Lz) {
    irrep_toric3d_params_t p;
    irrep_toric3d_init(&p, Lx, Ly, Lz);
    irrep_css_code_t css;
    if (irrep_toric3d_build(&p, &css) != IRREP_OK) return 1;
    int N = Lx * Ly * Lz;
    irrep_single_shot_code_t ss;
    irrep_single_shot_code_new(&ss, &css, 1, N);
    build_meta_for_toric3d(&ss, Lx, Ly, Lz);

    int n_syndrome_words = (3 * N + 63) / 64;
    int n_meta_words = (N + 63) / 64;

    /* Case 1: fake measurement-error pattern — flip a single syndrome bit
     * (e.g., bit 0 → corresponds to the xy-face stab at vertex (0,0,0)).
     * A real data error CANNOT produce this single-bit pattern (since it
     * would have to anti-commute with exactly one face-stab; but every
     * single-edge X error anti-commutes with the 2 cube faces meeting on
     * that edge — at least 2). So this is a measurement-only syndrome. */
    uint64_t *syndrome = (uint64_t *)calloc(n_syndrome_words, sizeof(uint64_t));
    uint64_t *meta = (uint64_t *)calloc(n_meta_words, sizeof(uint64_t));
    syndrome[0] |= 1; /* flip face-stab #0 */
    irrep_single_shot_meta_syndrome_Z(&ss, syndrome, meta);
    int meta_nonzero_meas = 0;
    for (int w = 0; w < n_meta_words; ++w) if (meta[w] != 0) meta_nonzero_meas = 1;

    /* Case 2: true syndrome from a single-edge X data error.
     * Edge (0, 0, 0, x) is x-edge between vertices (0,0,0) and (1,0,0).
     * Cube faces containing this edge: 2 (its bottom-xy face of cube
     * (0,0,0) AND front-xz face of cube (0,0,0)).
     * Wait, every edge belongs to exactly 4 faces in the lattice, but
     * here we only care about 2 stabilizers in a single cube. Let me
     * compute manually: the bottom-xy face of cube (0,0,0) = face index 0;
     * the top-xy face of cube (0,0,Lz-1) = top of that wrapped cube;
     * the front-xz face of cube (0,0,0) = face 2N + 0;
     * the back-xz face of cube (0,Ly-1,0) (wrapped).
     * Let's just compute the actual H_Z column for this edge.
     */
    memset(syndrome, 0, n_syndrome_words * sizeof(uint64_t));
    memset(meta, 0, n_meta_words * sizeof(uint64_t));
    /* Build true syndrome: for each Z-stab row, set bit if it includes
     * edge 0 (= (0,0,0,x)). */
    for (int row = 0; row < ss.css.H_Z.n_rows; ++row) {
        if (irrep_parity_matrix_get(&ss.css.H_Z, row, 0)) {
            syndrome[row / 64] |= ((uint64_t)1 << (row % 64));
        }
    }
    irrep_single_shot_meta_syndrome_Z(&ss, syndrome, meta);
    int meta_nonzero_data = 0;
    for (int w = 0; w < n_meta_words; ++w) if (meta[w] != 0) meta_nonzero_data = 1;

    free(syndrome);
    free(meta);
    irrep_single_shot_code_free(&ss);

    /* meas-only meta should be NON-zero; data-only meta should be ZERO. */
    if (!meta_nonzero_meas) return 1;
    if (meta_nonzero_data) return 1;
    return 0;
}

/* Auto-lifter: irrep_single_shot_lift computes meta-checks via F₂
 * left-nullspace. Verifies:
 *   (a) verify_meta passes on the lifted code (so M_X · H_X = 0 and
 *       M_Z · H_Z = 0 by construction).
 *   (b) For 3D toric code at (Lx, Ly, Lz), the lifted M_X has at least
 *       1 row (the all-vertex redundancy) and M_Z has at least
 *       Lx·Ly·Lz - 1 rows (the cube redundancies modulo a global
 *       redundancy on T³). Exact rank check is L-dependent; the lower
 *       bound is robust. */
static int test_3d_toric_auto_lift(int Lx, int Ly, int Lz) {
    irrep_toric3d_params_t p;
    irrep_toric3d_init(&p, Lx, Ly, Lz);
    irrep_css_code_t css;
    if (irrep_toric3d_build(&p, &css) != IRREP_OK) return 1;

    irrep_single_shot_code_t ss;
    if (irrep_single_shot_lift(&css, &ss) != IRREP_OK) {
        irrep_css_code_free(&css);
        return 1;
    }
    int rc = 0;
    if (irrep_single_shot_verify_meta(&ss) != IRREP_OK) rc = 1;

    /* X-side: at least 1 redundancy on T³ (sum of all vertex stabs = 0). */
    if (ss.M_X.n_rows < 1) {
        fprintf(stderr, "FAIL auto-lift X: expected ≥1 meta row, got %d\n",
                ss.M_X.n_rows);
        rc = 1;
    }
    /* Z-side: cube redundancies, expected count = Lx·Ly·Lz - 1 (one
     * cube redundancy per cube, minus one global redundancy on T³). */
    int expected_Z_min = Lx * Ly * Lz - 1;
    if (ss.M_Z.n_rows < expected_Z_min) {
        fprintf(stderr, "FAIL auto-lift Z: expected ≥%d meta rows, got %d\n",
                expected_Z_min, ss.M_Z.n_rows);
        rc = 1;
    }
    irrep_single_shot_code_free(&ss);
    irrep_css_code_free(&css);
    return rc;
}

/* The Steane [[7, 1, 3]] code has no inter-stabilizer redundancy
 * (rank(H_X) = m_X = 3, rank(H_Z) = m_Z = 3), so the lift returns
 * M_X and M_Z each with 0 rows. */
static int test_steane_no_redundancy(void) {
    irrep_css_code_t css;
    if (irrep_css_code_new(&css, 7, 3, 3) != IRREP_OK) return 1;
    /* Steane stabilizers (H_X = H_Z = parity matrix of [7,4,3] Hamming):
     *   row 0: {0, 2, 4, 6}
     *   row 1: {1, 2, 5, 6}
     *   row 2: {3, 4, 5, 6}  */
    const int row0[4] = { 0, 2, 4, 6 };
    const int row1[4] = { 1, 2, 5, 6 };
    const int row2[4] = { 3, 4, 5, 6 };
    for (int i = 0; i < 4; ++i) {
        irrep_parity_matrix_set(&css.H_X, 0, row0[i]);
        irrep_parity_matrix_set(&css.H_Z, 0, row0[i]);
        irrep_parity_matrix_set(&css.H_X, 1, row1[i]);
        irrep_parity_matrix_set(&css.H_Z, 1, row1[i]);
        irrep_parity_matrix_set(&css.H_X, 2, row2[i]);
        irrep_parity_matrix_set(&css.H_Z, 2, row2[i]);
    }
    if (irrep_css_code_verify(&css) != IRREP_OK) {
        irrep_css_code_free(&css);
        return 1;
    }
    irrep_single_shot_code_t ss;
    if (irrep_single_shot_lift(&css, &ss) != IRREP_OK) {
        irrep_css_code_free(&css);
        return 1;
    }
    int rc = 0;
    if (ss.M_X.n_rows != 0) {
        fprintf(stderr, "FAIL Steane M_X: expected 0 rows, got %d\n",
                ss.M_X.n_rows);
        rc = 1;
    }
    if (ss.M_Z.n_rows != 0) {
        fprintf(stderr, "FAIL Steane M_Z: expected 0 rows, got %d\n",
                ss.M_Z.n_rows);
        rc = 1;
    }
    /* Trivially M_X · H_X = 0 and M_Z · H_Z = 0 when M has 0 rows. */
    if (irrep_single_shot_verify_meta(&ss) != IRREP_OK) rc = 1;
    irrep_single_shot_code_free(&ss);
    irrep_css_code_free(&css);
    return rc;
}

int main(void) {
    int rc = 0;
    int sizes[][3] = { {2, 2, 2}, {3, 2, 2}, {3, 3, 3} };
    for (size_t k = 0; k < sizeof(sizes)/sizeof(sizes[0]); ++k) {
        int Lx = sizes[k][0], Ly = sizes[k][1], Lz = sizes[k][2];
        if (test_3d_toric_meta_verify(Lx, Ly, Lz))
            { fprintf(stderr, "FAIL test_3d_toric_meta_verify(%d,%d,%d)\n", Lx, Ly, Lz); rc = 1; }
        if (test_meta_syndrome_zero(Lx, Ly, Lz))
            { fprintf(stderr, "FAIL test_meta_syndrome_zero(%d,%d,%d)\n", Lx, Ly, Lz); rc = 1; }
        if (test_meta_syndrome_decoupling(Lx, Ly, Lz))
            { fprintf(stderr, "FAIL test_meta_syndrome_decoupling(%d,%d,%d)\n", Lx, Ly, Lz); rc = 1; }
        if (test_3d_toric_auto_lift(Lx, Ly, Lz))
            { fprintf(stderr, "FAIL test_3d_toric_auto_lift(%d,%d,%d)\n", Lx, Ly, Lz); rc = 1; }
    }
    if (test_steane_no_redundancy())
        { fprintf(stderr, "FAIL test_steane_no_redundancy\n"); rc = 1; }
    return rc;
}
