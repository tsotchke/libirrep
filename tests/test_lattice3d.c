/* SPDX-License-Identifier: MIT */
/* 3D Bravais-lattice regression tests.
 *
 * Asserted invariants (per family):
 *   - sites_per_cell, num_sites, num_cells
 *   - NN / NNN distances match the closed-form geometric values
 *   - Coordination number: NN bond list yields the expected per-site
 *     coordination (6 SC, 8 BCC, 12 FCC, 4 Diamond)
 *   - Bond endpoints are valid sites with `i < j` (canonical)
 *   - PBC translation: translate(site, Lx, 0, 0) == site, etc.
 *   - Reciprocal lattice satisfies a_i · b_j = 2π δ_ij
 *   - Site position round-trips through site_index → site_position → cell_of */

#include <irrep/lattice3d.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
        double da_ = (a), db_ = (b);                                                               \
        if (fabs(da_ - db_) > (tol)) {                                                             \
            fprintf(stderr, "  FAIL  %s:%d  %s  got %.12g vs %.12g\n", __FILE__, __LINE__, msg,    \
                    da_, db_);                                                                     \
            ++failed;                                                                              \
        }                                                                                          \
    } while (0)

typedef struct {
    irrep_lattice3d_kind_t kind;
    const char            *name;
    int                    spc;
    int                    coord_nn;
    int                    coord_nnn;
    double                 nn;
    double                 nnn;
} expected_t;

static const expected_t TABLE[] = {
    {IRREP_LATTICE3D_SC, "SC", 1, 6, 12, 1.0, 1.41421356237309504880},
    {IRREP_LATTICE3D_BCC, "BCC", 2, 8, 6, 0.86602540378443864676, 1.0},
    {IRREP_LATTICE3D_FCC, "FCC", 4, 12, 6, 0.70710678118654752440, 1.0},
    {IRREP_LATTICE3D_DIAMOND, "Diamond", 8, 4, 12, 0.43301270189221932338,
     0.70710678118654752440},
    {IRREP_LATTICE3D_PYROCHLORE, "Pyrochlore", 16, 6, 12, 0.35355339059327376220,
     0.61237243569579452455}};

static void test_family(const expected_t *E, int Lx, int Ly, int Lz) {
    fprintf(stderr, "  %s (%d×%d×%d)\n", E->name, Lx, Ly, Lz);
    irrep_lattice3d_t *L = irrep_lattice3d_build(E->kind, Lx, Ly, Lz);
    ASSERT(L != NULL, "build");
    if (!L)
        return;

    int spc = irrep_lattice3d_sites_per_cell(L);
    int ncell = irrep_lattice3d_num_cells(L);
    int nsite = irrep_lattice3d_num_sites(L);
    ASSERT(spc == E->spc, "sites_per_cell");
    ASSERT(ncell == Lx * Ly * Lz, "num_cells");
    ASSERT(nsite == ncell * spc, "num_sites");
    ASSERT_NEAR(irrep_lattice3d_nn_distance(L), E->nn, 1e-12, "nn_distance");
    ASSERT_NEAR(irrep_lattice3d_nnn_distance(L), E->nnn, 1e-12, "nnn_distance");

    /* Coordination number from NN bond list. */
    int n_nn = irrep_lattice3d_num_bonds_nn(L);
    int n_nnn = irrep_lattice3d_num_bonds_nnn(L);
    int *bi = malloc((size_t)(n_nn > n_nnn ? n_nn : n_nnn) * sizeof(int));
    int *bj = malloc((size_t)(n_nn > n_nnn ? n_nn : n_nnn) * sizeof(int));
    int *deg = calloc((size_t)nsite, sizeof(int));

    irrep_lattice3d_fill_bonds_nn(L, bi, bj);
    for (int b = 0; b < n_nn; ++b) {
        ASSERT(bi[b] >= 0 && bi[b] < nsite, "nn endpoint i in range");
        ASSERT(bj[b] >= 0 && bj[b] < nsite, "nn endpoint j in range");
        ASSERT(bi[b] < bj[b], "nn canonical i < j");
        deg[bi[b]]++;
        deg[bj[b]]++;
    }
    int min_deg = deg[0], max_deg = deg[0];
    for (int s = 1; s < nsite; ++s) {
        if (deg[s] < min_deg)
            min_deg = deg[s];
        if (deg[s] > max_deg)
            max_deg = deg[s];
    }
    ASSERT(min_deg == E->coord_nn && max_deg == E->coord_nn, "NN coordination uniform");

    /* NNN coordination — only meaningful when not aliased by PBC; for our
     * cluster sizes we expect uniform coordination. */
    for (int s = 0; s < nsite; ++s)
        deg[s] = 0;
    irrep_lattice3d_fill_bonds_nnn(L, bi, bj);
    for (int b = 0; b < n_nnn; ++b) {
        ASSERT(bi[b] < bj[b], "nnn canonical i < j");
        deg[bi[b]]++;
        deg[bj[b]]++;
    }
    min_deg = deg[0];
    max_deg = deg[0];
    for (int s = 1; s < nsite; ++s) {
        if (deg[s] < min_deg)
            min_deg = deg[s];
        if (deg[s] > max_deg)
            max_deg = deg[s];
    }
    ASSERT(min_deg == E->coord_nnn && max_deg == E->coord_nnn, "NNN coordination uniform");

    free(bi);
    free(bj);
    free(deg);

    /* Bond geometric distance check on a sample. */
    {
        int ns = (n_nn < 32 ? n_nn : 32);
        int *si = malloc((size_t)n_nn * sizeof(int));
        int *sj = malloc((size_t)n_nn * sizeof(int));
        irrep_lattice3d_fill_bonds_nn(L, si, sj);
        for (int b = 0; b < ns; ++b) {
            double pi[3], pj[3];
            irrep_lattice3d_site_position(L, si[b], pi);
            irrep_lattice3d_site_position(L, sj[b], pj);
            /* PBC may extend the bond across the cluster — take the minimum
             * image distance. */
            double best = 1e300;
            double a1[3], a2[3], a3[3];
            irrep_lattice3d_primitive_vectors(L, a1, a2, a3);
            for (int wx = -1; wx <= 1; ++wx)
                for (int wy = -1; wy <= 1; ++wy)
                    for (int wz = -1; wz <= 1; ++wz) {
                        double dx = pj[0] - pi[0] + wx * Lx * a1[0] + wy * Ly * a2[0] +
                                    wz * Lz * a3[0];
                        double dy = pj[1] - pi[1] + wx * Lx * a1[1] + wy * Ly * a2[1] +
                                    wz * Lz * a3[1];
                        double dz = pj[2] - pi[2] + wx * Lx * a1[2] + wy * Ly * a2[2] +
                                    wz * Lz * a3[2];
                        double r = sqrt(dx * dx + dy * dy + dz * dz);
                        if (r < best)
                            best = r;
                    }
            ASSERT_NEAR(best, E->nn, 1e-9, "bond cartesian length matches NN distance");
        }
        free(si);
        free(sj);
    }

    /* PBC translation closure. */
    for (int s = 0; s < nsite && s < 8; ++s) {
        int t = irrep_lattice3d_translate(L, s, Lx, 0, 0);
        ASSERT(t == s, "translate(Lx, 0, 0) closes");
        t = irrep_lattice3d_translate(L, s, 0, Ly, 0);
        ASSERT(t == s, "translate(0, Ly, 0) closes");
        t = irrep_lattice3d_translate(L, s, 0, 0, Lz);
        ASSERT(t == s, "translate(0, 0, Lz) closes");
    }

    /* site_index round-trip. */
    for (int s = 0; s < nsite; s += (nsite / 8 > 0 ? nsite / 8 : 1)) {
        int ix, iy, iz;
        irrep_lattice3d_cell_of(L, s, &ix, &iy, &iz);
        int sub = irrep_lattice3d_sublattice_of(L, s);
        int t = irrep_lattice3d_site_index(L, ix, iy, iz, sub);
        ASSERT(t == s, "site_index round-trip");
    }

    /* Reciprocal lattice. */
    {
        double a1[3], a2[3], a3[3], b1[3], b2[3], b3[3];
        irrep_lattice3d_primitive_vectors(L, a1, a2, a3);
        irrep_lattice3d_reciprocal_vectors(L, b1, b2, b3);
        double *as[3] = {a1, a2, a3};
        double *bs[3] = {b1, b2, b3};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double dot = as[i][0] * bs[j][0] + as[i][1] * bs[j][1] + as[i][2] * bs[j][2];
                double exp = (i == j) ? (2.0 * M_PI) : 0.0;
                ASSERT_NEAR(dot, exp, 1e-12, "a_i · b_j = 2π δ_ij");
            }
    }

    /* k-grid sanity: contains origin. */
    {
        double *kx = malloc((size_t)ncell * sizeof(double));
        double *ky = malloc((size_t)ncell * sizeof(double));
        double *kz = malloc((size_t)ncell * sizeof(double));
        irrep_lattice3d_k_grid(L, kx, ky, kz);
        ASSERT_NEAR(kx[0], 0.0, 1e-12, "k_grid[0].x = 0");
        ASSERT_NEAR(ky[0], 0.0, 1e-12, "k_grid[0].y = 0");
        ASSERT_NEAR(kz[0], 0.0, 1e-12, "k_grid[0].z = 0");
        free(kx);
        free(ky);
        free(kz);
    }

    irrep_lattice3d_free(L);
}

int main(void) {
    fprintf(stderr, "test_lattice3d:\n");

    /* Sizes chosen to keep PBC artefacts away while staying fast:
     *   SC needs ≥ 3³ to keep NN bond count = 3·N_sites
     *   BCC, FCC, Diamond fine at 2³ but use 3³ for NNN robustness on FCC. */
    test_family(&TABLE[0], 3, 3, 3); /* SC */
    test_family(&TABLE[1], 3, 3, 3); /* BCC */
    test_family(&TABLE[2], 3, 3, 3); /* FCC */
    test_family(&TABLE[3], 3, 3, 3); /* Diamond */
    test_family(&TABLE[4], 2, 2, 2); /* Pyrochlore — 2³ cell already 128 sites */

    fprintf(stderr, "  %d / %d assertions passed\n", total - failed, total);
    return failed == 0 ? 0 : 1;
}
