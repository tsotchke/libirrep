/* SPDX-License-Identifier: MIT */
/* Heisenberg AFM ground state of the 16-site pyrochlore cluster.
 *
 * Pyrochlore is the canonical 3D frustrated antiferromagnet — the
 * sublattice of corner-sharing tetrahedra in the spinel B-site family
 * (Tb₂Ti₂O₇, Yb₂Ti₂O₇, Dy₂Ti₂O₇, ZnCr₂O₄). It hosts spin-ice and U(1)
 * quantum-spin-liquid candidates depending on anisotropy and longer-
 * range exchange. The pure NN Heisenberg case has no closed-form ground
 * state; small-cluster ED is one of the few benchmarks available.
 *
 * Cluster: pyrochlore 1×1×1 = 16 sites, 6 NN per site, 48 NN bonds. PBC
 * along all three cubic axes — cell offsets at L=1 wrap to a single
 * image, but the pyrochlore basis preserves NN coordination because
 * neighbours of a site live on different sublattices within the same
 * conventional cubic cell.
 *
 * Hilbert space: 2¹⁶ = 65536. Lanczos with full reorthogonalisation
 * picks out the few lowest eigenvalues (the "tower" structure expected
 * for a frustrated 3D system). */

#include <irrep/hamiltonian.h>
#include <irrep/lattice3d.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define K_WANTED 6 /* low-lying eigenvalues — useful for the gap structure */

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(void) {
    printf("=== libirrep — pyrochlore 16-site Heisenberg AFM ===\n");
    printf("    H = J · Σ_<ij> S_i · S_j ,  S = ½,  J = 1\n\n");

    irrep_lattice3d_t *L = irrep_lattice3d_build(IRREP_LATTICE3D_PYROCHLORE, 1, 1, 1);
    if (!L) {
        fprintf(stderr, "  pyrochlore 1×1×1 build failed\n");
        return 1;
    }
    int N = irrep_lattice3d_num_sites(L);
    int nb = irrep_lattice3d_num_bonds_nn(L);
    int *bi = malloc((size_t)nb * sizeof(int));
    int *bj = malloc((size_t)nb * sizeof(int));
    irrep_lattice3d_fill_bonds_nn(L, bi, bj);

    printf("  Geometry: pyrochlore 1×1×1\n");
    printf("    N sites          = %d\n", N);
    printf("    NN bonds         = %d   (%.1f per site, expected 6.0)\n", nb, 2.0 * nb / N);
    printf("    NN distance      = %.6f   (= √2/4 — tetrahedral edge)\n",
           irrep_lattice3d_nn_distance(L));

    /* Print one tetrahedron: sites 0..3 form the "up" tetrahedron of FCC
     * sublattice 0 — pairwise NN at √2/4. */
    double r0[3], r1[3], r2[3], r3[3];
    irrep_lattice3d_site_position(L, 0, r0);
    irrep_lattice3d_site_position(L, 1, r1);
    irrep_lattice3d_site_position(L, 2, r2);
    irrep_lattice3d_site_position(L, 3, r3);
    printf("    'Up' tetrahedron (sites 0..3):\n");
    printf("      site 0 r = (%+.3f, %+.3f, %+.3f)\n", r0[0], r0[1], r0[2]);
    printf("      site 1 r = (%+.3f, %+.3f, %+.3f)\n", r1[0], r1[1], r1[2]);
    printf("      site 2 r = (%+.3f, %+.3f, %+.3f)\n", r2[0], r2[1], r2[2]);
    printf("      site 3 r = (%+.3f, %+.3f, %+.3f)\n", r3[0], r3[1], r3[2]);

    irrep_heisenberg_t *H = irrep_heisenberg_new(N, nb, bi, bj, /*J=*/1.0);
    long long           dim = irrep_heisenberg_dim(H);

    /* Sz=0 spread seed — the GS of the AFM is a singlet, so it lives in
     * Sz=0. A seed concentrated on Sz=0 configurations is non-pathological. */
    double _Complex *seed = malloc((size_t)dim * sizeof(*seed));
    double           sn = 0;
    int              target_pop = N / 2;
    for (long long s = 0; s < dim; ++s) {
        if (__builtin_popcountll((unsigned long long)s) == target_pop) {
            seed[s] = 0.1 * sin(0.37 * s) + I * 0.05 * cos(0.23 * s);
        } else {
            seed[s] = 0;
        }
        sn += creal(seed[s]) * creal(seed[s]) + cimag(seed[s]) * cimag(seed[s]);
    }
    sn = sqrt(sn);
    for (long long s = 0; s < dim; ++s)
        seed[s] /= sn;

    printf("\n  Lanczos: dim=%lld, k_wanted=%d, max_iters=300 (with reorth)\n", dim, K_WANTED);
    fflush(stdout);

    double t0 = now_sec();
    double eigs[K_WANTED];
    irrep_status_t st = irrep_lanczos_eigvals_reorth(irrep_heisenberg_apply, H, dim, K_WANTED, 300,
                                                     seed, eigs);
    double dt = now_sec() - t0;
    if (st != IRREP_OK) {
        fprintf(stderr, "  Lanczos failed, status=%d\n", (int)st);
        return 1;
    }

    printf("\n  Lowest %d eigenvalues (J units):\n", K_WANTED);
    for (int k = 0; k < K_WANTED; ++k) {
        printf("    E_%d = %+12.6f   (E_%d / N_sites = %+8.6f, "
               "E_%d / N_bonds = %+8.6f)\n",
               k, eigs[k], k, eigs[k] / N, k, eigs[k] / nb);
    }
    double gap_01 = eigs[1] - eigs[0];
    double gap_02 = eigs[2] - eigs[0];
    printf("\n  Lowest gaps:\n");
    printf("    Δ_01 = E_1 - E_0 = %.6f J\n", gap_01);
    printf("    Δ_02 = E_2 - E_0 = %.6f J\n", gap_02);
    printf("\n  Lanczos wall-clock = %.2fs\n", dt);

    /* For comparison: 4-site single tetrahedron (S=½ Heisenberg AFM on
     * K_4) has GS = -3/2 J = -0.5 J / bond. With 6 corner-sharing
     * tetrahedra in our cluster (16 sites / 4 sites per tetra ≈ 4
     * corner-share-counted tetrahedra → 16 · 6 / 4 ≈ 24 — actually
     * exactly 8 tetrahedra in the 16-site cluster, each sharing 4 corners
     * × 1 bond per shared corner — the geometry is non-trivial). The
     * "naive tetrahedron" lower bound on E₀ would be N_tetra · (-3/2 J)
     * = 8 · (-1.5) = -12 J. Frustration prevents saturating this. */
    printf("\n  Reference: lower bound from independent-tetrahedron decomposition\n");
    printf("    (8 tetrahedra × E_tetra = -3/2 J each) = -12 J\n");
    printf("    Actual E_0 = %+.6f J (saturation ratio: %.2f%%)\n", eigs[0],
           100.0 * eigs[0] / -12.0);

    free(seed);
    free(bi);
    free(bj);
    irrep_heisenberg_free(H);
    irrep_lattice3d_free(L);
    return 0;
}
