/* SPDX-License-Identifier: MIT */
/* Demonstration of the 3D lattice module.
 *
 * Builds each supported family at a modest cluster size and prints:
 *   - sites/cell, total sites, NN distance, NN bond count
 *   - cartesian positions of the first few sites
 *   - one sample NN bond with its end-point cartesian coordinates
 *
 * Bond lists from this module plug directly into the existing Heisenberg
 * builder (`irrep_heisenberg_new(N, n_bonds, bi, bj, J)`); the geometry is
 * the only piece that was missing for 3D physics work. */

#include <irrep/lattice3d.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const char *kind_name(irrep_lattice3d_kind_t k) {
    switch (k) {
    case IRREP_LATTICE3D_SC:
        return "Simple cubic";
    case IRREP_LATTICE3D_BCC:
        return "Body-centered cubic";
    case IRREP_LATTICE3D_FCC:
        return "Face-centered cubic";
    case IRREP_LATTICE3D_DIAMOND:
        return "Diamond";
    case IRREP_LATTICE3D_PYROCHLORE:
        return "Pyrochlore";
    }
    return "unknown";
}

static void demo(irrep_lattice3d_kind_t kind, int Lx, int Ly, int Lz) {
    irrep_lattice3d_t *L = irrep_lattice3d_build(kind, Lx, Ly, Lz);
    if (!L) {
        fprintf(stderr, "build failed for kind %d\n", (int)kind);
        return;
    }

    int spc = irrep_lattice3d_sites_per_cell(L);
    int N = irrep_lattice3d_num_sites(L);
    int nb = irrep_lattice3d_num_bonds_nn(L);
    int nb2 = irrep_lattice3d_num_bonds_nnn(L);

    printf("\n━ %s, %d×%d×%d  (sites/cell=%d, N=%d) ━\n", kind_name(kind), Lx, Ly, Lz, spc, N);
    printf("  NN  distance = %.6f, %d bonds (%.1f per site)\n",
           irrep_lattice3d_nn_distance(L), nb, 2.0 * nb / N);
    printf("  NNN distance = %.6f, %d bonds (%.1f per site)\n",
           irrep_lattice3d_nnn_distance(L), nb2, 2.0 * nb2 / N);

    /* Sample positions. */
    int n_show = spc < 4 ? spc : 4;
    printf("  Sublattice positions (cell 0):\n");
    for (int s = 0; s < n_show; ++s) {
        double xyz[3];
        irrep_lattice3d_site_position(L, s, xyz);
        printf("    sub %d   r = (%+0.3f, %+0.3f, %+0.3f)\n", s, xyz[0], xyz[1], xyz[2]);
    }

    /* Sample NN bond. */
    int *bi = malloc(nb * sizeof(int));
    int *bj = malloc(nb * sizeof(int));
    irrep_lattice3d_fill_bonds_nn(L, bi, bj);
    if (nb > 0) {
        double ri[3], rj[3];
        irrep_lattice3d_site_position(L, bi[0], ri);
        irrep_lattice3d_site_position(L, bj[0], rj);
        printf("  NN bond[0]: site %d at (%+0.3f,%+0.3f,%+0.3f) -- site %d at "
               "(%+0.3f,%+0.3f,%+0.3f)\n",
               bi[0], ri[0], ri[1], ri[2], bj[0], rj[0], rj[1], rj[2]);
    }
    free(bi);
    free(bj);
    irrep_lattice3d_free(L);
}

int main(void) {
    printf("=== libirrep 3D lattice module — geometry demonstration ===\n");
    demo(IRREP_LATTICE3D_SC, 3, 3, 3);
    demo(IRREP_LATTICE3D_BCC, 3, 3, 3);
    demo(IRREP_LATTICE3D_FCC, 3, 3, 3);
    demo(IRREP_LATTICE3D_DIAMOND, 2, 2, 2);
    demo(IRREP_LATTICE3D_PYROCHLORE, 2, 2, 2);
    return 0;
}
