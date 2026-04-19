/* SPDX-License-Identifier: MIT */
/* libFuzzer entry point for the 1.3 lattice + space_group builders.
 *
 * Reads two bytes of input as (kind, L) and one more byte as a bit-packed
 * field for Lx/Ly aspect ratio and wallpaper group. Exercises the site-
 * count, bond-list, permutation, and free paths on every valid combination
 * plus the rejection paths for invalid ones.
 *
 * Build with:
 *   clang -fsanitize=fuzzer,address -O1 -g -Iinclude \
 *         tests/fuzz/fuzz_lattice_build.c build/lib/liblibirrep.a \
 *         -lm -o fuzz_lattice_build
 *
 * Run:  ./fuzz_lattice_build -max_total_time=3600
 */

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/lattice.h>
#include <irrep/space_group.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 3) return 0;
    int kind_byte = data[0] % 5;   /* covers [0, 4] — one invalid value at 4 */
    int L         = (data[1] % 9) + 1;   /* 1..9 — 1 triggers the Lx<2 reject */
    int asp       = data[2] & 0x3;       /* 0..3: 00 square, 01/10 rect, 11 square */
    int wp        = (data[2] >> 2) & 0x3; /* 0..3 — includes one invalid */
    int Lx = L, Ly = L;
    if (asp == 1)      Ly = L + 1;
    else if (asp == 2) Lx = L + 1;

    irrep_lattice_t *lat = irrep_lattice_build((irrep_lattice_kind_t)kind_byte, Lx, Ly);
    if (!lat) return 0;

    /* Exercise query surface (no asserts; fuzzer catches UB / leaks / OOB). */
    (void)irrep_lattice_num_sites(lat);
    (void)irrep_lattice_num_bonds_nn(lat);
    (void)irrep_lattice_num_bonds_nnn(lat);
    double a1[2], a2[2], b1[2], b2[2];
    irrep_lattice_primitive_vectors (lat, a1, a2);
    irrep_lattice_reciprocal_vectors(lat, b1, b2);

    /* Bond list fill */
    int nb = irrep_lattice_num_bonds_nn(lat);
    if (nb > 0) {
        int *ii = malloc(sizeof(int) * nb);
        int *jj = malloc(sizeof(int) * nb);
        if (ii && jj) irrep_lattice_fill_bonds_nn(lat, ii, jj);
        free(ii); free(jj);
    }

    /* Site position lookup for a sample of sites */
    int n = irrep_lattice_num_sites(lat);
    for (int s = 0; s < n; s += (n / 4 + 1)) {
        double xy[2];
        (void)irrep_lattice_site_position(lat, s, xy);
        int ix, iy;
        (void)irrep_lattice_cell_of(lat, s, &ix, &iy);
        (void)irrep_lattice_sublattice_of(lat, s);
    }

    /* Space-group build + application (if the wallpaper group matches the
     * lattice; the library internally rejects incompatible pairs). */
    irrep_space_group_t *G = irrep_space_group_build(lat, (irrep_wallpaper_t)wp);
    if (G) {
        int order = irrep_space_group_order(G);
        int ns    = irrep_space_group_num_sites(G);
        int *perm = malloc(sizeof(int) * ns);
        if (perm) {
            for (int g = 0; g < order; g += (order / 8 + 1)) {
                irrep_space_group_permutation(G, g, perm);
                irrep_space_group_permutation_inverse(G, g, perm);
                (void)irrep_space_group_apply(G, g, 0);
            }
            free(perm);
        }
        irrep_space_group_free(G);
    }

    irrep_lattice_free(lat);
    return 0;
}
