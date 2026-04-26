/* SPDX-License-Identifier: MIT */
/* 3D Heisenberg AFM ground state via Lanczos on small cubic clusters.
 *
 * Demonstrates the new `lattice3d` module flowing into the existing
 * Hamiltonian and Lanczos primitives:
 *
 *   irrep_lattice3d_build         -> geometry
 *   irrep_lattice3d_fill_bonds_nn -> NN bond list
 *   irrep_heisenberg_new          -> H = J · Σ_<ij> S_i · S_j
 *   irrep_lanczos_eigvals         -> ground-state eigenvalue
 *
 * Cluster choices:
 *   SC  2×2×2  (N=8, dim=256)   — cube graph, bipartite, classic AFM benchmark
 *   BCC 2×2×2  (N=16, dim=65536) — bipartite (A vs body-centered B), 4-fold
 *                                    coordination after PBC dedup on 2-cell axes
 *
 * On a 2-cell-per-axis cluster, periodic-boundary-condition wrap places
 * "+x" and "-x" on the same image; the canonical bond list deduplicates
 * to a single bond per pair, so the resulting graph is a finite undirected
 * graph (cube for SC, the 16-vertex graph for BCC) rather than the
 * thermodynamic-limit lattice. Use 3³ or larger when you want the
 * thermodynamic-limit bond multiplicities. */

#include <irrep/hamiltonian.h>
#include <irrep/lattice3d.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static const char *kind_name(irrep_lattice3d_kind_t k) {
    switch (k) {
    case IRREP_LATTICE3D_SC:
        return "SC";
    case IRREP_LATTICE3D_BCC:
        return "BCC";
    case IRREP_LATTICE3D_FCC:
        return "FCC";
    case IRREP_LATTICE3D_DIAMOND:
        return "Diamond";
    case IRREP_LATTICE3D_PYROCHLORE:
        return "Pyrochlore";
    }
    return "?";
}

static void run(irrep_lattice3d_kind_t kind, int Lx, int Ly, int Lz, int max_iters) {
    irrep_lattice3d_t *L = irrep_lattice3d_build(kind, Lx, Ly, Lz);
    if (!L) {
        fprintf(stderr, "  build failed for %s\n", kind_name(kind));
        return;
    }
    int N = irrep_lattice3d_num_sites(L);
    int nb = irrep_lattice3d_num_bonds_nn(L);

    if (N > 24) {
        fprintf(stderr, "  %s %d×%d×%d: N=%d sites — full ED beyond 2^24 dim, skipping\n",
                kind_name(kind), Lx, Ly, Lz, N);
        irrep_lattice3d_free(L);
        return;
    }

    int *bi = malloc((size_t)nb * sizeof(int));
    int *bj = malloc((size_t)nb * sizeof(int));
    irrep_lattice3d_fill_bonds_nn(L, bi, bj);

    irrep_heisenberg_t *H = irrep_heisenberg_new(N, nb, bi, bj, /*J=*/1.0);
    long long dim = irrep_heisenberg_dim(H);

    /* Deterministic Sz=0 spread seed — non-pathological overlap with the GS
     * for any bipartite cluster (and most non-bipartite ones). */
    double _Complex *seed = malloc((size_t)dim * sizeof(*seed));
    double snorm = 0;
    for (long long s = 0; s < dim; ++s) {
        seed[s] = 0.1 * sin(0.37 * s) + I * 0.05 * cos(0.23 * s);
        snorm += creal(seed[s]) * creal(seed[s]) + cimag(seed[s]) * cimag(seed[s]);
    }
    snorm = sqrt(snorm);
    for (long long s = 0; s < dim; ++s)
        seed[s] /= snorm;

    double t0 = now_sec();
    double E0 = 0;
    irrep_status_t st =
        irrep_lanczos_eigvals(irrep_heisenberg_apply, H, dim, 1, max_iters, seed, &E0);
    double dt = now_sec() - t0;

    if (st != IRREP_OK) {
        fprintf(stderr, "  Lanczos returned status %d\n", (int)st);
    } else {
        printf("  %-10s %d×%d×%d  N=%-3d  bonds=%-4d  dim=%-9lld  E_0 = %+.6f J  "
               "(E_0/N_bonds = %+.6f, %.2fs)\n",
               kind_name(kind), Lx, Ly, Lz, N, nb, dim, E0, E0 / nb, dt);
    }

    free(seed);
    free(bi);
    free(bj);
    irrep_heisenberg_free(H);
    irrep_lattice3d_free(L);
}

int main(void) {
    printf("=== libirrep — 3D Heisenberg AFM ground state via Lanczos ===\n");
    printf("    H = J · Σ_<ij> S_i · S_j ,  S = ½,  J = 1\n\n");

    /* SC 2³ = 8 sites: cube graph (each site has 3 distinct NN after PBC
     * dedup on the 2-cell axes). Heisenberg AFM on the cube has GS at the
     * standard small-cluster value. */
    run(IRREP_LATTICE3D_SC, 2, 2, 2, 80);

    /* BCC 2³ = 16 sites: dim = 65536. BCC is bipartite (A sites vs. body-
     * centered B sites), so the AFM ground state is an exact singlet on
     * a finite cluster up to finite-size corrections. */
    run(IRREP_LATTICE3D_BCC, 2, 2, 2, 200);

    /* SC 4×2×2 = 16 sites: same dim as BCC 2³ but a different graph
     * (rectangular cube prism). Useful as a cross-check. */
    run(IRREP_LATTICE3D_SC, 4, 2, 2, 200);

    return 0;
}
