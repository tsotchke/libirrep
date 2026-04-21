/* SPDX-License-Identifier: MIT */
/* CLOCK_MONOTONIC is POSIX.1b, not strict C11. */
#define _POSIX_C_SOURCE 199309L
/* End-to-end configuration-space A₁ projection on a 6×6 × 3 = 108-site
 * kagome cluster — the target lattice for the symmetric neural-
 * quantum-state pipeline (kagome Heisenberg S = ½ ground-state nature,
 * open since Yan–Huse–White 2011).
 *
 * The example demonstrates the full 1.3 pipeline as a downstream VMC
 * driver would assemble it, end to end in ≲ 100 LOC:
 *
 *   1. Build the 6×6 kagome lattice (216 NN bonds).
 *   2. Attach its p6mm space group: 432 permutation elements on 108 sites.
 *   3. Declare a toy "reference wavefunction" — a classical Ising product
 *      state, specifically a √3×√3 Néel pattern (the order expected on a
 *      conventionally-ordered kagome).
 *   4. For a single sampled spin configuration σ, enumerate the full orbit
 *      {g·σ : g ∈ G} via the cached pullback permutations.
 *   5. Evaluate the reference amplitude at each orbit image.
 *   6. Reduce to the A₁-projected amplitude with the totally-symmetric
 *      character row (@ref irrep_sg_project_A1).
 *   7. Print a tiny report: cluster size, group order, projection time,
 *      projected amplitude.
 *
 * A real NQS driver substitutes step 5 with a neural-network forward pass;
 * everything else stays verbatim.
 *
 * Compile with `make examples` or
 *   cc -I include -L build/lib -llibirrep examples/kagome_a1_projection.c
 */

#include <irrep/config_project.h>
#include <irrep/lattice.h>
#include <irrep/space_group.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Reference "wavefunction": sign-structured Ising product state. For a
 * classical Néel state σ on kagome, the amplitude is ±1 depending on the
 * ordering pattern; for the sake of the example we use a simple linear
 * combination sensitive to the site-0 phase, so that most orbit images
 * contribute the same value and the A₁ projector delivers a non-trivial
 * scalar. */
static double _Complex toy_amplitude(int num_sites, const double *sigma) {
    double acc = 0.0;
    for (int s = 0; s < num_sites; ++s) acc += sigma[s];
    /* acc / √N as a simple magnetisation-like observable. */
    return (acc / sqrt((double)num_sites)) + 0.0 * I;
}

static double now_ms_(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e3 + (double)ts.tv_nsec * 1e-6;
}

int main(void) {
    /* 1. Lattice. */
    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 6, 6);
    if (!L) { fprintf(stderr, "lattice build failed\n"); return 1; }
    int num_sites = irrep_lattice_num_sites(L);

    /* 2. Space group. */
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
    if (!G) { fprintf(stderr, "space group build failed\n"); irrep_lattice_free(L); return 1; }
    int order = irrep_space_group_order(G);

    printf("kagome %d×%d = %d sites, p6mm order = %d\n",
           irrep_lattice_Lx(L), irrep_lattice_Ly(L), num_sites, order);
    printf("NN bonds: %d, NNN bonds: %d\n",
           irrep_lattice_num_bonds_nn(L),
           irrep_lattice_num_bonds_nnn(L));

    /* 3. Toy configuration: ±1 per site, sign determined by sublattice. */
    double *sigma = malloc((size_t)num_sites * sizeof(double));
    for (int s = 0; s < num_sites; ++s) {
        int sub = irrep_lattice_sublattice_of(L, s);
        sigma[s] = (sub == 0) ? +1.0 : (sub == 1 ? -1.0 : +1.0);
    }

    /* 4. Enumerate the orbit. */
    double t_orbit_0 = now_ms_();
    double *orbit = malloc((size_t)order * num_sites * sizeof(double));
    irrep_sg_enumerate_orbit(G, sigma, orbit);
    double t_orbit_1 = now_ms_();

    /* 5. Evaluate the reference amplitude on each image. */
    double t_eval_0 = now_ms_();
    double _Complex *amps = malloc((size_t)order * sizeof(double _Complex));
    for (int g = 0; g < order; ++g) {
        amps[g] = toy_amplitude(num_sites, orbit + (size_t)g * num_sites);
    }
    double t_eval_1 = now_ms_();

    /* 6. A₁ projection. */
    double t_proj_0 = now_ms_();
    double _Complex proj = irrep_sg_project_A1(G, amps);
    double t_proj_1 = now_ms_();

    /* 7. Report. */
    printf("orbit enumeration:   %.3f ms  (432 × 108 = %d doubles copied)\n",
           t_orbit_1 - t_orbit_0, order * num_sites);
    printf("amplitude evaluate:  %.3f ms  (432 calls to toy_amplitude)\n",
           t_eval_1 - t_eval_0);
    printf("A₁ projection:       %.3f ms  (432-term dot product)\n",
           t_proj_1 - t_proj_0);
    printf("P_{A₁} ψ(σ) = %+.12f  %+.12f i\n", creal(proj), cimag(proj));

    /* Sanity: the toy amplitude only depends on the sum Σσ, which is a
     * translation-invariant and point-group-invariant scalar. So every
     * orbit image has the same amplitude, and the A₁ projection returns
     * that amplitude unchanged. */
    double _Complex expected = toy_amplitude(num_sites, sigma);
    double drift = cabs(proj - expected);
    printf("drift from translation-invariant reference: %.3e\n", drift);

    free(amps);
    free(orbit);
    free(sigma);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
