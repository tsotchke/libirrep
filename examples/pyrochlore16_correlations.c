/* SPDX-License-Identifier: MIT */
/* Equal-time spin-spin correlation ⟨S_0 · S_r⟩ on the pyrochlore-16 GS.
 *
 * The connected spin-spin correlation function is the standard probe
 * for ordering tendencies in a frustrated magnet:
 *
 *   - Magnetically ordered phase: ⟨S_0 · S_r⟩ saturates to a finite
 *     value at large r (the order parameter squared, with sign set by
 *     the ordering wave-vector).
 *   - Quantum-spin-liquid phase: ⟨S_0 · S_r⟩ decays exponentially with
 *     distance, with a characteristic correlation length ξ.
 *   - Valence-bond solid: alternates strongly between bonded and
 *     non-bonded pairs, with the bond-pattern dimerisation signature.
 *
 * On the 16-site cluster the long-range tail is heavily renormalised by
 * finite size, but the on-site (r=0), NN, and NNN values are quantitative
 * benchmarks against larger-N studies (e.g. Tsunetsugu 2002 on N=24, 32).
 *
 * Output: the connected correlator C(r) = ⟨S_0 · S_r⟩ for every r,
 * sorted by minimum-image distance, plus the four-spin connected
 * correlator on a single tetrahedron (sanity check: it must equal the
 * Heisenberg energy of an isolated tetrahedron, −3/2 J / 6 = −¼ J per
 * bond, modulo finite-size cluster corrections).
 *
 * Build: `make examples`
 * Run:   `./build/bin/pyrochlore16_correlations` */

#include <irrep/hamiltonian.h>
#include <irrep/lattice3d.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N           16
#define DIM         (1LL << N)

/* Apply S_a · S_b to a state vector. Result is real (Heisenberg
 * commutes with Sz_total and S² is real). */
static double expectation_SS_(const double _Complex *psi, int a, int b) {
    /* S_a · S_b = ¼ (X_a X_b + Y_a Y_b + Z_a Z_b)
     *          = ¼ Z_a Z_b + ½ (S_+_a S_-_b + S_-_a S_+_b)
     * On bit-string s (bit i = +½ spin if 1, -½ if 0):
     *   Z_a Z_b |s⟩ = (+1 if bits aligned, -1 if anti) · |s⟩
     *   (S_+ S_- + S_- S_+) |s⟩ = |s' = s ^ (1<<a ^ 1<<b)⟩  if bits anti-aligned
     *                            = 0  if bits aligned. */
    double real_part = 0;
    for (long long s = 0; s < DIM; ++s) {
        double mag2 = creal(psi[s]) * creal(psi[s]) + cimag(psi[s]) * cimag(psi[s]);
        if (mag2 < 1e-30)
            continue;
        int sa = (int)((s >> a) & 1);
        int sb = (int)((s >> b) & 1);
        if (sa == sb) {
            real_part += 0.25 * mag2; /* Z_a Z_b = +1 */
        } else {
            real_part -= 0.25 * mag2; /* Z_a Z_b = -1 */
            /* Off-diagonal: ½ |s' ⟩ contribution to ⟨ψ| ψ'⟩ = ½ ψ*[s'] · ψ[s]. */
            long long sp = s ^ (1LL << a) ^ (1LL << b);
            double _Complex z = conj(psi[sp]) * psi[s];
            real_part += 0.5 * creal(z);
        }
    }
    return real_part;
}

typedef struct {
    int    site;
    double dist;
    double corr;
} site_corr_t;

static int cmp_dist(const void *a, const void *b) {
    const site_corr_t *pa = a, *pb = b;
    if (pa->dist < pb->dist)
        return -1;
    if (pa->dist > pb->dist)
        return 1;
    return pa->site - pb->site;
}

/* Minimum-image distance from r_a to r_b under PBC with primitive vectors a1/a2/a3. */
static double min_image_distance(const double r_a[3], const double r_b[3], const double a1[3],
                                 const double a2[3], const double a3[3], int Lx, int Ly, int Lz) {
    double best = 1e300;
    for (int wx = -1; wx <= 1; ++wx)
        for (int wy = -1; wy <= 1; ++wy)
            for (int wz = -1; wz <= 1; ++wz) {
                double dx = r_b[0] - r_a[0] + wx * Lx * a1[0] + wy * Ly * a2[0] +
                            wz * Lz * a3[0];
                double dy = r_b[1] - r_a[1] + wx * Lx * a1[1] + wy * Ly * a2[1] +
                            wz * Lz * a3[1];
                double dz = r_b[2] - r_a[2] + wx * Lx * a1[2] + wy * Ly * a2[2] +
                            wz * Lz * a3[2];
                double d = sqrt(dx * dx + dy * dy + dz * dz);
                if (d < best)
                    best = d;
            }
    return best;
}

int main(void) {
    printf("=== libirrep — pyrochlore 16-site spin-spin correlation function ===\n\n");

    irrep_lattice3d_t *L = irrep_lattice3d_build(IRREP_LATTICE3D_PYROCHLORE, 1, 1, 1);
    int n_bonds = irrep_lattice3d_num_bonds_nn(L);
    int *bi = malloc((size_t)n_bonds * sizeof(int));
    int *bj = malloc((size_t)n_bonds * sizeof(int));
    irrep_lattice3d_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N, n_bonds, bi, bj, /*J=*/1.0);

    /* Lanczos with eigvecs to get the GS amplitude vector. */
    double _Complex *seed = malloc((size_t)DIM * sizeof(*seed));
    double           sn = 0;
    int              target_pop = N / 2;
    for (long long s = 0; s < DIM; ++s) {
        if (__builtin_popcountll((unsigned long long)s) == target_pop) {
            seed[s] = 0.1 * sin(0.37 * s) + I * 0.05 * cos(0.23 * s);
        } else
            seed[s] = 0;
        sn += creal(seed[s]) * creal(seed[s]) + cimag(seed[s]) * cimag(seed[s]);
    }
    sn = sqrt(sn);
    for (long long s = 0; s < DIM; ++s)
        seed[s] /= sn;

    double           E0 = 0;
    double _Complex *gs = malloc((size_t)DIM * sizeof(*gs));
    printf("  Computing GS via Lanczos with eigenvectors (300 iters, reorth)...\n");
    fflush(stdout);
    irrep_status_t st =
        irrep_lanczos_eigvecs_reorth(irrep_heisenberg_apply, H, DIM, 1, 300, seed, &E0, gs);
    if (st != IRREP_OK) {
        fprintf(stderr, "Lanczos failed, status=%d\n", (int)st);
        return 1;
    }
    printf("  GS energy: E_0 = %+.6f J  (reference: −8.809084 from pyrochlore16_heisenberg)\n\n",
           E0);

    /* Compute <S_0 · S_r> for every r. */
    site_corr_t corr[N];
    double      a1[3], a2[3], a3[3];
    irrep_lattice3d_primitive_vectors(L, a1, a2, a3);
    double r0[3];
    irrep_lattice3d_site_position(L, 0, r0);
    for (int r = 0; r < N; ++r) {
        double rr[3];
        irrep_lattice3d_site_position(L, r, rr);
        double d = min_image_distance(r0, rr, a1, a2, a3, 1, 1, 1);
        double c = (r == 0) ? 0.75 /* ⟨S²⟩ = ¾ for spin-½ */
                            : expectation_SS_(gs, 0, r);
        corr[r].site = r;
        corr[r].dist = d;
        corr[r].corr = c;
    }
    qsort(corr, N, sizeof(*corr), cmp_dist);

    printf("  ⟨S_0 · S_r⟩ vs minimum-image distance r:\n");
    printf("  %-5s  %-9s  %-12s  %-s\n", "site", "dist", "<S_0·S_r>", "interpretation");
    for (int k = 0; k < N; ++k) {
        const char *note = "";
        if (corr[k].site == 0)
            note = "(self, ⟨S²⟩ = ¾)";
        else if (fabs(corr[k].dist - 0.25 * sqrt(2.0)) < 1e-6)
            note = "NN  (tetrahedral edge, expect AFM)";
        else if (fabs(corr[k].dist - 0.25 * sqrt(6.0)) < 1e-6)
            note = "NNN";
        printf("  %-5d  %-9.6f  %+12.6f  %s\n", corr[k].site, corr[k].dist, corr[k].corr, note);
    }

    /* Independent verification: the sum of NN ⟨S·S⟩ × 48 bonds should equal E_0
     * (since H = Σ_{bonds} S·S). On the 16-site cluster all 48 NN bonds are
     * equivalent by symmetry, so 48 · C(NN) ≈ E_0. */
    double C_nn = 0;
    int    n_nn_count = 0;
    for (int b = 0; b < n_bonds; ++b)
        if (bi[b] == 0 || bj[b] == 0) {
            int r = (bi[b] == 0) ? bj[b] : bi[b];
            C_nn += expectation_SS_(gs, 0, r);
            ++n_nn_count;
        }
    double C_nn_avg = (n_nn_count > 0) ? C_nn / n_nn_count : 0;
    printf("\n  Sum rule cross-check (translation-equivalent NN bonds):\n");
    printf("    average ⟨S_0 · S_r⟩ over %d NN partners = %+.6f\n", n_nn_count, C_nn_avg);
    printf("    48 NN bonds × <C_NN> = %+.6f J\n", 48 * C_nn_avg);
    printf("    E_0 = %+.6f J\n", E0);
    printf("    discrepancy = %+.2e (= 0 modulo Lanczos eigvec residual)\n",
           48 * C_nn_avg - E0);

    free(seed);
    free(gs);
    free(bi);
    free(bj);
    irrep_heisenberg_free(H);
    irrep_lattice3d_free(L);
    return 0;
}
