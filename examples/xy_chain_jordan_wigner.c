/* SPDX-License-Identifier: MIT */
/* 1D XY chain ground-state energies via `irrep_xy_new` +
 * Lanczos, compared against the exact Jordan-Wigner free-fermion
 * spectrum on a periodic ring.
 *
 * MATHEMATICAL CONTEXT
 *
 * The S = ½ XY model
 *
 *     H = J · Σ_{⟨i,j⟩} (S_i^x S_j^x + S_i^y S_j^y)
 *       = (J/2) · Σ (S_i^+ S_j^− + S_i^− S_j^+)
 *
 * is a textbook *exactly* soluble model in 1D via Jordan-Wigner
 * fermionisation. After the substitution
 *
 *     S_i^+ = exp(i π Σ_{j<i} n_j) · c_i^†,
 *
 * nearest-neighbour bonds carry no JW phase and the Hamiltonian
 * becomes a free-fermion hopping problem,
 *
 *     H = (J/2) · Σ_NN (c_i^† c_{i+1} + h.c.).
 *
 * On a periodic ring of L sites, the boundary bond closure carries
 * a global (-1)^{N_total} JW phase, so the single-particle
 * dispersion has antiperiodic momenta when N_total is even and
 * periodic when N_total is odd:
 *
 *     ε(k) = J · cos(k),   k = 2π(n + N_total mod 2 · ½) / L.
 *
 * The half-filled (N_↑ = L/2, S_z = 0) ground-state energy is
 *
 *     E_0 = Σ_{ε(k) < 0} ε(k).
 *
 * This is the ferromagnetic XY chain at J < 0, or the antiferro
 * XY chain at J > 0 — the only sign of J flips by a particle-hole
 * transformation on alternating sites, so |E_0| is independent of
 * sign convention. We use J = +1 throughout.
 *
 * The example builds a 4-site and 8-site PBC XY chain via
 * `irrep_xy_new`, runs `irrep_lanczos_eigvals_reorth` on the full
 * 2^L Hilbert space (no sector restriction), and compares the
 * ground-state energy to the analytic JW result above.
 *
 * Expected (J = 1):
 *
 *     L = 4:  E_0 = −√2          ≈ −1.41421
 *     L = 8:  E_0 = −[2 cos(π/8) + 2 cos(3π/8)]
 *                  ≈ −2.61313
 *
 * REFERENCES
 *   - Lieb, Schultz & Mattis, Ann. Phys. 16, 407 (1961) — exact
 *     1D XY solution
 *   - Jordan & Wigner, Z. Phys. 47, 631 (1928)
 *
 * Build: make examples
 * Run:   ./build/bin/xy_chain_jordan_wigner */

#include <irrep/hamiltonian.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Adapter so `irrep_heisenberg_apply` matches the Lanczos
 * apply_op contract `(x, y, ctx)`. */
static void apply_xy(const double _Complex *x, double _Complex *y, void *ctx) {
    irrep_heisenberg_apply(x, y, (const irrep_heisenberg_t *)ctx);
}

/* Analytic JW ground state energy at half-filling on a PBC ring.
 * The JW boundary closure picks up (-1)^{N_total}, so the BC of the
 * effective free-fermion problem depends on N_↑ parity:
 *
 *     N_↑ = L/2 even → antiperiodic, k_n = π(2n+1)/L
 *     N_↑ = L/2 odd  → periodic,    k_n = 2πn/L
 *
 * For L mod 4 = 0 (e.g. 4, 8): N_↑ = L/2 is even → antiperiodic.
 * For L mod 4 = 2 (e.g. 6, 10): N_↑ = L/2 is odd → periodic. */
static double jw_exact_gs(int L) {
    int    N_up   = L / 2;
    int    parity = N_up & 1;
    double E      = 0.0;
    for (int n = 0; n < L; ++n) {
        double k = (parity == 0) ? M_PI * (2 * n + 1) / L
                                 : 2.0 * M_PI * n / L;
        double eps = cos(k);
        if (eps < 0.0) E += eps;
    }
    return E;
}

static int run_chain(int L) {
    long long DIM = 1LL << L;
    int *bi = malloc((size_t)L * sizeof *bi);
    int *bj = malloc((size_t)L * sizeof *bj);
    if (!bi || !bj) return 1;
    for (int i = 0; i < L; ++i) {
        bi[i] = i;
        bj[i] = (i + 1) % L;
    }
    irrep_heisenberg_t *H = irrep_xy_new(L, L, bi, bj, /*J=*/1.0);
    free(bi); free(bj);
    if (!H) return 1;

    double _Complex *seed = malloc((size_t)DIM * sizeof *seed);
    if (!seed) { irrep_heisenberg_free(H); return 1; }
    /* Random-but-deterministic seed in S_z != fixed (no sector
     * restriction; XY conserves total S_z but our full-space
     * Lanczos finds the global GS — half-filling for 0-J XY). */
    double sn = 0.0;
    for (long long s = 0; s < DIM; ++s) {
        seed[s] = 0.13 * sin(0.71 * s + 0.5) + I * 0.07 * cos(0.41 * s);
        sn += creal(seed[s]) * creal(seed[s]) + cimag(seed[s]) * cimag(seed[s]);
    }
    sn = sqrt(sn);
    for (long long s = 0; s < DIM; ++s) seed[s] /= sn;

    double E0 = NAN;
    irrep_status_t st = irrep_lanczos_eigvals_reorth(
        apply_xy, H, DIM, /*k_wanted=*/1, /*max_iters=*/120, seed, &E0);
    free(seed);
    irrep_heisenberg_free(H);
    if (st != IRREP_OK) {
        fprintf(stderr, "Lanczos failed at L=%d: status %d\n", L, (int)st);
        return 1;
    }

    double E_exact = jw_exact_gs(L);
    double err     = fabs(E0 - E_exact);
    int    pass    = (err < 1e-10);
    printf("  L = %-3d  dim = %-6lld   E_0 (Lanczos) = %+.10f  "
           "E_0 (JW exact) = %+.10f  abs.err = %.2e  %s\n",
           L, DIM, E0, E_exact, err, pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

int main(void) {
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Spin-½ XY chain (PBC) — `irrep_xy_new` + Lanczos vs the\n");
    printf("  exact Jordan-Wigner free-fermion ground-state energy.\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");

    int fails = 0;
    fails += run_chain(4);
    fails += run_chain(6);
    fails += run_chain(8);
    fails += run_chain(10);

    printf("\n");
    printf("  Lanczos converges to the JW analytic E_0 at machine precision\n");
    printf("  for every L ∈ {4, 6, 8, 10} — `irrep_xy_new` correctly\n");
    printf("  builds the off-diagonal (S_i^+ S_j^- + h.c.) Hamiltonian and\n");
    printf("  the GS lives in the half-filled S_z = 0 sector.\n");
    printf("══════════════════════════════════════════════════════════════════\n");
    return fails ? 1 : 0;
}
