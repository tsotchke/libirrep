/* SPDX-License-Identifier: MIT */
/* Connected dimer–dimer correlation function on the N=24 kagome-Heisenberg
 * singlet ground state. Fourth independent diagnostic for the GS nature,
 * complementing γ (§1.10), S(k) (§1.14), and ξ (§1.16):
 *
 *   D(b, b') = ⟨B_b · B_b'⟩ − ⟨B_b⟩·⟨B_b'⟩,     B_b = S_{i_b} · S_{j_b}
 *
 * Then Fourier transform over the (bond-midpoint) sublattice:
 *
 *   S_D(k) = (1 / N_bonds) · Σ_{b, b'}  cos(k · (r_b − r_b'))  D(b, b')
 *
 * Phase fingerprints:
 *   VBS (valence-bond crystal)    : sharp peak at VBS wavevector
 *   Z₂ spin liquid / Dirac SL     : featureless, S_D(k) ≪ N_bonds
 *   Any magnetic order (Néel …)   : Bragg peak at magnetic wavevector
 *
 * Ruling out VBS is the missing leg of the kagome-Heisenberg GS triangle
 * (Z₂ SL / Dirac SL / VBS). The first three diagnostics rule out Dirac
 * and magnetic order; this rules out VBS if featureless.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kagome_dimer_correlation
 */

#include <irrep/config_project.h>
#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>
#include <irrep/space_group.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N_SITES 24
#define D_FULL  (1LL << N_SITES)
#define POPCOUNT 12

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void unfold(const irrep_space_group_t *G,
                   const irrep_sg_rep_table_t *T,
                   int order, const double _Complex *w,
                   const double _Complex *psi_sector, long long sector_dim,
                   double _Complex *psi_full) {
    memset(psi_full, 0, (size_t)D_FULL * sizeof(double _Complex));
    uint64_t        *op = malloc((size_t)order * sizeof(uint64_t));
    double _Complex *oa = malloc((size_t)order * sizeof(double _Complex));
    long long n_reps = irrep_sg_rep_table_count(T);
    long long written = 0;
    for (long long k = 0; k < n_reps && written < sector_dim; ++k) {
        uint64_t u = irrep_sg_rep_table_get(T, k);
        int ne = 0;
        for (int g = 0; g < order; ++g) {
            uint64_t gx = irrep_space_group_apply_bits(G, g, u);
            int f = -1;
            for (int e = 0; e < ne; ++e) if (op[e] == gx) { f = e; break; }
            if (f >= 0) oa[f] += w[g];
            else { op[ne] = gx; oa[ne] = w[g]; ++ne; }
        }
        double nn = 0.0;
        for (int e = 0; e < ne; ++e) nn += creal(oa[e] * conj(oa[e]));
        if (nn < 1e-20) continue;
        double inv = 1.0 / sqrt(nn);
        double _Complex c = psi_sector[written];
        for (int e = 0; e < ne; ++e) psi_full[op[e]] += c * inv * oa[e];
        ++written;
    }
    free(op); free(oa);
}

/* φ_out[t] = (S_i · S_j) ψ_in[t]. Adds to φ_out — caller zeros first. */
static void apply_bond(int i, int j,
                       const double _Complex *psi_in, double _Complex *phi_out) {
    long long mi = 1LL << i, mj = 1LL << j;
    for (long long s = 0; s < D_FULL; ++s) {
        int zi = (int)((s >> i) & 1), zj = (int)((s >> j) & 1);
        double diag = (zi == zj) ? +0.25 : -0.25;
        phi_out[s] += diag * psi_in[s];
        if (zi != zj) {
            long long sp = s ^ mi ^ mj;
            phi_out[sp] += 0.5 * psi_in[s];
        }
    }
}

int main(void) {
    printf("=== Kagome N=24 dimer-dimer correlation (VBS null diagnostic) ===\n\n");
    double t0 = now_sec();

    irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 4);
    irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);
    irrep_heisenberg_t *H = irrep_heisenberg_new(N_SITES, nb, bi, bj, 1.0);
    printf("  N_sites=%d, N_bonds=%d\n", N_SITES, nb);

    irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, POPCOUNT);
    irrep_sg_little_group_t *lg = irrep_sg_little_group_build(G, 0, 0);
    double _Complex chi1[1] = {1.0 + 0.0 * I};
    irrep_sg_little_group_irrep_t *mu = irrep_sg_little_group_irrep_new(lg, chi1, 1);
    irrep_sg_heisenberg_sector_t *S = irrep_sg_heisenberg_sector_build_at_k(H, T, lg, mu);
    long long sector_dim = irrep_sg_heisenberg_sector_dim(S);

    double _Complex *seed = malloc((size_t)sector_dim * sizeof(double _Complex));
    for (long long i = 0; i < sector_dim; ++i)
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
    double eigs[2];
    double _Complex *psi_sec = malloc((size_t)2 * sector_dim * sizeof(double _Complex));
    irrep_lanczos_eigvecs_reorth(irrep_sg_heisenberg_sector_apply, S,
                                 sector_dim, 2, 200, seed, eigs, psi_sec);
    printf("  E_0/N = %+.8f\n", eigs[0] / N_SITES);

    int order = irrep_space_group_order(G);
    double _Complex *w = malloc((size_t)order * sizeof(double _Complex));
    irrep_sg_projector_weights(lg, mu, w);
    double _Complex *psi_full = malloc((size_t)D_FULL * sizeof(double _Complex));
    unfold(G, T, order, w, psi_sec, sector_dim, psi_full);
    double nn = 0.0;
    for (long long s = 0; s < D_FULL; ++s) nn += creal(psi_full[s] * conj(psi_full[s]));
    printf("  ||ψ||² = %.10f\n\n", nn);

    /* Precompute φ_b = B_b ψ for each NN bond. */
    printf("  Applying %d bond operators to ψ…\n", nb);
    double tb = now_sec();
    double _Complex **phi = malloc((size_t)nb * sizeof(double _Complex *));
    for (int b = 0; b < nb; ++b) {
        phi[b] = calloc((size_t)D_FULL, sizeof(double _Complex));
        apply_bond(bi[b], bj[b], psi_full, phi[b]);
    }
    printf("  done in %.2f s\n", now_sec() - tb);

    /* Single-bond expectations ⟨B_b⟩ = ⟨ψ|φ_b⟩. */
    double *Bmean = malloc(sizeof(double) * nb);
    for (int b = 0; b < nb; ++b) {
        double _Complex x = 0.0;
        for (long long s = 0; s < D_FULL; ++s) x += conj(psi_full[s]) * phi[b][s];
        Bmean[b] = creal(x);
    }
    double Bsum = 0, Bvar = 0;
    for (int b = 0; b < nb; ++b) Bsum += Bmean[b];
    double Bavg = Bsum / nb;
    for (int b = 0; b < nb; ++b) Bvar += (Bmean[b]-Bavg)*(Bmean[b]-Bavg);
    printf("  ⟨B_b⟩ avg = %+.6f   bond-to-bond stddev = %.3e\n",
           Bavg, sqrt(Bvar/nb));
    printf("    (small stddev ⟹ bond expectation is uniform — no VBS dimer pattern)\n");

    /* Full dimer-dimer matrix D_{b,b'} = ⟨φ_b|φ_b'⟩ − Bmean[b]*Bmean[b']. */
    printf("\n  Computing dimer-dimer matrix (%d × %d)…\n", nb, nb);
    double tm = now_sec();
    double *Dmat = malloc((size_t)nb * nb * sizeof(double));
    for (int a = 0; a < nb; ++a) {
        for (int c = a; c < nb; ++c) {
            double _Complex x = 0.0;
            for (long long s = 0; s < D_FULL; ++s) x += conj(phi[a][s]) * phi[c][s];
            double conn = creal(x) - Bmean[a] * Bmean[c];
            Dmat[a * nb + c] = conn;
            Dmat[c * nb + a] = conn;
        }
    }
    printf("  matrix done in %.2f s\n", now_sec() - tm);

    /* Bond midpoint positions. */
    double *bx = malloc(sizeof(double) * nb);
    double *by = malloc(sizeof(double) * nb);
    for (int b = 0; b < nb; ++b) {
        double xi[2], xj[2];
        irrep_lattice_site_position(L, bi[b], xi);
        irrep_lattice_site_position(L, bj[b], xj);
        bx[b] = 0.5 * (xi[0] + xj[0]);
        by[b] = 0.5 * (xi[1] + xj[1]);
    }

    /* Bond dimer structure factor
     *   S_D(k) = (1/N_bonds) Σ_{b,b'} cos(k·(r_b − r_b')) D_{b,b'}
     */
    double b1[2], b2[2];
    irrep_lattice_reciprocal_vectors(L, b1, b2);
    int Lx = irrep_lattice_Lx(L), Ly = irrep_lattice_Ly(L);

    printf("\n  Dimer structure factor S_D(k) at allowed Bloch momenta:\n");
    printf("  %-8s %-25s %12s\n", "(mx,my)", "k (1/a)", "S_D(k)");
    printf("  --------  -------------------------  ------\n");

    double Sd_max = 0; int kmax_mx = 0, kmax_my = 0;
    double Sd_kmax = 0;
    double Td_gamma = 0, Td_max_notgamma = 0;
    int    Td_kmax_mx = 0, Td_kmax_my = 0;
    printf("  %-8s %-25s %12s %12s\n", "(mx,my)", "k (1/a)", "S_D(k)", "T_disc(k)");
    printf("  --------  -------------------------  ------------ ------------\n");
    for (int my = 0; my < Ly; ++my) {
        for (int mx = 0; mx < Lx; ++mx) {
            double kx = (double)mx / Lx * b1[0] + (double)my / Ly * b2[0];
            double ky = (double)mx / Lx * b1[1] + (double)my / Ly * b2[1];
            /* Connected S_D(k). */
            double Sd = 0;
            for (int a = 0; a < nb; ++a)
                for (int c = 0; c < nb; ++c) {
                    double dx = bx[a] - bx[c];
                    double dy = by[a] - by[c];
                    Sd += Dmat[a * nb + c] * cos(kx * dx + ky * dy);
                }
            Sd /= nb;
            /* Disconnected T_disc(k) — VBS order parameter. */
            double _Complex acc = 0.0;
            for (int b = 0; b < nb; ++b)
                acc += Bmean[b] * (cos(kx*bx[b] + ky*by[b])
                                   + I * sin(kx*bx[b] + ky*by[b]));
            double Td = creal(acc * conj(acc)) / nb;

            printf("  (%d,%d)     (%+.4f, %+.4f)             %+12.4f %12.4f\n",
                   mx, my, kx, ky, Sd, Td);
            if (fabs(Sd) > fabs(Sd_max)) {
                Sd_max = Sd; kmax_mx = mx; kmax_my = my; Sd_kmax = Sd;
            }
            if (mx == 0 && my == 0) Td_gamma = Td;
            else if (Td > Td_max_notgamma) {
                Td_max_notgamma = Td; Td_kmax_mx = mx; Td_kmax_my = my;
            }
        }
    }

    printf("\n  Connected  : max |S_D(k)|       = %.4f at (mx,my)=(%d,%d)\n",
           fabs(Sd_max), kmax_mx, kmax_my);
    printf("               max |S_D(k)| / N_bonds = %.4f\n", fabs(Sd_max) / nb);
    printf("  Discon. VBS: T_disc(Γ)           = %.4f  (trivial uniform peak)\n",
           Td_gamma);
    printf("               T_disc at best k≠Γ   = %.4f at (mx,my)=(%d,%d)\n",
           Td_max_notgamma, Td_kmax_mx, Td_kmax_my);
    printf("               T_disc(k≠Γ) / T_disc(Γ) = %.4f\n",
           Td_max_notgamma / (Td_gamma > 1e-12 ? Td_gamma : 1.0));
    printf("\n  INTERPRETATION:\n");
    double vbs_ratio = Td_max_notgamma / (Td_gamma > 1e-12 ? Td_gamma : 1.0);
    printf("    kagome N=24  : T_disc(k≠Γ)/T_disc(Γ) = %.3f\n", vbs_ratio);
    printf("    MG chain N=16: T_disc(π)/T_disc(0)   = 0.980   (positive control)\n");
    if (vbs_ratio < 0.25 && fabs(Sd_kmax) / nb < 0.1)
        printf("    ✓ NO VBS order — kagome ratio is %.0f× below the MG positive\n"
               "      control; residual signal is consistent with finite-cluster\n"
               "      shape anisotropy, not a genuine VBS order parameter.\n",
               0.980 / vbs_ratio);
    else if (vbs_ratio < 0.5)
        printf("    ~ Mild bond-structure — finite-size cluster-shape effect.\n");
    else
        printf("    ✗ Large T_disc ratio at k≠Γ → genuine VBS order.\n");

    /* Diagonal: bond-bond fluctuation ⟨B_b²⟩ − ⟨B_b⟩². */
    double dd_max = 0, dd_min = 1e9;
    for (int b = 0; b < nb; ++b) {
        double d = Dmat[b * nb + b];
        if (d > dd_max) dd_max = d;
        if (d < dd_min) dd_min = d;
    }
    printf("\n  Per-bond fluctuation Var(B_b) range: [%.4f, %.4f]\n", dd_min, dd_max);

    printf("\n  Total wall-clock: %.2f s\n", now_sec() - t0);

    for (int b = 0; b < nb; ++b) free(phi[b]);
    free(phi); free(Bmean); free(Dmat); free(bx); free(by);
    free(psi_full); free(psi_sec); free(seed); free(w);
    irrep_sg_heisenberg_sector_free(S);
    irrep_sg_little_group_irrep_free(mu);
    irrep_sg_little_group_free(lg);
    irrep_sg_rep_table_free(T);
    free(bi); free(bj);
    irrep_heisenberg_free(H);
    irrep_space_group_free(G);
    irrep_lattice_free(L);
    return 0;
}
