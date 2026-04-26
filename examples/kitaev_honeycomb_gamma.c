/* SPDX-License-Identifier: MIT */
/* Positive-control calibration for γ-extraction methodology on an
 * exactly-solvable topological ground state.
 *
 * The Kitaev honeycomb model,
 *
 *   H = − Σ_xx K_x σ^x_i σ^x_j − Σ_yy K_y σ^y_i σ^y_j − Σ_zz K_z σ^z_i σ^z_j,
 *
 * is exactly solvable via Majorana factorization (Kitaev 2006). In the
 * anisotropic A-phase (K_z ≫ K_x, K_y) the effective low-energy theory
 * is the toric code with EXACT γ = log 2 at any finite N. This is the
 * positive control we need: if our area-law / KP γ extraction
 * methodology cannot find γ = log 2 on this state, the methodology is
 * finite-N-limited regardless of system. If it DOES find it, the
 * kagome result earns credit for that methodology.
 *
 * Cluster: 3×4 honeycomb torus (N = 24, same size as kagome 2×4).
 * Bond-type labeling: x-bond has Δy ≈ 0; y-bond has Δy > 0; z-bond
 * has Δy < 0. Anisotropy: K_z = 1, K_x = K_y = 0.1 (deep A-phase).
 *
 * Diagnostics: E₀, ⟨H²⟩ residual, Lanczos convergence, then area-law
 * γ on 3 region families and KP γ on 3 region geometries.
 *
 *   make USE_OPENMP=1 examples
 *   ./build/bin/kitaev_honeycomb_gamma
 */

#include <irrep/hamiltonian.h>
#include <irrep/lattice.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define LX 3
#define LY 4
#define N_SITES (2 * LX * LY)          /* 2 sites / unit cell */
#define D_FULL  (1LL << N_SITES)

static double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Kitaev context: bond list with per-bond type (0 = x, 1 = y, 2 = z)
 * and 3-component coupling vector {K_x, K_y, K_z}. */
typedef struct {
    int     nb;
    int    *bi, *bj, *bt;
    double  Kx, Ky, Kz;
} kitaev_ctx_t;

/* Kitaev matrix-vector: out = H · psi.
 *
 *   σ^z σ^z: diagonal, +1 if bits aligned, −1 if anti-aligned
 *   σ^x σ^x: flips both bits, always coefficient +1
 *   σ^y σ^y: flips both bits, coefficient +1 if anti-aligned, −1 if aligned
 *
 * H has a leading minus sign, so contributions are negated. */
static void kitaev_apply(const double _Complex *psi,
                         double _Complex *out, void *opaque) {
    kitaev_ctx_t *K = (kitaev_ctx_t *)opaque;
    memset(out, 0, (size_t)D_FULL * sizeof(double _Complex));
    for (int b = 0; b < K->nb; ++b) {
        int i = K->bi[b], j = K->bj[b], t = K->bt[b];
        long long mi = 1LL << i, mj = 1LL << j;
        double coup = (t == 0) ? K->Kx : (t == 1) ? K->Ky : K->Kz;
        /* H has -coup * Σ σ^α_i σ^α_j → apply with minus. */
        for (long long s = 0; s < D_FULL; ++s) {
            int zi = (int)((s >> i) & 1), zj = (int)((s >> j) & 1);
            if (t == 2) {
                /* σ^z σ^z: diagonal */
                double sgn = (zi == zj) ? +1.0 : -1.0;
                out[s] += -coup * sgn * psi[s];
            } else if (t == 0) {
                /* σ^x σ^x: flip both bits, coefficient +1 */
                long long sp = s ^ mi ^ mj;
                out[sp] += -coup * psi[s];
            } else {
                /* σ^y σ^y: flip both bits, +1 if anti-aligned, −1 if aligned */
                long long sp = s ^ mi ^ mj;
                double sgn = (zi == zj) ? -1.0 : +1.0;
                out[sp] += -coup * sgn * psi[s];
            }
        }
    }
}

static int classify_bond(double dx, double dy) {
    if (fabs(dy) < 0.1)          return 0;  /* horizontal → x */
    if (dy > 0)                  return 1;  /* upward slant → y */
    return 2;                                /* downward slant → z */
}

static int boundary_bonds(const int *A, int nA, const int *bi, const int *bj, int nb) {
    char in_A[64] = {0};
    for (int k = 0; k < nA; ++k) in_A[A[k]] = 1;
    int boundary = 0;
    for (int b = 0; b < nb; ++b)
        if (in_A[bi[b]] != in_A[bj[b]]) ++boundary;
    return boundary;
}

static double entropy_of_region(const double _Complex *psi, const int *X, int nX) {
    long long dX = 1LL << nX;
    double _Complex *rho = malloc((size_t)dX * dX * sizeof(double _Complex));
    irrep_partial_trace(N_SITES, 2, psi, X, nX, rho);
    double S = irrep_entropy_vonneumann(rho, (int)dX);
    free(rho);
    return S;
}

int main(void) {
    printf("=== Kitaev honeycomb A-phase γ POSITIVE CONTROL ===\n");
    printf("    3×4 torus (N=24), K_z=1, K_x=K_y=0.1\n");
    printf("    Exact prediction: γ = log 2 = %.4f (toric-code A-phase)\n\n", log(2.0));
    double t0 = now_sec();

    irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_HONEYCOMB, LX, LY);
    int nb = irrep_lattice_num_bonds_nn(L);
    int *bi = malloc(sizeof(int) * nb);
    int *bj = malloc(sizeof(int) * nb);
    int *bt = malloc(sizeof(int) * nb);
    irrep_lattice_fill_bonds_nn(L, bi, bj);

    /* Classify bonds by geometric direction. */
    double a1[2], a2[2];
    irrep_lattice_primitive_vectors(L, a1, a2);
    int nx_type[3] = {0, 0, 0};
    for (int b = 0; b < nb; ++b) {
        double pi_[2], pj_[2];
        irrep_lattice_site_position(L, bi[b], pi_);
        irrep_lattice_site_position(L, bj[b], pj_);
        double dx0 = pj_[0] - pi_[0], dy0 = pj_[1] - pi_[1];
        double best_r = 1e300, bdx = 0, bdy = 0;
        for (int nx = -1; nx <= 1; ++nx)
            for (int ny = -1; ny <= 1; ++ny) {
                double dx = dx0 + nx * LX * a1[0] + ny * LY * a2[0];
                double dy = dy0 + nx * LX * a1[1] + ny * LY * a2[1];
                double r = sqrt(dx*dx + dy*dy);
                if (r < best_r) { best_r = r; bdx = dx; bdy = dy; }
            }
        bt[b] = classify_bond(bdx, bdy);
        nx_type[bt[b]]++;
    }
    printf("  bonds: %d total (x=%d, y=%d, z=%d)\n",
           nb, nx_type[0], nx_type[1], nx_type[2]);

    kitaev_ctx_t K = { nb, bi, bj, bt, 0.1, 0.1, 1.0 };

    /* Seed + Lanczos on full 2^N. */
    double _Complex *seed = malloc((size_t)D_FULL * sizeof(double _Complex));
    for (long long s = 0; s < D_FULL; ++s)
        seed[s] = 0.1 * sin(0.37 * s) + I * 0.05 * cos(0.23 * s);
    double sn = 0;
    for (long long s = 0; s < D_FULL; ++s) sn += creal(seed[s] * conj(seed[s]));
    sn = sqrt(sn);
    for (long long s = 0; s < D_FULL; ++s) seed[s] /= sn;

    double eigs[4];
    double _Complex *psi = malloc((size_t)4 * D_FULL * sizeof(double _Complex));
    double tL = now_sec();
    irrep_lanczos_eigvecs_reorth(kitaev_apply, &K, D_FULL, 4, 80, seed, eigs, psi);
    printf("  Lanczos (200 iters, dim=%lld): %.1f s\n", D_FULL, now_sec() - tL);
    printf("  E_0 = %+.8f   E_1 = %+.8f   ΔE₀₁ = %+.6e\n",
           eigs[0], eigs[1], eigs[1] - eigs[0]);
    printf("  E_2 = %+.8f   E_3 = %+.8f\n", eigs[2], eigs[3]);
    fflush(stdout);

    /* Residual. */
    double _Complex *Hpsi = malloc((size_t)D_FULL * sizeof(double _Complex));
    kitaev_apply(psi, Hpsi, &K);
    double res2 = 0;
    for (long long i = 0; i < D_FULL; ++i) {
        double _Complex r = Hpsi[i] - eigs[0] * psi[i];
        res2 += creal(r * conj(r));
    }
    printf("  ‖H·ψ − E·ψ‖ = %.3e\n\n", sqrt(res2));
    fflush(stdout);
    free(Hpsi); free(seed);

    /* ═══ Area-law γ across 3 region families ═══ */
    printf("━ Area-law γ across 3 region families (nested [0, |A|)) ━\n");
    const int fam[3][12] = {
        {0,1,2,3,4,5,6,7,8,9,10,11},
        {0,2,4,6,8,10,12,14,16,18,20,22},
        {0,3,6,9,12,15,18,21, 1,4,7,10},
    };
    const char *fam_name[3] = {"contiguous    ", "stride-2      ", "mixed-sublatt "};
    for (int f = 0; f < 3; ++f) {
        double xs[12], ys[12];
        int nA_max = 9;
        for (int nA = 1; nA <= nA_max; ++nA) {
            int A[16];
            for (int i = 0; i < nA; ++i) A[i] = fam[f][i];
            int bdry = boundary_bonds(A, nA, bi, bj, nb);
            double S_A = entropy_of_region(psi, A, nA);
            xs[nA - 1] = bdry; ys[nA - 1] = S_A;
        }
        int npts = nA_max;
        double mxv=0, myv=0, mxxv=0, mxyv=0;
        for (int i = 0; i < npts; ++i) {
            mxv+=xs[i]; myv+=ys[i]; mxxv+=xs[i]*xs[i]; mxyv+=xs[i]*ys[i];
        }
        mxv/=npts; myv/=npts; mxxv/=npts; mxyv/=npts;
        double alpha = (mxyv - mxv*myv)/(mxxv - mxv*mxv);
        double beta  = myv - alpha*mxv;
        double ss_res=0, ss_tot=0;
        for (int i = 0; i < npts; ++i) {
            double pred = alpha*xs[i] + beta;
            ss_res += (ys[i]-pred)*(ys[i]-pred);
            ss_tot += (ys[i]-myv)*(ys[i]-myv);
        }
        double R2 = 1.0 - ss_res/ss_tot;
        printf("  %s α=%+.4f β=%+.4f γ=%+.4f (R²=%.4f)\n",
               fam_name[f], alpha, beta, -beta, R2);
        fflush(stdout);
    }

    /* ═══ KP γ across 3 geometries (3-site regions → |ABC|=9, fast) ═══ */
    printf("\n━ KP γ (S_A+S_B+S_C−S_AB−S_AC−S_BC+S_ABC) across 3 geometries ━\n");
    int geoms[3][3][8] = {
        {{0,1,2}, {3,4,5}, {6,7,8}},
        {{0,2,4}, {6,8,10}, {12,14,16}},
        {{0,3,6}, {1,4,7}, {2,5,8}},
    };
    const char *geom_name[3] = {
        "compact strips",
        "separated strips (stride-2)",
        "sublattice mix",
    };
    int reg_size = 3;
    double gammas[3];
    for (int g_i = 0; g_i < 3; ++g_i) {
        double tg = now_sec();
        int *A = geoms[g_i][0], *B = geoms[g_i][1], *C = geoms[g_i][2];
        int AB[8], AC[8], BC[8], ABC[12];
        int nAB=0, nAC=0, nBC=0, nABC=0;
        for (int i = 0; i < reg_size; ++i) { AB[nAB++]=A[i]; AC[nAC++]=A[i]; ABC[nABC++]=A[i]; }
        for (int i = 0; i < reg_size; ++i) { AB[nAB++]=B[i]; BC[nBC++]=B[i]; ABC[nABC++]=B[i]; }
        for (int i = 0; i < reg_size; ++i) { AC[nAC++]=C[i]; BC[nBC++]=C[i]; ABC[nABC++]=C[i]; }
        double S_A  = entropy_of_region(psi, A, reg_size);
        double S_B  = entropy_of_region(psi, B, reg_size);
        double S_C  = entropy_of_region(psi, C, reg_size);
        double S_AB = entropy_of_region(psi, AB, nAB);
        double S_AC = entropy_of_region(psi, AC, nAC);
        double S_BC = entropy_of_region(psi, BC, nBC);
        double S_ABC = entropy_of_region(psi, ABC, nABC);
        double g = S_A + S_B + S_C - S_AB - S_AC - S_BC + S_ABC;
        gammas[g_i] = g;
        printf("  %-35s γ_KP = %+.4f  (%.1fs)\n",
               geom_name[g_i], g, now_sec() - tg);
        fflush(stdout);
    }
    double g_mean = (gammas[0]+gammas[1]+gammas[2]) / 3.0;
    double g_min = gammas[0], g_max = gammas[0];
    for (int i = 0; i < 3; ++i) {
        if (gammas[i] < g_min) g_min = gammas[i];
        if (gammas[i] > g_max) g_max = gammas[i];
    }
    printf("\n  KP γ: mean = %+.4f  range = [%+.4f, %+.4f]  spread = %.4f\n",
           g_mean, g_min, g_max, g_max - g_min);
    printf("  EXACT prediction for A-phase Kitaev: γ = log 2 = %+.4f\n", log(2.0));
    printf("  Deviation (mean − log 2) = %+.4f (%.0f%%)\n",
           g_mean - log(2.0), 100.0 * (g_mean - log(2.0)) / log(2.0));
    if (fabs(g_mean - log(2.0)) < 0.15 && (g_max - g_min) < 0.3)
        printf("  ✓ POSITIVE CONTROL PASSED: methodology finds γ ≈ log 2 on Kitaev.\n"
               "    Kagome KP γ can be trusted to identify topological phase.\n");
    else
        printf("  ✗ POSITIVE CONTROL FAILED: KP γ cannot robustly identify log 2\n"
               "    on exactly-solvable Kitaev at this cluster size. Need bigger N\n"
               "    or different diagnostic.\n");

    free(psi); free(bi); free(bj); free(bt);
    irrep_lattice_free(L);

    printf("\n━━━ Total wall-clock: %.1f s ━━━\n", now_sec() - t0);
    return 0;
}
