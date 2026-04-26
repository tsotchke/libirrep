/* SPDX-License-Identifier: MIT */
/* Momentum-resolved Heisenberg ED on a 3D Bravais lattice.
 *
 * Generalises `lattice3d_sector_ed.c` from the Γ-momentum sector to every
 * allowed k = (kx, ky, kz) on the cluster's Brillouin-zone grid. The
 * ground-state energy E₀(k) per momentum sector is the standard tower-
 * of-states diagnostic: degeneracy of E₀ across multiple k sectors
 * signals symmetry-broken phases (e.g. Anderson tower for Néel / valence-
 * bond crystals); a clean E₀(Γ) gap to all k ≠ Γ signals a featureless
 * symmetric phase.
 *
 * Sector basis (general k, 1D irreps of the translation group):
 *
 *   |k, u⟩ = (1/√(|G|·σ_u)) Σ_g e^{-i k·t_g} |T_g · u⟩
 *
 * where Stab_u = {g : T_g · u = u} and σ_u = Σ_{g ∈ Stab_u} e^{-i k·t_g}.
 * Orbits with σ_u = 0 are annihilated by the projector and filtered out.
 *
 * Off-diagonal Heisenberg matrix element (for anti-aligned bond (a,b) in
 * canonical u, flipped to canonical(u_ab) = v with canonicalising
 * translation t_R, i.e. T_{t_R} · u_ab = v):
 *
 *   ⟨k, v | H | k, u⟩ += ½ J · e^{-i k·t_{t_R}} · √(σ_v / σ_u)
 *
 * which reduces to the Γ formula `½ J · √(N_u/N_v)` (since σ_u = |Stab_u|
 * = |G|/N_u and the phase factor is 1) when k = 0.
 *
 * Cluster: SC 2×2×2 (8 sites, dim 256). At this size we can verify the
 * momentum-resolved spectrum against full ED — the minimum E₀(k) across
 * all 8 k sectors must equal the full-Hilbert-space ground state. */

#include <irrep/lattice3d.h>
#include <irrep/rdm.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* -------------------------------------------------------------------------- *
 * Sector context (per-cluster, momentum-independent)                         *
 * -------------------------------------------------------------------------- */

typedef struct {
    int      N;
    int      n_bonds;
    int     *bi, *bj;
    double   J;

    /* Translation table — same for every k. */
    int      n_trans;
    int    **trans_perm;     /* [n_trans][N] */
    int     *trans_dx;       /* [n_trans] integer cell offsets in the conventional cell */
    int     *trans_dy;
    int     *trans_dz;

    double a1[3], a2[3], a3[3];
    double b1[3], b2[3], b3[3];

    /* Canonical orbit table — momentum-independent. */
    int       n_reps;
    uint64_t *reps;          /* canonical rep per orbit, sorted ascending */
    int      *orbit_size;    /* |orb_u| per orbit */
    int      *stab_size;     /* |Stab_u| = n_trans / orbit_size */

    /* Per-orbit stabiliser membership: stab_g[i*n_trans + g] = 1 if T_g
     * stabilises canonical rep i. Used to compute σ_u for each k. */
    char     *stab_g;
} cluster_t;

/* k-dependent sector data: σ_u, filtered rep list, complex apply. */
typedef struct {
    const cluster_t *C;
    double           k[3]; /* cartesian momentum */

    int              dim;            /* number of non-filtered orbits at this k */
    int             *rep_index;      /* [dim] → cluster's canonical rep index */
    double _Complex *sigma;          /* [dim] σ_u value */
    double          *sigma_abs;      /* [dim] |σ_u|, precomputed */

    /* Mapping: cluster rep index → dim-index, or -1 if filtered out. */
    int *full_to_dim;
} ksector_t;

/* -------------------------------------------------------------------------- *
 * Bit-permutation helper                                                     *
 * -------------------------------------------------------------------------- */

static uint64_t apply_perm_bits(uint64_t s, const int *perm, int N) {
    uint64_t r = 0;
    for (int i = 0; i < N; ++i)
        if ((s >> i) & 1ULL)
            r |= 1ULL << perm[i];
    return r;
}

/* Canonical rep of `s` and the index of the translation that achieves it. */
static uint64_t canonicalise2(uint64_t s, const cluster_t *C, int *out_t) {
    uint64_t best = s;
    int      best_t = 0;
    for (int t = 0; t < C->n_trans; ++t) {
        uint64_t r = apply_perm_bits(s, C->trans_perm[t], C->N);
        if (r < best) {
            best = r;
            best_t = t;
        }
    }
    *out_t = best_t;
    return best;
}

static int cmp_u64(const void *a, const void *b) {
    uint64_t ua = *(const uint64_t *)a, ub = *(const uint64_t *)b;
    return ua < ub ? -1 : ua > ub ? 1 : 0;
}

static int rep_index(uint64_t v, const cluster_t *C) {
    int lo = 0, hi = C->n_reps - 1;
    while (lo <= hi) {
        int      m = (lo + hi) / 2;
        uint64_t mv = C->reps[m];
        if (mv == v)
            return m;
        if (mv < v)
            lo = m + 1;
        else
            hi = m - 1;
    }
    return -1;
}

/* -------------------------------------------------------------------------- *
 * Cluster build (Sz=0 sector, all canonical reps)                            *
 * -------------------------------------------------------------------------- */

static void build_cluster_(cluster_t *C, irrep_lattice3d_kind_t kind, int Lx, int Ly, int Lz,
                           int Sz_target, double J) {
    irrep_lattice3d_t *L = irrep_lattice3d_build(kind, Lx, Ly, Lz);
    C->N = irrep_lattice3d_num_sites(L);
    C->n_bonds = irrep_lattice3d_num_bonds_nn(L);
    C->bi = malloc((size_t)C->n_bonds * sizeof(int));
    C->bj = malloc((size_t)C->n_bonds * sizeof(int));
    irrep_lattice3d_fill_bonds_nn(L, C->bi, C->bj);
    C->J = J;

    irrep_lattice3d_primitive_vectors(L, C->a1, C->a2, C->a3);
    irrep_lattice3d_reciprocal_vectors(L, C->b1, C->b2, C->b3);

    C->n_trans = Lx * Ly * Lz;
    C->trans_perm = malloc((size_t)C->n_trans * sizeof(*C->trans_perm));
    C->trans_dx = malloc((size_t)C->n_trans * sizeof(int));
    C->trans_dy = malloc((size_t)C->n_trans * sizeof(int));
    C->trans_dz = malloc((size_t)C->n_trans * sizeof(int));
    int t = 0;
    for (int dz = 0; dz < Lz; ++dz)
        for (int dy = 0; dy < Ly; ++dy)
            for (int dx = 0; dx < Lx; ++dx) {
                C->trans_perm[t] = malloc((size_t)C->N * sizeof(int));
                for (int s = 0; s < C->N; ++s)
                    C->trans_perm[t][s] = irrep_lattice3d_translate(L, s, dx, dy, dz);
                C->trans_dx[t] = dx;
                C->trans_dy[t] = dy;
                C->trans_dz[t] = dz;
                ++t;
            }

    /* Enumerate Sz-target canonical reps. */
    int       target_pop = C->N / 2 + Sz_target;
    long long total_dim = 1LL << C->N;
    uint64_t *buf = malloc((size_t)total_dim * sizeof(uint64_t));
    int       n_can = 0;
    for (long long s = 0; s < total_dim; ++s) {
        if (__builtin_popcountll((unsigned long long)s) != target_pop)
            continue;
        int      t_dummy;
        uint64_t c = canonicalise2((uint64_t)s, C, &t_dummy);
        if ((uint64_t)s == c)
            buf[n_can++] = c;
    }
    qsort(buf, (size_t)n_can, sizeof(uint64_t), cmp_u64);
    C->n_reps = n_can;
    C->reps = malloc((size_t)n_can * sizeof(uint64_t));
    memcpy(C->reps, buf, (size_t)n_can * sizeof(uint64_t));
    free(buf);

    /* Stabiliser membership and orbit / stab sizes. */
    C->orbit_size = malloc((size_t)n_can * sizeof(int));
    C->stab_size = malloc((size_t)n_can * sizeof(int));
    C->stab_g = calloc((size_t)n_can * (size_t)C->n_trans, sizeof(char));
    for (int i = 0; i < n_can; ++i) {
        uint64_t u = C->reps[i];
        int      n_stab = 0;
        for (int g = 0; g < C->n_trans; ++g) {
            uint64_t r = apply_perm_bits(u, C->trans_perm[g], C->N);
            if (r == u) {
                C->stab_g[i * C->n_trans + g] = 1;
                ++n_stab;
            }
        }
        C->stab_size[i] = n_stab;
        C->orbit_size[i] = C->n_trans / n_stab;
    }

    irrep_lattice3d_free(L);
}

static void free_cluster_(cluster_t *C) {
    free(C->bi);
    free(C->bj);
    for (int t = 0; t < C->n_trans; ++t)
        free(C->trans_perm[t]);
    free(C->trans_perm);
    free(C->trans_dx);
    free(C->trans_dy);
    free(C->trans_dz);
    free(C->reps);
    free(C->orbit_size);
    free(C->stab_size);
    free(C->stab_g);
}

/* -------------------------------------------------------------------------- *
 * Per-momentum sector build                                                  *
 * -------------------------------------------------------------------------- */

/* Cartesian translation t_g = dx · a1 + dy · a2 + dz · a3. */
static double k_dot_t_(const cluster_t *C, const double k[3], int g) {
    double tx = C->trans_dx[g] * C->a1[0] + C->trans_dy[g] * C->a2[0] +
                C->trans_dz[g] * C->a3[0];
    double ty = C->trans_dx[g] * C->a1[1] + C->trans_dy[g] * C->a2[1] +
                C->trans_dz[g] * C->a3[1];
    double tz = C->trans_dx[g] * C->a1[2] + C->trans_dy[g] * C->a2[2] +
                C->trans_dz[g] * C->a3[2];
    return k[0] * tx + k[1] * ty + k[2] * tz;
}

static void build_ksector_(ksector_t *K, const cluster_t *C, const double k[3]) {
    K->C = C;
    memcpy(K->k, k, 3 * sizeof(double));

    K->full_to_dim = malloc((size_t)C->n_reps * sizeof(int));
    K->rep_index = malloc((size_t)C->n_reps * sizeof(int));
    K->sigma = malloc((size_t)C->n_reps * sizeof(double _Complex));
    K->sigma_abs = malloc((size_t)C->n_reps * sizeof(double));

    int dim = 0;
    for (int i = 0; i < C->n_reps; ++i) {
        double _Complex sig = 0;
        for (int g = 0; g < C->n_trans; ++g) {
            if (C->stab_g[i * C->n_trans + g]) {
                double phase = -k_dot_t_(C, k, g);
                sig += cexp(I * phase);
            }
        }
        if (cabs(sig) < 1e-9) {
            K->full_to_dim[i] = -1; /* filtered out */
        } else {
            K->full_to_dim[i] = dim;
            K->rep_index[dim] = i;
            K->sigma[dim] = sig;
            K->sigma_abs[dim] = cabs(sig);
            ++dim;
        }
    }
    K->dim = dim;
}

static void free_ksector_(ksector_t *K) {
    free(K->full_to_dim);
    free(K->rep_index);
    free(K->sigma);
    free(K->sigma_abs);
}

/* -------------------------------------------------------------------------- *
 * Heisenberg apply at fixed k                                                *
 * -------------------------------------------------------------------------- */

static void H_apply_k(const double _Complex *x, double _Complex *y, void *ctx) {
    ksector_t       *K = ctx;
    const cluster_t *C = K->C;
    memset(y, 0, (size_t)K->dim * sizeof(double _Complex));
    for (int i_dim = 0; i_dim < K->dim; ++i_dim) {
        int      i = K->rep_index[i_dim];
        uint64_t u = C->reps[i];
        if (x[i_dim] == 0.0)
            continue;
        double diag = 0.0;
        for (int b = 0; b < C->n_bonds; ++b) {
            int a = C->bi[b], q = C->bj[b];
            int sa = (int)((u >> a) & 1ULL);
            int sq = (int)((u >> q) & 1ULL);
            if (sa == sq) {
                diag += 0.25 * C->J;
            } else {
                diag -= 0.25 * C->J;
                uint64_t up = u ^ (1ULL << a) ^ (1ULL << q);
                int      t_R;
                uint64_t v = canonicalise2(up, C, &t_R);
                int      j = rep_index(v, C);
                if (j < 0)
                    continue;
                int j_dim = K->full_to_dim[j];
                if (j_dim < 0)
                    continue; /* v's orbit annihilated at this k */
                /* Phase factor e^{-i k·t_{t_R}} from the canonicalising
                 * translation T_{t_R} · u_ab = v. */
                double          phase = -k_dot_t_(C, K->k, t_R);
                double _Complex e_phase = cexp(I * phase);
                /* Coefficient √(σ_v / σ_u) — generalises √(N_u/N_v) to k≠Γ. */
                double _Complex coeff = csqrt(K->sigma[j_dim] / K->sigma[i_dim]);
                y[j_dim] += 0.5 * C->J * e_phase * coeff * x[i_dim];
            }
        }
        y[i_dim] += diag * x[i_dim];
    }
}

/* -------------------------------------------------------------------------- *
 * Main                                                                       *
 * -------------------------------------------------------------------------- */

int main(void) {
    printf("=== libirrep — momentum-resolved 3D Heisenberg ED ===\n");
    printf("    Cluster: SC 2×2×2  (N=8, dim=256)\n");
    printf("    Scanning all 8 allowed k-sectors at Sz=0.\n\n");

    cluster_t C = {0};
    build_cluster_(&C, IRREP_LATTICE3D_SC, 2, 2, 2, /*Sz=*/0, /*J=*/1.0);

    int Lx = 2, Ly = 2, Lz = 2;

    printf("  Sz=0 canonical reps: %d  (full Sz=0 dim = %d)\n", C.n_reps,
           1 << (C.N - 1) /* binomial guess */);
    printf("\n  k-sector summary:\n");
    printf("  %-26s  %-5s  %-12s  %-9s\n", "k = (kx, ky, kz)", "dim", "E_0", "time");

    double                gs_min = 1e300;
    int                   gs_min_k_index = -1;
    for (int n3 = 0; n3 < Lz; ++n3)
        for (int n2 = 0; n2 < Ly; ++n2)
            for (int n1 = 0; n1 < Lx; ++n1) {
                double f1 = (double)n1 / Lx;
                double f2 = (double)n2 / Ly;
                double f3 = (double)n3 / Lz;
                double k[3] = {f1 * C.b1[0] + f2 * C.b2[0] + f3 * C.b3[0],
                               f1 * C.b1[1] + f2 * C.b2[1] + f3 * C.b3[1],
                               f1 * C.b1[2] + f2 * C.b2[2] + f3 * C.b3[2]};

                ksector_t K = {0};
                build_ksector_(&K, &C, k);

                double E0 = 0;
                double t_ed = 0;
                if (K.dim > 0) {
                    double _Complex *seed = malloc((size_t)K.dim * sizeof(*seed));
                    double           sn = 0;
                    for (int i = 0; i < K.dim; ++i) {
                        seed[i] = 0.1 * sin(0.37 * (i + 1)) + I * 0.05 * cos(0.23 * (i + 1));
                        sn += creal(seed[i]) * creal(seed[i]) +
                              cimag(seed[i]) * cimag(seed[i]);
                    }
                    sn = sqrt(sn);
                    for (int i = 0; i < K.dim; ++i)
                        seed[i] /= sn;

                    double t1 = now_sec();
                    int    iters = K.dim < 80 ? K.dim : 80;
                    irrep_status_t st =
                        irrep_lanczos_eigvals(H_apply_k, &K, K.dim, 1, iters, seed, &E0);
                    t_ed = now_sec() - t1;
                    if (st != IRREP_OK) {
                        fprintf(stderr, "  Lanczos failed at k_index=%d\n",
                                (n3 * Ly + n2) * Lx + n1);
                        free(seed);
                        free_ksector_(&K);
                        continue;
                    }
                    free(seed);

                    if (E0 < gs_min) {
                        gs_min = E0;
                        gs_min_k_index = (n3 * Ly + n2) * Lx + n1;
                    }
                }

                /* Pretty-print k as fractions of the BZ (n_i / L_i units). */
                char k_label[64];
                snprintf(k_label, sizeof k_label, "(%.3f, %.3f, %.3f)·π", 2 * f1, 2 * f2,
                         2 * f3);
                if (K.dim > 0)
                    printf("  %-26s  %-5d  %+12.6f  %.3fs\n", k_label, K.dim, E0, t_ed);
                else
                    printf("  %-26s  %-5d  %-12s  --\n", k_label, K.dim,
                           "(annihilated)");

                free_ksector_(&K);
            }

    printf("\n  Minimum E_0(k) over all sectors = %+.6f J  at k_index=%d\n", gs_min,
           gs_min_k_index);
    printf("  Reference (full-space Lanczos, lattice3d_heisenberg_ed): E_0 = -4.820089 J\n");
    printf("  Difference = %+.2e\n", gs_min - (-4.820089));

    free_cluster_(&C);
    return 0;
}
