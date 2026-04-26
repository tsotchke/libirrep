/* SPDX-License-Identifier: MIT */
/* Translation-symmetry-resolved Heisenberg ED on a 3D Bravais lattice.
 *
 * Demonstrates the standard Γ-momentum sector reduction on a small 3D
 * cluster: orbit-canonicalise spin configurations under translation,
 * build the symmetric basis |Γ, u⟩ = (1/√|orb_u|) Σ_g |T_g · u⟩, apply
 * the Heisenberg operator in that basis, then Lanczos for the ground
 * state. The same technique drives 2D kagome / square ED at N≥27 in
 * libirrep's `irrep_heisenberg_apply_in_sector`; this example shows the
 * pattern using only `lattice3d_translate` so it works on every 3D
 * family without waiting for a full 3D space-group module.
 *
 * Cluster: BCC 2×2×2 = 16 sites, Sz=0 sector (C(16,8)=12870 states),
 * Γ momentum (k=0). Expected reduction: full dim 65536 → Sz=0
 * 12870 → Γ ≈ 1600. GS energy must match the full-space result
 * (E₀ = −20.000 J exactly — the K₈,₈ closed form). */

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

/* -------------------------------------------------------------------------- *
 * Sector context                                                             *
 * -------------------------------------------------------------------------- */

typedef struct {
    int      N;            /* number of sites */
    int      n_bonds;      /* NN bond count */
    int     *bi;           /* bond endpoints (length n_bonds) */
    int     *bj;
    double   J;            /* Heisenberg coupling */
    int      n_trans;      /* number of translations (= num_cells) */
    int    **trans_perm;   /* trans_perm[t][site] = translated site */

    /* Canonical reps in the (Sz=0, Γ=0) sector. */
    int       n_reps;      /* number of distinct orbits */
    uint64_t *reps;        /* canonical bitstring per orbit (sorted asc) */
    int      *orbit_size;  /* |orb_u| per orbit */

    /* Hash from any bitstring → canonical rep index (we use a sorted
     * array + bsearch on canonical reps, so we lookup by canonicalising
     * any state and binary-searching). */
} sector_t;

static int cmp_u64(const void *a, const void *b) {
    uint64_t ua = *(const uint64_t *)a, ub = *(const uint64_t *)b;
    return ua < ub ? -1 : ua > ub ? 1 : 0;
}

/* Apply permutation `perm` to a bitstring: bit i of input maps to bit perm[i]
 * of output. */
static uint64_t apply_perm_bits(uint64_t s, const int *perm, int N) {
    uint64_t r = 0;
    for (int i = 0; i < N; ++i)
        if ((s >> i) & 1ULL)
            r |= 1ULL << perm[i];
    return r;
}

/* Canonical rep of `s` = minimum over its translation orbit. */
static uint64_t canonicalise(uint64_t s, const sector_t *S) {
    uint64_t best = s;
    for (int t = 0; t < S->n_trans; ++t) {
        uint64_t r = apply_perm_bits(s, S->trans_perm[t], S->N);
        if (r < best)
            best = r;
    }
    return best;
}

/* Given a canonical rep, return its index in the sorted reps array.
 * Returns -1 if not found (caller must canonicalise first). */
static int rep_index(uint64_t v, const sector_t *S) {
    /* Binary search. */
    int lo = 0, hi = S->n_reps - 1;
    while (lo <= hi) {
        int      m = (lo + hi) / 2;
        uint64_t mv = S->reps[m];
        if (mv == v)
            return m;
        if (mv < v)
            lo = m + 1;
        else
            hi = m - 1;
    }
    return -1;
}

/* Heisenberg apply in (Sz=0, Γ) sector. Signature matches the Lanczos
 * apply-callback contract so we can pass `S` as ctx. */
static void H_apply_sector(const double _Complex *x, double _Complex *y, void *ctx) {
    sector_t *S = ctx;
    memset(y, 0, (size_t)S->n_reps * sizeof(double _Complex));
    for (int i = 0; i < S->n_reps; ++i) {
        if (x[i] == 0.0)
            continue;
        uint64_t u = S->reps[i];
        double   diag = 0.0;
        for (int b = 0; b < S->n_bonds; ++b) {
            int a = S->bi[b], q = S->bj[b];
            int sa = (int)((u >> a) & 1ULL);
            int sq = (int)((u >> q) & 1ULL);
            if (sa == sq) {
                diag += 0.25 * S->J;
            } else {
                diag -= 0.25 * S->J;
                /* Off-diagonal: flip and canonicalise. */
                uint64_t up = u ^ (1ULL << a) ^ (1ULL << q);
                uint64_t v = canonicalise(up, S);
                int      j = rep_index(v, S);
                if (j >= 0) {
                    /* Matrix element ⟨Γ,v|H|Γ,u⟩ = ½J · √(N_u/N_v) · k_uv,
                     * where N_u, N_v are orbit sizes — NB the source orbit is
                     * in the numerator (cf. hamiltonian.c:234). */
                    double scale = sqrt((double)S->orbit_size[i] / (double)S->orbit_size[j]);
                    y[j] += 0.5 * S->J * scale * x[i];
                }
            }
        }
        y[i] += diag * x[i];
    }
}

/* -------------------------------------------------------------------------- *
 * Sector build                                                               *
 * -------------------------------------------------------------------------- */

static void build_sector_(sector_t *S, irrep_lattice3d_kind_t kind, int Lx, int Ly, int Lz,
                          int Sz_target) {
    irrep_lattice3d_t *L = irrep_lattice3d_build(kind, Lx, Ly, Lz);
    S->N = irrep_lattice3d_num_sites(L);
    S->n_bonds = irrep_lattice3d_num_bonds_nn(L);
    S->bi = malloc((size_t)S->n_bonds * sizeof(int));
    S->bj = malloc((size_t)S->n_bonds * sizeof(int));
    irrep_lattice3d_fill_bonds_nn(L, S->bi, S->bj);
    S->J = 1.0;

    /* Translation table: every (dx, dy, dz) ∈ [0,Lx) × [0,Ly) × [0,Lz). */
    S->n_trans = Lx * Ly * Lz;
    S->trans_perm = malloc((size_t)S->n_trans * sizeof(*S->trans_perm));
    int t = 0;
    for (int dz = 0; dz < Lz; ++dz)
        for (int dy = 0; dy < Ly; ++dy)
            for (int dx = 0; dx < Lx; ++dx) {
                S->trans_perm[t] = malloc((size_t)S->N * sizeof(int));
                for (int s = 0; s < S->N; ++s)
                    S->trans_perm[t][s] = irrep_lattice3d_translate(L, s, dx, dy, dz);
                ++t;
            }

    /* Enumerate Sz-target configurations (popcount = N/2 + Sz_target).
     * Find canonical rep for each; deduplicate; sort ascending. */
    int      target_pop = S->N / 2 + Sz_target;
    long long total_dim = 1LL << S->N;
    uint64_t *canonical_buf = malloc((size_t)total_dim * sizeof(uint64_t));
    int       n_can = 0;
    for (long long s = 0; s < total_dim; ++s) {
        if (__builtin_popcountll((unsigned long long)s) != target_pop)
            continue;
        uint64_t c = canonicalise((uint64_t)s, S);
        if ((uint64_t)s == c)
            canonical_buf[n_can++] = c; /* `s` is its own canonical */
    }
    qsort(canonical_buf, (size_t)n_can, sizeof(uint64_t), cmp_u64);
    S->n_reps = n_can;
    S->reps = malloc((size_t)n_can * sizeof(uint64_t));
    memcpy(S->reps, canonical_buf, (size_t)n_can * sizeof(uint64_t));
    free(canonical_buf);

    /* Orbit size for each canonical rep. */
    S->orbit_size = malloc((size_t)n_can * sizeof(int));
    for (int i = 0; i < n_can; ++i) {
        uint64_t u = S->reps[i];
        int      seen[256];
        int      n_seen = 0;
        (void)seen;
        (void)n_seen;
        /* Orbit size = # of distinct translates. Track set of distinct images. */
        uint64_t images[1024];
        int      n_images = 0;
        for (int q = 0; q < S->n_trans; ++q) {
            uint64_t r = apply_perm_bits(u, S->trans_perm[q], S->N);
            int      seen_q = 0;
            for (int k = 0; k < n_images; ++k)
                if (images[k] == r) {
                    seen_q = 1;
                    break;
                }
            if (!seen_q)
                images[n_images++] = r;
        }
        S->orbit_size[i] = n_images;
    }

    irrep_lattice3d_free(L);
}

static void free_sector_(sector_t *S) {
    free(S->bi);
    free(S->bj);
    for (int t = 0; t < S->n_trans; ++t)
        free(S->trans_perm[t]);
    free(S->trans_perm);
    free(S->reps);
    free(S->orbit_size);
}

/* -------------------------------------------------------------------------- *
 * Main                                                                       *
 * -------------------------------------------------------------------------- */

int main(void) {
    printf("=== libirrep — 3D translation-sector Heisenberg ED ===\n\n");

    sector_t S = {0};
    double   t0 = now_sec();
    build_sector_(&S, IRREP_LATTICE3D_BCC, 2, 2, 2, /*Sz=*/0);
    double t_build = now_sec() - t0;

    long long full_dim = 1LL << S.N;
    long long sz0_dim = 0;
    for (int i = 0; i < S.n_reps; ++i)
        sz0_dim += S.orbit_size[i];

    printf("  Cluster: BCC 2×2×2  N=%d  bonds=%d  translations=%d\n", S.N, S.n_bonds, S.n_trans);
    printf("  Full Hilbert dim       = %lld\n", full_dim);
    printf("  Sz=0 sector dim        = %lld   (= Σ |orb_u|)\n", sz0_dim);
    printf("  Γ-momentum sector dim  = %d   (canonical reps; %.1f× shrink vs full)\n",
           S.n_reps, (double)full_dim / (double)S.n_reps);
    printf("  Build wall-clock       = %.2fs\n\n", t_build);

    /* Lanczos seed: deterministic spread across canonical reps. */
    double _Complex *seed = malloc((size_t)S.n_reps * sizeof(*seed));
    double           sn = 0;
    for (int i = 0; i < S.n_reps; ++i) {
        seed[i] = 0.1 * sin(0.37 * i) + I * 0.05 * cos(0.23 * i);
        sn += creal(seed[i]) * creal(seed[i]) + cimag(seed[i]) * cimag(seed[i]);
    }
    sn = sqrt(sn);
    for (int i = 0; i < S.n_reps; ++i)
        seed[i] /= sn;

    double t1 = now_sec();
    double E0 = 0;
    irrep_status_t st =
        irrep_lanczos_eigvals(H_apply_sector, &S, S.n_reps, 1, 200, seed, &E0);
    double t_ed = now_sec() - t1;
    if (st != IRREP_OK) {
        fprintf(stderr, "Lanczos failed, status=%d\n", (int)st);
        return 1;
    }

    printf("  Γ-sector ground state\n");
    printf("    E_0 = %+.6f J   (full-space reference: −20.000000 J — K₈,₈ closed form)\n", E0);
    printf("    E_0/N_bonds = %+.6f J\n", E0 / S.n_bonds);
    printf("    Lanczos wall-clock = %.2fs\n", t_ed);

    free(seed);
    free_sector_(&S);
    return 0;
}
