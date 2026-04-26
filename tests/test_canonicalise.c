/* SPDX-License-Identifier: MIT */
/* Tests for the orbit-representative canonical form primitives:
 *   - irrep_space_group_apply_bits (bitstring site permutation)
 *   - irrep_sg_canonicalise        (orbit minimum + reaching g)
 *   - irrep_sg_orbit_size          (|orbit| = |G| / |stabiliser|)
 *
 * These are the building blocks for sparse-Lanczos-in-sector at N ≥ 27
 * where representatives replace full-Hilbert vectors. Correctness here
 * is load-bearing for every downstream ED / NQS / DMRG handshake. */
#include "harness.h"
#include <irrep/config_project.h>
#include <irrep/lattice.h>
#include <irrep/space_group.h>
#include <stdlib.h>
#include <string.h>

/* Count set bits in a uint64_t (popcount). */
static int popcount_u64(uint64_t x) {
    int c = 0;
    while (x) {
        x &= x - 1;
        ++c;
    }
    return c;
}

int main(void) {
    IRREP_TEST_START("canonicalise");

    /* ---- Identity action: g = 0 leaves every config unchanged ---------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        IRREP_ASSERT(G != NULL);
        int N = irrep_space_group_num_sites(G); /* 12 */
        IRREP_ASSERT(N == 12);

        /* Every config is fixed by g = 0 (identity element). */
        for (uint64_t c = 0; c < 200; ++c) {
            uint64_t out = irrep_space_group_apply_bits(G, 0, c);
            IRREP_ASSERT(out == c);
        }

        /* Out-of-range g returns input unchanged (soft error). */
        IRREP_ASSERT(irrep_space_group_apply_bits(G, -1,     0x123) == 0x123);
        IRREP_ASSERT(irrep_space_group_apply_bits(G, 9999,   0x123) == 0x123);
        IRREP_ASSERT(irrep_space_group_apply_bits(NULL, 0,   0x123) == 0x123);

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- apply_bits preserves popcount (site-permutation invariant) ---- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int order = irrep_space_group_order(G); /* 48 */

        uint64_t probes[] = {0, 1, 0xAAAULL, 0xE38ULL, 0xFC0ULL};
        for (unsigned p = 0; p < sizeof(probes) / sizeof(probes[0]); ++p) {
            int pc_in = popcount_u64(probes[p]);
            for (int g = 0; g < order; ++g) {
                uint64_t out = irrep_space_group_apply_bits(G, g, probes[p]);
                IRREP_ASSERT(popcount_u64(out) == pc_in);
            }
        }

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Canonicalise idempotence + invariance under group action ------
     * canonicalise(rep) == rep for any rep; canonicalise(g · c) == canonicalise(c)
     * for any c and any g. */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int order = irrep_space_group_order(G);

        uint64_t probes[] = {0, 1, 0xABC, 0x555, 0xF0F, 0x333};
        for (unsigned p = 0; p < sizeof(probes) / sizeof(probes[0]); ++p) {
            uint64_t c = probes[p] & 0xFFFULL; /* clamp to 12 sites */
            uint64_t rep; int g_idx;
            irrep_sg_canonicalise(G, c, &rep, &g_idx);

            /* (1) Idempotence: canonicalise(rep) == rep. */
            uint64_t rep2; int g2;
            irrep_sg_canonicalise(G, rep, &rep2, &g2);
            IRREP_ASSERT(rep2 == rep);

            /* (2) g_idx · c == rep. */
            uint64_t check = irrep_space_group_apply_bits(G, g_idx, c);
            IRREP_ASSERT(check == rep);

            /* (3) Orbit invariance: canonicalise(h · c) == rep for every h. */
            for (int h = 0; h < order; ++h) {
                uint64_t hc = irrep_space_group_apply_bits(G, h, c);
                uint64_t r2; int g2b;
                irrep_sg_canonicalise(G, hc, &r2, &g2b);
                IRREP_ASSERT(r2 == rep);
            }
        }

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Orbit-size · (# reps) sum rule on kagome 12 at Sz = 0 ---------
     * Enumerate every 12-bit config with popcount = 6 (C(12,6) = 924).
     * Canonicalise each, keep only those that are their own rep. Sum of
     * orbit sizes over representatives must equal 924. This is the
     * Burnside / orbit-counting theorem in computational form. */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int N = irrep_space_group_num_sites(G); /* 12 */
        int order = irrep_space_group_order(G); /* 48 */
        uint64_t D = 1ULL << N;

        int total_configs = 0;
        int num_reps      = 0;
        int sum_orbit     = 0;
        for (uint64_t c = 0; c < D; ++c) {
            if (popcount_u64(c) != 6)
                continue;
            ++total_configs;

            uint64_t rep; int g_idx;
            irrep_sg_canonicalise(G, c, &rep, &g_idx);
            if (rep == c) {
                ++num_reps;
                int osz = irrep_sg_orbit_size(G, c);
                IRREP_ASSERT(osz > 0);
                IRREP_ASSERT(order % osz == 0); /* orbit size divides |G| */
                sum_orbit += osz;
            }
        }
        IRREP_ASSERT(total_configs == 924); /* C(12, 6) */
        IRREP_ASSERT(sum_orbit == total_configs);
        /* The count of representatives must be less than 924 / order plus
         * a handful (non-trivial stabilisers), and more than 924 / order minus
         * (no reps possible below that). For Sz = 0 kagome 12 on p6mm: 30. */
        IRREP_ASSERT(num_reps == 30);

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Fully-polarised states: unique orbit, full stabiliser --------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int order = irrep_space_group_order(G);

        /* |0⟩ (all down): stabiliser = G, orbit = {|0⟩}. */
        uint64_t rep; int g_idx;
        irrep_sg_canonicalise(G, 0, &rep, &g_idx);
        IRREP_ASSERT(rep == 0);
        IRREP_ASSERT(irrep_sg_orbit_size(G, 0) == 1);

        /* |1…1⟩ (all up on 12 sites): same. */
        uint64_t all_up = (1ULL << 12) - 1;
        irrep_sg_canonicalise(G, all_up, &rep, &g_idx);
        IRREP_ASSERT(rep == all_up);
        IRREP_ASSERT(irrep_sg_orbit_size(G, all_up) == 1);

        /* Single-up state (12 sites, popcount 1): orbit size divides order. */
        uint64_t single = 1ULL;
        IRREP_ASSERT(order % irrep_sg_orbit_size(G, single) == 0);

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Canonicalise reaches the minimum g on ties ---------------------
     * When a config is in the stabiliser orbit of multiple g's that map
     * to the same rep, the API must return the smallest such g. */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int order = irrep_space_group_order(G);

        /* |0⟩ is fixed by every g, so canonicalise returns g = 0. */
        int g_idx;
        irrep_sg_canonicalise(G, 0, NULL, &g_idx);
        IRREP_ASSERT(g_idx == 0);

        /* Scan all configs; verify the returned g_idx really maps c → rep
         * and is the minimum such. */
        for (uint64_t c = 0; c < 512; ++c) {
            uint64_t rep_c; int g_c;
            irrep_sg_canonicalise(G, c, &rep_c, &g_c);
            IRREP_ASSERT(irrep_space_group_apply_bits(G, g_c, c) == rep_c);
            /* No earlier g reaches rep_c from c. */
            for (int h = 0; h < g_c; ++h)
                IRREP_ASSERT(irrep_space_group_apply_bits(G, h, c) != rep_c);
            (void)order;
        }

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- rep_out = NULL skipping: caller can request just g_idx -------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);

        int g_only;
        irrep_sg_canonicalise(G, 0x55, NULL, &g_only);
        IRREP_ASSERT(g_only >= 0);

        uint64_t rep_only;
        irrep_sg_canonicalise(G, 0x55, &rep_only, NULL);

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Stabiliser: full group on trivial configs, trivial on generic
     * configs, orbit_size · stabiliser_size = |G| always. -------------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int order = irrep_space_group_order(G); /* 48 */
        int *stab = malloc(sizeof(int) * order);

        /* |0⟩ (all spins down) is fixed by every g. */
        int n_stab = irrep_sg_stabiliser(G, 0, stab);
        IRREP_ASSERT(n_stab == order);
        for (int g = 1; g < order; ++g)
            IRREP_ASSERT(stab[g] > stab[g - 1]); /* ascending */
        IRREP_ASSERT(stab[0] == 0);

        /* |all-up⟩ (0xFFF) same. */
        n_stab = irrep_sg_stabiliser(G, 0xFFFULL, stab);
        IRREP_ASSERT(n_stab == order);

        /* Generic config: orbit-stabiliser theorem. */
        uint64_t probes[] = {0x1ULL, 0x3ULL, 0x55ULL, 0xAFULL};
        for (unsigned p = 0; p < sizeof(probes)/sizeof(probes[0]); ++p) {
            int ns = irrep_sg_stabiliser(G, probes[p], stab);
            int os = irrep_sg_orbit_size(G, probes[p]);
            IRREP_ASSERT(ns > 0);
            IRREP_ASSERT(os > 0);
            IRREP_ASSERT(ns * os == order);
            /* Each returned g genuinely fixes config. */
            for (int k = 0; k < ns; ++k)
                IRREP_ASSERT(irrep_space_group_apply_bits(G, stab[k], probes[p]) == probes[p]);
        }

        /* Error paths. */
        IRREP_ASSERT(irrep_sg_stabiliser(NULL, 0, stab) == 0);
        IRREP_ASSERT(irrep_sg_stabiliser(G, 0, NULL) == 0);

        free(stab);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- p1 (translations only): orbit = {c, T·c, T²·c, …} ------------
     * Sum rule sanity on a 3×3 square with local dim 2, popcount 4:
     * C(9, 4) = 126 total configs; p1 order = 9; expected rep count is
     * close to 14 (with a few stabiliser-boosted orbits). */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 3, 3);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
        int N = irrep_space_group_num_sites(G); /* 9 */
        int order = irrep_space_group_order(G); /* 9 */

        int total = 0, sum_osz = 0;
        for (uint64_t c = 0; c < (1ULL << N); ++c) {
            if (popcount_u64(c) != 4)
                continue;
            ++total;
            uint64_t rep; int g_idx;
            irrep_sg_canonicalise(G, c, &rep, &g_idx);
            if (rep == c)
                sum_osz += irrep_sg_orbit_size(G, c);
        }
        IRREP_ASSERT(total == 126);
        IRREP_ASSERT(sum_osz == 126);
        IRREP_ASSERT(order == 9);

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    return IRREP_TEST_END();
}
