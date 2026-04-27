/* SPDX-License-Identifier: MIT */
/* Tests for irrep_sg_rep_table_* — orbit-representative enumeration and
 * rep → linear-index lookup.
 *
 * Covers: Burnside sum rule, per-rep canonical-form property, sorted
 * order, binary-search lookup round-trip, rejection of non-reps and
 * wrong-popcount queries, trivial-sector edge cases (popcount = 0, N). */
#include "harness.h"
#include <irrep/config_project.h>
#include <irrep/lattice.h>
#include <irrep/space_group.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Portable temp-file path. Works regardless of CWD (make runs from
 * project root; ctest runs from build/). Uses TMPDIR / TMP / TEMP
 * environment vars and falls back to "." on toolchains that expose
 * none of them. */
static const char *tmp_path_(const char *suffix) {
    static char buf[256];
    const char *tmp = getenv("TMPDIR");
    if (!tmp) tmp = getenv("TMP");
    if (!tmp) tmp = getenv("TEMP");
    if (!tmp) tmp = ".";
    snprintf(buf, sizeof buf, "%s/libirrep_%s.bin", tmp, suffix);
    return buf;
}

static int popcount_u64(uint64_t x) {
    int c = 0;
    while (x) { x &= x - 1; ++c; }
    return c;
}

static long long binom(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    long long r = 1;
    for (int i = 0; i < k; ++i)
        r = r * (n - i) / (i + 1);
    return r;
}

int main(void) {
    IRREP_TEST_START("rep_table");

    /* ---- Burnside sum rule on kagome 2×2 (N = 12, p6mm, |G| = 48) ------
     * Build table at every popcount 0..12 and verify
     *   Σ_k orbit_size(k) == C(N, p).
     * Also verify: every stored rep is canonical (== own orbit minimum)
     * and count matches count-of-canonical-forms from direct scan. */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int N = irrep_space_group_num_sites(G);
        IRREP_ASSERT(N == 12);

        for (int p = 0; p <= N; ++p) {
            irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, p);
            IRREP_ASSERT(T != NULL);
            long long n = irrep_sg_rep_table_count(T);

            /* Burnside: Σ orbit_size = C(N, p). */
            long long sum = 0;
            for (long long k = 0; k < n; ++k) {
                uint64_t rep = irrep_sg_rep_table_get(T, k);
                int      osz = irrep_sg_rep_table_orbit_size(T, k);
                IRREP_ASSERT(osz > 0);
                sum += osz;

                /* Every stored rep matches its own canonical form. */
                uint64_t rep_canon; int g_idx;
                irrep_sg_canonicalise(G, rep, &rep_canon, &g_idx);
                IRREP_ASSERT(rep_canon == rep);

                /* Popcount matches the table's fixed popcount. */
                IRREP_ASSERT(popcount_u64(rep) == p);
            }
            IRREP_ASSERT(sum == binom(N, p));

            /* Sorted ascending. */
            for (long long k = 1; k < n; ++k)
                IRREP_ASSERT(irrep_sg_rep_table_get(T, k) >
                             irrep_sg_rep_table_get(T, k - 1));

            irrep_sg_rep_table_free(T);
        }

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- rep → index round-trip on kagome 2×2 Sz = 0 --------------------
     * For every config c with popcount = 6, canonicalise(c) → rep; then
     * rep_table_index(rep) must return a valid k ∈ [0, count), and
     * rep_table_get(T, k) must equal rep. */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        int N = irrep_space_group_num_sites(G);

        irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 6);
        IRREP_ASSERT(T != NULL);
        long long n = irrep_sg_rep_table_count(T);
        IRREP_ASSERT(n == 30); /* established in test_canonicalise */

        int verified = 0;
        for (uint64_t c = 0; c < (1ULL << N); ++c) {
            if (popcount_u64(c) != 6)
                continue;
            uint64_t rep; int g_idx;
            irrep_sg_canonicalise(G, c, &rep, &g_idx);
            long long k = irrep_sg_rep_table_index(T, rep);
            IRREP_ASSERT(k >= 0 && k < n);
            IRREP_ASSERT(irrep_sg_rep_table_get(T, k) == rep);
            ++verified;
        }
        IRREP_ASSERT(verified == (int)binom(N, 6)); /* 924 */

        irrep_sg_rep_table_free(T);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Non-rep / wrong-popcount rejection ----------------------------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 6);
        IRREP_ASSERT(T != NULL);

        /* A popcount-5 config isn't in the popcount-6 table. */
        IRREP_ASSERT(irrep_sg_rep_table_index(T, 0x1F) == -1);

        /* A non-canonical config at the right popcount isn't in the table.
         * Pick a random popcount-6 config and test: unless it's a rep it
         * must NOT be in the table, even if popcount matches. */
        int non_rep_count = 0;
        for (uint64_t c = 0; c < (1ULL << 12); ++c) {
            if (popcount_u64(c) != 6)
                continue;
            uint64_t rep; int g_idx;
            irrep_sg_canonicalise(G, c, &rep, &g_idx);
            if (c != rep) {
                IRREP_ASSERT(irrep_sg_rep_table_index(T, c) == -1);
                if (++non_rep_count >= 50) break;
            }
        }
        IRREP_ASSERT(non_rep_count > 0);

        /* NULL table returns 0/-1 without crashing. */
        IRREP_ASSERT(irrep_sg_rep_table_count(NULL) == 0);
        IRREP_ASSERT(irrep_sg_rep_table_get(NULL, 0) == 0);
        IRREP_ASSERT(irrep_sg_rep_table_orbit_size(NULL, 0) == 0);
        IRREP_ASSERT(irrep_sg_rep_table_index(NULL, 0) == -1);

        irrep_sg_rep_table_free(T);
        irrep_sg_rep_table_free(NULL); /* no-op */
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Popcount edge cases: 0 and N each yield a single orbit -------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);

        irrep_sg_rep_table_t *T0 = irrep_sg_rep_table_build(G, 0);
        IRREP_ASSERT(T0 != NULL);
        IRREP_ASSERT(irrep_sg_rep_table_count(T0) == 1);
        IRREP_ASSERT(irrep_sg_rep_table_get(T0, 0) == 0);
        IRREP_ASSERT(irrep_sg_rep_table_orbit_size(T0, 0) == 1); /* G fixes |0⟩ */
        IRREP_ASSERT(irrep_sg_rep_table_index(T0, 0) == 0);
        irrep_sg_rep_table_free(T0);

        irrep_sg_rep_table_t *TN = irrep_sg_rep_table_build(G, 12);
        IRREP_ASSERT(TN != NULL);
        IRREP_ASSERT(irrep_sg_rep_table_count(TN) == 1);
        uint64_t all_up = (1ULL << 12) - 1;
        IRREP_ASSERT(irrep_sg_rep_table_get(TN, 0) == all_up);
        IRREP_ASSERT(irrep_sg_rep_table_orbit_size(TN, 0) == 1);
        IRREP_ASSERT(irrep_sg_rep_table_index(TN, all_up) == 0);
        irrep_sg_rep_table_free(TN);

        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- p1 on 4×4 square, popcount 4: translation-only representatives
     * N = 16 sites, order = 16 (translations). C(16,4) = 1820. With only
     * translations, stabilisers come from the lattice periodicity. */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_SQUARE, 4, 4);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P1);
        int order = irrep_space_group_order(G); /* 16 */
        IRREP_ASSERT(order == 16);

        irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 4);
        IRREP_ASSERT(T != NULL);
        long long n = irrep_sg_rep_table_count(T);

        long long sum = 0;
        for (long long k = 0; k < n; ++k)
            sum += irrep_sg_rep_table_orbit_size(T, k);
        IRREP_ASSERT(sum == binom(16, 4)); /* 1820 */

        irrep_sg_rep_table_free(T);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    /* ---- Save/load round-trip ------------------------------------------- */
    {
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);

        irrep_sg_rep_table_t *T = irrep_sg_rep_table_build(G, 6);
        IRREP_ASSERT(T != NULL);
        long long cnt = irrep_sg_rep_table_count(T);

        /* Save to /tmp. */
        const char *path = tmp_path_("rep_table_test");
        int rc = irrep_sg_rep_table_save(T, path);
        IRREP_ASSERT(rc == 0);

        /* Load back with matching metadata → reproduces the original. */
        irrep_sg_rep_table_t *T2 = irrep_sg_rep_table_load(G, 6, path);
        IRREP_ASSERT(T2 != NULL);
        IRREP_ASSERT(irrep_sg_rep_table_count(T2) == cnt);
        for (long long k = 0; k < cnt; ++k) {
            IRREP_ASSERT(irrep_sg_rep_table_get(T2, k) ==
                         irrep_sg_rep_table_get(T, k));
            IRREP_ASSERT(irrep_sg_rep_table_orbit_size(T2, k) ==
                         irrep_sg_rep_table_orbit_size(T, k));
        }

        /* Popcount mismatch → load returns NULL. */
        IRREP_ASSERT(irrep_sg_rep_table_load(G, 5, path) == NULL);
        /* Bad path → NULL. */
        IRREP_ASSERT(irrep_sg_rep_table_load(G, 6, "/nonexistent/path") == NULL);
        /* NULL arguments → NULL. */
        IRREP_ASSERT(irrep_sg_rep_table_load(NULL, 6, path) == NULL);
        IRREP_ASSERT(irrep_sg_rep_table_load(G, 6, NULL) == NULL);
        IRREP_ASSERT(irrep_sg_rep_table_save(NULL, path) == -1);
        IRREP_ASSERT(irrep_sg_rep_table_save(T, NULL) == -1);

        /* Group mismatch: save p6mm, try to load with p1 G'. */
        irrep_lattice_t     *L2 = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *Gp = irrep_space_group_build(L2, IRREP_WALLPAPER_P1);
        IRREP_ASSERT(irrep_sg_rep_table_load(Gp, 6, path) == NULL);
        irrep_space_group_free(Gp);
        irrep_lattice_free(L2);

        irrep_sg_rep_table_free(T2);

        /* Corrupted-file rejection: tamper with the magic byte and
         * verify load cleanly returns NULL rather than segfaulting or
         * returning garbage. */
        {
            FILE *fp = fopen(path, "r+b");
            IRREP_ASSERT(fp != NULL);
            char bad = 'X';
            fwrite(&bad, 1, 1, fp); /* overwrite first magic byte */
            fclose(fp);
            IRREP_ASSERT(irrep_sg_rep_table_load(G, 6, path) == NULL);
        }

        /* Wrong-version rejection: write a valid-magic file with an
         * unknown version field — load must reject. */
        {
            const char *vp = tmp_path_("rep_table_badver");
            FILE *fp = fopen(vp, "wb");
            IRREP_ASSERT(fp != NULL);
            char magic[8] = {'I','R','R','E','P','_','R','T'};
            uint32_t bad_version = 9999;
            int dummy_int = 0;
            long long dummy_ll = 0;
            fwrite(magic, 1, 8, fp);
            fwrite(&bad_version, sizeof(uint32_t), 1, fp);
            fwrite(&dummy_int, sizeof(int), 1, fp);
            fwrite(&dummy_int, sizeof(int), 1, fp);
            fwrite(&dummy_int, sizeof(int), 1, fp);
            fwrite(&dummy_ll, sizeof(long long), 1, fp);
            fclose(fp);
            IRREP_ASSERT(irrep_sg_rep_table_load(G, 6, vp) == NULL);
            remove(vp);
        }

        /* Truncated file: zero-length → rejection. */
        {
            const char *tp = tmp_path_("rep_table_trunc");
            FILE *fp = fopen(tp, "wb");
            IRREP_ASSERT(fp != NULL);
            fclose(fp);
            IRREP_ASSERT(irrep_sg_rep_table_load(G, 6, tp) == NULL);
            remove(tp);
        }

        irrep_sg_rep_table_free(T);
        irrep_space_group_free(G);
        irrep_lattice_free(L);
        remove(path);
    }

    /* ---- Invalid builds return NULL ------------------------------------- */
    {
        IRREP_ASSERT(irrep_sg_rep_table_build(NULL, 0) == NULL);
        irrep_lattice_t     *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 2, 2);
        irrep_space_group_t *G = irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
        IRREP_ASSERT(irrep_sg_rep_table_build(G, -1) == NULL);
        IRREP_ASSERT(irrep_sg_rep_table_build(G, 13) == NULL); /* > N */
        irrep_space_group_free(G);
        irrep_lattice_free(L);
    }

    return IRREP_TEST_END();
}
