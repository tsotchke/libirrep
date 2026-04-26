/* SPDX-License-Identifier: MIT */
/* Orbit-representative table for sparse Lanczos in symmetry-projected
 * Hilbert-space sectors. See include/irrep/config_project.h for the
 * public API.
 *
 * Enumeration strategy: iterate every bitstring with popcount = p using
 * Gosper's next-lex-greater-with-same-popcount formula, canonicalise
 * each, keep only those that are their own canonical form. Because
 * canonicalise returns the orbit's lex-minimum and iteration is
 * monotone increasing, representatives are accepted in sorted order for
 * free. Binary search backs rep → index lookup.
 *
 * Build cost: O(C(N, p) · |G| · N). At kagome 4×4, N = 48, p = 24:
 * infeasible (C(48,24) ≈ 2.4e13). At N = 30, p = 15: C ≈ 1.55e8; a few
 * minutes on a single CPU. The table is built once per sector and
 * reused across every Lanczos iteration. */
#include <irrep/config_project.h>
#include <irrep/space_group.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct irrep_sg_rep_table {
    const irrep_space_group_t *G;        /* borrowed, not freed */
    int                        popcount;
    long long                  count;
    uint64_t                  *reps;       /* sorted ascending */
    int                       *orbit_sizes; /* parallel to reps */
};

static int popcount_u64_(uint64_t x) {
    int c = 0;
    while (x) {
        x &= x - 1;
        ++c;
    }
    return c;
}

/* Gosper's hack: given x with popcount p, return the next larger uint64_t
 * with the same popcount, or 0 if no such value fits in 64 bits. */
static uint64_t next_same_popcount_(uint64_t x) {
    if (x == 0)
        return 0;
    uint64_t c = x & (uint64_t)(-(int64_t)x);
    uint64_t r = x + c;
    /* r == 0 means we overflowed — no next value in 64 bits. */
    if (r == 0)
        return 0;
    return (((r ^ x) >> 2) / c) | r;
}

/* Early-exit canonicalise: returns 1 iff `config` is its own orbit rep.
 * Stops the loop as soon as any g produces a strictly smaller image. */
static int is_canonical_(const irrep_space_group_t *G, uint64_t config, int order) {
    for (int g = 1; g < order; ++g) { /* skip g = 0 (identity) */
        uint64_t image = irrep_space_group_apply_bits(G, g, config);
        if (image < config)
            return 0;
    }
    return 1;
}

irrep_sg_rep_table_t *irrep_sg_rep_table_build(const irrep_space_group_t *G, int popcount) {
    if (!G || popcount < 0)
        return NULL;
    int N = irrep_space_group_num_sites(G);
    if (N <= 0 || N > 64 || popcount > N)
        return NULL;

    int order = irrep_space_group_order(G);

    irrep_sg_rep_table_t *T = calloc(1, sizeof(*T));
    if (!T)
        return NULL;
    T->G        = G;
    T->popcount = popcount;

    /* Reserve an initial capacity; grow by doubling. Upper bound is
     * C(N, p) / |G| (tight up to stabiliser corrections). */
    long long cap = 64;
    T->reps         = malloc((size_t)cap * sizeof(uint64_t));
    T->orbit_sizes  = malloc((size_t)cap * sizeof(int));
    if (!T->reps || !T->orbit_sizes) {
        irrep_sg_rep_table_free(T);
        return NULL;
    }

    long long count = 0;
    /* Starting config: the lex-smallest uint64_t with popcount p is
     * 0b0...0 1...1 (popcount ones at the low bits). */
    uint64_t c;
    if (popcount == 0) {
        c = 0;
        if (is_canonical_(G, c, order)) {
            T->reps[count]         = c;
            T->orbit_sizes[count]  = irrep_sg_orbit_size(G, c);
            ++count;
        }
    } else if (popcount == 64) {
        c = ~(uint64_t)0;
        if (is_canonical_(G, c, order)) {
            T->reps[count]         = c;
            T->orbit_sizes[count]  = irrep_sg_orbit_size(G, c);
            ++count;
        }
    } else {
        /* Upper bound: the largest bitstring occupying only sites [0, N).
         * For N == 64 that's ~0; for N < 64 it's (1 << N) − 1. Gosper must
         * terminate at this cap — otherwise at popcount = N the iteration
         * ticks forever through the exponentially-many 64-bit popcount-N
         * values with high bits set. */
        uint64_t upper_bound = (N == 64) ? ~(uint64_t)0
                                         : ((uint64_t)1 << N) - 1;
        c = ((uint64_t)1 << popcount) - 1;
        while (c != 0 && c <= upper_bound) {
            if (is_canonical_(G, c, order)) {
                if (count == cap) {
                    cap *= 2;
                    uint64_t *new_reps = realloc(T->reps, (size_t)cap * sizeof(uint64_t));
                    int      *new_osz  = realloc(T->orbit_sizes, (size_t)cap * sizeof(int));
                    if (!new_reps || !new_osz) {
                        free(new_reps ? new_reps : T->reps);
                        free(new_osz ? new_osz : T->orbit_sizes);
                        T->reps        = NULL;
                        T->orbit_sizes = NULL;
                        free(T);
                        return NULL;
                    }
                    T->reps        = new_reps;
                    T->orbit_sizes = new_osz;
                }
                T->reps[count]        = c;
                T->orbit_sizes[count] = irrep_sg_orbit_size(G, c);
                ++count;
            }
            uint64_t nxt = next_same_popcount_(c);
            if (nxt == 0 || nxt > upper_bound)
                break;
            c = nxt;
        }
    }

    T->count = count;

    /* Trim overallocation to exact size. */
    if (count > 0) {
        uint64_t *trim_r = realloc(T->reps, (size_t)count * sizeof(uint64_t));
        int      *trim_o = realloc(T->orbit_sizes, (size_t)count * sizeof(int));
        if (trim_r)
            T->reps = trim_r;
        if (trim_o)
            T->orbit_sizes = trim_o;
    }

    /* Sanity: reps are already sorted (Gosper iteration is monotone). */
    return T;
}

void irrep_sg_rep_table_free(irrep_sg_rep_table_t *T) {
    if (!T)
        return;
    free(T->reps);
    free(T->orbit_sizes);
    free(T);
}

long long irrep_sg_rep_table_count(const irrep_sg_rep_table_t *T) {
    return T ? T->count : 0;
}

uint64_t irrep_sg_rep_table_get(const irrep_sg_rep_table_t *T, long long k) {
    if (!T || k < 0 || k >= T->count)
        return 0;
    return T->reps[k];
}

int irrep_sg_rep_table_orbit_size(const irrep_sg_rep_table_t *T, long long k) {
    if (!T || k < 0 || k >= T->count)
        return 0;
    return T->orbit_sizes[k];
}

const irrep_space_group_t *irrep_sg_rep_table_space_group(const irrep_sg_rep_table_t *T) {
    return T ? T->G : NULL;
}

/* On-disk format:
 *   magic       : 8 bytes "IRREP_RT"
 *   version     : uint32 (currently 1)
 *   popcount    : int32
 *   num_sites   : int32 (cross-check against G at load time)
 *   group_order : int32 (cross-check)
 *   count       : int64
 *   reps        : count × uint64
 *   orbit_sizes : count × int32 */

static const char REP_TABLE_MAGIC[8] = {'I','R','R','E','P','_','R','T'};
static const uint32_t REP_TABLE_VERSION = 1;

int irrep_sg_rep_table_save(const irrep_sg_rep_table_t *T, const char *path) {
    if (!T || !path)
        return -1;
    FILE *f = fopen(path, "wb");
    if (!f)
        return -1;

    int num_sites = irrep_space_group_num_sites(T->G);
    int order     = irrep_space_group_order(T->G);
    int pc        = T->popcount;
    long long cnt = T->count;

    if (fwrite(REP_TABLE_MAGIC, 1, 8, f) != 8) goto fail;
    if (fwrite(&REP_TABLE_VERSION, sizeof(uint32_t), 1, f) != 1) goto fail;
    if (fwrite(&pc,        sizeof(int),   1, f) != 1) goto fail;
    if (fwrite(&num_sites, sizeof(int),   1, f) != 1) goto fail;
    if (fwrite(&order,     sizeof(int),   1, f) != 1) goto fail;
    if (fwrite(&cnt,       sizeof(long long), 1, f) != 1) goto fail;
    if (cnt > 0) {
        if ((long long)fwrite(T->reps,        sizeof(uint64_t), (size_t)cnt, f) != cnt) goto fail;
        if ((long long)fwrite(T->orbit_sizes, sizeof(int),      (size_t)cnt, f) != cnt) goto fail;
    }
    fclose(f);
    return 0;
fail:
    fclose(f);
    return -1;
}

irrep_sg_rep_table_t *irrep_sg_rep_table_load(const irrep_space_group_t *G,
                                              int popcount, const char *path) {
    if (!G || !path)
        return NULL;
    FILE *f = fopen(path, "rb");
    if (!f)
        return NULL;

    char     magic[8];
    uint32_t ver;
    int      file_pc, file_N, file_order;
    long long cnt;

    if (fread(magic, 1, 8, f) != 8) goto fail;
    if (memcmp(magic, REP_TABLE_MAGIC, 8) != 0) goto fail;
    if (fread(&ver, sizeof(uint32_t), 1, f) != 1) goto fail;
    if (ver != REP_TABLE_VERSION) goto fail;
    if (fread(&file_pc,    sizeof(int), 1, f) != 1) goto fail;
    if (fread(&file_N,     sizeof(int), 1, f) != 1) goto fail;
    if (fread(&file_order, sizeof(int), 1, f) != 1) goto fail;
    if (fread(&cnt,        sizeof(long long), 1, f) != 1) goto fail;

    /* Cross-check against caller-supplied group. */
    if (file_pc    != popcount)                           goto fail;
    if (file_N     != irrep_space_group_num_sites(G))     goto fail;
    if (file_order != irrep_space_group_order(G))         goto fail;
    if (cnt < 0)                                          goto fail;

    irrep_sg_rep_table_t *T = calloc(1, sizeof(*T));
    if (!T) goto fail;
    T->G        = G;
    T->popcount = popcount;
    T->count    = cnt;
    if (cnt > 0) {
        T->reps        = malloc((size_t)cnt * sizeof(uint64_t));
        T->orbit_sizes = malloc((size_t)cnt * sizeof(int));
        if (!T->reps || !T->orbit_sizes) {
            irrep_sg_rep_table_free(T);
            goto fail;
        }
        if ((long long)fread(T->reps,        sizeof(uint64_t), (size_t)cnt, f) != cnt) {
            irrep_sg_rep_table_free(T);
            goto fail;
        }
        if ((long long)fread(T->orbit_sizes, sizeof(int),      (size_t)cnt, f) != cnt) {
            irrep_sg_rep_table_free(T);
            goto fail;
        }
    }
    fclose(f);
    return T;
fail:
    fclose(f);
    return NULL;
}

long long irrep_sg_rep_table_index(const irrep_sg_rep_table_t *T, uint64_t rep) {
    if (!T || T->count == 0)
        return -1;
    /* Popcount gate — non-matching popcount cannot be in the table. */
    if (popcount_u64_(rep) != T->popcount)
        return -1;
    /* Binary search on sorted reps[]. */
    long long lo = 0, hi = T->count - 1;
    while (lo <= hi) {
        long long mid  = lo + (hi - lo) / 2;
        uint64_t  midv = T->reps[mid];
        if (midv == rep)
            return mid;
        if (midv < rep)
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return -1;
}
