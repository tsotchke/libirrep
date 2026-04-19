/* SPDX-License-Identifier: MIT */
/* libFuzzer entry point for point-group projectors.
 *
 * Picks a group + irrep index from the first bytes of the input, then
 * interprets the remaining bytes as a feature vector of the fixed
 * multiset "2x0e + 1x1o + 1x2e" (dim = 2 + 3 + 5 = 10 doubles = 80 B).
 * Runs `irrep_pg_project` and asserts every output is finite. A NaN
 * through the character-weighted D^l(R) chain would indicate a real
 * defect (e.g., divide-by-sin(β) at gimbal lock, or a malformed
 * Wigner-D extraction).
 *
 * The multiset is fixed at init so the fuzzer doesn't waste cycles
 * rebuilding tables; each iteration only allocates the output buffer. */

#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include <irrep/multiset.h>
#include <irrep/point_group.h>

static irrep_pg_table_t *g_tables[4] = { NULL, NULL, NULL, NULL };
static irrep_multiset_t *g_spec       = NULL;
static int               g_num_irreps[4];

static void init_once_(void) {
    if (g_spec) return;
    g_spec = irrep_multiset_parse("2x0e + 1x1o + 1x2e");
    g_tables[0] = irrep_pg_table_build(IRREP_PG_C4V);
    g_tables[1] = irrep_pg_table_build(IRREP_PG_D6);
    g_tables[2] = irrep_pg_table_build(IRREP_PG_C3V);
    g_tables[3] = irrep_pg_table_build(IRREP_PG_D3);
    for (int i = 0; i < 4; ++i) {
        g_num_irreps[i] = g_tables[i] ? irrep_pg_num_irreps(g_tables[i]) : 0;
    }
}

#define IN_DIM  10  /* 2·1 + 1·3 + 1·5 */
#define MIN_INPUT (2u + IN_DIM * sizeof(double))

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    init_once_();
    if (!g_spec || size < MIN_INPUT) return 0;

    int group_idx = (int)(data[0] & 0x3);                       /* 0..3 */
    irrep_pg_table_t *t = g_tables[group_idx];
    if (!t) return 0;
    int mu = (int)(data[1] % (unsigned)g_num_irreps[group_idx]);

    double in[IN_DIM];
    memcpy(in, data + 2, IN_DIM * sizeof(double));
    for (int i = 0; i < IN_DIM; ++i) {
        if (!isfinite(in[i])) in[i] = 0.0;
    }

    double out[IN_DIM] = { 0 };
    irrep_pg_project(t, mu, g_spec, in, out);
    for (int i = 0; i < IN_DIM; ++i) {
        if (!isfinite(out[i])) __builtin_trap();
    }
    return 0;
}
