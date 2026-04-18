/* SPDX-License-Identifier: MIT */
/* libFuzzer entry point for irrep_tp_apply on a fixed descriptor:
 *   (1x0e + 1x1o) ⊗ (1x0e + 1x1o) → (1x0e + 1x1o + 1x2e)
 * Inputs are reinterpreted as the 4-element a and b real-basis features. */

#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#include <irrep/multiset.h>
#include <irrep/tensor_product.h>

static tp_descriptor_t *g_desc   = NULL;
static irrep_multiset_t *g_a_ms  = NULL;
static irrep_multiset_t *g_b_ms  = NULL;
static irrep_multiset_t *g_c_ms  = NULL;

static void init_once_(void) {
    if (g_desc) return;
    g_a_ms = irrep_multiset_parse("1x0e + 1x1o");
    g_b_ms = irrep_multiset_parse("1x0e + 1x1o");
    g_c_ms = irrep_multiset_parse("1x0e + 1x1o + 1x2e");
    int paths[32];
    int n = irrep_tp_enumerate_paths(g_a_ms, g_b_ms, g_c_ms, paths, 32);
    g_desc = irrep_tp_build(g_a_ms, g_b_ms, g_c_ms, paths, n);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    init_once_();
    if (!g_desc || size < 64) return 0;
    double a[4], b[4], c[9];
    memcpy(a, data,      sizeof(a));
    memcpy(b, data + 32, sizeof(b));
    for (int i = 0; i < 4; ++i) {
        if (!isfinite(a[i])) a[i] = 0.0;
        if (!isfinite(b[i])) b[i] = 0.0;
    }
    irrep_tp_apply(g_desc, a, b, c);
    for (int i = 0; i < 9; ++i) {
        if (!isfinite(c[i])) __builtin_trap();
    }
    return 0;
}
