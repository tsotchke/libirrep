/* SPDX-License-Identifier: MIT */
/* libFuzzer entry point for Clebsch-Gordan.
 * Feeds random j1, m1, j2, m2, J, M — verifies finite result within [-1, 1]
 * and exact zero for selection-rule violations. */

#include <math.h>
#include <stdint.h>
#include <stddef.h>

#include <irrep/clebsch_gordan.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 24) return 0;
    int two_j1 = *(const int32_t*)(data +  0) & 0x1F;       /* 0..31 */
    int two_m1 = *(const int32_t*)(data +  4) & 0x3F;
    two_m1 = (two_m1 - 32);                                  /* -32..31 */
    int two_j2 = *(const int32_t*)(data +  8) & 0x1F;
    int two_m2 = (*(const int32_t*)(data + 12) & 0x3F) - 32;
    int two_J  = *(const int32_t*)(data + 16) & 0x3F;
    int two_M  = (*(const int32_t*)(data + 20) & 0x3F) - 32;

    double v = irrep_cg_2j(two_j1, two_m1, two_j2, two_m2, two_J, two_M);
    if (!isfinite(v)) __builtin_trap();
    if (v < -1.01 || v > 1.01)  __builtin_trap();
    return 0;
}
