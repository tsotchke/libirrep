/* SPDX-License-Identifier: MIT */
/* Non-libFuzzer driver for the LLVMFuzzerTestOneInput harnesses in this
 * directory. Apple's shipped clang strips libclang_rt.fuzzer_osx.a, and
 * mixing a Homebrew clang front-end with Apple's ld produces a linker
 * mismatch on aarch64. This driver feeds the same harnesses with a
 * deterministic PCG64-style stream so they can run under plain toolchains
 * and in CI without requiring libFuzzer.
 *
 * Correctness trade-off versus libFuzzer: we do not get coverage-guided
 * mutation; we get plain uniform-random bytes. Effective for catching
 * crashes and assertion violations in small-input harnesses (all of ours
 * read <= 256 bytes). Not effective at reaching paths gated by specific
 * magic values — but no libirrep harness does that.
 *
 * Usage:  driver [seed] [iters] [max_bytes]
 * Defaults: seed=42, iters=1000000, max_bytes=256.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size);

int main(int argc, char **argv) {
    uint64_t seed      = (argc > 1) ? strtoull(argv[1], NULL, 10) : 42ULL;
    uint64_t iters     = (argc > 2) ? strtoull(argv[2], NULL, 10) : 1000000ULL;
    size_t   max_bytes = (argc > 3) ? (size_t)strtoul(argv[3], NULL, 10) : 256;
    if (max_bytes == 0 || max_bytes > 4096) max_bytes = 256;

    /* 64-bit LCG (MMIX parameters); state becomes the random stream. */
    uint64_t state = seed ? seed : 0xDEADBEEFCAFEBABEULL;
    uint8_t *buf = (uint8_t*)malloc(max_bytes);
    if (!buf) { fprintf(stderr, "driver: out of memory\n"); return 1; }

    for (uint64_t it = 0; it < iters; ++it) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        size_t sz = 1 + (size_t)(state & (max_bytes - 1));
        for (size_t i = 0; i < sz; ++i) {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            buf[i] = (uint8_t)(state >> 56);
        }
        LLVMFuzzerTestOneInput(buf, sz);

        if ((it & 0xFFFFF) == 0 && it > 0) {
            /* progress every ~1M iters, to stderr so stdout stays clean. */
            fprintf(stderr, "  .. %llu iters\n",
                    (unsigned long long)it);
        }
    }

    free(buf);
    printf("fuzz-driver: %llu iterations ok\n",
           (unsigned long long)iters);
    return 0;
}
