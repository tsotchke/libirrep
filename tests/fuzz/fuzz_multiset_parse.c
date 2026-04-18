/* SPDX-License-Identifier: MIT */
/* libFuzzer entry point for irrep_multiset_parse.
 *
 * Build with:
 *   clang -fsanitize=fuzzer,address -O1 -g -Iinclude \
 *         tests/fuzz/fuzz_multiset_parse.c build/lib/liblibirrep.a \
 *         -lm -o fuzz_multiset_parse
 *
 * Then run:  ./fuzz_multiset_parse -max_total_time=3600
 */

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/multiset.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size > 4096) return 0;
    char *s = malloc(size + 1);
    if (!s) return 0;
    memcpy(s, data, size);
    s[size] = '\0';

    irrep_multiset_t *m = irrep_multiset_parse(s);
    if (m) {
        char buf[1024];
        (void)irrep_multiset_format(m, buf, sizeof(buf));
        irrep_multiset_simplify(m);
        (void)irrep_multiset_dim(m);
        irrep_multiset_free(m);
    }
    free(s);
    return 0;
}
