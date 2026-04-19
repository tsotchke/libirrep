/* SPDX-License-Identifier: MIT */
/* libFuzzer entry point for the nequip spec-string parser.
 *
 * Feeds the raw input as a NUL-terminated C string directly to
 * `irrep_nequip_layer_from_spec`. The parser must either return a live
 * layer (which we free) or NULL + last-error text — no crash, no leak,
 * no read-past-NUL. A small stack buffer with a terminator guarantees
 * `strlen(spec)` stays bounded by the libFuzzer input.
 *
 * The fuzzer cap is intentionally modest (1 KB); the parser is linear in
 * spec length and longer inputs don't exercise new code paths. */

#include <stdint.h>
#include <stddef.h>
#include <string.h>

#include <irrep/nequip.h>

#define MAX_SPEC_LEN 1024u

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size > MAX_SPEC_LEN) size = MAX_SPEC_LEN;

    char spec[MAX_SPEC_LEN + 1];
    memcpy(spec, data, size);
    spec[size] = '\0';

    irrep_nequip_layer_t *layer = irrep_nequip_layer_from_spec(spec);
    if (layer) irrep_nequip_layer_free(layer);
    return 0;
}
