/* SPDX-License-Identifier: MIT */
#include <irrep/version.h>

int irrep_version_major(void) {
    return IRREP_VERSION_MAJOR;
}
int irrep_version_minor(void) {
    return IRREP_VERSION_MINOR;
}
int irrep_version_patch(void) {
    return IRREP_VERSION_PATCH;
}
const char *irrep_version_string(void) {
    return IRREP_VERSION_STRING;
}

/* ABI hash placeholder; real value is baked in by the release script
 * (scripts/generate_abi_hash.sh), which replaces the sentinel below. */
#ifndef IRREP_BAKED_ABI_HASH
#define IRREP_BAKED_ABI_HASH "unset-abi-hash-dev-build"
#endif

const char *irrep_abi_hash(void) {
    return IRREP_BAKED_ABI_HASH;
}
