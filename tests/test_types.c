/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/types.h>
#include <irrep/version.h>

int main(void) {
    IRREP_TEST_START("types");

    /* Status strings exist. */
    IRREP_ASSERT(irrep_strerror(IRREP_OK) != NULL);
    IRREP_ASSERT(strcmp(irrep_strerror(IRREP_OK), "ok") == 0);

    /* Limits sensible. */
    IRREP_ASSERT(IRREP_L_MAX >= 8);
    IRREP_ASSERT(IRREP_TWO_J_MAX >= 2 * IRREP_L_MAX);

    /* Parity constants. */
    IRREP_ASSERT(IRREP_EVEN == +1);
    IRREP_ASSERT(IRREP_ODD  == -1);

    /* Version accessors match macros. */
    IRREP_ASSERT(irrep_version_major() == IRREP_VERSION_MAJOR);
    IRREP_ASSERT(irrep_version_minor() == IRREP_VERSION_MINOR);
    IRREP_ASSERT(irrep_version_patch() == IRREP_VERSION_PATCH);
    IRREP_ASSERT(strcmp(irrep_version_string(), IRREP_VERSION_STRING) == 0);

    return IRREP_TEST_END();
}
