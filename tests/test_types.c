/* SPDX-License-Identifier: MIT */
#include "harness.h"
#include <irrep/types.h>
#include <irrep/version.h>
#include <irrep/multiset.h>
#include <string.h>

int main(void) {
    IRREP_TEST_START("types");

    /* Status strings exist and cover every enum value. */
    IRREP_ASSERT(irrep_strerror(IRREP_OK) != NULL);
    IRREP_ASSERT(strcmp(irrep_strerror(IRREP_OK), "ok") == 0);
    IRREP_ASSERT(strcmp(irrep_strerror(IRREP_ERR_INVALID_ARG),    "invalid argument") == 0);
    IRREP_ASSERT(strcmp(irrep_strerror(IRREP_ERR_OUT_OF_MEMORY),  "out of memory") == 0);
    IRREP_ASSERT(strcmp(irrep_strerror(IRREP_ERR_SELECTION_RULE), "selection rule violation") == 0);
    IRREP_ASSERT(strcmp(irrep_strerror(IRREP_ERR_NOT_IMPLEMENTED),"not implemented") == 0);
    IRREP_ASSERT(strcmp(irrep_strerror(IRREP_ERR_PRECONDITION),   "precondition failed") == 0);
    IRREP_ASSERT(strcmp(irrep_strerror(IRREP_ERR_PARSE),          "parse error") == 0);
    /* Out-of-range enum value falls through to the "unknown status" branch. */
    IRREP_ASSERT(strcmp(irrep_strerror((irrep_status_t)9999),     "unknown status") == 0);

    /* Thread-local last-error buffer. Before any library call it should be
     * the empty string; a failing constructor call populates it. */
    IRREP_ASSERT(irrep_last_error() != NULL);
    /* (can't assume empty since prior tests in the same process may have set it) */

    /* Trigger a failing multiset parse to populate last-error. */
    IRREP_ASSERT(irrep_multiset_parse("not a valid multiset spec") == NULL);
    const char *msg = irrep_last_error();
    IRREP_ASSERT(msg != NULL && msg[0] != '\0');
    /* Reading twice returns the same pointer/content. */
    IRREP_ASSERT(strcmp(irrep_last_error(), msg) == 0);

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
