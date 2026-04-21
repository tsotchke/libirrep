/* SPDX-License-Identifier: MIT */
/*
 * Thread-local error channel for the library.
 *
 * Every failing entry point returns an `irrep_status_t` enum; `irrep_strerror`
 * names the enumerator, and `irrep_set_error_` (internal) records a formatted
 * message consulted via `irrep_last_error`. Storage is `_Thread_local` so
 * concurrent callers do not race on each other's messages, and so a caller
 * observing a non-OK status on this thread sees the string that was written
 * on this thread.
 */
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include <irrep/types.h>

#if defined(_MSC_VER)
#define IRREP_THREAD __declspec(thread)
#else
#define IRREP_THREAD _Thread_local
#endif

static IRREP_THREAD char g_last_error[256];

/* Internal helper used by module implementations. */
void irrep_set_error_(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    (void)vsnprintf(g_last_error, sizeof(g_last_error), fmt, ap);
    va_end(ap);
}

const char *irrep_last_error(void) {
    return g_last_error[0] ? g_last_error : "";
}

const char *irrep_strerror(irrep_status_t status) {
    switch (status) {
    case IRREP_OK:
        return "ok";
    case IRREP_ERR_INVALID_ARG:
        return "invalid argument";
    case IRREP_ERR_OUT_OF_MEMORY:
        return "out of memory";
    case IRREP_ERR_SELECTION_RULE:
        return "selection rule violation";
    case IRREP_ERR_NOT_IMPLEMENTED:
        return "not implemented";
    case IRREP_ERR_PRECONDITION:
        return "precondition failed";
    case IRREP_ERR_PARSE:
        return "parse error";
    }
    return "unknown status";
}
