/* SPDX-License-Identifier: MIT */
#ifndef IRREP_VERSION_H
#define IRREP_VERSION_H

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

#define IRREP_VERSION_MAJOR 1
#define IRREP_VERSION_MINOR 0
#define IRREP_VERSION_PATCH 0
#define IRREP_VERSION_STRING "1.0.0-dev"
#define IRREP_VERSION_HEX   ((IRREP_VERSION_MAJOR << 16) | \
                             (IRREP_VERSION_MINOR <<  8) | \
                             (IRREP_VERSION_PATCH <<  0))

IRREP_API int         irrep_version_major (void);
IRREP_API int         irrep_version_minor (void);
IRREP_API int         irrep_version_patch (void);
IRREP_API const char *irrep_version_string(void);

/* ABI hash: sha256 of sorted exported-symbol list + public struct layouts,
 * baked in at build time. Consumers can compare against the value reported
 * by the installed header to catch binary/header drift. */
IRREP_API const char *irrep_abi_hash(void);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_VERSION_H */
