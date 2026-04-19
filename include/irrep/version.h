/* SPDX-License-Identifier: MIT */
/** @file version.h
 *  @brief Compile-time version macros and runtime accessors, plus the ABI
 *         hash for header/binary drift detection.
 */
#ifndef IRREP_VERSION_H
#define IRREP_VERSION_H

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

#define IRREP_VERSION_MAJOR 1                  /**< Major: API-incompatible changes. */
#define IRREP_VERSION_MINOR 3                  /**< Minor: additive API. */
#define IRREP_VERSION_PATCH 0                  /**< Patch: bug fixes / internals. */
#define IRREP_VERSION_STRING "1.3.0-alpha"     /**< Human-readable version string. */
/** Packed integer version, e.g. `0x010203` for 1.2.3 — useful for preprocessor
 *  checks like `#if IRREP_VERSION_HEX >= 0x010100`. */
#define IRREP_VERSION_HEX   ((IRREP_VERSION_MAJOR << 16) | \
                             (IRREP_VERSION_MINOR <<  8) | \
                             (IRREP_VERSION_PATCH <<  0))

/** @brief Runtime major version (may differ from header if mis-linked). */
IRREP_API int         irrep_version_major (void);
/** @brief Runtime minor version. */
IRREP_API int         irrep_version_minor (void);
/** @brief Runtime patch version. */
IRREP_API int         irrep_version_patch (void);
/** @brief Runtime version string. */
IRREP_API const char *irrep_version_string(void);

/** @brief SHA-256 of the sorted exported-symbol list concatenated with public
 *         struct layouts, baked in at build time. A mismatch between this and
 *         the value recorded at header-install time indicates ABI drift. */
IRREP_API const char *irrep_abi_hash(void);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_VERSION_H */
