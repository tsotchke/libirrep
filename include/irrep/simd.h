/* SPDX-License-Identifier: MIT */
/** @file simd.h
 *  @brief Runtime CPU feature detection. The library's batched kernels route
 *         through an internal dispatch table that reads these flags once at
 *         load and picks the widest available kernel per operation.
 *
 *  Populated lazily on first access, thread-safe via a C11 atomic guard.
 */
#ifndef IRREP_SIMD_H
#define IRREP_SIMD_H

#include <stdbool.h>

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Host CPU feature mask reflecting what the library's dispatch table
 *         will use. All fields are one-way: they don't reset within a process. */
typedef struct {
    bool neon;       /**< aarch64 NEON (always true on arm64). */
    bool sse42;      /**< x86 SSE4.2. */
    bool avx2;       /**< x86 AVX2. */
    bool avx512f;    /**< x86 AVX-512 foundation (reserved, no kernels yet). */
    bool avx512dq;   /**< x86 AVX-512 doubleword/quadword (reserved). */
    bool fma;        /**< x86 FMA3. */
} irrep_cpu_features_t;

/** @brief Pointer to the lazy-initialised feature struct. */
IRREP_API const irrep_cpu_features_t *irrep_cpu_features(void);

/** @brief Convenience: NEON present. */
IRREP_API bool irrep_cpu_has_neon   (void);
/** @brief Convenience: SSE4.2 present. */
IRREP_API bool irrep_cpu_has_sse42  (void);
/** @brief Convenience: AVX2 present. */
IRREP_API bool irrep_cpu_has_avx2   (void);
/** @brief Convenience: AVX-512F present. */
IRREP_API bool irrep_cpu_has_avx512f(void);
/** @brief Convenience: FMA present. */
IRREP_API bool irrep_cpu_has_fma    (void);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SIMD_H */
