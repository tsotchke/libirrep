/* SPDX-License-Identifier: MIT */
#ifndef IRREP_SIMD_H
#define IRREP_SIMD_H

#include <stdbool.h>

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Populated once at library init (or first call). All fields reflect the
 * host CPU; dispatch tables select kernels based on these. */
typedef struct {
    bool neon;       /* aarch64 / arm64 only */
    bool sse42;
    bool avx2;
    bool avx512f;    /* reserved */
    bool avx512dq;   /* reserved */
    bool fma;
} irrep_cpu_features_t;

IRREP_API const irrep_cpu_features_t *irrep_cpu_features(void);

IRREP_API bool irrep_cpu_has_neon   (void);
IRREP_API bool irrep_cpu_has_sse42  (void);
IRREP_API bool irrep_cpu_has_avx2   (void);
IRREP_API bool irrep_cpu_has_avx512f(void);
IRREP_API bool irrep_cpu_has_fma    (void);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_SIMD_H */
