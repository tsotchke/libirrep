/* SPDX-License-Identifier: MIT */
#include <stdbool.h>

#include <irrep/simd.h>

static irrep_cpu_features_t g_features;
static bool                 g_features_initialized = false;

static void detect_once(void) {
    if (g_features_initialized) return;

#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
    g_features.neon = true;
#endif

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#  if defined(__GNUC__) || defined(__clang__)
    __builtin_cpu_init();
    g_features.sse42    = __builtin_cpu_supports("sse4.2") ? true : false;
    g_features.avx2     = __builtin_cpu_supports("avx2")   ? true : false;
    g_features.fma      = __builtin_cpu_supports("fma")    ? true : false;
    g_features.avx512f  = __builtin_cpu_supports("avx512f")  ? true : false;
    g_features.avx512dq = __builtin_cpu_supports("avx512dq") ? true : false;
#  endif
#endif

    g_features_initialized = true;
}

const irrep_cpu_features_t *irrep_cpu_features(void) {
    detect_once();
    return &g_features;
}

bool irrep_cpu_has_neon   (void) { detect_once(); return g_features.neon;    }
bool irrep_cpu_has_sse42  (void) { detect_once(); return g_features.sse42;   }
bool irrep_cpu_has_avx2   (void) { detect_once(); return g_features.avx2;    }
bool irrep_cpu_has_avx512f(void) { detect_once(); return g_features.avx512f; }
bool irrep_cpu_has_fma    (void) { detect_once(); return g_features.fma;     }
