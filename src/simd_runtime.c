/* SPDX-License-Identifier: MIT */
/* Runtime CPU feature detection + dispatch table population.
 *
 * Features are detected once (lazy, thread-safe via C11 atomic flag). The
 * dispatch table is populated from the detected features on first access —
 * widest available kernel wins, falls back to the scalar reference.
 *
 * NEON is implied by aarch64 (mandatory in ARMv8-A). x86 features use the
 * gcc/clang `__builtin_cpu_supports` protocol (portable across gcc ≥ 4.8
 * and clang ≥ 3.8). Windows MSVC fallback left for a future pass.
 */
#include <stdatomic.h>
#include <stdbool.h>
#include <stddef.h>

#include <irrep/simd.h>

#include "internal/dispatch.h"

/* Init synchronisation:
 *
 *   `g_initialized` (acquire/release) is the fast-path read. `g_init_lock`
 *   (a C11 atomic_flag used as a spin-lock) serialises the slow-path
 *   initialisation so there is no data race on the non-atomic writes to
 *   `g_features` and `g_dispatch`. A prior revision used only the atomic
 *   flag with a "benign race" comment; under the C11 memory model that was
 *   formally UB (concurrent writes to non-atomic memory are racy regardless
 *   of whether the values agree).
 *
 *   The spin is bounded by the init cost (microseconds), only hit by the
 *   very first set of concurrent callers at library load, and there is no
 *   link-time dependency on pthread or <threads.h>. */
static irrep_cpu_features_t g_features;
static irrep_dispatch_t     g_dispatch;
static atomic_bool          g_initialized = false;
static atomic_flag          g_init_lock   = ATOMIC_FLAG_INIT;

static void detect_features(void) {
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
}

static void populate_dispatch(void) {
    /* Start with scalar fallbacks — always correct. */
    g_dispatch.cutoff_cosine_batch        = irrep_cutoff_cosine_batch_scalar;
    g_dispatch.cutoff_cosine_d_batch      = irrep_cutoff_cosine_d_batch_scalar;
    g_dispatch.cutoff_polynomial_batch    = irrep_cutoff_polynomial_batch_scalar;
    g_dispatch.cutoff_polynomial_d_batch  = irrep_cutoff_polynomial_d_batch_scalar;
    g_dispatch.rbf_bessel_batch           = irrep_rbf_bessel_batch_scalar;
    g_dispatch.sph_harm_cart_all_batch    = irrep_sph_harm_cart_all_batch_scalar;

#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
    if (g_features.neon) {
        /* Pure multiply/add — vectorizes cleanly with float64x2_t. The cosine
         * variants still route through the scalar path: vectorized cos() would
         * require either Accelerate (non-portable) or a polynomial kernel that
         * isn't worth the precision trade-off at this size. */
        g_dispatch.cutoff_polynomial_batch   = irrep_cutoff_polynomial_batch_neon;
        g_dispatch.cutoff_polynomial_d_batch = irrep_cutoff_polynomial_d_batch_neon;
        /* SH batch: lane-paired iteration over edges. Bit-exact against the
         * scalar kernel on odd-N-included fixtures. */
        g_dispatch.sph_harm_cart_all_batch   = irrep_sph_harm_cart_all_batch_neon;
    }
#endif

#if (defined(__x86_64__) || defined(_M_X64)) && defined(__AVX2__) && defined(__FMA__)
    extern void irrep_sph_harm_cart_all_batch_avx2(int, size_t, const double *, double *);
    if (g_features.avx2 && g_features.fma) {
        /* SH batch: 4 edges per __m256d, Chebyshev trig recurrence
         * vectorised across lanes; per-lane Legendre grid. Bit-exact
         * against the scalar kernel on tail-included fixtures. */
        g_dispatch.sph_harm_cart_all_batch = irrep_sph_harm_cart_all_batch_avx2;
    }
#endif
}

static void init_once(void) {
    /* Fast path: already initialised. */
    if (atomic_load_explicit(&g_initialized, memory_order_acquire)) return;

    /* Slow path: acquire the init lock. Losers spin until the winner
     * publishes `g_initialized`. */
    while (atomic_flag_test_and_set_explicit(&g_init_lock, memory_order_acquire)) {
        if (atomic_load_explicit(&g_initialized, memory_order_acquire)) return;
    }

    /* Double-check after locking — another thread may have finished between
     * our initial fast-path load and our lock acquisition. */
    if (!atomic_load_explicit(&g_initialized, memory_order_relaxed)) {
        detect_features();
        populate_dispatch();
        atomic_store_explicit(&g_initialized, true, memory_order_release);
    }
    atomic_flag_clear_explicit(&g_init_lock, memory_order_release);
}

const irrep_cpu_features_t *irrep_cpu_features(void) {
    init_once();
    return &g_features;
}

const irrep_dispatch_t *irrep_dispatch_get(void) {
    init_once();
    return &g_dispatch;
}

bool irrep_cpu_has_neon   (void) { init_once(); return g_features.neon;    }
bool irrep_cpu_has_sse42  (void) { init_once(); return g_features.sse42;   }
bool irrep_cpu_has_avx2   (void) { init_once(); return g_features.avx2;    }
bool irrep_cpu_has_avx512f(void) { init_once(); return g_features.avx512f; }
bool irrep_cpu_has_fma    (void) { init_once(); return g_features.fma;     }
