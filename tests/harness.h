/* SPDX-License-Identifier: MIT */
/* Minimal TAP-style C11 test harness — no framework dependency.
 *
 * Usage:
 *   #include "harness.h"
 *   int main(void) {
 *       IRREP_TEST_START("my-module");
 *       IRREP_ASSERT(some_condition);
 *       IRREP_ASSERT_NEAR(actual, expected, 1e-12);
 *       IRREP_ASSERT_NEAR_REL(actual, expected, 1e-10);
 *       return IRREP_TEST_END();
 *   }
 */
#ifndef LIBIRREP_TEST_HARNESS_H
#define LIBIRREP_TEST_HARNESS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Per-file counters — not thread-safe; tests run single-threaded. */
static int g_irrep_test_count     = 0;
static int g_irrep_test_fail      = 0;
static const char *g_irrep_test_name = "";

#define IRREP_TEST_START(name) \
    do { \
        g_irrep_test_name  = (name); \
        g_irrep_test_count = 0; \
        g_irrep_test_fail  = 0; \
        printf("# %s\n", g_irrep_test_name); \
    } while (0)

#define IRREP_TEST_END() \
    (printf("1..%d\n# %s: %d passed, %d failed\n", \
            g_irrep_test_count, g_irrep_test_name, \
            g_irrep_test_count - g_irrep_test_fail, g_irrep_test_fail), \
     g_irrep_test_fail == 0 ? 0 : 1)

#define IRREP_ASSERT(cond) \
    do { \
        g_irrep_test_count++; \
        if (cond) { \
            printf("ok %d - %s\n", g_irrep_test_count, #cond); \
        } else { \
            g_irrep_test_fail++; \
            printf("not ok %d - %s (at %s:%d)\n", \
                   g_irrep_test_count, #cond, __FILE__, __LINE__); \
        } \
    } while (0)

#define IRREP_ASSERT_NEAR(actual, expected, tol) \
    do { \
        g_irrep_test_count++; \
        double _a = (double)(actual); \
        double _e = (double)(expected); \
        double _t = (double)(tol); \
        double _d = fabs(_a - _e); \
        if (_d <= _t) { \
            printf("ok %d - |%.15g - %.15g| = %.3e <= %.3e\n", \
                   g_irrep_test_count, _a, _e, _d, _t); \
        } else { \
            g_irrep_test_fail++; \
            printf("not ok %d - |%.15g - %.15g| = %.3e > %.3e (at %s:%d)\n", \
                   g_irrep_test_count, _a, _e, _d, _t, __FILE__, __LINE__); \
        } \
    } while (0)

#define IRREP_ASSERT_NEAR_REL(actual, expected, rel_tol) \
    do { \
        g_irrep_test_count++; \
        double _a = (double)(actual); \
        double _e = (double)(expected); \
        double _t = (double)(rel_tol); \
        double _scale = fabs(_e) > 1e-300 ? fabs(_e) : 1.0; \
        double _d = fabs(_a - _e) / _scale; \
        if (_d <= _t) { \
            printf("ok %d - rel |%.15g - %.15g|/%.3g = %.3e <= %.3e\n", \
                   g_irrep_test_count, _a, _e, _scale, _d, _t); \
        } else { \
            g_irrep_test_fail++; \
            printf("not ok %d - rel |%.15g - %.15g|/%.3g = %.3e > %.3e (at %s:%d)\n", \
                   g_irrep_test_count, _a, _e, _scale, _d, _t, __FILE__, __LINE__); \
        } \
    } while (0)

#define IRREP_SKIP(reason) \
    do { \
        g_irrep_test_count++; \
        printf("ok %d - # SKIP %s\n", g_irrep_test_count, (reason)); \
    } while (0)

#endif /* LIBIRREP_TEST_HARNESS_H */
