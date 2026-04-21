/* SPDX-License-Identifier: MIT */
/** @file export.h
 *  @brief Symbol-visibility macro.
 *
 *  The library builds with `-fvisibility=hidden` on gcc/clang; public
 *  symbols are tagged @c IRREP_API to opt back in to default visibility.
 *  On Windows the macro resolves to `__declspec(dllexport)` while the DLL
 *  is being built (define @c IRREP_BUILDING_DLL) and `__declspec(dllimport)`
 *  for consumers.
 *
 *  This header deliberately omits `extern "C"` guards: it contains no
 *  declarations, only a preprocessor macro. The guards would be inert
 *  and would give a false signal of C/C++ cross-language linkage content.
 */
#ifndef IRREP_EXPORT_H
#define IRREP_EXPORT_H

#if defined(_WIN32) || defined(__CYGWIN__)
#ifdef IRREP_BUILDING_DLL
#define IRREP_API __declspec(dllexport)
#else
#define IRREP_API __declspec(dllimport)
#endif
#elif (defined(__GNUC__) && __GNUC__ >= 4) || defined(__clang__)
#define IRREP_API __attribute__((visibility("default")))
#else
#define IRREP_API
#endif

#endif /* IRREP_EXPORT_H */
