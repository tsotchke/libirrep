/* SPDX-License-Identifier: MIT */
#ifndef IRREP_EXPORT_H
#define IRREP_EXPORT_H

#if defined(_WIN32) || defined(__CYGWIN__)
#  ifdef IRREP_BUILDING_DLL
#    define IRREP_API __declspec(dllexport)
#  else
#    define IRREP_API __declspec(dllimport)
#  endif
#elif (defined(__GNUC__) && __GNUC__ >= 4) || defined(__clang__)
#  define IRREP_API __attribute__((visibility("default")))
#else
#  define IRREP_API
#endif

#endif /* IRREP_EXPORT_H */
