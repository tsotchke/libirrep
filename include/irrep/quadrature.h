/* SPDX-License-Identifier: MIT */
#ifndef IRREP_QUADRATURE_H
#define IRREP_QUADRATURE_H

#include <stdbool.h>

#include <irrep/export.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Lebedev rules on S^2 (available orders:
 *   3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41). */
IRREP_API int  irrep_lebedev_size(int order);                      /* points; 0 if unknown order */
IRREP_API bool irrep_lebedev_fill(int order, double *xyz_weights); /* writes size * 4 doubles */

/* Gauss-Legendre nodes & weights on [-1, 1]. */
IRREP_API bool irrep_gauss_legendre(int n, double *nodes, double *weights);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_QUADRATURE_H */
