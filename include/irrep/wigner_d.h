/* SPDX-License-Identifier: MIT */
#ifndef IRREP_WIGNER_D_H
#define IRREP_WIGNER_D_H

#include <complex.h>

#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Small Wigner-d scalar (real, depends on beta only). */
IRREP_API double irrep_wigner_d_small   (int j, int mp, int m, double beta);
IRREP_API double irrep_wigner_d_small_2j(int two_j, int two_mp, int two_m, double beta);

/* Full Wigner-D scalar: D^j_{m'm}(alpha, beta, gamma). */
IRREP_API double _Complex irrep_wigner_D   (int j, int mp, int m,
                                            double alpha, double beta, double gamma);
IRREP_API double _Complex irrep_wigner_D_2j(int two_j, int two_mp, int two_m,
                                            double alpha, double beta, double gamma);

/* Full matrix (2j+1) x (2j+1), row-major, with m' as row, m as column. */
IRREP_API void irrep_wigner_D_matrix(int j, double _Complex *out,
                                     double alpha, double beta, double gamma);
IRREP_API void irrep_wigner_d_matrix(int j, double          *out, double beta);

/* Block-diagonal D on an irrep_multiset_t; writes `total_dim x total_dim` complex. */
IRREP_API void irrep_wigner_D_multiset(const irrep_multiset_t *m, double _Complex *out,
                                       double alpha, double beta, double gamma);

/* Derivative d/dβ of small-d (for force evaluation). */
IRREP_API double irrep_wigner_d_small_dbeta(int j, int mp, int m, double beta);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_WIGNER_D_H */
