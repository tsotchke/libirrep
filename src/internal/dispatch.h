/* SPDX-License-Identifier: MIT */
/* Private dispatch table: function pointers swapped at load time to the best
 * available SIMD kernel for the host. Not a public header — do not ship.
 *
 * Each public batched kernel has a `_scalar` reference implementation and
 * (optionally) `_neon` / `_avx2` / `_avx512` variants. At library load the
 * pointer is set to the widest kernel whose CPU features are present, else
 * the scalar fallback. The public entry point is a trampoline through the
 * pointer — one indirect call per batched invocation, negligible overhead.
 */
#ifndef IRREP_INTERNAL_DISPATCH_H
#define IRREP_INTERNAL_DISPATCH_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Radial / cutoff batched kernels. */
typedef void (*irrep_fn_cutoff_cosine_batch)     (size_t n, const double *r, double r_cut, double *out);
typedef void (*irrep_fn_cutoff_cosine_d_batch)   (size_t n, const double *r, double r_cut, double *out);
typedef void (*irrep_fn_cutoff_polynomial_batch) (size_t n, const double *r, double r_cut, int p, double *out);
typedef void (*irrep_fn_cutoff_polynomial_d_batch)(size_t n, const double *r, double r_cut, int p, double *out);
typedef void (*irrep_fn_rbf_bessel_batch)        (int n_deg, size_t n, const double *r, double r_cut, double *out);

/* Spherical-harmonic batched kernel (hot path for NequIP / MACE messages). */
typedef void (*irrep_fn_sph_harm_cart_all_batch)(int l_max, size_t n, const double *r_hats, double *out);

typedef struct {
    irrep_fn_cutoff_cosine_batch        cutoff_cosine_batch;
    irrep_fn_cutoff_cosine_d_batch      cutoff_cosine_d_batch;
    irrep_fn_cutoff_polynomial_batch    cutoff_polynomial_batch;
    irrep_fn_cutoff_polynomial_d_batch  cutoff_polynomial_d_batch;
    irrep_fn_rbf_bessel_batch           rbf_bessel_batch;
    irrep_fn_sph_harm_cart_all_batch    sph_harm_cart_all_batch;
} irrep_dispatch_t;

/* Returns the populated dispatch table (thread-safe, lazy init). */
const irrep_dispatch_t *irrep_dispatch_get(void);

/* Scalar reference implementations (always available, always correct). */
void irrep_cutoff_cosine_batch_scalar        (size_t n, const double *r, double r_cut, double *out);
void irrep_cutoff_cosine_d_batch_scalar      (size_t n, const double *r, double r_cut, double *out);
void irrep_cutoff_polynomial_batch_scalar    (size_t n, const double *r, double r_cut, int p, double *out);
void irrep_cutoff_polynomial_d_batch_scalar  (size_t n, const double *r, double r_cut, int p, double *out);
void irrep_rbf_bessel_batch_scalar           (int n_deg, size_t n, const double *r, double r_cut, double *out);
void irrep_sph_harm_cart_all_batch_scalar    (int l_max, size_t n, const double *r_hats, double *out);

/* NEON kernels (aarch64 only; prototypes always visible for the dispatch
 * table, but the compiled-in definitions live in src/radial_neon.c which is
 * only part of the build on aarch64 hosts). */
#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
void irrep_cutoff_polynomial_batch_neon    (size_t n, const double *r, double r_cut, int p, double *out);
void irrep_cutoff_polynomial_d_batch_neon  (size_t n, const double *r, double r_cut, int p, double *out);
void irrep_sph_harm_cart_all_batch_neon    (int l_max, size_t n, const double *r_hats, double *out);
#endif

#ifdef __cplusplus
}
#endif

#endif /* IRREP_INTERNAL_DISPATCH_H */
