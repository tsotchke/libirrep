/* SPDX-License-Identifier: MIT */
/* One NequIP-style message-passing step:
 *
 *   m_ij = W(r_ij) · ( h_j ⊗ Y(r̂_ij) ) · φ_cut(r_ij)
 *
 * where Y(r̂) are cartesian real spherical harmonics, φ_cut is a NequIP
 * polynomial cutoff, and the tensor product is weighted path-by-path.
 *
 * Full NequIP would multiply this by a learned radial function before the
 * tensor product; we bake the cutoff plus a single Bessel radial into one
 * per-path scalar weight to keep the example focused. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <irrep/multiset.h>
#include <irrep/spherical_harmonics.h>
#include <irrep/radial.h>
#include <irrep/tensor_product.h>

int main(void) {
    const double r_cut = 3.0;
    const int    p_cut = 6;

    /* Neighbour feature h_j: scalar + vector = 1x0e + 1x1o (4 real values). */
    irrep_multiset_t *h_mset = irrep_multiset_parse("1x0e + 1x1o");
    double h_j[4] = { 0.75,   0.1, 0.2, -0.3 };

    /* Edge geometry: vector from node i to node j. */
    double r_vec[3] = { 0.4, -0.6, 1.2 };
    double r_mag    = sqrt(r_vec[0]*r_vec[0] + r_vec[1]*r_vec[1] + r_vec[2]*r_vec[2]);
    double r_hat[3] = { r_vec[0] / r_mag, r_vec[1] / r_mag, r_vec[2] / r_mag };

    /* Y(r̂) up to l=2: 1 + 3 + 5 = 9 real values. */
    int sh_lmax = 2;
    irrep_multiset_t *y_mset = irrep_multiset_parse("1x0e + 1x1o + 1x2e");
    double y[9];
    irrep_sph_harm_cart_all(sh_lmax, y, r_hat);

    /* Output features: same irrep structure as h_j, promoted by the TP. */
    /* Output = (parity-natural couplings only) of h_j ⊗ Y(r̂). */
    irrep_multiset_t *m_mset = irrep_multiset_parse("1x0e + 1x1o + 1x2e");

    /* Enumerate valid paths. */
    int max_paths = 32;
    int *paths = malloc((size_t)max_paths * 3 * sizeof(int));
    int num_paths = irrep_tp_enumerate_paths(h_mset, y_mset, m_mset,
                                              paths, max_paths);

    tp_descriptor_t *desc = irrep_tp_build(h_mset, y_mset, m_mset,
                                            paths, num_paths);

    /* Per-path weight: Bessel(n=1) radial × polynomial cutoff. */
    double weights[32];
    for (int k = 0; k < num_paths; ++k) {
        weights[k] = irrep_rbf_bessel(1, r_mag, r_cut)
                   * irrep_cutoff_polynomial(r_mag, r_cut, p_cut);
    }

    double m_out[10] = { 0 };
    irrep_tp_apply_weighted(desc, weights, h_j, y, m_out);

    printf("nequip_message\n");
    printf("  r_mag        = %.6f,   r_cut = %.2f\n", r_mag, r_cut);
    printf("  Bessel(1, r) = %.6f,   cutoff_poly(r, p=6) = %.6f\n",
           irrep_rbf_bessel(1, r_mag, r_cut),
           irrep_cutoff_polynomial(r_mag, r_cut, p_cut));
    printf("  num_paths    = %d\n", num_paths);
    printf("  output m_ij  (dim = %d):\n", m_mset->total_dim);
    int dim = m_mset->total_dim;
    for (int i = 0; i < dim; ++i) printf("    m[%2d] = % .8f\n", i, m_out[i]);

    irrep_tp_free(desc);
    irrep_multiset_free(h_mset);
    irrep_multiset_free(y_mset);
    irrep_multiset_free(m_mset);
    free(paths);
    return 0;
}
