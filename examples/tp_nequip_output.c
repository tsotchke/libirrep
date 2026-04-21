/* SPDX-License-Identifier: MIT */
/* Emits the output of `irrep_tp_apply_weighted` on a fixed NequIP-like
 * descriptor as a JSON array. Consumed by `scripts/bench_vs_e3nn.py`
 * to verify that libirrep's output agrees numerically with e3nn's
 * FullyConnectedTensorProduct on the same inputs.
 *
 * Not benched, not in `make bench` — only useful paired with the
 * Python script. Built by `make build/bin/tp_emit_output` or the
 * script does it automatically. */

#include <irrep/tensor_product.h>
#include <irrep/multiset.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    irrep_multiset_t *a = irrep_multiset_parse("1x0e + 1x1o + 1x2e");
    irrep_multiset_t *b = irrep_multiset_parse("1x0e + 1x1o + 1x2e + 1x3o");
    irrep_multiset_t *c = irrep_multiset_parse("1x0e + 1x1o + 1x2e");
    int paths[64];
    int np = irrep_tp_enumerate_paths(a, b, c, paths, 64);
    tp_descriptor_t *d = irrep_tp_build(a, b, c, paths, np);

    /* Fixed inputs matching bench_tensor_product.c and the e3nn side. */
    double A[9]  = { 1, 0.5, 0.3, -0.2, 0.1, 0.4, -0.1, 0.6, -0.3 };
    double B[16] = { 1, 0.2, -0.3, 0.4, 0.5, -0.1, 0.2, 0.3, -0.4,
                     0.1, 0.2, -0.3, 0.5, -0.6, 0.7, -0.1 };
    double C[9]  = { 0 };
    /* Unweighted (raw CG-only) apply — matches the mathematical
     * object that e3nn's o3.FullTensorProduct computes, without
     * e3nn's per-path normalisation. */
    irrep_tp_apply(d, A, B, C);

    /* Emit as JSON array on one line. */
    printf("[");
    for (int i = 0; i < 9; ++i) {
        printf("%.17g", C[i]);
        if (i + 1 < 9) printf(", ");
    }
    printf("]\n");

    /* Also emit the weight vector (e3nn's internal weights arrange
     * paths in a different order; the Python side multiplies through
     * its own weight layout to compare output, not internals). */

    irrep_tp_free(d);
    irrep_multiset_free(a);
    irrep_multiset_free(b);
    irrep_multiset_free(c);
    return 0;
}
