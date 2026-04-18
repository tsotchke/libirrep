/* SPDX-License-Identifier: MIT */
/* M6: O(3) parity helpers and path-list filter. */

#include <irrep/parity.h>

#define UNUSED(x) ((void)(x))

int irrep_parity(irrep_label_t label) {
    return label.parity == IRREP_ODD ? -1 : +1;
}

int irrep_parity_product(irrep_label_t a, irrep_label_t b) {
    return irrep_parity(a) * irrep_parity(b);
}

/* Compact `paths` in place, keeping only triplets (i_a, i_b, i_c) whose
 * parity_a · parity_b equals parity_c. Returns the new count. */
int irrep_parity_filter_paths(const irrep_multiset_t *a,
                              const irrep_multiset_t *b,
                              const irrep_multiset_t *c,
                              int *paths, int num_paths) {
    if (!a || !b || !c || !paths || num_paths <= 0) return 0;
    int write = 0;
    for (int r = 0; r < num_paths; ++r) {
        int ia = paths[r * 3 + 0];
        int ib = paths[r * 3 + 1];
        int ic = paths[r * 3 + 2];
        if (ia < 0 || ia >= a->num_terms) continue;
        if (ib < 0 || ib >= b->num_terms) continue;
        if (ic < 0 || ic >= c->num_terms) continue;

        int pa = irrep_parity(a->labels[ia]);
        int pb = irrep_parity(b->labels[ib]);
        int pc = irrep_parity(c->labels[ic]);
        if (pa * pb != pc) continue;

        paths[write * 3 + 0] = ia;
        paths[write * 3 + 1] = ib;
        paths[write * 3 + 2] = ic;
        write++;
    }
    return write;
}
