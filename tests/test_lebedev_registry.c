/* SPDX-License-Identifier: MIT */
/* Runtime Lebedev-rule registry tests. The data tables at orders 9..41
 * are not bundled in the library tree (fetched by
 * `scripts/fetch_lebedev_tables.sh`). These tests exercise the registry
 * contract with a small synthetic rule that satisfies the validation
 * invariants — the higher-order Burkardt-parsed rules are exercised
 * end-to-end by `examples/register_lebedev.c`.
 */

#include "harness.h"
#include <irrep/quadrature.h>
#include <irrep/types.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    IRREP_TEST_START("lebedev_registry");

    /* Before any registration: orders 9..41 return 0 / false. */
    irrep_lebedev_clear_registry();
    IRREP_ASSERT(irrep_lebedev_size(9) == 0);
    IRREP_ASSERT(irrep_lebedev_size(11) == 0);
    IRREP_ASSERT(irrep_lebedev_size(41) == 0);
    double dummy[4];
    IRREP_ASSERT(irrep_lebedev_fill(9, dummy) == false);

    /* Hard-coded orders always resolve, pre and post any registry work. */
    IRREP_ASSERT(irrep_lebedev_size(3) == 6);
    IRREP_ASSERT(irrep_lebedev_size(5) == 14);
    IRREP_ASSERT(irrep_lebedev_size(7) == 26);

    /* Register a synthetic "rule" at order 9 — a 6-point axis orbit with
     * weight 1/6 per point. Polynomial exactness of this specific rule
     * is only degree 3 (it IS the Lebedev order-3 rule); we register it
     * at order-9 here solely to exercise the registry contract without
     * depending on the fetched data files. */
    double rule[6 * 4] = {
        +1, 0, 0, 1.0 / 6.0,
        -1, 0, 0, 1.0 / 6.0,
        0, +1, 0, 1.0 / 6.0,
        0, -1, 0, 1.0 / 6.0,
        0, 0, +1, 1.0 / 6.0,
        0, 0, -1, 1.0 / 6.0,
    };
    IRREP_ASSERT(irrep_lebedev_register_rule(9, 6, rule) == IRREP_OK);
    IRREP_ASSERT(irrep_lebedev_size(9) == 6);

    /* Fill reproduces what we registered. */
    double filled[6 * 4];
    IRREP_ASSERT(irrep_lebedev_fill(9, filled) == true);
    for (int k = 0; k < 6 * 4; ++k)
        IRREP_ASSERT(filled[k] == rule[k]);

    /* Replace-in-place: register a different rule at the same order. */
    double rule2[6 * 4];
    memcpy(rule2, rule, sizeof(rule));
    for (int k = 0; k < 6; ++k)
        rule2[4 * k + 3] = 1.0 / 6.0; /* same weights — keep it valid */
    rule2[0] = -1;
    rule2[4] = +1; /* swap axis-0 pair */
    IRREP_ASSERT(irrep_lebedev_register_rule(9, 6, rule2) == IRREP_OK);
    IRREP_ASSERT(irrep_lebedev_size(9) == 6);

    /* Hard-coded orders cannot be overridden. */
    IRREP_ASSERT(irrep_lebedev_register_rule(3, 6, rule) == IRREP_ERR_PRECONDITION);
    IRREP_ASSERT(irrep_lebedev_register_rule(5, 6, rule) == IRREP_ERR_PRECONDITION);
    IRREP_ASSERT(irrep_lebedev_register_rule(7, 6, rule) == IRREP_ERR_PRECONDITION);

    /* Validation: non-unit point → INVALID_ARG. */
    double bad_point[6 * 4];
    memcpy(bad_point, rule, sizeof(rule));
    bad_point[0] = 2.0; /* |r|² = 4 ≠ 1 */
    IRREP_ASSERT(irrep_lebedev_register_rule(11, 6, bad_point) == IRREP_ERR_INVALID_ARG);

    /* Validation: weights don't sum to 1 → INVALID_ARG. */
    double bad_weights[6 * 4];
    memcpy(bad_weights, rule, sizeof(rule));
    bad_weights[3] = 2.0; /* wrecks the weight sum */
    IRREP_ASSERT(irrep_lebedev_register_rule(11, 6, bad_weights) == IRREP_ERR_INVALID_ARG);

    /* Validation: even order / negative order / NULL → INVALID_ARG. */
    IRREP_ASSERT(irrep_lebedev_register_rule(8, 6, rule) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_lebedev_register_rule(-1, 6, rule) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_lebedev_register_rule(11, 0, rule) == IRREP_ERR_INVALID_ARG);
    IRREP_ASSERT(irrep_lebedev_register_rule(11, 6, NULL) == IRREP_ERR_INVALID_ARG);

    /* Clear removes registered orders; hard-coded remain. */
    irrep_lebedev_clear_registry();
    IRREP_ASSERT(irrep_lebedev_size(9) == 0);
    IRREP_ASSERT(irrep_lebedev_size(3) == 6);

    return IRREP_TEST_END();
}
