/* SPDX-License-Identifier: MIT */
/* Logical-operator demo: surface code [[9,1,3]] and 2×2 toric code [[8,2,2]].
 *
 * Demonstrates the v1.5.0-era extensions to the QEC stack:
 * `irrep_surface_logical_X` / `_Z` and `irrep_toric_logical_X1` / `_Z1`
 * / `_X2` / `_Z2`. Prints the supports + verifies all expected
 * (anti-)commutation relations against the stabilizer group.
 *
 * Build / run:
 *   make examples
 *   ./build/bin/qec_logical_operators_demo
 */
#include <irrep/css_code.h>
#include <irrep/stabilizer_group.h>
#include <irrep/surface_code.h>
#include <irrep/toric_code.h>

#include <stdio.h>
#include <stdlib.h>

static void print_pauli(const char *label, const irrep_pauli_t *p) {
    printf("  %-6s  [w=%d]  ", label, irrep_pauli_weight(p));
    for (int q = 0; q < p->n; ++q) {
        irrep_pauli_letter_t L = irrep_pauli_get(p, q);
        char c = (L == IRREP_PAULI_LETTER_X) ? 'X'
                : (L == IRREP_PAULI_LETTER_Y) ? 'Y'
                : (L == IRREP_PAULI_LETTER_Z) ? 'Z' : '.';
        putchar(c);
    }
    putchar('\n');
}

static int surface_demo(void) {
    printf("Rotated surface code [[9, 1, 3]]:\n");
    irrep_surface_params_t p;
    irrep_surface_init(&p, 3);
    irrep_css_code_t cs;
    if (irrep_surface_build(&p, &cs) != IRREP_OK) return 1;
    irrep_pauli_t Lx, Lz;
    irrep_surface_logical_X(&p, &Lx);
    irrep_surface_logical_Z(&p, &Lz);
    print_pauli("X̄", &Lx);
    print_pauli("Z̄", &Lz);
    printf("  X̄ Z̄ inner product (1 = anti-commute) = %d\n",
           irrep_pauli_symp_inner(&Lx, &Lz));

    /* Verify both commute with every stabilizer. */
    irrep_stabilizer_group_t g;
    irrep_css_code_to_stabilizer_group(&cs, &g);
    int anti_x = 0, anti_z = 0;
    for (int i = 0; i < g.n_generators; ++i) {
        if (!irrep_pauli_commute(&g.gens[i], &Lx)) ++anti_x;
        if (!irrep_pauli_commute(&g.gens[i], &Lz)) ++anti_z;
    }
    printf("  Stab anti-commutes:  X̄ = %d,  Z̄ = %d  (both must be 0)\n",
           anti_x, anti_z);

    irrep_pauli_free(&Lx);
    irrep_pauli_free(&Lz);
    irrep_stabilizer_group_free(&g);
    irrep_css_code_free(&cs);
    return (anti_x == 0 && anti_z == 0 &&
            irrep_pauli_symp_inner(&Lx, &Lz) >= 0) ? 0 : 1;
}

static int toric_demo(int Lx_dim, int Ly_dim) {
    printf("\nToric code (Lx=%d, Ly=%d) — [[%d, 2, %d]]:\n",
           Lx_dim, Ly_dim, 2 * Lx_dim * Ly_dim,
           Lx_dim < Ly_dim ? Lx_dim : Ly_dim);
    irrep_toric_params_t p;
    irrep_toric_init(&p, Lx_dim, Ly_dim);

    irrep_pauli_t Lx1, Lz1, Lx2, Lz2;
    irrep_toric_logical_X1(&p, &Lx1);
    irrep_toric_logical_Z1(&p, &Lz1);
    irrep_toric_logical_X2(&p, &Lx2);
    irrep_toric_logical_Z2(&p, &Lz2);
    print_pauli("X̄_1", &Lx1);
    print_pauli("Z̄_1", &Lz1);
    print_pauli("X̄_2", &Lx2);
    print_pauli("Z̄_2", &Lz2);

    printf("\n  Commutation matrix (1 = anti-commute):\n");
    irrep_pauli_t *ops[4] = { &Lx1, &Lz1, &Lx2, &Lz2 };
    const char *names[4] = { "X̄_1", "Z̄_1", "X̄_2", "Z̄_2" };
    printf("           ");
    for (int j = 0; j < 4; ++j) printf("  %-5s", names[j]);
    printf("\n");
    for (int i = 0; i < 4; ++i) {
        printf("  %-7s ", names[i]);
        for (int j = 0; j < 4; ++j) {
            printf("    %d  ", irrep_pauli_symp_inner(ops[i], ops[j]));
        }
        printf("\n");
    }
    printf("  Expected: 1 at (X̄_1, Z̄_1) and (X̄_2, Z̄_2) pairs (anti-commute);\n");
    printf("            0 elsewhere (cross-logical pairs commute).\n");

    irrep_pauli_free(&Lx1);
    irrep_pauli_free(&Lz1);
    irrep_pauli_free(&Lx2);
    irrep_pauli_free(&Lz2);
    return 0;
}

int main(void) {
    int rc = 0;
    rc |= surface_demo();
    rc |= toric_demo(3, 3);
    return rc;
}
