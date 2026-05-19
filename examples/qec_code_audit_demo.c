/* SPDX-License-Identifier: MIT */
/* QEC code audit demo: walks the library's named CSS codes, prints
 * the auto-derived [[n, k, d]] tuple via runtime F₂-rank and brute-
 * force distance primitives.
 *
 * # What this demonstrates
 *
 * Concrete runtime audit of every named QEC code shipped by libirrep,
 * using only the auto-discovery primitives (no hand-coded analytical
 * shortcuts). Output for each code:
 *
 *   - n_qubits        from the CSS struct
 *   - m_X, m_Z        check-matrix row counts
 *   - k               via irrep_css_code_logical_qubits
 *                     (= n - rank(H_X) - rank(H_Z))
 *   - d_X, d_Z        via irrep_css_code_compute_logical_X / _Z
 *                     (auto-discovers minimum-weight side logicals)
 *   - d = min(d_X, d_Z)  the CSS code distance
 *
 * # Output format (one row per code)
 *
 *   code-name    [[n, k, d]]   n_qubits  m_X  m_Z  k  d_X  d_Z  d
 *
 * Codes shown: Steane [[7,1,3]], hex [[19,1,5]], [[37,1,7]] from the
 * parametric hex generator, [[15,7,3]] Hamming-CSS, [[15,1,3]]
 * Reed-Muller (triorthogonal — d_X ≠ d_Z), [[8,3,2]] cubic 3D color
 * code, [[13,1,3]] HGP(rep-3, rep-3), [[25,1,4]] HGP(rep-4, rep-4),
 * Steane⊗Steane [[49,1,9]].
 *
 * Each row's d should match the published distance for that code;
 * runtime mismatches would expose a code-construction bug. */
#include <irrep/color_code.h>
#include <irrep/color_code_3d.h>
#include <irrep/color_codes_2d.h>
#include <irrep/concatenated_code.h>
#include <irrep/css_code.h>
#include <irrep/hypergraph_product.h>
#include <irrep/stabilizer_group.h>

#include <stdio.h>
#include <stdlib.h>

static void
audit_code(const char *label, irrep_css_code_t *cs, int max_weight)
{
    if (cs == NULL) return;
    int k = irrep_css_code_logical_qubits(cs);
    irrep_pauli_t Lx, Lz;
    int d_X = irrep_css_code_compute_logical_X(cs, max_weight, &Lx);
    int d_Z = irrep_css_code_compute_logical_Z(cs, max_weight, &Lz);
    int d = (d_X < d_Z) ? d_X : d_Z;
    /* If either side hit max_weight+1, the actual d_X/d_Z exceeds the
     * search window — flag with a '*' to distinguish from a confirmed
     * value. */
    int d_X_bounded = (d_X <= max_weight);
    int d_Z_bounded = (d_Z <= max_weight);
    printf("  %-30s  [[%2d, %d, %d%s]]   m_X=%2d m_Z=%2d  d_X=%d%s d_Z=%d%s\n",
           label, cs->n, k, d,
           (!d_X_bounded || !d_Z_bounded) ? "*" : "",
           cs->H_X.n_rows, cs->H_Z.n_rows,
           d_X, d_X_bounded ? "" : "+",
           d_Z, d_Z_bounded ? "" : "+");
    if (d_X <= max_weight) irrep_pauli_free(&Lx);
    if (d_Z <= max_weight) irrep_pauli_free(&Lz);
}

int main(void) {
    printf("libirrep QEC code audit (auto-derived [[n, k, d]] via F₂-rank + brute distance)\n");
    printf("================================================================================\n");

    /* Steane. */
    { irrep_css_code_t cs;
      if (irrep_color_steane(&cs) == IRREP_OK) {
          audit_code("Steane [[7, 1, 3]]", &cs, 4);
          irrep_css_code_free(&cs);
      } }
    /* [[15, 7, 3]] Hamming-CSS. */
    { irrep_css_code_t cs;
      if (irrep_color_hamming_15_7_3(&cs) == IRREP_OK) {
          audit_code("[[15, 7, 3]] Hamming-CSS", &cs, 4);
          irrep_css_code_free(&cs);
      } }
    /* [[19, 1, 5]] hex. */
    { irrep_css_code_t cs;
      if (irrep_color_hex_19_1_5(&cs) == IRREP_OK) {
          audit_code("[[19, 1, 5]] triangular hex", &cs, 6);
          irrep_css_code_free(&cs);
      } }
    /* [[37, 1, 7]] from parametric hex generator. */
    { irrep_css_code_t cs;
      if (irrep_color_hex_triangular_build(7, &cs) == IRREP_OK) {
          /* d_X / d_Z for n=37 at weight 7 takes a few seconds; cap
           * search at weight 5 (= proves d_X, d_Z ≥ 6 → d ≥ 6). For
           * concrete d = 7 verification, use a larger max_weight. */
          audit_code("[[37, 1, 7]] hex(d=7)", &cs, 5);
          irrep_css_code_free(&cs);
      } }
    /* [[8, 3, 2]] cubic 3D color code. */
    { irrep_css_code_t cs;
      if (irrep_color_3d_cube_8_3_2(&cs) == IRREP_OK) {
          audit_code("[[8, 3, 2]] 3D cubic color", &cs, 3);
          irrep_css_code_free(&cs);
      } }
    /* [[15, 1, 3]] Reed-Muller — triorthogonal, d_X ≠ d_Z. */
    { irrep_css_code_t cs;
      if (irrep_color_3d_rm_15_1_3(&cs) == IRREP_OK) {
          audit_code("[[15, 1, 3]] RM (triorth)", &cs, 8);
          irrep_css_code_free(&cs);
      } }
    /* [[13, 1, 3]] HGP. */
    { irrep_css_code_t cs;
      if (irrep_hgp_repetition_3_13_1_3(&cs) == IRREP_OK) {
          audit_code("[[13, 1, 3]] HGP(rep-3)", &cs, 4);
          irrep_css_code_free(&cs);
      } }
    /* [[25, 1, 4]] HGP. */
    { irrep_css_code_t cs;
      if (irrep_hgp_repetition_4_25_1_4(&cs) == IRREP_OK) {
          audit_code("[[25, 1, 4]] HGP(rep-4)", &cs, 5);
          irrep_css_code_free(&cs);
      } }

    printf("================================================================================\n");
    printf("All distances are runtime-computed via F_2-rank + brute-force enumeration.\n");
    printf("'*' / '+' annotations indicate the search window was capped below the true distance\n");
    printf("(e.g., [[37, 1, 7]] needs max_weight >= 7 to confirm d = 7).\n");
    return 0;
}
