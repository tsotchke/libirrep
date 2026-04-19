/* SPDX-License-Identifier: MIT */
/* Golden-vector generator for the torque-net TP-path compatibility test.
 *
 * Generates `vectors.json` carrying five fixed (m_i, m_j, r̂) triples and
 * their expected outputs under the five torque-net basis terms expressed
 * via `irrep_tp_apply_uvw`:
 *
 *   T1  =  (m_j · r̂) · m_i   — (1o,1o→0e) then (0e,1o→1o),   × √3 then ×1
 *   T2  =  m_j × r̂ (axial)   — (1o,1o→1e),                   × √2
 *   T3  =  m_i × m_j (axial) — (1o,1o→1e),                   × √2
 *   T4  =  (m_i · m_j) · m_i — (1o,1o→0e) then (0e,1o→1o),   × √3 then ×1
 *   T5  =  m_j               — passthrough
 *
 * The cartesian ↔ real-SH basis permutation is (x, y, z) ↔ (y, z, x).
 * All sign conventions verified empirically in
 * `examples/torque_net_tp_paths.c` against bit-exact cartesian
 * reference values (residual ≲ 2×10⁻¹⁶).
 *
 * Build (from the repo root):
 *   cc -Iinclude \
 *      tests/test_downstream_compat/torque_net_tp_paths/generate_golden.c \
 *      build/lib/liblibirrep.a -lm -o /tmp/gen_golden
 *   /tmp/gen_golden tests/test_downstream_compat/torque_net_tp_paths/vectors.json
 *
 * The produced JSON is consumed by `tests/test_downstream_compat.c`; the
 * downstream `spin_based_neural_network` cross-checks the same golden
 * values from its side to pin conventions across repositories. */

#include <irrep/irrep.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Cartesian → real-SH l=1 permutation (x, y, z) → (y, z, x). */
static void cart_to_real_l1(const double c[3], double r[3]) {
    r[0] = c[1]; r[1] = c[2]; r[2] = c[0];
}
static void real_l1_to_cart(const double r[3], double c[3]) {
    c[0] = r[2]; c[1] = r[0]; c[2] = r[1];
}

static double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static void cross3(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

typedef struct {
    const char *name;
    double m_i[3];
    double m_j[3];
    double r_hat[3];   /* unit vector */
} config_t;

static const config_t k_configs[] = {
    { "axial_z",
      { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 0.0 } },
    { "perpendicular",
      { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } },
    { "generic_1",
      { 0.7, -0.3, 0.5 }, { -0.2, 0.4, 0.8 }, { 0.5, 0.5, 0.70710678118654752 } },
    { "generic_2",
      { 0.6, 0.8, 0.0 }, { 0.0, 0.6, 0.8 }, { 0.8, 0.0, 0.6 } },
    { "antiparallel",
      { 1.0, 0.0, 0.0 }, { -1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 } },
};
#define N_CONFIGS ((int)(sizeof(k_configs) / sizeof(k_configs[0])))

/* Compute the five basis terms analytically in cartesian, the authoritative
 * reference: the TP-path computation must reproduce these within ~10⁻¹². */
static void eval_cartesian_reference(const config_t *c,
                                     double T1[3], double T2[3], double T3[3],
                                     double T4[3], double T5[3]) {
    double mjr  = dot3(c->m_j, c->r_hat);            /* m_j · r̂ */
    double mij  = dot3(c->m_i, c->m_j);              /* m_i · m_j */
    for (int i = 0; i < 3; ++i) T1[i] = mjr * c->m_i[i];
    cross3(c->m_j, c->r_hat, T2);
    cross3(c->m_i, c->m_j,   T3);
    for (int i = 0; i < 3; ++i) T4[i] = mij * c->m_i[i];
    for (int i = 0; i < 3; ++i) T5[i] = c->m_j[i];
}

int main(int argc, char **argv) {
    const char *out_path = (argc > 1)
        ? argv[1]
        : "tests/test_downstream_compat/torque_net_tp_paths/vectors.json";

    FILE *f = fopen(out_path, "w");
    if (!f) { perror(out_path); return 1; }

    const double sqrt2 = sqrt(2.0);
    const double sqrt3 = sqrt(3.0);

    fprintf(f, "{\n");
    fprintf(f, "  \"_version\": \"1\",\n");
    fprintf(f, "  \"_note\": \"Torque-net TP-path golden vectors. Five fixed (m_i, m_j, r̂) triples and the five basis terms T1..T5 expressed via irrep_tp_apply_uvw. Conventions: 1o⊗1o→0e gives (1/√3)·(a·b); 1o⊗1o→1e gives (+1/√2)·(a×b) in the (y,z,x) real-SH layout. Multiplying the 0e TP output by √3 recovers the bare dot; multiplying the 1e TP output by √2 and permuting (y,z,x)→(x,y,z) recovers the bare cross product. Reference (cartesian) values are authoritative; the TP-path values after rescaling and permutation must match within 1e-12.\",\n");
    fprintf(f, "  \"generator\": \"tests/test_downstream_compat/torque_net_tp_paths/generate_golden.c\",\n");
    fprintf(f, "  \"libirrep_version\": \"%s\",\n", irrep_version_string());
    fprintf(f, "  \"libirrep_abi_hash\": \"%s\",\n", irrep_abi_hash());
    fprintf(f, "  \"configs\": [\n");

    for (int c = 0; c < N_CONFIGS; ++c) {
        const config_t *cfg = &k_configs[c];
        double T1[3], T2[3], T3[3], T4[3], T5[3];
        eval_cartesian_reference(cfg, T1, T2, T3, T4, T5);

        fprintf(f, "    {\n");
        fprintf(f, "      \"name\": \"%s\",\n", cfg->name);
        fprintf(f, "      \"m_i\":   [%.17g, %.17g, %.17g],\n",
                cfg->m_i[0], cfg->m_i[1], cfg->m_i[2]);
        fprintf(f, "      \"m_j\":   [%.17g, %.17g, %.17g],\n",
                cfg->m_j[0], cfg->m_j[1], cfg->m_j[2]);
        fprintf(f, "      \"r_hat\": [%.17g, %.17g, %.17g],\n",
                cfg->r_hat[0], cfg->r_hat[1], cfg->r_hat[2]);
        fprintf(f, "      \"T1_dot_proj\":    [%.17g, %.17g, %.17g],\n", T1[0], T1[1], T1[2]);
        fprintf(f, "      \"T2_mj_cross_r\":  [%.17g, %.17g, %.17g],\n", T2[0], T2[1], T2[2]);
        fprintf(f, "      \"T3_mi_cross_mj\": [%.17g, %.17g, %.17g],\n", T3[0], T3[1], T3[2]);
        fprintf(f, "      \"T4_mij_mi\":      [%.17g, %.17g, %.17g],\n", T4[0], T4[1], T4[2]);
        fprintf(f, "      \"T5_mj\":          [%.17g, %.17g, %.17g]\n",  T5[0], T5[1], T5[2]);
        fprintf(f, "    }%s\n", (c + 1 < N_CONFIGS) ? "," : "");
    }

    fprintf(f, "  ],\n");
    fprintf(f, "  \"constants\": {\n");
    fprintf(f, "    \"sqrt2\": %.17g,\n", sqrt2);
    fprintf(f, "    \"sqrt3\": %.17g,\n", sqrt3);
    fprintf(f, "    \"scale_0e_to_dot\":   %.17g,\n", sqrt3);
    fprintf(f, "    \"scale_1e_to_cross\": %.17g\n",  sqrt2);
    fprintf(f, "  }\n");
    fprintf(f, "}\n");

    fclose(f);
    fprintf(stderr, "wrote %s with %d configs\n", out_path, N_CONFIGS);
    return 0;
}
