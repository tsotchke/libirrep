/* SPDX-License-Identifier: MIT */
/* Downstream-compat golden-vector test.
 *
 * Runs five fixed (m_i, m_j, r̂) configurations through the UVW tensor-
 * product paths that reproduce the spin-torque-net's cartesian basis, and
 * asserts bit-exact agreement with the golden JSON under
 *   tests/test_downstream_compat/torque_net_tp_paths/vectors.json
 *
 * The JSON is authoritative; the test compares TP output (after the
 * documented √3 / √2 prefactor rescaling and the (y,z,x)→(x,y,z) layout
 * permutation) against the pre-computed reference cartesian values.
 * Residual threshold 1e-12 relative.
 *
 * When the downstream spin_based_neural_network vendors its own golden-
 * forces file (`expected_forces.json`), this test will extend to assert
 * bit-equality there too. Until then a `# SKIP:` line surfaces the gap in
 * TAP output without silently passing. */

#include "harness.h"
#include "test_downstream_compat/lattice_loader.h"
#include <irrep/multiset.h>
#include <irrep/tensor_product.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static int file_exists(const char *path) {
    struct stat st;
    return stat(path, &st) == 0 && S_ISREG(st.st_mode);
}

/* Minimal ad-hoc JSON reader for the fixed schema of vectors.json. Not a
 * general parser — consumes only the exact field order the generator emits
 * (name, m_i, m_j, r_hat, T1..T5). */

static const char *find_(const char *p, const char *needle) {
    return strstr(p, needle);
}

/* Parse a JSON number-triple `[a, b, c]` starting at or before `*pp`.
 * Advances *pp past the closing bracket on success. */
static int parse_vec3_(const char **pp, double out[3]) {
    const char *p = *pp;
    while (*p && *p != '[')
        ++p;
    if (*p != '[')
        return -1;
    ++p;
    char *end = NULL;
    for (int i = 0; i < 3; ++i) {
        out[i] = strtod(p, &end);
        if (end == p)
            return -1;
        p = end;
        while (*p == ' ' || *p == ',' || *p == '\t' || *p == '\n')
            ++p;
    }
    while (*p && *p != ']')
        ++p;
    if (*p != ']')
        return -1;
    *pp = p + 1;
    return 0;
}

/* Load one config from the vectors.json content; returns pointer advanced
 * past the config block, or NULL on parse failure. */
static const char *load_config_(const char *p, double m_i[3], double m_j[3], double r_hat[3],
                                double T1[3], double T2[3], double T3[3], double T4[3],
                                double T5[3]) {
    const char *q = p;
    q = find_(q, "\"m_i\":");
    if (!q)
        return NULL;
    q += strlen("\"m_i\":");
    if (parse_vec3_(&q, m_i) != 0)
        return NULL;
    q = find_(q, "\"m_j\":");
    if (!q)
        return NULL;
    q += strlen("\"m_j\":");
    if (parse_vec3_(&q, m_j) != 0)
        return NULL;
    q = find_(q, "\"r_hat\":");
    if (!q)
        return NULL;
    q += strlen("\"r_hat\":");
    if (parse_vec3_(&q, r_hat) != 0)
        return NULL;
    /* parse_vec3_ already scans forward to the next '[', so the exact
     * advance past the key doesn't matter — leave it at 1 so find_
     * starts past the current key and can't match it again. */
    q = find_(q, "\"T1_dot_proj\":");
    if (!q)
        return NULL;
    ++q;
    if (parse_vec3_(&q, T1) != 0)
        return NULL;
    q = find_(q, "\"T2_mj_cross_r\":");
    if (!q)
        return NULL;
    ++q;
    if (parse_vec3_(&q, T2) != 0)
        return NULL;
    q = find_(q, "\"T3_mi_cross_mj\":");
    if (!q)
        return NULL;
    ++q;
    if (parse_vec3_(&q, T3) != 0)
        return NULL;
    q = find_(q, "\"T4_mij_mi\":");
    if (!q)
        return NULL;
    ++q;
    if (parse_vec3_(&q, T4) != 0)
        return NULL;
    q = find_(q, "\"T5_mj\":");
    if (!q)
        return NULL;
    ++q;
    if (parse_vec3_(&q, T5) != 0)
        return NULL;
    return q;
}

/* Slurp whole file. */
static char *slurp_(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f)
        return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (n < 0) {
        fclose(f);
        return NULL;
    }
    char *buf = malloc((size_t)n + 1);
    if (!buf) {
        fclose(f);
        return NULL;
    }
    size_t got = fread(buf, 1, (size_t)n, f);
    fclose(f);
    buf[got] = '\0';
    return buf;
}

/* (x, y, z) → (y, z, x) real-SH l=1 permutation. */
static void cart_to_real_(const double c[3], double r[3]) {
    r[0] = c[1];
    r[1] = c[2];
    r[2] = c[0];
}
static void real_to_cart_(const double r[3], double c[3]) {
    c[0] = r[2];
    c[1] = r[0];
    c[2] = r[1];
}

/* Build the (1o ⊗ 1o → 0e + 1e) tensor-product descriptor once and reuse.
 * Returns a newly-allocated descriptor; caller frees with irrep_tp_free. */
static tp_descriptor_t *build_tp_(void) {
    irrep_multiset_t *a = irrep_multiset_parse("1x1o");
    irrep_multiset_t *b = irrep_multiset_parse("1x1o");
    irrep_multiset_t *c = irrep_multiset_parse("1x0e + 1x1e");
    int               np = irrep_tp_enumerate_paths(a, b, c, NULL, 0);
    int              *paths = malloc((size_t)np * 3 * sizeof(int));
    irrep_tp_enumerate_paths(a, b, c, paths, np);
    tp_descriptor_t *tp = irrep_tp_build_uvw(a, b, c, paths, np);
    free(paths);
    irrep_multiset_free(a);
    irrep_multiset_free(b);
    irrep_multiset_free(c);
    return tp;
}

/* Evaluate (a · b) via the 1o⊗1o→0e path, scaled by √3. */
static double tp_dot_(tp_descriptor_t *tp, const double *w, const double *aw,
                      const double cart_a[3], const double cart_b[3]) {
    double a_real[3], b_real[3];
    cart_to_real_(cart_a, a_real);
    cart_to_real_(cart_b, b_real);
    double out[4] = {0}; /* 0e (1) + 1e (3) */
    irrep_tp_apply_uvw(tp, w, a_real, b_real, out);
    (void)aw;
    return out[0] * sqrt(3.0);
}

/* Evaluate (a × b) via the 1o⊗1o→1e path, scaled by √2 and permuted. */
static void tp_cross_(tp_descriptor_t *tp, const double *w, const double cart_a[3],
                      const double cart_b[3], double cart_out[3]) {
    double a_real[3], b_real[3];
    cart_to_real_(cart_a, a_real);
    cart_to_real_(cart_b, b_real);
    double out[4] = {0};
    irrep_tp_apply_uvw(tp, w, a_real, b_real, out);
    double cross_real[3] = {
        out[1] * sqrt(2.0),
        out[2] * sqrt(2.0),
        out[3] * sqrt(2.0),
    };
    real_to_cart_(cross_real, cart_out);
}

static double max_abs_diff_(const double a[3], const double b[3]) {
    double m = 0.0;
    for (int i = 0; i < 3; ++i) {
        double d = fabs(a[i] - b[i]);
        if (d > m)
            m = d;
    }
    return m;
}

int main(void) {
    IRREP_TEST_START("downstream_compat");

    /* ---- lattice connectivity JSON (shared with benchmarks) ---- */
    const char *conn = "tests/test_downstream_compat/lattice_connectivity.json";
    IRREP_ASSERT(file_exists(conn));

    lattice_shape_t shape;
    IRREP_ASSERT(lattice_shape_load(conn, "4x4_periodic_64_edges", &shape) == 0);
    IRREP_ASSERT(shape.num_nodes == 16);
    IRREP_ASSERT(shape.num_edges == 64);
    IRREP_ASSERT(shape.edge_src != NULL);
    IRREP_ASSERT(shape.edge_dst != NULL);
    IRREP_ASSERT(shape.edge_vec != NULL);
    lattice_shape_free(&shape);

    /* ---- Golden-vector TP-path round-trip ---- */
    const char *vectors = "tests/test_downstream_compat/torque_net_tp_paths/vectors.json";
    IRREP_ASSERT(file_exists(vectors));

    char *text = slurp_(vectors);
    IRREP_ASSERT(text != NULL);

    tp_descriptor_t *tp = build_tp_();
    IRREP_ASSERT(tp != NULL);
    int     nw = irrep_tp_num_weights_uvw(tp);
    double *w = calloc((size_t)nw, sizeof(double));
    for (int i = 0; i < nw; ++i)
        w[i] = 1.0; /* identity weights */

    /* Find the first config block, then iterate. */
    const char *p = find_(text, "\"configs\":");
    IRREP_ASSERT(p != NULL);

    const double kTol = 1e-12;
    int          n_loaded = 0;
    for (;;) {
        const char *block = find_(p, "\"name\":");
        if (!block)
            break;
        /* Advance p past the name so the next iteration finds the next block. */
        p = block + 7;

        double m_i[3], m_j[3], r_hat[3];
        double T1_ref[3], T2_ref[3], T3_ref[3], T4_ref[3], T5_ref[3];
        if (load_config_(block, m_i, m_j, r_hat, T1_ref, T2_ref, T3_ref, T4_ref, T5_ref) == NULL)
            break;

        /* T1 = (m_j · r̂) · m_i  via (1o⊗1o→0e) × √3, then scalar × vector. */
        double mjr_tp = tp_dot_(tp, w, NULL, m_j, r_hat);
        double T1_got[3];
        for (int i = 0; i < 3; ++i)
            T1_got[i] = mjr_tp * m_i[i];
        IRREP_ASSERT(max_abs_diff_(T1_got, T1_ref) < kTol);

        /* T2 = m_j × r̂  via (1o⊗1o→1e) × √2 + permute. */
        double T2_got[3];
        tp_cross_(tp, w, m_j, r_hat, T2_got);
        IRREP_ASSERT(max_abs_diff_(T2_got, T2_ref) < kTol);

        /* T3 = m_i × m_j  via the same path. */
        double T3_got[3];
        tp_cross_(tp, w, m_i, m_j, T3_got);
        IRREP_ASSERT(max_abs_diff_(T3_got, T3_ref) < kTol);

        /* T4 = (m_i · m_j) · m_i — dot path × √3, then scale m_i. */
        double mij_tp = tp_dot_(tp, w, NULL, m_i, m_j);
        double T4_got[3];
        for (int i = 0; i < 3; ++i)
            T4_got[i] = mij_tp * m_i[i];
        IRREP_ASSERT(max_abs_diff_(T4_got, T4_ref) < kTol);

        /* T5 = m_j (passthrough; trivially exact). */
        IRREP_ASSERT(max_abs_diff_(m_j, T5_ref) < kTol);

        ++n_loaded;
    }
    IRREP_ASSERT(n_loaded == 5); /* all five configs exercised */

    free(w);
    irrep_tp_free(tp);
    free(text);

    /* ---- Forces golden gated on downstream ship ---- */
    const char *forces = "tests/test_downstream_compat/torque_net_tp_paths/expected_forces.json";
    if (!file_exists(forces)) {
        printf("# SKIP: %s — awaiting spin_based_neural_network _apply_forces golden values\n",
               forces);
    }

    return IRREP_TEST_END();
}
