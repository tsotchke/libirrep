/* SPDX-License-Identifier: MIT */
/* Focused reader for lattice_connectivity.json — knows only the exact
 * schema the generator emits, nothing else. See lattice_loader.h. */

#include "lattice_loader.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Diagnostic emitter. Mirrors fprintf(stderr, …) style so CI log output
 * pinpoints which field of which shape failed to parse. Avoids pulling in
 * the library's irrep_set_error_ hook — this loader is a test-tree utility
 * and shouldn't depend on the public error channel. */
static int load_fail_(const char *json_path, const char *shape_name,
                      const char *detail) {
    fprintf(stderr, "lattice_loader: %s: shape '%s': %s\n",
            json_path, shape_name, detail);
    return -1;
}
#define FAIL(msg) return load_fail_(json_path, shape_name, (msg))

/* --- whole-file slurp ---------------------------------------------------- */
static char *slurp_file_(const char *path, size_t *len_out) {
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (n < 0) { fclose(f); return NULL; }
    char *buf = malloc((size_t)n + 1);
    if (!buf) { fclose(f); return NULL; }
    size_t got = fread(buf, 1, (size_t)n, f);
    fclose(f);
    buf[got] = '\0';
    if (len_out) *len_out = got;
    return buf;
}

/* Find the first occurrence of `needle` starting at `start`; return pointer
 * to it within the buffer or NULL. */
static const char *find_after_(const char *haystack, const char *needle) {
    return strstr(haystack, needle);
}

/* Skip whitespace. */
static const char *skip_ws_(const char *p) {
    while (*p && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')) ++p;
    return p;
}

/* Parse a non-negative integer at `p`, advance `p`. */
static int parse_int_(const char **pp) {
    const char *p = skip_ws_(*pp);
    int neg = 0;
    if (*p == '-') { neg = 1; ++p; }
    int v = 0;
    while (isdigit((unsigned char)*p)) { v = v * 10 + (*p - '0'); ++p; }
    *pp = p;
    return neg ? -v : v;
}

/* Parse a double at `p`, advance `p`. */
static double parse_dbl_(const char **pp) {
    const char *p = skip_ws_(*pp);
    char *end = NULL;
    double v = strtod(p, &end);
    *pp = end ? end : p;
    return v;
}

/* Expect a character, skipping leading whitespace. */
static int expect_(const char **pp, char c) {
    const char *p = skip_ws_(*pp);
    if (*p != c) return -1;
    *pp = p + 1;
    return 0;
}

/* Parse an array of ints into a freshly-allocated buffer sized @p expected.
 * `*pp` must point at or before the opening `[`. */
static int parse_int_array_(const char **pp, int expected, int **out) {
    if (expect_(pp, '[') != 0) return -1;
    int *buf = malloc((size_t)expected * sizeof(int));
    if (!buf) return -1;
    for (int i = 0; i < expected; ++i) {
        buf[i] = parse_int_(pp);
        if (i + 1 < expected) {
            if (expect_(pp, ',') != 0) { free(buf); return -1; }
        }
    }
    if (expect_(pp, ']') != 0) { free(buf); return -1; }
    *out = buf;
    return 0;
}

/* Parse an array of 3-tuples of doubles. */
static int parse_vec3_array_(const char **pp, int expected, double **out) {
    if (expect_(pp, '[') != 0) return -1;
    double *buf = malloc((size_t)expected * 3 * sizeof(double));
    if (!buf) return -1;
    for (int i = 0; i < expected; ++i) {
        if (expect_(pp, '[') != 0) { free(buf); return -1; }
        buf[3*i + 0] = parse_dbl_(pp);
        if (expect_(pp, ',') != 0) { free(buf); return -1; }
        buf[3*i + 1] = parse_dbl_(pp);
        if (expect_(pp, ',') != 0) { free(buf); return -1; }
        buf[3*i + 2] = parse_dbl_(pp);
        if (expect_(pp, ']') != 0) { free(buf); return -1; }
        if (i + 1 < expected) {
            if (expect_(pp, ',') != 0) { free(buf); return -1; }
        }
    }
    if (expect_(pp, ']') != 0) { free(buf); return -1; }
    *out = buf;
    return 0;
}

int lattice_shape_load(const char *json_path, const char *shape_name,
                       lattice_shape_t *out) {
    if (!json_path || !shape_name || !out) return -1;
    memset(out, 0, sizeof(*out));

    size_t flen = 0;
    char *text = slurp_file_(json_path, &flen);
    if (!text) {
        fprintf(stderr, "lattice_loader: %s: cannot open for reading\n", json_path);
        return -1;
    }

    /* Locate the shape block by its quoted key. */
    char key[256];
    int kn = snprintf(key, sizeof(key), "\"%s\"", shape_name);
    if (kn < 0 || kn >= (int)sizeof(key)) { free(text); FAIL("shape name too long"); }
    const char *p = find_after_(text, key);
    if (!p) { free(text); FAIL("shape key not found in JSON"); }

    /* Parse fields strictly in the order the generator emits. */
    const char *q;
    int num_nodes = 0, num_edges = 0;

    q = find_after_(p, "\"num_nodes\":");
    if (!q) { free(text); FAIL("missing 'num_nodes' field"); }
    q += strlen("\"num_nodes\":");
    num_nodes = parse_int_(&q);
    if (num_nodes < 0) { free(text); FAIL("num_nodes must be non-negative"); }

    q = find_after_(q, "\"num_edges\":");
    if (!q) { free(text); FAIL("missing 'num_edges' field"); }
    q += strlen("\"num_edges\":");
    num_edges = parse_int_(&q);
    /* Explicit validation — without this a malicious or malformed JSON with
     * e.g. `"num_edges": -1` would reach malloc((size_t)-1 ...) and rely on
     * the allocator to fail the request. Safer to bail with a diagnostic. */
    if (num_edges < 0) { free(text); FAIL("num_edges must be non-negative"); }

    q = find_after_(q, "\"edge_src\":");
    if (!q) { free(text); FAIL("missing 'edge_src' field"); }
    q += strlen("\"edge_src\":");
    int *src = NULL;
    if (parse_int_array_(&q, num_edges, &src) != 0) {
        free(text);
        FAIL("failed to parse 'edge_src' array (check num_edges agrees with length)");
    }

    q = find_after_(q, "\"edge_dst\":");
    if (!q) { free(src); free(text); FAIL("missing 'edge_dst' field"); }
    q += strlen("\"edge_dst\":");
    int *dst = NULL;
    if (parse_int_array_(&q, num_edges, &dst) != 0) {
        free(src); free(text);
        FAIL("failed to parse 'edge_dst' array");
    }

    q = find_after_(q, "\"edge_vec\":");
    if (!q) { free(src); free(dst); free(text); FAIL("missing 'edge_vec' field"); }
    q += strlen("\"edge_vec\":");
    double *vec = NULL;
    if (parse_vec3_array_(&q, num_edges, &vec) != 0) {
        free(src); free(dst); free(text);
        FAIL("failed to parse 'edge_vec' array of [x, y, z] triples");
    }

    free(text);
    out->num_nodes = num_nodes;
    out->num_edges = num_edges;
    out->edge_src  = src;
    out->edge_dst  = dst;
    out->edge_vec  = vec;
    return 0;
}

void lattice_shape_free(lattice_shape_t *shape) {
    if (!shape) return;
    free(shape->edge_src);
    free(shape->edge_dst);
    free(shape->edge_vec);
    memset(shape, 0, sizeof(*shape));
}
