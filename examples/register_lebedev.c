/* SPDX-License-Identifier: MIT */
/* Parse Lebedev-Laikov 1999 quadrature tables (fetched via
 * `scripts/fetch_lebedev_tables.sh`) and register them at runtime so
 * `irrep_lebedev_size` / `irrep_lebedev_fill` return the high-order
 * rules without the library bundling any data.
 *
 * File format (per-line, Burkardt's convention):
 *
 *     phi_deg   theta_deg   weight
 *
 * `phi` is the azimuth (typically in [−180°, 360°]; sign-wrapped to
 * [0°, 360°) in-parser), `theta` is the polar angle from the +z axis
 * in [0°, 180°]. Converts to cartesian
 * `(x, y, z) = (sin θ cos φ, sin θ sin φ, cos θ)`. Weights as supplied
 * (summing to 1 on the unit sphere).
 *
 *   make examples
 *   scripts/fetch_lebedev_tables.sh
 *   ./build/bin/register_lebedev data/lebedev
 *
 * Exits non-zero on parse failure; prints a summary on success. */

#include <irrep/quadrature.h>
#include <irrep/types.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const int orders[] = {9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41};
#define N_ORDERS (sizeof(orders) / sizeof(orders[0]))

static int parse_and_register(const char *dir, int order) {
    char path[512];
    snprintf(path, sizeof(path), "%s/lebedev_%d.txt", dir, order);
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "lebedev order %d: cannot open %s\n", order, path);
        return -1;
    }
    size_t  cap   = 64;
    size_t  n     = 0;
    double *buf   = malloc(cap * 4 * sizeof(double));
    double  line_phi, line_theta, line_w;
    while (fscanf(f, "%lf %lf %lf", &line_phi, &line_theta, &line_w) == 3) {
        if (n + 1 > cap) {
            cap *= 2;
            buf = realloc(buf, cap * 4 * sizeof(double));
        }
        double th = line_theta * (M_PI / 180.0);
        double ph = line_phi   * (M_PI / 180.0);
        buf[4 * n + 0] = sin(th) * cos(ph);
        buf[4 * n + 1] = sin(th) * sin(ph);
        buf[4 * n + 2] = cos(th);
        buf[4 * n + 3] = line_w;
        ++n;
    }
    fclose(f);

    irrep_status_t rc = irrep_lebedev_register_rule(order, (int)n, buf);
    free(buf);
    if (rc != IRREP_OK) {
        fprintf(stderr, "lebedev order %d: register_rule failed (rc=%d)\n", order, rc);
        return -1;
    }
    printf("  order %2d: %4zu points registered\n", order, n);
    return 0;
}

int main(int argc, char **argv) {
    const char *dir = (argc > 1) ? argv[1] : "data/lebedev";
    printf("Registering Lebedev-Laikov 1999 tables from %s/:\n", dir);
    int ok = 0, fail = 0;
    for (size_t i = 0; i < N_ORDERS; ++i) {
        if (parse_and_register(dir, orders[i]) == 0)
            ++ok;
        else
            ++fail;
    }
    printf("\n%d / %zu orders registered. Library now serves\n", ok, N_ORDERS);
    printf("`irrep_lebedev_size` and `irrep_lebedev_fill` at each of them.\n");
    return fail > 0 ? 1 : 0;
}
