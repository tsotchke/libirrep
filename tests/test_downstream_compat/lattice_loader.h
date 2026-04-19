/* SPDX-License-Identifier: MIT */
/* Minimal reader for tests/test_downstream_compat/lattice_connectivity.json.
 *
 * Loads a single named shape ("4x4_periodic_64_edges", …) into caller-owned
 * buffers. The JSON is the shared source of truth with
 * `spin_based_neural_network`'s tree; the generator
 * (generate_lattice_connectivity.c) lives downstream and emits both copies
 * from the same `torque_net_build_grid` enumeration, so any edge-ordering
 * drift surfaces in both trees simultaneously.
 *
 * This reader is deliberately minimal — it understands exactly the schema
 * the generator produces and nothing else. Not a general JSON parser.
 */
#ifndef LIBIRREP_LATTICE_LOADER_H
#define LIBIRREP_LATTICE_LOADER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int     num_nodes;
    int     num_edges;
    int    *edge_src;       /* num_edges  */
    int    *edge_dst;       /* num_edges  */
    double *edge_vec;       /* num_edges × 3 */
} lattice_shape_t;

/* Load named shape from @p json_path into @p out. Caller frees via
 * `lattice_shape_free`. Returns 0 on success, non-zero on error. */
int  lattice_shape_load(const char *json_path, const char *shape_name,
                        lattice_shape_t *out);

void lattice_shape_free(lattice_shape_t *shape);

#ifdef __cplusplus
}
#endif

#endif /* LIBIRREP_LATTICE_LOADER_H */
