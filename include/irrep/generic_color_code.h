/* SPDX-License-Identifier: MIT */
/** @file generic_color_code.h
 *  @brief Face-list-driven 2D color code construction.
 *
 *  Subsumes any 2D color code on a 3-colourable trivalent lattice:
 *
 *    - Steane [[7, 1, 3]] (6.6.6 / hex, smallest, L=1)
 *    - [[19, 1, 5]] hex (6.6.6, L=2)
 *    - [[17, 1, 5]] square-octagon (4.8.8)
 *    - Triangular `[[37, 1, 7]]` (L=3) and beyond, by formula
 *
 *  ## Construction
 *
 *  A 2D color code is specified by:
 *    - `n_qubits`: number of vertices (data qubits) of the lattice.
 *    - `n_faces`: number of faces of the (3-colourable) tiling.
 *    - For each face `f`, a list of the qubits on its boundary, of
 *      length `face_size[f]`.
 *
 *  Each face hosts BOTH an X-stabilizer and a Z-stabilizer on its
 *  boundary qubits — the defining property of 2D color codes.
 *
 *  ## CSS orthogonality
 *
 *  CSS orthogonality `H_X · H_Z^T = 0 (mod 2)` reduces (since
 *  `H_X = H_Z` for color codes) to the requirement that every pair of
 *  face boundaries shares an *even* number of qubits. On a 3-colourable
 *  trivalent lattice this is automatic: two faces share either an edge
 *  (2 qubits) or no boundary (0 qubits).
 *
 *  This module verifies CSS orthogonality at construction time so
 *  caller errors in the face list are caught early.
 *
 *  ## Face colours (optional)
 *
 *  For decoding and analysis, the user may supply a colour `∈ {R, G, B}`
 *  for each face. If present, the global redundancy `Σ_R H_R + Σ_G H_G + Σ_B H_B = 0`
 *  (each vertex appears in exactly one face of each colour) provides
 *  additional structure.
 *
 *  ## Primary references
 *
 *  - Bombín-Martín-Delgado, *Topological quantum distillation*, PRL 97
 *    (2006) 180501.
 *  - Landahl-Anderson-Rice, *Fault-tolerant quantum computing with color
 *    codes*, arXiv:1108.5738 (2011).
 *  - Kubica, *The ABCs of the color code*, PhD thesis Caltech (2018).
 */
#ifndef IRREP_GENERIC_COLOR_CODE_H
#define IRREP_GENERIC_COLOR_CODE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <irrep/css_code.h>
#include <irrep/export.h>
#include <irrep/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Optional face colours for a 2D color code. */
typedef enum {
    IRREP_COLOR_NONE = 0,
    IRREP_COLOR_R    = 1,
    IRREP_COLOR_G    = 2,
    IRREP_COLOR_B    = 3,
} irrep_color_t;

/** @brief Specification of a 2D color-code lattice as a face list.
 *
 *  Caller owns all arrays. The struct does NOT take ownership; the build
 *  function reads from it. */
typedef struct {
    int n_qubits;          /**< Number of vertices / data qubits. */
    int n_faces;           /**< Number of faces. */
    const int *face_sizes; /**< face_sizes[f] = weight of face f, length n_faces. */
    const int * const *face_qubits; /**< face_qubits[f] = array of qubit indices (length face_sizes[f]). */
    const irrep_color_t *face_color; /**< Optional, length n_faces (may be NULL). */
} irrep_color_lattice_t;

/** @brief Build a 2D color code from its face list.
 *
 *  Allocates `out` as a CSS code with `m_X = m_Z = n_faces` stabilizers.
 *  Each face contributes one X-stab AND one Z-stab on the same boundary
 *  qubits (the defining color-code property).
 *
 *  Returns `IRREP_ERR_PRECONDITION` if any face references an out-of-range
 *  qubit or the resulting CSS-orthogonality check fails. */
IRREP_API irrep_status_t
irrep_generic_color_build(const irrep_color_lattice_t *lattice,
                          irrep_css_code_t *out);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_GENERIC_COLOR_CODE_H */
