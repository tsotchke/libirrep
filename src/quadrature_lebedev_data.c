/* SPDX-License-Identifier: MIT */
/* Lebedev quadrature tables (V.I. Lebedev & D.N. Laikov, 1999).
 *
 * Tables are public-domain; the original Fortran-source node & weight values
 * are imported in M8. For M1 we ship the file empty of data so `make` links
 * but `irrep_lebedev_fill()` returns false for every order. */

/* Typedef sentinel — makes this a valid ISO C translation unit without
 * producing a symbol or a mutable static. Replaced by real tables in M8. */
typedef int irrep_lebedev_data_placeholder_t;
