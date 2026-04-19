/* SPDX-License-Identifier: MIT */
/** @file irrep.h
 *  @brief Umbrella header — includes every public module. Use this when you
 *         don't care about a narrow dependency graph.
 *
 *  Individual headers are also self-contained and can be included directly
 *  if you'd rather keep your translation-unit dependencies tight.
 *
 *  This header has no `extern "C"` guards of its own — each included
 *  header supplies its own. The umbrella is a pure composition, and
 *  adding guards here would nest redundantly.
 */
#ifndef IRREP_H
#define IRREP_H

#include <irrep/export.h>
#include <irrep/version.h>
#include <irrep/types.h>
#include <irrep/simd.h>
#include <irrep/so3.h>
#include <irrep/su2.h>
#include <irrep/clebsch_gordan.h>
#include <irrep/spherical_harmonics.h>
#include <irrep/solid_harmonics.h>
#include <irrep/wigner_d.h>
#include <irrep/multiset.h>
#include <irrep/multiset_2j.h>
#include <irrep/tensor_product.h>
#include <irrep/recoupling.h>
#include <irrep/radial.h>
#include <irrep/quadrature.h>
#include <irrep/time_reversal.h>
#include <irrep/parity.h>
#include <irrep/equivariant_layers.h>
#include <irrep/nequip.h>
#include <irrep/point_group.h>
#include <irrep/lattice.h>
#include <irrep/space_group.h>
#include <irrep/config_project.h>
#include <irrep/rdm.h>
#include <irrep/spin_project.h>
#include <irrep/sym_group.h>

#endif /* IRREP_H */
