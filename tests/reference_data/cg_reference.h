/* SPDX-License-Identifier: MIT */
/* Reference Clebsch-Gordan coefficients, hand-tabulated from authoritative
 * sources and independently verified against sympy's
 *
 *     sympy.physics.quantum.cg.CG(j1, m1, j2, m2, J, M).doit().simplify()
 *
 * before committing. Each entry stores the signed rational square
 * (`sign · sqrt(num / den)`) as two integers — this lets the test re-derive
 * the reference value at run time via libm `sqrt()`, so our comparison is
 * bit-aligned with how the library itself rounds intermediate results.
 *
 * Sources:
 *   Sakurai, Modern Quantum Mechanics (2e), Appendix A, Tables A.1–A.4.
 *   Varshalovich, Quantum Theory of Angular Momentum, §8.5.
 *
 * Indices are doubled (two_j, two_m) so half-integer spin is representable. */

#ifndef IRREP_TEST_CG_REFERENCE_H
#define IRREP_TEST_CG_REFERENCE_H

typedef struct {
    int two_j1, two_m1;
    int two_j2, two_m2;
    int two_J,  two_M;
    int sign;                      /* ±1 or 0 */
    int num, den;                  /* reference = sign · sqrt(num / den) */
    const char *description;
} irrep_cg_reference_t;

static const irrep_cg_reference_t IRREP_CG_REFERENCES[] = {
    /* (½, ½) coupling */
    { 1,  1,  1,  1,  2,  2,  +1,  1, 1, "⟨½, ½; ½, ½ | 1, 1⟩ = 1"              },
    { 1,  1,  1, -1,  2,  0,  +1,  1, 2, "⟨½, ½; ½, -½ | 1, 0⟩ = 1/√2"          },
    { 1, -1,  1,  1,  2,  0,  +1,  1, 2, "⟨½, -½; ½, ½ | 1, 0⟩ = 1/√2"          },
    { 1, -1,  1, -1,  2, -2,  +1,  1, 1, "⟨½, -½; ½, -½ | 1, -1⟩ = 1"           },
    { 1,  1,  1, -1,  0,  0,  +1,  1, 2, "⟨½, ½; ½, -½ | 0, 0⟩ = +1/√2"         },
    { 1, -1,  1,  1,  0,  0,  -1,  1, 2, "⟨½, -½; ½, ½ | 0, 0⟩ = -1/√2"         },

    /* (1, 1) coupling */
    { 2,  2,  2,  2,  4,  4,  +1,  1, 1, "⟨1, 1; 1, 1 | 2, 2⟩ = 1"              },
    { 2,  2,  2,  0,  4,  2,  +1,  1, 2, "⟨1, 1; 1, 0 | 2, 1⟩ = 1/√2"           },
    { 2,  0,  2,  2,  4,  2,  +1,  1, 2, "⟨1, 0; 1, 1 | 2, 1⟩ = 1/√2"           },
    { 2,  2,  2, -2,  4,  0,  +1,  1, 6, "⟨1, 1; 1, -1 | 2, 0⟩ = 1/√6"          },
    { 2,  0,  2,  0,  4,  0,  +1,  2, 3, "⟨1, 0; 1, 0 | 2, 0⟩ = √(2/3)"         },
    { 2, -2,  2,  2,  4,  0,  +1,  1, 6, "⟨1, -1; 1, 1 | 2, 0⟩ = 1/√6"          },
    { 2,  2,  2, -2,  0,  0,  +1,  1, 3, "⟨1, 1; 1, -1 | 0, 0⟩ = 1/√3"          },
    { 2,  0,  2,  0,  0,  0,  -1,  1, 3, "⟨1, 0; 1, 0 | 0, 0⟩ = -1/√3"          },
    { 2, -2,  2,  2,  0,  0,  +1,  1, 3, "⟨1, -1; 1, 1 | 0, 0⟩ = 1/√3"          },

    /* (1, ½) coupling */
    { 2,  2,  1, -1,  3,  1,  +1,  1, 3, "⟨1, 1; ½, -½ | 3/2, ½⟩ = 1/√3"        },
    { 2,  0,  1,  1,  3,  1,  +1,  2, 3, "⟨1, 0; ½,  ½ | 3/2, ½⟩ = √(2/3)"      },
    { 2,  2,  1, -1,  1,  1,  +1,  2, 3, "⟨1, 1; ½, -½ | ½, ½⟩ = √(2/3)"        },
    { 2,  0,  1,  1,  1,  1,  -1,  1, 3, "⟨1, 0; ½,  ½ | ½, ½⟩ = -1/√3"         },

    /* Selection-rule zeros — must return exactly 0.0. */
    { 2,  2,  2,  0,  4,  0,   0,  0, 1, "m-sum mismatch"                       },
    { 2,  0,  2,  0,  6,  0,   0,  0, 1, "J outside triangle"                   },
    { 2,  2,  2,  2,  0,  0,   0,  0, 1, "|M| > J"                              },
};

enum { IRREP_CG_NUM_REFERENCES =
       (int)(sizeof(IRREP_CG_REFERENCES) / sizeof(IRREP_CG_REFERENCES[0])) };

#endif /* IRREP_TEST_CG_REFERENCE_H */
