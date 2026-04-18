# libirrep — References

Every non-trivial formula or data table shipped by the library is cited here.
If you add a new formula, add it to this file in the same commit.

## Clebsch-Gordan coefficients

- G. Racah, "Theory of Complex Spectra. II," *Phys. Rev.* **62**, 438 (1942).
  Supplies the single-sum "Racah formula" used in `src/clebsch_gordan.c`.
- D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii, *Quantum Theory of
  Angular Momentum*, World Scientific (1988), §8. Selection rules and
  tabulated hand values against which our tests are checked.
- J.J. Sakurai, *Modern Quantum Mechanics*, Appendix A. Source of the hand
  values in `tests/reference_data/cg_reference.json` for small j.

## Wigner 3j / 6j / 9j

- Varshalovich §9 and Edmonds §6 for 6j / 9j definitions and the single-sum
  formulas implemented in `src/recoupling.c`.
- A.R. Edmonds, *Angular Momentum in Quantum Mechanics*, Princeton (1957).

## Spherical harmonics

- Sakurai Appendix A and Jackson, *Classical Electrodynamics*, §3, for the
  Condon-Shortley convention used throughout.
- T. Limpanuparb & J. Milthorpe, "Associated Legendre Polynomials and
  Spherical Harmonics for Efficient Computation," arXiv:1410.1748 (2014).
  Basis for the stable cartesian recurrence in `src/spherical_harmonics.c`.
- W.H. Press et al., *Numerical Recipes 3e*, §6.7, for the associated
  Legendre three-term recurrence.

## Wigner-D matrices

- Varshalovich §4.3.4 Eq. (10): Jacobi-polynomial form of the small-d matrix,
  numerically stable at large j.
- Varshalovich §4.16: tabulated closed-form small-d at j ∈ {½, 1, 3/2, 2}.
- L.C. Biedenharn & J.D. Louck, *Angular Momentum in Quantum Physics*,
  Addison-Wesley (1981).

## Rotation conversions

- F.L. Markley, "Unit Quaternion from Rotation Matrix," *J. Guidance, Control,
  and Dynamics* **31** (2), 440-442 (2008). Branch-switching log (`rot_log`).
- S.W. Shepperd, "Quaternion from Rotation Matrix," *J. Guidance and Control*
  **1** (3), 223-224 (1978). Quaternion-from-matrix extraction.
- K. Shoemake, "Uniform Random Rotations," *Graphics Gems III*, p. 124
  (1992). `irrep_quat_random`.

## Quadrature

- V.I. Lebedev & D.N. Laikov, "A quadrature formula for the sphere of the
  131st algebraic order of accuracy," *Doklady Mathematics* **59** (3),
  477-481 (1999). Public-domain node and weight tables; shipped in
  `src/quadrature_lebedev_data.c`.

## Radial basis / cutoffs

- J. Klicpera, J. Groß, S. Günnemann, "Directional Message Passing for
  Molecular Graphs," ICLR (2020). Bessel radial basis.
- S. Batzner et al., "E(3)-Equivariant Graph Neural Networks for Data-Efficient
  and Accurate Interatomic Potentials," *Nature Communications* **13**, 2453
  (2022). Cutoff envelopes, NequIP message structure.

## Tensor products / equivariance

- N. Thomas et al., "Tensor Field Networks," arXiv:1802.08219 (2018).
- M. Geiger et al., "e3nn: Euclidean Neural Networks," arXiv:2207.09453 (2022).
  Path-indexed tensor-product decomposition that libirrep mirrors.
- I. Batatia et al., "MACE: Higher Order Equivariant Message Passing Neural
  Networks for Fast and Accurate Force Fields," NeurIPS (2022).

## Time reversal

- Sakurai §4.4 for the antiunitary T operator, Kramers degeneracy, and the
  T² = ±1 dichotomy.
- F. Haake, *Quantum Signatures of Chaos*, §2 for the explicit matrix
  representation used in `src/time_reversal.c`.
