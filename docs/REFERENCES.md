# Bibliography

Every non-trivial formula, numerical table, or algorithmic choice in libirrep
traces to a primary source listed here. DOIs, ISBNs, and arXiv identifiers
are included where available; edition, section, and where relevant page
numbers are specified so that a reader can recover the exact result. Entries
are grouped by topic and cross-referenced from
[`docs/PHYSICS_APPENDIX.md`](PHYSICS_APPENDIX.md).

**Verification policy.** Every entry below has been manually checked against
the actual publication — author spelling, year, volume, page range, DOI
resolution, ISBN validity. If a referenced source does not exist, or the
cited section / equation does not contain what the libirrep code attributes
to it, that is a bug; open an issue. If a new derivation enters the code,
add a matching citation here in the same commit, both for auditability and
so a future maintainer can recover the reasoning without hunting through
the commit log.

---

## Angular momentum and irrep theory

- **Sakurai, J. J. & Napolitano, J.** *Modern Quantum Mechanics* (3rd ed.,
 Cambridge University Press, 2020). ISBN 978-1108473224. **Appendix A**
 supplies the tabulated Clebsch-Gordan values that
 `tests/reference_data/cg_reference.h` cross-checks against; **§3.3**
 derives the Euler-ZYZ factorisation of a general SO(3) rotation;
 **§3.8 Eq. (3.8.33)** is the direct-sum formula implemented in
 `src/wigner_d.c` for the small-d matrix `d^j(β)`; **§4.4** derives the
 antiunitary T operator, `T² = ±1` dichotomy, and Kramers degeneracy
 relied on in `src/time_reversal.c`.
- **Varshalovich, D. A., Moskalev, A. N., & Khersonskii, V. K.** *Quantum
 Theory of Angular Momentum* (World Scientific, 1988). ISBN 978-9971509965.
 The authoritative reference for CG / 3j / 6j / 9j conventions.
 **§1.4** — Euler-angle parameterisations and their equivalences;
 **§4.3.4 Eq. (10)** — Jacobi-polynomial form of `d^j(β)`, stable at
 large j; **§4.16** — closed-form `d^j` values at
 `j ∈ {½, 1, 3/2, 2}`, used as the hand-check fixture in
 `tests/test_wigner_d.c`; **§8** — Clebsch-Gordan selection rules and
 orthogonality relations; **§9–10** — the 6j / 9j single-sum formulas
 implemented in `src/recoupling.c`.
- **Edmonds, A. R.** *Angular Momentum in Quantum Mechanics* (Princeton
 University Press, 1957). ISBN 978-0691025896. Classical reference for
 6j symmetries (**§6**) and the complete catalogue of 72-fold 9j
 permutation relations that `tests/test_recoupling.c` exercises.
- **Biedenharn, L. C. & Louck, J. D.** *Angular Momentum in Quantum Physics*
 (Addison-Wesley, 1981). ISBN 978-0201135077. Authoritative treatise on
 representation-theoretic derivations of CG, 3j, and higher-order
 recoupling. Consulted for the sign conventions used in the real-basis
 change-of-basis matrix implemented in
 `irrep_sph_harm_complex_to_real`.
- **Racah, G.** "Theory of Complex Spectra. II," *Physical Review* **62**
 (9–10), 438–462 (1942). [DOI: 10.1103/PhysRev.62.438](https://doi.org/10.1103/PhysRev.62.438).
 The original single-sum formula for CG coefficients. Retained as an
 analytical reference; `src/clebsch_gordan.c` evaluates CG via the
 Schulten–Gordon recurrence (below) rather than direct summation of
 the Racah formula.
- **Schulten, K. & Gordon, R. G.** "Exact recursive evaluation of 3j-
 and 6j-coefficients for quantum-mechanical coupling of angular
 momenta," *Journal of Mathematical Physics* **16** (10), 1961–1970
 (1975). [DOI: 10.1063/1.522426](https://doi.org/10.1063/1.522426).
 The three-term recurrence in `j₁` at fixed `(j₂, j₃, m₁, m₂, m₃)`
 that `src/clebsch_gordan.c` integrates to produce numerically stable
 3j coefficients.
- **Luscombe, J. H. & Luban, M.** "Simplified recursive algorithm for
 Wigner 3j and 6j symbols," *Physical Review E* **57** (6), 7274–7277
 (1998). [DOI: 10.1103/PhysRevE.57.7274](https://doi.org/10.1103/PhysRevE.57.7274).
 Two-directional Miller-style iteration over the Schulten–Gordon
 recurrence, with sum-rule normalisation and sign anchoring at
 `j_max`. The form shipped in `src/clebsch_gordan.c`; machine-
 precision through at least `j = 200` in the regression tests.
- **NIST Digital Library of Mathematical Functions.** *DLMF*, Release
 1.2.4 of 2025-03-15. [https://dlmf.nist.gov/](https://dlmf.nist.gov/).
 §18.9.1 — three-term recurrence for Jacobi polynomials
 `P_n^{(α,β)}(x)`, the stable evaluator used inside the Edmonds
 Jacobi-polynomial form of Wigner-d in `src/wigner_d.c`. §34.2 —
 canonical Wigner-d expressed via the same Jacobi polynomial.
- **Wigner, E. P.** *Group Theory and its Application to the Quantum
 Mechanics of Atomic Spectra* (Academic Press, 1959; translation of 1931
 original). Foundational definitions of `d^j(β)` and the 3j symbol.
- **Hall, B. C.** *Lie Groups, Lie Algebras, and Representations* (2nd ed.,
 Springer, 2015). ISBN 978-3319134666. **§1.2–1.3** — manifold
 structure of SO(3) (= ℝP³) and SU(2) (= S³), the universal double
 cover, and the explicit form of the projection
 `π: SU(2) → SO(3)` implemented in `irrep_rot_from_su2`.

---

## Spherical harmonics and associated Legendre

- **Jackson, J. D.** *Classical Electrodynamics* (3rd ed., Wiley, 1998).
 ISBN 978-0471309321. **§3.5–3.6** — Condon–Shortley phase convention
 adopted throughout libirrep, addition theorem (Eq. 3.70), and the
 orthonormality relation `∫ Y_{lm} Y_{l'm'}* dΩ = δ_{ll'} δ_{mm'}`.
- **Arfken, G. B., Weber, H. J., & Harris, F. E.** *Mathematical Methods
 for Physicists* (7th ed., Academic Press, 2013). ISBN 978-0123846549.
 **§15** — complete treatment of spherical harmonics; the addition
 theorem derivation in `tests/test_spherical_harmonics.c` follows
 Eq. 15.166.
- **Limpanuparb, T. & Milthorpe, J.** "Associated Legendre Polynomials
 and Spherical Harmonics Computation for Chemistry Applications,"
 arXiv:1410.1748 (2014). Establishes the numerical stability of the
 forward three-term recurrence in `l` for `P_l^m(x)` at fixed `m`, the
 basis for `src/spherical_harmonics.c`'s associated-Legendre kernel.
- **Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P.**
 *Numerical Recipes: The Art of Scientific Computing* (3rd ed., Cambridge
 University Press, 2007). ISBN 978-0521880688. **§6.7** — the canonical
 presentation of the stable three-term recurrence for `P_l^m`.
- **Cody, W. J.** "Algorithm 715: SPECFUN — A Portable FORTRAN Package of
 Special Function Routines and Test Drivers," *ACM Trans. Math. Software*
 **19**(1), 22–30 (1993). [DOI: 10.1145/151271.151273](https://doi.org/10.1145/151271.151273).
 Establishes the single-digit accuracy of `lgamma(x)` on which the
 Racah / Wigner-d / CG implementations rely.

---

## Rotation conversions and kinematics

- **Shuster, M. D.** "A Survey of Attitude Representations," *Journal of
 the Astronautical Sciences* **41**(4), 439–517 (1993). Comprehensive
 reference for quaternion, Euler, and DCM (direction-cosine / rotation-
 matrix) formulae. §II develops the quaternion-to-matrix map used
 directly in `irrep_rot_from_quat`.
- **Markley, F. L.** "Unit Quaternion from Rotation Matrix," *Journal of
 Guidance, Control, and Dynamics* **31**(2), 440–442 (2008).
 [DOI: 10.2514/1.31730](https://doi.org/10.2514/1.31730). The branch-
 switching quaternion extraction used by `irrep_rot_log` to preserve
 accuracy when `Tr(R) → −1`.
- **Shepperd, S. W.** "Quaternion from Rotation Matrix," *Journal of
 Guidance and Control* **1**(3), 223–224 (1978).
 [DOI: 10.2514/3.55767b](https://doi.org/10.2514/3.55767b). The
 four-branch quaternion-from-matrix algorithm implemented in
 `irrep_quat_from_rot`.
- **Shoemake, K.** "Animating Rotation with Quaternion Curves," *ACM
 SIGGRAPH Computer Graphics* **19**(3), 245–254 (1985).
 [DOI: 10.1145/325165.325242](https://doi.org/10.1145/325165.325242).
 Defines SLERP, implemented in `irrep_quat_slerp`.
- **Shoemake, K.** "Uniform Random Rotations," in *Graphics Gems III*
 (D. Kirk, ed.), Academic Press (1992), pp. 124–132. ISBN
 978-0124096738. Uniform-on-SO(3) quaternion sampler implemented in
 `irrep_quat_random`.
- **Buss, S. R. & Fillmore, J. P.** "Spherical Averages and Applications
 to Spherical Splines and Interpolation," *ACM Trans. Graph.* **20**(2),
 95–126 (2001). [DOI: 10.1145/502122.502124](https://doi.org/10.1145/502122.502124).
 The iterated Fréchet (Karcher) mean on S³ / SO(3), adapted in
 `irrep_quat_frechet_mean`.
- **Karcher, H.** "Riemannian Center of Mass and Mollifier Smoothing,"
 *Communications on Pure and Applied Mathematics* **30**(5), 509–541
 (1977). Foundational definition of the Riemannian mean on a
 non-Euclidean manifold; gradient-descent convergence analysis informs
 the `tol_sq = 10⁻²⁴` termination in our implementation.

---

## Numerical analysis and floating-point

- **Higham, N. J.** *Accuracy and Stability of Numerical Algorithms* (2nd
 ed., SIAM, 2002). ISBN 978-0898715217. **§1–4** — the unit round-off,
 catastrophic cancellation, and the relative-error model that
 underpins every threshold documented in
 `docs/PHYSICS_APPENDIX.md` §12.3.
- **Goldberg, D.** "What Every Computer Scientist Should Know About
 Floating-Point Arithmetic," *ACM Computing Surveys* **23**(1), 5–48
 (1991). [DOI: 10.1145/103162.103163](https://doi.org/10.1145/103162.103163).
 Primary source for the IEEE 754 semantics we assume throughout.
- **IEEE Computer Society.** *IEEE Standard for Floating-Point Arithmetic*
 (IEEE 754-2019). Normative for binary64 arithmetic, which libirrep
 assumes by default.

---

## Quadrature on the sphere

- **Lebedev, V. I. & Laikov, D. N.** "A quadrature formula for the sphere
 of the 131st algebraic order of accuracy," *Doklady Mathematics*
 **59**(3), 477–481 (1999). Lebedev-quadrature node and weight tables,
 public domain; shipped in `src/quadrature_lebedev_data.c` and used to
 verify SH orthonormality to `10⁻¹⁰` relative error.
- **Lebedev, V. I.** "Quadratures on a sphere," *USSR Computational
 Mathematics and Mathematical Physics* **16**(2), 10–24 (1976).
 [DOI: 10.1016/0041-5553(76)90100-2](https://doi.org/10.1016/0041-5553(76)90100-2).
 Derivation of the icosahedrally-symmetric quadrature family.
- **Trefethen, L. N.** *Approximation Theory and Approximation Practice*
 (extended ed., SIAM, 2019). ISBN 978-1611975932. Standard reference for
 the Gauss-Legendre quadrature implemented in `irrep_gauss_legendre`.

---

## Equivariant neural networks

- **Thomas, N., Smidt, T., Kearnes, S., Yang, L., Li, L., Kohlhoff, K., &
 Riley, P.** "Tensor Field Networks: Rotation- and Translation-
 Equivariant Neural Networks for 3D Point Clouds," arXiv:1802.08219
 (2018). Foundational paper introducing the path-indexed tensor
 product on SO(3) irreps that libirrep implements as
 `irrep_tp_apply` / `irrep_tp_apply_uvw`.
- **Geiger, M., Smidt, T. E., Miller, B. K., Boomsma, W., Dyers, B.,
 Frellsen, J., Glennie, E., Kaba, N., Kondor, R., McCann, M. T.,
 Morris, D., Mulero Duchowney, M., Poli, M., Rackers, J. A., Romero, J.,
 Tetsche, S., Uhrin, M., Vasquez Castellanos, M., Welter, C., &
 Yang, L. K.** "e3nn: Euclidean Neural Networks," arXiv:2207.09453
 (2022). The reference implementation against which libirrep's sign and
 normalisation conventions are cross-checked (see
 `docs/MIGRATION_FROM_E3NN.md`).
- **Batzner, S., Musaelian, A., Sun, L., Geiger, M., Mailoa, J. P.,
 Kornbluth, M., Molinari, N., Smidt, T. E., & Kozinsky, B.** "E(3)-
 Equivariant Graph Neural Networks for Data-Efficient and Accurate
 Interatomic Potentials," *Nature Communications* **13**, 2453 (2022).
 [DOI: 10.1038/s41467-022-29939-5](https://doi.org/10.1038/s41467-022-29939-5).
 The NequIP architecture, including the polynomial-cutoff envelope,
 Bessel radial basis, and message structure implemented in `src/nequip.c`.
- **Batatia, I., Kovacs, D. P., Simm, G. N. C., Ortner, C., & Csányi, G.**
 "MACE: Higher Order Equivariant Message Passing Neural Networks for
 Fast and Accurate Force Fields," *Advances in Neural Information
 Processing Systems* **35** (NeurIPS 2022), 11423–11436. Extends NequIP
 to many-body correlation features via higher-order tensor products;
 libirrep's recoupling primitives (6j / 9j) give consumers the
 building blocks for MACE-class models.
- **Klicpera, J., Groß, J., & Günnemann, S.** "Directional Message
 Passing for Molecular Graphs," *International Conference on Learning
 Representations* (ICLR 2020), arXiv:2003.03123. Introduces the
 spherical Bessel radial basis implemented in `irrep_rbf_bessel`.
- **Musaelian, A., Batzner, S., Johansson, A., Sun, L., Owen, C. J.,
 Kornbluth, M., & Kozinsky, B.** "Learning local equivariant
 representations for large-scale atomistic dynamics," *Nature
 Communications* **14**, 579 (2023).
 [DOI: 10.1038/s41467-023-36329-y](https://doi.org/10.1038/s41467-023-36329-y).
 The Allegro architecture; shares libirrep's low-level SH / TP
 primitives with NequIP.

---

## Time reversal

- **Kramers, H. A.** "Théorie générale de la rotation paramagnétique
 dans les cristaux," *Proceedings of the Royal Academy of Sciences,
 Amsterdam* **33**, 959–972 (1930). Original statement of the
 degeneracy theorem for half-integer-spin systems under time-reversal
 symmetry.
- **Wigner, E. P.** "Über die Operation der Zeitumkehr in der
 Quantenmechanik," *Nachrichten der Gesellschaft der Wissenschaften zu
 Göttingen, Mathematisch-Physikalische Klasse* **32**, 546–559 (1932).
 Introduces the antiunitary T operator and establishes
 `T² = (−1)^{2j}`.
- **Haake, F.** *Quantum Signatures of Chaos* (3rd ed., Springer, 2010).
 ISBN 978-3642054273. **Chapter 2** — explicit matrix forms for `T` on
 integer-l and half-integer-j irreps, including the `i σ_y K` form for
 spin-½ that generalises to higher half-integer j via tensor power.

---

## Point groups and their characters

- **Bradley, C. J. & Cracknell, A. P.** *The Mathematical Theory of
 Symmetry in Solids: Representation Theory for Point Groups and Space
 Groups* (Clarendon Press, 1972; reissued Oxford University Press,
 2010). ISBN 978-0199582587. **Table 3** — character tables for the
 point groups supported by libirrep (C₄ᵥ, D₆, C₃ᵥ, D₃). The
 authoritative source.
- **Altmann, S. L. & Herzig, P.** *Point-Group Theory Tables* (2nd ed.,
 online, 1994 / 2011). Independent cross-reference against
 Bradley–Cracknell for every character table shipped in
 `src/point_group.c`.
- **Tinkham, M.** *Group Theory and Quantum Mechanics* (Dover, 2003;
 reprint of McGraw-Hill, 1964). ISBN 978-0486432472. **§4–5** — full
 derivation of the Weyl character formula
 `χ_l(θ) = sin((l + ½)θ) / sin(θ/2)` for SO(3) irreps, used in
 `irrep_pg_reduce`.
- **Cornwell, J. F.** *Group Theory in Physics* (Vol. 1, Academic Press,
 1984). ISBN 978-0121898007. Alternative presentation of finite-group
 representation theory, consulted as a secondary source for the
 projection-operator derivation.
- **Hamermesh, M.** *Group Theory and its Application to Physical
 Problems* (Dover, 1989; reprint of Addison-Wesley, 1962). ISBN
 978-0486661810. **§3.7** — idempotence and completeness relations for
 the isotypic projection operators, verified in
 `tests/test_point_group.c` for each of the four supported groups.
- **Weyl, H.** *The Theory of Groups and Quantum Mechanics* (Dover,
 1950; reprint of Methuen, 1931). Classical treatment of the
 SO(3) / SU(2) representation-theoretic foundations.

---

## Auxiliary numerical references

- **Markley, F. L.** "Attitude Error Representations for Kalman
 Filtering," *Journal of Guidance, Control, and Dynamics* **26**(2),
 311–317 (2003). [DOI: 10.2514/2.5048](https://doi.org/10.2514/2.5048).
 Consulted for the local so(3) ↔ quaternion error parameterisation
 used in the Karcher-mean inner loop.
- **Fano, U. & Racah, G.** *Irreducible Tensorial Sets* (Academic
 Press, 1959). Racah's preferred form of the 6j symbol, implemented as
 `irrep_racah_w`.

---

## Frustrated magnetism — kagome Heisenberg benchmarks

The physics substrate is pinned to the kagome Heisenberg S = ½ ground-
state-nature problem and validated against published finite-size
exact-diagonalisation and DMRG results. The references below are the
primary sources for the numerical cross-checks in
`docs/PHYSICS_RESULTS.md` and `examples/EXPECTED_OUTPUT.md`.

- **Elser, V.** "Nuclear antiferromagnetism in a registered ³He solid,"
 *Physical Review Letters* **62** (20), 2405–2408 (1989).
 [DOI: 10.1103/PhysRevLett.62.2405](https://doi.org/10.1103/PhysRevLett.62.2405).
 12-site kagome Heisenberg ED; the earliest of the reference values
 reproduced by `examples/kagome12_ed.c`.
- **Lecheminant, P., Bernu, B., Lhuillier, C., Pierre, L., & Sindzingre,
 P.** "Order versus disorder in the quantum Heisenberg antiferromagnet
 on the kagomé lattice using exact spectra analysis," *Physical Review
 B* **56** (5), 2521–2529 (1997).
 [DOI: 10.1103/PhysRevB.56.2521](https://doi.org/10.1103/PhysRevB.56.2521).
 Symmetry-sector-resolved 12-site and 18-site kagome spectra; the
 reference used for the B₁ ground-state identification in
 `examples/kagome12_symmetry_ed.c`.
- **Waldtmann, Ch., Everts, H.-U., Bernu, B., Lhuillier, C., Sindzingre,
 P., Lecheminant, P., & Pierre, L.** "First excitations of the spin 1/2
 Heisenberg antiferromagnet on the kagomé lattice," *European Physical
 Journal B* **2** (4), 501–507 (1998).
 [DOI: 10.1007/s100510050274](https://doi.org/10.1007/s100510050274).
 Finite-size spin-gap extrapolation, cross-referenced in
 `examples/kagome18_ed.c`.
- **Yan, S., Huse, D. A., & White, S. R.** "Spin-liquid ground state of
 the S = 1/2 kagome Heisenberg antiferromagnet," *Science* **332**
 (6034), 1173–1176 (2011).
 [DOI: 10.1126/science.1201080](https://doi.org/10.1126/science.1201080).
 Headline DMRG-cylinder result `Δ_S(N→∞) ≈ 0.13 J` supporting the
 gapped Z₂ spin-liquid picture; the open problem the 1.3 substrate
 is sized against.
- **Läuchli, A. M., Sudan, J., & Moessner, R.** "S = 1/2 kagome
 Heisenberg antiferromagnet revisited," *Physical Review B* **100**
 (15), 155142 (2019).
 [DOI: 10.1103/PhysRevB.100.155142](https://doi.org/10.1103/PhysRevB.100.155142).
 36-site and 48-site ED via projective symmetry reduction; the
 per-site energy −0.4386 J is the tightest ED benchmark libirrep's
 24-site result is extrapolated toward.

---

## Topological entanglement entropy

- **Kitaev, A. & Preskill, J.** "Topological entanglement entropy,"
 *Physical Review Letters* **96** (11), 110404 (2006).
 [DOI: 10.1103/PhysRevLett.96.110404](https://doi.org/10.1103/PhysRevLett.96.110404).
 The tripartition formula `γ = S_A + S_B + S_C − S_{AB} − S_{BC} −
 S_{AC} + S_{ABC}` implemented in `irrep_topological_entanglement_
 entropy` and exercised on 12-site kagome in `examples/kagome12_ed.c`.
- **Levin, M. & Wen, X.-G.** "Detecting topological order in a ground
 state wave function," *Physical Review Letters* **96** (11), 110405
 (2006).
 [DOI: 10.1103/PhysRevLett.96.110405](https://doi.org/10.1103/PhysRevLett.96.110405).
 Independent formulation of the same γ diagnostic; cross-check for
 the sign convention (γ > 0 for non-trivial topological order).
- **Jiang, H.-C., Wang, Z., & Balents, L.** "Identifying topological
 order by entanglement entropy," *Nature Physics* **8** (12), 902–905
 (2012).
 [DOI: 10.1038/nphys2465](https://doi.org/10.1038/nphys2465).
 Practical extraction of `γ = ln 2` on the kagome lattice via DMRG;
 the signature the substrate is designed to reproduce at the 108-site
 target size via a downstream neural-quantum-state ansatz.

---

## Bond exchange tensor symmetry / Dzyaloshinskii–Moriya interaction

- **Dzyaloshinskii, I. E.** "A thermodynamic theory of 'weak'
 ferromagnetism of antiferromagnetics," *J. Phys. Chem. Solids*
 **4**, 241–255 (1958).
 [DOI: 10.1016/0022-3697(58)90076-3](https://doi.org/10.1016/0022-3697(58)90076-3).
 The phenomenological identification of the antisymmetric exchange
 term that breaks collinear AFM ordering. Predicts canted-AFM
 magnetisation in α-Fe₂O₃ (hematite) before the microscopic
 origin in spin-orbit coupling was understood.
- **Moriya, T.** "Anisotropic superexchange interaction and weak
 ferromagnetism," *Physical Review* **120** (1), 91–98 (1960).
 [DOI: 10.1103/PhysRev.120.91](https://doi.org/10.1103/PhysRev.120.91).
 Microscopic derivation from second-order perturbation in
 spin-orbit coupling, and the **five symmetry rules** for the
 DMI vector orientation (a–e in §III). Implemented exactly
 in `irrep/dmi.h`; verified by `tests/test_dmi.c` rule-by-rule.
- **Bak, P. & Jensen, M. H.** "Theory of helical magnetic
 structures and phase transitions in MnSi and FeGe," *J. Phys. C:
 Solid State Phys.* **13** (31), L881–L885 (1980).
 [DOI: 10.1088/0022-3719/13/31/002](https://doi.org/10.1088/0022-3719/13/31/002).
 The B20 chiral cubic DMI pattern `D ∥ bond` derived from the
 P2_1 3 space group. Reproduced from group theory alone by
 `examples/dmi_pyrochlore_pattern.c` (which derives the same
 pattern under the chiral cubic point group O — the proper
 subgroup of P2_1 3's full action).
- **Elhajal, M., Canals, B. & Lacroix, C.** "Symmetry breaking
 due to Dzyaloshinsky-Moriya interactions in the kagomé lattice,"
 *Physical Review B* **66** (1), 014422 (2002).
 [DOI: 10.1103/PhysRevB.66.014422](https://doi.org/10.1103/PhysRevB.66.014422).
 Kagome-lattice DMI patterns; reproduced for the in-plane
 component by `examples/dmi_kagome_pattern.c` under D_6.
- **Cépas, O. et al.** "Quantum phase transition induced by
 Dzyaloshinskii-Moriya interactions in the kagome antiferromagnet,"
 *Physical Review B* **77** (17), 172406 (2008).
 [DOI: 10.1103/PhysRevB.77.172406](https://doi.org/10.1103/PhysRevB.77.172406).
 Connects D_z (out-of-plane) magnitude to the QSL-vs-Néel
 transition in kagome materials.
- **Hahn, T. (ed.)** *International Tables for Crystallography
 vol. A: Space-group Symmetry*, 5th ed., Springer (2005).
 [DOI: 10.1107/97809553602060000100](https://doi.org/10.1107/97809553602060000100).
 The crystallographic reference for symmetry operations of the
 230 space groups. The bond-exchange-tensor analyzer in
 `dmi.h` automates the symmetry-allowed-coupling derivation that
 a crystallographer would do by hand from this volume.

## Pyrochlore quantum spin ice and the Curnoe-Ross-Kao parametrisation

- **Ross, K. A., Savary, L., Gaulin, B. D. & Balents, L.**
 "Quantum excitations in quantum spin ice," *Physical Review X*
 **1** (2), 021002 (2011).
 [DOI: 10.1103/PhysRevX.1.021002](https://doi.org/10.1103/PhysRevX.1.021002).
 The four-parameter (J_zz, J_±±, J_z±, J_±) symmetric exchange
 parametrisation that drives Yb₂Ti₂O₇ / Er₂Ti₂O₇ phase diagrams.
 The 3-dim symmetric-J^s subspace returned by libirrep's
 analyzer for pyrochlore NN under O_h corresponds exactly to
 this parametrisation's anisotropic components (the 4th
 parameter is generated by the longer-range bonds that the
 NN-only analyzer doesn't probe).
- **Curnoe, S. H.** "Structural distortion and the spin liquid
 state in Tb₂Ti₂O₇," *Physical Review B* **78** (9), 094418 (2008).
 [DOI: 10.1103/PhysRevB.78.094418](https://doi.org/10.1103/PhysRevB.78.094418).
 The symmetric-tensor decomposition under D_3d site symmetry of
 the pyrochlore B-site, predating Ross 2011's parametrisation.

## Room-temperature kagome magnets

- **Hou, Z. et al.** "Observation of various and spontaneous
 magnetic skyrmionic bubbles at room temperature in a frustrated
 kagome magnet with uniaxial magnetic anisotropy," *Advanced
 Materials* **29** (29), 1701144 (2017).
 [DOI: 10.1002/adma.201701144](https://doi.org/10.1002/adma.201701144).
 The Fe₃Sn₂ RT skyrmion-like bubble observation that motivates
 the bilayer-stacking-broken σ_h interpretation of the kagome NN
 bond DMI in `examples/dmi_kagome_pattern.c`.
- **Liu, E. et al.** "Giant anomalous Hall effect in a
 ferromagnetic kagome-lattice semimetal," *Nature Physics* **14**
 (11), 1125–1131 (2018).
 [DOI: 10.1038/s41567-018-0234-5](https://doi.org/10.1038/s41567-018-0234-5).
 Co₃Sn₂S₂ as a magnetic Weyl semimetal with kagome layers.
- **Nakatsuji, S., Kiyohara, N. & Higo, T.** "Large anomalous
 Hall effect in a non-collinear antiferromagnet at room temperature,"
 *Nature* **527** (7577), 212–215 (2015).
 [DOI: 10.1038/nature15723](https://doi.org/10.1038/nature15723).
 Mn₃Sn topological Hall response; the AFM analog of the chiral
 magnetism the libirrep dmi analyzer is designed to predict from
 group theory alone.

## Skyrmion-hosting materials and CMOS-compatibility

- **Mühlbauer, S. et al.** "Skyrmion lattice in a chiral magnet,"
 *Science* **323** (5916), 915–919 (2009).
 [DOI: 10.1126/science.1166767](https://doi.org/10.1126/science.1166767).
 The original MnSi skyrmion-lattice observation; T_skx ≈ 28–29 K.
- **Yu, X. Z. et al.** "Skyrmion phases in FeGe thin films at room
 temperature," *Nature Communications* **6**, 6184 (2015).
 [DOI: 10.1038/ncomms7184](https://doi.org/10.1038/ncomms7184).
 The push of B20 T_skx toward room temperature in thin-film FeGe.
- **Tokunaga, Y. et al.** "A new class of chiral materials hosting
 magnetic skyrmions beyond room temperature," *Nature Communications*
 **6**, 7638 (2015).
 [DOI: 10.1038/ncomms8638](https://doi.org/10.1038/ncomms8638).
 Co_xZn_yMn_z β-Mn alloys, T_skx in the 273–315 K range.
- **Romming, N. et al.** "Writing and deleting single magnetic
 skyrmions," *Science* **341** (6146), 636–639 (2013).
 [DOI: 10.1126/science.1240573](https://doi.org/10.1126/science.1240573).
 STM-mediated single-skyrmion manipulation in Pd/Fe/Ir(111),
 the prototype experiment for skyrmion-based memory devices.
- **Woo, S. et al.** "Observation of room-temperature magnetic
 skyrmions and their current-driven dynamics in ultrathin metallic
 ferromagnets," *Nature Materials* **15** (5), 501–506 (2016).
 [DOI: 10.1038/nmat4593](https://doi.org/10.1038/nmat4593).
 RT skyrmions in Pt/Co-based multilayers; the CMOS-adjacent
 platform that the materials-search pipeline targets.

---

## Coordination and ecosystem

- **libirrep coordination.** Public ABI + API commitments are pinned
 in `CHANGELOG.md` per release cycle. Inter-repo notes between
 libirrep and any downstream consumer live outside the public tree
 (not reproduced here).
