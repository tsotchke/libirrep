/* SPDX-License-Identifier: MIT */
/** @file hamiltonian.h
 *  @brief On-the-fly Hamiltonian apply-operators for ED / Lanczos.
 *
 *  Every example under `examples/` that does exact diagonalisation
 *  previously re-implemented the same Heisenberg `apply_op` callback
 *  by hand. This header promotes the canonical forms to library
 *  primitives, so downstream consumers plug them directly into
 *  `irrep_lanczos_eigvals` without writing their own bit-twiddling.
 *
 *  Convention for spin-½ computational basis: integer `s ∈ [0, 2^N)`
 *  represents the `N`-bit Ising configuration with bit `i = (s >> i) & 1`
 *  = {1 → up, 0 → down}. `S_z_total = (popcount(s) − N/2)/2` is exactly
 *  conserved by every model here; callers who want fixed-magnetisation
 *  ED should seed the Lanczos starting vector on a single `S_z` sector
 *  (the apply functions preserve it automatically).
 *
 *  All apply functions have the signature required by the sparse
 *  eigensolver `irrep_lanczos_eigvals`: `(const double _Complex *psi,
 *  double _Complex *out, void *opaque)`. Pass the returned opaque
 *  handle as the `ctx` argument to Lanczos.
 */
#ifndef IRREP_HAMILTONIAN_H
#define IRREP_HAMILTONIAN_H

#include <complex.h>
#include <stddef.h>

#include <irrep/export.h>
#include <irrep/types.h>

/* Forward decl — avoid circular include with config_project.h */
struct irrep_sg_rep_table;
struct irrep_sg_little_group;
struct irrep_sg_little_group_irrep;

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Opaque handle for a spin-½ Heisenberg Hamiltonian
 *  `H = J · Σ_{⟨i,j⟩ ∈ bonds} S_i · S_j`.
 *
 *  Built from an explicit NN-bond list (as produced by
 *  @ref irrep_lattice_fill_bonds_nn) plus a coupling constant. */
typedef struct irrep_heisenberg irrep_heisenberg_t;

/** @brief Construct a Heisenberg operator on `num_sites` spin-½ sites.
 *
 *  The bond list `{(bi[b], bj[b]) : b ∈ [0, num_bonds)}` is copied into
 *  the handle — the caller is free to free its arrays after this call.
 *
 *  @param num_sites  Number of sites (≤ 30 for typical 2^N Hilbert spaces).
 *  @param num_bonds  Length of the bond arrays.
 *  @param bi, bj     NN bond endpoints; `bi[b] < bj[b]` is not required
 *                    but duplicates are not deduplicated.
 *  @param J          Coupling constant (positive = antiferromagnetic).
 *  @return Owned handle, or NULL on OOM / invalid input. */
IRREP_API irrep_heisenberg_t *irrep_heisenberg_new(int num_sites, int num_bonds, const int *bi,
                                                   const int *bj, double J);

/** @brief Release a handle returned by @ref irrep_heisenberg_new. */
IRREP_API void irrep_heisenberg_free(irrep_heisenberg_t *H);

/** @brief Apply callback with the signature required by
 *         @ref irrep_lanczos_eigvals.
 *
 *  Computes `out = H · psi` on length-`2^num_sites` complex vectors.
 *  For each bond `(i, j)`:
 *    - `S_i^z · S_j^z |s⟩ = ±¼ |s⟩`   (sign = +1 if aligned, −1 if anti)
 *    - `½(S_i^+ S_j^- + S_i^- S_j^+) |s⟩ = ½ |s ⊕ 2^i ⊕ 2^j⟩`
 *      when the two spins are anti-aligned, zero otherwise.
 *
 *  Time complexity: `O(num_bonds · 2^num_sites)` per call. Memory:
 *  writes into `out`; reads `psi`. Thread-safe on disjoint `(psi, out)`
 *  pairs with a shared `opaque`. */
IRREP_API void irrep_heisenberg_apply(const double _Complex *psi, double _Complex *out,
                                      void *opaque);

/** @brief Number of sites the handle was built with. */
IRREP_API int irrep_heisenberg_num_sites(const irrep_heisenberg_t *H);

/** @brief Hilbert-space dimension (`2^num_sites`). Convenience for
 *         sizing Lanczos buffers. */
IRREP_API long long irrep_heisenberg_dim(const irrep_heisenberg_t *H);

/** @brief Build a J₁-J₂ Heisenberg on two bond sets (NN + NNN):
 *
 *    H = J₁ · Σ_{⟨i,j⟩ ∈ nn}  S_i · S_j
 *      + J₂ · Σ_{⟨i,j⟩ ∈ nnn} S_i · S_j
 *
 *  Produced as an `irrep_heisenberg_t` so the apply callback is shared
 *  with the pure-NN case; internally stored as a concatenated bond list
 *  with per-bond effective coupling. Callers fill the nn and nnn bond
 *  arrays from @ref irrep_lattice_fill_bonds_nn and
 *  @ref irrep_lattice_fill_bonds_nnn.
 *
 *  The J₁-J₂ square Heisenberg is one of the standard frustrated-
 *  magnetism benchmarks; J₂/J₁ ≈ 0.5 sits near the proposed spin-liquid
 *  window and is the secondary target the 1.3 substrate supports. */
IRREP_API irrep_heisenberg_t *irrep_heisenberg_j1j2_new(int num_sites, int num_bonds_nn,
                                                        const int *nn_i, const int *nn_j, double J1,
                                                        int num_bonds_nnn, const int *nnn_i,
                                                        const int *nnn_j, double J2);

/** @brief Build a spin-½ XY chain / lattice:
 *
 *    H = J · Σ_{⟨i,j⟩}  (S_i^x S_j^x + S_i^y S_j^y)  =  (J/2) · Σ  (S_i^+ S_j^− + h.c.)
 *
 *  No S^z-S^z term. S_z_total is still conserved, but the eigenvalue
 *  spectrum is continuous-like for large N (Bethe-ansatz soluble in 1D).
 *  Same apply callback as Heisenberg. */
IRREP_API irrep_heisenberg_t *irrep_xy_new(int num_sites, int num_bonds, const int *bi,
                                           const int *bj, double J);

/** @brief Sparse Hamiltonian apply on the orbit-representative basis in
 *         the totally-symmetric sector of the space group @p T was built on.
 *
 *  Interpretation depends on the space group used to build the rep table:
 *    - `IRREP_WALLPAPER_P1`  : pure-momentum Γ sector (pure Sz, translation-invariant)
 *    - `IRREP_WALLPAPER_P4MM`: full C_4v + translations trivial rep (Γ, A_1)
 *    - `IRREP_WALLPAPER_P6MM`: full C_6v + translations trivial rep (Γ, A_1)
 *
 *  For each representative `u` with orbit size `N_u`, the basis vector is
 *  `|ũ⟩ = (1 / √N_u) · Σ_{x ∈ orbit(u)} |x⟩`. Applying a Heisenberg bond
 *  with coefficients `(c_zz, c_pm)`:
 *  - Diagonal: `y_u += c_zz · (aligned? +1 : −1) · x_u` for every bond.
 *  - Off-diagonal: if bits `(i, j)` of `u` are anti-aligned, flip to get
 *    `u'`; canonicalise `u' → v`; then `y_v += c_pm · √(N_u / N_v) · x_u`
 *    — source-orbit-size in numerator (cf. derivation:
 *    `⟨Γ,v|H|Γ,u⟩ = ½J · √(N_u/N_v) · k_uv` for K _uv = number of
 *    anti-aligned bonds in canonical `u` flipping into `v`'s orbit).
 *
 *  Complexity: `O(num_bonds · |reps|)`. At N = 36 kagome Sz = 0, the dense
 *  path is `O(27 · 6.9e10)` while this path is `O(27 · 2.1e7)` — a 3300×
 *  reduction. This is the primitive that unlocks workstation-scale ED past
 *  N = 24 and the torus side of the KH cross-validation benchmark.
 *
 *  Matches the signature `irrep_lanczos_eigvals` / `_reorth` expect, so
 *  callers pass `H` as `opaque` and this function as the apply callback.
 *  Thread-safe on disjoint `(psi_in, psi_out)` pairs with shared context. */
IRREP_API void irrep_heisenberg_apply_in_sector(const irrep_heisenberg_t *H,
                                                const struct irrep_sg_rep_table *T,
                                                const double _Complex *psi_in,
                                                double _Complex *psi_out);

/** @brief Precomputed sector-binding: a (Hamiltonian, rep_table) pair with
 *  all `(rep, bond) → target_rep_index + coefficient` transitions cached.
 *
 *  The uncached @ref irrep_heisenberg_apply_in_sector spends most of its
 *  time canonicalising flipped configurations inside the matvec (one
 *  `irrep_space_group_apply_bits` per group element, per anti-aligned
 *  bond, per rep — O(|G| · |bonds| · |reps|) canonicalisation work per
 *  iteration). This binding canonicalises once, upfront, and the
 *  matvec becomes pure O(|bonds| · |reps|) memory-bound arithmetic.
 *
 *  Measured on kagome 3×3 (N = 27, 186616 reps, 27 bonds): ~60× speedup
 *  on Apple M2 Ultra vs the uncached path.
 *
 *  Memory: `|reps| · (num_bonds · sizeof(int64 + double) + sizeof(double))`
 *  ≈ 40 MB at the N = 27 scale. */
typedef struct irrep_sg_heisenberg_sector irrep_sg_heisenberg_sector_t;

/** @brief Build the cached sector-binding. Canonicalises every reachable
 *  `(rep, bond)` transition upfront; returns NULL on OOM. */
IRREP_API irrep_sg_heisenberg_sector_t *irrep_sg_heisenberg_sector_build(
    const irrep_heisenberg_t *H, const struct irrep_sg_rep_table *T);

/** @brief Build a cached binding at a non-trivial `(k, μ_k)` sector.
 *
 *  1D-irrep path. For each rep `u` in `T`, computes the sector-basis norm
 *  `σ_u = Σ_{s ∈ Stab(u)} conj(w_s)` where `w_g` is the composite projector
 *  weight from @ref irrep_sg_projector_weights. Reps with `σ_u = 0` are
 *  annihilated by the projector and filtered out of the sector basis.
 *
 *  Off-diagonal CSR coefficient per bond:
 *    `c_b · |G| · conj(w_{g_idx}) · √(σ_v / σ_u)`
 *  where `g_idx` is the canonicalising element mapping `flip_b(u)` to `v`.
 *  At `(Γ, A_1)` this reduces to `c_b · √(orbit_u / orbit_v)` — the
 *  trivial-sector formula from @ref irrep_sg_heisenberg_sector_build.
 *
 *  Caller passes vectors of length @ref irrep_sg_heisenberg_sector_dim
 *  (the count of non-zero-σ reps), NOT `irrep_sg_rep_table_count(T)`.
 *
 *  @param H     Heisenberg apply-operator handle from
 *               @ref irrep_heisenberg_new (carries the bond list and the
 *               coupling magnitude).
 *  @param T     symmetry-adapted rep table from
 *               @ref irrep_sg_rep_table_build (selects the popcount /
 *               Sz sector and the parent space group).
 *  @param lg    little group at the target `k` (from @ref irrep_sg_little_group_build)
 *  @param mu_k  target little-group irrep (from @ref irrep_sg_little_group_irrep_named
 *               or @ref irrep_sg_little_group_irrep_new). Must be 1D (dim = 1);
 *               2D irreps need a different code path and return NULL here. */
IRREP_API irrep_sg_heisenberg_sector_t *irrep_sg_heisenberg_sector_build_at_k(
    const irrep_heisenberg_t *H, const struct irrep_sg_rep_table *T,
    const struct irrep_sg_little_group *lg,
    const struct irrep_sg_little_group_irrep *mu_k);

IRREP_API void irrep_sg_heisenberg_sector_free(irrep_sg_heisenberg_sector_t *S);

/** @brief Dimension of the sector basis (number of non-filtered reps). */
IRREP_API long long irrep_sg_heisenberg_sector_dim(const irrep_sg_heisenberg_sector_t *S);

/** @brief Build the dense sector Hamiltonian matrix via a first-principles
 *  construction — intended as a debugging / cross-validation ORACLE for the
 *  sparse sector builders. For each representative `u` in the rep table,
 *  build `|ũ⟩ = P_μ|u⟩ / √σ_u` as an explicit 2^N-dim complex vector
 *  (using the projector weights), then sandwich `H` through pairs of basis
 *  vectors to fill `H_out[row j, col i] = ⟨ũ_j | H | ũ_i⟩`.
 *
 *  For 2D irreps, the character-projector basis has one vector per orbit
 *  (representing ONE copy of the d_μ-fold-degenerate isotypic block).
 *  Eigenvalues from this matrix match H's spectrum on V_μ (each eigenvalue
 *  appearing m_μ(u) times instead of d_μ·m_μ times).
 *
 *  Complexity: `O(|reps| · 2^N · |G|)` time, `O(|reps| · 2^N)` memory.
 *  Practical up to N ≈ 16; beyond that the sparse builder is the only
 *  option. The slow-but-direct construction makes this the definitive
 *  reference for sparse-formula cross-checks.
 *
 *  @param H      Heisenberg handle
 *  @param T      rep table (fixed popcount)
 *  @param lg     little group at k
 *  @param mu_k   irrep of the little point group (1D or 2D)
 *  @param H_out  caller-allocated `n_max * n_max` complex-double matrix
 *                (row-major); on success the top-left `dim * dim` block
 *                is filled
 *  @param n_max  row/column capacity of `H_out`
 *  @return actual dim on success, -1 on invalid input / OOM. */
IRREP_API int irrep_sg_heisenberg_sector_build_dense(const irrep_heisenberg_t *H,
                                                     const struct irrep_sg_rep_table *T,
                                                     const struct irrep_sg_little_group *lg,
                                                     const struct irrep_sg_little_group_irrep *mu_k,
                                                     double _Complex *H_out, int n_max);

/** @brief Serialise @p S to @p path (binary format). The file encodes
 *  whether the binding was built via the trivial-sector path (real
 *  coefficients) or the (k, μ_k) path (complex coefficients). Returns
 *  0 on success, -1 on I/O error.
 *
 *  The space group and rep table are NOT serialised — the caller must
 *  reconstruct them separately and pass them in on load for cross-
 *  validation. This keeps serialisation format-independent of ABI
 *  changes to those structures. */
IRREP_API int irrep_sg_heisenberg_sector_save(const irrep_sg_heisenberg_sector_t *S,
                                              const char *path);

/** @brief Load a sector binding previously saved via
 *  @ref irrep_sg_heisenberg_sector_save. Validates that the on-disk
 *  `n_reps` matches the caller-supplied expected dimension (pass -1 to
 *  skip validation). Returns NULL on any mismatch or I/O error. */
IRREP_API irrep_sg_heisenberg_sector_t *irrep_sg_heisenberg_sector_load(
    const char *path, long long expected_dim);

/** @brief Matvec using the cached binding. Signature matches the
 *  Lanczos apply-callback contract — pass `S` as `ctx`. */
IRREP_API void irrep_sg_heisenberg_sector_apply(const double _Complex *psi_in,
                                                double _Complex *psi_out, void *S);

#ifdef __cplusplus
}
#endif

#endif /* IRREP_HAMILTONIAN_H */
