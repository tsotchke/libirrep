# Tutorial 07 — The 1.3 Kagome NQS substrate

This tutorial walks through the libirrep 1.3 primitives as they compose
into a symmetric neural-quantum-state (NQS) ansatz on the kagome lattice.
The motivation — and the 1.3 scope lock — is set by
the 1.3 CHANGELOG:
the Kagome Heisenberg S = ½ ground-state-nature problem (open since Yan,
Huse & White, *Science* **332**, 1173, 2011), and its immediate
follow-ups (J₁–J₂ square Heisenberg, 2D Hubbard at finite doping).

No library ships the full NQS driver — that lives in
`spin_based_neural_network` downstream. What libirrep 1.3 contributes is
the group-theoretic bookkeeping and the entanglement diagnostics, so the
driver can focus on the wavefunction ansatz, the MCMC sampler, and the
optimiser.

## 1. The pipeline in one picture

```
 (Lx, Ly, kind) lattice.h
 ───────────────► irrep_lattice_t
 │ sites, bonds, k-grid
 ▼
 wallpaper kind space_group.h
 ───────────────► irrep_space_group_t (p4mm / p6mm, order up to 432)
 │ site permutations π_g
 ▼
 configuration σ config_project.h
 ───────────────► orbit {g·σ : g ∈ G}
 │ wavefunction evaluation per orbit image
 ▼
 amplitudes ψ_g
 │ character-weighted reducer
 ▼
 P_μ ψ(σ) = (d_μ/|G|) Σ_g χ_μ*(g) ψ_g
```

Separately, once the variational optimisation has converged to a candidate
ground state, the Z₂-vs-trivial diagnosis runs through:

```
 learned |ψ⟩ rdm.h
 ─────────► partial trace ρ_A = Tr_B |ψ⟩⟨ψ|
 │
 ▼
 Jacobi eigvals → spectrum
 │
 ▼
 S_VN, Rényi_α, Kitaev–Preskill γ

 γ = ln 2 → Z₂ spin liquid
 γ = 0 → no topological order
```

All the pieces live on the stable C ABI and compose without allocations in
the hot loop.

## 2. Building the lattice

```c
#include <irrep/lattice.h>

irrep_lattice_t *L = irrep_lattice_build(IRREP_LATTICE_KAGOME, 6, 6);
assert(irrep_lattice_num_sites(L) == 108); /* 6·6·3 */
assert(irrep_lattice_num_bonds_nn(L) == 216); /* 4 NN per site · 108 / 2 */
assert(irrep_lattice_num_bonds_nnn(L) == 216); /* across-hex mixed sublattice */
```

The builder uses canonical direction sets per lattice (sqare, triangular,
honeycomb, kagome) and deduplicates bond candidates after PBC wrap, so the
bond list is a clean `i < j` tuple set even on small clusters where PBC
causes self-images to collapse. NN bond length is 1 on every supported
lattice, so Heisenberg J and Hubbard t transfer between lattices without
a rescale.

## 3. Attaching the wallpaper group

```c
#include <irrep/space_group.h>

irrep_space_group_t *G =
 irrep_space_group_build(L, IRREP_WALLPAPER_P6MM);
assert(irrep_space_group_order(G) == 432); /* 36 cells × 12 point ops */
assert(irrep_space_group_num_sites(G) == 108);
```

The `build` call validates that the torus spanned by `(Lx·a₁, Ly·a₂)` is
invariant under the full point group. On 6×6 kagome that holds; on
`Lx ≠ Ly` triangular or kagome it does not, and the builder returns
`NULL` with an explanatory `irrep_last_error`.

Inside, every group element is realised as a site permutation
`π_g : [0, num_sites) → [0, num_sites)`. Forward and inverse tables are
cached so a single application is a pointer read plus an array index:

```c
int image = irrep_space_group_apply(G, g, site);
```

Benchmark on M2 Ultra:

| operation | time |
| --------------------------------------------------- | ---- |
| `irrep_space_group_apply(g, s)` (one site) | 1.35 ns |
| `irrep_space_group_permutation(g, out)` (108 ints) | 7.4 ns |
| `irrep_space_group_apply_config(g, in, out)` | 40.6 ns |
| full 432-element orbit sum | 17.6 µs |

Subordinate to any realistic neural-network forward pass (target 10 ms).

## 4. Projecting onto a target irrep

The pullback action on a real-valued configuration `σ ∈ {±1}^{num_sites}`
is `(g·σ)(s) = σ(g⁻¹·s)`. The config-projection header wraps the orbit
enumeration so callers don't have to manage the inverse permutations:

```c
#include <irrep/config_project.h>

double *orbit = malloc(sizeof(double) * order * num_sites);
irrep_sg_enumerate_orbit(G, sigma, orbit);
```

For each orbit image, the caller evaluates their wavefunction amplitude
(the one expensive step: typically a neural-network forward pass):

```c
double _Complex *amps = malloc(sizeof(double _Complex) * order);
for (int g = 0; g < order; ++g) {
 amps[g] = neural_net_forward(orbit + g * num_sites);
}
```

And finally the character-weighted reduce:

```c
/* Totally-symmetric (A₁) projection: (1/|G|) Σ_g ψ(g·σ) */
double _Complex proj = irrep_sg_project_A1(G, amps);

/* Generic irrep μ with a caller-supplied character row: */
irrep_sg_irrep_t *mu = irrep_sg_irrep_new(G, chi_mu, d_mu);
double _Complex projected = irrep_sg_project_amplitude(mu, amps);
```

At the 6×6 kagome scale, step 1 (enumerate) and step 3 (reduce) together
add < 20 µs per MCMC sample, dominated by the 432 neural-network forward
passes — exactly the posture the 1.3 scope commits to.

See [`examples/kagome_a1_projection.c`](../../examples/kagome_a1_projection.c)
for a runnable demo exercising the whole chain with a toy amplitude.

## 5. Closing the physics: entanglement diagnostics

Once the variational optimisation has landed on a candidate ground state
(or a sample from its posterior), the Z₂ topological-order diagnosis runs
through the `rdm.h` primitives. On an ED-scale validation cluster:

```c
#include <irrep/rdm.h>

/* |ψ⟩ in the computational basis, length 2^N. */
/* sites_A partitions the cluster; Tr_B is implicit in the partial_trace. */
double _Complex rho_A[D * D]; /* D = 2^|A| */
irrep_partial_trace(N, 2, psi, sites_A, nA, rho_A);

double S_A = irrep_entropy_vonneumann(rho_A, D);

/* Kitaev–Preskill for a tripartition (A, B, C) of a large contractible */
/* region. γ distinguishes Z₂ (ln 2) from trivial (0) topological order. */
double gamma = irrep_topological_entanglement_entropy(
 S_A, S_B, S_C, S_AB, S_BC, S_AC, S_ABC);
```

The Hermitian eigendecomposition under the hood is a phase-reduction +
real-Givens cyclic-Jacobi sweep — no LAPACK dependency. It converges to
`10⁻¹⁴` relative precision in ~10 sweeps for dense
full-rank ρ_A.

The ED validation path is demonstrated end-to-end in
[`examples/heisenberg4_ed.c`](../../examples/heisenberg4_ed.c), which
solves the 2×2 square Heisenberg cluster (4 sites, 16-dim Hilbert space),
verifies `E_0 = −2 J`, confirms the ground state is a singlet via the
`spin_project.h` total-J projector, and reports the bipartite
entanglement entropy at a 2-vs-2 cut.

## 6. Total-J restriction (optional)

The `spin_project.h` module restricts the NQS ansatz to a fixed-total-J
sector — useful at half filling on an even number of sites, where the
Heisenberg ground state lives in `J = 0` (Marshall sign rule, Lieb–Mattis
1962):

```c
#include <irrep/spin_project.h>

/* 16·6·16 Euler grid, N sites, two_J_target = 0 for the singlet sector. */
irrep_spin_project_spin_half(
 /*two_J_target=*/0, N,
 /*n_alpha=*/16, /*n_beta=*/6, /*n_gamma=*/16,
 psi_in, psi_out);
```

Downstream, the projected `psi_out` is renormalised and becomes the
variational ansatz the sampler draws from. The singlet weight
`‖P_{J=0}|ψ⟩‖²` is a direct quality metric during training — when it drops
the learning rate is too aggressive.

## 7. Scaling beyond ED

The entanglement-entropy path above assumes `|ψ⟩` is stored explicitly —
feasible only up to ≲ 24 sites on a single machine. For the 108-site
6×6 target, ρ_A is assembled from MCMC or matrix-product-state
contractions outside libirrep; downstream passes the resulting dense `ρ_A`
directly into the entropy routines, which scale as O(n³) in the subsystem
dimension `n = 2^|A|` — comfortable up to `|A| = 12` or so (n = 4096).

The lattice, space-group, and configuration-projection primitives stay
valid at all scales: they operate on site indices, not on the full
Hilbert space.

## 8. Where to go next

- [`examples/kagome_a1_projection.c`](../../examples/kagome_a1_projection.c)
 — end-to-end 6×6 kagome A₁ projection on a reference configuration.
- [`examples/heisenberg4_ed.c`](../../examples/heisenberg4_ed.c) —
 2×2 Heisenberg ED + J=0 verification + entropy; a minimal physics
 regression test.
- the 1.3 CHANGELOG
 — scope lock for the cycle; enumerates the Kagome Heisenberg target,
 the secondary targets, and the exit criteria.
- [`docs/DESIGN.md`](../DESIGN.md) §2.9 — module responsibilities for
 the seven 1.3 headers in the layering context of the existing library.
