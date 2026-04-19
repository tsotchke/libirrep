# test_downstream_compat/

Golden-vector test fixtures vendored from downstream consumers so that any
convention drift on the libirrep side fails here before the downstream CI
picks it up.

## `torque_net_tp_paths/`

From `spin_based_neural_network`. Five fixed `(h_in, edge_vec)` configurations
and their expected tensor-product outputs under the UVW TP formulation of
the torque-net's five Cartesian basis terms.

Parity is multiplicative under the TP path filter, so `1o ⊗ 1o` can only
land in even-parity output irreps. The cross product of two polar vectors
is therefore built via the `(1, 1, 1)` path into an axial `1e` output, not
`1o` — see `docs/tutorials/05_tensor_products.md` §7 for the full sign /
normalisation derivation and `examples/torque_net_tp_paths.c` for the
bit-exact (~2e-16) cartesian round trip.

```
 T1 = (m_j · r̂) · m_i — (1,1,0) → scalar in 0e, then (0,1,1) → m_i
 T2 = m_j × r̂ (axial) — (1,1,1) into 1e, prefactor √2
 T3 = m_i × m_j (axial) — (1,1,1) into 1e, prefactor √2
 T4 = (m_i · m_j) · m_i — (1,1,0) → scalar in 0e, then (0,1,1) → m_i
 T5 = m_j — passthrough
```

Prefactor convention:
- `1o ⊗ 1o → 0e` produces `(1/√3) · (a · b)` — multiply by `√3` for bare dot.
- `1o ⊗ 1o → 1e` produces `(1/√2) · (a × b)` in the `(y, z, x)` real-SH
 layout — multiply by `√2` and permute to `(x, y, z)` for bare cross.

Downstream consumers whose torque-net treats `m` as axial (pseudovector)
should use `1x1e` inputs and swap the output parities accordingly
(`1e ⊗ 1o → 1o` for polar cross product, etc.).

### How this test runs (wired in as `tests/test_downstream_compat.c` once
`torque_net_tp_paths/vectors.json` lands):

1. Read each `(h_in, edge_vec)` input + expected output from
 `torque_net_tp_paths/vectors.json`.
2. Build the UVW tensor-product descriptor that matches each term.
3. Compute `irrep_tp_apply_uvw` and assert bit-equality with the expected
 output (modulo the pinned prefactor).
4. Exit non-zero on any mismatch; the `# FAIL` line names which of T1..T5
 broke, making convention drift trivial to localise.

### Pending ship from downstream

- `vectors.json` with five `(h_in, edge_vec, expected)` triples.
- Bridge-side test `spin_based_neural_network/tests/test_torque_net_irrep.c`
 asserts the same golden vectors from the other direction.

Once the JSON lands this README becomes obsolete — replaced by the comment
block atop the test harness.
