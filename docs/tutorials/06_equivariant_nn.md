# Tutorial 06 — Building an Equivariant MLP

Full content lands in M14. Topics covered:

- Embedding an edge vector: cartesian SH × radial basis × cutoff.
- Composing a NequIP-style message: `weighted_tp(h_j, Y(r̂_ij)) · φ(r_ij)`.
- Stacking layers: linear-on-irreps, RMS norm, gate activation.
- Training loop glue: where autograd lives (caller's responsibility).
- End-to-end equivariance check across a whole network.
