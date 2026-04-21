# Benchmark baselines

Each file is a JSON array of bench records produced by
`make bench` (or, for cross-platform, a manually-driven run of the
`build/bin/bench_*` binaries). One file per (OS, arch) pair.

Baselines are intentionally single-run snapshots: they catch
regressions (a kernel that was 100 ns/op yesterday is now 1000 ns/op),
not 1% variability. `scripts/perf_compare.sh` flags any regression > 5%
against the matching baseline.

## Files

| file                          | host              | native?    | notes |
|-------------------------------|-------------------|------------|-------|
| `darwin-arm64.json`           | Apple M2 Ultra    | native     | primary dev box; NEON kernels active |
| `darwin-x86_64-rosetta.json`  | Apple M2 Ultra    | Rosetta 2  | AVX2 kernels active via Rosetta's AVX2 translation (macOS Sequoia 15.0+) |

## Rosetta caveat

`darwin-x86_64-rosetta.json` numbers are about **2–5× slower** than
`darwin-arm64.json` for identical arithmetic. This is Rosetta's dynamic-
translation overhead, not the AVX2 path being slow. The Rosetta
baseline is valuable for:

- Correctness: bit-exact agreement between the x86_64 build and the
  arm64 build on every kernel that the test suite covers.
- Regression detection on the AVX2 code path specifically: if an AVX2
  kernel ever stops running faster than the scalar x86_64 fallback,
  this baseline moves substantially and the CI should notice.

Native x86_64 performance (Intel Skylake-X, AMD Zen 2+) is **not**
measured here. First-push CI on ubuntu-22.04 x86_64 runners will
establish a third baseline; at that point this Rosetta file becomes
a correctness oracle rather than a perf reference.

## Measured hot-path comparison (apples to Rosetta-apples)

| kernel                              | arm64 native  | x86_64 Rosetta  |
|-------------------------------------|--------------:|----------------:|
| `sph_harm_cart_all_l4`              |   3.0 ns/op   |   6.8 ns/op     |
| `sph_harm_cart_all_batch_l4_4096E`  |  38.0 ns/edge | 175.2 ns/edge   |
| `tp_apply_nequip`                   | 144.0 ns/call | 229.6 ns/call   |
| `tp_apply_nequip_weighted`          | 144.8 ns/call | 229.6 ns/call   |
| `cg_110_110_110`                    |  22.1 ns/call |  94.7 ns/call   |
| `rbf_bessel`                        |   6.2 ns/op   |  12.5 ns/op     |

Both paths validated bit-exact on the full 28-test suite.
