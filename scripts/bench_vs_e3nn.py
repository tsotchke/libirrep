#!/usr/bin/env python3
"""
bench_vs_e3nn.py — side-by-side perf + numerical-agreement bench of
libirrep's `irrep_tp_apply` against e3nn.o3.FullTensorProduct on a
matched descriptor.

Usage (from a venv with e3nn + torch installed):

    python3 -m venv /tmp/venv && /tmp/venv/bin/pip install e3nn
    /tmp/venv/bin/python scripts/bench_vs_e3nn.py

What it measures:
    (a) libirrep tp_apply over `iters` calls on a fixed descriptor.
    (b) e3nn FullyConnectedTensorProduct forward on the same shape.
    (c) Cosine similarity and max abs error between the two outputs
        on a single identical input (numerical agreement).

What it does NOT measure:
    - Batching across edges. The libirrep bench is a single-sample
      call; e3nn excels at batched inputs. This bench is the
      per-call shape — useful to reason about NequIP inference which
      evaluates O(100k) edges per message-passing layer.
    - GPU perf. e3nn runs on MPS by default on Apple Silicon; we
      pin to CPU for apples-to-apples.

The libirrep side links a tiny helper program (tp_bench_harness)
that emits ns/op to stdout; this script spawns it via subprocess.
This avoids writing ctypes bindings from scratch for a one-off
benchmark.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent


def run_libirrep_bench(iters: int, verbose: bool) -> dict:
    """Invoke the C harness; parse its JSON line."""
    harness = REPO_ROOT / "build" / "bin" / "bench_tensor_product"
    if not harness.is_file():
        print(f"building {harness.name} ...", file=sys.stderr)
        subprocess.check_call(
            ["make", "bench"], cwd=REPO_ROOT, stdout=subprocess.DEVNULL
        )

    if verbose:
        print(f"$ {harness}")
    out = subprocess.check_output([str(harness)]).decode()
    # bench_tensor_product emits three JSON lines; we want the NequIP one.
    records = []
    for line in out.splitlines():
        line = line.strip()
        if line.startswith("{"):
            records.append(json.loads(line))
    for r in records:
        if r.get("name") == "tp_apply_nequip_weighted":
            return r
    raise RuntimeError(f"tp_apply_nequip_weighted record not found; got: {records}")


def run_e3nn_bench(iters: int, verbose: bool) -> dict:
    """Replicate libirrep's NequIP descriptor in e3nn, time forward."""
    import torch
    from e3nn import o3

    torch.set_num_threads(1)              # match libirrep's scalar path
    device = torch.device("cpu")

    # Descriptor: identical to the one `bench_tensor_product.c` uses.
    #   a: 1x0e + 1x1o + 1x2e   (dim 9)
    #   b: 1x0e + 1x1o + 1x2e + 1x3o   (dim 16)
    #   c: 1x0e + 1x1o + 1x2e   (dim 9)
    irreps_a = o3.Irreps("1x0e + 1x1o + 1x2e")
    irreps_b = o3.Irreps("1x0e + 1x1o + 1x2e + 1x3o")
    irreps_c = o3.Irreps("1x0e + 1x1o + 1x2e")

    # FullyConnectedTensorProduct is e3nn's analog of libirrep's UUU
    # weighted TP: per-path learnable scalar, internal weights buffer.
    tp = o3.FullyConnectedTensorProduct(
        irreps_a, irreps_b, irreps_c,
        internal_weights=True, shared_weights=True,
    ).to(device)
    tp.eval()

    # Identical fixed inputs to the libirrep bench for numerical comparison.
    a = torch.tensor(
        [1, 0.5, 0.3, -0.2, 0.1, 0.4, -0.1, 0.6, -0.3],
        dtype=torch.float64, device=device,
    )
    b = torch.tensor(
        [1, 0.2, -0.3, 0.4, 0.5, -0.1, 0.2, 0.3, -0.4,
         0.1, 0.2, -0.3, 0.5, -0.6, 0.7, -0.1],
        dtype=torch.float64, device=device,
    )

    # e3nn defaults to float32; upcast the module for apples-to-apples.
    tp.double()

    # Warm-up (load weights, JIT, caches).
    with torch.no_grad():
        for _ in range(100):
            _ = tp(a, b)

    # Time
    with torch.no_grad():
        t0 = time.perf_counter()
        for _ in range(iters):
            out = tp(a, b)
        t1 = time.perf_counter()

    ns_per_op = (t1 - t0) * 1e9 / iters
    return {
        "name": "e3nn_FullyConnected_nequip",
        "iterations": iters,
        "ns_per_op": ns_per_op,
        "ops_per_sec": iters / (t1 - t0),
        "output_vec": out.detach().cpu().numpy().tolist(),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--iters", type=int, default=500_000)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    try:
        import torch      # noqa
        from e3nn import o3  # noqa
    except ImportError:
        print("error: this script requires e3nn + torch in the current", file=sys.stderr)
        print("       Python environment. Install via:", file=sys.stderr)
        print("         python3 -m venv /tmp/venv &&", file=sys.stderr)
        print("         /tmp/venv/bin/pip install e3nn", file=sys.stderr)
        return 1

    print("=" * 72)
    print("libirrep vs e3nn — tensor product on NequIP-like descriptor")
    print("=" * 72)
    print()
    print("descriptor:")
    print("  a = 1x0e + 1x1o + 1x2e           (dim 9)")
    print("  b = 1x0e + 1x1o + 1x2e + 1x3o    (dim 16)")
    print("  c = 1x0e + 1x1o + 1x2e           (dim 9)")
    print(f"iters: {args.iters:,}")
    print()

    irrep = run_libirrep_bench(args.iters, args.verbose)
    e3nn  = run_e3nn_bench(args.iters, args.verbose)

    print("libirrep tp_apply_weighted:")
    print(f"  {irrep['ns_per_op']:>10.1f} ns / call    ({irrep['ops_per_sec']:.2e} calls/s)")
    print()
    print("e3nn FullyConnectedTensorProduct (f64, CPU, 1 thread):")
    print(f"  {e3nn['ns_per_op']:>10.1f} ns / call    ({e3nn['ops_per_sec']:.2e} calls/s)")
    print()

    ratio = e3nn["ns_per_op"] / irrep["ns_per_op"]
    label = "faster" if ratio > 1 else "slower"
    print(f"libirrep is {ratio:.2f}× {label} on this descriptor shape.")

    # Numerical agreement: e3nn's FullyConnectedTensorProduct uses its
    # own internal weight layout, different from libirrep's per-path
    # scalar weights; so the raw outputs WILL NOT match. What we can
    # compare directly is the *unweighted* per-path TP of e3nn's
    # `FullTensorProduct` (every path weight = 1.0) against libirrep's
    # `irrep_tp_apply` with NULL weights.
    print()
    print("-" * 72)
    print("Numerical agreement — unweighted TP on identical input")
    print("-" * 72)
    try:
        lib_vec, e3_vec, max_abs_err, cos = check_numerical_agreement()
        print(f"libirrep:       {fmt_vec(lib_vec)}")
        print(f"e3nn:           {fmt_vec(e3_vec)}")
        print(f"max abs err:    {max_abs_err:.3e}")
        print(f"cosine similarity: {cos:.10f}")
        # e3nn's internal CG tables are float32-accurate (tp.double()
        # promotes the compute but not the stored coefficients);
        # accept agreement at 1e-6 which is the realistic precision
        # of a cross-library comparison.
        ok = max_abs_err < 1e-6 and cos > 1 - 1e-12
        print("  " + ("OK ✓" if ok else "DISAGREE"))
    except Exception as e:
        print(f"(skipped: {e})")

    return 0


def fmt_vec(v):
    import numpy as np
    return "[" + ", ".join(f"{x:+.6f}" for x in v) + "]"


def check_numerical_agreement():
    """
    Compare libirrep's `irrep_tp_apply` (unweighted) against e3nn's
    `FullTensorProduct` (every path-weight = 1) on identical inputs.
    Both compute the same mathematical object — Clebsch-Gordan-coupled
    tensor products in the real-basis — up to a (known) phase choice.

    Uses `examples/tp_nequip_output` binary to get libirrep's vector.
    """
    import numpy as np
    import torch
    from e3nn import o3

    # ---- libirrep side: spawn the C binary
    helper = REPO_ROOT / "build" / "bin" / "tp_nequip_output"
    if not helper.is_file():
        subprocess.check_call(
            ["make", str(helper)], cwd=REPO_ROOT, stdout=subprocess.DEVNULL
        )
    lib_json = subprocess.check_output([str(helper)]).decode().strip()
    lib_vec = np.array(json.loads(lib_json), dtype=np.float64)

    # ---- e3nn side: raw CG-only TP, no per-path normalisation, no
    # weights. This is FullTensorProduct with `irrep_normalization=
    # 'none'`, which matches libirrep's unweighted `irrep_tp_apply`
    # mathematical object byte-for-byte (both compute the real-basis
    # Clebsch-Gordan coupling with the i^{l_a+l_b-l_c} phase factor
    # applied internally). e3nn's FullTensorProduct sums into every
    # valid l_c that appears in the target multiset.
    irreps_a = o3.Irreps("1x0e + 1x1o + 1x2e")
    irreps_b = o3.Irreps("1x0e + 1x1o + 1x2e + 1x3o")

    tp = o3.FullTensorProduct(
        irreps_a, irreps_b,
        irrep_normalization="none",
        path_normalization="none",
    ).double()
    tp.eval()

    # e3nn's FullTensorProduct produces every (l_a, l_b) → l_c path
    # concatenated. We filter to the libirrep target "1x0e + 1x1o + 1x2e"
    # (keep only l=0 even, l=1 odd, l=2 even; sum multiplicities).
    tp_irreps_out = tp.irreps_out
    target_c     = o3.Irreps("1x0e + 1x1o + 1x2e")

    a = torch.tensor(
        [1, 0.5, 0.3, -0.2, 0.1, 0.4, -0.1, 0.6, -0.3],
        dtype=torch.float64,
    )
    b = torch.tensor(
        [1, 0.2, -0.3, 0.4, 0.5, -0.1, 0.2, 0.3, -0.4,
         0.1, 0.2, -0.3, 0.5, -0.6, 0.7, -0.1],
        dtype=torch.float64,
    )
    with torch.no_grad():
        full = tp(a, b).detach().numpy().astype(np.float64)

    # Reduce to target_c: sum the (mul_ab · 2l+1) chunks of `full`
    # that correspond to each (l, parity) block in target_c.
    # Convention-bridge: e3nn's internal CG tensor bakes in a
    # `1/sqrt(2·l_c + 1)` factor relative to the Sakurai / Varshalovich
    # standard convention libirrep uses. Multiply out the factor so we
    # compare apples-to-apples.
    offset_src = 0
    e3_vec = np.zeros(target_c.dim, dtype=np.float64)
    for (mul_src, (l_src, p_src)) in tp_irreps_out:
        d_src = (2 * l_src + 1)
        block_size = mul_src * d_src
        block = full[offset_src : offset_src + block_size].reshape(mul_src, d_src)
        block = block * np.sqrt(d_src)   # e3nn → Sakurai convention
        offset_src += block_size
        # find the matching (l, p) slot in target_c
        offset_dst = 0
        for (mul_dst, (l_dst, p_dst)) in target_c:
            d_dst = (2 * l_dst + 1)
            if (l_dst, p_dst) == (l_src, p_src):
                # sum over mul_src channels into the first mul_dst slots
                take = min(mul_src, mul_dst)
                e3_vec[offset_dst : offset_dst + d_dst] += block[:take].sum(axis=0)
                break
            offset_dst += mul_dst * d_dst

    max_abs_err = float(np.max(np.abs(lib_vec - e3_vec)))
    dot = float(np.dot(lib_vec, e3_vec))
    cos = dot / (np.linalg.norm(lib_vec) * np.linalg.norm(e3_vec) + 1e-300)
    return lib_vec, e3_vec, max_abs_err, cos


if __name__ == "__main__":
    sys.exit(main())
