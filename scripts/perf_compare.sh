#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
#
# Compare two benchmark JSON files and flag regressions past a threshold.
# Benchmark files are arrays of `{ "name": string, "ns_per_op": number, ... }`
# produced by benchmarks/run_benchmarks.sh.
#
# Usage:
#   scripts/perf_compare.sh baseline.json candidate.json [pct_threshold]
#
# Exit status:
#   0 — no regressions past threshold
#   1 — at least one benchmark regressed; offenders listed on stdout
#   2 — bad input (missing file, malformed JSON, no overlapping names)
set -euo pipefail

BASELINE="${1:?usage: scripts/perf_compare.sh baseline.json candidate.json [pct] [min_ns]}"
CANDIDATE="${2:?usage: scripts/perf_compare.sh baseline.json candidate.json [pct] [min_ns]}"
PCT="${3:-5}"            # default: fail on ≥5% regression
MIN_ABS_NS="${4:-1.0}"   # absolute floor: ignore deltas below this many ns — keeps
                         # system noise on sub-10-ns benchmarks (rbf_bessel, tp_apply_empty,
                         # ..) from tripping the gate. Tune per-machine if needed.

if [ ! -f "$BASELINE" ]  || [ ! -f "$CANDIDATE" ]; then
    echo "perf_compare: missing input" >&2
    exit 2
fi

command -v python3 >/dev/null 2>&1 || {
    echo "perf_compare: python3 required for JSON parsing" >&2
    exit 2
}

python3 - "$BASELINE" "$CANDIDATE" "$PCT" "$MIN_ABS_NS" <<'PY'
import json, sys

baseline_path, candidate_path, pct_str, min_abs_str = sys.argv[1:5]
threshold    = float(pct_str)
min_abs_ns   = float(min_abs_str)

def load(p):
    with open(p) as f:
        data = json.load(f)
    if isinstance(data, dict) and "benchmarks" in data:
        data = data["benchmarks"]
    if not isinstance(data, list):
        sys.stderr.write(f"perf_compare: {p}: not a list of benchmarks\n")
        sys.exit(2)
    return { b["name"]: b for b in data if "name" in b }

b_all = load(baseline_path)
c_all = load(candidate_path)

# Extract and compare the _meta sentinel so mismatched hardware produces a
# visible warning. Benchmarks whose name starts with "_" are skipped in the
# diff itself — they're metadata, not measurements.
meta_b = b_all.get("_meta", {})
meta_c = c_all.get("_meta", {})
cpu_b  = meta_b.get("cpu_model", "unknown")
cpu_c  = meta_c.get("cpu_model", "unknown")
print(f"baseline cpu : {cpu_b}")
print(f"candidate cpu: {cpu_c}")
if cpu_b != "unknown" and cpu_c != "unknown" and cpu_b != cpu_c:
    print(f"WARNING: hardware mismatch — comparison may generate false positives")
print()

b = { k: v for k, v in b_all.items() if not k.startswith("_") }
c = { k: v for k, v in c_all.items() if not k.startswith("_") }
common = sorted(set(b) & set(c))

if not common:
    sys.stderr.write("perf_compare: no overlapping benchmark names\n")
    sys.exit(2)

regressed = []
improved  = []
unchanged = []

print(f"{'benchmark':<40}  {'baseline':>12}  {'candidate':>12}  {'delta':>8}")
print("-" * 80)
for name in common:
    base_ns = float(b[name].get("ns_per_op", 0.0))
    cand_ns = float(c[name].get("ns_per_op", 0.0))
    if base_ns <= 0:
        continue
    delta_ns  = cand_ns - base_ns
    delta_pct = 100.0 * delta_ns / base_ns
    marker = " "
    # Suppress noise on very fast benchmarks: the gate requires both the
    # relative delta exceed `threshold` AND the absolute delta exceed
    # `min_abs_ns`. Printed delta is unchanged.
    over_pct = delta_pct >  threshold and abs(delta_ns) >= min_abs_ns
    und_pct  = delta_pct < -threshold and abs(delta_ns) >= min_abs_ns
    if   over_pct: marker = "!"; regressed.append((name, delta_pct))
    elif und_pct:  marker = "+"; improved.append((name, delta_pct))
    else: unchanged.append(name)
    print(f"{marker} {name:<38}  {base_ns:>12.2f}  {cand_ns:>12.2f}  {delta_pct:>7.1f}%")

print("-" * 80)
print(f"regressed: {len(regressed)}    improved: {len(improved)}    unchanged: {len(unchanged)}")
print(f"(threshold: ±{threshold:.1f}%   min abs delta: {min_abs_ns:.2f} ns)")

if regressed:
    print()
    print("regressions past threshold:")
    for name, pct in regressed:
        print(f"  {name}: {pct:+.1f}%")
    sys.exit(1)
sys.exit(0)
PY
