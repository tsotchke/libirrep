#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Run all compiled benchmarks, concatenating JSON output into a timestamped file.
set -euo pipefail

BUILD_DIR="${BUILD_DIR:-build}"
BIN_DIR="${BUILD_DIR}/bin"
STAMP="$(date -u +'%Y-%m-%dT%H-%M-%SZ')"
TRIPLE="$(uname -s | tr '[:upper:]' '[:lower:]')-$(uname -m)"
OUT_DIR="benchmarks/results/${STAMP}"
mkdir -p "${OUT_DIR}"
OUT_FILE="${OUT_DIR}/${TRIPLE}.json"

# CPU fingerprint — different cores at the same arch (e.g. M1 Max vs M3 Pro)
# produce enough variance to trip the perf_compare gate. Baking the model
# into the results JSON lets the diff tool flag incompatible baselines.
case "$(uname -s)" in
    Darwin) CPU_MODEL=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo "unknown") ;;
    Linux)  CPU_MODEL=$(awk -F': ' '/^model name/ {print $2; exit}' /proc/cpuinfo 2>/dev/null || echo "unknown") ;;
    *)      CPU_MODEL="unknown" ;;
esac

echo "# benchmark results → ${OUT_FILE}"
echo "# cpu model        → ${CPU_MODEL}"

# Each bench binary emits one or more `{"name":..., ...}` records per line.
# We collect them all, drop blank lines, and join with commas so the final
# file is a valid JSON array (consumable by scripts/perf_compare.sh).
#
# Iterate over benchmarks/bench_*.c sources rather than build/bin/bench_*
# binaries so we pick up exactly the benchmark suite (not example-folder
# demos that happen to share the `bench_` prefix and emit human-readable
# banners on stdout — those would otherwise corrupt the JSON aggregate).
TMP="$(mktemp)"
for src in benchmarks/bench_*.c; do
    [ -f "${src}" ] || continue
    name=$(basename "${src}" .c)
    bin="${BIN_DIR}/${name}"
    [ -x "${bin}" ] || continue
    "${bin}" >> "${TMP}"
done

{
    echo "["
    # First record is a metadata sentinel — perf_compare skips records whose
    # name starts with "_meta" when diffing, but still reads them to extract
    # cpu_model / stamp / triple for the compatibility check.
    printf '{"name":"_meta","cpu_model":"%s","triple":"%s","stamp":"%s"}' \
        "${CPU_MODEL}" "${TRIPLE}" "${STAMP}"
    # Remaining records, comma-prefixed.
    awk 'NF > 0' "${TMP}" | awk '{ print ","; print $0; }'
    echo
    echo "]"
} > "${OUT_FILE}"

rm -f "${TMP}"

echo "# done"
