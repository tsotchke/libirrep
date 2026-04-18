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

echo "# benchmark results → ${OUT_FILE}"
{
    echo "["
    first=1
    for b in "${BIN_DIR}"/bench_*; do
        [ -x "${b}" ] || continue
        if [ ${first} -eq 0 ]; then echo ","; fi
        "${b}"
        first=0
    done
    echo "]"
} > "${OUT_FILE}"

echo "# done"
